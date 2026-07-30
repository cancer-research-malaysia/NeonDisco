#!/usr/bin/env python3
"""
validate-checkpoint-manifest.py

Batch-validates that S3 objects required to re-enter a Nextflow pipeline
mid-DAG (e.g. rerunning PVACFUSE from a saved AGFusion checkpoint) actually
exist, are non-empty, and are readable without a Glacier/Archive restore.

INPUT MANIFEST (TSV, header required):
    sample_name    file_role         s3_uri
    MB0001         agfusion_dir      s3://bucket/work/de/6734ac.../MB0001_validated_agfusion
    MB0001         hla_tsv           s3://bucket/.../cohort_wide_hla_list.tsv
    MB0001         metadata_dir      s3://bucket/refdata/neondisco-metadata

Rows sharing the same file_role + s3_uri (e.g. hla_tsv, metadata_dir, which
are global/shared across all samples) are deduplicated automatically so
they're only checked once.

`file_role` doesn't need to match any fixed vocabulary -- it's just a label
carried through to the report so you know which take: input a failure
belongs to. Directory-type roles (any s3_uri ending in '/', or a role name
containing "dir") are checked via list_objects_v2; everything else is
checked via head_object.

OUTPUT:
    - <report>.tsv        every row checked, with status + detail
    - <ready_manifest>.tsv  only sample_names where every required role is OK
                             (this is what you feed into neondisco's
                             --startFrom seed channel)

Usage:
    pip install boto3 --break-system-packages
    python3 validate-checkpoint-manifest.py \\
        --manifest checkpoint_manifest.tsv \\
        --report validation_report.tsv \\
        --ready-manifest ready_manifest.tsv \\
        --max-workers 32
"""

import argparse
import csv
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import urlparse

try:
    import boto3
    from botocore.config import Config
    from botocore.exceptions import ClientError
except ImportError:
    sys.exit("boto3 is required: pip install boto3 --break-system-packages")

# Storage classes that require a restore request before they're readable.
# GLACIER_IR is deliberately excluded -- despite the name, Glacier Instant
# Retrieval objects are readable directly, no restore needed.
FROZEN_STORAGE_CLASSES = {"GLACIER", "DEEP_ARCHIVE"}

# Intelligent-Tiering objects that have aged into an archive tier report
# this via ArchiveStatus on head_object, independent of StorageClass.
FROZEN_ARCHIVE_STATUSES = {"ARCHIVE_ACCESS", "DEEP_ARCHIVE_ACCESS"}


def parse_s3_uri(uri):
    p = urlparse(uri)
    if p.scheme != "s3":
        raise ValueError(f"not an s3:// URI: {uri}")
    return p.netloc, p.path.lstrip("/")


def is_dir_like(role, uri):
    return uri.endswith("/") or "dir" in role.lower()


def check_object(s3, uri):
    """Check a single S3 object (file). Returns (status, detail, storage_class)."""
    bucket, key = parse_s3_uri(uri)
    try:
        resp = s3.head_object(Bucket=bucket, Key=key)
    except ClientError as e:
        code = e.response.get("Error", {}).get("Code", "")
        if code in ("404", "NoSuchKey", "NotFound"):
            return "MISSING", "object not found", None
        return "ERROR", f"{code}: {e}", None

    size = resp.get("ContentLength", 0)
    storage_class = resp.get("StorageClass", "STANDARD")
    archive_status = resp.get("ArchiveStatus")
    restore = resp.get("Restore")  # e.g. 'ongoing-request="false", expiry-date="..."'

    if archive_status in FROZEN_ARCHIVE_STATUSES:
        return "FROZEN", f"intelligent-tiering archive status: {archive_status}", storage_class

    if storage_class in FROZEN_STORAGE_CLASSES:
        if restore and 'ongoing-request="false"' in restore:
            return "OK", f"restored copy available ({restore})", storage_class
        return "FROZEN", f"storage class {storage_class}, not restored", storage_class

    if size == 0:
        return "EMPTY", "zero-byte object", storage_class

    return "OK", f"{size} bytes", storage_class


def check_prefix(s3, uri):
    """Check a directory-like S3 prefix. Returns (status, detail, storage_class)."""
    bucket, key = parse_s3_uri(uri)
    if key and not key.endswith("/"):
        key += "/"

    paginator = s3.get_paginator("list_objects_v2")
    objects = []
    try:
        for page in paginator.paginate(Bucket=bucket, Prefix=key):
            objects.extend(page.get("Contents", []))
    except ClientError as e:
        code = e.response.get("Error", {}).get("Code", "")
        return "ERROR", f"{code}: {e}", None

    if not objects:
        return "MISSING", "prefix has no objects", None

    frozen = [o for o in objects if o.get("StorageClass") in FROZEN_STORAGE_CLASSES]
    empty = [o for o in objects if o.get("Size", 0) == 0]
    zero_size_total = sum(o.get("Size", 0) for o in objects) == 0

    if frozen:
        classes = sorted({o["StorageClass"] for o in frozen})
        return "FROZEN", f"{len(frozen)}/{len(objects)} objects in {classes}", ",".join(classes)

    if zero_size_total:
        return "EMPTY", f"all {len(objects)} objects are zero-byte", None

    worst_class = sorted({o.get("StorageClass", "STANDARD") for o in objects})
    return "OK", f"{len(objects)} objects, {sum(o.get('Size', 0) for o in objects)} bytes total", ",".join(worst_class)


def check_row(s3, role, uri):
    if is_dir_like(role, uri):
        return check_prefix(s3, uri)
    return check_object(s3, uri)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="TSV: sample_name, file_role, s3_uri")
    ap.add_argument("--report", default="validation_report.tsv")
    ap.add_argument("--ready-manifest", default="ready_manifest.tsv")
    ap.add_argument("--max-workers", type=int, default=32)
    ap.add_argument("--profile", default=None, help="AWS profile name (optional)")
    ap.add_argument("--region", default=None, help="AWS region (optional)")
    args = ap.parse_args()

    session = boto3.Session(profile_name=args.profile, region_name=args.region)
    s3 = session.client("s3", config=Config(max_pool_connections=args.max_workers))

    with open(args.manifest, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    required_cols = {"sample_name", "file_role", "s3_uri"}
    if not rows or not required_cols.issubset(rows[0].keys()):
        sys.exit(f"manifest must be a TSV with columns: {sorted(required_cols)}")

    # Dedupe identical (file_role, s3_uri) pairs -- shared/global inputs
    # (e.g. hla_tsv, metadata_dir) would otherwise get checked once per sample.
    unique_checks = {}
    row_refs = defaultdict(list)
    for row in rows:
        key = (row["file_role"], row["s3_uri"])
        unique_checks[key] = row
        row_refs[key].append(row["sample_name"])

    print(f"{len(rows)} manifest rows -> {len(unique_checks)} unique S3 locations to check "
          f"(shared inputs deduplicated)", file=sys.stderr)

    results = {}
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futures = {
            pool.submit(check_row, s3, role, uri): (role, uri)
            for (role, uri) in unique_checks
        }
        done = 0
        for fut in as_completed(futures):
            role, uri = futures[fut]
            try:
                status, detail, storage_class = fut.result()
            except Exception as e:
                status, detail, storage_class = "ERROR", str(e), None
            results[(role, uri)] = (status, detail, storage_class)
            done += 1
            if done % 100 == 0 or done == len(unique_checks):
                print(f"  checked {done}/{len(unique_checks)}", file=sys.stderr)

    # Write full report (one line per original manifest row, dedup expanded back out)
    with open(args.report, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_name", "file_role", "s3_uri", "status", "detail", "storage_class"])
        for row in rows:
            key = (row["file_role"], row["s3_uri"])
            status, detail, storage_class = results[key]
            w.writerow([row["sample_name"], row["file_role"], row["s3_uri"], status, detail, storage_class])

    # A sample is "ready" only if every role it needs came back OK
    sample_roles_needed = defaultdict(set)
    sample_roles_ok = defaultdict(set)
    for row in rows:
        key = (row["file_role"], row["s3_uri"])
        status = results[key][0]
        sample_roles_needed[row["sample_name"]].add(row["file_role"])
        if status == "OK":
            sample_roles_ok[row["sample_name"]].add(row["file_role"])

    ready_samples = {
        s for s in sample_roles_needed
        if sample_roles_needed[s] == sample_roles_ok[s]
    }

    with open(args.ready_manifest, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_name", "file_role", "s3_uri"])
        for row in rows:
            if row["sample_name"] in ready_samples:
                w.writerow([row["sample_name"], row["file_role"], row["s3_uri"]])

    # Summary
    status_counts = defaultdict(int)
    for status, _, _ in results.values():
        status_counts[status] += 1

    print("\n=== summary (unique locations) ===", file=sys.stderr)
    for status in ("OK", "MISSING", "FROZEN", "EMPTY", "ERROR"):
        if status_counts.get(status):
            print(f"  {status}: {status_counts[status]}", file=sys.stderr)

    all_samples = set(sample_roles_needed)
    not_ready = all_samples - ready_samples
    print(f"\n{len(ready_samples)}/{len(all_samples)} samples fully ready to recover from checkpoint.",
          file=sys.stderr)
    if not_ready:
        preview = ", ".join(sorted(not_ready)[:10])
        more = f" (+{len(not_ready) - 10} more)" if len(not_ready) > 10 else ""
        print(f"{len(not_ready)} samples need recompute: {preview}{more}", file=sys.stderr)

    print(f"\nWrote {args.report} (full detail) and {args.ready_manifest} "
          f"(seed manifest for --startFrom).", file=sys.stderr)


if __name__ == "__main__":
    main()
