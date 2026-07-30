#!/usr/bin/env python3
"""
build-checkpoint-manifest.py
 
Builds a checkpoint manifest for re-entering NeonDisco at the pvacFuse stage,
by listing the pipeline's own published S3 output structure
 
This works specifically because these two upstream processes already
publishDir with sampleName baked into the path:
 
    KEEP_VALIDATED_FUSIONS_PYENV
        -> ${outputDir}/${sampleName}/POST-IN-SILICO-VALIDATION-TSV-out/
               validated-agfusion-dirs/                       (agfusion_dir)
 
    REFORMAT_AND_COLLATE_HLA_RESULTS_PYENV
        -> ${outputDir}/Cohortwide-HLA-typing-OUT/
               Cohortwide_HLA_types.tsv                       (hla_tsv, global)
 
Which agfusion source directory is correct depends on params.recurrentFusionsNeoPredOnly
AS ACTUALLY RESOLVED FOR THE RUN BEING RECOVERED (profile params override base.config):
    - false (e.g. -profile dualNeoPredMode, sampleHLANeoPredMode) ->
          KEEP_VALIDATED_FUSIONS_PYENV.out.validatedAgfusionDir      [default below]
    - true  (e.g. -profile sharedHLANeoPredMode, or base.config default with
             no mode profile overriding it) ->
          FILTER_SAMPLE_LEVEL_VALIDATED_FUSIONS_FOR_RECURRENT_PYENV.out.validatedRecurrentAgfusionDir
          i.e. ${sampleName}/RECURRENT-VALIDATED-FUSIONS-out/FI-validated-recurrent-agfusion-dirs/
    Use --recurrent-only to switch.
 
Usage:
    python3 build-checkpoint-manifest.py \\
        --output-dir s3://crmy-gb-main/neondisco/main-analyses/MSA/running/NeonDisco-outputs/MYBRCA-n990-dataset-awsbatch-dualMode-ND-v0-2-2-RERUN-FINAL-v5/ \\
        --manifest-out checkpoint_manifest.tsv
        # add --recurrent-only if the run used sharedHLANeoPredMode (or otherwise had recurrentFusionsNeoPredOnly=true)
"""
 
import argparse
import csv
import sys
from urllib.parse import urlparse
 
try:
    import boto3
except ImportError:
    sys.exit("boto3 is required: pip install boto3 --break-system-packages")
 
AGFUSION_SUBPATH_ALL_VALIDATED = "POST-IN-SILICO-VALIDATION-TSV-out/validated-agfusion-dirs/"
AGFUSION_SUBPATH_RECURRENT_ONLY = "RECURRENT-VALIDATED-FUSIONS-out/FI-validated-recurrent-agfusion-dirs/"
HLA_TSV_SUBPATH = "Cohortwide-HLA-typing-OUT/Cohortwide_HLA_types.tsv"
 
# Top-level prefixes under outputDir that are NOT sample directories.
NON_SAMPLE_PREFIXES = {"Cohortwide-HLA-typing-OUT", "Cohortwide-Fusions-OUT"}
 
 
def parse_s3_uri(uri):
    p = urlparse(uri)
    return p.netloc, p.path.lstrip("/")
 
 
def list_sample_names(s3, output_dir_uri):
    bucket, prefix = parse_s3_uri(output_dir_uri)
    if prefix and not prefix.endswith("/"):
        prefix += "/"
 
    paginator = s3.get_paginator("list_objects_v2")
    samples = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix, Delimiter="/"):
        for cp in page.get("CommonPrefixes", []):
            name = cp["Prefix"][len(prefix):].rstrip("/")
            if name and name not in NON_SAMPLE_PREFIXES:
                samples.append(name)
    return sorted(samples)
 
 
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output-dir", required=True, help="s3:// URI of params.outputDir for the run being recovered")
    ap.add_argument("--recurrent-only", action="store_true",
                     help="Set if the run being recovered had recurrentFusionsNeoPredOnly=true "
                          "(e.g. -profile sharedHLANeoPredMode). Default assumes false "
                          "(dualNeoPredMode / sampleHLANeoPredMode), matching KEEP_VALIDATED_FUSIONS_PYENV output.")
    ap.add_argument("--manifest-out", default="checkpoint_manifest.tsv")
    ap.add_argument("--profile", default=None, help="AWS profile name (unrelated to Nextflow -profile)")
    ap.add_argument("--region", default=None)
    args = ap.parse_args()
 
    session = boto3.Session(profile_name=args.profile, region_name=args.region)
    s3 = session.client("s3")
    sts = session.client("sts")
    identity = sts.get_caller_identity()
    resolved_region = session.region_name or "unset"
    print(f"Using AWS identity: {identity['Arn']} (account {identity['Account']}), " f"region: {resolved_region}", file=sys.stderr)
    if resolved_region != "ap-southeast-5":
        print(f"WARNING: resolved region '{resolved_region}' does not match " f"NeonDisco's configured aws.region 'ap-southeast-5'", file=sys.stderr)


    output_dir = args.output_dir.rstrip("/")
    print(f"Listing sample directories under {output_dir} ...", file=sys.stderr)
    samples = list_sample_names(s3, output_dir)
    print(f"Found {len(samples)} candidate sample directories.", file=sys.stderr)
 
    if not samples:
        sys.exit("No sample directories found -- check --output-dir is correct and matches the run you're recovering.")
 
    agfusion_subpath = AGFUSION_SUBPATH_RECURRENT_ONLY if args.recurrent_only else AGFUSION_SUBPATH_ALL_VALIDATED
    print(f"Using agfusion source: {agfusion_subpath} "
          f"({'recurrent-only' if args.recurrent_only else 'all validated fusions'})", file=sys.stderr)
 
    hla_tsv_uri = f"{output_dir}/{HLA_TSV_SUBPATH}"
 
    with open(args.manifest_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_name", "file_role", "s3_uri"])
        for sample in samples:
            agfusion_uri = f"{output_dir}/{sample}/{agfusion_subpath}"
            w.writerow([sample, "agfusion_dir", agfusion_uri])
            w.writerow([sample, "hla_tsv", hla_tsv_uri])
 
    print(f"Wrote {args.manifest_out}: {len(samples)} samples x 2 roles "
          f"({len(samples) * 2} rows, hla_tsv deduplicated during validation). "
          f"metaDataDir intentionally excluded (local path, not an S3 checkpoint concern).",
          file=sys.stderr)
    print("Next: python3 validate_checkpoint_manifest.py --manifest "
          f"{args.manifest_out} --report validation_report.tsv --ready-manifest ready_manifest.tsv",
          file=sys.stderr)
 
 
if __name__ == "__main__":
    main()
