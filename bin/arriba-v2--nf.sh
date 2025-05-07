#!/usr/bin/bash

# Set env variables
ARRIBA_PKG="/opt/conda/var/lib/arriba"

BAM=$1
SAMPLE_ID=$2
ARRIBA_DB=$3
CORE=$4

ARRIBA_ASSEMBLY="${ARRIBA_DB}/GRCh38viral.fa"
ARRIBA_ANNOT="${ARRIBA_DB}/ENSEMBL113.gtf"

echo "Environment variables set! Listing fastq files..."
echo "BAM File: $BAM"
echo "Sample ID: $SAMPLE_ID"
echo "Number of cores: $CORE"

run_arriba() {
  local BAM=$1
  local SAMPLE_ID=$2
  local ARRIBA_ASSEMBLY=$3
  local ARRIBA_ANNOT=$4
  local ARRIBA_PKG=$5
  local CORE=$6

  echo "Running Arriba..."
  arriba \
      -x "${BAM}" \
      -o "${SAMPLE_ID}-arriba-fusions.tsv" \
      -O "${SAMPLE_ID}-fusions.discarded.tsv" \
      -a "${ARRIBA_ASSEMBLY}" \
      -g "${ARRIBA_ANNOT}" \
      -b "${ARRIBA_PKG}/blacklist_hg38_GRCh38_v2.3.0.tsv.gz" \
      -k "${ARRIBA_PKG}/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz" \
      -t "${ARRIBA_PKG}/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz" \
      -p "${ARRIBA_PKG}/protein_domains_hg38_GRCh38_v2.3.0.gff3" 2>&1 | tee "arriba-run-nf.log-$(date +%Y%m%d_%H-%M-%S).txt"
}

echo "Starting Arriba..."
# measure execution time
STARTTIME=$(date +%s)
if run_arriba "${BAM}" "${SAMPLE_ID}" "${ARRIBA_ASSEMBLY}" "${ARRIBA_ANNOT}" "${ARRIBA_PKG}" "${CORE}"; then
    ENDTIME=$(date +%s)
    ELAP=$(( ENDTIME - STARTTIME ))
    echo "Arriba run of ${SAMPLE_ID} completed successfully. Time taken: ${ELAP} seconds. Check log file for run details."
else
    echo "Something went wrong during Arriba run of ${SAMPLE_ID}. Check log file."
fi
