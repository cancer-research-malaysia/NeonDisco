#!/usr/bin/bash

# Set env variables
CTAT_LIB="/work/libs"
ARR_OUTDIR="/work/nf_work"
ARRIBA_PKG="/opt/conda/var/lib/arriba"

READ1=$1
READ2=$2
SAMPLE_ID="${READ1%%_*}"
CORE=$3

echo "Environment variables set! Listing fastq files..."
echo "READ1 File: $READ1"
echo "READ2 File: $READ2"
echo "Sample ID: $SAMPLE_ID"
echo "Number of cores: $CORE"

###### TESTING ########
# touch "${ARR_OUTDIR}/${SAMPLE_ID}-sample.txt"
# echo "Sample file created!"
# echo "${READ1} & ${READ2}: The sample ID is ${SAMPLE_ID}" > "${ARR_OUTDIR}/${SAMPLE_ID}-sample.txt"
# echo "Done!"
########################

run_star_and_arriba() {
  local READ1=$1
  local READ2=$2
  local SAMPLE_ID=$3
  local CTAT_LIB=$4
  local ARRIBA_PKG=$5
  local ARR_OUTDIR=$6
  local CORES=$7

  if STAR --runThreadN "${CORES}" \
  --genomeDir "${CTAT_LIB}/ref_genome.fa.star.idx" \
  --genomeLoad LoadAndRemove \
  --readFilesIn "${READ1}" "${READ2}" \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped Within \
  --outBAMcompression 0 \
  --outFilterMultimapNmax 50 \
  --peOverlapNbasesMin 10 \
  --alignSplicedMateMapLminOverLmate 0.5 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --chimSegmentMin 10 \
  --chimOutType WithinBAM HardClip \
  --chimJunctionOverhangMin 10 \
  --chimScoreDropMax 30 \
  --chimScoreJunctionNonGTAG 0 \
  --chimScoreSeparation 1 \
  --chimSegmentReadGapMax 3 \
  --chimMultimapNmax 50 \
  --outFileNamePrefix "${ARR_OUTDIR}/${SAMPLE_ID}-STAR_"; then
  echo "STAR run of ${SAMPLE_ID} completed successfully. Running Arriba..."
  arriba \
      -x "${ARR_OUTDIR}/${SAMPLE_ID}-STAR_Aligned.out.bam" \
      -o "${ARR_OUTDIR}/${SAMPLE_ID}-arriba-fusions.tsv" \
      -O "${ARR_OUTDIR}/${SAMPLE_ID}-fusions.discarded.tsv" \
      -a "${CTAT_LIB}/ref_genome.fa" \
      -g "${CTAT_LIB}/ref_annot.gtf" \
      -b "${ARRIBA_PKG}/blacklist_hg38_GRCh38_v2.3.0.tsv.gz" \
      -k "${ARRIBA_PKG}/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz" \
      -t "${ARRIBA_PKG}/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz" \
      -p "${ARRIBA_PKG}/protein_domains_hg38_GRCh38_v2.3.0.gff3" 2>&1 | tee "${ARR_OUTDIR}/arriba-run-nf.log-$(date +%Y%m%d_%H-%M-%S).txt"
  fi
}

echo "Starting STAR and then Arriba..."
# measure execution time
STARTTIME=$(date +%s)
if run_star_and_arriba "${READ1}" "${READ2}" "${SAMPLE_ID}" "${CTAT_LIB}" "${ARRIBA_PKG}" "${ARR_OUTDIR}" "${CORE}"; then
    ENDTIME=$(date +%s)
    ELAP=$(( ENDTIME - STARTTIME ))
    echo "Arriba run of ${SAMPLE_ID} completed successfully. Time taken: ${ELAP} seconds. Check log file for run details."
else
    echo "Something went wrong during Arriba run of ${SAMPLE_ID}. Check log file."
fi
