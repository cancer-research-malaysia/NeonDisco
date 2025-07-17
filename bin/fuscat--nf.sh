#!/usr/bin/env bash

# get script arguments
READ1=$1
READ2=$2
SAMPLE_ID="${READ1%%_*}"
DB=$3
CORE=$4

run_fuscat() {
  local FASTQ_R1_FILE=$1
  local FASTQ_R2_FILE=$2
  local DB=$3
  local CORES=$4
  local OUTDIR=$5

  fusioncatcher -d ${DB} --input ${FASTQ_R1_FILE},${FASTQ_R2_FILE} --output ${OUTDIR} -p ${CORES} 2>&1 | tee "${OUTDIR}/fuscat-run.log-$(date +%Y%m%d_%H-%M-%S).txt"
}

# Set variables
OUTDIR="${SAMPLE_ID}"

# measure execution time
STARTTIME=$(date +%s)
if run_fuscat "${READ1}" "${READ2}" "${DB}" "${CORE}" "${OUTDIR}"; then
    ENDTIME=$(date +%s)
    ELAP=$(( ENDTIME - STARTTIME ))
    echo "FusionCatcher run completed successfully. Time taken: ${ELAP}. Check log file for run details."
else
    echo "Something went wrong during FusionCatcher run. Check log file."
fi
