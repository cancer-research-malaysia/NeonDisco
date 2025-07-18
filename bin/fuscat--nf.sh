#!/usr/bin/env bash

# get script arguments
READ1=$1
READ2=$2
DB=$3
CORE=$4
SAMPLE_ID=$5

run_fuscat() {
  local FASTQ_R1_FILE=$1
  local FASTQ_R2_FILE=$2
  local DB=$3
  local CORES=$4
  local SAMPLE=$5

  fusioncatcher -d ${DB} --input ${FASTQ_R1_FILE},${FASTQ_R2_FILE} --output ${SAMPLE} -p ${CORES} 2>&1 | tee "${SAMPLE}-fuscat-run.log-$(date +%Y%m%d_%H-%M-%S).txt"
}

# measure execution time
STARTTIME=$(date +%s)
if run_fuscat "${READ1}" "${READ2}" "${DB}" "${CORE}" "${SAMPLE_ID}"; then
    ENDTIME=$(date +%s)
    ELAP=$(( ENDTIME - STARTTIME ))
    echo "FusionCatcher run completed successfully. Time taken: ${ELAP}. Check log file for run details."
else
    echo "Something went wrong during FusionCatcher run. Check log file."
fi
