#!/usr/bin/bash
# get script arguments
READ1=$1
READ2=$2
SAMPLE_ID=$3
CORE=$4

run_starfusion() {
  local FASTQ_R1_FILE=$1
  local FASTQ_R2_FILE=$2
  local DB=$3
  local OUTDIR=$4
  local CORES=$5

  fusioncatcher -d ${DB} --input ${FASTQ_R1_FILE},${FASTQ_R2_FILE} --output ${OUTDIR} -p $5
}

# Set env variables
export DB="/home/app/libs"
export OUTDIR_PREFIX="/home/app/nf_work"
export OUTDIR="${OUTDIR_PREFIX}/${SAMPLE_ID}"
echo "Environment variables set and exported!"

# measure execution time
STARTTIME=$(date +%s)
if run_fuscat "${READ1}" "${READ2}" "${DB}" "${OUTDIR_PREFIX}" "${CORE}"; then
    ENDTIME=$(date +%s)
    ELAP=$(( ENDTIME - STARTTIME ))
    echo "FusionCatcher run completed successfully. Time taken: ${ELAP}. Check log file for run details."
else
    echo "Something went wrong during FusionCatcher run. Check log file."
fi