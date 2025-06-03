#!/usr/bin/bash
# get script arguments
READ1=$1
READ2=$2
SAMPLE_ID=$3
CORE=$4
CTATDB=$5

run_starfusion() {
  local FASTQ_R1_FILE=$1
  local FASTQ_R2_FILE=$2
  local SAMPLE_ID=$3
  local CORES=$4
  local DB=$5

  STAR-Fusion --genome_lib_dir /path/to/your/CTAT_resource_lib --left_fq $FASTQ_R1_FILE --right_fq $FASTQ_R2_FILE --output_dir "${SAMPLE_ID}-STARFusion-out"
}

# measure execution time
STARTTIME=$(date +%s)
if run_starfusion "${READ1}" "${READ2}" "${SAMPLE_ID}" "${CORE}" "${CTATDB}"; then
    ENDTIME=$(date +%s)
    ELAP=$(( ENDTIME - STARTTIME ))
    echo "STARFusion run completed successfully. Time taken: ${ELAP}. Check log file for run details."
else
    echo "Something went wrong during STARFusion run. Check log file."
fi