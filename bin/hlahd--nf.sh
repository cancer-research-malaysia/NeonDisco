#!/usr/bin/bash 

HLAHD_DIR="/tmp/hlahd.1.7.1"

SAMPLE_ID=$1
READ1=$2
READ2=$3
CORES=$4

echo "Running HLA-HD for sample: ${SAMPLE_ID}"
echo "FastQ input files: ${READ1} & ${READ2}"
echo "Number of cores for HLA-HD specified: ${CORES}"

# measure execution time
STARTTIME=$(date +%s)

# run HLA-HD (discard reads <30bp)
mkdir -p HLAHD-out/
if hlahd.sh -t ${CORES} -m 30 -f ${HLAHD_DIR}/freq_data/ "${READ1}" "${READ2}" "${HLAHD_DIR}/HLA_gene.split.3.50.0.txt" "${HLAHD_DIR}/dictionary/" "${SAMPLE_ID}" "HLAHD-out"; then
    echo "HLA-HD run successfully for sample ${SAMPLE_ID}."
    ENDTIME=$(date +%s)
    ELAP=$(( ENDTIME - STARTTIME ))
    echo "Time taken: ${ELAP}. Check log file for run details."
else
    echo "HLA-HD run failed for sample ${SAMPLE_ID}."
    exit 1
fi
