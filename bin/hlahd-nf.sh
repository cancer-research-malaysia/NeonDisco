#!/usr/bin/bash 

echo $(id)

OUTDIR="/work/nf_work"
HLAHD_DIR="/tmp/hlahd.1.7.1"

SAMPLE_ID=$1
READ1=$2
READ2=$3
CORES=$4

echo "Running HLA-HD for sample: ${SAMPLE_ID}"
echo "FastQ input files: ${READ1} & ${READ2}"
echo "Number of cores for HLA-HD specified: ${CORES}"

mkdir -p "${OUTDIR}/hla-hd-out"

# measure execution time
STARTTIME=$(date +%s)

# run HLA-HD (discard reads <30bp)
if hlahd.sh -t ${CORES} -m 30 -f ${HLAHD_DIR}/freq_data/ "${READ1}" "${READ2}" "${HLAHD_DIR}/HLA_gene.split.3.50.0.txt" "${HLAHD_DIR}/dictionary/" "${SAMPLE_ID}_HLAHD" "${OUTDIR}" 2>&1 | tee "${SAMPLE_ID}-HLAHD.log"; then
    echo "HLA-HD run successfully for sample ${SAMPLE_ID}." >> "${SAMPLE_ID}-HLAHD.log"
    ENDTIME=$(date +%s)
    ELAP=$(( ENDTIME - STARTTIME ))
    echo "Time taken: ${ELAP}. Check log file for run details." >> "${SAMPLE_ID}-HLAHD.log"
else
    echo "HLA-HD run failed for sample ${SAMPLE_ID}."
    exit 1
fi