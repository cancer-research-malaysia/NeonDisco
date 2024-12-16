#!/usr/bin/bash 

echo $(id)

OUTDIR="/work/nf_work"
HLAHD_DIR="/tmp/hlahd.1.7.0"

SAMPLE_ID=$1
READ1=$2
READ2=$3
CORES=$4

echo "Running HLA-HD for sample: ${SAMPLE_ID}"
echo "FastQ input files: ${READ1} & ${READ2}"
echo "Number of cores for HLA-HD specified: ${CORES}"

mkdir -p "${OUTDIR}/logs/"
mkdir -p "${OUTDIR}/hla-hd-out"

# measure execution time
STARTTIME=$(date +%s)

# run HLA-HD (discard reads <30bp)
# if hlahd.sh -t ${CORES} -m 30 -f ${HLAHD_DIR}/freq_data/ "${OUTDIR}/merged-${}-MHC-WES-b2f-R1.fq" "${OUTDIR}/merged-${SAMPLE_NAME}-MHC-WES-b2f-R2.fq" "${HLAHD_DIR}/HLA_gene.split.3.50.0.txt" "${HLAHD_DIR}/dictionary/" "${SAMPLE_NAME}_WES-MHC-bam" "${OUTDIR}/hla-hd-out" 2>&1 | tee "${OUTDIR}/logs/${SAMPLE_NAME}-HLAHD-run.log-$(date +%Y%m%d_%H-%M-%S).txt"; then
#     echo "HLA-HD run successfully for sample ${SAMPLE_NAME}."
#     ENDTIME=$(date +%s)
#     ELAP=$(( ENDTIME - STARTTIME ))
#     echo "Time taken: ${ELAP}. Check log file for run details."
# else
#     echo "HLA-HD run failed for sample ${SAMPLE_NAME}."
#     exit 1
# fi