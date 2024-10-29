#!/usr/bin/bash 

echo $(id)  

# ENSURE THAT THE fq.gz HAVE BEEN PREPROCESSED FIRST BEFORE YOU RUN HLA-HD! 
# get script arguments

READ_ONE=$1
READ_TWO=$2
HLAHD_DIR="/tmp/hlahd.1.7.0"

# check R1 and R2 fq file existence 
if [ -f "${READ_ONE}" ] && [ -f "${READ_TWO}" ]; then
    echo "Running HLA-HD for sample ${SAMPLE_NAME}."
    mkdir -p "${OUTDIR}/logs/"
    mkdir -p "${OUTDIR}/hla-hd-out"
    # measure execution time
    STARTTIME=$(date +%s)
    # run HLA-HD (8 threads, discard reads <30bp)
    if hlahd.sh -t 8 -m 30 -f ${HLAHD_DIR}/freq_data/ "${OUTDIR}/merged-${SAMPLE_NAME}-MHC-WES-b2f-R1.fq" "${OUTDIR}/merged-${SAMPLE_NAME}-MHC-WES-b2f-R2.fq" "${HLAHD_DIR}/HLA_gene.split.3.50.0.txt" "${HLAHD_DIR}/dictionary/" "${SAMPLE_NAME}_WES-MHC-bam" "${OUTDIR}/hla-hd-out" 2>&1 | tee "${OUTDIR}/logs/${SAMPLE_NAME}-HLAHD-run.log-$(date +%Y%m%d_%H-%M-%S).txt"; then
        echo "HLA-HD run successfully for sample ${SAMPLE_NAME}."
        ENDTIME=$(date +%s)
        ELAP=$(( ENDTIME - STARTTIME ))
        echo "Time taken: ${ELAP}. Check log file for run details."
    else
        echo "HLA-HD run failed for sample ${SAMPLE_NAME}."
        exit 1
    fi
else
    echo "One or both FASTQ files are missing for sample ${SAMPLE_NAME}."
    exit 1
fi
