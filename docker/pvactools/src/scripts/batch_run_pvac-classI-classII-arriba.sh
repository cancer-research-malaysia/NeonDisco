#!/usr/bin/bash

INDIR="/work/input"
OUTDIR="/work/output"
TOOL="arriba"
ALLELES='A*11:01,HLA-A*24:02,HLA-A*33:03,HLA-A*02:07,HLA-A*11:353,HLA-A*02:01,HLA-A*02:03,HLA-A*11:02,HLA-A*11:12,HLA-B*40:01,HLA-B*46:01,HLA-B*58:01,HLA-B*13:01,HLA-B*15:02,HLA-B*38:02,HLA-B*51:01,HLA-B*54:01,HLA-B*27:04,HLA-B*15:01,HLA-C*07:02,HLA-C*01:02,HLA-C*03:04,HLA-C*08:01,HLA-C*03:02,HLA-C*15:02,HLA-C*12:02,HLA-C*04:01,HLA-C*14:02,HLA-C*06:02,HLA-C*04:03,DRB1*09:01,DRB1*12:02,DRB1*15:01,DRB1*03:01,DRB1*11:01,DRB1*08:03,DRB1*04:05,DRB1*07:01,DRB1*16:02,DRB1*15:02,DQA1*01:02,DQA1*06:01,DQA1*03:02,DQA1*01:04,DQA1*01:03,DQA1*05:01,DQA1*05:05,DQA1*03:01,DQA1*03:03,DQA1*02:01,DQB1*03:01,DQB1*03:03,DQB1*05:02,DQB1*06:01,DQB1*02:01,DQB1*05:03,DQB1*03:02,DQB1*04:01,DQB1*05:01,DQB1*02:02,DPA1*02:02,DPA1*01:03,DPA1*02:01,DPA1*04:01,DPB1*05:01,DPB1*02:01,DPB1*04:01,DPB1*02:02,DPB1*13:01,DPB1*03:01,DRA*01:01,DRA*01:02,DRB2*01:01,DRB3*02:02,DRB3*03:01,DRB3*01:01,DRB3*03:06,DRB3*03:10,DRB4*01:03,DRB5*01:01,DRB5*01:11,DRB5*01:22,DRB6*02:01,DRB6*02:02,DRB6*01:01,DRB7*01:01,DRB8*01:01,DRB9*01:01,DRB9*01:02'
#################################################

# first search the agfusion dataspace for Arriba output folders
export AGF_DIRS=($(find "${INDIR}" -mindepth 1 -maxdepth 1 -type d -name 'agf_output*'))
echo "Found ${#AGF_DIRS[@]} AGFusion folders in the input directory."

# Loop through each BAM file
for DIR in "${AGF_DIRS[@]}"; do
    # Extract the sample dirname from the file path
    SAMPLE_DIRNAME=$(basename "${DIR}")
    SAMPLE_NAME="${SAMPLE_DIRNAME##agf_output_}"
    SAMPLE_NAME="${SAMPLE_NAME%%_arriba}"
    echo "Processing sample: ${SAMPLE_NAME}"
    mkdir -p ${OUTDIR}/${SAMPLE_NAME}
    
    # measure execution time
    STARTTIME=$(date +%s)

    # run pvacfuse
    # if pvacfuse run -e1 8,9,10,11,12,13,14 -e2 11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 --iedb-install-directory /tmp/ -k -t 8 ${INDIR} pvf_${SAMPLE_NAME}_${TOOL} ${ALLELES} NetMHCpan NetMHCIIpan ${OUTDIR}/${SAMPLE_NAME} 2>&1 | tee "${OUTDIR}/${SAMPLE_NAME}/pvf-${SAMPLE_NAME}-arriba-run-$(date +%Y%m%d_%H-%M-%S).txt"; then
    #     echo "pVacFuse finished successfully for sample ${SAMPLE_NAME}."
    #     ENDTIME=$(date +%s)
    #     ELAP=$(( ENDTIME - STARTTIME ))
    #     echo "Time taken: ${ELAP}. Check log file for run details."
    # else
    #     echo "pVacFuse run failed for sample ${SAMPLE_NAME}."
    #     exit 1
    # fi
done

