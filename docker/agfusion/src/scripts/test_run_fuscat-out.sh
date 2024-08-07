#!/usr/bin/bash

# Set env variables
TSVINPUT="/work/data"
TOOL="fusioncatcher"
AGF_DB="/work/db"
OUTDIR="/work/out"

# This script is used to run AGFusion on a single sample
# if agfusion batch -f "${TSVINPUT}/final-list_candidate-fusion-genes.txt" -a ${TOOL} -db "${AGF_DB}/agfusion.homo_sapiens.95.db" -o "${OUTDIR}/agf_output_${TOOL}" --middlestar --noncanonical; then
#     echo "AGFusion run successfully."
# else
#     echo "AGFusion run failed."
# fi

#################################################

# This script is used to run AGFusion on multiple samples in a loop
# first search the data space for unique sample directories
export SAMPLE_DIRS=($(find "${TSVINPUT}" -mindepth 1 -maxdepth 1 -type d))
# this command expects all directories under ${TSVINPUT} to be sample directories
for SAMPLE in "${SAMPLE_DIRS[@]}"; do
    SAMPLE_NAME=$(basename "${SAMPLE}")
    echo "Processing sample: ${SAMPLE_NAME}"
    if [ -f "${SAMPLE}/final-list_candidate-fusion-genes.txt" ]; then
        echo "Fusion file found for sample ${SAMPLE_NAME}."
        # Run AGFusion on the sample
        echo "Running AGFusion for sample ${SAMPLE_NAME}..."
        # if agfusion batch -f "${SAMPLE}/final-list_candidate-fusion-genes.txt" -a ${TOOL} -db "${AGF_DB}/agfusion.homo_sapiens.95.db" -o "${OUTDIR}/agf_out_${SAMPLE_NAME}_${TOOL}" --middlestar --noncanonical; then
        #     echo "AGFusion run successfully for sample ${SAMPLE_NAME}."
        # else
        #     echo "AGFusion run failed for sample ${SAMPLE_NAME}."
        # fi
    else 
        echo "Fusion file not found for sample ${SAMPLE_NAME}. Skipping."
        continue
    fi
done