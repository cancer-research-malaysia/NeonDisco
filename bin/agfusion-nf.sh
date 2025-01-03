#!/usr/bin/bash

# Set env variables
INPUT_FILE=$1
SAMPLE_ID=$2
AGF_DB="/work/libs"
OUTDIR="/work/nf_work"

# Extract tool identifier from filename
TOOL_ID=${INPUT_FILE##*_}  # Get part after last _
TOOL_ID=${TOOL_ID%%.tsv}   # Remove .tsv extension

# Set TOOL based on identifier
if [ "$TOOL_ID" = "arr" ]; then
    TOOL="arriba"
elif [ "$TOOL_ID" = "fc" ]; then
    TOOL="fusioncatcher"
else
    echo "Unknown tool identifier in filename"
    exit 1
fi

if agfusion batch -f "${INPUT_FILE}" -a ${TOOL} -db "${AGF_DB}/agfusion.homo_sapiens.95.db" -o "${OUTDIR}/${SAMPLE_ID}_agfusion_${TOOL}" --middlestar --noncanonical; then
    echo "AGFusion run successfully."
else
    echo "AGFusion run failed."
    exit 1
fi

