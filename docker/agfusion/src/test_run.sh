#!/usr/bin/bash

# Set env variables
TSVINPUT="/data/test.tsv"
TOOL="arriba"
AGF_DB="/db"
OUTDIR="/out"

# This script is used to run AGFusion on a single sample
if agfusion batch -f ${TSVINPUT} -a ${TOOL} -db "${AGF_DB}/agfusion.homo_sapiens.95.db" -o ${OUTDIR}/agf_output_${TOOL} --middlestar --noncanonical; then
    echo "AGFusion run successfully."
else
    echo "AGFusion run failed."
fi

