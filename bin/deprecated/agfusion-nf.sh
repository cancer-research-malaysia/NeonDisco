#!/usr/bin/bash

# Set env variables
INPUT_FILE=$1
SAMPLE_ID=$2
AGF_DB="/work/libs"

if agfusion annotate -g5 MAPK13 -g3 C1QL1 -j5 36132629 -j3 44965446 -db "${AGF_DB}/agfusion.homo_sapiens.95.db" -o "${SAMPLE_ID}-MAPK13-C1QL1" --middlestar --noncanonical; then
    echo "AGFusion run successfully."
else
    echo "AGFusion run failed."
    exit 1
fi

