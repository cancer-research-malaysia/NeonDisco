#!/usr/bin/env bash

echo $(id)

READ1=$1
READ2=$2
SAMPLE_ID=$3
CORES=$4

fastp -i "${READ1}" -I "${READ2}" -o "${SAMPLE_ID}_trimmed.R1.fq.gz" -O "${SAMPLE_ID}_trimmed.R2.fq.gz" --overrepresentation_analysis --detect_adapter_for_pe --thread $"${CORES}"
