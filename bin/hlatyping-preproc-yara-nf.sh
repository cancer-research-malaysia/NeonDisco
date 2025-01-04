#!/usr/bin/bash 

echo $(id)

OUTDIR="/work/nf_work"
SAMPLE_ID=$1
BAM_INPUT=$2

#################################################
# samtools --help
# picard SamToFastq --help

# first preprocess WES data
# use "chr6:28477797-33448354" if using hg19

echo "${SAMPLE_ID}"
echo "${BAM_INPUT}"


    