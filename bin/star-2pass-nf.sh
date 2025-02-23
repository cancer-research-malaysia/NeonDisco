#!/usr/bin/bash

echo $(id)

READ1=$1
READ2=$2
SAMPLE_ID=$3
CORES=$4
INDEX=$5

STAR --runThreadN "${CORES}" \
        --genomeDir "${INDEX}" \
        --readFilesIn "${READ1}" "${READ2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE_ID}-STAR_2pass_" \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outFilterMultimapNmax 50 \
        --chimOutType WithinBAM HardClip \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreSeparation 1 \
        --chimScoreJunctionNonGTAG 0 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 50 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --peOverlapNbasesMin 10 \
        --twopassMode Basic
