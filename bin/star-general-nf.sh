#!/usr/bin/bash

READ1=$1
READ2=$2
SAMPLE_ID=$3
CORES=$4
INDEX=$5

# check if the input reads are gzipped or not, important for STAR command
case $READ1 in
    *.gz) 
    # run STAR with zcat flag
        echo "Running STAR with zcat flag"
        STAR --runThreadN "${CORES}" \
        --genomeDir "${INDEX}" \
        --readFilesIn "${READ1}" "${READ2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE_ID}-STAR-ARR_" \
        --outFilterMultimapNmax 1000 \
        --outSAMmultNmax 1 \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --chimOutJunctionFormat 1 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 200000 \
        --alignIntronMax 200000 \
        --chimSegmentReadGapMax 3 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --seedSearchStartLmax 20 \
        --winAnchorMultimapNmax 50 \
        --outSAMtype BAM Unsorted \
        --chimOutType Junctions WithinBAM \
        --outSAMunmapped Within KeepPairs
        ;;
    *.fastq | *.fq)
    # run STAR without zcat flag
        echo "Running STAR without zcat flag"
        STAR --runThreadN "${CORES}" \
        --genomeDir "${INDEX}" \
        --readFilesIn "${READ1}" "${READ2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE_ID}-STAR-ARR_" \
        --outFilterMultimapNmax 1000 \
        --outSAMmultNmax 1 \
        --chimSegmentMin 10 \
        --chimJunctionOverhangMin 10 \
        --chimOutJunctionFormat 1 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 200000 \
        --alignIntronMax 200000 \
        --chimSegmentReadGapMax 3 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --seedSearchStartLmax 20 \
        --winAnchorMultimapNmax 50 \
        --outSAMtype BAM Unsorted \
        --chimOutType Junctions WithinBAM \
        --outSAMunmapped Within KeepPairs
        ;;
    *)
        echo "Input reads are not in fastq or fastq.gz format"
        exit 1
        ;;
esac



