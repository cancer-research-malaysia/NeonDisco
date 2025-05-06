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
        --genomeLoad LoadAndRemove \
        --readFilesIn "${READ1}" "${READ2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE_ID}-STAR_1pass_" \
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
        --peOverlapNbasesMin 10
        ;;
    *.fastq | *.fq)
    # run STAR without zcat flag
        echo "Running STAR without zcat flag"
        STAR --runThreadN "${CORES}" \
        --genomeDir "${INDEX}" \
        --genomeLoad LoadAndRemove \
        --readFilesIn "${READ1}" "${READ2}" \
        --outFileNamePrefix "${SAMPLE_ID}-STAR_1pass_" \
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
        --peOverlapNbasesMin 10
        ;;
    *)
        echo "Input reads are not in fastq or fastq.gz format"
        exit 1
        ;;
esac

