#!/usr/bin/bash

SAMPLE_ID=$1
BAM_INPUT=$2

#################################################
# this strategy is used for either input data in aligned bam formats, or RNA-based input data in aligned bam format
# 	The genomic region that harbors the MHC: In GRCh37, it corresponds to chr6:28,477,797-33,448,354 (6p22.1-21.3).In GRCh38, it corresponds to chr6:28,510,020-33,480,577 (Sakaue et al., 2023)

echo "${SAMPLE_ID}"
echo "${BAM_INPUT}"

# run samtools filtering
samtools view -b -h "${BAM_INPUT}" "6:28,510,020-33,480,577" > "MHC-${SAMPLE_ID}.bam"
samtools view -b -f 4 "${BAM_INPUT}" > "unmapped-${SAMPLE_ID}.bam"
samtools merge -o "merged-${SAMPLE_ID}.bam" "MHC-${SAMPLE_ID}.bam" "unmapped-${SAMPLE_ID}.bam"

# run Picard conversion
if picard SamToFastq -I "merged-${SAMPLE_ID}.bam" -F "${SAMPLE_ID}_UNMAP_MERGED_R1.fastq" -F2 "${SAMPLE_ID}_UNMAP_MERGED_R2.fastq" -VALIDATION_STRINGENCY SILENT -QUIET TRUE; then
    echo "Picard conversion successful for sample ${SAMPLE_ID}."
    # relabel records with awk
    awk '{if(NR%4 == 1){O=$0; gsub("/1","1",O); print O}else{print $0}}' "${SAMPLE_ID}_UNMAP_MERGED_R1.fastq" > "${SAMPLE_ID}_Bam2Fq_R1.fastq" && \
    awk '{if(NR%4 == 1){O=$0; gsub("/2","2",O); print O}else{print $0}}' "${SAMPLE_ID}_UNMAP_MERGED_R2.fastq" > "${SAMPLE_ID}_Bam2Fq_R2.fastq"
else
    echo "Picard conversion failed for sample ${SAMPLE_ID}."
    exit 1
fi
    