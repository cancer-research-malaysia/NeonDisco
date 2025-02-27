#!/usr/bin/bash

SAMPLE_ID=$1
FQ_INPUT1=$2
FQ_INPUT2=$3
CORES=$4
BT_IDX=$5

#################################################

echo "${SAMPLE_ID}"
echo "READ 1:${FQ_INPUT1} & READ 2:${FQ_INPUT2}"

# map fastq files to HLA loci
mkdir -p "bowtie2-o"
if bowtie2 -p ${CORES} -x "${BT_IDX}" -1 "${FQ_INPUT1}" -2 "${FQ_INPUT2}" -S "${SAMPLE_ID}.hlamap.sam"; then
    echo "HLA-loci read mapping for sample ${SAMPLE_ID} is complete."
    # extract mapped reads
    samtools view -h -F 4 "${SAMPLE_ID}.hlamap.sam" > "${SAMPLE_ID}.MAPPED.sam" && echo "Mapped reads extracted." || echo "Mapped reads extraction failed." && exit 1
else
    echo "HLA-loci read mapping for sample ${SAMPLE_ID} failed."
    exit 1
fi

# convert mapped to fastq
mkdir -p "picard-o"
if picard SamToFastq -I "bowtie2-o/${SAMPLE_ID}.MAPPED.sam" -F "picard-o/${SAMPLE_ID}_HLA-UNMAP_R1.fastq" -F2 "picard-o/${SAMPLE_ID}_HLA-UNMAP_R2.fastq" -VALIDATION_STRINGENCY SILENT -QUIET TRUE; then
    echo "Picard conversion successful for sample ${SAMPLE_ID}."
    # relabel records with awk
    awk '{if(NR%4 == 1){O=$0; gsub("/1","1",O); print O}else{print $0}}' "picard-o/${SAMPLE_ID}_HLA-UNMAP_R1.fastq" > "${SAMPLE_ID}_Bam2Fq_R1.fastq" && \
    awk '{if(NR%4 == 1){O=$0; gsub("/2","2",O); print O}else{print $0}}' "picard-o/${SAMPLE_ID}_HLA-UNMAP_R2.fastq" > "${SAMPLE_ID}_Bam2Fq_R2.fastq"
else
    echo "Picard conversion failed for sample ${SAMPLE_ID}."
    exit 1
fi