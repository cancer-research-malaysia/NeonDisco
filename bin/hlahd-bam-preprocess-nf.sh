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

# run samtools filtering
mkdir -p "${OUTDIR}/samtools-o"
samtools view -b -h "${BAM_INPUT}" "6:28510120-33480577" > "${OUTDIR}/samtools-o/MHC-${SAMPLE_ID}.bam"
samtools view -b -f 4 "${BAM_INPUT}" > "${OUTDIR}/samtools-o/unmapped-${SAMPLE_ID}.bam"
samtools merge -o "${OUTDIR}/samtools-o/merged-${SAMPLE_ID}.bam" "${OUTDIR}/samtools-o/MHC-${SAMPLE_ID}.bam" "${OUTDIR}/samtools-o/unmapped-${SAMPLE_ID}.bam"

# run Picard conversion
mkdir -p "${OUTDIR}/picard-o"
if picard SamToFastq -I "${OUTDIR}/samtools-o/merged-${SAMPLE_ID}.bam" -F "${OUTDIR}/picard-o/${SAMPLE_ID}_UNMAP_MERGED_R1.fastq" -F2 "${OUTDIR}/picard-o/${SAMPLE_ID}_UNMAP_MERGED_R2.fastq" -VALIDATION_STRINGENCY SILENT -QUIET TRUE; then
    echo "Picard conversion successful for sample ${SAMPLE_ID}."
    # relabel records with awk
    awk '{if(NR%4 == 1){O=$0; gsub("/1","1",O); print O}else{print $0}}' "${OUTDIR}/picard-o/${SAMPLE_ID}_UNMAP_MERGED_R1.fastq" > "${SAMPLE_ID}_Bam2Fq_R1.fastq" && \
    awk '{if(NR%4 == 1){O=$0; gsub("/2","2",O); print O}else{print $0}}' "${OUTDIR}/picard-o/${SAMPLE_ID}_UNMAP_MERGED_R2.fastq" > "${SAMPLE_ID}_Bam2Fq_R2.fastq"
fi
    