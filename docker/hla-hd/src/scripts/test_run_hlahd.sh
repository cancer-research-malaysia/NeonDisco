#!/usr/bin/bash

# Set env variables
TSVINPUT="/work/data"
OUTDIR="/work/out"

# # first preprocess WES data
# mkdir "${OUTDIR}/samtools-out"
# # use "chr6:28477797-33448354" for hg19
# samtools view -b -h in.bam "chr6:28510120-33480577" > "${OUTDIR}/samtools-out/MHC.bam" 
# samtools view -b -f 4 in.bam > unmapped.bam
# samtools merge -o merged.bam MHC.bam unmapped.bam
# SamToFastq I=merged.bam F=merged_R1.fq F2=merged_R2.fq
# awk '{if(NR%4 == 1){O=$0; gsub("/1","1",O); print O}else{print $0}}' merged_R1.fq > ${sample-id}-MHC_WES-b2f-R1.fq
# awk '{if(NR%4 == 1){O=$0; gsub("/2","2",O); print O}else{print $0}}' merged_R2.fq > ${sample-id}-MHC_WES-b2f-R2.fq
# # execution of HLA-HD (4 threads, discard reads <30bp)
# hlahd.sh -t 4 -m 30 -f ${hlahd-dir}/freq_data/ ${sample-id}-MHC_WES-b2f-R1.fq ${sample-id}-MHC_WES-b2f-R2.fq ${hlahd-
# dir}/HLA_gene.split.3.50.0.txt ${hlahd-dir}/dictionary/ ${sample-id}_WES-MHC-bam ${outdir} >> log.txt 2>&1

#################################################


# This script is used to run HLA-HD for multiple samples in a loop
# first search the data space for unique bam files
export BAMS=($(find "${TSVINPUT}" -mindepth 1 -maxdepth 1 -type f -name '*.bam'))

for SAMPLE in "${BAMS[@]}"; do
    # Extract the sample name from the file path
    SAMPLE_NAME=$(basename "${SAMPLE}" .bam)
    echo "Processing sample: ${SAMPLE_NAME}"
    samtools --help
    SamToFastq --help
#     if [ -f "${SAMPLE}/arriba-fusions.tsv" ]; then
#         echo "Fusion file found for sample ${SAMPLE_NAME}."
#         # Run AGFusion on the sample
#         echo "Running AGFusion for sample ${SAMPLE_NAME}..."
#         if agfusion batch -f "${SAMPLE}/arriba-fusions.tsv" -a ${TOOL} -db "${AGF_DB}/agfusion.homo_sapiens.95.db" -o "${OUTDIR}/agf_output_${SAMPLE_NAME}_${TOOL}" --middlestar --noncanonical; then
#             echo "AGFusion run successfully for sample ${SAMPLE_NAME}."
#         else
#             echo "AGFusion run failed for sample ${SAMPLE_NAME}."
#         fi
#     else 
#         echo "Fusion file not found for sample ${SAMPLE_NAME}. Skipping."
#         continue
#     fi
done


