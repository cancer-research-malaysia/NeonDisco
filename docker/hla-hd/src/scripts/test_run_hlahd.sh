#!/usr/bin/bash

# Set env variables
TSVINPUT="/work/data"
OUTDIR="/work/out"
SAMT-OUT="${OUTDIR}/samtools-out"
mkdir ${SAMT-OUT}

#################################################
samtools --help
picard SamToFastq --help
hlahd.sh --help

# This script is used to run HLA-HD for multiple samples in a loop
# first search the data space for unique bam files
export BAMS=($(find "${TSVINPUT}" -mindepth 1 -maxdepth 1 -type f -name '*.bam'))
echo "Found ${#BAMS[@]} BAM files in the input directory."

# Check if the bam array is empty or not
if [ ${#BAMS[@]} -eq 0 ]; then
    echo "No BAM files found in the input directory. Exiting."
    exit 1
else
    # Loop through each BAM file
    for SAMPLE in "${BAMS[@]}"; do
        # Extract the sample name from the file path
        SAMPLE_NAME=$(basename "${SAMPLE}" .bam)
        echo "Processing sample: ${SAMPLE_NAME}"
        # # first preprocess WES data
        # # use "chr6:28477797-33448354" for hg19
        # samtools view -b -h "$SAMPLE" "chr6:28510120-33480577" > "${SAMT-OUT}/MHC-${SAMPLE_NAME}.bam" 
        # samtools view -b -f 4 "$SAMPLE" > "${SAMT-OUT}/unmapped-${SAMPLE_NAME}.bam"
        # samtools merge -o "${SAMT-OUT}/merged-${SAMPLE_NAME}.bam" "${SAMT-OUT}/MHC-${SAMPLE_NAME}.bam" "${SAMT-OUT}/unmapped-${SAMPLE_NAME}.bam"

        # # run Picard conversion
        # if picard SamToFastq I="${SAMT-OUT}/merged-${SAMPLE_NAME}.bam" F="${SAMT-OUT}/merged-${SAMPLE_NAME}_R1.fq" F2="${SAMT-OUT}/merged-${SAMPLE_NAME}_R2.fq"; then
        #     echo "Picard conversion successful for sample ${SAMPLE_NAME}."
        #     # relabel records with awk
        #     awk '{if(NR%4 == 1){O=$0; gsub("/1","1",O); print O}else{print $0}}' "${SAMT-OUT}/merged-${SAMPLE_NAME}_R1.fq" > "${OUTDIR}/merged-${SAMPLE_NAME}-MHC-WES-b2f-R1.fq" && \
        #     awk '{if(NR%4 == 1){O=$0; gsub("/2","2",O); print O}else{print $0}}' "${SAMT-OUT}/merged-${SAMPLE_NAME}_R2.fq" > "${OUTDIR}/merged-${SAMPLE_NAME}-MHC-WES-b2f-R2.fq"
            
        #     HLAHD_DIR="/tmp/hlahd.1.7.0"
        #     # check R1 and R2 fq file existence 
        #     if [ -f "${OUTDIR}/merged-${SAMPLE_NAME}-MHC-WES-b2f-R1.fq" ] && [ -f "${OUTDIR}/merged-${SAMPLE_NAME}-MHC-WES-b2f-R2.fq" ]; then
        #         echo "Running HLA-HD for sample ${SAMPLE_NAME}."
                
        #         mkdir -p "${OUTDIR}/logs/"

        #         # measure execution time
        #         STARTTIME=$(date +%s)

        #         # run HLA-HD (8 threads, discard reads <30bp)
        #         if hlahd.sh -t 8 -m 30 -f ${HLAHD_DIR}/freq_data/ "${OUTDIR}/merged-${SAMPLE_NAME}-MHC-WES-b2f-R1.fq" "${OUTDIR}/merged-${SAMPLE_NAME}-MHC-WES-b2f-R2.fq" "${HLAHD_DIR}/HLA_gene.split.3.50.0.txt" "${HLAHD_DIR}/dictionary/" "${OUTDIR}/${SAMPLE_NAME}_WES-MHC-bam" ${OUTDIR} 2>&1 | tee "${OUTDIR}/logs/${SAMPLE_NAME}-HLAHD-run.log-$(date +%Y%m%d_%H-%M-%S).txt"; then
        #             echo "HLA-HD run successfully for sample ${SAMPLE_NAME}."
        #             ENDTIME=$(date +%s)
        #             ELAP=$(( ENDTIME - STARTTIME ))
        #             echo "Time taken: ${ELAP}. Check log file for run details."
        #         else
        #             echo "HLA-HD run failed for sample ${SAMPLE_NAME}."
        #             exit 1
        #         fi
        #     else
        #         echo "One or both FASTQ files are missing for sample ${SAMPLE_NAME}."
        #         exit 1
        #     fi
        # else
        #     echo "Picard conversion failed for sample ${SAMPLE_NAME}."
        #     exit 1
        # fi
    done
fi
