#!/usr/bin/bash

echo $(id)

run_fusionInspector() {
  local CTAT_LIB=$1
  local FASTQS=$2
  local ID=$3
  local THREADS=$4
  local OUTDIR=$5

  FusionInspector --annotate --examine_coding_effect --predict_cosmic_like --include_Trinity --vis \
  --CPU ${THREADS} --cleanup --only_fusion_reads --extract_fusion_reads_file sample_${ID} \
  --fusions "${FUSIONLIST_DIR}/fusionlist_${ID}.txt" \
  --genome_lib_dir "${CTAT_LIB}" \
  --left_fq "${FASTQS}/${ID}_r1.fq.gz" --right_fq "${FASTQS}/${ID}_r2.fq.gz" \
  --out_prefix sample_${ID} -O "${OUTDIR}/" 2>&1 | tee "${OUTDIR}/sample_${ID}.arriba.log.$(date +%Y%m%d_%H-%M-%S).txt"
}

# Set env variables
TOOL="arriba"
CTAT_LIB="/work/libs"
FASTQS="/work/data"
FUSIONLIST_DIR="/work/fusionlist"
THREADS="8"
OUTDIR="/work/out"

echo "Environment variables set! Listing fastq files..."

# Check if the subdirectory contains a .fastq.gz or fq.gz file
if [[ $(find "$FASTQS" \( -name '*.fastq.gz' -o -name '*.fq.gz' \) -type f) ]]; then
    echo "Fastq files are present. Listing..."
    find "$FASTQS" -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \)
    echo "Total number of files: $(find "$FASTQS" -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | wc -l)"
    # extract read1 group number of files
    R1_COUNT=$(find "$FASTQS" -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | xargs -n 1 basename | awk -F'_' '/r1/ {print $1}' | sort | uniq | wc -l) && echo $R1_COUNT
    echo "Total number of R1 files: $R1_COUNT"
    # extract read2 group number of files
    R2_COUNT=$(find "$FASTQS" -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | xargs -n 1 basename | awk -F'_' '/r2/ {print $1}' | sort | uniq | wc -l) && echo $R2_COUNT
    echo "Total number of R2 files: $R2_COUNT"
    # extract sample IDs into an array
    readarray -t SAMPLE_ID < <(find "$FASTQS" -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | xargs -n 1 basename | awk -F'_' '/r2/ {print $1}' | sort | uniq)
    # print out the sample IDs
    echo "Sample IDs: ${SAMPLE_ID[@]}"

    # Check if R1_COUNT == R2_COUNT (implying paired reads data)
    if (( R1_COUNT == R2_COUNT )); then
        if (( R1_COUNT == 0 && R2_COUNT == 0 )); then 
            echo "The counts are zeroes. Something is wrong. Check your input files."
            exit 1
        else
            echo "The input files appear to be paired. Looping through the sample ID array..."
            for prefix in "${SAMPLE_ID[@]}"; do
                echo "Sample ID: ${prefix}"
                echo "Running FusionInspector..."
                mkdir -p "${OUTDIR}/from-${TOOL}/${prefix}"
                # measure execution time
                STARTTIME=$(date +%s)
                if run_fusionInspector "${CTAT_LIB}" "${FASTQS}" "${prefix}" "${THREADS}" "${OUTDIR}/from-${TOOL}/${prefix}"; then
                    ENDTIME=$(date +%s)
                    ELAP=$(( ENDTIME - STARTTIME ))
                    echo "FusionInspector run on ${TOOL} outputs completed successfully. Time taken: ${ELAP}. Check log file for run details."
                else
                    echo "Something went wrong. Check log file."
                fi
            done
        fi
    else
        echo "The number of files is odd. Make sure the input files are paired before proceeding."
        exit 1
    fi
else
    echo "Fastq input files are not found. Double check your input path."
    exit 1
fi
