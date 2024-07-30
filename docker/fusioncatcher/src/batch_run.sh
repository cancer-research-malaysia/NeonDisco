#!/usr/bin/bash

run_fuscat() {
  local DB=$1
  local FASTQ_R1_FILE=$2
  local FASTQ_R2_FILE=$3
  local OUTDIR=$4

  fusioncatcher -d ${DB} --input ${FASTQ_R1_FILE},${FASTQ_R2_FILE} --output ${OUTDIR}
}

# Set env variables
export DB="/db"
export INP_DIR="/in"
export OUTDIR_PREFIX="/out"
echo "Environment variables set and exported! Finding fastq files..."

# Check if the subdirectory contains a .fastq.gz or fq.gz file
FILES=$(find "$INP_DIR" \( -name '*.fastq.gz' -o -name '*.fq.gz' \) -type f > found.txt && cat found.txt | wc -l)

# check if FILES is zero or not
if [[ ${FILES} -eq 0 ]]; then
    echo "Fastq files are not present. Skipping..."
    exit 1
else
    echo "Fastq files found!"
    echo "Total number of files: ${FILES}"
    # extract read1 group number of files
    R1_COUNT=$(find "$INP_DIR" -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | xargs -n 1 basename | awk -F'_' '/r1/ {print $1}' | sort | uniq | wc -l) && echo $R1_COUNT
    echo "Total number of R1 files: $R1_COUNT"
    # extract read2 group number of files
    R2_COUNT=$(find "$INP_DIR" -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | xargs -n 1 basename | awk -F'_' '/r2/ {print $1}' | sort | uniq | wc -l) && echo $R2_COUNT
    echo "Total number of R2 files: $R2_COUNT"
    # extract sample IDs into an array
    readarray -t SAMPLE_ID < <(find "$INP_DIR" -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | xargs -n 1 basename | awk -F'_' '/r2/ {print $1}' | sort | uniq)
    # print out the sample IDs
    echo "Sample IDs: ${SAMPLE_ID[@]}"

    # Check if R1_COUNT == R2_COUNT (implying paired reads data)
    if (( R1_COUNT == R2_COUNT )); then
        if (( R1_COUNT == 0 && R2_COUNT == 0 )); then 
            echo "The sample counts are zeroes. Something is wrong. Check your input files."
            exit 1
        else
            echo "The input files appear to be paired. Looping through the sample ID array..."
            for prefix in "${SAMPLE_ID[@]}"; do
                export prefix
                echo "Sample ID: ${prefix}"
                # find the paired fastq files for the sample ID and assign to variables
                FASTQS=$(find "$INP_DIR" -mindepth 1 -type f \( -name "*${prefix}*fastq.gz" -o -name "*${prefix}*.fq.gz" \))
                echo "Paired fastq files: ${FASTQS}"
                export FASTQS
                # Extract the first two lines
                FILE_R1=$(echo "$FASTQS" | sed -n '1p')
                FILE_R2=$(echo "$FASTQS" | sed -n '2p')

                echo "First file: $FILE_R1"
                echo "Second file: $FILE_R2"

                if [[ -n "$FILE_R1" && -n "$FILE_R2" ]]; then
                    # make output directory
                    mkdir -p "${OUTDIR_PREFIX}/${prefix}"
                    # measure execution time
                    STARTTIME=$(date +%s)
                    # if run_fuscat "${DB}" "${FILE_R1}" "${FILE_R2}" "${OUTDIR_PREFIX}/${prefix}"; then
                    #     ENDTIME=$(date +%s)
                    #     ELAP=$(( ENDTIME - STARTTIME ))
                    #     echo "FusionCatcher run completed successfully. Time taken: ${ELAP}. Check log file for run details."
                    # else
                    #     echo "Something went wrong during FusionCatcher run. Check log file."
                    # fi
                else
                    echo "Paired fastq files not found for ${prefix}. Skipping..."
                    continue
                fi
            
            done
        fi
    else
        echo "The number of files is odd. Make sure the input files are paired before proceeding."
        exit 1
    fi
fi
