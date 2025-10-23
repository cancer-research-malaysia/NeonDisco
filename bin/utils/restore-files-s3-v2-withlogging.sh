#!/bin/bash
# Script to restore RNA-seq files based on sample list
set -e  # Exit on any error

FILELIST=$1
S3_BUCKET=""
S3_PATH=""
LOG_FILE="restore_$(date +%Y%m%d_%H%M%S).log"

# Check if file list is provided
if [[ -z "$FILELIST" ]]; then
    echo "Usage: $0 <sample_list_file>"
    exit 1
fi

# Check if file exists
if [[ ! -f "$FILELIST" ]]; then
    echo "Error: File $FILELIST not found"
    exit 1
fi

echo "Starting restoration process at $(date)" | tee "$LOG_FILE"

# Read each line from the file
while IFS= read -r sample_num; do
    # Skip empty lines
    [[ -z "$sample_num" ]] && continue
    
    # Construct filenames
    READ1="${S3_PATH}/${sample_num}T_r1.fq.gz"
    READ2="${S3_PATH}/${sample_num}T_r2.fq.gz"
    
    echo "Processing sample: $sample_num" | tee -a "$LOG_FILE"
    
    # Restore READ1
    if aws s3api restore-object --bucket "${S3_BUCKET}" --key "${READ1}" --restore-request '{"Days":180,"GlacierJobParameters":{"Tier":"Standard"}}' 2>> "$LOG_FILE"; then
        echo "✓ Initiated restore: $READ1" | tee -a "$LOG_FILE"
    else
        echo "✗ Failed to restore: $READ1" | tee -a "$LOG_FILE"
    fi
    
    # Restore READ2
    if aws s3api restore-object --bucket "${S3_BUCKET}" --key "${READ2}" --restore-request '{"Days":180,"GlacierJobParameters":{"Tier":"Standard"}}' 2>> "$LOG_FILE"; then
        echo "✓ Initiated restore: $READ2" | tee -a "$LOG_FILE"
    else
        echo "✗ Failed to restore: $READ2" | tee -a "$LOG_FILE"
    fi
    
    echo "---" | tee -a "$LOG_FILE"
    
done < "${FILELIST}"

echo "Restoration process completed at $(date)" | tee -a "$LOG_FILE"
echo "Check $LOG_FILE for detailed results"
