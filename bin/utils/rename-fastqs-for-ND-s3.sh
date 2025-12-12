#!/bin/bash

# Configuration
SOURCE_BUCKET="s3://crm.ntuh"
S3_PREFIX="sequencing-data/RNAseq/fastq_raw" 
DEST_BUCKET="s3://crmy-gb-main/neondisco/dataset-ntuh/RNAseq-fastq-raw"
 
METADATA_FILE="NTUH_sample_renamed_metadata.txt"
DRY_RUN=true  # Set to false to actually execute

# Limit results in dry-run mode for faster testing
DRY_RUN_LIMIT=()
if [ "$DRY_RUN" = true ]; then
  DRY_RUN_LIMIT=(awk 'NR<=10')
else
  DRY_RUN_LIMIT=(cat)
fi

# Initialize metadata file
if [ "$DRY_RUN" = false ]; then
  echo -e "sample_id\toriginal_R1\toriginal_R2\tnew_R1\tnew_R2" > "$METADATA_FILE"
fi

# Counter for naming
counter=1

echo "========================================="
echo "S3 Paired-End File Rename Script"
echo "========================================="
echo "SOURCE BUCKET: ${SOURCE_BUCKET}/${S3_PREFIX}"
echo "DEST BUCKET:   ${DEST_BUCKET}/"
echo "DRY RUN MODE:  ${DRY_RUN}"
echo "========================================="
echo ""

# Get all fastq.gz files and extract unique sample prefixes
# This assumes files are named like: sample_R1_001.fastq.gz and sample_R2_001.fastq.gz
declare -A processed_samples

while read -r file_path; do

    filename=$(basename "$file_path")
    
    # Extract sample prefix by removing _R1_001 or _R2_001 suffix
    # This handles patterns like: A010001RSS00_S126_R1_001.fastq.gz
    sample_prefix=$(echo "$filename" | sed -E 's/_R[12]_001\.fastq\.gz$//')
    
    # Skip if we've already processed this sample
    if [[ -v "processed_samples[$sample_prefix]" ]]; then
      continue
    fi
    processed_samples[$sample_prefix]=1
    
    # Find both R1 and R2 files for this sample
    r1_file=""
    r2_file=""
    r1_path=""
    r2_path=""
    
    # Search for R1 and R2 variants (load list into array for faster iteration)
    mapfile -t search_paths < <(aws s3 ls "${SOURCE_BUCKET}/${S3_PREFIX}" --recursive | grep '\.fastq\.gz$' | awk '{print $4}')
    for search_path in "${search_paths[@]}"; do
      search_file=$(basename "$search_path")
      if [[ "$search_file" == "${sample_prefix}_R1_001.fastq.gz" ]]; then
        r1_file="$search_file"
        r1_path="$search_path"
      elif [[ "$search_file" == "${sample_prefix}_R2_001.fastq.gz" ]]; then
        r2_file="$search_file"
        r2_path="$search_path"
      fi
    done
    unset search_paths
    
    # Check if we found both pairs
    if [[ -z "$r1_file" || -z "$r2_file" ]]; then
      echo "⚠ WARNING: Incomplete pair for sample prefix: ${sample_prefix}"
      echo "  Found R1: ${r1_file:-MISSING}"
      echo "  Found R2: ${r2_file:-MISSING}"
      echo "  Skipping this sample..."
      echo ""
      continue
    fi
    
    # Generate new names
    new_r1="${counter}T_R1.fastq.gz"
    new_r2="${counter}T_R2.fastq.gz"
    
    # Construct full S3 paths
    original_r1_s3="${SOURCE_BUCKET}/${r1_path}"
    original_r2_s3="${SOURCE_BUCKET}/${r2_path}"
    new_r1_s3="${DEST_BUCKET}/${new_r1}"
    new_r2_s3="${DEST_BUCKET}/${new_r2}"
    
    # Print detailed information
    echo "========================================"
    echo "PAIRED SAMPLE #${counter}: ${sample_prefix}"
    echo "========================================"
    echo "R1 (Forward Read):"
    echo "  Original path:      ${r1_path}"
    echo "  Original filename:  ${r1_file}"
    echo "  New filename:       ${new_r1}"
    echo "  Source S3 URI:      ${original_r1_s3}"
    echo "  Destination S3 URI: ${new_r1_s3}"
    echo ""
    echo "R2 (Reverse Read):"
    echo "  Original path:      ${r2_path}"
    echo "  Original filename:  ${r2_file}"
    echo "  New filename:       ${new_r2}"
    echo "  Source S3 URI:      ${original_r2_s3}"
    echo "  Destination S3 URI: ${new_r2_s3}"
    echo ""
    
    if [ "$DRY_RUN" = true ]; then
      echo "  [DRY RUN] Would execute:"
      echo "    aws s3 cp \"${original_r1_s3}\" \"${new_r1_s3}\""
      echo "    aws s3 cp \"${original_r2_s3}\" \"${new_r2_s3}\""
    else
      echo "  Executing copies..."
      
      # Copy R1
      aws s3 cp "$original_r1_s3" "$new_r1_s3"
      r1_status=$?
      
      # Copy R2
      aws s3 cp "$original_r2_s3" "$new_r2_s3"
      r2_status=$?
      
      if [ $r1_status -eq 0 ] && [ $r2_status -eq 0 ]; then
        echo "  ✓ Both files copied successfully"
        # Record mapping
        echo -e "${counter}T\t${r1_file}\t${r2_file}\t${new_r1}\t${new_r2}" >> "$METADATA_FILE"
      else
        echo "  ✗ Copy failed!"
        [ $r1_status -ne 0 ] && echo "    R1 copy failed"
        [ $r2_status -ne 0 ] && echo "    R2 copy failed"
      fi
    fi
    
    echo ""
    ((counter++))
  
  done < <(aws s3 ls "${SOURCE_BUCKET}/${S3_PREFIX}" --recursive | grep '\.fastq\.gz$' | "${DRY_RUN_LIMIT[@]}" | awk '{print $4}')

echo "========================================="
if [ "$DRY_RUN" = true ]; then
  echo "DRY RUN COMPLETE - No changes made"
  echo "Set DRY_RUN=false to execute"
else
  echo "Renaming complete!"
  echo "Metadata saved to: ${METADATA_FILE}"
fi
echo "Total sample pairs processed: $((counter - 1))"
echo "========================================="
