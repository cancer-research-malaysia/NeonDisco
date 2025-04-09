#!/bin/bash

# This utility script generates a manifest file from a list of sample names

# check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <SOURCEPATH> <SAMPLE_LIST> <OUTPUT_MANIFEST>"
    echo "    SOURCEPATH: Source path (either POSIX paths for local directory or S3 bucket remote path) where the input files are stored"
    echo "    SAMPLE_LIST: File containing sample names (this is ideally the prefix of the input read files; e.g., for a file named sample1_r1.fq.gz, the prefix would be sample1)"
    echo "    OUTPUT_MANIFEST: Output manifest filename"

    echo "Example: $0 s3://bucket/path/to/files/ sample_list.txt manifest.tsv"
    echo "Please provide all the required arguments at runtime."
    exit 1
fi

# Set up variables
SOURCEPATH=$1  # S3 bucket path where the input files are stored
SAMPLE_LIST=$2  # File containing sample names (this is ideally the prefix of the input read files; e.g., for a file named sample1_r1.fq.gz, the prefix would be sample1)
OUTPUT_MANIFEST=$3 # Output manifest filename

# Validate input file exists
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "Error: $SAMPLE_LIST file not found. Please create a file with sample names."
    exit 1
fi

# Create manifest file with header
echo -e "sampleName\tread1Path\tread2Path" > "$OUTPUT_MANIFEST"

# Process each sample
while read -r sample_base; do
    # Skip empty lines or comments
    [[ -z "$sample_base" || "$sample_base" =~ ^# ]] && continue
    
    # Define the R1 and R2 file paths on S3
    r1_path="${SOURCEPATH}/${sample_base}_r1.fq.gz"
    r2_path="${SOURCEPATH}/${sample_base}_r2.fq.gz"
    
    # Write to manifest file
    echo -e "${sample_base}\t${r1_path}\t${r2_path}" >> "$OUTPUT_MANIFEST"
    
    echo "Added sample: $sample_base"
done < "$SAMPLE_LIST"

# Show manifest file contents
echo "Manifest file created: $OUTPUT_MANIFEST"
cat "$OUTPUT_MANIFEST"

# Print total number of samples processed
total_samples=$(grep -v "^sampleName" "$OUTPUT_MANIFEST" | wc -l)
echo "Total samples processed: $total_samples"
