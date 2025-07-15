#!/usr/bin/env bash

# Usage: ./fusins-preproc--nf.sh fusins-genepair.txt filtered-agfusion-dir output_file.txt
# If output_file is not provided, it will overwrite the input file

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_file> <target_directory> [output_file]"
    echo "Example: $0 fusins.txt filtered-agfusion-dir filtered_fusins.txt"
    exit 1
fi

INPUT_FILE="$1"
TARGET_DIR="$2"
OUTPUT_FILE="${3:-"${INPUT_FILE%.*}-filtered.txt"}"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found"
    exit 1
fi

# Check if target directory exists
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Target directory '$TARGET_DIR' not found"
    exit 1
fi

# Create a temporary file for the filtered results
TEMP_FILE=$(mktemp)

# Read each line from the input file
while IFS= read -r line; do
    # Skip empty lines
    if [ -z "$line" ]; then
        continue
    fi
    
    # Handle quoted lines by removing quotes, preserve everything else as a unit
    if [[ "$line" == \"*\" ]]; then
        # For quoted lines, remove only the quotes but keep everything else including commas
        clean_line=$(echo "$line" | sed 's/^"//;s/"$//')
    else
        # For non-quoted lines, keep everything as is (preserve commas too)
        clean_line="$line"
    fi
    
    # Check if any subdirectory contains this fusion pair as a substring
    found=false
    for subdir in "$TARGET_DIR"/*; do
        if [ -d "$subdir" ]; then
            subdir_name=$(basename "$subdir")
            if [[ "$subdir_name" == *"$clean_line"* ]]; then
                found=true
                break
            fi
        fi
    done
    
    # If found, add the original line to temp file
    if [ "$found" = true ]; then
        echo "$line" >> "$TEMP_FILE"
    fi
done < "$INPUT_FILE"

# Move the temp file to the output location
mv "$TEMP_FILE" "$OUTPUT_FILE"

echo "Filtering complete. Results saved to: $OUTPUT_FILE"
echo "Lines retained: $(wc -l < "$OUTPUT_FILE")"