#!/usr/bin/bash
# get script arguments

SAMPLENAME=$1
FI_VALIDATED_TSV=$2
RECURRENT_FUSIONS_TSV=$3

echo "Filtering validated fusions for recurrent ones only..."
echo "Sample: ${SAMPLENAME}"
echo "Validated fusions: ${FI_VALIDATED_TSV}"
echo "Recurrent fusions TSV: ${RECURRENT_FUSIONS_TSV}"

# Initialize report
REPORT="${SAMPLENAME}_recurrent_filter_report.txt"
echo "Recurrent Fusion Filtering Report for ${SAMPLENAME}" > ${REPORT}
echo "===================================================" >> ${REPORT}
echo "Date: $(date)" >> ${REPORT}
echo "" >> ${REPORT}

# output file for recurrent validated fusions
OUTPUT_FILE="${SAMPLENAME}-validated-recurrent-fusions-only.tsv"
echo "Output file for recurrent validated fusions: ${OUTPUT_FILE}" >> ${REPORT}


# Check input files
if [ $(tail -n +2 ${FI_VALIDATED_TSV} | wc -l) -eq 0 ]; then
    echo "No validated fusions for this sample."
    echo "STATUS: No validated fusions!" >> ${REPORT}
    head -1 ${FI_VALIDATED_TSV} > "${OUTPUT_FILE}"
    mkdir -p validated-recurrent-agfusion-outdir
    exit 0
fi

if [ $(tail -n +2 "${RECURRENT_FUSIONS_TSV}" | wc -l) -eq 0 ]; then
    echo "No recurrent fusions found in this cohort."
    echo "STATUS: No recurrent fusions in this cohort!" >> ${REPORT}
    head -1 ${FI_VALIDATED_TSV} > "${OUTPUT_FILE}"
    mkdir -p validated-recurrent-agfusion-outdir
    exit 0
fi

# Extract recurrent gene pairs for this sample
GENEPAIRS="${SAMPLENAME}_recurrent_gene_pairs.txt"
tail -n +2 "${RECURRENT_FUSIONS_TSV}" | awk -F'\t' -v sample="${SAMPLENAME}" '$8 == sample {print $2}' | sort -u > ${GENEPAIRS}

sample_recurrent_count=$(wc -l < "${GENEPAIRS}")
echo "Recurrent fusions for this sample ${SAMPLENAME}: $sample_recurrent_count" >> ${REPORT}

# Get validated fusion count
validated_count=$(tail -n +2 "${FI_VALIDATED_TSV}" | wc -l)
echo "FusionInspector-validated fusions for this sample ${SAMPLENAME}: $validated_count" >> ${REPORT}

# Filter validated fusions for recurrent ones
head -1 "${FI_VALIDATED_TSV}" > "${OUTPUT_FILE}"

if [ $sample_recurrent_count -gt 0 ]; then
    while IFS= read -r gene_pair_row; do
        if [ -n "$gene_pair_row" ]; then  # Check if gene pair is not empty
            tail -n +2 "${FI_VALIDATED_TSV}" | awk -F'\t' -v gene_pair="$gene_pair_row" '$2 == gene_pair' >> "${OUTPUT_FILE}"
        else
            echo "WARNING: Empty recurrent gene pair found in ${GENEPAIRS}" >> ${REPORT}
        fi
    done < ${GENEPAIRS}
fi

final_count=$(tail -n +2 "${OUTPUT_FILE}" | wc -l)
echo "Final recurrent FusionInspector-validated fusions for this sample ${SAMPLENAME}: $final_count" >> ${REPORT}

# Filter AGFusion directories
mkdir -p validated-recurrent-agfusion-outdir

if [ $final_count -gt 0 ]; then
    # Ensure destination directory exists
    mkdir -p validated-recurrent-agfusion-outdir
    
    while IFS= read -r recurrent_gene_pair; do
        # Transform the gene pair format to match directory naming
        # PARG::BMS1__10:49885203-10:42791627 -> 26T_PARG--BMS1__10-49885203--10-42791627
        echo "Processing recurrent gene pair: $recurrent_gene_pair"
        transformed_pair=$(echo "$recurrent_gene_pair" | sed 's/-/--/g; s/::/--/g; s/:/-/g')
        echo "Transformed recurrent gene pair: $transformed_pair"
        expected_dir_pattern="${SAMPLENAME}_${transformed_pair}"
        
        for agf_dir in validated-agfusion-outdir/*; do
            if [ -d "$agf_dir" ]; then
                dir_basename=$(basename "$agf_dir")
                echo "Checking directory: $dir_basename"
                # Match the transformed pattern
                if [[ "$dir_basename" == "$expected_dir_pattern" ]]; then
                    cp -r "$agf_dir" validated-recurrent-agfusion-outdir/
                    echo "COPIED: $dir_basename" >> ${REPORT}
                    echo "COPIED: $dir_basename"
                else
                    echo "SKIPPED: $dir_basename does not match expected pattern $expected_dir_pattern" >> ${REPORT}
                    echo "SKIPPED: $dir_basename does not match expected pattern $expected_dir_pattern"
                fi
            fi
        done
    done < <(tail -n +2 "${OUTPUT_FILE}" | cut -f1)

    # Handle case where no directories exist
    if [ ! -d "validated-agfusion-outdir" ] || [ -z "$(ls -A validated-agfusion-outdir 2>/dev/null)" ]; then
      echo "WARNING: Input AGFusion directories are empty or do not exist." >> ${REPORT}
    fi
fi

# Final report
echo "Recurrent Fusion Filtering Completed for ${SAMPLENAME}" >> ${REPORT}
echo "===================================================" >> ${REPORT}
echo "Final recurrent validated fusions file: ${OUTPUT_FILE}" >> ${REPORT}
if [ $validated_count -gt 0 ]; then
    reduction_pct=$(( (validated_count - final_count) * 100 / validated_count ))
    echo "Reduction: ${reduction_pct}% filtered out" >> ${REPORT}
else
    echo "Reduction: N/A (no input fusions)" >> ${REPORT}
fi

