#!/usr/bin/env bash

# Usage: ./fusins--nf.sh input_genepair.txt ctat_genome_lib_path read1.fq read2.fq sample_id cores
FUSIONLIST=$1
CTATDB=$2
READ1=$3
READ2=$4
SAMPLEID=$5
CORES=$6

# Function to create empty output file with proper header
create_empty_output() {
    local sample_id=$1
    local reason=$2
    
    echo "INFO: ${reason}...creating empty output file"
    mkdir -p FI
    
    # FusionInspector abridged TSV header
    cat > "FI/${sample_id}.FusionInspector.fusions.abridged.tsv" <<-'EOF'
	#FusionName	JunctionReadCount	SpanningFragCount	est_J	est_S	LeftGene	LeftLocalBreakpoint	LeftBreakpoint	RightGene	RightLocalBreakpoint	RightBreakpoint	SpliceType	LargeAnchorSupport	NumCounterFusionLeft	NumCounterFusionRight	FAR_left	FAR_right	LeftBreakDinuc	LeftBreakEntropy	RightBreakDinuc	RightBreakEntropy	FFPM	microh_brkpt_dist	num_microh_near_brkpt	annot_splice	consensus_splice	left_counter_ffpm	right_counter_ffpm	pred_cluster	fusion_cluster_att	annots
	EOF
}

# Run FusionInspector and capture output
output=$(FusionInspector --fusions "$FUSIONLIST" \
                --genome_lib "${CTATDB}" \
                --left_fq "${READ1}" --right_fq "${READ2}" \
                --out_prefix "${SAMPLEID}" \
                --predict_cosmic_like \
                --cleanup --CPU "$CORES" 2>&1)

exit_code=$?

# Print the output
echo "$output"

# Handle different exit scenarios
if [[ "$output" == *"No fusions listed in input file"* ]]; then
    create_empty_output "$SAMPLEID" "No fusions found in input file"
    exit 0
elif [ $exit_code -eq 0 ]; then
    echo "FusionInspector completed successfully."
    exit 0
else
    echo "ERROR: FusionInspector failed with exit code: $exit_code"
    create_empty_output "$SAMPLEID" "Irrecoverable error"
    exit $exit_code
fi
