#!/usr/bin/bash
# Usage: ./fusins--nf.sh input_genepair.txt ctat_genome_lib_path

FUSIONLIST=$1
CTATDB=$2
READ1=$3
READ2=$4
SAMPLEID=$5
# Example usage:
# ./fusins--nf.sh input_genepair.txt /path/to/ctat_genome_lib/ rnaseq_1.fq rnaseq_2.fq 124T

# check arguments
if [ $# -ne 5 ]; then
	echo "Usage: $0 <input_genepair.txt> <ctat_genome_lib_path> <read1.fq> <read2.fq> <sample_id>"
	exit 1
fi

# Capture both stdout/stderr and exit code
output=$(FusionInspector --fusions "$FUSIONLIST" \
                --genome_lib "${CTATDB}" \
                --left_fq "${READ1}" --right_fq "${READ2}" \
                --out_prefix "${SAMPLEID}" \
                --predict_cosmic_like \
                --cleanup 2>&1)

exit_code=$?

# Print the output
echo "$output"

# Check for specific messages and exit codes
if [[ "$output" == "No fusions listed in input file"* ]]; then
    echo "INFO: No fusions found in input file...creating empty output file"
    mkdir -p FI && touch "FI/${SAMPLEID}.FusionInspector.fusions.abridged.tsv" && printf "#FusionName\tJunctionReadCount\tSpanningFragCount\test_J\test_S\tLeftGene\tLeftLocalBreakpoint\tLeftBreakpoint\tRightGene\tRightLocalBreakpoint\tRightBreakpoint\tSpliceType\tLargeAnchorSupport\tNumCounterFusionLeft\tNumCounterFusionRight\tFAR_left\tFAR_right\tLeftBreakDinuc\tLeftBreakEntropy\tRightBreakDinuc\tRightBreakEntropy\tFFPM\tmicroh_brkpt_dist\tnum_microh_near_brkpt\tannot_splice\tconsensus_splice\tleft_counter_ffpm\tright_counter_ffpm\tpred_cluster\tfusion_cluster_att\tannots\n" > "FI/${SAMPLEID}.FusionInspector.fusions.abridged.tsv"
    exit 0
elif [ $exit_code -eq 0 ]; then
    echo "FusionInspector completed successfully"
    exit 0
else
    echo "ERROR: FusionInspector failed with exit code: $exit_code"
    exit $exit_code
fi
