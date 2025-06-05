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

FusionInspector --fusions "$FUSIONLIST" \
                --genome_lib "${CTATDB}" \
                --left_fq "${READ1}" --right_fq "${READ2}" \
                --out_prefix "${SAMPLEID}" \
                --predict_cosmic_like \
                --cleanup
				