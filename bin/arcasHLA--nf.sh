#!/usr/bin/env bash

SAMPLE_ID=$1
BAM=$2
CORES=$3

echo "Running arcasHLA for sample: ${SAMPLE_ID}"
echo "Input BAM file: ${BAM}"
echo "Number of cores to use: ${CORES}"

# measure execution time
STARTTIME=$(date +%s)

# run arcasHLA extract command
# arcasHLA extract ${BAM} -t ${CORES} -v && arcasHLA genotype *.1.fq.gz *.2.fq.gz -t ${CORES} -v
# MHCI only run 
# then rename output json file to match sample ID
arcasHLA extract ${BAM} -t ${CORES} -v && arcasHLA genotype *.1.fq* *.2.fq* --genes A,B,C -t ${CORES} -v && mv ${SAMPLE_ID}*.genotype.json ${SAMPLE_ID}.genotype.json
if [ $? -eq 0 ]; then
    echo "arcasHLA run successfully for sample ${SAMPLE_ID}."
    ENDTIME=$(date +%s)
    echo "arcasHLA run time: $(($ENDTIME - $STARTTIME)) seconds"
else
    echo "arcasHLA failed for sample ${SAMPLE_ID}."
    exit 1
fi

