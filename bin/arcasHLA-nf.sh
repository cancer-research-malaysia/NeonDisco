#!/usr/bin/bash 

SAMPLE_ID=$1
BAM=$2
CORES=$3

echo "Running HLA-HD for sample: ${SAMPLE_ID}"
echo "Input BAM file: ${BAM}"
echo "Number of cores to use: ${CORES}"

# measure execution time
STARTTIME=$(date +%s)

# run arcasHLA extract command
# arcasHLA extract ${BAM} -t ${CORES} -v && arcasHLA genotype *.1.fq.gz *.2.fq.gz -t ${CORES} -v
# MHCI only run
arcasHLA extract ${BAM} -t ${CORES} -v && arcasHLA genotype *.1.fq.gz *.2.fq.gz --genes A,B,C -t ${CORES} -v

if [ $? -eq 0 ]; then
    echo "arcasHLA run successfully for sample ${SAMPLE_ID}."
    ENDTIME=$(date +%s)
    echo "arcasHLA run time: $(($ENDTIME - $STARTTIME)) seconds"
else
    echo "arcasHLA failed for sample ${SAMPLE_ID}."
    exit 1
fi

