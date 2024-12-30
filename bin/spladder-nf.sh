#!/usr/bin/bash

echo $(id)

BAMFILE=$1
THREADS=$2
CTAT_LIB="/work/libs"
OUTDIR=$3

if spladder build --parallel ${THREADS} --bams ${BAMFILE} --annotation ${CTAT_LIB}/ref_annot.gtf --outdir ${OUTDIR} --set-mm-tag nM -v 2>&1 | tee "${OUTDIR}/spladder-run-nf.log-$(date +%Y%m%d_%H-%M-%S).txt"; then
    echo "Spladder build done!"
else
    echo "Spladder build failed."
    exit 1
fi

