#!/usr/bin/bash

echo $(id)

READ1=$1
READ2=$2
CORES=$3
SAMPLE_ID="${READ1%%_*}"
CTAT_LIB="/work/libs"
OUTDIR="/work/nf_work"

if STAR \
--runMode alignReads \
--runThreadN "${CORES}" \
--genomeDir "${CTAT_LIB}/ref_genome.fa.star.idx" \
--readFilesIn "${READ1}" "${READ2}" \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--limitBAMsortRAM 20000000000 \
--outFileNamePrefix "${OUTDIR}/${SAMPLE_ID}-STAR_" \
--genomeLoad LoadAndRemove >> "${OUTDIR}/mapping.log.txt" 2>&1; then
    samtools index "${OUTDIR}/${SAMPLE_ID}-STAR_Aligned.sortedByCoord.out.bam" -@ "${CORES}" -o "${OUTDIR}/${SAMPLE_ID}-STAR_Aligned.sortedByCoord.out.bam.bai" && echo "STAR alignment complete."
fi

