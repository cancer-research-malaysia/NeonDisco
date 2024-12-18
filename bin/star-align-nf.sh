#!/usr/bin/bash

SAMPLE_ID="${READ1%%_*}"
CTAT_LIB="/work/libs"
INPUTDIR="/work/input"
OUTDIR="/work/nf_work"

STAR \
--runMode alignReads \
--runThreadN 16 \
--genomeDir "${CTAT_LIB}/ref_genome.fa.star.idx" \
--readFilesIn "${INPUTDIR}"/12T_R1.fq "${INPUTDIR}"/12T_R2.fq \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--limitBAMsortRAM 20000000000 \
--outFileNamePrefix "${OUTDIR}/${SAMPLE_ID}-STAR_" \
--genomeLoad NoSharedMemory >> mapping.log.txt 2>&1

# samtools index Aligned.sortedByCoord.out.bam -@ 16 -o Aligned.sortedByCoord.out.bai

#/home/suffian/libs/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx

