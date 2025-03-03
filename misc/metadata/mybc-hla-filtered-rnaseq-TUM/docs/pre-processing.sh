workdir="/shared_2/minhui/mybrca/RNAseqHLAHD/STAR_mapping_hs37d5_2023" 
S3="s3://crm.tumorstudy.analysis"
NUM=$1
ID=$(sed -n ${NUM}p allMyBrCa_ExceptThoseWithoutHLAHD.txt)

#download fq files
mkdir ${workdir}/sample_${ID}
cd ${workdir}/sample_${ID}
if [[ $ID -lt 173 ]]; then ID_s3="${ID}T"; else ID_s3=$ID; fi  #account for inconsistent file naming
parallel -j 2 aws s3 cp $S3/mybrca-sequencing-data/RNAseq/${ID_s3}_r{}.fq.gz ${ID}_r{}.fq.gz ::: 1 2
parallel -j 2 gunzip ${ID}_r{}.fq.gz ::: 1 2

#START MAPPING + PREPROCESSING
echo "LINE_${NUM} SAMPLE_${ID} processing : START $(date)" >> ${workdir}/pre_process.progress.log.txt
STAR \
--runMode alignReads \
--runThreadN 16 \
--genomeDir /shared_2/minhui/ref/STAR_GenomeIndex_hs37d5_2023 \
--readFilesIn ${ID}_r1.fq ${ID}_r2.fq \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--limitBAMsortRAM 20000000000 \
--genomeLoad NoSharedMemory >> mapping.log.txt 2>&1

samtools index Aligned.sortedByCoord.out.bam -@ 16
samtools view -b -h Aligned.sortedByCoord.out.bam "6:28477797-33448354" -@ 16 > MHC.bam
samtools view -b -f 4 Aligned.sortedByCoord.out.bam -@ 16 > unmapped.bam
samtools merge -o merged.bam MHC.bam unmapped.bam -@ 16
samtools sort -n -@ 16 -m 2G -o merged_sorted.bam merged.bam
samtools fastq -@ 16 -1 merged_R1.fq -2 merged_R2.fq merged_sorted.bam

if [[ -e merged_R1.fq && -e merged_R2.fq ]]; then
	echo "LINE_${NUM} SAMPLE_${ID} processing : DONE  $(date)" >> ${workdir}/pre_process.progress.log.txt
	aws s3 cp merged_R1.fq $S3/mybrca-sequencing-data/RNAseq_MHCplusUnmapped_forHLAHD/${ID}T_R1.fq
	aws s3 cp merged_R2.fq $S3/mybrca-sequencing-data/RNAseq_MHCplusUnmapped_forHLAHD/${ID}T_R2.fq
	aws s3 cp mapping.log.txt $S3/mybrca-sequencing-data/RNAseq_MHCplusUnmapped_forHLAHD/${ID}T_STAR.log.txt
else
	echo "LINE_${NUM} SAMPLE_${ID} processing : FAILED $(date)" >> ${workdir}/pre_process.failed.txt
fi

cd ${workdir}
rm -rd sample_${ID}
