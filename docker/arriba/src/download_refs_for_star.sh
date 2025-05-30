#!/bin/bash
# this script is a modified version of the script from Arriba package `download_references.sh`
set -o pipefail
set -e -u

declare -A ASSEMBLIES
ASSEMBLIES[hs37d5]="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
ASSEMBLIES[hg19]="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz"
ASSEMBLIES[GRCh37]="http://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
ASSEMBLIES[hg38]="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chromFa.tar.gz"
ASSEMBLIES[GRCh38]="http://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
ASSEMBLIES[mm10]="http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/chromFa.tar.gz"
ASSEMBLIES[GRCm38]="http://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
ASSEMBLIES[mm39]="http://hgdownload.cse.ucsc.edu/goldenpath/mm39/bigZips/mm39.chromFa.tar.gz"
ASSEMBLIES[GRCm39]="http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

declare -A ANNOTATIONS
ANNOTATIONS[GENCODE19]="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
ANNOTATIONS[RefSeq_hg19]="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz"
ANNOTATIONS[ENSEMBL87]="http://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz"
ANNOTATIONS[GENCODE38]="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"
ANNOTATIONS[RefSeq_hg38]="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refGene.txt.gz"
ANNOTATIONS[ENSEMBL104]="http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz"
ANNOTATIONS[ENSEMBL113]="http://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.chr.gtf.gz"
ANNOTATIONS[GENCODEM25]="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
ANNOTATIONS[RefSeq_mm10]="http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/refGene.txt.gz"
ANNOTATIONS[GENCODEM27]="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.annotation.gtf.gz"
ANNOTATIONS[RefSeq_mm39]="http://hgdownload.cse.ucsc.edu/goldenpath/mm39/database/refGene.txt.gz"

declare -A COMBINATIONS
COMBINATIONS["hs37d5+GENCODE19"]="hs37d5+GENCODE19"
COMBINATIONS["hs37d5+RefSeq"]="hs37d5+RefSeq_hg19"
COMBINATIONS["hs37d5+ENSEMBL87"]="hs37d5+ENSEMBL87"
COMBINATIONS["hg19+GENCODE19"]="hg19+GENCODE19"
COMBINATIONS["hg19+RefSeq"]="hg19+RefSeq_hg19"
COMBINATIONS["hg19+ENSEMBL87"]="hg19+ENSEMBL87"
COMBINATIONS["GRCh37+GENCODE19"]="GRCh37+GENCODE19"
COMBINATIONS["GRCh37+RefSeq"]="GRCh37+RefSeq_hg19"
COMBINATIONS["GRCh37+ENSEMBL87"]="GRCh37+ENSEMBL87"
COMBINATIONS["hg38+GENCODE38"]="hg38+GENCODE38"
COMBINATIONS["hg38+RefSeq"]="hg38+RefSeq_hg38"
COMBINATIONS["hg38+ENSEMBL104"]="hg38+ENSEMBL104"
COMBINATIONS["GRCh38+GENCODE38"]="GRCh38+GENCODE38"
COMBINATIONS["GRCh38+RefSeq"]="GRCh38+RefSeq_hg38"
COMBINATIONS["GRCh38+ENSEMBL104"]="GRCh38+ENSEMBL104"
COMBINATIONS["GRCh38+ENSEMBL113"]="GRCh38+ENSEMBL113"
COMBINATIONS["GRCm38+GENCODEM25"]="GRCm38+GENCODEM25"
COMBINATIONS["GRCm38+RefSeq"]="GRCm38+RefSeq_mm10"
COMBINATIONS["mm10+GENCODEM25"]="mm10+GENCODEM25"
COMBINATIONS["mm10+RefSeq"]="mm10+RefSeq_mm10"
COMBINATIONS["GRCm39+GENCODEM27"]="GRCm39+GENCODEM27"
COMBINATIONS["GRCm39+RefSeq"]="GRCm39+RefSeq_mm39"
COMBINATIONS["mm39+GENCODEM27"]="mm39+GENCODEM27"
COMBINATIONS["mm39+RefSeq"]="mm39+RefSeq_mm39"
for COMBINATION in ${!COMBINATIONS[@]}; do
	COMBINATIONS["${COMBINATION%+*}viral+${COMBINATION#*+}"]="${COMBINATIONS[$COMBINATION]%+*}viral+${COMBINATIONS[$COMBINATION]#*+}"
done

if [ -z "$1" ] || [ -z "${COMBINATIONS[$1]}" ]; then
	echo "Usage: $(basename $0) ASSEMBLY+ANNOTATION" 1>&2
	echo "Available assemblies and annotations:" 1>&2
	tr ' ' '\n' <<<"${!COMBINATIONS[@]}" | sort 1>&2
	exit 0
fi

DESIRED_REF=$1
ASSEMBLY=$2
ANNOTATION=$3

# set defaults
THREADS="${THREADS-8}"
SJDBOVERHANG="${SJDBOVERHANG-250}"
WGET=$(which wget 2> /dev/null && echo " -q -O -" || echo "curl -L -s -S")

# now check if the assembly and annotation are provided locally
if [ ! -z $ASSEMBLY ] && [ ! -z $ANNOTATION ]; then
	echo "Using provided assembly and annotation"
	
	ASSEMBLY_NAME="${COMBINATIONS[${DESIRED_REF}]%+*}"
	ANNOTATION_NAME="${COMBINATIONS[${DESIRED_REF}]#*+}"
	VIRAL=""
	# check for viral parameter
	if [[ $DESIRED_REF =~ viral ]]; then
		VIRAL="viral"
	fi

	# check zipped
	if [[ ${ASSEMBLY} =~ \.tar\.gz$ ]]; then
		tar -x -O -z ${ASSEMBLY}
	elif [[ ${ASSEMBLY} =~ \.gz$ ]]; then
		gunzip -c ${ASSEMBLY}
	else
		cat ${ASSEMBLY}
	fi | 
	if [ "$VIRAL" = "viral" ]; then
		# drop viral contigs from assembly
		awk '/^>/{ contig=$1 } contig!~/^>NC_|^>AC_/{ print }'
	else
		cat
	fi > "/home/app/refs/${ASSEMBLY_NAME}.fa"

	echo "Using provided annotation"
	if [[ ${ANNOTATION} =~ \.gz$ ]]; then
		gunzip -c ${ANNOTATION}
	else
		cat ${ANNOTATION}
	fi |
	if [[ $ANNOTATION_NAME =~ RefSeq ]]; then
		# convert genePred to GTF
		awk -F '\t' -v OFS='\t' '
		function min(x, y) { return (x>y) ? y : x }
		function max(x, y) { return (x<y) ? y : x }
		{
			split($10, start, ",")
			split($11, end, ",")
			split($16, frame, ",")
			# remove stop codon from left end for coding genes on the minus strand
			if ($4=="-" && $14=="cmpl" && (start[1]!=$7 || (min(end[1],$8)-start[1]+frame[1])%3==0)) {
				$7+=3
				for (i in end)
					if ($7>=end[i] && $7<=end[i]+2)
						$7+=start[i+1]-end[i]
			}
			# remove stop codon from right end for coding genes on the plus strand
			if ($4=="+" && $15=="cmpl" && (end[$9]!=$8 || (end[$9]-max(start[$9],$7)+frame[$9])%3==0)) {
				$8-=3
				for (i in start)
					if ($8<=start[i] && $8>=start[i]-2)
						$8-=start[i]-end[i-1]
			}
			# append running number to duplicate uses of the same transcript ID
			gene_id=$13
			if (transcripts[$2]++) {
				gene_id=$13"_"transcripts[$2]
				$2=$2"_"transcripts[$2]
			}
			# print one line for each exon
			for (i=1; i<=$9; i++) {
				exon=($4=="+") ? i : $9-i+1
				attributes="gene_id \""gene_id"\"; transcript_id \""$2"\"; exon_number \""exon"\"; exon_id \""$2"."exon"\"; gene_name \""$13"\";"
				print $3,"RefSeq","exon",start[i]+1,end[i],".",$4,".",attributes
				# print one line for each coding region
				if ($14~/cmpl/ && $7<=end[i] && $8>=start[i])
					print $3,"RefSeq","CDS",max($7,start[i])+1,min($8,end[i]),".",$4,frame[i],attributes
			}
		}' | sort -k1,1V -k4,4n -k5,5n -k3,3 -S4G 
	else
		cat
	fi |
	if ! grep -q '^>chr' "/home/app/refs/${ASSEMBLY_NAME}.fa"; then
		sed -e 's/^chrM/MT/' -e 's/^chr//'
	else
		sed -e 's/^MT/chrM/' -e 's/^\([1-9XY]\|[12][0-9]\)/chr\1/'
	fi > "/home/app/refs/${ANNOTATION_NAME}.gtf"

else
	ASSEMBLY="${COMBINATIONS[$1]%+*}"
	ANNOTATION="${COMBINATIONS[$1]#*+}"
	VIRAL=""
	# check for viral parameter
	if [[ $ASSEMBLY =~ viral ]]; then
		ASSEMBLY="${ASSEMBLY%viral}"
		VIRAL="viral"
	fi

	echo "Downloading assembly: ${ASSEMBLIES[$ASSEMBLY]}"
	$WGET "${ASSEMBLIES[$ASSEMBLY]}" |
	if [[ ${ASSEMBLIES[$ASSEMBLY]} =~ \.tar\.gz$ ]]; then
		tar -x -O -z
	elif [[ ${ASSEMBLIES[$ASSEMBLY]} =~ \.gz$ ]]; then
		gunzip -c
	else
		cat
	fi |
	if [ "$VIRAL" = "viral" ]; then
		# drop viral contigs from assembly
		awk '/^>/{ contig=$1 } contig!~/^>NC_|^>AC_/{ print }'
	else
		cat
	fi > "/home/app/refs/$ASSEMBLY.fa"

	echo "Downloading annotation: ${ANNOTATIONS[$ANNOTATION]}"
	$WGET "${ANNOTATIONS[$ANNOTATION]}" |
	if [[ ${ANNOTATIONS[$ANNOTATION]} =~ \.gz$ ]]; then
		gunzip -c
	else
		cat
	fi |
	if [[ $ANNOTATION =~ RefSeq ]]; then
		# convert genePred to GTF
		awk -F '\t' -v OFS='\t' '
		function min(x, y) { return (x>y) ? y : x }
		function max(x, y) { return (x<y) ? y : x }
		{
			split($10, start, ",")
			split($11, end, ",")
			split($16, frame, ",")
			# remove stop codon from left end for coding genes on the minus strand
			if ($4=="-" && $14=="cmpl" && (start[1]!=$7 || (min(end[1],$8)-start[1]+frame[1])%3==0)) {
				$7+=3
				for (i in end)
					if ($7>=end[i] && $7<=end[i]+2)
						$7+=start[i+1]-end[i]
			}
			# remove stop codon from right end for coding genes on the plus strand
			if ($4=="+" && $15=="cmpl" && (end[$9]!=$8 || (end[$9]-max(start[$9],$7)+frame[$9])%3==0)) {
				$8-=3
				for (i in start)
					if ($8<=start[i] && $8>=start[i]-2)
						$8-=start[i]-end[i-1]
			}
			# append running number to duplicate uses of the same transcript ID
			gene_id=$13
			if (transcripts[$2]++) {
				gene_id=$13"_"transcripts[$2]
				$2=$2"_"transcripts[$2]
			}
			# print one line for each exon
			for (i=1; i<=$9; i++) {
				exon=($4=="+") ? i : $9-i+1
				attributes="gene_id \""gene_id"\"; transcript_id \""$2"\"; exon_number \""exon"\"; exon_id \""$2"."exon"\"; gene_name \""$13"\";"
				print $3,"RefSeq","exon",start[i]+1,end[i],".",$4,".",attributes
				# print one line for each coding region
				if ($14~/cmpl/ && $7<=end[i] && $8>=start[i])
					print $3,"RefSeq","CDS",max($7,start[i])+1,min($8,end[i]),".",$4,frame[i],attributes
			}
		}' | sort -k1,1V -k4,4n -k5,5n -k3,3 -S4G
	else
		cat
	fi |
	if ! grep -q '^>chr' "/home/app/refs/$ASSEMBLY.fa"; then
		sed -e 's/^chrM/MT/' -e 's/^chr//'
	else
		sed -e 's/^MT/chrM/' -e 's/^\([1-9XY]\|[12][0-9]\)/chr\1/'
	fi > "/home/app/refs/$ANNOTATION.gtf"

fi

if [ "$VIRAL" = "viral" ]; then
	echo "Appending RefSeq viral genomes"
	REFSEQ_VIRAL_GENOMES=$(dirname "$0")/RefSeq_viral_genomes_v2.3.0.fa.gz
	if [ ! -e "$REFSEQ_VIRAL_GENOMES" ]; then
		echo "Viral genomes from RefSeq not found in the expected location: $REFSEQ_VIRAL_GENOMES"
		exit 1
	fi
	if [ ! -z ${ASSEMBLY_NAME} ]; then
		gunzip -c "$REFSEQ_VIRAL_GENOMES" >> "/home/app/refs/${ASSEMBLY_NAME}.fa"
	else
		gunzip -c "$REFSEQ_VIRAL_GENOMES" >> "/home/app/refs/${ASSEMBLY}.fa"
	fi
fi

if [[ $(samtools --version-only 2> /dev/null) =~ ^1\. ]]; then
	echo "Indexing assembly"
	if [ ! -z "$ASSEMBLY_NAME" ]; then
		samtools faidx "/home/app/refs/$ASSEMBLY_NAME.fa"
	else
		samtools faidx "/home/app/refs/$ASSEMBLY.fa"
	fi
fi

if [[ ! -z $2 ]]; then
	mkdir -p /home/app/refs/STAR_index_${ASSEMBLY_NAME}_${ANNOTATION_NAME}
	if STAR --runMode genomeGenerate --genomeDir /home/app/refs/STAR_index_${ASSEMBLY_NAME}_${ANNOTATION_NAME} --genomeFastaFiles "/home/app/refs/${ASSEMBLY_NAME}.fa" --sjdbGTFfile "/home/app/refs/$ANNOTATION_NAME.gtf" --runThreadN "$THREADS" --sjdbOverhang "$SJDBOVERHANG"; then
		echo "STAR index created successfully"
	else
		echo "STAR index creation failed"
		exit 1
	fi
else
	mkdir -p /home/app/refs/STAR_index_${ASSEMBLY}_${ANNOTATION}
	if STAR --runMode genomeGenerate --genomeDir /home/app/refs/STAR_index_${ASSEMBLY}_${ANNOTATION} --genomeFastaFiles "/home/app/refs/$ASSEMBLY.fa" --sjdbGTFfile "/home/app/refs/$ANNOTATION.gtf" --runThreadN "$THREADS" --sjdbOverhang "$SJDBOVERHANG"; then
		echo "STAR index created successfully"
	else
		echo "STAR index creation failed"
		exit 1
	fi
fi

echo "Reference for STAR and Arriba downloaded and indexed successfully."
