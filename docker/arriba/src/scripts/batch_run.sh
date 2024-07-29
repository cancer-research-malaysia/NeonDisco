#!/usr/bin/bash

run_star_and_arriba() {
  local CTAT_LIB=$1
  local FASTQS=$2
  local SAMPLE=$3
  local ARRIBA_PKG=$4
  local ARR_OUTDIR=$5
  local STAR_TMPDIR=$6

  STAR --runThreadN 8 \
  --genomeDir "${CTAT_LIB}/ref_genome.fa.star.idx" \
  --genomeLoad NoSharedMemory \
  --readFilesIn "${FASTQS}/${SAMPLE}_r1.fq.gz" "${FASTQS}/${SAMPLE}_r2.fq.gz" \
  --readFilesCommand zcat \
  --outStd BAM_Unsorted \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped Within \
  --outBAMcompression 0 \
  --outFilterMultimapNmax 50 \
  --peOverlapNbasesMin 10 \
  --alignSplicedMateMapLminOverLmate 0.5 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --chimSegmentMin 10 \
  --chimOutType WithinBAM HardClip \
  --chimJunctionOverhangMin 10 \
  --chimScoreDropMax 30 \
  --chimScoreJunctionNonGTAG 0 \
  --chimScoreSeparation 1 \
  --chimSegmentReadGapMax 3 \
  --chimMultimapNmax 50 \
  --outTmpDir "${STAR_TMPDIR}" | tee "${ARR_OUTDIR}/${SAMPLE}-Aligned.out.bam" | arriba \
  -x /dev/stdin \
  -o "${ARR_OUTDIR}/${SAMPLE}/arriba-fusions.tsv" \
  -O "${ARR_OUTDIR}/${SAMPLE}/fusions.discarded.tsv" \
  -a "${CTAT_LIB}/ref_genome.fa" \
  -g "${CTAT_LIB}/ref_annot.gtf" \
  -b "${ARRIBA_PKG}/blacklist_hg38_GRCh38_v2.3.0.tsv.gz" \
  -k "${ARRIBA_PKG}/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz" \
  -t "${ARRIBA_PKG}/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz" \
  -p "${ARRIBA_PKG}/protein_domains_hg38_GRCh38_v2.3.0.gff3" 2>&1 | tee "${ARR_OUTDIR}/${SAMPLE}/arriba-test.log.txt"
}

# Example usage:
# run_star_and_arriba "/path/to/CTAT_LIB" "/path/to/FASTQS" "/path/to/ARRIBA_PKG" "/path/to/ARR_OUTDIR" "/path/to/STAR_TMPDIR"

# Set env variables
CTAT_LIB="/tmp/libs"
FASTQS="/tmp/data"
ARR_OUTDIR="/tmp/out"
STAR_TMPDIR="/tmp/out/star"
ARRIBA_PKG="/opt/conda/var/lib/arriba"

echo "Environment variables set! Listing fastq files..."

# Check if the subdirectory contains a .fastq.gz or fq.gz file
if [[ $(find "$FASTQS" \( -name '*.fastq.gz' -o -name '*.fq.gz' \) -type f) ]]; then
    echo "Fastq files are present. Listing..."
    find "$FASTQS" -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \)
    echo "Total number of files: $(find "$FASTQS" -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | wc -l)"
    # extract read1 group number of files
    R1_COUNT=$(find ./ -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | xargs -n 1 basename | awk -F'_' '/r1/ {print $1}' | sort | uniq | wc -l)
    # extract read2 group number of files
    R2_COUNT=$(find ./ -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | xargs -n 1 basename | awk -F'_' '/r2/ {print $1}' | sort | uniq | wc -l)
    # extract sample IDs into an array
    SAMPLE_ID=($(find ./ -mindepth 1 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' \) | xargs -n 1 basename | awk -F'_' '/r2/ {print $1}' | sort | uniq))

    # Check if R1_COUNT == R2_COUNT (implying paired reads data)
    if (( R1_COUNT == R2_COUNT )); then
        echo "The input files appear to be paired. Looping through the sample ID array..."
        for prefix in "${SAMPLE_ID[@]}"; do
            echo "Sample ID: ${prefix}"
            echo "Running STAR while piping to Arriba..."
            mkdir -p "${ARR_OUTDIR}/${prefix}"
            # Measure execution time
            STARTTIME=$(date +%s)
            if run_star_and_arriba "${CTAT_LIB}" "${FASTQS}" "${prefix}" "${ARRIBA_PKG}" "${ARR_OUTDIR}" "${STAR_TMPDIR}"; then
                ENDTIME=$(date +%s)
                ELAP=$(( ENDTIME - STARTTIME ))
                echo "Arriba run completed successfully. Time taken: ${ELAP}. Check log file for run details."
            else
                echo "Something went wrong during Arriba run. Check log file."
            fi
        done
    else
        echo "The number of files is odd. Make sure the input files are paired before proceeding."
        exit 1
    fi
else
    echo "Fastq input files are not found. Double check your input path."
    exit 1
fi
