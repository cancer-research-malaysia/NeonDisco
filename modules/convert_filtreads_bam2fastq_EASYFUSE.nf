//
process CONVERT_FILTREADS_BAM2FASTQ_EASYFUSE {
    
    label 'convertFilteredReads'

    container "${params.container__preproc}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FILTERED-BAM-2-FASTQ"

    input:
        tuple val(sampleName), path(filteredBam)

    output:
        tuple val(sampleName), path("*-filtered.read?.fastq.gz", arity: '2'), emit: filtered_fastqs

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}
    BAM=${filteredBam}

    echo "Processing files of sample \${SAMPLE_ID}"
    echo "The filtered bam file: \${BAM}"
    
    samtools fastq \
    -0  "\${SAMPLE_ID}-filtered.other.fastq.gz" \
    -1  "\${SAMPLE_ID}-filtered.read1.fastq.gz" \
    -2 "\${SAMPLE_ID}-filtered.read2.fastq.gz" \
    --threads ${params.numCores} \
    \${BAM}

    """
    stub:
    """
    touch ${sampleName}-filtered.read1.fastq.gz
    touch ${sampleName}-filtered.read2.fastq.gz
    echo "Stub run finished!"
    """
}
