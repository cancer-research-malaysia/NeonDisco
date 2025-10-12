//
process CONVERT_FILTREADS_BAM2FASTQ_EASYFUSE {
    errorStrategy 'retry'
    maxRetries 3
    cpus params.numCores
    
    label 'convertFilteredReads'

    container "${params.container__preproc}"

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
    
    samtools fastq -1 "\${SAMPLE_ID}-filtered.read1.fastq.gz" -2 "\${SAMPLE_ID}-filtered.read2.fastq.gz" -0 /dev/null -@ ${task.cpus} \${BAM}

    """
    stub:
    """
    touch ${sampleName}-filtered.read1.fastq.gz
    touch ${sampleName}-filtered.read2.fastq.gz
    echo "Stub run finished!"
    """
}
