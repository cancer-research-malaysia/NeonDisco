//
process ALIGN_READS_STAR_ARRIBA {
    errorStrategy 'retry'
    maxRetries 3
    maxForks 1
    
    label 'alignReadsArriba'

    container "${params.container__preproc}"

    input:
        tuple val(sampleName), path(trimmedReads)
        path starIndex
    output:
        tuple val(sampleName), path("*-STAR-ARR_Aligned.out.bam", arity: '1'), emit: aligned_bam

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}
    CORES=${params.numCores}
    STAR_INDEX=${starIndex}
    READ1=${trimmedReads[0]}  # First file in the nested list will be read 1 file
    READ2=${trimmedReads[1]}

    echo "Processing files of sample \${SAMPLE_ID}"
    echo "Number of cores to use: ${params.numCores}"
    echo "The index path: \${STAR_INDEX}"
    
    # STAR normal alignment (1-pass) for Arriba
    if star-for-arriba--nf.sh "\${READ1}" "\${READ2}" "\${SAMPLE_ID}" ${params.numCores} "\${STAR_INDEX}"; then
        echo "STAR normal 1-pass alignment for Arriba is complete!"
    else
        echo "STAR alignment failed. Check logs. Exiting..."
        exit 1
    fi

    """
    stub:
    """
    touch ${sampleName}-STAR-ARR_Aligned.out.bam
    echo "Stub run finished!" > test_stub_STAR-Arriba-align.log
    """
}
