//
process ALIGN_READS_STAR_ARRIBA {
    cpus params.numCores
    
    label 'alignReadsArriba'

    container "${params.container__preproc}"

    input:
        tuple val(sampleName), path(filtFastqs)
        path starIndex
    output:
        tuple val(sampleName), path("*-STAR-ARR_Aligned.out.bam", arity: '1'), emit: aligned_bam

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}
    CORES=${task.cpus}
    STAR_INDEX=${starIndex}
    READ1=${filtFastqs[0]}  # First file in the nested list will be read 1 file
    READ2=${filtFastqs[1]}

    echo "Processing files of sample \${SAMPLE_ID}"
    echo "Number of cores to use: ${task.cpus}"
    echo "The index path: \${STAR_INDEX}"
    
    # STAR normal alignment (1-pass) for Arriba
    if star-for-arriba--nf.sh "\${READ1}" "\${READ2}" "\${SAMPLE_ID}" ${task.cpus} "\${STAR_INDEX}"; then
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
