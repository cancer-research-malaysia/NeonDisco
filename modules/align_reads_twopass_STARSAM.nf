//
process ALIGN_READS_TWOPASS_STARSAM {
    errorStrategy 'retry'
    maxRetries 3
    maxForks 5
    cpus params.numCores

    label 'alignReads2Pass'
    
    container "${params.container__preproc}"

    input:
        tuple val(sampleName), path(trimmedReads)
        path starIndex

    output:
        tuple val(sampleName), path("*-STAR*Aligned.out.bam", arity: '1'), emit: aligned_bam

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}
    CORES=${task.cpus}
    STAR_INDEX=${starIndex}
    READ1=${trimmedReads[0]}  # First file in the nested list will be read 1 file
    READ2=${trimmedReads[1]}

    echo "Processing files of sample \${SAMPLE_ID}"
    echo "Number of cores to use: ${task.cpus}"
    echo "The index path: \${STAR_INDEX}"

    # STAR 2-pass alignment
    if star-2pass--nf.sh "\${READ1}" "\${READ2}" "\${SAMPLE_ID}" ${task.cpus} "\${STAR_INDEX}"; then
       echo "STAR sample-level 2-pass alignment is complete!"
    else
       echo "STAR alignment failed. Check logs. Exiting..."
       exit 1
    fi

    """
    stub:
    """
    touch ${sampleName}-STAR_2pass_Aligned.out.bam
    touch ${sampleName}-STAR_2pass_Log.final.out
    touch ${sampleName}-STAR_2pass_SJ.out.tab
    echo "Stub run finished!" > test_stub_STAR-2pass-align.log
    """
}
