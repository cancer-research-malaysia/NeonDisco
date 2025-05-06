//
process ALIGN_READS_STAR_GENERAL {
    maxForks 1
    
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name ALIGNMENT-STARGENERAL -v ${params.starIndex}:/home/app/starIdx -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"

    input:
        tuple val(sampleName), path(trimmedReads)

    output:
        tuple val(sampleName), path("*-STAR-GEN_Aligned.out.bam", arity: '1'), emit: aligned_bams
        tuple val(sampleName), path("*-STAR-GEN_Chimeric.out.junction"), emit: chimeric_reads
        tuple val(sampleName), path("*-STAR-GEN_Log.final.out"), emit: read_stats

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}
    CORES=${params.numCores}
    STAR_INDEX="/home/app/starIdx"
    READ1=${trimmedReads[0]}  # First file in the nested list will be read 1 file
    READ2=${trimmedReads[1]}

    echo "Processing files of sample \${SAMPLE_ID}"
    echo "Number of cores to use: ${params.numCores}"
    echo "The index path: \${STAR_INDEX}"
    
    # STAR normal alignment for general tools
    if bash /home/app/scripts/star-general-nf.sh "\${READ1}" "\${READ2}" "\${SAMPLE_ID}" ${params.numCores} "\${STAR_INDEX}"; then
        echo "STAR general alignment is complete!"
    else
        echo "STAR alignment failed. Check logs. Exiting..."
        exit 1
    fi

    """
    stub:
    """
    touch ${sampleName}-STAR-GEN_Aligned.out.bam
    touch ${sampleName}-STAR-GEN_Log.final.out
    touch ${sampleName}-STAR-GEN_Chimeric.out.junction
    echo "Stub run finished!" > test_stub_STAR-Arriba-align.log
    """
}
