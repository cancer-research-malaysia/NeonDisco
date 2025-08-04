//
process ALIGN_READS_STAR_GENERAL {
    maxForks 1
    label 'alignReadsGeneral'
    
    container "${params.container__preproc}"

    input:
        tuple val(sampleName), path(trimmedReads)
        path starIndex

    output:
        tuple val(sampleName), path("*-STAR-GEN_Aligned.out.bam"), path("*-STAR-GEN_Aligned.out.bam.bai"), emit: final_bam

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
    
    # STAR normal alignment for general tools
    if star-general--nf.sh "\${READ1}" "\${READ2}" "\${SAMPLE_ID}" ${params.numCores} "\${STAR_INDEX}"; then
        echo "STAR general alignment is complete!"
        # sort and index
        samtools sort -@ ${params.numCores} -m 4G -O bam \${SAMPLE_ID}-STAR-GEN_Aligned.out.bam | samtools index -@ ${params.numCores} - && echo "Indexing complete!"
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
