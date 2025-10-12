//
process ALIGN_READS_STAR_GENERAL {
    errorStrategy 'retry'
    maxRetries 3
    //maxForks 10
    cpus params.numCores
    
    label 'alignReadsGeneral'
    
    container "${params.container__preproc}"

    input:
        tuple val(sampleName), path(trimmedReads)
        path starIndex

    output:
        tuple val(sampleName), path("${sampleName}-STAR-GEN_Aligned.out.sorted.bam"), path("${sampleName}-STAR-GEN_Aligned.out.sorted.bam.bai"), path("${sampleName}-STAR-GEN_Aligned.out.sorted.byqueryname.bam"), emit: final_bam

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
    
    # STAR normal alignment for general tools
    if star-general--nf.sh "\${READ1}" "\${READ2}" "\${SAMPLE_ID}" ${task.cpus} "\${STAR_INDEX}"; then
        echo "STAR general alignment is complete!"
        # sort and index - create files with the names Nextflow expects
        samtools sort -@ ${task.cpus} -m 4G -O bam -o \${SAMPLE_ID}-STAR-GEN_Aligned.out.sorted.bam \${SAMPLE_ID}-STAR-GEN_Aligned.out.bam
        samtools index -@ ${task.cpus} \${SAMPLE_ID}-STAR-GEN_Aligned.out.sorted.bam

        # now sort by queryname for fusion detection
        samtools sort -n -@ ${task.cpus} -o \${SAMPLE_ID}-STAR-GEN_Aligned.out.sorted.byqueryname.bam \${SAMPLE_ID}-STAR-GEN_Aligned.out.sorted.bam
    else
        echo "STAR alignment failed. Check logs. Exiting..."
        exit 1
    fi

    """
    stub:
    """
    touch ${sampleName}-STAR-GEN_Aligned.out.bam
    touch ${sampleName}-STAR-GEN_Aligned.out.sorted.bam
    touch ${sampleName}-STAR-GEN_Aligned.out.sorted.bam.bai
    touch ${sampleName}-STAR-GEN_Aligned.out.sorted.byqueryname.bam
    touch ${sampleName}-STAR-GEN_Log.final.out
    touch ${sampleName}-STAR-GEN_Chimeric.out.junction
    echo "Stub run finished!" > test_stub_STAR-Arriba-align.log
    """
}
