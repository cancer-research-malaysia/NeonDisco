// 
process FIXMATES_MARKDUPES_SAMTOOLS {
    maxForks 2
    
    // Apply publishDir conditionally based on whether outputDir starts with 's3://'
    if (params.outputDir.toString().startsWith('s3://')) {
        publishDir "${params.outputDir}/${sampleName}/2PASS-ALIGNMENT-out", 
            mode: 'copy',
            saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    }
    
    afterScript params.deleteIntMedFiles ? "find ./ -name \"${sampleName}*STAR_2pass_Aligned.out.bam\" -type l -exec sh -c 'rm -f \$(readlink -f \"{}\")' \\; -delete" : "echo 'Skipping intermediate file cleanup...'"
    
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name BAM-POSTPROC -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(bamFile)
    
    output:
        tuple val(sampleName), path("*_fixmates_markdupes.bam"), path("*_fixmates_markdupes.bam.bai"), emit: final_bam
    
    script:
    """
    samtools sort -n -@ ${params.numCores} -m 4G -O bam ${bamFile} | \
    samtools fixmate -pcmu -O bam - ${sampleName}_fixmates.bam

    if [ -f ${sampleName}_fixmates.bam ]; then
        echo "Fixmate information for ${sampleName} is complete!"
        # mark duplicates
        samtools sort -@ ${params.numCores} -m 4G -O bam ${sampleName}_fixmates.bam | \
        samtools markdup -@ ${params.numCores} - ${sampleName}_fixmates_markdupes.bam

        # then index
        samtools index ${sampleName}_fixmates_markdupes.bam && echo "Indexing complete!"

        # clean up intermediate bam
        rm -f ${sampleName}_fixmates.bam
    else
        echo "Grokking fixmate information for ${sampleName} failed. Check logs. Exiting..."
        exit 1
    fi
    """
    stub:
    """
    touch ${sampleName}_fixmates_markdupes.bam
    touch ${sampleName}_fixmates_markdupes.bam.bai
    echo "Stub run finished!" > test_stub_fixmates_markdupes.log
    """

}
