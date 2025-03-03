// Fix mate information
process FIXMATES_MARKDUPES_BAMS_SAMTOOLS {
    errorStrategy 'finish'
    maxForks 2
    publishDir "${params.output_dir}/${sampleName}/SAMTOOLS-postproc-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name samtools-postproc-bams -v \$(pwd):/home/app/nf_work -v ${params.bin_dir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(bam)
    
    output:
        tuple val(sampleName), path("*_fixmates_markdupes.ba*", arity: '2'), emit: final_bams
    
    script:
    """
    samtools sort -n -@ ${params.num_cores} -m 4G -O bam ${bam} | \
    samtools fixmate -pcmu -O bam - ${sampleName}_fixmates.bam

    if [ -f ${sampleName}_fixmates.bam ]; then
        echo "Fixmate information for ${sampleName} is complete!"
        # mark duplicates
        samtools sort -@ ${params.num_cores} -m 4G -O bam ${sampleName}_fixmates.bam | \
        samtools markdup -@ ${params.num_cores} - ${sampleName}_fixmates_markdupes.bam

        # then index
        samtools index ${sampleName}_fixmates_markdupes.bam && echo "Indexing complete!"
    else
        echo "Fixmate information for ${sampleName} failed. Check logs. Exiting..."
        exit 1
    fi
    """
    stub:
    """
    touch ${sampleName}_trimmed.R1.fq.gz
    touch ${sampleName}_trimmed.R2.fq.gz
    echo "Stub run finished!" > test_stub_fastp-trim.log
    """

}
