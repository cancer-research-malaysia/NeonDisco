// Fix mate information
process FIXMATES_MARKDUPES_SAMTOOLS_S3LOCAL {
    maxForks 4
    // publishDir "${params.output_dir}/${sampleName}/SAMTOOLS-postproc-out", mode: 'copy',
    //     saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
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

        # clean up original bam
        rm -f ${bam}
        rm -f ${sampleName}_fixmates.bam

        # transfer to s3
        aws s3 cp ${sampleName}_fixmates_markdupes.bam s3://${params.s3_bucket}/main-outputs/HLA-typing-arcasHLA/MyBrCa-RNA-seq/${sampleName}/
        aws s3 cp ${sampleName}_fixmates_markdupes.bam.bai s3://${params.s3_bucket}/main-outputs/HLA-typing-arcasHLA/MyBrCa-RNA-seq/${sampleName}/
        if [ \$? -eq 0 ]; then
            echo "Transfer to S3 complete!"
        else
            echo "Transfer to S3 failed. Check logs. Exiting..."
            exit 1
        fi

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
