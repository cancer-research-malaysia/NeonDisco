// generic module to upload files to S3
process UPLOAD_S3FILES {
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name star-align-reads -v ${params.arriba_db}:/home/app/libs -v \$(pwd):/home/app/nf_work -v ${params.bin_dir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(arcasHlaOutPath)

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}

    echo "Uploading arcasHLA output of \${SAMPLE_ID}..."

    # Uploading to S3
    if aws s3 cp ${arcasHlaOutPath}/ "s3://${params.s3_bucket}/main-outputs/HLA-typing-arcasHLA/MyBrCa-RNA-seq/\${SAMPLE_ID}/arcasHLA-out/" --recursive --region ap-southeast-5; then
        echo "Upload finished successfully for \${SAMPLE_ID}!"
    else
        echo "Upload failed. Check logs. Exiting..."
        exit 1
    fi
    """
}
