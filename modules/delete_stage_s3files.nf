// align trimmed reads
process DELETE_STAGE_S3FILES {
    maxForks 4
    afterScript "find ./ -name '${sampleName}*.f*.gz' -type l -exec sh -c 'rm -f \$(readlink -f \"{}\")' \\; -delete"
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name star-align-reads -v ${params.arriba_db}:/home/app/libs -v \$(pwd):/home/app/nf_work -v ${params.bin_dir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(readFile1), path(readFile2)

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}
    CORES=${params.num_cores}

    echo "Processing files of sample \${SAMPLE_ID}"
    echo "Number of cores to use: ${params.num_cores}"

    # test
    if touch S3_${sampleName}.txt; then
        echo "${readFile1} and ${readFile2} are staged at \$(realpath ${readFile1}) and \$(realpath ${readFile2})" > S3_${sampleName}.txt
        cp ${readFile1} ./read1-test.fq.gz
        cp ${readFile2} ./read2-test.fq.gz
        echo "Test files copied to current directory"
    else
        echo "Error. Exiting..."
        exit 1
    fi

    """
}
