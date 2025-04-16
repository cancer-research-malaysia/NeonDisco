// Trim raw reads
process TRIM_READS_FASTP {
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name TRIM-READS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(readFiles)

    output:
        tuple val(sampleName), path("*_trimmed.R?.fq.gz", arity: '2'), emit: trimmed_reads
    
    script:
    """
    # Initialize variables
    READ1=${readFiles[0]}  # First file in the list will be our main input
    READ2=${readFiles[1]}
    SAMPLE_ID=${sampleName}
    echo "Processing files: \${READ1} & \${READ2} of sample \${SAMPLE_ID}"
    echo "Starting FASTP trimming..."

    # Running FASTP
    if bash /home/app/scripts/fastp-nf.sh "\${READ1}" "\${READ2}" "\${SAMPLE_ID}"; then
        echo "FASTP trimming finished successfully!"
    else
        echo "FASTP trimming failed. Check logs. Exiting..."
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
