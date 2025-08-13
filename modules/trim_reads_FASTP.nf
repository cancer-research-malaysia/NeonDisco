// Trim raw reads
process TRIM_READS_FASTP {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'trimReads'
    
    afterScript params.deleteStagedFiles ? "find ./ -name \"${sampleName}*_*.f*q*\" -type l -exec sh -c 'rm -f \$(readlink -f \"{}\")' \\; -delete" : "echo 'Skipping staged file cleanup...'"

    container "${params.container__preproc}"
    
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
    if fastp--nf.sh "\${READ1}" "\${READ2}" "\${SAMPLE_ID}"; then
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
