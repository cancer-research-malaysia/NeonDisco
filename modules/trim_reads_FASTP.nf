// Trim raw reads
process TRIM_READS_FASTP {
    publishDir "${params.output_dir}/${sampleName}/FASTP-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name fastp-trimming -v \$(pwd):/home/app/nf_work -v ${params.bin_dir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(readFiles)

    script:
    """
    # Initialize variables
    READ1=${readFiles[0]}  # First file in the list will be our main input
    READ2=${readFiles[1]}
    echo "Processing files: \${READ1} & \${READ2}"

    echo "Starting FASTP trimming..."
    # Running FASTP
    #if fastp -i "\${READ1}" -I "\${READ2}" -o trimmed.R1.fq.gz -O trimmed.R2.fq.gz -p; then
    #    echo "FASTP trimming finished successfully."
    #else
    #    echo "FASTP trimming failed. Exiting..."
    #    exit 1
    #fi
    """
    stub:
    """
    touch trimmed.R1.fq.gz
    touch trimmed.R2.fq.gz
    echo "Stub run finished!" > test_stub.log
    """
}
