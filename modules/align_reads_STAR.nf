// align trimmed reads
process ALIGN_READS_STAR {
    publishDir "${params.output_dir}/${sampleName}/STAR-alignment", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name star-align-reads -v ${params.arriba_db}:/home/app/libs -v \$(pwd):/home/app/nf_work -v ${params.bin_dir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(readFiles)

    output:
        tuple val(sampleName), path("*Aligned.sortedByCoord.out.bam"), path("*Aligned.sortedByCoord.out.bam.bai"), emit: aligned_reads

    script:
    """
    # Initialize variables
    READ1=${readFiles[0]}  # First file in the list will be our main input
    READ2=${readFiles[1]}
    echo "Processing files: \${READ1} & \${READ2}"
    echo "Number of cores to use: ${params.num_cores}"

    echo "Starting STAR alignment..."
    # Running STAR on read files
    if bash /work/scripts/star-align-nf.sh "\${READ1}" "\${READ2}" ${params.num_cores}; then
        echo "STAR alignment complete. Check outputs for run status."
    else
        echo "STAR alignment failed."
        exit 1
    fi
    """
    stub:
    """
    touch test_stub_Aligned.sortedByCoord.out.bam
    touch test_stub_Aligned.sortedByCoord.out.bam.bai
    echo "stub run finished!" > test_stub_Aligned.sortedByCoord.out.bam
    """
}
