// Run HLA typing module
process ALIGN_READS_STAR {
    publishDir "${params.output_dir}/${sampleName}/EIS/STAR-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__arriba}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name align-reads -v ${params.arriba_db}:/work/libs -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(readFiles)
        val(numCores)

    output:
        tuple val(sampleName), path("*Aligned.sortedByCoord.out.bam"), emit: aligned_reads

    script:
    """
    # Initialize variables
    READ1=${readFiles[0]}  # First file in the list will be our main input
    READ2=${readFiles[1]}
    echo "Processing files: \${READ1} & \${READ2}"
    echo "Number of cores to use: ${numCores}"

    echo "Starting STAR alignment..."
    # Running STAR on read files
    if bash /work/scripts/star-align-nf.sh "\${READ1}" "\${READ2}" ${numCores}; then
        echo "STAR alignment complete. Check outputs for run status."
    else
        echo "STAR alignment failed."
        exit 1
    fi
    """
    stub:
    """
    touch test_stub_Alignment.sortedByCoord.out.bam
    echo "stub run finished!" > test_stub_Alignment.sortedByCoord.out.bam
    """
}
