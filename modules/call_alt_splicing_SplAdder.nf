// Run HLA typing module
process CALL_ALT_SPLICING_SPLADDER {
    publishDir "${params.output_dir}/${sampleName}/AS/SplAdder-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__spladder}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name call-eis-spladder -v ${params.arriba_db}:/work/libs -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(bam), path(bai)
        val(numCores)

    output:
        tuple val(sampleName), path("SplAdder_${sampleName}")

    script:
    """
    # Initialize variables
    BAM=${bam}  # First file in the list will be our main input
    BAI=${bai}
    echo "Processing bam file \${BAM} & index file \${BAI}"
    echo "Number of cores to use for SplAdder: ${numCores}"

    # create output folder
    mkdir -p /work/nf_work/SplAdder_${sampleName}
    touch \${BAI}

    echo "Starting SplAdder..."
    # Running SPLADDER
    if bash /work/scripts/spladder-nf.sh "\${BAM}" ${numCores} /work/nf_work/SplAdder_${sampleName}; then
        echo "Alternative splicing event calling using SplAdder is complete. Check outputs for run status."
    else
        echo "Alternative splicing event calling failed."
        exit 1
    fi
    """
    stub:
    """
    mkdir -p /work/nf_work/SplAdder_${sampleName}
    echo "stub run finished!" > /work/nf_work/SplAdder_${sampleName}/test.stub
    """
}
