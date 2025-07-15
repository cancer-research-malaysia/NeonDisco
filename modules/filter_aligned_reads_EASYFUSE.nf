//
process FILTER_ALIGNED_READS_EASYFUSE {
    
    label 'filterReads'
    
    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name ALIGNED-READS-FILTERING"

    input:
        tuple val(sampleName), path(bamFile)

    output:
        tuple val(sampleName), path("*.filtered.bam", arity: '1'), emit: filtered_bam

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}
    BAM=${bamFile}

    echo "Processing files of sample \${SAMPLE_ID}"
    echo "The bam file: \${BAM}"
    
    # Use the read filtering script from EasyFuse program
    if easyfuse-fusionreadfilter--nf.py --input "\${BAM}" --output "\${SAMPLE_ID}.filtered.bam"; then
        echo "EasyFuse read filtering is complete!"
    else
        echo "EasyFuse read filtering failed. Check logs. Exiting..."
        exit 1
    fi

    """
    stub:
    """
    touch ${sampleName}.filtered.bam
    echo "Stub run finished!"
    """
}
