//
process FILTER_ALIGNED_READS_EASYFUSE {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'filterReads'
    
    container "${params.container__pyenv}"

    input:
        tuple val(sampleName), path(bamFile), path(bamIdx), path(querynameBamFile)

    output:
        tuple val(sampleName), path("*.filtered.bam", arity: '1'), emit: filtered_bam

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}
    QUERYNAMEBAM=${querynameBamFile}

    echo "Processing files of sample \${SAMPLE_ID}"
    echo "The bam file: \${QUERYNAMEBAM}"
    
    # Use the read filtering script from EasyFuse program
    if easyfuse-fusionreadfilter--nf.py --input "\${QUERYNAMEBAM}" --output "\${SAMPLE_ID}.filtered.bam"; then
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
