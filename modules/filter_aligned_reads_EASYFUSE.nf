//
process FILTER_ALIGNED_READS_EASYFUSE {
    cpus params.numCores / 2
    
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

    # Validate BAM integrity
    if ! samtools quickcheck "\${QUERYNAMEBAM}"; then
        echo "ERROR: BAM file \${QUERYNAMEBAM} failed integrity check" >&2
        exit 1
    fi
    echo "BAM validation successful."

    # Check for pairing issues
    echo "Checking BAM pairing statistics..."
    samtools flagstat "\${QUERYNAMEBAM}" > flagstat.txt
    cat flagstat.txt

    PROPERLY_PAIRED=\$(grep "properly paired" flagstat.txt | awk '{print \$1}')
    TOTAL_READS=\$(grep "in total" flagstat.txt | awk '{print \$1}')

    if [ "\${PROPERLY_PAIRED}" -eq 0 ]; then
        echo "WARNING: No properly paired reads found in BAM file" >&2
        echo "This may cause issues with EasyFuse. Consider checking alignment parameters." >&2
    fi

    echo "Proceeding with EasyFuse read filtering..."
    if easyfuse-fusionreadfilter--nf.py --input "\${QUERYNAMEBAM}" --output "\${SAMPLE_ID}.filtered.bam"; then
        echo "EasyFuse read filtering is complete!"
    else
        echo "ERROR: EasyFuse read filtering failed for sample \${SAMPLE_ID}" >&2
        echo "This may indicate pairing or BAM format issues. Check the logs above." >&2
        exit 1
    fi
    """
    stub:
    """
    touch ${sampleName}.filtered.bam
    echo "Stub run finished!"
    """
}
