process TESTING_MODULE {
    tag "${sample_id}"
    
    input:
        tuple val(sample_id), path(input_file)
    
    output:
        tuple val(sample_id), path("${sample_id}_output.txt"), emit: test_output
    
    script:
    """
    # This just creates an empty output file
    touch ${sample_id}_output.txt
    
    # Optional: Add some dummy content
    echo "Processing sample: ${sample_id}" > ${sample_id}_output.txt
    """
}
