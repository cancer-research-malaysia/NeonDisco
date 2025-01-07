// Run HLA typing module
process TYPE_HLA_ALLELES_HLAHD {
    publishDir "${params.output_dir}/${sampleName}/HLAHD-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__hlahd}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name hla-typing -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(inputFiles)

    output:
        tuple val(sampleName), path("**/*_final.result.txt"), emit: hla_combined_result

    script:
    """
    # Initialize variables
    input_file1=${inputFiles[0]}  # First file in the list will be our main input
    input_file2=${inputFiles[1]}
    echo "Processing input file: \${input_file1}"
    echo "Number of cores to use: ${params.num_cores}"
    
    echo "Running HLA-HD..."
    # Running HLA-HD on read files
    if bash /work/scripts/hlahd-nf.sh "${sampleName}" "\${read1_file}" "\${read2_file}" ${params.num_cores}; then
        echo "HLA-HD completed its run. Check outputs for run status."
    else
        echo "HLA-HD failed to complete."
        exit 1
    fi
    """
    stub:
    """
    mkdir -p test-stub-dir
    touch test-stub-dir/test_stub_final.result.txt
    echo "stub run finished!" > test-stub-dir/test_stub_final.result.txt
    """
}
