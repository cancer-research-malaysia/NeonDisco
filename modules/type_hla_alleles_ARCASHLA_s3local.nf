// type HLA allotypes using arcasHLA (s3local profile version)
process TYPE_HLA_ALLELES_ARCASHLA_S3LOCAL {
    publishDir "${params.output_dir}/${sampleName}/arcasHLA-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__arcashla}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name arcashla-typing -v \$(pwd):/home/app/nf_work -v ${params.bin_dir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(bambaiPair)

    output:
        tuple val(sampleName), path("*.genotype.json"), emit: allotype_json 
        path("*.genotype.log")
        tuple val(sampleName), path("${params.output_dir}/${sampleName}/arcasHLA-out"), emit: arcasHLA_out
    
    script:
    """
    # Initialize variables
    SAMPLE_ID=${sampleName}
    BAM=${bambaiPair[0]}
    BAM_IDX=${bambaiPair[1]}
    
    echo "Processing files: \${BAM} of sample \${SAMPLE_ID}"
    echo "The index file is: \${BAM_IDX}"
    echo "Number of cores to use: ${params.num_cores}"
    echo "Starting arcasHLA typing..."

    # Running arcasHLA
    if bash /home/app/scripts/arcasHLA-nf.sh "\${SAMPLE_ID}" "\${BAM}" "${params.num_cores * 2}"; then
        echo "arcasHLA typing finished successfully!"
        # remove bam file
        rm -f \${BAM}
        rm -f \${BAM_IDX}
    else
        echo "arcasHLA typing failed. Check logs. Exiting..."
        exit 1
    fi
    """
    stub:
    """
    touch ${sampleName}_test-stub.genotype.json
    touch ${sampleName}_test-stub.genotype.log
    echo "Stub run finished!" > test_stub_arcashla-type.log
    """
}
