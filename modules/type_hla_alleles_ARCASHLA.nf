// type HLA allotypes using arcasHLA
process TYPE_HLA_ALLELES_ARCASHLA {
    publishDir "${params.outputDir}/${sampleName}/HLA-TYPING-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    container "${params.container__arcashla}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name TYPE-HLAS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(bam), path(bamIdx)

    output:
        tuple val(sampleName), path("*.genotype.json"), emit: allotype_json
    
    script:
    """
    # Initialize variables
    SAMPLE_ID=${sampleName}
    BAM=${bam}
    BAM_IDX=${bamIdx}
    
    echo "Processing files: \${BAM} of sample \${SAMPLE_ID}"
    echo "The index file is: \${BAM_IDX}"
    echo "Number of cores to use: ${params.numCores}"
    echo "Starting arcasHLA typing..."

    # Running arcasHLA
    if bash /home/app/scripts/arcasHLA--nf.sh "\${SAMPLE_ID}" "\${BAM}" "${params.numCores}"; then
        echo "arcasHLA typing finished successfully!"
    else
        echo "arcasHLA typing failed. Check logs. Exiting..."
        exit 1
    fi
    """
    stub:
    """
    touch ${sampleName}_test-stub.genotype.json

    echo "Stub run finished!" > test_stub_type-hlas.log
    """
}
