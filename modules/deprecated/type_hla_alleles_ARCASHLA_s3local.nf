// 
process TYPE_HLA_ALLELES_ARCASHLA_S3LOCAL {
    afterScript "find ./ -name ${bambaiPair[0]} -type l -exec sh -c 'rm -f \$(readlink -f \"{}\")' \\; -delete"
    publishDir "${params.outputDir}/${sampleName}/arcasHLA-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__arcashla}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name HLATYPING -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(bambaiPair)

    output:
        tuple val(sampleName), path("*.genotype.json"), emit: allotype_json 
        path("*.genotype.log")
    
    script:
    """
    # Initialize variables
    SAMPLE_ID=${sampleName}
    BAM=${bambaiPair[0]}
    BAM_IDX=${bambaiPair[1]}
    
    echo "Processing files: \${BAM} of sample \${SAMPLE_ID}"
    echo "The index file is: \${BAM_IDX}"
    echo "Number of cores to use: ${params.numCores}"
    echo "Starting arcasHLA typing..."

    # Running arcasHLA
    if bash /home/app/scripts/arcasHLA-nf.sh "\${SAMPLE_ID}" "\${BAM}" "${params.numCores * 2}"; then
        echo "arcasHLA typing finished successfully!"
        # remove fq files
        rm -f *.f*q.*
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
