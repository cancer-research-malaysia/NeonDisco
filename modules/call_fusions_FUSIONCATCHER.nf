// Run FT calling module
process CALL_FUSIONS_FUSIONCATCHER {
    maxForks 2
    
    publishDir "${params.outputDir}/${sampleName}/FUSIONCATCHER-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__fuscat}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FT-CALLING-FUSCAT -v ${params.fuscatDB}:/home/app/libs -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(filtFastqs)

    output:
        tuple val(sampleName), path("${sampleName}_fc.tsv"), emit: fuscat_fusion_tuple

    script:
    """
    echo "Path to input read file 1: ${filtFastqs[0]}"
    echo "Path to input read file 2: ${filtFastqs[1]}"
    if bash /home/app/scripts/fuscat-nf.sh ${filtFastqs[0]} ${filtFastqs[1]} ${params.numCores}; then
        echo "FusionCatcher has finished running on ${sampleName}. Copying main output file..."
        cp final-list_candidate-fusion-genes.txt ${sampleName}_fc.tsv
    fi
    """
    stub:
    """
    touch ${sampleName}_fc.tsv
    echo "stub run finished!\thello my world!" > ${sampleName}_fc.tsv
    """
}
