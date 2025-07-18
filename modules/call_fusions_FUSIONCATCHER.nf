// Run FT calling module
process CALL_FUSIONS_FUSIONCATCHER {
    maxForks 4
    label 'callFusionsFC'
    
    container "${params.container__fuscat}"

    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out/FUSIONCATCHER-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    
    input:
        tuple val(sampleName), path(filtFastqs)

    output:
        tuple val(sampleName), path("${sampleName}_fc.tsv"), emit: fuscat_fusion_tuple

    script:
    """
    echo "Path to input read file 1: ${filtFastqs[0]}"
    echo "Path to input read file 2: ${filtFastqs[1]}"
    if fuscat--nf.sh ${filtFastqs[0]} ${filtFastqs[1]} /tmp/fuscat-db ${params.numCores} ${sampleName}; then
        echo "FusionCatcher has finished running on ${sampleName}. Copying main output file..."
        cp ${sampleName}/final-list_candidate-fusion-genes.txt ${sampleName}_fc.tsv
    fi
    """
    stub:
    """
    touch ${sampleName}_fc.tsv
    echo "stub run finished!\thello my world!" > ${sampleName}_fc.tsv
    """
}
