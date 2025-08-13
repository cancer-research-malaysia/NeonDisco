// Run FT calling module
process CALL_FUSIONS_STARFUSION {
    errorStrategy 'retry'
    maxRetries 3
    maxForks 2
    
    label 'callFusionsSF'
    
    container "${params.container__starfusion}"

    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out/STARFUSION-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(filtFastqs)
        path ctatDB

    output:
        tuple val(sampleName), path("${sampleName}_sf.tsv"), emit: starfus_fusion_tuple

    script:
    """
    echo "Path to input read file 1: ${filtFastqs[0]}"
    echo "Path to input read file 2: ${filtFastqs[1]}"
    if starfusion--nf.sh ${filtFastqs[0]} ${filtFastqs[1]} ${sampleName} ${params.numCores} ${ctatDB}; then
        echo "STARFusion has finished running on ${sampleName}. Copying main output file..."
        cp ${sampleName}-STARFusion-out/star-fusion.fusion_predictions.abridged.tsv ${sampleName}_sf.tsv
    fi
    """
    stub:
    """
    touch ${sampleName}_sf.tsv
    echo "stub run finished!\thello my world!" > ${sampleName}_sf.tsv
    """
}
