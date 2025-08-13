// Run FT calling module
process CALL_FUSIONS_ARRIBA {
    errorStrategy 'retry'
    maxRetries 3
    maxForks 1
    
    label 'callFusionsAR'
    
    container "${params.container__arriba}"
    
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out/ARRIBA-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }


    input:
        tuple val(sampleName), path(bamFile)
        path arribaDB

    output:
        tuple val(sampleName), path("${sampleName}_arr.tsv"), emit: arriba_fusion_tuple

    script:
    """
    echo "Path to input bam for Arriba: ${bamFile}"
    if arriba-v2--nf.sh ${bamFile} ${sampleName} ${arribaDB} ${params.numCores}; then
        echo "Arriba has finished running on ${sampleName}. Copying main output file..."
        mv ${sampleName}-arriba-fusions.tsv ${sampleName}_arr.tsv
    fi
    """
    stub:
    """
    touch ${sampleName}_arr.tsv
    echo "stub run finished!\thello my world!" > ${sampleName}_arr.tsv
    """
}