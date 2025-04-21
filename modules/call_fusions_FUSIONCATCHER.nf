// Run FT calling module
process CALL_FUSIONS_FUSIONCATCHER {
    maxForks 2
    publishDir "${params.outputDir}/${sampleName}/FUSIONCATCHER-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__fuscat}"
    
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FT-CALLING-FUSCAT -v ${params.fuscatDB}:/work/libs -v \$(pwd):/work/nf_work -v ${params.binDir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(trimmedReads)

    output:
        tuple val(sampleName), path("${sampleName}_fc.tsv"), emit: fuscat_fusion_tuple

    script:
    """
    echo "Path to input read file 1: ${trimmedReads[0]}"
    echo "Path to input read file 2: ${trimmedReads[1]}"
    if bash /work/scripts/fuscat-nf.sh ${trimmedReads[0]} ${trimmedReads[1]} ${params.numCores}; then
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
