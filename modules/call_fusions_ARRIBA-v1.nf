// Run FT calling module
process CALL_FUSIONS_ARRIBA_V1 {
    maxForks 1
    
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out/ARRIBA-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__arriba}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FT-CALLING-ARRIBA -v ${params.arribaDB}:/home/app/arriba-db -v ${params.starIndex}:/home/app/starIdx -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(trimmedFastqs)

    output:
        tuple val(sampleName), path("${sampleName}_arr.tsv"), emit: arriba_fusion_tuple

    script:
    """
    if bash /home/app/scripts/arriba-v1-nf.sh ${trimmedFastqs[0]} ${trimmedFastqs[1]} ${params.numCores} "/home/app/starIdx" "/home/app/arriba-db"; then
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