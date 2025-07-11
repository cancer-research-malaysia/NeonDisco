// Run FT calling module
process CALL_FUSIONS_ARRIBA {
    maxForks 1
    label 'callFusionsArriba'
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out/ARRIBA-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__arriba}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FT-CALLING-ARRIBA -v ${params.arribaDB}:/home/app/arriba-db -v ${params.starIndex}:/home/app/starIdx -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(bamFile)

    output:
        tuple val(sampleName), path("${sampleName}_arr.tsv"), emit: arriba_fusion_tuple

    script:
    """
    echo "Path to input bam for Arriba: ${bamFile}"
    if bash /home/app/scripts/arriba-v2--nf.sh ${bamFile} ${sampleName} "/home/app/arriba-db" ${params.numCores}; then
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