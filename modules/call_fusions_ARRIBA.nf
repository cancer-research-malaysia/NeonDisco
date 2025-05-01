// Run FT calling module
process CALL_FUSIONS_ARRIBA {
    maxForks 1
    
    publishDir "${params.outputDir}/${sampleName}/ARRIBA-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__arriba}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FT-CALLING-ARRIBA -v ${params.arribaDB}:/home/app/arriba-db -v ${params.starIndex}:/home/app/starIdx -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(trimmedReads)

    output:
        tuple val(sampleName), path("${sampleName}_arr.tsv"), emit: arriba_fusion_tuple

    script:
    """
    echo "Path to input read file 1: ${trimmedReads[0]}"
    echo "Path to input read file 2: ${trimmedReads[1]}"
    if bash /home/app/scripts/arriba-nf.sh ${trimmedReads[0]} ${trimmedReads[1]} ${params.numCores} "/home/app/starIdx" "/home/app/arriba-db"; then
        echo "Arriba has finished running on ${sampleName}. Copying main output file..."
        cp ${sampleName}-arriba-fusions.tsv ${sampleName}_arr.tsv
    fi
    """
    stub:
    """
    touch ${sampleName}_arr.tsv
    echo "stub run finished!\thello my world!" > ${sampleName}_arr.tsv
    """
}