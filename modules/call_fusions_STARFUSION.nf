// Run FT calling module
process CALL_FUSIONS_STARFUSION {
    maxForks 1
    
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out/STARFUSION-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__starfusion}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FT-CALLING-STARFUSION -v ${params.ctatDB}:/home/app/libs -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(filtFastqs)

    output:
        tuple val(sampleName), path("${sampleName}_sf.tsv"), emit: starfus_fusion_tuple

    script:
    """
    echo "Path to input read file 1: ${filtFastqs[0]}"
    echo "Path to input read file 2: ${filtFastqs[1]}"
    if bash /home/app/scripts/starfusion--nf.sh ${filtFastqs[0]} ${filtFastqs[1]} ${sampleName} ${params.numCores} /home/app/libs; then
        echo "STARFusion has finished running on ${sampleName}. Copying main output file..."
        #cp final-list_candidate-fusion-genes.txt ${sampleName}_sf.tsv
        touch ${sampleName}_sf.tsv
    fi
    """
    stub:
    """
    touch ${sampleName}_sf.tsv
    echo "stub run finished!\thello my world!" > ${sampleName}_sf.tsv
    """
}
