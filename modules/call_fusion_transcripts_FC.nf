// Run calling module
process CALL_FUSION_TRANSCRIPTS_FC {
    // maxForks 10
    publishDir "${params.output_dir}/${sampleName}", mode: 'copy'
    container "${params.container__fuscat}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name fuscat-ftcall -v ${params.fuscat_db}:/work/libs -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(readFiles)

    output:
        tuple val(sampleName), path("${sampleName}_fc.tsv"), emit: fuscat_fusion_tuple

    script:
    """
    echo "Path to input read file 1: ${readFiles[0]}"
    echo "Path to input read file 2: ${readFiles[1]}"
    if bash /work/scripts/fuscat-nf.sh ${readFiles[0]} ${readFiles[1]} 8; then
        echo "FusionCatcher has finished running on ${sampleName}. Copying main output file..."
        cp final-list_candidate-fusion-genes.txt ${sampleName}_fc.tsv
    fi
    """
}
