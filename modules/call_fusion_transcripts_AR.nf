// Run FT calling module
process CALL_FUSION_TRANSCRIPTS_AR {
    // maxForks 10
    publishDir "${params.output_dir}/${sampleName}/FT/ARRIBA-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__arriba}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name arriba-ftcall -v ${params.arriba_db}:/work/libs -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(readFiles)
        val(numCores)

    output:
        tuple val(sampleName), path("${sampleName}_arr.tsv"), emit: arriba_fusion_tuple

    script:
    """
    echo "Path to input read file 1: ${readFiles[0]}"
    echo "Path to input read file 2: ${readFiles[1]}"
    if bash /work/scripts/arriba-nf.sh ${readFiles[0]} ${readFiles[1]} ${numCores}; then
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