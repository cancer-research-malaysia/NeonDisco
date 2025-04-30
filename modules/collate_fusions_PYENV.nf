// Run POLARS post-processing module
process COLLATE_FUSIONS_PYENV {
    publishDir "${params.outputDir}/${sampleName}/Fusion-Call-Combined-UNFILTERED-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FT-POSTPROCESS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(arFile), path(fcFile)

    output:
        tuple val(sampleName), path("${sampleName}-collated-FT-UNFILTERED.tsv"), emit: collatedFTList
        tuple val(sampleName), path("${sampleName}-collated-FT-UNFILTERED.parquet"), emit: collatedFTParquet

    script:
    """
    echo "Path to Arriba tsv file of raw fusion transcripts: ${arFile}"
    echo "Path to FusionCatcher tsv file of raw fusion transcripts: ${fcFile}"
    if python /home/app/scripts/collate-ft-nf.py ${sampleName} ${arFile} arr ${fcFile} fc; then
        echo "Wrangling completed."
    fi
    """
    stub:
    """
    touch ${sampleName}-collated-FT-UNFILTERED.tsv
    touch ${sampleName}-collated-FT-UNFILTERED.parquet
    echo "stub run finished!\thello my world!" > ${sampleName}-collated-FT-UNFILTERED.tsv
    """
}
