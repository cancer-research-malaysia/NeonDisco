// Run POLARS post-processing module
process COMBINE_FUSIONS_PYENV {
    publishDir "${params.outputDir}/${sampleName}/Combined-FT-UNFILTERED-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FT-POSTPROCESS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(arFile), path(fcFile)

    output:
        tuple val(sampleName), path("${sampleName}-combined-tool-FT-UNFILTERED.tsv"), emit: combinedFTList
        tuple val(sampleName), path("${sampleName}-combined-tool-FT-UNFILTERED.parquet"), emit: combinedFTParquet

    script:
    """
    echo "Path to Arriba tsv file of raw fusion transcripts: ${arFile}"
    echo "Path to FusionCatcher tsv file of raw fusion transcripts: ${fcFile}"
    if python /home/app/scripts/combine-ft-nf.py ${sampleName} ${arFile} arr ${fcFile} fc; then
        echo "Wrangling completed."
    fi
    """
    stub:
    """
    touch ${sampleName}-combined-tool-FT-UNFILTERED.tsv
    touch ${sampleName}-combined-tool-FT-UNFILTERED.parquet
    echo "stub run finished!\thello my world!" > ${sampleName}-combined-tool-FT-UNFILTERED.tsv
    """
}
