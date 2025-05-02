// 
process FILTER_FUSIONS_PYENV {
    
    publishDir "${params.outputDir}/${sampleName}/FILTERED-FTs-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name APPLY-FILTERS-FTS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(combinedFTParquet)

    output:
        tuple val(sampleName), path("${sampleName}-combined-tool-FT-FILTERED.tsv"), emit: filteredFusionList
        tuple val(sampleName), path("${sampleName}-combined-tool-FT-FILTERED.parquet"), emit: filteredFusionParquet

    script:
    """
    echo "Path to collated fusion transcripts: ${combinedFTParquet}"
    echo "Running filtering script to filter FT lists based on multiple criteria..."
    if python /home/app/scripts/filter-ft-nf.py ${sampleName} ${combinedFTParquet}; then
        echo "Filtering completed."
    fi

    """
    stub:
    """
    echo "stub run finished!" > ${sampleName}-combined-tool-FT-FILTERED.tsv
    echo "stub run finished!" > ${sampleName}-combined-tool-FT-FILTERED.parquet
    """
}
