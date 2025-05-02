// 
process FILTER_FUSIONS_PYENV {
    
    publishDir "${params.outputDir}/${sampleName}/FILTERED-FTs-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name APPLY-FILTERS-FTS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(collatedFTParquet)

    // output:
    //     tuple val(sampleName), path("${sampleName}_agfusion_*"), emit: agfusion_outdir

    script:
    """
    echo "Path to collated fusion transcripts: ${collatedFTParquet}"
    echo "Running filtering script to filter FT lists based on multiple criteria..."
    if python /home/app/scripts/filter-ft-nf.py ${sampleName} ${collatedFTParquet}; then
        echo "Filtering completed."
    fi

    """
    stub:
    """
    echo "stub run finished!" > ${sampleName}
    """
}
