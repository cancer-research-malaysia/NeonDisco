// 
process FILTER_FUSIONS_PYENV {
    
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FILTER-FUSIONS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(collatedFTParquet)

    output:
        tuple val(sampleName), path("${sampleName}-collated-FT-normFiltered.tsv"), emit: filteredFusions

    script:
    """
    echo "Path to collated fusion transcripts: ${collatedFTParquet}"
    echo "Sample name: ${sampleName}"
    echo "Generating output file name: ${sampleName}-collated-FT-normFiltered"
    OUTPUT_NAME="${sampleName}-collated-FT-normFiltered"
    echo "Running filtering script to filter for tumor-specific FTs..."
    if python /home/app/scripts/filter-FTs--nf.py ${sampleName} ${combinedFTParquet} ${params.panelOfNormalsTsv} \${OUTPUT_NAME}; then
        echo "Filtering completed."
    fi

    """
    stub:
    """
    echo "stub run finished!" > ${sampleName}-collated-FT-normFiltered.tsv
    """
}
