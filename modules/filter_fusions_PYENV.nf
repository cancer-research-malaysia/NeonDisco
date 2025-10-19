process FILTER_FUSIONS_PYENV {
    cpus 4
    
    label 'filterFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(collatedFTParquet)
        path metaDataDir

    output:
        tuple val(sampleName), path("${sampleName}-collated-FT-normFiltered.tsv"), emit: filteredFusions
        tuple val(sampleName), path("${sampleName}-collated-FT-normFiltered-unique-genePairs-for-FusInspector.txt"), emit: uniqueFiltFusionPairsForFusIns

    script:
    """
    echo "Path to collated fusion transcripts: ${collatedFTParquet}"
    echo "Sample name: ${sampleName}"
    echo "Generating output file name: ${sampleName}-collated-FT-normFiltered"
    OUTPUT_NAME="${sampleName}-collated-FT-normFiltered"

    echo "Running filtering script to filter for tumor-specific FTs..."
    if wrangle-and-filter-FTs--nf.py ${sampleName} \
    ${collatedFTParquet} \
    ${metaDataDir}/${params.panelOfNormalsPq} \
    ${metaDataDir}/${params.babiNormalsPq} \
    ${metaDataDir}/${params.panelOfCCLEInternalsPq} \
    ${metaDataDir}/${params.gaoFusionsPq} \
    ${metaDataDir}/${params.mitelmanFusionsPq} \
    ${metaDataDir}/${params.klijnFusionsPq} \
    \${OUTPUT_NAME}; then
        echo "Filtering completed."
    fi
    """
    stub:
    """
    echo "stub run finished!" > ${sampleName}-collated-FT-normFiltered.tsv
    echo "stub run finished!" > ${sampleName}-collated-FT-normFiltered-unique-genePairs-for-FusInspector.txt
    echo "Stub run for FILTER_FUSIONS_PYENV completed successfully."
    """
}
