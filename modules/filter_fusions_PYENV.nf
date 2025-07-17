// 
process FILTER_FUSIONS_PYENV {
    
    label 'filterFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(collatedFTParquet)

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
    if wrangle-and-filter-FTs--nf.py ${sampleName} ${collatedFTParquet} /tmp/metadata/${params.panelOfNormalsPq} /tmp/metadata/${params.babiNormalsPq} /tmp/metadata/${params.panelOfCCLEInternalsPq} /tmp/metadata/${params.gaoFusionsPq} /tmp/metadata/${params.mitelmanFusionsPq} /tmp/metadata/${params.klijnFusionsPq} \${OUTPUT_NAME}; then
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
