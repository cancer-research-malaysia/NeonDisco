// 
process FILTER_FUSIONS_PYENV {
    
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FILTER-FUSIONS -v ${params.metaDataDir}:/home/app/metadata"
    
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
    if wrangle-and-filter-FTs--nf.py ${sampleName} ${collatedFTParquet} /home/app/metadata/${params.panelOfNormalsPq} /home/app/metadata/${params.babiNormalsPq} /home/app/metadata/${params.panelOfCCLEInternalsPq} /home/app/metadata/${params.gaoFusionsPq} /home/app/metadata/${params.mitelmanFusionsPq} /home/app/metadata/${params.klijnFusionsPq} \${OUTPUT_NAME}; then
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
