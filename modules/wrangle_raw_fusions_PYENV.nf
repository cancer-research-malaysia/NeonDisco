// 
process WRANGLE_RAW_FUSIONS_PYENV {
    cpus 4
    
    label 'wrangleRawFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(collatedFTParquet)

    output:
        tuple val(sampleName), path("${sampleName}-wrangled-unfiltered-fusions.tsv"), emit: wrangledUnfilteredFusionsTsv

    script:
    """
    echo "Path to collated fusion transcripts: ${collatedFTParquet}"
    echo "Sample name: ${sampleName}"
    echo "Generating output file name: ${sampleName}-wrangled-unfiltered-fusions"
    OUTPUT_NAME="${sampleName}-wrangled-unfiltered-fusions"
    export POLARS_MAX_THREADS=${task.cpus}
    echo "Setting POLARS_MAX_THREADS to \${POLARS_MAX_THREADS}..."

    echo "Running wrangling script to process raw FTs..."
    if wrangle-FTs-only--nf.py ${sampleName} ${collatedFTParquet} \${OUTPUT_NAME}; then
        echo "Wrangling completed."
    fi
    """
    stub:
    """
    echo "stub run finished!" > ${sampleName}-wrangled-unfiltered-fusions.tsv
    echo "Stub run for WRANGLE_RAW_FUSIONS completed successfully."
    """
}
