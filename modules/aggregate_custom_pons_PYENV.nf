
process AGGREGATE_CUSTOM_PONS_PYENV {
    errorStrategy 'retry'
    maxRetries 3

    label 'aggregateCustomPons'
    
    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
        

    input:
    path(ponFiles)
    
    output:
    path("${params.ponsOutputName}.parquet"), emit: finalCustomPONs
    path("${params.ponsOutputName}.tsv")
    
    script:
    """
    echo "Path to collated PON files: ${ponFiles}"
    echo "Output file name: ${params.ponsOutputName}.parquet"
    echo "Running Python script to aggregate custom PON files..."

    if concatenate-parquet-pons--nf.py ${params.ponsOutputName}.parquet ${ponFiles}; then
        echo "Aggregation completed successfully."
    else
        echo "Error during aggregation. Please check the input files and script."
        exit 1
    fi

    """


}