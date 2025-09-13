// 
process COLLATE_FUSIONS_PYENV {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'collateFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(arFile), path(fcFile), path(sfFile)

    output:
        tuple val(sampleName), path("${sampleName}-collated-FT-raw.parquet"), emit: collatedFTParquet
        tuple val(sampleName), path("${sampleName}-collated-FT-raw.tsv"), emit: collatedFTFile

    script:
    """
    echo "Path to Arriba tsv file of raw fusion transcripts: ${arFile}"
    echo "Path to FusionCatcher tsv file of raw fusion transcripts: ${fcFile}"
    echo "Path to StarFusion tsv file of raw fusion transcripts: ${sfFile}"
    echo "Sample name: ${sampleName}"

    OUTPUT_NAME="${sampleName}-collated-FT-raw"
    echo "Generated output filename using sample name: \${OUTPUT_NAME}"

    echo "Running Python script to collate FT files..."
    if collate-FTs--nf.py -s ${sampleName} -o \${OUTPUT_NAME} --arriba ${arFile} --fusioncatcher ${fcFile} --starfusion ${sfFile}; then
        echo "Collation completed successfully."
    else
        echo "Collation failed. Please check the logs for errors."
        exit 1
    fi
    """
    stub:
    """
    touch ${sampleName}-collated-FT-raw.tsv
    touch ${sampleName}-collated-FT-raw.parquet
    echo "stub run finished!\thello my world!" > ${sampleName}-collated-FT-raw.tsv
    """
}
