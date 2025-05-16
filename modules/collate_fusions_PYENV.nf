// 
process COLLATE_FUSIONS_PYENV {
    
    publishDir "${params.outputDir}/${sampleName}/AGGREGATE-FUSION-CALLING-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name COLLATE-FUSIONS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(arFile), path(fcFile)

    output:
        tuple val(sampleName), path("${sampleName}-collated-FT-raw.parquet"), emit: collatedFTParquet

    script:
    """
    echo "Path to Arriba tsv file of raw fusion transcripts: ${arFile}"
    echo "Path to FusionCatcher tsv file of raw fusion transcripts: ${fcFile}"

    OUTPUT_NAME="${sampleName}-collated-FT-raw"
    echo "Generated output filename using sample name: \${OUTPUT_NAME}"

    echo "Running Python script to collate FT files..."
    if python /home/app/scripts/collate-FTs--nf.py -s ${sampleName} -o \${OUTPUT_NAME} --arriba ${arFile} --fusioncatcher ${fcFile}; then
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
