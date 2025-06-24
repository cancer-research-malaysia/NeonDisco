// 
process COLLATE_CUSTOM_PONS_PYENV {
    
    publishDir "${params.outputDir}/${sampleName}/CUSTOM-FUSION-PONS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name COLLATE-FUSION-PONS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(arFile), path(fcFile), path(sfFile)

    output:
        path("${sampleName}-collated-Normals-FTs.parquet"), emit: collatedFTParquet
        path("${sampleName}-collated-Normals-FTs.tsv")

    script:
    """
    echo "Path to Arriba tsv file of raw fusion transcripts: ${arFile}"
    echo "Path to FusionCatcher tsv file of raw fusion transcripts: ${fcFile}"
    echo "Path to StarFusion tsv file of raw fusion transcripts: ${sfFile}"
    echo "Sample name: ${sampleName}"

    OUTPUT_NAME="${sampleName}-collated-Normals-FTs"
    echo "Generated output filename using sample name: \${OUTPUT_NAME}"

    echo "Running Python script to collate FT files..."
    if python /home/app/scripts/collate-FTs--nf.py -s ${sampleName} -o \${OUTPUT_NAME} --arriba ${arFile} --fusioncatcher ${fcFile} --starfusion ${sfFile}; then
        echo "Collation completed successfully."
    else
        echo "Collation failed. Please check the logs for errors."
        exit 1
    fi
    """
    stub:
    """
    touch ${sampleName}-collated-Normals-FTs.tsv
    touch ${sampleName}-collated-Normals-FTs.parquet
    echo "stub run finished!\thello my world!" > ${sampleName}-collated-Normals-FTs.tsv
    """
}
