
process AGGREGATE_CUSTOM_PONS_PYENV {
    
    publishDir "${params.outputDir}", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name AGGREGATE-PONS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    

    input:
    path(ponFiles)
    
    output:
    path("${params.ponsOutputName}.parquet"), emit: finalCustomPONs
    
    script:
    """
    echo "Path to collated PON files: ${ponFiles}"
    echo "Output file name: ${params.ponsOutputName}.parquet"
    echo "Running Python script to aggregate custom PON files..."

    if python /home/app/scripts/concatenate-parquet-pons--nf.py ${params.ponsOutputName}.parquet ${ponFiles}; then
        echo "Aggregation completed successfully."
    else
        echo "Error during aggregation. Please check the input files and script."
        exit 1
    fi

    """


}