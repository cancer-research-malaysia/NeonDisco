// Run intermediate processing module
process COLLATE_FUSIONS {
    publishDir "${params.output_dir}/${sampleName}", mode: 'copy'
    container "${params.container__pypolars}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name polars-postprocess -v ${params.input_dir}:/work/data -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(arFile), path(fcFile)

    output:
        path "*-fusiontranscript-raw-list.*"

    script:
    """
    echo "Path to Arriba tsv file of raw fusion transcripts: ${arFile}"
    echo "Path to FusionCatcher tsv file of raw fusion transcripts: ${fcFile}"
    if python /work/scripts/collate-ft-nf.py ${sampleName} ${arFile} arr ${fcFile} fc; then
        echo "Wrangling completed."
    fi
    """
    stub:
    """
    mkdir -p stub-o
    touch stub-o/${sampleName}-collated-FT-raw-list.tsv
    """
}
