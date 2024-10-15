#!/usr/bin/env nextflow

// Run first module
process collateFusionTranscripts {
    //publishDir "${params.output_dir}/${sampleName}", mode: 'copy'
    container "${params.container__pypolars}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name polars-postprocess -v ${params.input_dir}:/work/data -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        tuple path(ar_file), path(fc_file)

    //output:
        //

    script:
    """
    echo "Path to input read file 1: ${readFiles[0]}"
    echo "Path to input read file 2: ${readFiles[1]}"
    bash /work/scripts/fuscat-nf.sh ${readFiles[0]} ${readFiles[1]}
    """
}
