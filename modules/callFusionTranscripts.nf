#!/usr/bin/env nextflow

// Run first module
process callFusionTranscripts {
    
    if ( params.ft_caller == 'arriba' ){
        container "${params.container__arriba}"
    } else if ( params.ft_caller == 'fusioncatcher' ){
        container "${params.container__fuscat}"
    }
    
    // Set this for cluster run
    // clusterOptions '-l select=1:ncpus=1:mem=16GB -l walltime=4:00:00 -P 12003580 -q normal'
    // maxForks 40
    // publishDir "${params.output_dir}/sorted_beds/", mode: 'copy'
    
    input:
        tuple val(sampleName), path(readFiles)

    script:
    """
    echo "Path to input read file 1: ${readFiles[0]}"
    echo "Path to input read file 2: ${readFiles[1]}"
    docker --version
    docker run 
    """
}