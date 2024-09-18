#!/usr/bin/env nextflow

// Run first module
process callFusionTranscripts {
    
    if ( params.ft_caller == 'arriba' ){
        container "${params.container__arriba}"
    } else if ( params.ft_caler == 'fusioncatcher' ){
        container "${params.container__fuscat}"
    }
    
    // Set this for cluster run
    // clusterOptions '-l select=1:ncpus=1:mem=16GB -l walltime=4:00:00 -P 12003580 -q normal'
    // maxForks 40
    // publishDir "${params.output_dir}/sorted_beds/", mode: 'copy'
    
    input:
        path readFile

    script:
    """
    echo "Path to input file: $readFile"
    
    """
}