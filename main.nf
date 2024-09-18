#!/usr/bin/env nextflow

// All of the default parameters are being set in `nextflow.config`
// Import sub-workflows
include { callFusionTranscripts } from './modules/callFusionTranscripts'
//include { generateBAMPaths } from './modules/generateBAMPaths'
//include { callVariants } from './modules/callVariants'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run main.nf <--OPTION NAME> <ARGUMENT>

Required Arguments:
---------------
  Input paths:
    --input_dir                     Path to base directory where directories of datasets containing BAM files [MANDATORY]
    --genome_fasta                Path to the reference genome fasta file [MANDATORY]
    --fpscore_matrix              Globbed path to base directory (path/to/dir/*.ext) where the merged footprint score matrix files of all studied motifs are located (saved as parquet files) [MANDATORY]

  Input manifest:
    --dataset_id_list             Path to a list of dataset IDs to work on. Required if [--run_mode] is set to `subset`

  Output path:
    --output_dir                  Directory path for output VCF files [MANDATORY]

---------------
    --help                        Print this help message and exit
    
    """.stripIndent()
}

// Main workflow
workflow {
    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    } else if ( !params.input_dir ) {
        log.error "The input_dir parameter is not set. Please provide a valid input directory path."
        exit 1
    } else {
        log.info "Initializing workflow..."
        log.info "Checking whether the directory of input files is valid..."

        // Check if input file dir exists or not
        if ( file(params.input_dir).exists() == true && file(params.input_dir).isDirectory() == true ) {
            log.info "The input file directory provided is valid."

            // Create input file channel
            read_pairs_ch = Channel
                    .fromFilePairs("${params.input_dir}/*{R,r}{1,2}*.{fastq,fq}{,.gz}", checkIfExists: true)
                    .ifEmpty { 
                            log.error "No valid input FASTQ read files found in ${params.input_dir}. Please check your input directory."
                            exit 1
                        }

            // Count and display the number of file pairs
            read_pairs_ch.count().view { count -> 
                    "Total number of valid input file pairs found: $count" 
                    }
            
            // View 
            read_pairs_ch.view { sample_id, files -> 
                "Sample: $sample_id, READ 1: ${files[0].name}, READ 2: ${files[1].name}"
                }
            
            // view channel raw
            //read_pairs_ch.view()

            /////////// BEGIN WORKFLOW ////////////////

            // call fusion transcripts

            callFusionTranscripts(read_pairs_ch)

        }
        else {
            log.error "The input file directory does not exist. Please provide a valid directory path."
            exit 1
        }
    }
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}