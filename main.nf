#!/usr/bin/env nextflow

// All of the default parameters are being set in `nextflow.config`
// Import sub-workflows
include { callFusionTranscriptsAR } from './modules/callFusionTranscriptsAR'
include { callFusionTranscriptsFC } from './modules/callFusionTranscriptsFC'
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
    --input_dir                   Path to base directory where directories of datasets containing BAM files [MANDATORY]
    --ftcaller                    Name of the fusion caller to be run. Either <arriba> or <fusioncatcher> [MANDATORY]

  Output path:
    --output_dir                  Directory path for output files [MANDATORY]

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
    } else if ( !params.input_dir || !params.ftcaller ){
        log.error "The input_dir or the fusion caller parameter is not set. Please provide a valid input directory path or fusion caller parameter."
        exit 1
    } else {
        log.info "Initializing workflow..."
        log.info "Checking whether the directory of input files is valid..."

        // Check if input file dir exists or not
        if ( file(params.input_dir).exists() == true && file(params.input_dir).isDirectory() == true ) {
            log.info "The input file directory provided is valid."

            // check if fusion caller is valid; ftcaller should either be 'arriba' or 'fusioncatcher' or 'both'
            if (!(params.ftcaller in ['arriba', 'fusioncatcher', 'both'])) {
                log.info "The fusion caller specified does not exist. Please specify either <arriba> or <fusioncatcher> or <both> only. <both> is set as default."
                exit 1
            } else {
                log.info "Fusion caller specified. Getting input files..."
                
                // Create input file channel
                read_pairs_ch = Channel
                    .fromFilePairs("${params.input_dir}/*{R,r}{1,2}*.{fastq,fq}{,.gz}", checkIfExists: true)
                    .ifEmpty { 
                            log.error "No valid input FASTQ read files found in ${params.input_dir}. Please check your input directory."
                            exit 1
                        }
                    // Sort the channel elements based on the first object of each tuple,
                    // that is, the sample name, and convert to a channel with a single
                    // element which is a list of tuples
                    .toSortedList( { a, b -> a[0] <=> b[0] } ) // <=> is an operator for comparison
                    // flatten the single-element channel to a channel with as many elements
                    // as there are samples, which is the original structure provided by
                    // fromFilePairs
                    .flatMap()
                    // View the channel elements by printing it to the screen
                    .view()
                
                // Count and display the number of file pairs
                read_pairs_ch.count().view { count -> 
                    "Total number of valid input file pairs found: $count" 
                    }
                // View 
                // read_pairs_ch.view { sample_id, files -> 
                //     "Sample: $sample_id, READ 1: ${files[0].name}, READ 2: ${files[1].name}"
                //     }

                // view channel raw
                //read_pairs_ch.view()

                /////////// BEGIN WORKFLOW ////////////////

                log.info "Starting fusion transcript calls..."
                
                if ( params.ftcaller == 'both' ){
                    // call fusion transcripts
                    log.info "<both> is specified so running Arriba and Fusioncatcher callers..."
                    callFusionTranscriptsAR(read_pairs_ch)
                    callFusionTranscriptsFC(read_pairs_ch)
                } else if ( params.ftcaller == 'arriba' ){
                    log.info "Running Arriba only..."
                    callFusionTranscriptsAR(read_pairs_ch)
                } else if ( params.ftcaller == 'fusioncatcher' ){
                    log.info "Running Fusioncatcher only..."
                    callFusionTranscriptsFC(read_pairs_ch)
                }
            }
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