#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// All of the default parameters are being set in `nextflow.config`
// Import sub-workflows
include { CALL_FUSION_TRANSCRIPTS_AR } from './modules/call_fusion_transcripts_AR'
include { CALL_FUSION_TRANSCRIPTS_FC } from './modules/call_fusion_transcripts_FC'
include { COLLATE_FUSIONS } from './modules/collate_fusions'
include { TYPE_HLA_ALLELES } from './modules/type_hla_alleles'
include { TESTING_MODULE } from './modules/testing_module'


// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run main.nf -profile <local/awsbatch> <--OPTION NAME> <ARGUMENT>

Required Arguments:
---------------
    -profile                      Either <local> for testing, or <awsbatch> for AWS Batch cluster [MANDATORY]
    --input_dir                   Path to base directory where directories of datasets containing raw fastq or fq files. [MANDATORY]
    --output_dir                  Directory path for output files [MANDATORY]

---------------
    --ftcaller                    Name of the fusion caller to be run. Either <arriba> or <fusioncatcher> or <both>. Defaults to <both> if not specified
    --hla_typing                  Set to <true> if HLA typing subworkflow is required. Defaults to <false> when not specified
    --hla_typing_dir              Directory path to WES or RNA-seq sequencing data to run HLA typing on. Required when --hla_typing is set to <true>
    --hla-only                    Setting to run just the HLA typing subpipeline. Default is set to <false>
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
    }
    log.info "[STATUS] Initializing workflow..."
    log.info "[STATUS] The HLA typing subworkflow is set to <${params.hla_typing}>."
    if ( params.hla_typing && !params.hla_typing_dir ) {
        log.error "The parameter 'hla_typing' is set to <true> (default is <false>) but WES or RNA-seq input directory for HLA typing is not provided. Please provide a valid input directory path."
        log.info "[STATUS] Aborting..."
        exit 1
    } else if ( params.hla_typing && params.hla_typing_dir ) {
        log.info "[STATUS] Checking whether the directory of HLA input files is valid..."
        // Check if input file dir exists or not
        if ( file(params.hla_typing_dir).exists() == true && file(params.hla_typing_dir).isDirectory() ){
            log.info "[STATUS] The input file directory provided is valid."
            if ( params.hla_only ) {
                log.info "[STATUS] HLA typing 'only' argument is set to <true>. Running HLA typing as the only subworkflow to be run..."
                log.info "[STATUS] Getting HLA input files..."
                // Create input file channel
                hla_reads_ch = Channel.fromFilePairs("${params.hla_typing_dir}/*.{fastq,fq,bam}{,.gz}", checkIfExists: true).ifEmpty { 
                    log.error "No valid input WES or RNA-seq read files found in ${params.hla_typing_dir}. Please check your input directory."
                    log.info "[STATUS] Aborting..."
                    exit 1
                }.toSortedList( { a, b -> a[0] <=> b[0] } )  // Sort the channel elements based on the first object of each tuple (sample name) and convert to a channel with a single element which is a list of tuples (NOTE: <=> is an operator for comparison)
    
                .flatMap() // Flatten the single-element channel to a channel with as many elements as there are samples, which is the original structure provided by fromFilePairs
                .view()
            } else {
                log.info "[STATUS] Running HLA typing subworkflow in parallel..."
            }
        } else {
            log.error "The HLA typing input directory is not valid. Please double-check the validity of the provided pathname."
            exit 0
        }
    }
    if ( !params.input_dir || !params.output_dir ) {
        log.error "The input directory path (--input_dir) or the output directory path (--output_dir) is not provided. Please provide a valid input path for all of these parameters."
        log.info "[STATUS] Aborting..."
        exit 1
    }
    log.info "[STATUS] Checking if an alternative fusion caller parameter is set..."
    if ( params.ftcaller != "both" ) {
        if ( !(params.ftcaller in ['arriba', 'fusioncatcher']) ) {
            log.error "Fusion caller is set but the value set <${params.ftcaller}> is unknown. Please specify either <arriba> or <fusioncatcher> or <both> only. <both> is set by default."
            log.info "[STATUS] Aborting..."
            exit 1
        } else {
            log.info "[STATUS] Fusion caller is set to <${params.ftcaller}>."
        }
    } else {
        log.info "[STATUS] Fusion caller parameter is set to <both>."
    }
    log.info "[STATUS] Checking whether the directory of input files is valid..."
    // Check if input file dir exists or not
    if ( file(params.input_dir).exists() == true && file(params.input_dir).isDirectory() ) {
        log.info "[STATUS] The input file directory provided is valid."
        log.info "[STATUS] Fusion caller specified. Getting input files..."
        // Create input file channel
        read_pairs_ch = Channel.fromFilePairs("${params.input_dir}/*{R,r}{1,2}*.{fastq,fq}{,.gz}", checkIfExists: true).ifEmpty { 
                log.error "No valid input FASTQ read files found in ${params.input_dir}. Please check your input directory."
                log.info "[STATUS] Aborting..."
                exit 1
        }.toSortedList( { a, b -> a[0] <=> b[0] } )  // Sort the channel elements based on the first object of each tuple (sample name) and convert to a channel with a single element which is a list of tuples (NOTE: <=> is an operator for comparison)
    
        .flatMap() // Flatten the single-element channel to a channel with as many elements as there are samples, which is the original structure provided by fromFilePairs
        .view()

        count_ch = read_pairs_ch.count() // Get the count into a value channel
        count_ch.view { count -> 
                "Total number of valid input file pairs found: $count" 
            } // View the count
        count_ch.map { count ->
            if (count > 0) {
                log.info "[STATUS] Starting fusion transcript calls..."
            } else {
                log.error "No valid input files found."
                log.info "[STATUS] Aborting..."
                exit 1
            } // Use count in conditional
        }
        /////////// BEGIN WORKFLOW ////////////////

        // Call fusion transcripts
        if (params.ftcaller == 'both' || params.ftcaller == 'arriba') {
            log.info "[STATUS] Running Arriba..."
            arResultTuple = CALL_FUSION_TRANSCRIPTS_AR(read_pairs_ch)
        }
        if (params.ftcaller == 'both' || params.ftcaller == 'fusioncatcher') {
            log.info "[STATUS] Running FusionCatcher..."
            fcResultTuple = CALL_FUSION_TRANSCRIPTS_FC(read_pairs_ch)
        }
        if (params.ftcaller == 'both') {
            log.info "[STATUS] Joining Arriba and FusionCatcher raw output data..."
                
            combinedResultFiles = arResultTuple.arriba_fusion_tuple.join(fcResultTuple.fuscat_fusion_tuple, by: 0).map { sampleName, arFile, fcFile ->
                        tuple(sampleName, arFile, fcFile)
                } // Join result files based on sample name and creates a channel of list of lists with each list containing the sample name at index 0 and then a tuple of 

            combinedResultFiles.view()

            // Wrangle raw tsv to get fusion transcripts called by both Arriba and FusionCatcher
            COLLATE_FUSIONS(combinedResultFiles)
        }
        
    } else {
        log.error "Either the input directory does not exist or it is not a directory. Please double-check the path."
        log.info "[STATUS] Aborting..."
        exit 1
    }     

    // log.info "The path to FusionCatcher DB: ${params.fuscat_db}"
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}





