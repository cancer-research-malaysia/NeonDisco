#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// All of the default parameters are being set in `nextflow.config`
// Import sub-workflows
include { CALL_FUSION_TRANSCRIPTS_AR } from './modules/call_fusion_transcripts_AR.nf'
include { CALL_FUSION_TRANSCRIPTS_FC } from './modules/call_fusion_transcripts_FC.nf'
include { PREDICT_CODING_SEQ_AGFUSION } from './modules/predict_coding_seq_AGFUSION.nf'
include { COLLATE_FUSIONS_POLARS } from './modules/collate_fusions_POLARS.nf'
include { TYPE_HLA_ALLELES_HLAHD } from './modules/type_hla_alleles_HLAHD.nf'
include { ALIGN_READS_STAR } from './modules/align_reads_STAR.nf'
include { CALL_ALT_SPLICING_SPLADDER } from './modules/call_alt_splicing_SPLADDER.nf'

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

def create_hla_reads_channel(dir_path) {
    // find bam files if any
    bam_files_ch = Channel.fromPath("${dir_path}/*.{bam,bai}", checkIfExists: true)
        .ifEmpty { 
            log.error "No valid input BAM files found in ${dir_path}. Please check your input directory."
            log.info "[STATUS] Aborting..."
            exit 1
        }
        .map { file -> 
            // Get filename without extension
            def name = file.name.replaceAll(/\.(bam|bai)$/, '')
            tuple(name, file)
        }
        .groupTuple()  // Group by sample name
        .map { sample_name, files -> 
            // Sort to ensure BAM comes before BAI
            def sorted_files = files.sort { a, b -> 
                a.name.endsWith("bam") ? -1 : 1  // BAM files come first
            }
            tuple(sample_name, sorted_files)
        }
    
    // find fastq files if any
    // Check for fastq files existence first
    def fastq_files = file("${dir_path}/*{R,r}{1,2}*.{fastq,fq}{,.gz}")
    def has_fastq = fastq_files.size() > 0

    if (has_fastq) {
        log.info "[STATUS] Found FASTQ files. Mixing with BAM files..."
        reads_files_ch = Channel.fromFilePairs("${dir_path}/*{R,r}{1,2}*.{fastq,fq}{,.gz}")
            .toSortedList( { a, b -> a[0] <=> b[0] } )
            .flatMap()
        return bam_files_ch.mix(reads_files_ch)
    } else {
        log.warn "No FASTQ files found. Returning BAM files only..."
        return bam_files_ch
    }
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

    // HLA-only workflow branch
    if (params.hla_only) {
        log.info "[STATUS] HLA typing 'hla_only' argument is set to <true>."
        log.info "[STATUS] Running HLA typing workflow only..."

        if (!params.hla_typing_dir) {
            log.error "HLA typing directory must be specified when running HLA-only mode."
            log.info "[STATUS] Aborting..."
            exit 1
        }
        
        if (!file(params.hla_typing_dir).exists() || !file(params.hla_typing_dir).isDirectory()) {
            log.error "The HLA typing input directory is not valid."
            log.info "[STATUS] Aborting..."
            exit 1
        }

        hla_reads_ch = create_hla_reads_channel(params.hla_typing_dir)
        // hla_reads_ch.view()
        TYPE_HLA_ALLELES_HLAHD(hla_reads_ch, params.num_cores)
    }

    // Full pipeline branch
    else {
        if ( !params.input_dir || !params.output_dir ) {
            log.error "The input directory path (--input_dir) or the output directory path (--output_dir) is not provided. Please provide a valid input path for all of these parameters."
            log.info "[STATUS] Aborting..."
            exit 1
        }

        if ( params.hla_typing ) {
            if ( !params.hla_typing_dir ) {
                log.error "HLA typing directory must be specified when HLA typing workflow is set to <true>."
                log.info "[STATUS] Aborting..."
                exit 1
            }

            if (file(params.hla_typing_dir).exists() && file(params.hla_typing_dir).isDirectory()) {
                hla_reads_ch = create_hla_reads_channel(params.hla_typing_dir)
                // hla_reads_ch.view()
                TYPE_HLA_ALLELES_HLAHD(hla_reads_ch, params.num_cores)
            }
        }

        /////////// BEGIN FUSION TRANSCRIPT CALLING WORKFLOW ////////////////
        log.info "[STATUS] Checking whether the directory of input files is valid..."
        // Check if input file dir exists or not
        if ( file(params.input_dir).exists() == true && file(params.input_dir).isDirectory() ) {
            log.info "[STATUS] The input file directory provided is valid."
            log.info "[STATUS] Checking if an alternative fusion caller parameter is set..."
            
            if ( params.ftcaller != "both" ) {
                if ( !(params.ftcaller in ['arriba', 'fusioncatcher']) ) {
                    log.error "Fusion caller is set but the value set <${params.ftcaller}> is unknown. Please specify either <arriba> or <fusioncatcher> or <both> only. <both> is set by default."
                    log.info "[STATUS] Aborting..."
                    exit 1
                } else {
                    log.info "[STATUS] Fusion caller is set to <${params.ftcaller}>."
                }
            }
            else {
                log.info "[STATUS] Fusion caller parameter is set to <both>."
            }

            log.info "[STATUS] Getting input files..."

            // Collect inputs into channels
            read_pairs_ch = Channel.fromFilePairs("${params.input_dir}/*{R,r}{1,2}*.{fastq,fq}{,.gz}", checkIfExists: true)
            .ifEmpty { 
                log.error "No valid input FASTQ read files found in ${params.input_dir}. Please check your input directory."
                log.info "[STATUS] Aborting..."
                exit 1
            }
            .toSortedList( { a, b -> a[0] <=> b[0] } )
            .flatMap()
            //.view()

            count_ch = read_pairs_ch.count() // Get the count into a value channel
            
            count_ch.view { count -> 
                "[STATUS] Total number of valid input file pairs for fusion transcript calling found: $count" 
            } // View the count
            
            count_ch.map { count ->
                if (count == 0) {
                    log.error "No valid input files found."
                    log.info "[STATUS] Aborting..."
                    exit 1
                }
            }

            // Call fusion transcripts
            if (params.ftcaller == 'arriba') {
                log.info "[STATUS] Running Arriba only..."
                arResultTuple = CALL_FUSION_TRANSCRIPTS_AR(read_pairs_ch, params.num_cores)
                // arResultTuple.view()

                // Create dummy file reference
                def DUMMY_FILE = file("${projectDir}/assets/DUMMY_FILE", checkIfExists: false)

                // Use map to add a dummy file to the tuple
                input_with_dummy_ch = arResultTuple
                    .map { sampleName, ftFile -> 
                        [sampleName, ftFile, DUMMY_FILE] 
                    }
                
                //input_with_dummy_ch.view()
                
                // predict coding sequences
                PREDICT_CODING_SEQ_AGFUSION(input_with_dummy_ch, params.num_cores)
            }

            if (params.ftcaller == 'fusioncatcher') {
                log.info "[STATUS] Running FusionCatcher only..."
                fcResultTuple = CALL_FUSION_TRANSCRIPTS_FC(read_pairs_ch, params.num_cores)
                //fcResultTuple.view()

                // Create dummy file reference
                def DUMMY_FILE = file("${projectDir}/assets/DUMMY_FILE", checkIfExists: false)

                // Use map to add a dummy file to the tuple
                input_with_dummy_ch = fcResultTuple
                    .map { sampleName, ftFile -> 
                        [sampleName, ftFile, DUMMY_FILE] 
                    }
                
                //input_with_dummy_ch.view()
                
                // predict coding sequences
                PREDICT_CODING_SEQ_AGFUSION(input_with_dummy_ch, params.num_cores)

            }

            if (params.ftcaller == 'both') {
                log.info "[STATUS] Running both Arriba and FusionCatcher asynchronously..."
                arResultTuple = CALL_FUSION_TRANSCRIPTS_AR(read_pairs_ch, params.num_cores)
                fcResultTuple = CALL_FUSION_TRANSCRIPTS_FC(read_pairs_ch, params.num_cores)

                combinedResults_ch = arResultTuple.arriba_fusion_tuple
                .join(fcResultTuple.fuscat_fusion_tuple, by: 0)
                .map { sampleName, arFile, fcFile -> tuple(sampleName, arFile, fcFile) }
                .view { sampleName, arFile, fcFile ->
                """
                \n===============================================
Sample: ${sampleName}
-----------------------------------------------
Arriba file    : ${arFile}
FusionCatcher  : ${fcFile}
===============================================
                """
            }

                // predict coding sequences
                PREDICT_CODING_SEQ_AGFUSION(combinedResults_ch, params.num_cores)

                // Wrangle raw TSVs to get fusion transcripts called by both Arriba and FusionCatcher
                COLLATE_FUSIONS_POLARS(combinedResults_ch)
            }

            // first align reads to bam then index
            alignedBams_ch = ALIGN_READS_STAR(read_pairs_ch, params.num_cores)
            // alignedBams_ch.view()

            // Call alt spliced events
            CALL_ALT_SPLICING_SPLADDER(alignedBams_ch, params.num_cores)
        } 
        else {
            log.error "Either the input directory does not exist or it is not a directory. Please double-check the path."
            log.info "[STATUS] Aborting..."
            exit 1
        }
    }
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}





