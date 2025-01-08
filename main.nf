#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import submodules
include { PREPROC_HLA_TYPING_INPUT_SAMPICARD } from './modules/preproc_hla_typing_input_SAMPICARD.nf'
include { TYPE_HLA_ALLELES_HLAHD } from './modules/type_hla_alleles_HLAHD.nf'
include { CALL_FUSION_TRANSCRIPTS_AR } from './modules/call_fusion_transcripts_AR.nf'
include { CALL_FUSION_TRANSCRIPTS_FC } from './modules/call_fusion_transcripts_FC.nf'
include { PREDICT_CODING_SEQ_AGFUSION } from './modules/predict_coding_seq_AGFUSION.nf'
include { COLLATE_FUSIONS_POLARS } from './modules/collate_fusions_POLARS.nf'
include { ALIGN_READS_STAR } from './modules/align_reads_STAR.nf'
include { CALL_ALT_SPLICING_SPLADDER } from './modules/call_alt_splicing_SPLADDER.nf'

// Function definitions
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
    def bam_files_ch = Channel.fromPath("${dir_path}/*.{bam,bai}", checkIfExists: true)
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
            def sorted_files = files.sort { a, _b -> 
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
        def reads_files_ch = Channel.fromFilePairs("${dir_path}/*{R,r}{1,2}*.{fastq,fq}{,.gz}")
            .toSortedList( { a, b -> a[0] <=> b[0] } )
            .flatMap()
        return bam_files_ch.mix(reads_files_ch)
    } else {
        log.warn "No FASTQ files found. Returning BAM files only..."
        return bam_files_ch
    }
}
///////////////////////////////////////////////////////////////////////////

// Read Trimming and Alignment Workflow
workflow TRIM_AND_ALIGN_READS {
    take:
        input_dir
    
    main:
        read_pairs_ch = Channel.fromFilePairs("${input_dir}/*{R,r}{1,2}*.{fastq,fq}{,.gz}", checkIfExists: true)
            .ifEmpty { error "No valid input FASTQ read files found in ${input_dir}" }
            .toSortedList( { a, b -> a[0] <=> b[0] } )
            .flatMap()

        alignedBams_ch = ALIGN_READS_STAR(read_pairs_ch)

    emit:
        reads = read_pairs_ch
        bams = alignedBams_ch
}


// HLA Typing Workflow
workflow HLA_TYPING {
    take:
        hla_typing_dir
    
    main:
        if (!file(hla_typing_dir).exists() || !file(hla_typing_dir).isDirectory()) {
            log.error "The HLA typing input directory is not valid."
            exit 1
        }

        hla_reads_ch = create_hla_reads_channel(hla_typing_dir)
        preprocessed_ch = PREPROC_HLA_TYPING_INPUT_SAMPICARD(hla_reads_ch)
        preprocessed_ch.view()
        TYPE_HLA_ALLELES_HLAHD(preprocessed_ch)
}

// Fusion Analysis Workflow
workflow FUSION_CALLING {
    take:
        read_pairs_ch
    
    main:

        // Handle different fusion caller scenarios
        if (params.ftcaller == 'arriba') {
            arResultTuple = CALL_FUSION_TRANSCRIPTS_AR(read_pairs_ch)
            def DUMMY_FILE = file("${projectDir}/assets/DUMMY_FILE", checkIfExists: false)
            input_with_dummy_ch = arResultTuple.map { sampleName, ftFile -> 
                [sampleName, ftFile, DUMMY_FILE] 
            }
            PREDICT_CODING_SEQ_AGFUSION(input_with_dummy_ch)
        }
        else if (params.ftcaller == 'fusioncatcher') {
            fcResultTuple = CALL_FUSION_TRANSCRIPTS_FC(read_pairs_ch)
            def DUMMY_FILE = file("${projectDir}/assets/DUMMY_FILE", checkIfExists: false)
            input_with_dummy_ch = fcResultTuple.map { sampleName, ftFile -> 
                [sampleName, ftFile, DUMMY_FILE] 
            }
            PREDICT_CODING_SEQ_AGFUSION(input_with_dummy_ch)
        }
        else if (params.ftcaller == 'both') {
            arResultTuple = CALL_FUSION_TRANSCRIPTS_AR(read_pairs_ch)
            fcResultTuple = CALL_FUSION_TRANSCRIPTS_FC(read_pairs_ch)
            
            combinedResults_ch = arResultTuple.arriba_fusion_tuple
                .join(fcResultTuple.fuscat_fusion_tuple, by: 0)
                .map { sampleName, arFile, fcFile -> tuple(sampleName, arFile, fcFile) }
            
            PREDICT_CODING_SEQ_AGFUSION(combinedResults_ch)
            COLLATE_FUSIONS_POLARS(combinedResults_ch)
        }
}

// Main workflow
workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 1
    }

    // Input validation
    if (params.hla_only && !params.hla_typing_dir) {
        error "HLA typing directory must be specified when running in HLA-only mode."
    }

    if (!params.hla_only && !params.input_dir) {
        error "The input directory path (--input_dir) is not provided."
    }

    // Execute workflows based on parameters
    if (params.hla_only) {
        log.info "[STATUS] Running HLA typing workflow only..."
        HLA_TYPING(params.hla_typing_dir)
    }
    else {
        if (params.hla_typing) {
            log.info "[STATUS] Running HLA typing workflow..."
            HLA_TYPING(params.hla_typing_dir)
        }
        
        //log.info "[STATUS] Running fusion analysis workflow..."
        //FUSION_ANALYSIS(params.input_dir, params.ftcaller)
    }

    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }
}


