#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import submodules
include { TRIM_READS_FASTP } from './modules/trim_reads_FASTP.nf'
include { ALIGN_READS_1PASS_STARSAM } from './modules/align_reads_1pass_STARSAM.nf'
include { ALIGN_READS_2PASS_STARSAM } from './modules/align_reads_2pass_STARSAM.nf'
include { FIXMATES_MARKDUPES_SAMTOOLS } from './modules/fixmates_markdupes_SAMTOOLS.nf'
include { FISH_HLA_READS_SAMPBOWT } from './modules/fish_hla_reads_SAMPBOWT.nf'
include { TYPE_HLA_ALLELES_HLAHD } from './modules/type_hla_alleles_HLAHD.nf'
include { TYPE_HLA_ALLELES_ARCASHLA } from './modules/type_hla_alleles_arcasHLA.nf'
include { CALL_FUSIONS_ARRIBA } from './modules/call_fusions_ARRIBA.nf'
include { CALL_FUSIONS_FUSIONCATCHER } from './modules/call_fusions_FUSIONCATCHER.nf'
include { COLLATE_FUSIONS_PYENV } from './modules/collate_fusions_PYENV.nf'

// Function to print help message
def helpMessage() {
    log.info"""
Usage:

nextflow run main.nf -profile <local/awsbatch> <--OPTION NAME> <ARGUMENT>

Required Arguments:
---------------
    -profile            Either <local> for testing, or <awsbatch> for AWS Batch cluster [MANDATORY]
    --input_dir         Path to directory containing RNA-seq FASTQ/BAM files [MANDATORY]
    --output_dir        Directory path for output files [MANDATORY]

Optional Arguments:
---------------
    --trimming          Set to <true> to perform read trimming on input files (only works with FASTQ inputs)
    --hla_only         Set to <true> to run only HLA typing workflow. Default: false
    --help             Print this help message and exit
    """.stripIndent()
}

// Function to validate input directory
def validateInputDir(dir_path) {
    if (!file(dir_path).exists() || !file(dir_path).isDirectory()) {
        log.error "The input directory '${dir_path}' does not exist or is not valid."
        return false
    }
    return true
}

// Function to create input channel
def createInputChannel(dir_path) {
    // Check for existence of different file types
    def bam_files = file("${dir_path}/*.bam")
    def fastq_files = file("${dir_path}/*{R,r}{1,2}*.{fastq,fq}{,.gz}")
    
    // Initialize empty channels
    def bam_Ch = Channel.empty()
    def fastq_Ch = Channel.empty()

    // Process BAM files if they exist
    if (bam_files) {
        log.info "[STATUS] Found BAM files in ${dir_path}"
        bam_Ch = Channel.fromPath("${dir_path}/*.{bam,bai}")
            .map { file -> 
                def name = file.name.replaceAll(/\.(bam|bai)$/, '')
                tuple(name, file)
            }
            .groupTuple()
            .map { sample_name, files -> 
                def sorted_files = files.sort { a, _b -> 
                    a.name.endsWith("bam") ? -1 : 1
                }
                tuple(sample_name, sorted_files)
            }
    }
    
    // Process FASTQ files if they exist
    if (fastq_files) {
        log.info "[STATUS] Found FASTQ files in ${dir_path}"
        fastq_Ch = Channel.fromFilePairs("${dir_path}/*{R,r}{1,2}*.{fastq,fq}{,.gz}", checkIfExists: true)
            .toSortedList( { a, b -> a[0] <=> b[0] } )
            .flatMap()
    }
    
    // Check if we found any files
    if (!bam_files && !fastq_files) {
        log.error "No BAM or FASTQ files found in ${dir_path}"
        exit 1
    }
    
    return bam_Ch.mix(fastq_Ch)
}

// Subworkflow definitions
workflow TRIM_READS {
    take:
        reads_Ch
    main:
        TRIM_READS_FASTP(reads_Ch)
    emit:
        trimmed_reads = TRIM_READS_FASTP.out.trimmed_reads
}

workflow ALIGN_READS_1PASS {
    take:
        reads_Ch
    main:
        FIXMATES_MARKDUPES_SAMTOOLS(ALIGN_READS_1PASS_STARSAM(reads_Ch).aligned_bams)
    emit:
        aligned_Bams = FIXMATES_MARKDUPES_SAMTOOLS.out.final_bams
}

workflow ALIGN_READS_2PASS {
    take:
        reads_Ch
    main:
        FIXMATES_MARKDUPES_SAMTOOLS(ALIGN_READS_2PASS_STARSAM(reads_Ch).aligned_bams)
    emit:
        aligned_Bams = FIXMATES_MARKDUPES_SAMTOOLS.out.final_bams
}

workflow HLA_TYPING_HLAHD {
    take:
        input_Ch
    main:
        // First fish for HLA reads
        FISH_HLA_READS_SAMPBOWT(input_Ch)
        
        // Then type the HLA alleles using the fished reads
        TYPE_HLA_ALLELES_HLAHD(FISH_HLA_READS_SAMPBOWT.out.fished_files)
    emit:
        hla_types = TYPE_HLA_ALLELES_HLAHD.out.hla_combined_result
}

workflow HLA_TYPING_ARCASHLA {
    take:
        input_Ch
    main:
        // Then type the HLA alleles
        TYPE_HLA_ALLELES_ARCASHLA(input_Ch)
    emit:
        hla_types = TYPE_HLA_ALLELES_ARCASHLA.out.allotype_json
}

// Main workflow
workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate required parameters
    if (!params.input_dir || !params.output_dir) {
        log.error "Input and output directories must be specified."
        exit 1
    }

    // Validate input directory
    if (!validateInputDir(params.input_dir)) {
        exit 1
    }

    // Create input channel
    def input_Ch = createInputChannel(params.input_dir)

    // Branch input channel based on file type (just peek at the 1st element of the second element [file list] of the input tuple)
    def branched = input_Ch.branch {
        fastq: it[1][0].name =~ /\.(fastq|fq)(\.gz)?$/
        bam: it[1][0].name =~ /\.bam$/
    }

    // Process FASTQ files if present and trimming is requested
    def procInput_Ch
    def trimmedFastqs = params.trimming ? TRIM_READS(branched.fastq).trimmed_reads : branched.fastq
    procInput_Ch = trimmedFastqs.mix(branched.bam)
    procInput_Ch.view()

    // Execute workflows based on hla_only parameter
    if (params.hla_only) {
        // Run only HLA typing using HLAHD
        // HLA_TYPING_HLAHD(procInput_Ch)

        // Run only HLA typing from fq files using arcasHLA
        aligned_Ch = ALIGN_READS_2PASS(procInput_Ch)
        HLA_TYPING_ARCASHLA(aligned_Ch)

    } else {
        // main pipeline
        aligned_Ch = ALIGN_READS_2PASS(procInput_Ch)
        HLA_TYPING_ARCASHLA(aligned_Ch)

        // gene fusion identification submodule
        CALL_FUSIONS_ARRIBA(procInput_Ch)
        CALL_FUSIONS_FUSIONCATCHER(procInput_Ch)

        // Join the outputs based on sample name
        CALL_FUSIONS_ARRIBA.out.arriba_fusion_tuple
            .join(CALL_FUSIONS_FUSIONCATCHER.out.fuscat_fusion_tuple)
            .set { combinedFcFiles_Ch }
    
        // Run the collation process with the joined output
        COLLATE_FUSIONS_PYENV(combinedFcFiles_Ch)
    }

	// Completion handler
	workflow.onComplete = {
    	println "Pipeline completed at: $workflow.complete"
    	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
	}
}
