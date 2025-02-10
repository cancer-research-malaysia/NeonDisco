#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import submodules
include { TRIM_READS_FASTP } from './modules/trim_reads_FASTP.nf'
include { TYPE_HLA_ALLELES_HLAHD } from './modules/type_hla_alleles_HLAHD.nf'

///////////// DEFINITIONS ////////////
// Function to print help message
def helpMessage() {
    log.info"""
Usage:

nextflow run main.nf -profile <local/awsbatch> <--OPTION NAME> <ARGUMENT>

Required Arguments:
---------------
    -profile            Either <local> for testing, or <awsbatch> for AWS Batch cluster [MANDATORY]
    --input_dir         Path to base directory where directories of datasets containing raw fastq or fq files. [MANDATORY]
    --output_dir        Directory path for output files [MANDATORY]

---------------
    --ftcaller          Name of the fusion caller to be run. Either <arriba> or <fusioncatcher> or <both>. Defaults to <both> if not specified
    --hla_typing        Set to <true> if HLA typing subworkflow is required. Defaults to <false> when not specified
    --hla_typing_dir    Directory path to WES or RNA-seq sequencing data to run HLA typing on. Required when --hla_typing is set to <true>
    --hla-only          Setting to run just the HLA typing subpipeline. Default is set to <false>
    --trimming          Set to <true> to perform read trimming on input files (only works with FASTQ inputs)
    --help              Print this help message and exit
    
    """.stripIndent()
}

// Function to create input channel for HLA typing
def createHlaChannel(dir_path) {
	// Check for existence of different file types
    def bam_files = file("${dir_path}/*.bam")
    def fastq_files = file("${dir_path}/*{R,r}{1,2}*.{fastq,fq}{,.gz}")
    
    // Initialize empty channels
    def bam_ch = Channel.empty()
    def fastq_ch = Channel.empty()

    // Process BAM files if they exist
    if (bam_files) {
        log.info "[STATUS] Found BAM files in ${dir_path}"
        bam_ch = Channel.fromPath("${dir_path}/*.{bam,bai}")
            .map { file -> 
                def name = file.name.replaceAll(/\.(bam|bai)$/, '')
                tuple(name, file)
            }
            .groupTuple()
            .map { sample_name, files -> 
                // Verify we have both BAM and BAI
                if (!files.any { it.name.endsWith('.bai') }) {
                    log.error "Missing index file (.bai) for BAM: ${sample_name}"
                }
                if (!files.any { it.name.endsWith('.bam') }) {
                    log.error "Missing BAM file for index: ${sample_name}"
                }
                
                // Sort to ensure BAM comes before BAI
                def sorted_files = files.sort { a, _b -> 
                    a.name.endsWith("bam") ? -1 : 1
                }
                tuple(sample_name, sorted_files)
            }
    }
    
    // Process FASTQ files if they exist
    if (fastq_files) {
        log.info "[STATUS] Found FASTQ files in ${dir_path}"
        fastq_ch = Channel.fromFilePairs("${dir_path}/*{R,r}{1,2}*.{fastq,fq}{,.gz}", checkIfExists: true)
            .toSortedList( { a, b -> a[0] <=> b[0] } )
            .flatMap()
    }
    
    // Check if we found any files at all
    if (!bam_files && !fastq_files) {
        log.error "No BAM or FASTQ files found in ${dir_path}"
		exit 1
    }
    
    // Mix the channels and return
    return bam_ch.mix(fastq_ch)
}

// Function to validate input directory
def validateInputDir(dir_path, dir_type) {
    if (!file(dir_path).exists() || !file(dir_path).isDirectory()) {
        log.error "The ${dir_type} directory '${dir_path}' is not valid or does not exist."
		exit 1
    }
}

workflow TRIM_READS {
    take:
        reads_ch
    main:
        TRIM_READS_FASTP(reads_ch)
    emit:
        trimmed_reads = TRIM_READS_FASTP.out.trimmed_reads
}

workflow TYPE_HLAS {
    take:
        input_ch
    main:
        TYPE_HLA_ALLELES_HLAHD(input_ch)
    emit:
        hla_types = TYPE_HLA_ALLELES_HLAHD.out.hla_combined_result
}

workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

	// Check output directory value
	if (!params.output_dir) {
		log.error "Output directory path or prefix must be specified."
		exit 1
	}

    // Validate HLA typing directory when in HLA-only mode
    if (params.hla_only) {
        if (!params.hla_typing_dir) {
            log.error "HLA typing directory must be specified when running in HLA-only mode."
			exit 1
        }
        validateInputDir(params.hla_typing_dir, "HLA TYPING")

        // Create input channel for all files
        input_ch = createHlaChannel(params.hla_typing_dir)

		// Debug view of initial channel
        log.info "[DEBUG] Initial mixed channel content:"
        input_ch.view({ id, files ->
            def fileType = files instanceof List ? "FASTQ" : "BAM"
            return "[Sample: $id, Type: $fileType, Files: $files]"
        })

        // Process FASTQ files if trimming is requested
        if (params.trimming) {
            input_ch.branch {
                fastq: (it[1][0].name.endsWith('.fastq.gz') || 
				it[1][0].name.endsWith('.fq.gz') || 
				it[1][0].name.endsWith('.fastq') || 
				it[1][0].name.endsWith('.fq'))

    			bam: (it[1][0].name.endsWith('.bam'))
            }
			.set { branched }

			// Debug views of branched channels
            log.info "[DEBUG] Branched FASTQ channel content:"
            branched.fastq.view({ id, files ->
                "[FASTQ Sample: $id, Files: $files]"
            })
			log.info "[DEBUG] Branched BAM channel content:"
            branched.bam.view({ id, files ->
                "[BAM Sample: $id, Files: $files]"
            })
			// Process FASTQ files through trimming
            TRIM_READS(branched.fastq)
            
            // Combine trimmed FASTQ with BAM files
            processed_input = TRIM_READS.out.trimmed_reads.mix(branched.bam)

			// Debug view of final mixed channel
            log.info "[DEBUG] Final mixed channel after trimming:"
            processed_input.view({ id, files ->
                def fileType = files instanceof List ? "Trimmed FASTQ" : "BAM"
                return "[Sample: $id, Type: $fileType, Files: $files]"
            })

            //// Process FASTQ files through trimming
            //TRIM_READS(branched_input.fastq)
            
            //// Combine trimmed FASTQ with BAM files
            //processed_input = TRIM_READS.out.trimmed_reads.mix(branched_input.bam)
        } else {
            processed_input = input_ch
        }

        // Run HLA typing
        TYPE_HLAS(processed_input)
    }

    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }
}
