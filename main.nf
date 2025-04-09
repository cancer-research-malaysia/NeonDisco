#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import submodules
include { TRIM_READS_FASTP } from './modules/trim_reads_FASTP.nf'
include { ALIGN_READS_1PASS_STARSAM } from './modules/align_reads_1pass_STARSAM.nf'
include { ALIGN_READS_2PASS_STARSAM } from './modules/align_reads_2pass_STARSAM.nf'
include { FIXMATES_MARKDUPES_SAMTOOLS } from './modules/fixmates_markdupes_SAMTOOLS.nf'
include { FISH_HLA_READS_SAMPBOWT } from './modules/fish_hla_reads_SAMPBOWT.nf'
include { TYPE_HLA_ALLELES_HLAHD } from './modules/type_hla_alleles_HLAHD.nf'
include { TYPE_HLA_ALLELES_ARCASHLA } from './modules/type_hla_alleles_ARCASHLA.nf'
include { CALL_FUSIONS_ARRIBA } from './modules/call_fusions_ARRIBA.nf'
include { CALL_FUSIONS_FUSIONCATCHER } from './modules/call_fusions_FUSIONCATCHER.nf'
include { COLLATE_FUSIONS_PYENV } from './modules/collate_fusions_PYENV.nf'
include { PREDICT_CODING_SEQ_AGFUSION } from './modules/predict_coding_seq_AGFUSION.nf'

// s3local testing submodules
include { ALIGN_READS_2PASS_STARSAM_S3LOCAL } from './modules/align_reads_2pass_STARSAM_s3local.nf'
include { FIXMATES_MARKDUPES_SAMTOOLS_S3LOCAL } from './modules/fixmates_markdupes_SAMTOOLS_s3local.nf'
include { TYPE_HLA_ALLELES_ARCASHLA_S3LOCAL } from './modules/type_hla_alleles_ARCASHLA_s3local.nf'

// Function to print help message
def helpMessage() {
    log.info"""
Usage:

nextflow run main.nf -profile <local/awsbatch/s3local> <--OPTION NAME> <ARGUMENT>

Required Arguments:
---------------
    -profile            Either <s3local> for testing, or <awsbatch> for Batch [MANDATORY]
    --manifestPath          Path to tab-delimited manifest file [MANDATORY]
                         – must contain sample ID and read1 and read2 local filepaths or remote s3 filepaths

Optional Arguments:
---------------
    --inputDir          Path to local directory containing BAM/FASTQ input files
    --outputDir         Local directory path for output
    --trimming          Set to <true> to perform read trimming on FASTQ input [DEFAULT: false]
    --HLATypingOnly     Set to <true> to exclusively run HLA typing workflow [DEFAULT: false]
    --help              Print this help message and exit
    """.stripIndent()
}

// Function to validate input directory
def validateInputDir(dirPath) {
    if (!file(dirPath).exists() || !file(dirPath).isDirectory()) {
        log.error "The local input directory '${dirPath}' does not exist or is not valid."
        return false
    }
    return true
}

// Function to create input channel if local directory is provided
def createInputChannelFromPOSIX(dirPath) {
    // Check for existence of different file types
    def bamFiles = file("${dirPath}/*.bam")
    def fastqFiles = file("${dirPath}/*{R,r}{1,2}*.{fastq,fq}{,.gz}")
    
    // Initialize empty channels
    def bamCh = Channel.empty()
    def fastqCh = Channel.empty()

    // Process BAM files if they exist
    if (bamFiles) {
        log.info "[STATUS] Found BAM files in ${dirPath}"
        bamCh = Channel.fromPath("${dirPath}/*.{bam,bai}")
            .map { file -> 
                def name = file.name.replaceAll(/\.(bam|bai)$/, '')
                tuple(name, file)
            }
            .groupTuple()
            .map { sample, files -> 
                def sortedFiles = files.sort { a, _b -> 
                    a.name.endsWith("bam") ? -1 : 1
                }
                tuple(sample, sortedFiles)
            }
    }
    
    // Process FASTQ files if they exist
    if (fastqFiles) {
        log.info "[STATUS] Found FASTQ files in ${dirPath}"
        fastqCh = Channel.fromFilePairs("${dirPath}/*{R,r}{1,2}*.{fastq,fq}{,.gz}", checkIfExists: true)
            .toSortedList( { a, b -> a[0] <=> b[0] } )
            .flatMap()
    }
    
    // Check if we found any files
    if (!bamFiles && !fastqFiles) {
        log.error "No BAM or FASTQ files found in ${dirPath}"
        exit 1
    }
    
    return bamCh.mix(fastqCh)
}

// Function to read manifest file and create input channels
def createInputChannelFromManifest(manifestPath) {
    // Validate that manifest file exists
    def manifestFile = file(manifestPath)
    if (!manifestFile.exists()) {
        log.error "Manifest file not found at path: ${manifestPath}"
        exit 1
    }

    // Create channel from manifest file
    def inputCh = Channel
        .fromPath(manifestPath)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            // Assuming the TSV has columns: sampleName, read1_S3Path, read2_S3Path
            def sampleName = row.sampleName
            def read1 = row.read1_S3Path
            def read2 = row.read2_S3Path
            
            // Validate required fields
            if (!sampleName || !read1 || !read2) {
                log.error "Invalid manifest row: ${row}. Ensure all fields are present."
                return null
            }
            
            // Return tuple of sample name and read file paths
            return tuple(sampleName, read1, read2)
        }
        .filter { it != null }

    return inputCh
}


// Subworkflow definitions
workflow TRIMMING {
    take:
        readsCh
    main:
        TRIM_READS_FASTP(readsCh)
    emit:
        trimmedCh = TRIM_READS_FASTP.out.trimmed_reads
}

workflow ALIGNMENT_1P {
    take:
        trimmedCh
    main:
        FIXMATES_MARKDUPES_SAMTOOLS(ALIGN_READS_1PASS_STARSAM(trimmedCh).aligned_bams)
    emit:
        alignedBamCh = FIXMATES_MARKDUPES_SAMTOOLS.out.final_bams
}

workflow ALIGNMENT_2P {
    take:
        trimmedCh
    main:
        FIXMATES_MARKDUPES_SAMTOOLS(ALIGN_READS_2PASS_STARSAM(trimmedCh).aligned_bams)
    emit:
        alignedBamCh = FIXMATES_MARKDUPES_SAMTOOLS.out.final_bams
}

workflow ALIGNMENT_2P_S3LOCAL {
    take:
        trimmedCh
    main:
        FIXMATES_MARKDUPES_SAMTOOLS_S3LOCAL(ALIGN_READS_2PASS_STARSAM_S3LOCAL(trimmedCh).aligned_bams)
    emit:
        alignedBamCh = FIXMATES_MARKDUPES_SAMTOOLS_S3LOCAL.out.final_bams
}

workflow HLA_TYPING_HLAHD {
    take:
        inputCh
    main:
        // First fish for HLA reads
        FISH_HLA_READS_SAMPBOWT(inputCh)
        
        // Then type the HLA alleles using the fished reads
        TYPE_HLA_ALLELES_HLAHD(FISH_HLA_READS_SAMPBOWT.out.fished_files)
    emit:
        hlaTypes = TYPE_HLA_ALLELES_HLAHD.out.hla_combined_result
}

workflow HLA_TYPING_ARCASHLA {
    take:
        inputCh
    main:
        // Then type the HLA alleles
        TYPE_HLA_ALLELES_ARCASHLA(inputCh)
    emit:
        hlaTypes = TYPE_HLA_ALLELES_ARCASHLA.out.allotype_json
}

workflow HLA_TYPING_ARCASHLA_S3LOCAL {
    take:
        inputCh
    main:
        // Then type the HLA alleles
        TYPE_HLA_ALLELES_ARCASHLA_S3LOCAL(inputCh)
    emit:
        hlaTypes = TYPE_HLA_ALLELES_ARCASHLA_S3LOCAL.out.allotype_json
}

// Main workflow
// workflow {
//     // Show help message if requested
//     if (params.help) {
//         helpMessage()
//         exit 0
//     }

//     // Validate required parameters
//     if (!params.input_dir || !params.output_dir) {
//         log.error "Input and output directories must be specified."
//         exit 1
//     }

//     // Validate input directory
//     if (!validateInputDir(params.input_dir)) {
//         exit 1
//     }

//     // Create input channel
//     def input_Ch = createInputChannel(params.input_dir)

//     // Branch input channel based on file type (just peek at the 1st element of the second element [file list] of the input tuple)
//     def branched = input_Ch.branch {
//         fastq: it[1][0].name =~ /\.(fastq|fq)(\.gz)?$/
//         bam: it[1][0].name =~ /\.bam$/
//     }

//     // Process FASTQ files if present and trimming is requested
//     def procInput_Ch
//     def trimmedFastqs = params.trimming ? TRIM_READS(branched.fastq).trimmed_reads : branched.fastq
//     procInput_Ch = trimmedFastqs.mix(branched.bam)
//     procInput_Ch.view()

//     // Execute workflows based on hla_only parameter
//     if (params.hla_only) {
//         // Run only HLA typing using HLAHD
//         HLA_TYPING_HLAHD(procInput_Ch)

//         // Run only HLA typing from fq files using arcasHLA
//         aligned_Ch = ALIGN_READS_2PASS(procInput_Ch)
//         HLA_TYPING_ARCASHLA(aligned_Ch)

//     } else {
//         // main pipeline
//         aligned_Ch = ALIGN_READS_2PASS(procInput_Ch)
//         HLA_TYPING_ARCASHLA(aligned_Ch)

//         // gene fusion identification submodule
//         CALL_FUSIONS_ARRIBA(procInput_Ch)
//         CALL_FUSIONS_FUSIONCATCHER(procInput_Ch)

//         // Join the outputs based on sample name
//         CALL_FUSIONS_ARRIBA.out.arriba_fusion_tuple
//             .join(CALL_FUSIONS_FUSIONCATCHER.out.fuscat_fusion_tuple)
//             .set { combinedFTFiles_Ch }
        
//         combinedFTFiles_Ch.view()
    
//         // Run the collation process with the joined output
//         COLLATE_FUSIONS_PYENV(combinedFTFiles_Ch)

//         // Run AGFusion to predict fusion protein sequences
//         PREDICT_CODING_SEQ_AGFUSION(COLLATE_FUSIONS_PYENV.out.collatedFTList)
//     }

// 	// Completion handler
// 	workflow.onComplete = {
//     	println "Pipeline completed at: $workflow.complete"
//     	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
// 	}
// }

// temp workflow for HLA retyping
workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Create input channel from manifest file
    // Add a parameter for manifest file path in your nextflow.config or pass via CLI
    def inputCh = createInputChannelFromManifest(params.manifestPath)

    inputCh.view()

    // Execute workflows based on hla_only parameter
    if (params.hlaOnly) {

        // Run only HLA typing from fq files using arcasHLA
        alignedCh = ALIGNMENT_2P_S3LOCAL(inputCh)
        alignedCh.view()
        HLA_TYPING_ARCASHLA_S3LOCAL(alignedCh)
        
    } else {
        // pass
    }

	// Completion handler
	workflow.onComplete = {
    	println "Pipeline completed at: $workflow.complete"
    	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        //workDir.resolve("stage-${workflow.sessionId}").deleteDir()
	}
}