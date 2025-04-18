#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import submodules
include { TRIM_READS_FASTP } from './modules/trim_reads_FASTP.nf'
include { ALIGN_READS_TWOPASS_STARSAM } from './modules/align_reads_twopass_STARSAM.nf'
include { FIXMATES_MARKDUPES_SAMTOOLS } from './modules/fixmates_markdupes_SAMTOOLS.nf'
include { FISH_HLA_READS_SAMPBOWT } from './modules/fish_hla_reads_SAMPBOWT.nf'
include { TYPE_HLA_ALLELES_HLAHD } from './modules/type_hla_alleles_HLAHD.nf'
include { TYPE_HLA_ALLELES_ARCASHLA } from './modules/type_hla_alleles_ARCASHLA.nf'
include { CALL_FUSIONS_ARRIBA } from './modules/call_fusions_ARRIBA.nf'
include { CALL_FUSIONS_FUSIONCATCHER } from './modules/call_fusions_FUSIONCATCHER.nf'
include { COLLATE_FUSIONS_PYENV } from './modules/collate_fusions_PYENV.nf'
include { PREDICT_CODING_SEQ_AGFUSION } from './modules/predict_coding_seq_AGFUSION.nf'

// Function to print help message
def helpMessage() {
    log.info"""
Usage:

nextflow run main.nf -profile <local | awsbatch> <--OPTION NAME> <ARGUMENT>

Required Arguments:
---------------
    -c <configFile>      Path to the config file. [REQUIRED]
    -profile             Either <local> for local runs or <awsbatch> for AWS Batch [REQUIRED]
    --manifestPath       Path to tab-delimited manifest file [REQUIRED if inputDir not provided]
                          – must contain sample ID and read1 and read2 local filepaths or remote s3 filepaths
    --inputDir           Path to local directory containing BAM/FASTQ input files [REQUIRED if manifestPath not provided]
    --inputSource        Input source type: <local> for local files, <s3> for S3 files [REQUIRED]
                          – if inputSource is set to <s3>, --inputDir cannot be used and --manifestPath must be provided
    

Optional Arguments:
---------------
    --outputDir          Directory path for output; can be s3 URIs [DEFAULT: ./outputs]
    --trimReads          Set to <false> to skip read trimming on FASTQ input [DEFAULT: true]
    --hlaTypingOnly      Set to <true> to exclusively run HLA typing workflow [DEFAULT: false]
    --deleteIntMedFiles  Set to <true> to delete intermediate files right after they are not needed. [DEFAULT: false]
    --help               Print this help message and exit
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
            // manifest TSV MUST HAVE THESE COLUMN NAMES AS HEADER: sampleName, read1Path, read2Path
            def sampleName = row.sampleName
            def read1 = row.read1Path
            def read2 = row.read2Path
            
            // Validate required fields
            if (!sampleName || !read1 || !read2) {
                log.error "Invalid manifest row: ${row}. Ensure all fields are present."
                return null
            }
            
            // Return tuple of sample name and read file paths
            return tuple(sampleName, [read1, read2])
        }
        .filter { it != null }

    return inputCh
}

// Subworkflow definitions
workflow TRIMMING_WF {
    take:
        readsCh
    main:
        TRIM_READS_FASTP(readsCh)
    emit:
        trimmedCh = TRIM_READS_FASTP.out.trimmed_reads
}

workflow TWOPASS_ALIGNMENT_WF {
    take:
        trimmedCh
    main:
        FIXMATES_MARKDUPES_SAMTOOLS(ALIGN_READS_TWOPASS_STARSAM(trimmedCh).aligned_bams)
    emit:
        alignedBamCh = FIXMATES_MARKDUPES_SAMTOOLS.out.final_bam
}

workflow HLA_TYPING_WF {
    take:
        alignedBamCh
    main:
        TYPE_HLA_ALLELES_ARCASHLA(alignedBamCh)
        // TYPE_HLA_ALLELES_HLAHD(alignedBamCh)
    emit:
        hlaTypingCh = TYPE_HLA_ALLELES_ARCASHLA.out.allotype_json
}

// Main workflow
workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Check that profile is set to one of the allowed values
    if (!workflow.profile) {
        log.error "No profile specified. Please specify -profile <local | awsbatch>"
        exit 1
    } else if (!['local', 'awsbatch'].contains(workflow.profile)) {
        log.error "Invalid profile: ${workflow.profile}. Must be one of: {local, awsbatch}"
        exit 1
    }
    
    // Check that either inputDir or manifestPath is provided
    if (!params.inputDir && !params.manifestPath) {
        log.error "Either --inputDir or --manifestPath must be specified"
        exit 1
    } else if (params.inputDir && params.manifestPath) {
        log.error "Both --inputDir and --manifestPath cannot be specified at the same time"
        exit 1
    }
    // Check the inputSource parameter
    if (!params.inputSource) {
        log.error "Input source must be specified with --inputSource <local | s3>"
        exit 1
    } else if (!['local', 's3'].contains(params.inputSource)) {
        log.error "Invalid input source: ${params.inputSource}. Must be one of: {local, s3}"
        exit 1
    }

    if (params.inputSource == 's3' && !params.manifestPath) {
        log.error "If inputSource is set to 's3', --manifestPath must be provided"
        exit 1
    } else if (params.inputSource == 's3' && params.inputDir) {
        log.error "If inputSource is set to 's3', --inputDir cannot be used"
        exit 1
    } else {
        log.info "Input source is set to ${params.inputSource}"
        // set deleteStagedFiles to false if inputSource is local
        if (params.inputSource == 's3') {
            params.deleteStagedFiles = true
        }
        log.info "deleteIntMedFiles parameter is set to ${params.deleteIntMedFiles}. Intermediate files will be deleted once dependent processes are complete..."
    }

    // Create input channel based on provided input method
    def inputCh
    if (params.manifestPath) {
        inputCh = createInputChannelFromManifest(params.manifestPath)
        log.info "Using manifest file: ${params.manifestPath}"
    } else {
        if (!validateInputDir(params.inputDir)) {
            exit 1
        }
        inputCh = createInputChannelFromPOSIX(params.inputDir)
        log.info "Input files are provided as local directory: ${params.inputDir}"
    }
    
    // Log the key parameters
    log.info "Output directory: ${params.outputDir}"
    log.info "Read trimming: ${params.trimReads}"
    log.info "HLA typing only mode: ${params.hlaTypingOnly}"
    
    // Process input based on trimReads parameter
    def processedInputCh = params.trimReads ? TRIMMING_WF(inputCh).trimmedCh : inputCh

    // Choose alignment workflow based on profile
    def alignedCh = TWOPASS_ALIGNMENT_WF(processedInputCh)
    alignedCh.view()
    
    // Execute workflows based on hlaTypingOnly parameter
    if (params.hlaTypingOnly) {
    // Run only HLA typing
        HLA_TYPING_WF(alignedCh)
    } else {
        // Run full workflow including HLA typing and fusion detection
        // Run HLA typing
        HLA_TYPING_WF(alignedCh)
        // Run fusion calling
        //FUSION_CALLING(processedInputCh)
    }

    // Completion handler
    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        workDir.resolve("stage-${workflow.sessionId}").deleteDir()
    }
}