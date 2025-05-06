#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import submodules
////// PREPROCESSING MODULES //////////
include { TRIM_READS_FASTP } from './modules/trim_reads_FASTP.nf'
include { ALIGN_READS_TWOPASS_STARSAM } from './modules/align_reads_twopass_STARSAM.nf'
include { FIXMATES_MARKDUPES_SAMTOOLS } from './modules/fixmates_markdupes_SAMTOOLS.nf'

////// FUSION ALIGNMENT MODULES //////////
include { ALIGN_READS_STAR_ARRIBA } from './modules/align_reads_STAR_ARRIBA.nf'
include { ALIGN_READS_STAR_GENERAL } from './modules/align_reads_STAR_GENERAL.nf'
include { FILTER_ALIGNED_READS_EASYFUSE } from './modules/filter_aligned_reads_EASYFUSE.nf'
include { CONVERT_FILT_BAMS2FASTQ_EASYFUSE } from './modules/convert_filt_bams2fastq_EASYFUSE.nf'

////// FUSION CALLING MODULES //////////
include { CALL_FUSIONS_ARRIBA } from './modules/call_fusions_ARRIBA.nf'
include { CALL_FUSIONS_FUSIONCATCHER } from './modules/call_fusions_FUSIONCATCHER.nf'
include { COMBINE_FUSIONS_PYENV } from './modules/combine_fusions_PYENV.nf'

////// FUSION FILTERING AND ANNOTATION MODULES //////////
include { FILTER_FUSIONS_PYENV } from './modules/filter_fusions_PYENV.nf'
include { PREDICT_CODING_SEQ_AGFUSION } from './modules/predict_coding_seq_AGFUSION.nf'

/////// HLA TYPING MODULES //////////
include { TYPE_HLA_ALLELES_ARCASHLA } from './modules/type_hla_alleles_ARCASHLA.nf'
include { EXTRACT_HLATYPING_JSONS_PYENV } from './modules/extract_hlatyping_jsons_PYENV.nf'
include { COMBINE_HLA_FILES_BASH } from './modules/combine_hla_files_BASH.nf'



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
        FIXMATES_MARKDUPES_SAMTOOLS(ALIGN_READS_TWOPASS_STARSAM(trimmedCh).aligned_bam)
    emit:
        alignedBamCh = FIXMATES_MARKDUPES_SAMTOOLS.out.final_bam
}

workflow GENERAL_ALIGNMENT_WF {
    take:
        trimmedCh
    main:
        FILTER_ALIGNED_READS_EASYFUSE(ALIGN_READS_STAR_GENERAL(trimmedCh).aligned_bam)
        CONVERT_FILT_BAMS2FASTQ_EASYFUSE(FILTER_ALIGNED_READS_EASYFUSE.out.filtered_bam)
    emit:
        filtFastqsCh = CONVERT_FILT_BAMS2FASTQ_EASYFUSE.out.filtered_fastqs
}

workflow AGGREGATE_FUSION_CALLING_WF {
    take:
        filtFastqsCh
    main:
        // gene fusion identification submodule
        CALL_FUSIONS_ARRIBA(ALIGN_READS_STAR_ARRIBA(filtFastqsCh).aligned_bam)
        CALL_FUSIONS_FUSIONCATCHER(filtFastqsCh)

        // Join the outputs based on sample name
        CALL_FUSIONS_ARRIBA.out.arriba_fusion_tuple
            .join(CALL_FUSIONS_FUSIONCATCHER.out.fuscat_fusion_tuple)
            .set { combinedFTFilesCh }
        
        //combinedFTFilesCh.view()

        // Run the combining process with the joined output then channel into filtering process
        FILTER_FUSIONS_PYENV(COMBINE_FUSIONS_PYENV(combinedFTFilesCh).combinedFTParquet)

    emit:
        filteredFusionCh = FILTER_FUSIONS_PYENV.out.filteredFusionList
        
}

workflow PREDICT_CODING_SEQ_WF {
    take:
        filteredFusionCh
    main:
        PREDICT_CODING_SEQ_AGFUSION(filteredFusionCh)

}

workflow HLA_TYPING_WF {
    take:
        alignedBamCh
    main:
        EXTRACT_HLATYPING_JSONS_PYENV(TYPE_HLA_ALLELES_ARCASHLA(alignedBamCh).allotype_json)

    emit:
        hlaTypingCh = EXTRACT_HLATYPING_JSONS_PYENV.out.hlaTypingTsv
}

workflow HLA_TYPE_COLLATION_WF {
    take:
        hlaTypingCh
    main:
        COMBINE_HLA_FILES_BASH(hlaTypingCh)
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
        log.info "deleteIntMedFiles parameter is set to ${params.deleteIntMedFiles}." 
        if (params.deleteIntMedFiles) {
            log.info "Intermediate files will be deleted once dependent processes are complete..."
        } else {
            log.info "Intermediate files will be kept for debugging purposes."
        }
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
    def qcProcInputCh = params.trimReads ? TRIMMING_WF(inputCh).trimmedCh : inputCh


    ///////// Two-pass STAR alignment workflow ////////////
    def alignedCh = TWOPASS_ALIGNMENT_WF(qcProcInputCh)
    ////////////////////////////////////////////////////////


    // Execute workflow branching based on hlaTypingOnly parameter
    if (params.hlaTypingOnly) {
    // Run only HLA typing
        def hlaFilesCh = HLA_TYPING_WF(alignedCh)
        // collate HLA typing results
        HLA_TYPE_COLLATION_WF(hlaFilesCh.collect(flat: false))

        
        // // combined HLA typing results
        // reformattedHlaFiles
        // .map { sampleId, file -> 
        // // Read file content and combine with sample ID
        //     def content = file.text.trim()
        //     return [sampleId, content]
        // }
        // .collectFile(
        //     name: 'combined_hla_types.tsv',
        //     newLine: true,
        //     seed: "SampleID\tHLA_Types",  // Header
        //     storeDir: "${params.outputDir}/combined-HLA-types"
        // ) { sampleId, content ->
        //     return "${sampleId}\t${content}"
        // }

    } else {
        // HLA typing
        def hlaFilesCh = HLA_TYPING_WF(alignedCh)
        // collate HLA typing results
        HLA_TYPE_COLLATION_WF(hlaFilesCh.collect(flat: false))
        
        // Run the general alignment workflow
        def filtFastqsCh = GENERAL_ALIGNMENT_WF(qcProcInputCh)
        
        // Fusion calling
        def filteredFusionCh = AGGREGATE_FUSION_CALLING_WF(filtFastqsCh)

        // run AGFUSION coding sequence prediction
        PREDICT_CODING_SEQ_WF(filteredFusionCh)
    }

    // Completion handler
    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        //workDir.resolve("stage-${workflow.sessionId}").deleteDir()
    }
}