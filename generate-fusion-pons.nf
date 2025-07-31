#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
// Set default parameter values
params.deleteStagedFiles = params.inputSource == 's3' ? true : false

// Import submodules - only what we need for PoN generation
////// PREPROCESSING MODULES //////////
include { TRIM_READS_FASTP } from './modules/trim_reads_FASTP.nf'

////// FUSION ALIGNMENT MODULES //////////
include { ALIGN_READS_STAR_GENERAL } from './modules/align_reads_STAR_GENERAL.nf'
include { ALIGN_READS_STAR_ARRIBA } from './modules/align_reads_STAR_ARRIBA.nf'
include { FILTER_ALIGNED_READS_EASYFUSE } from './modules/filter_aligned_reads_EASYFUSE.nf'
include { CONVERT_FILTREADS_BAM2FASTQ_EASYFUSE } from './modules/convert_filtreads_bam2fastq_EASYFUSE.nf'

////// FUSION CALLING MODULES //////////
include { CALL_FUSIONS_ARRIBA } from './modules/call_fusions_ARRIBA.nf'
include { CALL_FUSIONS_FUSIONCATCHER } from './modules/call_fusions_FUSIONCATCHER.nf'
include { CALL_FUSIONS_STARFUSION } from './modules/call_fusions_STARFUSION.nf'

// Include PoN-specific modules (you'll need to create these)
include { COLLATE_CUSTOM_PONS_PYENV } from './modules/collate_custom_pons_PYENV.nf'
include { AGGREGATE_CUSTOM_PONS_PYENV } from './modules/aggregate_custom_pons_PYENV.nf'

// Function to print help message
def helpMessage() {
    log.info"""
Usage:

nextflow run generate_fusion-pons.nf -profile <local | awsbatch> <--OPTION NAME> <ARGUMENT>

Description:
This pipeline processes normal samples to generate Panel of Normals (PoN) entries for filtering non-cancer specific fusions in tumor analysis.

Required Arguments:
---------------
    -c <configFile>      Path to the config file. [REQUIRED]
    -profile             Either <local> for local runs or <awsbatch> for AWS Batch [REQUIRED]
    --manifestPath       Path to tab-delimited manifest file [REQUIRED if inputDir not provided]
                          – must contain sample ID and read1 and read2 local filepaths or remote s3 filepaths
    --inputDir           Path to local directory containing BAM/FASTQ input files [REQUIRED if manifestPath not provided]
    --inputSource        Input source type: <local> for local files, <s3> for S3 files [REQUIRED]
                          – if inputSource is set to <s3>, --inputDir cannot be used and --manifestPath must be provided
    --ponsOutputName   Output filename for the generated Panel of Normals [REQUIRED]
                          – this will be used to name the final output file, e.g., <generated_pon.tsv>
    

Optional Arguments:
---------------
    --outputDir          Directory path for output; can be s3 URIs [DEFAULT: ./outputs]
    --existingPoN        Path to existing PoN file to append to [OPTIONAL]
    --trimReads          Set to <false> to skip read trimming on FASTQ input [DEFAULT: true]
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

// Function to read manifest file and create input channels with sample type branching
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
            // manifest TSV MUST HAVE THESE COLUMN NAMES AS HEADER: sampleName, read1Path, read2Path, sampleType
            def sampleName = row.sampleName
            def read1 = row.read1Path
            def read2 = row.read2Path
            def sampleType = row.sampleType
            
            // Validate required fields
            if (!sampleName || !read1 || !read2 || !sampleType) {
                log.error "Invalid manifest row: ${row}. Ensure all fields (sampleName, read1Path, read2Path, sampleType) are present."
                return null
            }
            
            // Return tuple of sample name, read file paths, and sample type
            return tuple(sampleName, [read1, read2], sampleType)
        }
        .filter { it != null }

    return inputCh
}

// Function to branch input channel by sample type
def branchInputChannelBySampleType(inputCh) {
    def branchedCh = inputCh.branch {
        tumor: it[2] == 'Tumor'
        normal: it[2] == 'Normal'
    }
    
    // Convert back to the original format (sampleName, [read1, read2]) for downstream compatibility
    def tumorCh = branchedCh.tumor.map { sampleName, reads, _sampleType -> tuple(sampleName, reads) }
    def normalCh = branchedCh.normal.map { sampleName, reads, _sampleType -> tuple(sampleName, reads) }
    
    return [tumorCh, normalCh]
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

workflow PON_FUSION_CALLING_WF {
    take:
        trimmedCh
        starIndex
        arribaDB
        fuscatDB
        ctatDB

    main:
        // Preprocess reads for fusion calling (same as your aggregate workflow)
        FILTER_ALIGNED_READS_EASYFUSE(ALIGN_READS_STAR_GENERAL(trimmedCh, starIndex).aligned_bam)
        CONVERT_FILTREADS_BAM2FASTQ_EASYFUSE(FILTER_ALIGNED_READS_EASYFUSE.out.filtered_bam)
        filtFastqsCh = CONVERT_FILTREADS_BAM2FASTQ_EASYFUSE.out.filtered_fastqs

        // Gene fusion identification submodule
        CALL_FUSIONS_ARRIBA(ALIGN_READS_STAR_ARRIBA(filtFastqsCh, starIndex).aligned_bam, arribaDB)
        CALL_FUSIONS_FUSIONCATCHER(filtFastqsCh, fuscatDB)
        CALL_FUSIONS_STARFUSION(filtFastqsCh, ctatDB)

        // Join the outputs based on sample name
        CALL_FUSIONS_ARRIBA.out.arriba_fusion_tuple
            .join(CALL_FUSIONS_FUSIONCATCHER.out.fuscat_fusion_tuple)
            .join(CALL_FUSIONS_STARFUSION.out.starfus_fusion_tuple)
            .set { combinedFTFilesCh }

        // Run the combining process with the joined output
        COLLATE_CUSTOM_PONS_PYENV(combinedFTFilesCh)

        // Generate PoN entries from collated fusions
        // Collect all PoN entries and aggregate them
        COLLATE_CUSTOM_PONS_PYENV.out.collatedFTParquet
            .collect()
            .set { allPonEntriesCh }

        // Aggregate all PoN entries into final PoN file
        AGGREGATE_CUSTOM_PONS_PYENV(allPonEntriesCh)

    emit:
        final_pon = AGGREGATE_CUSTOM_PONS_PYENV.out.finalCustomPONs
    
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
        log.error "No profile specified. Please specify -profile < persoMode-local | popMode-local | aws-batch >"
        exit 1
    } else if (!['persoMode-local', 'popMode-local', 'aws-batch'].contains(workflow.profile)) {
        log.error "Invalid profile: ${workflow.profile}. Must be one of: {persoMode-local, popMode-local, aws-batch}"
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

    // check that ponsOutputName is provided
    if (!params.ponsOutputName) {
        log.error "Output filename for the generated Panel of Normals must be specified with --ponsOutputName"
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
        if (params.inputSource == 's3') {
            log.info "deleteStagedFiles parameter is set to <<${params.deleteStagedFiles}>>"
            if (params.deleteStagedFiles) {
                log.info "Staged files will be deleted once dependent processes are complete..."
            } else {
                log.info "Staged files will be kept for debugging purposes."
            }
        } else {
            log.info "Input files will be read from local directory: ${params.inputDir}"
        }
        log.info "deleteIntMedFiles parameter is set to <<${params.deleteIntMedFiles}>>" 
        if (params.deleteIntMedFiles) {
            log.info "Intermediate files will be deleted once dependent processes are complete..."
        } else {
            log.info "Intermediate files will be kept for debugging purposes."
        }
    }

    // Create input channel based on provided input method
    def tumorCh = Channel.empty()
    def normalCh = Channel.empty()
    
    if (params.manifestPath) {
        def inputCh = createInputChannelFromManifest(params.manifestPath)
        log.info "Using manifest file: ${params.manifestPath}"
        
        // Branch the input channel by sample type
        def (branchedTumorCh, branchedNormalCh) = branchInputChannelBySampleType(inputCh)
        tumorCh = branchedTumorCh
        normalCh = branchedNormalCh
        
        // Log sample counts
        tumorCh.count().subscribe { count ->
            log.info "Found ${count} tumor samples!"
        }
        normalCh.count().subscribe { count ->
            log.info "Found ${count} normal samples!"
        }
        
    } else {
        if (!validateInputDir(params.inputDir)) {
            exit 1
        }
        normalCh = createInputChannelFromPOSIX(params.inputDir)
        log.info "Input files are provided as local directory: ${params.inputDir}"
        log.warn "All samples from input directory will be processed as NORMAL samples! If you have tumor samples mixed in, please provide a manifest file with sampleType column at runtime instead of --inputDir !!!"
    }
    
    // Log the key parameters
    log.info "Output directory: ${params.outputDir}"
    log.info "Read trimming: ${params.trimReads}"
    log.info "PoN output file: ${params.ponsOutputName}"
    
    // Process input based on trimReads parameter
    def qcProcInputCh = params.trimReads ? TRIMMING_WF(normalCh).trimmedCh : normalCh

    // Run PoN fusion calling workflow
    PON_FUSION_CALLING_WF(qcProcInputCh, 
        params.starIndex, 
        params.arribaDB, 
        params.fuscatDB, 
        params.ctatDB)

    // Completion handler
    workflow.onComplete = {
        println "Panel of Normals generation completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }
}
