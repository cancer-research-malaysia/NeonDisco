#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
// Set default parameter values
params.deleteStagedFiles = (params.inputSource == 's3')

// Import submodules
////// PREPROCESSING MODULES //////////
include { TRIM_READS_FASTP } from './modules/trim_reads_FASTP.nf'
include { ALIGN_READS_TWOPASS_STARSAM } from './modules/align_reads_twopass_STARSAM.nf'
include { FIXMATES_MARKDUPES_SAMTOOLS } from './modules/fixmates_markdupes_SAMTOOLS.nf'

////// FUSION ALIGNMENT MODULES //////////
include { ALIGN_READS_STAR_ARRIBA } from './modules/align_reads_STAR_ARRIBA.nf'
include { ALIGN_READS_STAR_GENERAL } from './modules/align_reads_STAR_GENERAL.nf'
include { FILTER_ALIGNED_READS_EASYFUSE } from './modules/filter_aligned_reads_EASYFUSE.nf'
include { CONVERT_FILTREADS_BAM2FASTQ_EASYFUSE } from './modules/convert_filtreads_bam2fastq_EASYFUSE.nf'

////// FUSION CALLING MODULES //////////
include { CALL_FUSIONS_ARRIBA } from './modules/call_fusions_ARRIBA.nf'
include { CALL_FUSIONS_FUSIONCATCHER } from './modules/call_fusions_FUSIONCATCHER.nf'
include { CALL_FUSIONS_STARFUSION } from './modules/call_fusions_STARFUSION.nf'
include { COLLATE_FUSIONS_PYENV } from './modules/collate_fusions_PYENV.nf'
include { WRANGLE_RAW_FUSIONS_PYENV } from './modules/wrangle_raw_fusions_PYENV.nf'

////// FUSION FILTERING AND ANNOTATION MODULES //////////
include { FILTER_FUSIONS_PYENV } from './modules/filter_fusions_PYENV.nf'
include { COLLECT_COHORTWIDE_UNFILTERED_FUSIONS_PYENV } from './modules/collect_cohortwide_unfiltered_fusions_PYENV.nf'
include { COLLECT_COHORTWIDE_NORMFILTERED_FUSIONS_PYENV } from './modules/collect_cohortwide_normfiltered_fusions_PYENV.nf'
include { GET_COHORT_RECURRENT_FUSIONS_PYENV } from './modules/get_cohort_recurrent_fusions_PYENV.nf'
include { TRANSLATE_IN_SILICO_AGFUSION } from './modules/translate_in_silico_AGFUSION.nf'
include { VALIDATE_IN_SILICO_FUSIONINSPECTOR } from './modules/validate_in_silico_FUSIONINSPECTOR.nf'
include { COLLECT_COHORTWIDE_VALIDATED_FUSIONS_PYENV } from './modules/collect_cohortwide_validated_fusions_PYENV.nf'

////// FUSION NEOEPITOPE PREDICTION MODULES //////////
include { KEEP_VALIDATED_FUSIONS_PYENV } from './modules/keep_validated_fusions_PYENV.nf'
include { FILTER_VALIDATED_FUSIONS_FOR_RECURRENT_PYENV } from './modules/filter_validated_fusions_for_recurrent_PYENV.nf'
include { PREDICT_NEOPEPTIDES_SAMPLE_LEVEL_HLAS_PVACFUSE } from './modules/predict_neopeptides_sample_level_hlas_PVACFUSE.nf'
include { PREDICT_NEOPEPTIDES_COHORT_LEVEL_HLAS_PVACFUSE } from './modules/predict_neopeptides_cohort_level_hlas_PVACFUSE.nf'

////// HLA TYPING MODULES //////////
include { TYPE_HLAS_WITH_FALLBACK_ARCASHLA } from './modules/type_hlas_with_fallback_ARCASHLA.nf'
include { REFORMAT_AND_COLLATE_HLA_RESULTS_PYENV } from './modules/reformat_and_collate_hla_results_PYENV.nf'
include { FILTER_HLA_ALLOTYPES_FREQ_PYENV } from './modules/filter_hla_allotypes_freq_PYENV.nf'


// Function to print help message
def helpMessage() {
    log.info"""
Usage:

nextflow run neondisco-main.nf -profile <EXECUTOR,MODE[,RESOURCE]> <--OPTION NAME> <ARGUMENT>

Profile Examples:
----------------
    -profile local,sampleNeoPredMode              # Run locally in sample-level neoantigen mode
    -profile local,sharedNeoPredMode              # Run locally in shared neoantigen mode  
    -profile local,dualNeoPredMode                # Run locally in both modes
    -profile awsbatch,sampleNeoPredMode           # Run on AWS Batch in sample-level neoantigen mode
    -profile awsbatch,sharedNeoPredMode           # Run on AWS Batch in shared neoantigen mode
    -profile awsbatch,dualNeoPredMode             # Run on AWS Batch in both modes

Required Arguments:
---------------
    -c <configFile>             Path to the config file. [REQUIRED]
    -profile                    Comma-separated list of profiles [REQUIRED]
                                    EXECUTOR:   local | awsbatch
                                    MODE:       sampleNeoPredMode | sharedNeoPredMode | dualNeoPredMode
    --inputDir                  Path to local directory containing BAM/FASTQ input files [REQUIRED if manifestPath not provided]
    --manifestPath              Path to tab-delimited manifest file [REQUIRED if inputDir not provided]
                                    – must contain sample ID and read1 and read2 local filepaths or remote s3 filepaths
                                    – must also contain sampleType column with values 'Tumor' or 'Normal'
    --inputSource               Input source type: <local> for local files, <s3> for S3 files [REQUIRED]
                                    – if inputSource is set to <s3>, --inputDir cannot be used and --manifestPath must be provided
                                    -if inputSource is set to <local> -profile 'awsbatch' is disallowed

Optional Arguments:
---------------
    --outputDir                 Directory path for output; can be s3 URIs [DEFAULT: ./outputs]
    --trimReads                 Skip read trimming on FASTQ input [DEFAULT: true]
    --hlaTypingOnly             Exclusively run HLA typing subworkflow [DEFAULT: false]
    --includeNeoPred            Run neopeptide prediction subworkflow [DEFAULT: true]
    --deleteIntMedFiles         Delete intermediate files right after they are not needed [DEFAULT: false]
    --deleteStagedFiles         Delete staged files after processing [DEFAULT: true if inputSource is 's3']
    --sampleHLANeoPred          Run neopeptide prediction using sample-level HLAs [DEFAULT: varies by mode]
    --cohortHLANeoPred          Run neopeptide prediction using cohort-level HLAs [DEFAULT: varies by mode]
    --recurrentFusionsOnly      Analyze recurrent fusions only [DEFAULT: varies by mode]
                                    – true: Only process cohort-wide recurrent fusions (skip if none found)
                                    – false: Process all validated fusions regardless of recurrence
    --recurrenceThreshold       Threshold for recurrent fusions; only applies if recurrentFusionsOnly is set to true
                                [DEFAULT: 0.005] (0.5% recurrence)
    --help                      Print this help message and exit
    """.stripIndent()
}

// Function to validate profiles
def validateProfiles() {
    if (!workflow.profile) {
        log.error "No profile specified. Please specify -profile <EXECUTOR,MODE[,RESOURCE]>"
        log.error "Examples: -profile local,sampleNeoPredMode or -profile awsbatch,sharedNeoPredMode"
        return false
    }
    
    def profiles = workflow.profile.split(',').collect { it.trim() }
    def executors = ['local', 'awsbatch']
    def modes = ['sampleNeoPredMode', 'sharedNeoPredMode', 'dualNeoPredMode']
    
    // Check for required executor and mode
    def hasExecutor = profiles.any { executors.contains(it) }
    def hasMode = profiles.any { modes.contains(it) }
    
    if (!hasExecutor) {
        log.error "You must specify an executor profile: ${executors.join(', ')}"
        return false
    }
    
    if (!hasMode) {
        log.error "You must specify a discovery mode profile: ${modes.join(', ')}"
        return false
    }
    
    // Check for multiple executors (still not allowed)
    def executorCount = profiles.count { executors.contains(it) }
    if (executorCount > 1) {
        log.error "Only one executor profile allowed. Found: ${profiles.findAll { executors.contains(it) }.join(', ')}"
        return false
    }
    
    // Check for mode combinations - now only allow one mode since we have dualNeoPredMode
    def modeCount = profiles.count { modes.contains(it) }
    if (modeCount > 1) {
        log.error "Only one neoantigen prediction mode allowed. Found: ${profiles.findAll { modes.contains(it) }.join(', ')}"
        log.error "Use 'dualNeoPredMode' to run both sample-level and shared predictions"
        return false
    }
    
    // Log which mode is being used
    if (profiles.contains('dualNeoPredMode')) {
        log.info "Running in dual mode: both sample-level and shared neoantigen discovery enabled"
    } else if (profiles.contains('sampleNeoPredMode')) {
        log.info "Running in sample-level neoantigen discovery mode"
    } else if (profiles.contains('sharedNeoPredMode')) {
        log.info "Running in shared neoantigen discovery mode"
    }
    
    // Validate unknown profiles
    def validProfiles = executors + modes
    def unknownProfiles = profiles.findAll { !validProfiles.contains(it) }
    if (unknownProfiles) {
        log.error "Unknown profile(s): ${unknownProfiles.join(', ')}"
        log.error "Valid profiles: ${validProfiles.join(', ')}"
        log.error "Please rerun by specifying a valid profile combination from the list above."
        return false
    }
    
    log.info "Using profiles: ${profiles.join(', ')}"
    return true
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


//////////////// Workflow for other discovery modules ////////////////
// workflow TWOPASS_READS_ALIGNMENT_WF {
//     take:
//         trimmedCh
//         starIndex
//     main:
//         FIXMATES_MARKDUPES_SAMTOOLS(ALIGN_READS_TWOPASS_STARSAM(trimmedCh, starIndex).aligned_bam)
//     emit:
//         alignedBam2PassCh = FIXMATES_MARKDUPES_SAMTOOLS.out.final_bam
// }

workflow GENERAL_READS_ALIGNMENT_WF {
    take:
        trimmedCh
        starIndex
    main:
        ALIGN_READS_STAR_GENERAL(trimmedCh, starIndex)
    emit:
        alignedBamCh = ALIGN_READS_STAR_GENERAL.out.final_bam.view()
}

/////////simplified with single HLA typing process with HLA-HD fallback
workflow HLA_TYPING_WITH_FALLBACK_WF {
    take:
        alignedBamCh  // [sampleName, bam, bamIdx]
    
    main:
        // Single process handles both arcasHLA and HLA-HD fallback
        TYPE_HLAS_WITH_FALLBACK_ARCASHLA(alignedBamCh)
    
        // Collect all JSON files
        all_json_files = TYPE_HLAS_WITH_FALLBACK_ARCASHLA.out.hla_json
            .map { _sampleName, json -> json }
            .collect()
    
        // Single process to reformat and collate everything
        REFORMAT_AND_COLLATE_HLA_RESULTS_PYENV(all_json_files)
    
    emit:
        sampleSpecificHLAsTsv = REFORMAT_AND_COLLATE_HLA_RESULTS_PYENV.out.cohortWideHLAList
        individualResultDir = REFORMAT_AND_COLLATE_HLA_RESULTS_PYENV.out.individualResults
}

//////////////////////////////////////////////////////////////////////////
workflow COLLECT_COHORTWIDE_UNFILTERED_FUSIONS {
    take:
        aggregatedTuplesCh
    main:

        // Collate fusion calls from different callers
        COLLATE_FUSIONS_PYENV(aggregatedTuplesCh)

        collatedFusionsforWrangling = COLLATE_FUSIONS_PYENV.out.collatedFTParquet

        // Wrangle raw fusions
        wrangledFusionsTsv = WRANGLE_RAW_FUSIONS_PYENV(collatedFusionsforWrangling).wrangledUnfilteredFusionsTsv

        // Collect cohort-wide unfiltered fusions
        COLLECT_COHORTWIDE_UNFILTERED_FUSIONS_PYENV(wrangledFusionsTsv
                            .collect { _meta, filepath -> filepath }
                        )
}

workflow COLLATE_FILTER_FUSIONS {
    take:
        aggregatedTuplesCh
        metaDataDir
    main:
        // Collate fusion calls from different callers
        collatedFusionsforFiltering = COLLATE_FUSIONS_PYENV(aggregatedTuplesCh).collatedFTParquet

        // Run the combining process with the joined output then channel into filtering process
        FILTER_FUSIONS_PYENV(collatedFusionsforFiltering, metaDataDir)

    emit:
        normFilteredFusionsCh = FILTER_FUSIONS_PYENV.out.filteredFusions
        uniqueFiltFusionPairsForFusInsCh = FILTER_FUSIONS_PYENV.out.uniqueFiltFusionPairsForFusIns
}

workflow COLLECT_COHORTWIDE_NORMFILTERED_FUSIONS {
    take:
        normFilteredFusionsCh
    main:
        // Collect cohort-wide normfiltered fusions
        COLLECT_COHORTWIDE_NORMFILTERED_FUSIONS_PYENV(normFilteredFusionsCh)
    emit:
        cohortwideFusionsFile = COLLECT_COHORTWIDE_NORMFILTERED_FUSIONS_PYENV.out.cohortwideFusionsFile
}

workflow AGGREGATE_FUSION_CALLING_WF {
    take:
        alignedBamCh
        starIndex
        arribaDB
        fuscatDB
        ctatDB
        metaDataDir
    main:
        ////// Preprocess reads for aggregate fusion calling
        FILTER_ALIGNED_READS_EASYFUSE(alignedBamCh)
        CONVERT_FILTREADS_BAM2FASTQ_EASYFUSE(FILTER_ALIGNED_READS_EASYFUSE.out.filtered_bam)
        filtFastqsCh = CONVERT_FILTREADS_BAM2FASTQ_EASYFUSE.out.filtered_fastqs

        ////// Run gene fusion identification submodules
        CALL_FUSIONS_ARRIBA(ALIGN_READS_STAR_ARRIBA(filtFastqsCh, starIndex).aligned_bam, arribaDB)
        CALL_FUSIONS_FUSIONCATCHER(filtFastqsCh, fuscatDB)
        CALL_FUSIONS_STARFUSION(filtFastqsCh, ctatDB)

        //////////////////////////////////////////////////////////////
        /////////////// no EasyFuse filtering variant ////////////////
        ////// gene fusion identification submodule
        // CALL_FUSIONS_ARRIBA(ALIGN_READS_STAR_ARRIBA(trimmedFqs, starIndex).aligned_bam, arribaDB)
        // CALL_FUSIONS_FUSIONCATCHER(trimmedFqs, fuscatDB)
        // CALL_FUSIONS_STARFUSION(trimmedFqs, ctatDB)
        //////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////

        // Join the outputs based on sample name
        CALL_FUSIONS_ARRIBA.out.arriba_fusion_tuple
           .join(CALL_FUSIONS_FUSIONCATCHER.out.fuscat_fusion_tuple)
           .join(CALL_FUSIONS_STARFUSION.out.starfus_fusion_tuple)
           .set { combinedFTFilesCh }

        // Collect cohort-wide unfiltered fusions
        COLLECT_COHORTWIDE_UNFILTERED_FUSIONS(combinedFTFilesCh)

        // Collate and filter fusions
        COLLATE_FILTER_FUSIONS(combinedFTFilesCh, metaDataDir)

        //// Collect cohort-wide normfiltered fusions
        COLLECT_COHORTWIDE_NORMFILTERED_FUSIONS(COLLATE_FILTER_FUSIONS.out.normFilteredFusionsCh)

    emit:
        normFilteredFusionsCh = COLLATE_FILTER_FUSIONS.out.normFilteredFusionsCh
        uniqueFiltFusionPairsForFusInsCh = COLLATE_FILTER_FUSIONS.out.uniqueFiltFusionPairsForFusInsCh
        cohortwideFusionsFile = COLLECT_COHORTWIDE_NORMFILTERED_FUSIONS.out.cohortwideFusionsFile
    
}

////////////////////////////////////////////////////////////////////////////////////
workflow IN_SILICO_TRANSCRIPT_VALIDATION_WF {
    take:
        normFilteredFusionsCh
        uniqueFiltFusionPairsForFusInsCh
        trimmedCh
        ctatDB
    main:
        agfusionOutput = TRANSLATE_IN_SILICO_AGFUSION(normFilteredFusionsCh)

        // Join channels and fix the tuple structure
        joinedInputs = uniqueFiltFusionPairsForFusInsCh
            .join(trimmedCh, by: 0)
            .join(agfusionOutput.filtered_agfusion_outdir, by: 0)
            .map { sampleName, uniqueFiltPairs, trimmedReads, agfusionDir -> 
                tuple(sampleName, agfusionDir, uniqueFiltPairs, trimmedReads)
            }
        
        VALIDATE_IN_SILICO_FUSIONINSPECTOR(joinedInputs, ctatDB)

    emit:
        fusInspectorTsv = VALIDATE_IN_SILICO_FUSIONINSPECTOR.out.fusInspectorTsv
        filteredAgfusionOutdir = TRANSLATE_IN_SILICO_AGFUSION.out.filtered_agfusion_outdir 
        
}

workflow RECURRENT_FUSION_FILTERING_WF {
    take:
        cohortwideFusionsFile
    main:
        // Get cohort recurrent fusions
        GET_COHORT_RECURRENT_FUSIONS_PYENV(cohortwideFusionsFile)
    emit:
        cohortRecurrentFusionsCh = GET_COHORT_RECURRENT_FUSIONS_PYENV.out.cohortRecurrentFusionTsv
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
workflow COLLECT_COHORTWIDE_FI_VALIDATED_FUSIONS {
    take:
        validatedFusionsTsvs
    main:
        // get cohort-wide validated fusions
        COLLECT_COHORTWIDE_VALIDATED_FUSIONS_PYENV(validatedFusionsTsvs)
    emit:
        cohortValidatedFusions = COLLECT_COHORTWIDE_VALIDATED_FUSIONS_PYENV.out.cohortwideValidatedFusionsFile

}

workflow SAMPLE_LEVEL_HLA_NEOANTIGENS {
    take:
        agfusionFinalDir
        sampleSpecificHLAsTsv
        metaDataDir
    main:
        PREDICT_NEOPEPTIDES_SAMPLE_LEVEL_HLAS_PVACFUSE(
            agfusionFinalDir, 
            sampleSpecificHLAsTsv,
            metaDataDir
        )
}

workflow COHORT_LEVEL_HLA_NEOANTIGENS {
    take:
        agfusionFinalDir
        sampleSpecificHLAsTsv
        metaDataDir
    main:
        // Filter HLA file based on sample count
        validHLAForCohort = sampleSpecificHLAsTsv
            .filter { file ->
                def sampleCount = file.countLines() - 1
                log.info "Found ${sampleCount} sample(s) in HLA file..."
                if (sampleCount >= 5) {
                    log.info "Neopeptide prediction using cohort-level HLA types will be executed with ${sampleCount} input sample(s)."
                    return true
                } else {
                    log.warn "Neopeptide prediction using cohort-level HLA types skipped: only ${sampleCount} input sample(s) available (minimum 5 required)..."
                    return false
                }
            }

        // Filter HLA allotypes by frequency
        FILTER_HLA_ALLOTYPES_FREQ_PYENV(validHLAForCohort)
        
        // Combine agfusion dirs with cohort HLA
        cohortInputs = agfusionFinalDir
                        .combine(FILTER_HLA_ALLOTYPES_FREQ_PYENV.out.hlaAllotypesFreqFilteredString)

        PREDICT_NEOPEPTIDES_COHORT_LEVEL_HLAS_PVACFUSE(
            cohortInputs.map { sampleName, agfusionDir, _cohortHLA -> 
                tuple(sampleName, agfusionDir) 
            },
            cohortInputs.map { _sampleName, _agfusionDir, cohortHLA -> 
                cohortHLA 
            }.first(),
            metaDataDir
        )
    
}

// Updated main neopeptide workflow - now just orchestrates the sub-workflows
workflow NEOANTIGEN_PREDICTION_WF {
    take:
        fusInspectorTsv
        filteredAgfusionOutdir
        normFilteredFusionsCh
        cohortRecurrentFusionsCh
        sampleSpecificHLAsTsv
        metaDataDir
    
    main:

        // join the fusInspectorTsv, filteredAgfusionOutdir, and normFilteredFusionsCh channel by sampleName
        joinedInputs = fusInspectorTsv
            .join(filteredAgfusionOutdir, by: 0)
            .join(normFilteredFusionsCh, by: 0)
            .map { sampleName, fusInspectorFile, agfusionDir, filteredFusions -> 
                tuple(sampleName, fusInspectorFile, agfusionDir, filteredFusions) 
            }
        
        // Preprocess agfusion output for neoepitope prediction
        KEEP_VALIDATED_FUSIONS_PYENV(joinedInputs)

        // Collect cohort-wide validated fusions
        COLLECT_COHORTWIDE_FI_VALIDATED_FUSIONS(KEEP_VALIDATED_FUSIONS_PYENV.out.validatedFusions
                            .collect { _meta, filepath -> filepath }
                        )

        validatedAgfusionDir = KEEP_VALIDATED_FUSIONS_PYENV.out.validatedAgfusionDir
        validatedFusions = KEEP_VALIDATED_FUSIONS_PYENV.out.validatedFusions

        // join the validatedFusions and validatedAgfusionDir channels by sampleName
        joinedValidatedFusionsDat = validatedFusions
            .join(validatedAgfusionDir, by: 0)
            .map { sampleName, validatedFusionsFile, validatedDir -> 
                tuple(sampleName, validatedFusionsFile, validatedDir) 
            }

        // Filter validated fusions for recurrent ones
        FILTER_VALIDATED_FUSIONS_FOR_RECURRENT_PYENV(
            joinedValidatedFusionsDat,
            cohortRecurrentFusionsCh
        )
        
        // Get the full validated fusions directory
        recurrentValidatedDir = FILTER_VALIDATED_FUSIONS_FOR_RECURRENT_PYENV.out.validatedRecurrentAgfusionDir
        
        // Logic based on recurrentFusionsOnly parameter
        def finalAgfusionDir = Channel.empty()

        if (params.recurrentFusionsOnly) {
            // Default mode: Only process recurrent fusions
            finalAgfusionDir = recurrentValidatedDir
                .ifEmpty { 
                    log.warn "No recurrent fusions found in this input cohort. Neoantigen prediction will be skipped."
                    log.info "Consider using [--recurrentFusionsOnly false] flag to process all validated fusions instead."
                    Channel.empty()
                }   
        } else {
            // Alternative mode: Process all validated fusions
            finalAgfusionDir = validatedAgfusionDir
            log.info "Processing all validated fusions [--recurrentFusionsOnly false]..."
        }
        
        // Run neoepitope prediction subworkflow
        if (params.includeNeoPred) {
            log.info "Running neoepitope prediction subworkflow..."
            // Run sample-specific if enabled
            if (params.sampleHLANeoPred) {
                SAMPLE_LEVEL_HLA_NEOANTIGENS(finalAgfusionDir, sampleSpecificHLAsTsv, metaDataDir)
            }

            // Run cohort-wide if enabled
            if (params.cohortHLANeoPred) {
                COHORT_LEVEL_HLA_NEOANTIGENS(finalAgfusionDir, sampleSpecificHLAsTsv, metaDataDir)
            }
        } else {
            log.info "Skipping neoepitope prediction subworkflow."
        }
}

// Main workflow
workflow {
    // Show help message if requested
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Validate profiles
    if (!validateProfiles()) {
        exit 1
    }

    // Check that either inputDir or manifestPath is provided
    if (!params.inputDir && !params.manifestPath) {
        log.error "Either --inputDir or --manifestPath must be specified."
        exit 1
    } else if (params.inputDir && params.manifestPath) {
        log.error "Both --inputDir and --manifestPath cannot be specified at the same time."
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

    // Validate inputSource and executor compatibility
    def profilesList = workflow.profile.split(',').collect { it.trim() }
    def isAwsBatch = profilesList.contains('awsbatch')
    def isLocal = profilesList.contains('local')
    
    if (isLocal && isAwsBatch) {
        log.error "AWS Batch executor is not compatible with local input source."
        log.error "Please use either:"
        log.error "  - Local executor with local input: -profile local,<MODE> --inputSource local"
        log.error "  - AWS Batch executor with S3 input: -profile awsbatch,<MODE> --inputSource s3"
        exit 1
    }

    if (params.inputSource == 's3' && !params.manifestPath) {
        log.error "If inputSource is set to 's3', --manifestPath must be provided."
        exit 1
    } else if (params.inputSource == 's3' && params.inputDir) {
        log.error "If inputSource is set to 's3', --inputDir cannot be used."
        exit 1
    } else {
        log.info "Input source is set to ${params.inputSource}"
        if (params.inputSource == 's3') {
            log.info "deleteStagedFiles parameter is set to << ${params.deleteStagedFiles} >>"
            if (params.deleteStagedFiles) {
                log.info "----Staged files will be deleted once dependent processes are complete..."
            } else {
                log.info "----Staged files will be kept for debugging purposes."
            }
        } else {
            log.info "Input files will be read from local directory: << ${params.inputDir} >>"
        }
        log.info "deleteIntMedFiles parameter is set to << ${params.deleteIntMedFiles} >>" 
        if (params.deleteIntMedFiles) {
            log.info "----Intermediate files will be deleted once dependent processes are complete..."
        } else {
            log.info "----Intermediate files will be kept for debugging purposes."
        }
    }

    // Create input channel based on provided input method
    def tumorCh = Channel.empty()
    def normalCh = Channel.empty()
    
    if (params.manifestPath) {
        def inputCh = createInputChannelFromManifest(params.manifestPath)
        log.info "Using manifest file: << ${params.manifestPath} >>"
        
        // Branch the input channel by sample type
        def (branchedTumorCh, branchedNormalCh) = branchInputChannelBySampleType(inputCh)
        tumorCh = branchedTumorCh.view()
        normalCh = branchedNormalCh.view()
        
        // Log sample counts
        normalCh.count().subscribe { count ->
            log.info "Found ${count} normal samples! ----Normal samples will not be processed in the main NeonDisco pipeline."
        }
        tumorCh.count().subscribe { count ->
            log.info "Found ${count} tumor samples! ----Initializing NeonDisco pipeline..."
            log.info ""
        }
        
    } else {
        if (!validateInputDir(params.inputDir)) {
            exit 1
        }
        tumorCh = createInputChannelFromPOSIX(params.inputDir)
        log.info "Input files are provided as local directory: << ${params.inputDir} >>"
        log.warn "----All samples from input directory will be processed as TUMOR samples! If you have normal samples mixed in, please provide a manifest file with sampleType column at runtime instead of --inputDir."
    }
    
    // Log the key parameters
    log.info "Output directory: << ${params.outputDir} >>"
    log.info "Read trimming: << ${params.trimReads} >>"
    log.info "HLA typing only? : << ${params.hlaTypingOnly} >>"
    log.info "Include neopeptide prediction? : << ${params.includeNeoPred} >>"
    // Log the fusion filtering mode
    def mode = params.recurrentFusionsOnly ? "Recurrent-only" : "All-validated"
    log.info "Fusion-derived neoantigen prediction input set: << ${mode} fusions >>"
    
    if (params.recurrentFusionsOnly) {
        log.info "----Recurrence threshold: << ${params.recurrenceThreshold * 100}% >>"
    } else {
        log.info "----Recurrent fusion–only parameter is disabled; all validated fusions will be processed"
    }
    log.info "Neopeptide prediction mode [Personalized (using sample-level HLA allotypes)]: << ${params.sampleHLANeoPred} >>"
    log.info "Neopeptide prediction mode [Cohort-based (using cohort-level HLA allotypes)]: << ${params.cohortHLANeoPred} >>"
    
    // Process tumor samples only (normal channel remains unused but available)
    def qcProcInputCh = params.trimReads ? TRIMMING_WF(tumorCh).trimmedCh : tumorCh

    ///////// STAR alignment workflow ////////////
    def alignedBamsCh = GENERAL_READS_ALIGNMENT_WF(qcProcInputCh, params.starIndex)
    ////////////////////////////////////////////////////////

    // Execute workflow branching based on hlaTypingOnly parameter
    if (params.hlaTypingOnly) {
        
        log.info "Running HLA typing only subworkflow..."
        // Run only HLA typing
        HLA_TYPING_WITH_FALLBACK_WF(alignedBamsCh)

    } else {
        
        //// HLA typing
        HLA_TYPING_WITH_FALLBACK_WF(alignedBamsCh)

        //// Fusion calling module
        AGGREGATE_FUSION_CALLING_WF(alignedBamsCh, params.starIndex, 
                params.arribaDB, params.fuscatDB, params.ctatDB, params.metaDataDir)

        /////////////////////////////////////////////////////////
        //// ALTERNATIVE MODE WITHOUT EASYFUSE PRE-FILTERING ////
        // AGGREGATE_FUSION_CALLING_WF(params.starIndex, 
        //         params.arribaDB, params.fuscatDB, params.ctatDB, params.metaDataDir, qcProcInputCh)
        ////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////
        // RNA-editing calling module
        //def _alignedBams2passCh = TWOPASS_READS_ALIGNMENT_WF(qcProcInputCh, params.starIndex)
        ////////////////////////////////////////////////////////


        //// AGFUSION coding sequence prediction
        IN_SILICO_TRANSCRIPT_VALIDATION_WF(
                AGGREGATE_FUSION_CALLING_WF.out.normFilteredFusionsCh,
                AGGREGATE_FUSION_CALLING_WF.out.uniqueFiltFusionPairsForFusInsCh,
                qcProcInputCh,
                params.ctatDB
            )

        //// recurrent fusion filtering
        RECURRENT_FUSION_FILTERING_WF(AGGREGATE_FUSION_CALLING_WF.out.cohortwideFusionsFile)
        
        //// Neoepitope prediction
        NEOANTIGEN_PREDICTION_WF(
            IN_SILICO_TRANSCRIPT_VALIDATION_WF.out.fusInspectorTsv,
            IN_SILICO_TRANSCRIPT_VALIDATION_WF.out.filteredAgfusionOutdir,
            AGGREGATE_FUSION_CALLING_WF.out.normFilteredFusionsCh,
            RECURRENT_FUSION_FILTERING_WF.out.cohortRecurrentFusionsCh,
            HLA_TYPING_WITH_FALLBACK_WF.out.sampleSpecificHLAsTsv,
            params.metaDataDir
        )
    }

    // Completion handler
    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        //workDir.resolve("stage-${workflow.sessionId}").deleteDir()
    }
}
