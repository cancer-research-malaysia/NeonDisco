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
include { TRANSLATE_IN_SILICO_AGFUSION } from './modules/translate_in_silico_AGFUSION.nf'
include { COLLECT_COHORTWIDE_PROTEIN_CODING_FUSIONS_PYENV } from './modules/collect_cohortwide_protein_coding_fusions_PYENV.nf'
include { VALIDATE_IN_SILICO_FUSIONINSPECTOR } from './modules/validate_in_silico_FUSIONINSPECTOR.nf'
include { COLLECT_COHORTWIDE_VALIDATED_FUSIONS_PYENV } from './modules/collect_cohortwide_validated_fusions_PYENV.nf'
include { GET_COHORTWIDE_RECURRENT_VALIDATED_FUSIONS_PYENV } from './modules/get_cohortwide_recurrent_validated_fusions_PYENV.nf'

////// FUSION NEOEPITOPE PREDICTION MODULES //////////
include { KEEP_VALIDATED_FUSIONS_PYENV } from './modules/keep_validated_fusions_PYENV.nf'
include { FILTER_SAMPLE_LEVEL_VALIDATED_FUSIONS_FOR_RECURRENT_PYENV } from './modules/filter_sample_level_validated_fusions_for_recurrent_PYENV.nf'
include { PREDICT_NEOPEPTIDES_SAMPLE_LEVEL_HLAS_PVACFUSE } from './modules/predict_neopeptides_sample_level_hlas_PVACFUSE.nf'
include { PREDICT_NEOPEPTIDES_COHORT_LEVEL_HLAS_PVACFUSE } from './modules/predict_neopeptides_cohort_level_hlas_PVACFUSE.nf'
include { COLLECT_COHORTWIDE_FUSION_NEOPEPTIDES_SAMPHLA_PYENV } from './modules/collect_cohortwide_fusion_neopeptides_sampHLA_PYENV.nf'
include { COLLECT_COHORTWIDE_FUSION_NEOPEPTIDES_COHOHLA_PYENV } from './modules/collect_cohortwide_fusion_neopeptides_cohoHLA_PYENV.nf'

////// HLA TYPING MODULES //////////
include { TYPE_HLAS_WITH_FALLBACK_ARCASHLA } from './modules/type_hlas_with_fallback_ARCASHLA.nf'
include { REFORMAT_AND_COLLATE_HLA_RESULTS_PYENV } from './modules/reformat_and_collate_hla_results_PYENV.nf'
include { FILTER_HLA_ALLOTYPES_FREQ_PYENV } from './modules/filter_hla_allotypes_freq_PYENV.nf'

////// RNA EDITING MODULES (NEW) //////////
// include { CALL_RNA_EDITS_JACUSA2 } from './modules/call_rna_edits_JACUSA2.nf'
// include { ANNOTATE_RNA_EDITS_ANNOVAR } from './modules/annotate_rna_edits_ANNOVAR.nf'
// include { FILTER_RNA_EDITS_PYENV } from './modules/filter_rna_edits_PYENV.nf'

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to print help message
def helpMessage() {
    log.info"""
Usage:

nextflow run neondisco-main.nf -profile <EXECUTOR,MODE[,RESOURCE]> <--OPTION NAME> <ARGUMENT>

Profile Examples:
----------------
    -profile local,sampleHLANeoPredMode              # Run locally and run sample-level HLA neoantigen prediction
    -profile local,sharedHLANeoPredMode              # Run locally and run shared, cohort-level HLA neoantigen prediction
    -profile local,dualNeoPredMode                   # Run locally in both modes
    -profile awsbatch,sampleHLANeoPredMode           # Run on AWS Batch in sample-level neoantigen mode
    -profile awsbatch,sharedHLANeoPredMode              # Run on AWS Batch in shared neoantigen mode
    -profile awsbatch,dualNeoPredMode                # Run on AWS Batch in both modes

Required Arguments:
---------------
    -c <configFile>                 Path to the config file. [REQUIRED]
    -profile                        Comma-separated list of profiles [REQUIRED]
                                        EXECUTOR:   local | awsbatch
                                        MODE:       sampleHLANeoPredMode | sharedHLANeoPredMode | dualNeoPredMode
    --inputDir                      Path to local directory containing BAM/FASTQ input files [REQUIRED if manifestPath not provided]
    --manifestPath                  Path to tab-delimited manifest file [REQUIRED if inputDir not provided]
                                        – must contain sample ID and read1 and read2 local filepaths or remote s3 filepaths
                                        – must also contain sampleType column with values 'Tumor' or 'Normal'
    --inputSource                   Input source type: <local> for local files, <s3> for S3 files [REQUIRED]
                                        – if inputSource is set to <s3>, --inputDir cannot be used and --manifestPath must be provided
                                        - if inputSource is set to <local> -profile 'awsbatch' is disallowed

Optional Arguments:
---------------
    --outputDir                     Directory path for output; can be s3 URIs [DEFAULT: ./outputs]
    --trimReads                     Skip read trimming on FASTQ input [DEFAULT: true]
    --hlaTypingOnly                 Exclusively run HLA typing subworkflow [DEFAULT: false]
    --includeNeoPred                Run neopeptide prediction subworkflow [DEFAULT: true]
    --deleteIntMedFiles             Delete intermediate files right after they are not needed [DEFAULT: false]
    --deleteStagedFiles             Delete staged files after processing [DEFAULT: true if inputSource is 's3']
    --sampleHLANeoPred              Run neopeptide prediction using sample-level HLAs [DEFAULT: varies by mode]
    --sharedHLANeoPred              Run neopeptide prediction using cohort-level HLAs [DEFAULT: varies by mode]
    --recurrentFusionsNeoPredOnly   Predict neopeptides for recurrent fusions only [DEFAULT: varies by mode]
                                        – true:  Only predict neopeptides for shared recurrent fusions
                                        – false: Predict neopeptides for all validated fusions regardless of recurrence
    --recurrenceThreshold           Threshold for recurrent fusions; only applies if --recurrentFusionsNeoPredOnly is set to true
                                        - [DEFAULT: 0.005] (0.5% recurrence)

Discovery Module Options (All enabled by default - use flags to disable):
--------------------------
    --disableFusionDiscovery        Disable fusion discovery module [DEFAULT: false]
        --disableArriba                 Disable Arriba fusion caller [DEFAULT: false]
        --disableFusioncatcher          Disable FusionCatcher fusion caller [DEFAULT: false]
        --disableStarFusion             Disable STAR-Fusion caller [DEFAULT: false]
    --disableRnaEditingDiscovery    Disable RNA editing discovery module [DEFAULT: true - not yet implemented]
        --rnaEditingMinDepth            Minimum depth for RNA editing calls [DEFAULT: 10]
        --rnaEditingMinAltFreq          Minimum alternate allele frequency for RNA editing [DEFAULT: 0.05]

    --help                          Print this help message and exit
    """.stripIndent()
}

// Function to validate profiles
def validateProfiles() {
    if (!workflow.profile) {
        log.error "No profile specified. Please specify -profile <EXECUTOR,MODE[,RESOURCE]>"
        log.error "Examples: -profile local,sampleHLANeoPredMode or -profile awsbatch,sharedHLANeoPredMode"
        return false
    }
    
    def profiles = workflow.profile.split(',').collect { profile -> profile.trim() }
    def executors = ['local', 'awsbatch']
    def modes = ['sampleHLANeoPredMode', 'sharedHLANeoPredMode', 'dualNeoPredMode']
    
    // Check for required executor and mode
    def hasExecutor = profiles.any { profile -> executors.contains(profile) }
    def hasMode = profiles.any { profile -> modes.contains(profile) }
    
    if (!hasExecutor) {
        log.error "You must specify an executor profile: ${executors.join(', ')}"
        return false
    }
    
    if (!hasMode) {
        if (params.includeNeoPred) {
            log.error "You must specify a neoantigen prediction mode profile: ${modes.join(', ')}, as --includeNeoPred is set to True"
            return false
        } else {
            log.warn "No neoantigen prediction mode profile specified, AND includeNeoPred is set to false. Proceeding without neoantigen prediction."
            return true
        }
    }
    
    // Check for multiple executors (still not allowed)
    def executorCount = profiles.count { profile -> executors.contains(profile) }
    if (executorCount > 1) {
        log.error "Only one executor profile allowed. Found: ${profiles.findAll { profile -> executors.contains(profile) }.join(', ')}"
        return false
    }
    
    // Check for mode combinations - now only allow one mode since we have dualNeoPredMode
    def modeCount = profiles.count { profile -> modes.contains(profile) }
    if (modeCount > 1) {
        log.error "Only one neoantigen prediction mode allowed. Found: ${profiles.findAll { profile -> modes.contains(profile) }.join(', ')}"
        log.error "Use dualNeoPredMode to run both sample-level and shared HLA neoantigen predictions"
        return false
    } else if (modeCount == 0 && params.includeNeoPred) {
        log.error "You must specify a neoantigen prediction mode profile: ${modes.join(', ')}, as --includeNeoPred is set to True"
        return false
    } else if (modeCount == 0 && !params.includeNeoPred) {
        log.warn "No neoantigen prediction mode profile specified, AND includeNeoPred is set to false. Proceeding without neoantigen prediction."
        return true
    }
    
    // Log which mode is being used
    if (profiles.contains('dualNeoPredMode')) {
        log.info "Running in dual mode: both sample-level and shared neoantigen discovery enabled"
    } else if (profiles.contains('sampleHLANeoPredMode')) {
        log.info "Running in sample-level HLA neoantigen prediction mode"
    } else if (profiles.contains('sharedHLANeoPredMode')) {
        log.info "Running in shared HLA neoantigen prediction mode"
    }
    
    // Validate unknown profiles
    def validProfiles = executors + modes
    def unknownProfiles = profiles.findAll { profile -> !validProfiles.contains(profile) }
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
    def bamCh = channel.empty()
    def fastqCh = channel.empty()

    // Process BAM files if they exist
    if (bamFiles) {
        log.info "[STATUS] Found BAM files in ${dirPath}"
        bamCh = channel.fromPath("${dirPath}/*.{bam,bai}")
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
        fastqCh = channel.fromFilePairs("${dirPath}/*{R,r}{1,2}*.{fastq,fq}{,.gz}", checkIfExists: true)
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
    def inputCh = channel
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
        .filter { rowData -> rowData != null }

    return inputCh
}

// Function to branch input channel by sample type
def branchInputChannelBySampleType(inputCh) {
    def branchedCh = inputCh.branch { sample ->
        tumor: sample[2] == 'Tumor'
        normal: sample[2] == 'Normal'
    }
    
    // Convert back to the original format (sampleName, [read1, read2])
    def tumorCh = branchedCh.tumor.map { sampleName, reads, _sampleType -> tuple(sampleName, reads) }
    def normalCh = branchedCh.normal.map { sampleName, reads, _sampleType -> tuple(sampleName, reads) }
    
    return [tumorCh, normalCh]
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// MAIN WORKFLOWS ////////////////////////////////////////////

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
        alignedBamCh = ALIGN_READS_STAR_GENERAL.out.final_bam
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

// ============================================================================
// REFACTORED MODULAR FUSION DISCOVERY WORKFLOW
// ============================================================================

workflow FUSION_DISCOVERY_MODULE {
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

        // Initialize empty channels for each caller
        def arribaResults = channel.empty()
        def fusioncatcherResults = channel.empty()
        def starfusionResults = channel.empty()

        // Conditionally run each fusion caller
        if (!params.disableArriba) {
            log.info "[FUSION DISCOVERY] Arriba enabled"
            CALL_FUSIONS_ARRIBA(
                ALIGN_READS_STAR_ARRIBA(filtFastqsCh, starIndex).aligned_bam, 
                arribaDB
            )
            arribaResults = CALL_FUSIONS_ARRIBA.out.arriba_fusion_tuple
        } else {
            log.info "[FUSION DISCOVERY] Arriba disabled"
        }

        if (!params.disableFusioncatcher) {
            log.info "[FUSION DISCOVERY] FusionCatcher enabled"
            CALL_FUSIONS_FUSIONCATCHER(filtFastqsCh, fuscatDB)
            fusioncatcherResults = CALL_FUSIONS_FUSIONCATCHER.out.fuscat_fusion_tuple
        } else {
            log.info "[FUSION DISCOVERY] FusionCatcher disabled"
        }

        if (!params.disableStarFusion) {
            log.info "[FUSION DISCOVERY] STAR-Fusion enabled"
            CALL_FUSIONS_STARFUSION(filtFastqsCh, ctatDB)
            starfusionResults = CALL_FUSIONS_STARFUSION.out.starfus_fusion_tuple
        } else {
            log.info "[FUSION DISCOVERY] STAR-Fusion disabled"
        }

        // Intelligently combine results based on what's enabled
        def combinedFTFilesCh = channel.empty()
        
        // Count enabled callers
        def enabledCount = [!params.disableArriba, !params.disableFusioncatcher, !params.disableStarFusion].count { caller -> caller }
        
        if (enabledCount == 0) {
            log.error "[FUSION DISCOVERY] All fusion callers disabled! At least one must be enabled."
            exit 1
        }
        
        if (!params.disableArriba && !params.disableFusioncatcher && !params.disableStarFusion) {
            combinedFTFilesCh = arribaResults
                .join(fusioncatcherResults)
                .join(starfusionResults)
        } else if (!params.disableArriba && !params.disableFusioncatcher) {
            combinedFTFilesCh = arribaResults
                .join(fusioncatcherResults)
                .map { sampleName, arribaFile, fuscatFile -> 
                    tuple(sampleName, arribaFile, fuscatFile, null) 
                }
        } else if (!params.disableArriba && !params.disableStarFusion) {
            combinedFTFilesCh = arribaResults
                .join(starfusionResults)
                .map { sampleName, arribaFile, starfusFile -> 
                    tuple(sampleName, arribaFile, null, starfusFile) 
                }
        } else if (!params.disableFusioncatcher && !params.disableStarFusion) {
            combinedFTFilesCh = fusioncatcherResults
                .join(starfusionResults)
                .map { sampleName, fuscatFile, starfusFile -> 
                    tuple(sampleName, null, fuscatFile, starfusFile) 
                }
        } else if (!params.disableArriba) {
            combinedFTFilesCh = arribaResults
                .map { sampleName, arribaFile -> 
                    tuple(sampleName, arribaFile, null, null) 
                }
        } else if (!params.disableFusioncatcher) {
            combinedFTFilesCh = fusioncatcherResults
                .map { sampleName, fuscatFile -> 
                    tuple(sampleName, null, fuscatFile, null) 
                }
        } else if (!params.disableStarFusion) {
            combinedFTFilesCh = starfusionResults
                .map { sampleName, starfusFile -> 
                    tuple(sampleName, null, null, starfusFile) 
                }
        }

        // Collate and filter fusions
        COLLATE_FUSIONS_PYENV(combinedFTFilesCh)
        collatedFusionsParquet = COLLATE_FUSIONS_PYENV.out.collatedFTParquet

        // Branch 1: Wrangle raw fusions for unfiltered collection
        WRANGLE_RAW_FUSIONS_PYENV(collatedFusionsParquet)
        COLLECT_COHORTWIDE_UNFILTERED_FUSIONS_PYENV(
            WRANGLE_RAW_FUSIONS_PYENV.out.wrangledUnfilteredFusionsTsv
                .collect { _sampleName, tsv -> tsv }
        )

        // Branch 2: Filter fusions for downstream processing
        FILTER_FUSIONS_PYENV(collatedFusionsParquet, metaDataDir)

        // Collect cohort-wide normfiltered fusions
        COLLECT_COHORTWIDE_NORMFILTERED_FUSIONS_PYENV(
            FILTER_FUSIONS_PYENV.out.filteredFusions
                .collect { _sampleName, tsv -> tsv }
        )

        // Format output to standardized schema for unified downstream processing
        standardizedOutput = FILTER_FUSIONS_PYENV.out.filteredFusions
            .map { sampleName, fusionFile -> 
                tuple(sampleName, 'fusion', fusionFile, null)
            }

    emit:
        // Standardized output for unified discovery processing
        discoveries = standardizedOutput
        
        // Module-specific outputs for fusion-specific downstream processing
        normFilteredFusionsCh = FILTER_FUSIONS_PYENV.out.filteredFusions
        uniqueFiltFusionPairsForFusInsCh = FILTER_FUSIONS_PYENV.out.uniqueFiltFusionPairsForFusIns
        cohortwideNormfilteredFusionsFile = COLLECT_COHORTWIDE_NORMFILTERED_FUSIONS_PYENV.out.cohortwideFusionsFile
}

// ============================================================================
// RNA EDITING DISCOVERY MODULE (NEW)
// ============================================================================

// workflow RNA_EDITING_DISCOVERY_MODULE {
//     take:
//         alignedBamCh        // [sampleName, bam, bamIdx]
//         referenceGenome
//         annovarDB
//         metaDataDir
    
//     main:
//         log.info "[RNA EDITING DISCOVERY] Starting RNA editing discovery module"
        
//         // Call RNA editing sites
//         CALL_RNA_EDITS_JACUSA2(alignedBamCh, referenceGenome)
        
//         // Annotate editing sites
//         ANNOTATE_RNA_EDITS_ANNOVAR(
//             CALL_RNA_EDITS_JACUSA2.out.editing_sites,
//             annovarDB
//         )
        
//         // Filter editing sites
//         FILTER_RNA_EDITS_PYENV(
//             ANNOTATE_RNA_EDITS_ANNOVAR.out.annotated_edits,
//             metaDataDir
//         )

//         // Format output to standardized schema
//         standardizedOutput = FILTER_RNA_EDITS_PYENV.out.filteredEdits
//             .map { sampleName, editFile -> 
//                 tuple(sampleName, 'rna_editing', editFile, null)
//             }

//     emit:
//         discoveries = standardizedOutput
//         filteredEdits = FILTER_RNA_EDITS_PYENV.out.filteredEdits
// }

// ============================================================================
// UNIFIED DISCOVERY ORCHESTRATOR
// ============================================================================

workflow UNIFIED_DISCOVERY_MODULE {
    take:
        alignedBamCh
        starIndex
        _referenceGenome
        arribaDB
        fuscatDB
        ctatDB
        _annovarDB
        metaDataDir
    
    main:
        log.info "=== STARTING UNIFIED DISCOVERY MODULE ==="
        
        // Initialize channel for all discoveries
        def allDiscoveries = channel.empty()
        def fusionModuleOutputs = null
        def _rnaEditingModuleOutputs = null

        // Conditionally run fusion discovery
        if (!params.disableFusionDiscovery) {
            log.info "=== FUSION DISCOVERY MODULE ENABLED ==="
            FUSION_DISCOVERY_MODULE(
                alignedBamCh,
                starIndex,
                arribaDB,
                fuscatDB,
                ctatDB,
                metaDataDir
            )
            allDiscoveries = allDiscoveries.mix(FUSION_DISCOVERY_MODULE.out.discoveries)
            fusionModuleOutputs = FUSION_DISCOVERY_MODULE.out
        } else {
            log.info "=== FUSION DISCOVERY MODULE DISABLED ==="
        }

        // Conditionally run RNA editing discovery
        // if (!params.disableRnaEditingDiscovery) {
        //     log.info "=== RNA EDITING DISCOVERY MODULE ENABLED ==="
        //     RNA_EDITING_DISCOVERY_MODULE(
        //         alignedBamCh,
        //         referenceGenome,
        //         annovarDB,
        //         metaDataDir
        //     )
        //     allDiscoveries = allDiscoveries.mix(RNA_EDITING_DISCOVERY_MODULE.out.discoveries)
        //     rnaEditingModuleOutputs = RNA_EDITING_DISCOVERY_MODULE.out
        // } else {
        //     log.info "=== RNA EDITING DISCOVERY MODULE DISABLED ==="
        // }

        // Branch discoveries by type for type-specific processing
        def branchedDiscoveries = allDiscoveries.branch { sample ->
            fusion: sample[1] == 'fusion'
            rna_editing: sample[1] == 'rna_editing'
            splicing: sample[1] == 'splicing'
        }

    emit:
        // Unified output for generic downstream processing
        allDiscoveries = allDiscoveries
        
        // Type-specific outputs
        fusionDiscoveries = branchedDiscoveries.fusion
        rnaEditingDiscoveries = branchedDiscoveries.rna_editing
        
        // Module-specific outputs (for specialized downstream processing)
        fusionNormFilteredFusions = !params.disableFusionDiscovery ? 
            fusionModuleOutputs.normFilteredFusionsCh : channel.empty()
        fusionUniquePairs = !params.disableFusionDiscovery ? 
            fusionModuleOutputs.uniqueFiltFusionPairsForFusInsCh : channel.empty()
        fusionCohortwideFile = !params.disableFusionDiscovery ?
            fusionModuleOutputs.cohortwideNormfilteredFusionsFile : channel.empty()
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
///////////////// IN-SILICO TRANSCRIPT VALIDATION WORKFLOW /////////////////////////
workflow PROTEIN_CODING_PREDICTION {
    take:
        normFilteredFusionsCh
        cohortwideNormfilteredFusionsFile
    main:
        // Translate filtered fusions in silico with AGFusion   
        TRANSLATE_IN_SILICO_AGFUSION(normFilteredFusionsCh)
        
        // Collect protein-coding fusion output manifests and concat into cohortwide file
        COLLECT_COHORTWIDE_PROTEIN_CODING_FUSIONS_PYENV(
            TRANSLATE_IN_SILICO_AGFUSION.out.protein_coding_fusions_manifest.collect { _sampleName, manifest -> manifest },
            cohortwideNormfilteredFusionsFile
        )
    emit:
        filteredAgfusionOutdir = TRANSLATE_IN_SILICO_AGFUSION.out.filtered_agfusion_outdir
        proteinCodingManifest = TRANSLATE_IN_SILICO_AGFUSION.out.protein_coding_fusions_manifest
}

workflow FUSION_INSPECTOR_VALIDATION {
    take:
        filteredAgfusionOutdir
        uniqueFiltFusionPairsForFusInsCh
        trimmedCh
        ctatDB
    main:
        // Join channels and fix the tuple structure
        joinedInputs = uniqueFiltFusionPairsForFusInsCh
            .join(trimmedCh, by: 0)
            .join(filteredAgfusionOutdir, by: 0)
            .map { sampleName, uniqueFiltPairs, trimmedReads, agfusionDir -> 
                tuple(sampleName, agfusionDir, uniqueFiltPairs, trimmedReads)
            }
        
        VALIDATE_IN_SILICO_FUSIONINSPECTOR(joinedInputs, ctatDB)

    emit:
        fusInspectorTsv = VALIDATE_IN_SILICO_FUSIONINSPECTOR.out.fusInspectorTsv
}

workflow IN_SILICO_FUSION_VALIDATION_WF {
    take:
        normFilteredFusionsCh
        cohortwideNormfilteredFusionsFile
        uniqueFiltFusionPairsForFusInsCh
        trimmedCh
        ctatDB
    main:
        PROTEIN_CODING_PREDICTION(normFilteredFusionsCh, cohortwideNormfilteredFusionsFile)
        
        FUSION_INSPECTOR_VALIDATION(
            PROTEIN_CODING_PREDICTION.out.filteredAgfusionOutdir,
            uniqueFiltFusionPairsForFusInsCh,
            trimmedCh,
            ctatDB
        )

    emit:
        fusInspectorTsv = FUSION_INSPECTOR_VALIDATION.out.fusInspectorTsv
        filteredAgfusionOutdir = PROTEIN_CODING_PREDICTION.out.filteredAgfusionOutdir
        proteinCodingManifest = PROTEIN_CODING_PREDICTION.out.proteinCodingManifest     
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
workflow GET_COHORTWIDE_FI_VALIDATED_FUSIONS {
    take:
        validatedFusionsTsvs
    main:
        // get cohort-wide validated fusions
        COLLECT_COHORTWIDE_VALIDATED_FUSIONS_PYENV(validatedFusionsTsvs)
    emit:
        cohortValidatedFusions = COLLECT_COHORTWIDE_VALIDATED_FUSIONS_PYENV.out.cohortwideValidatedFusionsFile
}

// Filter validated fusions
workflow VALIDATED_FUSION_FILTERING_WF {
    take:
        fusInspectorTsv
        filteredAgfusionOutdir
        normFilteredFusionsCh
        proteinCodingManifest
    main:
        // Join ALL inputs by sampleName, including the manifest
        joinedInputs = fusInspectorTsv
            .join(filteredAgfusionOutdir, by: 0)
            .join(normFilteredFusionsCh, by: 0)
            .join(proteinCodingManifest, by: 0)
            .map { sampleName, fusInspectorFile, agfusionDir, filteredFusions, proteinManifest -> 
                tuple(sampleName, fusInspectorFile, agfusionDir, filteredFusions, proteinManifest) 
            }
        
        // Preprocess agfusion output for neoepitope prediction
        KEEP_VALIDATED_FUSIONS_PYENV(joinedInputs)

        // Join the validatedFusions and validatedAgfusionDir channels by sampleName
        joinedValidatedFusionsDat = KEEP_VALIDATED_FUSIONS_PYENV.out.validatedFusions
            .join(KEEP_VALIDATED_FUSIONS_PYENV.out.validatedAgfusionDir, by: 0)
            .map { sampleName, validatedFusionsFile, validatedDir -> 
                tuple(sampleName, validatedFusionsFile, validatedDir) 
            }
        
        // Collect cohort-wide validated fusions
        GET_COHORTWIDE_FI_VALIDATED_FUSIONS(
            KEEP_VALIDATED_FUSIONS_PYENV.out.validatedFusions.collect { _sampleName, tsv -> tsv }
        )

    emit:
        validatedFusSampleData = joinedValidatedFusionsDat
        cohortValidatedFusions = GET_COHORTWIDE_FI_VALIDATED_FUSIONS.out.cohortValidatedFusions
}

/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////COHORTWIDE VALIDATED FUSION RECURRENCE WORKFLOW//////////////////

workflow COHORT_VALIDATED_FUSION_RECURRENCE_COMPILATION_WF {
    take:
        cohortValidatedFusions
        totalSampleCount
    main:
        // Get cohort recurrent fusions
        GET_COHORTWIDE_RECURRENT_VALIDATED_FUSIONS_PYENV(cohortValidatedFusions, totalSampleCount)
    emit:
        cohortRecurrentFusionsCh = GET_COHORTWIDE_RECURRENT_VALIDATED_FUSIONS_PYENV.out.cohortRecurrentFusionTsv
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////// NEOANTIGEN PREDICTION SUBWORKFLOWS ///////////////////////
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
    emit:
        predictedSampleNeopeptides = PREDICT_NEOPEPTIDES_SAMPLE_LEVEL_HLAS_PVACFUSE.out.predictedSampleSpecificNeopeptides
}

workflow GET_COHORTWIDE_SAMPHLA_NEOPEPTIDES {
    take:
        sampleLevelHlaPvacFuseOutTsvs
    main:
        // get cohort-wide fusion neopeptides from sample HLA pvacfuse output
        COLLECT_COHORTWIDE_FUSION_NEOPEPTIDES_SAMPHLA_PYENV(sampleLevelHlaPvacFuseOutTsvs)
    emit:
        sampleLevelHlaNeopeptides = COLLECT_COHORTWIDE_FUSION_NEOPEPTIDES_SAMPHLA_PYENV.out.cohortwideSampleHLAFusNeopeptides
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
    emit:
        predictedCohortNeopeptides = PREDICT_NEOPEPTIDES_COHORT_LEVEL_HLAS_PVACFUSE.out.predictedCohortNeopeptides
}
workflow GET_COHORTWIDE_COHORTHLA_NEOPEPTIDES {
    take:
        cohortLevelHlaPvacFuseOutTsvs
    main:
        // get cohort-wide fusion neopeptides from cohort HLA pvacfuse output
        COLLECT_COHORTWIDE_FUSION_NEOPEPTIDES_COHOHLA_PYENV(cohortLevelHlaPvacFuseOutTsvs)
    emit:
        cohortLevelHlaNeopeptides = COLLECT_COHORTWIDE_FUSION_NEOPEPTIDES_COHOHLA_PYENV.out.cohortwideCohortHLAFusNeopeptides
}

// Updated main neopeptide workflow - now just orchestrates the sub-workflows
workflow NEOANTIGEN_PREDICTION_WF {
    take:
        validatedFusSampleData
        cohortRecurrentFusionsCh
        sampleSpecificHLAsTsv
        metaDataDir
    main:

        // Filter sample-specific validated fusions for recurrent ones
        FILTER_SAMPLE_LEVEL_VALIDATED_FUSIONS_FOR_RECURRENT_PYENV(
            validatedFusSampleData.combine(
                cohortRecurrentFusionsCh.map { _label, file -> file }  // Drop the label
            )
        )
        
        // Get the full validated fusions directory of a sample
        recurrentValidatedDir = FILTER_SAMPLE_LEVEL_VALIDATED_FUSIONS_FOR_RECURRENT_PYENV.out.validatedRecurrentAgfusionDir
        
        // Logic based on recurrentFusionsNeoPredOnly parameter
        def finalAgfusionDir = channel.empty()

        if (params.recurrentFusionsNeoPredOnly) {
            // Default mode: Only process recurrent fusions
            finalAgfusionDir = recurrentValidatedDir
                .ifEmpty { 
                    log.warn "No recurrent fusions found in this input cohort. Neoantigen prediction will be skipped."
                    log.info "Consider using [--recurrentFusionsNeoPredOnly false] flag to process all validated fusions instead."
                    channel.empty()
                }   
        } else {
            // Alternative mode: Process all validated fusions
            finalAgfusionDir = validatedFusSampleData
                .map { sampleName, _validatedFusionsFile, validatedDir -> 
                    tuple(sampleName, validatedDir) 
                }
            log.info "Processing all validated fusions [--recurrentFusionsNeoPredOnly false]..."
        }
        
        // Run neoepitope prediction subworkflow
        if (params.includeNeoPred) {
            log.info "Running neoepitope prediction subworkflow..."
            // Run sample-specific if enabled
            if (params.sampleHLANeoPred) {
                SAMPLE_LEVEL_HLA_NEOANTIGENS(finalAgfusionDir, sampleSpecificHLAsTsv, metaDataDir)
                // Collect cohort-wide sample HLA neopeptides
                sampleLevelFusNeopeps = SAMPLE_LEVEL_HLA_NEOANTIGENS.out.predictedSampleNeopeptides
                GET_COHORTWIDE_SAMPHLA_NEOPEPTIDES(sampleLevelFusNeopeps.collect())
            }

            // Run cohort-wide if enabled
            if (params.sharedHLANeoPred) {
                COHORT_LEVEL_HLA_NEOANTIGENS(finalAgfusionDir, sampleSpecificHLAsTsv, metaDataDir)
                // Collect cohort-wide cohort HLA neopeptides
                cohortLevelFusNeopeps = COHORT_LEVEL_HLA_NEOANTIGENS.out.predictedCohortNeopeptides
                GET_COHORTWIDE_COHORTHLA_NEOPEPTIDES(cohortLevelFusNeopeps.collect())
            }

        } else {
            log.info "Skipping neoepitope prediction subworkflow."
        }
}
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// NEONDISCO MAIN WORKFLOW ////////////////////////////////////////////

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
    def profilesList = workflow.profile.split(',').collect { profile -> profile.trim() }
    def isAwsBatch = profilesList.contains('awsbatch')
    def isLocal = profilesList.contains('local')

    if (isLocal && isAwsBatch) {
        log.error "AWS Batch executor is not compatible with local input source."
        log.error "Please use either:"
        log.error "  - Local executor with local or S3 input: -profile local,<MODE> --inputSource local/s3"
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
    def tumorCh = channel.empty()
    def normalCh = channel.empty()
    def totalSampleCountCh = channel.value(0)

    if (params.manifestPath) {
        def inputCh = createInputChannelFromManifest(params.manifestPath)
        log.info "Using manifest file: << ${params.manifestPath} >>"
        
        // Branch the input channel by sample type
        def (branchedTumorCh, branchedNormalCh) = branchInputChannelBySampleType(inputCh)

        // Create a copy of tumor channel for counting (since channels can only be consumed once)
        def tumorChForCount = branchedTumorCh.map { sample -> sample }
        def tumorChForPipeline = branchedTumorCh.map { sample -> sample }
        
        tumorCh = tumorChForPipeline
        normalCh = branchedNormalCh

        // Capture tumor count as a CHANNEL VALUE (not a variable)
        totalSampleCountCh = tumorChForCount
                            .toList()
                            .map { sampleList -> sampleList.size() }
        
        // Log sample counts by type AND capture tumor count for recurrence calculation
        normalCh.count().subscribe { count ->
            log.info "Found ${count} normal sample(s)! ----Normal samples will not be processed in the main NeonDisco pipeline."
        }
        tumorCh.count().subscribe { count ->
            log.info "Found ${count} tumor sample(s)! ----Initializing NeonDisco pipeline..."
            log.info "Using ${count} tumor samples for recurrence frequency calculations."
            log.info ""
        }

    } else {

        if (!validateInputDir(params.inputDir)) {
            exit 1
        }
        
        def inputCh = createInputChannelFromPOSIX(params.inputDir)

        // Create a copy for counting
        def tumorChForCount = inputCh.map { sample -> sample }
        def tumorChForPipeline = inputCh.map { sample -> sample }
        
        tumorCh = tumorChForPipeline
        
        // Capture count as a CHANNEL VALUE
        totalSampleCountCh = tumorChForCount
            .toList()
            .map { sampleList -> sampleList.size() }
        
        totalSampleCountCh.subscribe { count ->
            log.info "Input files are provided as local directory: << ${params.inputDir} >>"
            log.info "Found ${count} sample(s) in directory! ----All samples will be processed as TUMOR samples!"
            log.warn "If you have normal samples mixed in, please provide a manifest file with sampleType column at runtime instead of --inputDir."
            log.info "Using ${count} samples for recurrence frequency calculations."
        }
    }

    
    // Log the key parameters
    log.info "Output directory: << ${params.outputDir} >>"
    log.info "Read trimming: << ${params.trimReads} >>"
    log.info "HLA typing only? : << ${params.hlaTypingOnly} >>"
    log.info "Include neopeptide prediction? : << ${params.includeNeoPred} >>"
    
    // Log discovery module settings
    log.info "=== DISCOVERY MODULE SETTINGS ==="
    log.info "Fusion discovery enabled: << ${!params.disableFusionDiscovery} >>"
    if (!params.disableFusionDiscovery) {
        log.info "  ├─ Arriba: << ${!params.disableArriba} >>"
        log.info "  ├─ FusionCatcher: << ${!params.disableFusioncatcher} >>"
        log.info "  └─ STAR-Fusion: << ${!params.disableStarFusion} >>"
    }
    log.info "RNA editing discovery enabled: << ${!params.disableRnaEditingDiscovery} >>"
    if (!params.disableRnaEditingDiscovery) {
        log.info "  ├─ Minimum depth: << ${params.rnaEditingMinDepth} >>"
        log.info "  └─ Minimum alt frequency: << ${params.rnaEditingMinAltFreq} >>"
    }
    
    // Log the fusion prediction mode
    def mode = params.recurrentFusionsNeoPredOnly ? "Recurrent-only" : "All-validated"
    log.info "Fusion-derived neoantigen prediction input set: << ${mode} fusions >>"
    
    if (params.recurrentFusionsNeoPredOnly) {
        log.info "----Recurrence threshold: << ${params.recurrenceThreshold * 100}% >>"
    } else {
        log.info "----Recurrent fusion–only parameter is disabled; all validated fusions will be predicted for neopeptides."
    }
    log.info "Neopeptide prediction mode [Personalized (using sample-level HLA allotypes)]: << ${params.sampleHLANeoPred} >>"
    log.info "Neopeptide prediction mode [Cohort-based (using cohort-level HLA allotypes)]: << ${params.sharedHLANeoPred} >>"
    
    // Process tumor samples only (normal channel remains unused but available)
    def qcProcInputCh = params.trimReads ? TRIMMING_WF(tumorCh).trimmedCh : tumorCh

    ///////// STAR alignment workflow ////////////
    def alignedBamsCh = GENERAL_READS_ALIGNMENT_WF(qcProcInputCh, params.starIndex).alignedBamCh
    ////////////////////////////////////////////////////////

    // Execute workflow branching based on hlaTypingOnly parameter
    if (params.hlaTypingOnly) {
        log.info "Running HLA typing only subworkflow..."
        // Run only HLA typing
        HLA_TYPING_WITH_FALLBACK_WF(alignedBamsCh)

    } else {
        //// HLA typing
        HLA_TYPING_WITH_FALLBACK_WF(alignedBamsCh)

        //// Unified discovery module - handles both fusion and RNA editing
        UNIFIED_DISCOVERY_MODULE(
            alignedBamsCh,
            params.starIndex,
            params.referenceGenome,
            params.arribaDB,
            params.fuscatDB,
            params.ctatDB,
            params.annovarDB,
            params.metaDataDir
        )

        // Continue with fusion-specific processing if enabled
        if (!params.disableFusionDiscovery) {
            //// AGFUSION coding sequence prediction
            IN_SILICO_FUSION_VALIDATION_WF(
                UNIFIED_DISCOVERY_MODULE.out.fusionNormFilteredFusions,
                UNIFIED_DISCOVERY_MODULE.out.fusionCohortwideFile,
                UNIFIED_DISCOVERY_MODULE.out.fusionUniquePairs,
                qcProcInputCh,
                params.ctatDB
            )
            
            //// Validated fusion filtering
            VALIDATED_FUSION_FILTERING_WF(
                IN_SILICO_FUSION_VALIDATION_WF.out.fusInspectorTsv,
                IN_SILICO_FUSION_VALIDATION_WF.out.filteredAgfusionOutdir,
                UNIFIED_DISCOVERY_MODULE.out.fusionNormFilteredFusions,
                IN_SILICO_FUSION_VALIDATION_WF.out.proteinCodingManifest
            )

            //// Cohortwide validated fusion recurrence
            COHORT_VALIDATED_FUSION_RECURRENCE_COMPILATION_WF(
                VALIDATED_FUSION_FILTERING_WF.out.cohortValidatedFusions,
                totalSampleCountCh
            )

            //// Neoepitope prediction
            NEOANTIGEN_PREDICTION_WF(
                VALIDATED_FUSION_FILTERING_WF.out.validatedFusSampleData,
                COHORT_VALIDATED_FUSION_RECURRENCE_COMPILATION_WF.out.cohortRecurrentFusionsCh,
                HLA_TYPING_WITH_FALLBACK_WF.out.sampleSpecificHLAsTsv,
                params.metaDataDir
            )
        }

        // RNA editing-specific processing would go here if enabled
        if (!params.disableRnaEditingDiscovery) {
            log.info "RNA editing discovery outputs available for downstream processing"
            // Add RNA editing-specific downstream workflows here
        }
    }

    // Completion handler
    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    }
}

