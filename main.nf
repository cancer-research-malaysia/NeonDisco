#!/usr/bin/env nextflow

// All of the default parameters are being set in `nextflow.config`
// Import sub-workflows
include { listInputFiles } from './modules/listInputFiles'
//include { generateBAMPaths } from './modules/generateBAMPaths'
//include { callVariants } from './modules/callVariants'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run main.nf <--OPTION NAME> <ARGUMENT>

Required Arguments:
---------------
  Input paths:
    --input_dir                     Path to base directory where directories of datasets containing BAM files [MANDATORY]
    --genome_fasta                Path to the reference genome fasta file [MANDATORY]
    --fpscore_matrix              Globbed path to base directory (path/to/dir/*.ext) where the merged footprint score matrix files of all studied motifs are located (saved as parquet files) [MANDATORY]

  Input manifest:
    --dataset_id_list             Path to a list of dataset IDs to work on. Required if [--run_mode] is set to `subset`

  Output path:
    --output_dir                  Directory path for output VCF files [MANDATORY]

---------------
    --help                        Print this help message and exit
    
    """.stripIndent()
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
    } else if ( !params.input_dir ) {
        log.error "The input_dir parameter is not set. Please provide a valid input directory path."
        exit 1
    } else {
        log.info "Initializing workflow..."

        // Preparing the input data
        log.info "Collating input files..."

        def inputFiles = []
        
        log.info "Checking whether the directory of input files is valid..."

        // Check if input file dir exists or not
        if ( file(params.input_dir).exists() == true && file(params.input_dir).isDirectory() == true ) {
            log.info "The input file directory provided is valid. Getting a list of the input files..."
            file(params.input_dir).eachFileMatch(~/.*\.fq$/) { file ->
                                                                        inputFiles << file
                                                                        }
            log.info "Total number of input files found: ${inputFiles.size}"
            if (inputFiles.size() == 0){
                log.error "No input files (fastqs) were found in the provided directory. Please provide a valid directory path containing raw input files with the .fastq extension."
                exit 1
            } else {
                log.info "Execution of workflow will proceed..."
            }
        }
        else {
            log.error "The input file directory does not exist. Please provide a valid directory path."
            exit 1
        }
    
        // Set up a channel to grab all the matrix files in the input folder
        in_files = Channel.fromPath(inputFiles).view()
        // Extract the prefix from the input files and return a tuple of the file and the prefix
        //motifMatrix_ch = in_matrix.map{ file -> [file.baseName.replaceAll("_tfbs_merged_matrix-full", ""), file] }//.view()
        // Extract the TF footprint regions (TFBS) from the input fps matrix files as bed files
        //bedFilesList = extractTFBSBeds(motifMatrix_ch).toList()//.view()


        // Check if run_mode is set
        // if (params.run_mode == "subset"){

        //     // Check if dataset_id_list is provided
        //     if (params.dataset_id_list == false){
        //         log.error "The [--dataset_id_list] option is required to run the workflow if [--run_mode] is set to 'subset'."
        //         exit 1
        //     } else {

        //         // Extract bam directory paths
        //         datasetIDs_ch = Channel.fromPath(params.dataset_id_list).splitText().map { it.trim() }//.view()
        //         datasetIDPaths_ch = datasetIDs_ch.map { id -> [id, file("${params.bam_dir}/${id}")] }//.view()
        //     }
        // }
        // else if ( params.run_mode == "all" ) {

        //     // Extract all unique dataset IDs in the input bam folder and the path to the dataset ID
        //     datasetIDPaths_ch = Channel.fromPath("${params.bam_dir}/*", type: 'dir').map { dir -> [dir.name, dir] }//.view()

        // }

        // // Set up a channel to grab all the bam files for each dataset ID
        // datasetIDBams_ch = datasetIDPaths_ch.map { id, path -> 
        //                                     def bams = files("${path}/**.bam")
        //                                     return [id, bams]
        //                                     }//.view()

        // // Generate a channel for the bam list of a dataset ID
        // bamPaths_ch = generateBAMPaths(datasetIDBams_ch)//.view()

        // // now we can run the variant-calling process
        // rawVCFs_ch = callVariants(bamPaths_ch, bedFilesList)//.view()
    }
}
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}