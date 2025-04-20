// preprocess files for HLA-HD typing using SAMPICARD (for BAMs) or BOWTIE2 (for FASTQs)
process FISH_HLA_READS_SAMPBOWT {
    maxForks 2
    publishDir "${params.output_dir}/${sampleName}/HLA-fishing-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name HLA-fishing -v ${params.bowtieIndex_db}:/work/libs/hla -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(inputFiles)

    output:
        tuple val(sampleName), path("*_R?.fastq"), emit: fished_files

    script:
    """
    # Initialize variables
    input_file1=${inputFiles[0]}  # First file in the list will be our main input
    input_file2=${inputFiles[1]}
    echo "Processing input file: \${input_file1}"
    echo "Number of cores to use: ${params.num_cores}"
    
    case "\${input_file1}" in
        *.fastq|*.fq)
            # Uncompressed fastq files
            echo "Processing already decompressed fastq files..."
            
            # map fastq files
            if bash /work/scripts/fish-bowtie2-nf.sh ${sampleName} "\${input_file1}" "\${input_file2}" "${params.num_cores}" "/work/libs/hla/hla_gen"; then
                echo "Fished fq files created."
            else
                echo "HLA fishing did not work."
                exit 1
            fi
            ;;
            
        *.fastq.gz|*.fq.gz)
            # Compressed fastq files
            gunzip -c "\${input_file1}" > "${sampleName}_R1.fastq"
            gunzip -c "\${input_file2}" > "${sampleName}_R2.fastq"

            # map fastq files
            if bash /work/scripts/fish-bowtie2-nf.sh ${sampleName} "${sampleName}_R1.fastq" "${sampleName}_R2.fastq" "${params.num_cores}" "/work/libs/hla/hla_gen"; then
                echo "Fished fq files created."
            else
                echo "HLA fishing did not work."
                exit 1
            fi
            ;;
            
        *.bam)
            # BAM file processing
            echo "Processing bam file: \${input_file1}"
            if bash /work/scripts/fish-sampicard-nf.sh "${sampleName}" "\${input_file1}"; then
                echo "File preprocessing has finished running on ${sampleName}."
            else
                echo "BAM preprocessing failed"
                exit 1
            fi
            ;;
            
        *)
            echo "Unsupported file format: \${input_file1} & \${input_file2}"
            exit 1
            ;;
    esac

    """
    stub:
    """
    touch test-stub_R1.fastq test-stub_R2.fastq
    """
}