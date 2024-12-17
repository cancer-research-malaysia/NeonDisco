// Run HLA typing module
process TYPE_HLA_ALLELES {
    publishDir "${params.output_dir}/${sampleName}-HLA", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__hlahd}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name hla-typing -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(inputFiles)
        val(numCores)

    output:
        tuple val(sampleName), path("**/*_final.result.txt"), emit: hla_combined_result
        path "*-HLAHD-complete.txt"

    script:
    """
    # Initialize variables
    input_file1=${inputFiles[0]}  # First file in the list will be our main input
    input_file2=${inputFiles[1]}
    echo "Processing input file: \${input_file1}"
    echo "Number of cores to use: ${numCores}"
    
    case "\${input_file1}" in
        *.fastq|*.fq)
            # Uncompressed fastq files
            cp "\${input_file1}" "${sampleName}_R1.fastq"
            cp "\${input_file2}" "${sampleName}_R2.fastq"
            read1_file="${sampleName}_R1.fastq"
            read2_file="${sampleName}_R2.fastq"
            ;;
            
        *.fastq.gz|*.fq.gz)
            # Compressed fastq files
            gunzip -c "\${input_file1}" > "${sampleName}_R1.fastq"
            gunzip -c "\${input_file2}" > "${sampleName}_R2.fastq"
            read1_file="${sampleName}_R1.fastq"
            read2_file="${sampleName}_R2.fastq"
            ;;
            
        *.bam)
            # BAM file processing
            echo "Processing bam file: \${input_file1}"
            if bash /work/scripts/hlahd-bam-preprocess-nf.sh "${sampleName}" "\${input_file1}"; then
                echo "File preprocessing has finished running on ${sampleName}."
                
                # assign output fastq to variables
                read1_file=\$(find . -maxdepth 1 -type f -name "*_Bam2Fq_R1.fastq")
                read2_file=\$(find . -maxdepth 1 -type f -name "*_Bam2Fq_R2.fastq")
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

    echo "Running HLA-HD..."
    # Running HLA-HD on read files
    if bash /work/scripts/hlahd-nf.sh "${sampleName}" "\${read1_file}" "\${read2_file}" ${numCores}; then
        echo "HLA-HD completed its run. Check outputs for run status."
        echo "YAY!" > ${sampleName}-HLAHD-complete.txt
    else
        echo "HLA-HD failed to complete."
        exit 1
    fi
    """
    stub:
    """
    touch ${sampleName}-HLAHD-complete.txt
    echo "stub run finished!" > ${sampleName}-HLAHD-complete.txt
    """
}
