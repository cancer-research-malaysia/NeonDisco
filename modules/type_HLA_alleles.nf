// Run HLA typing module
process TYPE_HLA_ALLELES {
    //publishDir "${params.output_dir}/${sampleName}", mode: 'copy'
    container "${params.container__hlahd}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name hla-typing -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
    tuple val(sampleName), path(inputFiles)
    val(numCores)

    script:

    if (inputFiles instanceof List && inputFiles.size() == 2) {
        """
        # The input has two files; assuming paired-read fastq files
        echo "Processing paired-end fastq files: ${inputFiles[0]} and ${inputFiles[1]}"
        echo "Number of cores to use: ${numCores}"

        read1_file="${inputFiles[0]}"
        read2_file="${inputFiles[1]}"

        # unpack fastq files if needed
        if [[ "\${read1_file}" == *.gz ]]; then
            gunzip -c \${read1_file} > ${sampleName}_R1.fastq
            read1_file=${sampleName}_R1.fastq
            echo \${read1_file}
        fi

        if [[ "\${read2_file}" == *.gz ]]; then
            gunzip -c \${read2_file} > ${sampleName}_R2.fastq
            read2_file=${sampleName}_R2.fastq
            echo \${read2_file}
        fi

        # Running HLA-HD on read files
        if bash /work/scripts/hlahd-nf.sh ${sampleName} \${read1_file} \${read2_file} ${numCores}; then
            echo "HLA-HD completed its run. Check outputs for run status."
            echo "Removing unzipped fastq files from work directory."
            # rm -rf *.fastq
        fi

        """
    } else {
        """
        # Input is solitary bam file
        echo "Processing bam file: ${inputFiles}"
        echo "Number of cores to use: ${numCores}"

        if bash /work/scripts/hlahd-bam-preprocess-nf.sh ${sampleName} ${inputFiles}; then
            echo "File prepropcessing has finished running on ${sampleName}."

            # assign output fastq to variables
            read1_file=\$(find . -maxdepth 1 -type f -name "*_Bam2Fq_R1.fastq")
            read2_file=\$(find . -maxdepth 1 -type f -name "*_Bam2Fq_R2.fastq")

            echo "Running HLA-HD..."
            # Running HLA-HD on read files
            if bash /work/scripts/hlahd-nf.sh ${sampleName} \${read1_file} \${read2_file} ${numCores}; then
                echo "HLA-HD completed its run. Check outputs for run status."
                # rm -rf *.fastq
            fi
        fi
        """
    }
}
