process TRANSLATE_IN_SILICO_AGFUSION {
    cpus params.numCores
    
    label 'translateInSilico'
    
    container "${params.container__agfusion}"

    publishDir "${params.outputDir}/${sampleName}/AGFUSION-out", mode: 'copy', overwrite: true
    
    input:
        tuple val(sampleName), path(filteredFusions)

    output:
        tuple val(sampleName), path("${sampleName}_filtered-agfusion-dirs.tar.gz"), emit: filtered_agfusion_outdir
        tuple val(sampleName), path("${sampleName}_agfusion-dirs.tar.gz"), emit: raw_agfusion_outdir
        tuple val(sampleName), path("${sampleName}_agfusion_filtered_manifest.txt"), emit: protein_coding_fusions_manifest
        tuple val(sampleName), path("${sampleName}_translate-fusions-with-agfusion-report.log")
        tuple val(sampleName), path("agfusion-logs/")

    script:
    """
    # Set up logging
    LOG_FILE="${sampleName}_translate-fusions-with-agfusion-report.log"
    
    # Function to log messages to both stdout and log file
    log_msg() {
        echo "\$1" | tee -a "\$LOG_FILE"
    }

    log_msg "=== TRANSLATE_IN_SILICO_AGFUSION Process Log ==="
    log_msg "Sample: ${sampleName}"
    log_msg "Start time: \$(date)"
    log_msg "================================================"
    log_msg ""
    
    log_msg "Path to filtered fusion transcripts: ${filteredFusions}"
    log_msg "Running parser to extract fusion transcripts from the filtered lists of fusion transcripts from the fusion calling tools..."

    if generate-agfusion-cmd--nf.py -i ${filteredFusions} -c; then
        log_msg "Parser has finished running the output of the selected FT calling tool of ${sampleName}."
    else
        log_msg "Failed to run generate-agfusion-cmd--nf.py. Please check the input file and the script."
        exit 1
    fi

    log_msg ""
    log_msg "Running AGFusion..."
    mkdir -p filtered-agfusion-dirs

    if bash agfusion-cmd.sh; then
        log_msg "AGFusion has finished."
        log_msg "Running AGFusion post-processing..."

        # Initialize manifest file
        touch ${sampleName}_agfusion_filtered_manifest.txt

        # now traverse the agfusion-dirs directory and do a conditional copy
        # for each directory in agfusion-dirs, if the directory contains a *_protein.fa file, copy the directory to the output directory
        for dir in agfusion-dirs/*; do
            if [ -d "\$dir" ]; then
                if ls "\$dir"/*_protein.fa 1> /dev/null 2>&1; then
                    log_msg "Copying \$dir to output directory..."
                    cp -r "\$dir" "filtered-agfusion-dirs/"
                    FTID=\$(basename "\$dir")
                    # wrangle the FTID to a specific format using awk and gsub
                    # first remove the sample prefix up to the first underscore
                    FTID="\${FTID#*_}"
                    # Use awk to handle each part separately
                    FTID=\$(echo "\$FTID" | awk -F'__' '{
                        # Part 1: Gene pair - replace -- with ::
                        gene_pair = \$1
                        gsub(/--/, "::", gene_pair)
    
                        # Part 2: Coordinates
                        coords = \$2
                        gsub(/-/, ":", coords)      # Replace ALL - with :
                        gsub(/::/, "-", coords)     # Replace :: back to -
                        print gene_pair "__" coords
                    }') && log_msg "Reformatted FT ID: \$FTID"
                    # Record the FTID in manifest
                    echo "\$FTID" >> ${sampleName}_agfusion_filtered_manifest.txt
                else
                    log_msg "No *_protein.fa file found in \$dir, skipping..."
                fi
            fi
        done
    else
        log_msg "AGFusion failed to run properly. Please check the logs for errors."
        exit 1
    fi

    # Ensure directory exists and is never empty for S3 compatibility
    if [ ! -d "agfusion-dirs" ]; then
        log_msg "agfusion-dirs does not exist. Creating it with placeholder."
        mkdir -p agfusion-dirs
        echo "No AGFusion directories (raw outputs) found for ${sampleName}" > agfusion-dirs/_placeholder.txt
    elif [ -z "\$(ls -A agfusion-dirs)" ]; then
        log_msg "agfusion-dirs is empty. Creating placeholder."
        echo "No AGFusion directories (raw outputs) found for ${sampleName}" > agfusion-dirs/_placeholder.txt
    else
        log_msg "agfusion-dirs has content. Proceeding with outputs..."
    fi

    # Ensure directory exists and is never empty for S3 compatibility
    if [ ! -d "filtered-agfusion-dirs" ]; then
        log_msg "filtered-agfusion-dirs does not exist. Creating it with placeholder."
        mkdir -p filtered-agfusion-dirs
        echo "No AGFusion directories with protein files found for ${sampleName}" > filtered-agfusion-dirs/_placeholder.txt
    elif [ -z "\$(ls -A filtered-agfusion-dirs)" ]; then
        log_msg "filtered-agfusion-dirs is empty. Creating placeholder."
        echo "No AGFusion directories with protein files found for ${sampleName}" > filtered-agfusion-dirs/_placeholder.txt
    else
        log_msg "filtered-agfusion-dirs has content. Proceeding with outputs..."
    fi

    log_msg ""
    log_msg "AGFusion post-processing completed."
    filtered_count=\$(cat ${sampleName}_agfusion_filtered_manifest.txt | wc -l)
    log_msg "Total filtered directories: \$filtered_count"
    
    # Package directories as tarballs for reliable S3 staging
    log_msg ""
    log_msg "Packaging directories as tarballs for S3 compatibility..."
    tar -czf ${sampleName}_filtered-agfusion-dirs.tar.gz filtered-agfusion-dirs/
    log_msg "Created ${sampleName}_filtered-agfusion-dirs.tar.gz"
    
    tar -czf ${sampleName}_agfusion-dirs.tar.gz agfusion-dirs/
    log_msg "Created ${sampleName}_agfusion-dirs.tar.gz"
    
    log_msg ""
    log_msg "================================================"
    log_msg "End time: \$(date)"
    log_msg "Process completed successfully"
    """

    stub:
    """
    LOG_FILE="${sampleName}_translate-fusions-with-agfusion-report.log"
    mkdir -p filtered-agfusion-dirs
    mkdir -p agfusion-dirs
    mkdir -p agfusion-logs
    echo "stub run finished!" > filtered-agfusion-dirs/stub.out
    echo "stub run finished!" > agfusion-dirs/stub.out
    echo "stub run finished!" > agfusion-logs/stub.out
    echo "# Stub run - no actual directories" > ${sampleName}_agfusion_filtered_manifest.txt
    echo "=== STUB RUN ===" > "\$LOG_FILE"
    echo "Sample: ${sampleName}" >> "\$LOG_FILE"
    echo "Stub execution completed at: \$(date)" >> "\$LOG_FILE"
    
    # Create stub tarballs
    tar -czf ${sampleName}_filtered-agfusion-dirs.tar.gz filtered-agfusion-dirs/
    tar -czf ${sampleName}_agfusion-dirs.tar.gz agfusion-dirs/
    """
}
