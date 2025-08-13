// 
process TRANSLATE_IN_SILICO_AGFUSION {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'translateInSilico'
    
    container "${params.container__agfusion}"

    publishDir "${params.outputDir}/${sampleName}/AGFUSION-out", mode: 'copy', overwrite: true
    
    input:
        tuple val(sampleName), path(filteredFusions)

    output:
        tuple val(sampleName), path("filtered-agfusion-dirs/"), emit: filtered_agfusion_outdir

    script:
    """
    echo "Path to filtered fusion transcripts: ${filteredFusions}"
    echo "Running parser to extract fusion transcripts from the filtered lists of fusion transcripts from the fusion calling tools..."
    
    if generate-agfusion-cmd--nf.py -i ${filteredFusions} -c; then
        echo "Parser has finished running the output of the selected FT calling tool of ${sampleName}."
    else
        echo "Failed to run generate-agfusion-cmd--nf.py. Please check the input file and the script."
        exit 1
    fi

    echo "Running AGFusion..."
    mkdir -p filtered-agfusion-dirs

    echo "Running AGFusion..."
    if bash agfusion-cmd.sh; then
        echo "AGFusion has finished."
        echo "Running AGFusion post-processing..."


        # now tranverse the agfusion-dirs directory and do a conditional copy
        #for each directory in agfusion-dirs, if the directory contains a *_protein.fa file, copy the directory to the output directory
        for dir in agfusion-dirs/*; do
            if [ -d "\$dir" ]; then
                if ls "\$dir"/*_protein.fa 1> /dev/null 2>&1; then
                    echo "Copying \$dir to output directory..."
                    cp -r "\$dir" "filtered-agfusion-dirs/"
                else
                    echo "No *_protein.fa file found in \$dir, skipping..."
                fi
            fi
        done
    else
        echo "AGFusion failed to run properly. Please check the logs for errors."
        exit 1
    fi

    # Ensure directory is never empty for S3 compatibility
    if find filtered-agfusion-dirs/ -mindepth 1 -maxdepth 1 -type d -quit | grep -q .; then
        echo "Has subdirectories. Proceeding with outputs..."
    else
        echo "No subdirectories found. Creating a dummy file to ensure directory is not empty."
        echo "No AGFusion directories with protein files found for ${sampleName}" > filtered-agfusion-dirs/_placeholder.txt
    fi
    
    
    echo "AGFusion post-processing completed."

    """
    stub:
    """
    mkdir -p filtered-agfusion-dirs
    echo "stub run finished!" > filtered-agfusion-dirs/stub.out
    """
}
