// 
process TRANSLATE_IN_SILICO_AGFUSION {
    
    publishDir "${params.outputDir}/${sampleName}/AGFUSION-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__agfusion}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name TRANSLATE-IN-SILICO -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(filteredFusions)

    output:
        tuple val(sampleName), path("filtered-agfusion-dirs/"), emit: agfusion_outdir

    script:
    """
    echo "Path to filtered fusion transcripts: ${filteredFusions}"
    echo "Running parser to extract fusion transcripts from the filtered lists of fusion transcripts from the fusion calling tools..."
    if python /home/app/scripts/generate-agfusion-cmd--nf.py -i ${filteredFusions} -c; then
        echo "Parser has finished running the output of the selected FT calling tool of ${sampleName}."
    fi

    echo "Running AGFusion..."
    if bash agfusion-cmd.sh; then
        echo "AGFusion has finished."
        echo "Running AGFusion post-processing..."
        mkdir -p filtered-agfusion-dirs
        # now tranverse the agfusion-dirs directory and do a conditional copy, for each directory in agfusion-dirs, if the directory contains a *_protein.fa file, copy the directory to the output directory
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
    fi

    
    """
    stub:
    """
    mkdir -p filtered-agfusion-dirs
    echo "stub run finished!" > filtered-agfusion-dirs/stub.out
    """
}
