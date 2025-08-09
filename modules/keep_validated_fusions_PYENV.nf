// 
process KEEP_VALIDATED_FUSIONS_PYENV {
    
    label 'keepValidatedFusions'
    
    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/${sampleName}/POST-IN-SILICO-VALIDATION-TSV-out", mode: 'copy', pattern: "*-collated-FT-normFiltered-FI-validated.tsv",
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(fusInspectorTsv), path(filtered_agfusion_outdir), path(filteredFusions)

    output:
        tuple val(sampleName), path("${sampleName}-collated-FT-normFiltered-FI-validated.tsv"), emit: validatedFusions
        tuple val(sampleName), path("validated-agfusion-outdir/"), emit: validatedAgfusionDir

    script:
    """
    echo "Path to Fusion Inspector validation output: ${fusInspectorTsv}"
    echo "Sample name: ${sampleName}"

    echo "Now retain only validated fusions from Fusion Inspector on the normal-filtered collated FT TSV file..."
    gawk 'NR==FNR {ref[\$1]=1; next} FNR==1 || (\$NF in ref)' ${fusInspectorTsv} ${filteredFusions} > ${sampleName}-collated-FT-normFiltered-FI-validated.tsv
    
    # Create validated agfusion output directory
    mkdir -p validated-agfusion-outdir

    # check if the output file is empty 
    if [ ! -s ${sampleName}-collated-FT-normFiltered-FI-validated.tsv ]; then
        echo "Output file empty. No validated fusions found in Fusion Inspector output. Aborting this sample run." | tee validated-agfusion-outdir/_empty.txt 
        exit 0
    fi
    
    echo "Copying validated AGFusion directories based on Fusion Inspector results..."
    
    # Extract unique gene pairs from Fusion Inspector output (skip header)
    tail -n +2 ${fusInspectorTsv} | cut -f1 | sort -u > validated_gene_pairs.txt
    
    # For each validated gene pair, find and copy matching AGFusion subdirectories
    while IFS= read -r gene_pair; do
        echo "Looking for AGFusion directories containing: \$gene_pair"
        
        # Find subdirectories in AGFusion output that contain the gene pair as substring
        # Use ls instead of find to debug what's actually in the directory
        echo "Contents of AGFusion directory:"
        
        for agf_dir in ${filtered_agfusion_outdir}/*; do
            if [ -d "\$agf_dir" ]; then
                dir_basename=\$(basename "\$agf_dir")
                echo "Checking directory: \$dir_basename"
                
                # Check if the gene pair is contained in the directory name
                if [[ "\$dir_basename" == *"\$gene_pair"* ]]; then
                    echo "Found matching directory: \$dir_basename"
                    echo "Copying \$agf_dir to validated-agfusion-outdir/\$dir_basename"
                    cp -r "\$agf_dir" validated-agfusion-outdir/
                else
                    echo "No match for \$gene_pair in \$dir_basename"
                fi
            fi
        done
        
    done < validated_gene_pairs.txt
    
    # Check if we copied any directories
    copied_dirs=\$(find validated-agfusion-outdir -maxdepth 1 -type d | wc -l)
    if [ \$copied_dirs -le 1 ]; then
        echo "Warning: No matching AGFusion directories found for validated fusions"
    else
        echo "Successfully copied \$((copied_dirs - 1)) AGFusion directories"
    fi
    
    # Ensure directory is never empty for S3 compatibility
    if find validated-agfusion-outdir/ -mindepth 1 -maxdepth 1 -type d -quit | grep -q .; then
        echo "Has subdirectories. Proceeding with outputs..."
    else
        echo "No subdirectories found. Creating a dummy file to ensure directory is not empty."
        echo "No AGFusion directories with protein files found for ${sampleName}" > validated-agfusion-outdir/_placeholder.txt
    fi

    # Clean up temporary file
    rm -f validated_gene_pairs.txt
    """
    
    stub:
    """
    mkdir -p validated-agfusion-outdir
    touch ${sampleName}-collated-FT-normFiltered-FI-validated.tsv
    echo "stub run finished!"
    """
}
