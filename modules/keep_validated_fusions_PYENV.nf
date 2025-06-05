
// 
process KEEP_VALIDATED_FUSIONS_PYENV {
    
    publishDir "${params.outputDir}/${sampleName}/POST-IN-SILICO-VALIDATION-TSV-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name KEEP-VALIDATED-FUSIONS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(fusInspectorTsv)
        path(filtered_agfusion_outdir)
        tuple val(dummyVar), path(filteredFusions)

    output:
        tuple val(sampleName), path("${sampleName}-collated-FT-normFiltered-FI-validated.tsv"), emit: validatedFusions
        tuple val(sampleName), path("validated-agfusion-outdir"), emit: validatedAgfusionDir

    script:
    """
    echo "Path to Fusion Inspector validation output: ${fusInspectorTsv}"
    echo "Sample name: ${sampleName}"

    echo "Now retain only validated fusions from Fusion Inspector on the normal-filtered collated FT TSV file..."
    gawk 'NR==FNR {ref[\$1]=1; next} FNR==1 || (\$NF in ref)' ${fusInspectorTsv} ${filteredFusions} > ${sampleName}-collated-FT-normFiltered-FI-validated.tsv
    
    # check if the output file is empty 
    if [ ! -s ${sampleName}-collated-FT-normFiltered-FI-validated.tsv ]; then
        echo "Output file empty. No validated fusions found in Fusion Inspector output. Aborting this sample run."
        exit 0
    fi

    # Create validated agfusion output directory
    mkdir -p validated-agfusion-outdir
    
    echo "Copying validated AGFusion directories based on Fusion Inspector results..."
    
    # Extract unique gene pairs from Fusion Inspector output (skip header)
    tail -n +2 ${fusInspectorTsv} | cut -f1 | sort -u > validated_gene_pairs.txt
    
    # For each validated gene pair, find and copy matching AGFusion subdirectories
    while IFS= read -r gene_pair; do
        echo "Looking for AGFusion directories containing: \$gene_pair"
        
        # Find subdirectories in AGFusion output that contain the gene pair as substring
        find ${filtered_agfusion_outdir} -maxdepth 1 -type d -name "*\${gene_pair}*" | while read -r agf_dir; do
            if [ -d "\$agf_dir" ]; then
                dir_basename=\$(basename "\$agf_dir")
                echo "Found matching directory: \$dir_basename"
                echo "Copying \$agf_dir to validated-agfusion-outdir/\$dir_basename"
                cp -r "\$agf_dir" validated-agfusion-outdir/
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
