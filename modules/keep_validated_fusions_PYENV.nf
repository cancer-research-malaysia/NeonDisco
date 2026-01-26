// 
process KEEP_VALIDATED_FUSIONS_PYENV {
    cpus 1
    
    label 'keepValidatedFusions'
    
    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/${sampleName}/POST-IN-SILICO-VALIDATION-TSV-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(fusInspectorTsv), path(filteredAgfusionTarball), path(filteredFusions), path(proteinCodingFusManifest)

    output:
        tuple val(sampleName), path("${sampleName}-collated-FT-normFiltered-protein-coding-only-FI-validated.tsv"), optional: true, emit: validatedFusions
        tuple val(sampleName), path("validated-agfusion-dirs/"), optional: true, emit: validatedAgfusionDir
        tuple val(sampleName), path("${sampleName}-validated-fusion-filtering.log")

    script:
    """
    echo "Path to Fusion Inspector validation output: ${fusInspectorTsv}" | tee "${sampleName}-validated-fusion-filtering.log"
    echo "Sample name: ${sampleName}" | tee -a "${sampleName}-validated-fusion-filtering.log"
    echo "Path to filtered AGFusion output tarball: ${filteredAgfusionTarball}" | tee -a "${sampleName}-validated-fusion-filtering.log"
    echo "Extracting filtered AGFusion output directory..." | tee -a "${sampleName}-validated-fusion-filtering.log"
    tar -xzf ${filteredAgfusionTarball}

    echo "First filter the collated FT TSV file to only include NORMAL-FILTERED, PROTEIN-CODING fusions, then pipe this to retain only those fusions validated by Fusion Inspector..."
    gawk 'NR==FNR {ref[\$1]=1; next} FNR==1 || (\$1 in ref)' ${proteinCodingFusManifest} ${filteredFusions} | gawk 'NR==FNR {ref[\$1]=1; next} FNR==1 || (\$NF in ref)' ${fusInspectorTsv} - > ${sampleName}-collated-FT-normFiltered-protein-coding-only-FI-validated.tsv

    # Check if the output file has data beyond the header
    data_rows=\$(tail -n +2 ${sampleName}-collated-FT-normFiltered-protein-coding-only-FI-validated.tsv | wc -l)
    if [ "\$data_rows" -eq 0 ]; then
        echo "No validated fusions found in Fusion Inspector output. Aborting this sample run." | tee -a "${sampleName}-validated-fusion-filtering.log"
        # delete the empty output tsv as well
        rm -f ${sampleName}-collated-FT-normFiltered-protein-coding-only-FI-validated.tsv
        exit 0
    fi
    
    # Create validated agfusion output directory
    mkdir -p validated-agfusion-dirs

    echo "Copying validated AGFusion directories based on Fusion Inspector results..." | tee -a "${sampleName}-validated-fusion-filtering.log"
    
    # Extract unique gene pairs from Fusion Inspector output (skip header)
    tail -n +2 ${fusInspectorTsv} | cut -f1 | sort -u > validated_gene_pairs.txt
    
    # For each validated gene pair, find and copy matching AGFusion subdirectories
    while IFS= read -r gene_pair; do
        echo "Looking for AGFusion directories containing: \$gene_pair" | tee -a "${sampleName}-validated-fusion-filtering.log"
        
        # Find subdirectories in AGFusion output that contain the gene pair as substring
        # Use ls instead of find to debug what's actually in the directory
        echo "Contents of AGFusion directory:" | tee -a "${sampleName}-validated-fusion-filtering.log"
        
        for agf_dir in filtered-agfusion-dirs/*; do
            if [ -d "\$agf_dir" ]; then
                dir_basename=\$(basename "\$agf_dir")
                echo "Checking directory: \$dir_basename" | tee -a "${sampleName}-validated-fusion-filtering.log"
                
                # Check if the gene pair is contained in the directory name
                if [[ "\$dir_basename" == *"\$gene_pair"* ]]; then
                    echo "Found matching directory: \$dir_basename" | tee -a "${sampleName}-validated-fusion-filtering.log"
                    echo "Copying \$agf_dir to validated-agfusion-dirs/\$dir_basename" | tee -a "${sampleName}-validated-fusion-filtering.log"
                    cp -r "\$agf_dir" validated-agfusion-dirs/
                else
                    echo "No match for \$gene_pair in \$dir_basename" | tee -a "${sampleName}-validated-fusion-filtering.log"
                fi
            fi
        done
        
    done < validated_gene_pairs.txt
    
    # Check if we copied any directories
    copied_dirs=\$(find validated-agfusion-dirs -mindepth 1 -maxdepth 1 -type d | wc -l)
    if [ \$copied_dirs -eq 0 ]; then
        echo "Warning: No matching AGFusion directories found for validated fusions (data mismatch)" | tee -a "${sampleName}-validated-fusion-filtering.log"
    else
        echo "Successfully copied \$copied_dirs AGFusion directories" | tee -a "${sampleName}-validated-fusion-filtering.log"
    fi

    # Clean up temporary file
    rm -f validated_gene_pairs.txt
    """
    
    stub:
    """
    mkdir -p validated-agfusion-dirs
    touch ${sampleName}-collated-FT-normFiltered-protein-coding-only-FI-validated.tsv
    touch ${sampleName}-validated-fusion-filtering.log
    echo "stub run finished!"
    """
}
