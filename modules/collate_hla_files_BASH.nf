// 
process COLLATE_HLA_FILES_BASH {
    
    label 'collateHLAFiles'
    
    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }


    input:
    path hlaFiles
    path mappingFile
    
    output:
    path "Cohort-wide_HLA_types.tsv", emit: cohortWideHLAList
    
    script:
    """
    touch "Cohort-wide_HLA_types.tsv"
    echo -e "sampleName\thlaAllotypes" > Cohort-wide_HLA_types.tsv

    # Process each file using the explicit mapping
    while IFS=\$'\t' read -r sample_name filename; do
        if [ -f "\$filename" ]; then
            # Extract only the HLA allotypes (everything after the first column)
            hla_content=\$(cut -f2- "\$filename")
            
            # Check if content is empty
            if [ -z "\$hla_content" ]; then
                echo "Warning: No HLA data found in \$filename for sample \$sample_name" >&2
            fi
            # Append to the combined file
            echo -e "\$sample_name\t\$hla_content" >> Cohort-wide_HLA_types.tsv
        else
            echo "Warning: File not found: \$filename for sample \$sample_name" >&2
            exit 1
        fi
    done < ${mappingFile}
    echo "Cohort-wide HLA types file created successfully."
    """

    stub:
    """
    touch "Cohort-wide_HLA_types.tsv"
    echo "Stub run finished!" > test_stub_combine-hla-files.log
    """
}