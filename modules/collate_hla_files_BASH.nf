// 
process COLLATE_HLA_FILES_BASH {
    
    label 'collateHLAFiles'
    
    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    
    input:
    val hlaData
    
    output:
    path "Cohort-wide_HLA_types.tsv", emit: cohortWideHLAList
    
    script:

    // Create a string that will be parsed correctly in bash
    def hlaDataStr = hlaData.collect { item -> 
        "${item[0]}:${item[1]}"
    }.join(' ')

    """
    touch "Cohort-wide_HLA_types.tsv"
    echo -e "sampleName\thlaAllotypes" > Cohort-wide_HLA_types.tsv
    
    # Process each sample file pair with a delimiter we can easily parse
    for item in ${hlaDataStr}; do
        # Split the item using the delimiter
        file_path=\$(echo "\$item" | cut -d':' -f2)
        
        # Get the content from the file if it exists
        if [ -f "\$file_path" ]; then
            hla_content=\$(cat "\$file_path")
            
            # Append to the combined file
            echo -e "\$hla_content" >> Cohort-wide_HLA_types.tsv
        else
            echo "Warning: File not found: \$file_path"
        fi
    done
    """

    stub:
    """
    touch "Cohort-wide_HLA_types.tsv"
    echo "Stub run finished!" > test_stub_combine-hla-files.log
    """
}