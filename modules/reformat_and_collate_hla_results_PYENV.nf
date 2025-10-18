// Streamlined reformat and collate process
process REFORMAT_AND_COLLATE_HLA_RESULTS_PYENV {
    cpus 1
    
    label 'reformatCollateHLAs'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-HLA-typing-OUT", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    input:
    path jsonFiles  // All JSON files collected
    
    output:
    path "Cohortwide_HLA_types.tsv", emit: cohortWideHLAList
    path "Sample-specific_HLA_results/*", emit: individualResults
    
    script:
    """
    mkdir -p Sample-specific_HLA_results
    
    # Initialize output files
    echo -e "sampleName\thlaAllotypes" > Cohortwide_HLA_types.tsv
    
    # Process each JSON file from the collected input
    for json_file in ${jsonFiles}; do
        if [[ -f "\$json_file" ]]; then
            sample_name=\$(basename "\$json_file" .genotype.json)
            
            # Reformat individual file using your existing script
            reformat-hlas--nf.py "\$sample_name" "\$json_file" > "Sample-specific_HLA_results/\${sample_name}-HLA-types-reformatted.tsv"
            
            # Extract HLA content for cohort file
            hla_content=\$(tail -n1 "Sample-specific_HLA_results/\${sample_name}-HLA-types-reformatted.tsv" | cut -f2-)
            
            # Add to cohort file
            echo -e "\$sample_name\t\$hla_content" >> Cohortwide_HLA_types.tsv
        fi
    done
    
    # Sort cohort file by sample name
    { head -n1 Cohortwide_HLA_types.tsv; tail -n+2 Cohortwide_HLA_types.tsv | sort -V; } > temp_sorted.tsv
    mv temp_sorted.tsv Cohortwide_HLA_types.tsv
    
    # Summary statistics
    total_samples=\$((\$(wc -l < Cohortwide_HLA_types.tsv) - 1))
    echo "Total samples processed: \$total_samples"
    echo "Cohortwide HLA types file created successfully"
    """
    
    stub:
    """
    mkdir -p Sample-specific_HLA_results
    touch Cohortwide_HLA_types.tsv
    echo "Stub run finished!" > test_stub_reformat_and_collate.log
    """
}
