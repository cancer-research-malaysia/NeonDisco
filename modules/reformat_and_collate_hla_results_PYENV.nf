// Streamlined reformat and collate process
process REFORMAT_AND_COLLATE_HLA_RESULTS_PYENV {
    
    label 'reformatCollateHLAs'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    input:
    path jsonFiles  // All JSON files collected
    
    output:
    path "Cohort-wide_HLA_types.tsv", emit: cohortWideHLAList
    path "individual_hla_results/*", emit: individualResults
    
    script:
    """
    mkdir -p individual_hla_results
    
    # Initialize output files
    echo -e "sampleName\thlaAllotypes" > Cohort-wide_HLA_types.tsv
    
    # Process each JSON file from the collected input
    for json_file in ${jsonFiles}; do
        if [[ -f "\$json_file" ]]; then
            sample_name=\$(basename "\$json_file" .genotype.json)
            
            # Reformat individual file using your existing script
            reformat-hlas--nf.py "\$sample_name" "\$json_file" > "individual_hla_results/\${sample_name}-HLA-types-reformatted.tsv"
            
            # Extract HLA content for cohort file
            hla_content=\$(tail -n1 "individual_hla_results/\${sample_name}-HLA-types-reformatted.tsv" | cut -f2-)
            
            # Add to cohort file
            echo -e "\$sample_name\t\$hla_content" >> Cohort-wide_HLA_types.tsv
        fi
    done
    
    # Sort cohort file by sample name
    { head -n1 Cohort-wide_HLA_types.tsv; tail -n+2 Cohort-wide_HLA_types.tsv | sort -V; } > temp_sorted.tsv
    mv temp_sorted.tsv Cohort-wide_HLA_types.tsv
    
    # Summary statistics
    total_samples=\$((\$(wc -l < Cohort-wide_HLA_types.tsv) - 1))
    echo "Total samples processed: \$total_samples"
    echo "Cohort-wide HLA types file created successfully"
    """
    
    stub:
    """
    mkdir -p individual_hla_results
    touch Cohort-wide_HLA_types.tsv
    echo "Stub run finished!" > test_stub_reformat_and_collate.log
    """
}
