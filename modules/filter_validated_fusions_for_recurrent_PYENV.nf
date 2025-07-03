// Alternative approach: Create a separate recurrent-aware filtering step
process FILTER_VALIDATED_FUSIONS_FOR_RECURRENT_PYENV {
    
    publishDir "${params.outputDir}/${sampleName}/POST-IN-SILICO-VALIDATION-TSV-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FILTER-FOR-RECURRENT -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
    tuple val(sampleName), path(validatedFusions)
    tuple val(dummyVar), path(validatedAgfusionDir)
    tuple val(cohortLabel), path(recurrentFusionsTsv)

    output:
    tuple val(sampleName), path("${sampleName}-validated-recurrent-only.tsv"), emit: validatedRecurrentFusions
    tuple val(sampleName), path("validated-recurrent-agfusion-outdir"), emit: validatedRecurrentAgfusionDir
    path("${sampleName}_recurrent_filter_report.txt"), emit: recurrentFilterReport

    script:
    """
    echo "Filtering validated fusions for recurrent ones only..."
    echo "Sample: ${sampleName}"
    echo "Validated fusions: ${validatedFusions}"
    echo "Recurrent fusions TSV: ${recurrentFusionsTsv}"

    # Initialize report
    echo "Recurrent Fusion Filtering Report for ${sampleName}" > ${sampleName}_recurrent_filter_report.txt
    echo "===================================================" >> ${sampleName}_recurrent_filter_report.txt
    echo "Date: \$(date)" >> ${sampleName}_recurrent_filter_report.txt
    echo "" >> ${sampleName}_recurrent_filter_report.txt

    # Check input files
    if [ $(tail -n +2 ${validatedFusions} | wc -l) -eq 0 ]; then
        echo "No validated fusions for this sample."
        echo "STATUS: No validated fusions" >> ${sampleName}_recurrent_filter_report.txt
        head -1 ${validatedFusions} > ${sampleName}-validated-recurrent-only.tsv
        mkdir -p validated-recurrent-agfusion-outdir
        exit 0
    fi

    if [ $(tail -n +2 ${recurrentFusionsTsv} | wc -l) -eq 0 ]; then
        echo "No recurrent fusions in cohort."
        echo "STATUS: No recurrent fusions in cohort" >> ${sampleName}_recurrent_filter_report.txt
        head -1 ${validatedFusions} > ${sampleName}-validated-recurrent-only.tsv
        mkdir -p validated-recurrent-agfusion-outdir
        exit 0
    fi

    # Get validated fusion count
    validated_count=\$(tail -n +2 ${validatedFusions} | wc -l)
    echo "Input validated fusions: \$validated_count" >> ${sampleName}_recurrent_filter_report.txt

    # Extract recurrent gene pairs for this sample
    tail -n +2 ${recurrentFusionsTsv} | awk -F'\\t' -v sample="${sampleName}" '\$8 == sample {print \$2}' | sort -u > sample_recurrent_gene_pairs.txt
    
    sample_recurrent_count=\$(wc -l < sample_recurrent_gene_pairs.txt)
    echo "Recurrent gene pairs for this sample: \$sample_recurrent_count" >> ${sampleName}_recurrent_filter_report.txt

    # Filter validated fusions for recurrent ones
    head -1 ${validatedFusions} > ${sampleName}-validated-recurrent-only.tsv
    
    if [ \$sample_recurrent_count -gt 0 ]; then
        while IFS= read -r recurrent_gene_pair; do
            tail -n +2 ${validatedFusions} | awk -F'\\t' -v gene_pair="\$recurrent_gene_pair" '\$2 == gene_pair' >> ${sampleName}-validated-recurrent-only.tsv
        done < sample_recurrent_gene_pairs.txt
    fi

    final_count=\$(tail -n +2 ${sampleName}-validated-recurrent-only.tsv | wc -l)
    echo "Final recurrent validated fusions: \$final_count" >> ${sampleName}_recurrent_filter_report.txt

    # Filter AGFusion directories
    mkdir -p validated-recurrent-agfusion-outdir
    
    if [ \$final_count -gt 0 ]; then
        while IFS= read -r recurrent_gene_pair; do
            for agf_dir in ${validatedAgfusionDir}/*; do
                if [ -d "\$agf_dir" ]; then
                    dir_basename=\$(basename "\$agf_dir")
                    if [[ "\$dir_basename" == *"\$recurrent_gene_pair"* ]]; then
                        cp -r "\$agf_dir" validated-recurrent-agfusion-outdir/
                        echo "COPIED: \$dir_basename" >> ${sampleName}_recurrent_filter_report.txt
                    fi
                fi
            done
        done < sample_recurrent_gene_pairs.txt
    fi

    # Summary
    echo "" >> ${sampleName}_recurrent_filter_report.txt
    echo "SUMMARY:" >> ${sampleName}_recurrent_filter_report.txt
    echo "Input validated fusions: \$validated_count" >> ${sampleName}_recurrent_filter_report.txt
    echo "Final recurrent fusions: \$final_count" >> ${sampleName}_recurrent_filter_report.txt
    echo "Reduction: \$(( (validated_count - final_count) * 100 / validated_count ))% filtered out" >> ${sampleName}_recurrent_filter_report.txt

    # Cleanup
    rm -f sample_recurrent_gene_pairs.txt
    """
    
    stub:
    """
    mkdir -p validated-recurrent-agfusion-outdir
    touch ${sampleName}-validated-recurrent-only.tsv
    touch ${sampleName}_recurrent_filter_report.txt
    echo "FILTER_VALIDATED_FOR_RECURRENT_PYENV: Stub run finished!" > ${sampleName}_recurrent_filter_report.txt
    """
}
