
// Process to filter for recurrent fusions
process GET_COHORTWIDE_RECURRENT_VALIDATED_FUSIONS_PYENV {
    cpus 2

    label 'getCohortRecurrentFusionsFiValidated'
    
    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Fusions-OUT/Cohortwide-Recurrent-Fusions-OUT", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    
    input:
    path cohortwideFusionsFile
    val totalSampleCount
    
    output:
    tuple val("Cohortwide-Recurrent-Validated-Fusions"), path("Cohortwide_normfiltered_FI-validated_recurrent_fusions.thresholded.tsv"), emit: cohortRecurrentFusionTsv
    tuple val("Cohortwide-Recurrent-Validated-Fusions"), path("Cohortwide_normfiltered_FI-validated_recurrent_fusions.thresholded.freq-report.*"), emit: fusionFrequencyReport
    tuple val("Cohortwide-Unthresholded-Recurrent-Validated-Fusions"), path("Cohortwide_normfiltered_FI-validated_recurrent_fusions.unthresh.tsv"), emit: cohortUnthresholdedRecurrentFusionTsv
    tuple val("Cohortwide-Unthresholded-Recurrent-Validated-Fusions"), path("Cohortwide_normfiltered_FI-validated_recurrent_fusions.unthresh.freq-report.*"), emit: unthresholdedFusionFrequencyReport
    tuple val("Cohortwide-Recurrent-Fusion-Visualization"), path("Cohortwide_recurrent_fusion_stats-vis_unthresh.html"), emit: fusionVisualizationUnthresholded, optional: true
    tuple val("Cohortwide-Recurrent-Fusion-Visualization"), path("Cohortwide_recurrent_fusion_stats-vis_thresholded.html"), emit: fusionVisualizationThresholded, optional: true
    tuple val("Cohortwide-Recurrent-Fusion-Heatmap"), path("Cohortwide_recurrent_fusion_cluster-heatmap-vis_unthresh.html"), emit: fusionHeatmapUnthresholded, optional: true
    tuple val("Cohortwide-Recurrent-Fusion-Heatmap"), path("Cohortwide_recurrent_fusion_cluster-heatmap-vis_thresholded.html"), emit: fusionHeatmapThresholded, optional: true
    
    script:
    """
    # Generate recurrent fusion files
    get-cohortwide-recurrent-fusions--nf.py \\
        --input ${cohortwideFusionsFile} \\
        --threshold ${params.recurrenceThreshold} \\
        --output Cohortwide_normfiltered_FI-validated_recurrent_fusions.thresholded.tsv \\
        --report Cohortwide_normfiltered_FI-validated_recurrent_fusions.thresholded.freq-report.txt \\
        --total-samples ${totalSampleCount} \\
        --unthresholded-output Cohortwide_normfiltered_FI-validated_recurrent_fusions.unthresh.tsv || echo "Recurrent fusion extraction failed"
    
    # Generate visualizations if cohort is large enough (>= 10 samples)
    if [ ${totalSampleCount} -ge 10 ]; then
        echo "Generating distribution visualizations..."
        
        # Unthresholded distribution (pie + bar charts)
        visualize-recurrent-fusion-stats--nf.py \\
            --input Cohortwide_normfiltered_FI-validated_recurrent_fusions.unthresh.freq-report.tsv \\
            --output Cohortwide_recurrent_fusion_stats-vis_unthresh.html \\
            --min-cohort-size 10 \\
            --report-type Unthresholded || echo "Unthresholded visualization generation failed or skipped"
        
        # Thresholded distribution (pie + bar charts)
        visualize-recurrent-fusion-stats--nf.py \\
            --input Cohortwide_normfiltered_FI-validated_recurrent_fusions.thresholded.freq-report.tsv \\
            --output Cohortwide_recurrent_fusion_stats-vis_thresholded.html \\
            --min-cohort-size 10 \\
            --report-type Thresholded || echo "Thresholded visualization generation failed or skipped"
        
        echo "Generating patient-fusion coverage heatmaps..."
        
        # Unthresholded heatmap (patient × fusion matrix)
        visualize-recurrent-fusion-cohort-coverage--nf.py \\
            --input Cohortwide_normfiltered_FI-validated_recurrent_fusions.unthresh.freq-report.tsv \\
            --output Cohortwide_recurrent_fusion_cluster-heatmap-vis_unthresh.html \\
            --report-type Unthresholded || echo "Unthresholded heatmap generation failed or skipped"
        
        # Thresholded heatmap (patient × fusion matrix)
        visualize-recurrent-fusion-cohort-coverage--nf.py \\
            --input Cohortwide_normfiltered_FI-validated_recurrent_fusions.thresholded.freq-report.tsv \\
            --output Cohortwide_recurrent_fusion_cluster-heatmap-vis_thresholded.html \\
            --report-type Thresholded || echo "Thresholded heatmap generation failed or skipped"
    else
        echo "Cohort size (${totalSampleCount}) below minimum (10) for visualization. Skipping."
    fi
    """
    stub:
    """
    touch Cohortwide_normfiltered_FI-validated_recurrent_fusions.thresholded.tsv
    touch Cohortwide_normfiltered_FI-validated_recurrent_fusions.thresholded.freq-report.txt
    touch Cohortwide_normfiltered_FI-validated_recurrent_fusions.thresholded.freq-report.tsv
    touch Cohortwide_normfiltered_FI-validated_recurrent_fusions.unthresh.tsv
    touch Cohortwide_normfiltered_FI-validated_recurrent_fusions.unthresh.freq-report.txt
    touch Cohortwide_normfiltered_FI-validated_recurrent_fusions.unthresh.freq-report.tsv

    echo "[GET_COHORTWIDE_RECURRENT_VALIDATED_FUSIONS_PYENV]: Stub run finished!"
    """
}
