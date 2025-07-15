// Process to collect all TSV files and concatenate them
process CONCAT_NORMFILTERED_FUSION_FILES_PYENV {
    
    label 'concatNormFilteredFusions'

    publishDir "${params.outputDir}/Cohortwide-Total-Fusions", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    
    input:
    path validatedFusionsTsvs
    
    output:
    path "Cohortwide_total_fusions.tsv", emit: cohortwideFusionsFile
    
    script:
    """
    concatenate-cohortwide-fusions--nf.py \\
        --input_files ${validatedFusionsTsvs} \\
        --output Cohortwide_total_fusions.tsv
    """
    stub:
    """
    touch Cohortwide_total_fusions.tsv
    echo "[CONCAT_NORMFILTERED_FUSION_FILES_PYENV]: Stub run finished!"
    """
}

// Process to filter for recurrent fusions
process GET_COHORT_RECURRENT_FUSIONS_PYENV {
    
    label 'getCohortRecurrentFusions'
    
    publishDir "${params.outputDir}/Cohortwide-Recurrent-Fusions", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    container "${params.container__pyenv}"
    
    input:
    path cohortwideFusionsFile
    
    output:
    tuple val("Cohortwide-Recurrent-Fusions"), path("Cohortwide_recurrent_fusions.tsv"), emit: cohortRecurrentFusionTsv
    tuple val("Cohort-Recurrent-Fusions"), path("Recurrent_fusion_frequency_report.txt"), emit: fusionFrequencyReport
    
    script:
    """
    get-cohortwide-recurrent-fusions--nf.py \\
        --input ${cohortwideFusionsFile} \\
        --threshold ${params.recurrenceThreshold} \\
        --output Cohortwide_recurrent_fusions.tsv \\
        --report Recurrent_fusion_frequency_report.txt
    """
    stub:
    """
    touch Cohortwide_recurrent_fusions.tsv
    touch Recurrent_fusion_frequency_report.txt
    echo "[GET_COHORT_RECURRENT_FUSIONS_PYENV]: Stub run finished!"
    """

}