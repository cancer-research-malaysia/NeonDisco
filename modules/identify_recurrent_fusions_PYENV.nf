// Process to collect all TSV files and concatenate them
process CONCAT_NORMFILTERED_FUSION_FILES_PYENV {
    
    label 'concatNormFilteredFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-NormFiltered-Fusions", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
    path normFilteredFusionsTsvs
    
    output:
    path "Cohortwide_norm-filtered_fusions.tsv", emit: cohortwideFusionsFile
    
    script:
    """
    concatenate-cohortwide-fusions--nf.py \\
        --input_files ${normFilteredFusionsTsvs} \\
        --output Cohortwide_norm-filtered_fusions.tsv
    """
    stub:
    """
    touch Cohortwide_norm-filtered_fusions.tsv
    echo "[CONCAT_NORMFILTERED_FUSION_FILES_PYENV]: Stub run finished!"
    """
}

// Process to filter for recurrent fusions
process GET_COHORT_RECURRENT_FUSIONS_PYENV {
    
    label 'getCohortRecurrentFusions'
    
    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Recurrent-NormFiltered-Fusions", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    
    input:
    path cohortwideFusionsFile
    
    output:
    tuple val("Cohortwide-Recurrent-Fusions"), path("Cohortwide_recurrent_normfiltered_fusions.tsv"), emit: cohortRecurrentFusionTsv
    tuple val("Cohort-Recurrent-Fusions"), path("Recurrent_normfiltered_fusion_frequency_report.txt"), emit: fusionFrequencyReport
    
    script:
    """
    get-cohortwide-recurrent-fusions--nf.py \\
        --input ${cohortwideFusionsFile} \\
        --threshold ${params.recurrenceThreshold} \\
        --output Cohortwide_recurrent_normfiltered_fusions.tsv \\
        --report Recurrent_normfiltered_fusion_frequency_report.txt
    """
    stub:
    """
    touch Cohortwide_recurrent_normfiltered_fusions.tsv
    touch Recurrent_normfiltered_fusion_frequency_report.txt
    echo "[GET_COHORT_RECURRENT_FUSIONS_PYENV]: Stub run finished!"
    """

}
