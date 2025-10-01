// Process to filter for recurrent fusions
process GET_COHORTWIDE_RECURRENT_VALIDATED_FUSIONS_PYENV {
    
    label 'getCohortRecurrentFusionsFiValidated'
    
    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Fusions-OUT/Recurrent-Normfiltered-FI-Validated-Fusions", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    
    input:
    path cohortwideFusionsFile
    
    output:
    tuple val("Cohortwide-Recurrent-Validated-Fusions"), path("Cohortwide_normfiltered_FI_validated_fusions_recurrent.tsv"), emit: cohortRecurrentFusionTsv
    tuple val("Cohort-Recurrent-Validated-Fusions"), path("Recurrent_normfiltered_FI_validated_fusion_frequency_report.*"), emit: fusionFrequencyReport
    
    script:
    """
    get-cohortwide-recurrent-fusions--nf.py \\
        --input ${cohortwideFusionsFile} \\
        --threshold ${params.recurrenceThreshold} \\
        --output Cohortwide_normfiltered_FI_validated_fusions_recurrent.tsv \\
        --report Recurrent_normfiltered_FI_validated_fusion_frequency_report.txt
    """
    stub:
    """
    touch Cohortwide_normfiltered_FI_validated_fusions_recurrent.tsv
    touch Recurrent_normfiltered_FI_validated_fusion_frequency_report.txt
    touch Recurrent_normfiltered_FI_validated_fusion_frequency_report.tsv

    echo "[GET_COHORT_RECURRENT_FUSIONS_FI_VALIDATED_PYENV]: Stub run finished!"
    """

}
