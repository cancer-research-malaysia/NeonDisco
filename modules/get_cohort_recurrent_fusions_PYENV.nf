// Process to filter for recurrent fusions
process GET_COHORT_RECURRENT_FUSIONS_PYENV {
    
    label 'getCohortRecurrentFusions'
    
    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Fusions-OUT/Recurrent-NormFiltered-Fusions", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    
    input:
    path cohortwideFusionsFile
    
    output:
    tuple val("Cohortwide-Recurrent-Fusions"), path("Cohortwide_normfiltered_fusions_RECURRENT.tsv"), emit: cohortRecurrentFusionTsv
    tuple val("Cohort-Recurrent-Fusions"), path("Recurrent_normfiltered_fusion_frequency_report.txt"), emit: fusionFrequencyReport
    
    script:
    """
    get-cohortwide-recurrent-fusions--nf.py \\
        --input ${cohortwideFusionsFile} \\
        --threshold ${params.recurrenceThreshold} \\
        --output Cohortwide_normfiltered_fusions_RECURRENT.tsv \\
        --report Recurrent_normfiltered_fusion_frequency_report.txt
    """
    stub:
    """
    touch Cohortwide_normfiltered_fusions_RECURRENT.tsv
    touch Recurrent_normfiltered_fusion_frequency_report.txt
    echo "[GET_COHORT_RECURRENT_FUSIONS_PYENV]: Stub run finished!"
    """

}
