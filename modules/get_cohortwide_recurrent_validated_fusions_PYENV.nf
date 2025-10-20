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
    tuple val("Cohortwide-Recurrent-Validated-Fusions"), path("Cohortwide_normfiltered_FI-validated_recurrent_fusions.tsv"), emit: cohortRecurrentFusionTsv
    tuple val("Cohortwide-Recurrent-Validated-Fusions"), path("Cohortwide_normfiltered_FI-validated_recurrent_fusion_freq-report.*"), emit: fusionFrequencyReport
    
    script:
    """
    get-cohortwide-recurrent-fusions--nf.py \\
        --input ${cohortwideFusionsFile} \\
        --threshold ${params.recurrenceThreshold} \\
        --output Cohortwide_normfiltered_FI-validated_recurrent_fusions.tsv \\
        --report Cohortwide_normfiltered_FI-validated_recurrent_fusion_freq-report.txt \\
        --total-samples ${totalSampleCount}
    """
    stub:
    """
    touch Cohortwide_normfiltered_FI-validated_recurrent_fusions.tsv
    touch Cohortwide_normfiltered_FI-validated_recurrent_fusion_freq-report.txt
    touch Cohortwide_normfiltered_FI-validated_recurrent_fusion_freq-report.tsv

    echo "[GET_COHORTWIDE_RECURRENT_VALIDATED_FUSIONS_PYENV]: Stub run finished!"
    """

}
