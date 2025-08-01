// Process to collect all TSV files and concatenate them
process COLLECT_COHORTWIDE_VALIDATED_FUSIONS_PYENV {
    
    label 'collectCohortValidatedFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-FI-Validated-Fusions", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
    path validatedFusionsTsvs
    
    output:
    path "Cohortwide_FI-validated_fusions.tsv", emit: cohortwideValidatedFusionsFile
    
    script:
    """
    concatenate-cohortwide-fusions--nf.py \\
        --input_files ${validatedFusionsTsvs} \\
        --output Cohortwide_FI-validated_fusions.tsv
    """
    stub:
    """
    touch Cohortwide_FI-validated_fusions.tsv
    echo "[COLLECT_COHORTWIDE_VALIDATED_FUSION_FILES_PYENV]: Stub run finished!"
    """
}
