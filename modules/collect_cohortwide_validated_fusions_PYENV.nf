// Process to collect all TSV files and concatenate them
process COLLECT_COHORTWIDE_VALIDATED_FUSIONS_PYENV {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'collectCohortValidatedFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Fusions-OUT", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
    path validatedFusionsTsvs
    
    output:
    path "Cohortwide_normfiltered_protein-coding_FI-validated_fusions.tsv", emit: cohortwideValidatedFusionsFile
    path "Unique-found-fusions__OUT/Cohortwide_normfiltered_protein-coding_FI-validated_fusions-UNIQUE.manifest.txt"

    script:
    """
    concatenate-cohortwide-fusions--nf.py \\
        --input_files ${validatedFusionsTsvs} \\
        --output Cohortwide_normfiltered_protein-coding_FI-validated_fusions.tsv && \\
        mkdir -p Unique-found-fusions__OUT && \\
        mv Cohortwide_normfiltered_protein-coding_FI-validated_fusions-UNIQUE.manifest.txt Unique-found-fusions__OUT/
    """
    stub:
    """
    touch Cohortwide_normfiltered_protein-coding_FI-validated_fusions.tsv
    mkdir -p Unique-found-fusions__OUT
    touch Unique-found-fusions__OUT/Cohortwide_normfiltered_protein-coding_FI-validated_fusions-UNIQUE.manifest.txt
    echo "[COLLECT_COHORTWIDE_VALIDATED_FUSION_FILES_PYENV]: Stub run finished!"
    """
}
