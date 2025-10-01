// Process to collect all TSV files and concatenate them
process COLLECT_COHORTWIDE_NORMFILTERED_FUSIONS_PYENV {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'collectNormFilteredFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Fusions-OUT", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
    path normFilteredFusionsTsvs
    
    output:
    path "Cohortwide_normfiltered_fusions.tsv", emit: cohortwideFusionsFile
    path "Cohortwide_normfiltered_fusions-UNIQUE.manifest.txt"
    
    script:
    """
    concatenate-cohortwide-fusions--nf.py \\
        --input_files ${normFilteredFusionsTsvs} \\
        --output Cohortwide_normfiltered_fusions.tsv
    """
    stub:
    """
    touch Cohortwide_normfiltered_fusions.tsv
    touch Cohortwide_normfiltered_fusions-UNIQUE.manifest.txt
    echo "[CONCAT_NORMFILTERED_FUSION_FILES_PYENV]: Stub run finished!"
    """
}

