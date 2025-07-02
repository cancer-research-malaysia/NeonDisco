// Process to collect all TSV files and concatenate them
process CONCAT_NORMFILTERED_FUSION_FILES_PYENV {
    
    publishDir "${params.outputDir}", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name CONCAT_FUSION_FILES -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
    path cohortwideFusionsTsvs
    
    output:
    path "cohortwide_concat_fusions.tsv", emit: cohortwideFusionsFile
    
    script:
    """
    python3 /home/app/scripts/concatenate-cohortwide-fusions--nf.py \\
        --input_files ${cohortwideFusionsTsvs} \\
        --output cohortwide_concat_fusions.tsv
    """
    stub:
    """
    touch cohortwide_concat_fusions.tsv
    echo "[CONCAT_NORMFILTERED_FUSION_FILES_PYENV]: Stub run finished!"
    """
}

// Process to filter for recurrent fusions
process GET_COHORT_RECURRENT_FUSIONS_PYENV {
    
    publishDir "${params.outputDir}", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name IDENTIFY_RECURRENT_FUSIONS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
    path cohortwideFusionsFile
    
    output:
    path "cohortwide_recurrent_fusions.tsv", emit: cohortRecurrentFusionTsv
    path "fusion_frequency_report.txt", emit: fusionFrequencyReport
    
    script:
    """
    python3 /home/app/scripts/get-cohortwide-recurrent-fusions--nf.py \\
        --input ${cohortwideFusionsFile} \\
        --threshold ${params.recurrenceThreshold} \\
        --output cohortwide_recurrent_fusions.tsv \\
        --report fusion_frequency_report.txt
    """
    stub:
    """
    touch cohortwide_recurrent_fusions.tsv
    touch fusion_frequency_report.txt
    echo "[GET_COHORT_RECURRENT_FUSIONS_PYENV]: Stub run finished!"
    """

}