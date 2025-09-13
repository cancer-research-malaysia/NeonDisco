// 
process COLLECT_COHORTWIDE_UNFILTERED_FUSIONS_PYENV {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'collectCohortUnfilteredFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Fusions-OUT", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
    path wrangledFusionsTsvs
    
    output:
    path "Cohortwide_raw_unfiltered_fusions.tsv", emit: cohortwideUnfilteredFusionsFile


    script:
    """
    echo "Running Python script to collect unfiltered cohort-wide FT files..."
    concatenate-cohortwide-fusions-RAW--nf.py \\
        --input_files ${wrangledFusionsTsvs} \\
        --output Cohortwide_raw_unfiltered_fusions.tsv
    """

    stub:
    """
    touch Cohortwide_raw_unfiltered_fusions.tsv
    echo "stub run finished!\thello my world!" > Cohortwide_raw_unfiltered_fusions.tsv
    """
}
