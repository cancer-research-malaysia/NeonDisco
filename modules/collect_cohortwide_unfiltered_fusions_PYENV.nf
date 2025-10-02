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
    path "Cohortwide_all_unfiltered_fusions.tsv", emit: cohortwideUnfilteredFusionsFile
    path "Unique-found-fusions__OUT/Cohortwide_all_unfiltered_fusions-UNIQUE.manifest.txt"

    script:
    """
    echo "Running Python script to collect unfiltered cohort-wide FT files..."
    concatenate-cohortwide-fusions-unfilt--nf.py \\
        --input_files ${wrangledFusionsTsvs} \\
        --output Cohortwide_all_unfiltered_fusions.tsv && \\
        mkdir -p Unique-found-fusions__OUT && \\
        mv Cohortwide_all_unfiltered_fusions-UNIQUE.manifest.txt Unique-found-fusions__OUT/
    """

    stub:
    """
    touch Cohortwide_all_unfiltered_fusions.tsv
    mkdir -p Unique-found-fusions__OUT
    touch Unique-found-fusions__OUT/Cohortwide_all_unfiltered_fusions-UNIQUE.manifest.txt
    echo "stub run finished!\thello my world!" > Cohortwide_all_unfiltered_fusions.tsv
    """
}
