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
    path "Unique-fusionTranscriptID-only/Cohortwide_all_unfiltered_fusions-UNIQUE.manifest.txt"

    script:
    """
    echo "Running Python script to collect unfiltered cohort-wide FT files..."
    mkdir -p Unique-fusionTranscriptID-only && concatenate-cohortwide-fusions-unfilt--nf.py --input_files ${wrangledFusionsTsvs} --output Unique-fusionTranscriptID-only/Cohortwide_all_unfiltered_fusions.tsv

    """
    stub:
    """
    touch Cohortwide_all_unfiltered_fusions.tsv
    mkdir -p Unique-fusionTranscriptID-only
    touch Unique-fusionTranscriptID-only/Cohortwide_all_unfiltered_fusions-UNIQUE.manifest.txt
    echo "stub run finished!"
    """
}
