// 
process COLLECT_COHORTWIDE_UNFILTERED_FUSIONS_PYENV {
    cpus 1
    
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
    concatenate-cohortwide-fusions-unfilt--nf.py --input_files ${wrangledFusionsTsvs} --output Cohortwide_all_unfiltered_fusions.tsv && mkdir -p Unique-fusionTranscriptID-only && mv Cohortwide_all_unfiltered_fusions-UNIQUE.manifest.txt Unique-fusionTranscriptID-only/ && echo "Cohort-wide unfiltered fusion collection complete!"

    """
    stub:
    """
    touch Cohortwide_all_unfiltered_fusions.tsv
    mkdir -p Unique-fusionTranscriptID-only
    touch Unique-fusionTranscriptID-only/Cohortwide_all_unfiltered_fusions-UNIQUE.manifest.txt
    echo "stub run finished!"
    """
}
