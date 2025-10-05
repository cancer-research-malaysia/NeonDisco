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
    path "Unique-fusionTranscriptID-only/Cohortwide_normfiltered_fusions-UNIQUE.manifest.txt"
    
    script:
    """
    mkdir -p Unique-fusionTranscriptID-only && concatenate-cohortwide-fusions--nf.py --input_files ${normFilteredFusionsTsvs} --output Unique-fusionTranscriptID-only/Cohortwide_normfiltered_fusions.tsv
    """
    stub:
    """
    touch Cohortwide_normfiltered_fusions.tsv
    mkdir -p Unique-fusionTranscriptID-only
    touch Unique-fusionTranscriptID-only/Cohortwide_normfiltered_fusions-UNIQUE.manifest.txt
    echo "Stub run finished!"
    """
}

