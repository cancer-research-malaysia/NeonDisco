// Process to collect all TSV files and concatenate them
process COLLECT_COHORTWIDE_VALIDATED_FUSIONS_PYENV {
    cpus 1
    
    label 'collectCohortValidatedFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Fusions-OUT", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
    path validatedFusionsTsvs
    
    output:
    path "Cohortwide_normfiltered_protein-coding_FI-validated_fusions.tsv", emit: cohortwideValidatedFusionsFile
    path "Unique-fusionTranscriptID-only/Cohortwide_normfiltered_protein-coding_FI-validated_fusions-UNIQUE.manifest.txt"

    script:
    """
    concatenate-cohortwide-fusions--nf.py --input_files ${validatedFusionsTsvs} --output Cohortwide_normfiltered_protein-coding_FI-validated_fusions.tsv && mkdir -p Unique-fusionTranscriptID-only && mv Cohortwide_normfiltered_protein-coding_FI-validated_fusions-UNIQUE.manifest.txt Unique-fusionTranscriptID-only/

    """
    stub:
    """
    touch Cohortwide_normfiltered_protein-coding_FI-validated_fusions.tsv
    mkdir -p Unique-fusionTranscriptID-only
    touch Unique-fusionTranscriptID-only/Cohortwide_normfiltered_protein-coding_FI-validated_fusions-UNIQUE.manifest.txt
    echo "Stub run finished!"
    """
}
