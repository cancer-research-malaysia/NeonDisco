// Process to collect all TSV files and concatenate them
process COLLECT_COHORTWIDE_PROTEIN_CODING_FUSIONS_PYENV {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'collectProteinCodingFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Fusions-OUT", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
    path proteinCodingFusionsTxts
    path cohortwideNormFilteredFusionsFile
    
    output:
    path "Cohortwide_normfiltered_protein-coding-only_fusions.tsv"
    path "Unique-fusionTranscriptID-only/Cohortwide_normfiltered_protein-coding-only_fusions-UNIQUE.manifest.txt"
    
    script:
    """
    mkdir -p Unique-fusionTranscriptID-only
    cat ${proteinCodingFusionsTxts} | sort -u > Unique-fusionTranscriptID-only/Cohortwide_normfiltered_protein-coding-only_fusions-UNIQUE.manifest.txt

    # Filter using awk (faster for many patterns)
    awk 'NR==FNR {ids[\$1]=1; next} FNR==1 {print; next} \$1 in ids' Unique-fusionTranscriptID-only/Cohortwide_normfiltered_protein-coding-only_fusions-UNIQUE.manifest.txt ${cohortwideNormFilteredFusionsFile} > Cohortwide_normfiltered_protein-coding-only_fusions.tsv

    """
    stub:
    """
    touch Cohortwide_normfiltered_protein-coding-only_fusions.tsv
    mkdir -p Unique-fusionTranscriptID-only
    touch Unique-fusionTranscriptID-only/Cohortwide_normfiltered_protein-coding-only_fusions-UNIQUE.manifest.txt
    echo "Stub run finished!"
    """
}

