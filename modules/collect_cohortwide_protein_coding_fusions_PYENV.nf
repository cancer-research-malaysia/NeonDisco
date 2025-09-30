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
    path "Cohortwide_protein-coding_fusions_manifest.txt"
    path "Cohortwide_normfiltered_protein-coding-only.tsv"
    
    script:
    """
    cat ${proteinCodingFusionsTxts} | sort -u > Cohortwide_protein-coding_fusions_manifest.txt

    # Filter using awk (faster for many patterns)
    awk 'NR==FNR {ids[\$1]=1; next} FNR==1 {print; next} \$1 in ids' Cohortwide_protein-coding_fusions_manifest.txt ${cohortwideNormFilteredFusionsFile} > Cohortwide_normfiltered_protein-coding-only.tsv

    """
    stub:
    """
    touch Cohortwide_normfiltered_protein-coding-only.tsv
    touch Cohortwide_protein-coding_fusions_manifest.txt
    echo "Stub run finished!"
    """
}

