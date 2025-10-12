// Process to collect all TSV files and concatenate them
process COLLECT_COHORTWIDE_FUSION_NEOPEPTIDES_SAMPHLA_PYENV {
    errorStrategy 'retry'
    maxRetries 3
    cpus 1
    
    label 'collectCohortFusionNeopeptidesSampHla'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Fusions-OUT/Cohortwide_Neopeptides-OUT", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
    path sampleLevelHlaPvacFuseOutTsvs
    
    output:
    path "Cohortwide_sample_level_HLA_fusion_neopeptides.tsv", emit: cohortwideSampleHLAFusNeopeptides

    script:
    """
    concatenate-cohortwide-neopeptides-pvacfuse--nf.py Cohortwide_sample_level_HLA_fusion_neopeptides.tsv ${sampleLevelHlaPvacFuseOutTsvs} && echo "Cohortwide HLA Fusion Neopeptides collection finished!"
    """
    stub:
    """
    touch Cohortwide_sample_level_HLA_fusion_neopeptides.tsv
    echo "Stub run finished!"
    """
}
