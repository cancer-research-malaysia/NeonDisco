// Process to collect all TSV files and concatenate them
process COLLECT_COHORTWIDE_FUSION_NEOPEPTIDES_COHOHLA_PYENV {
    cpus 1
    
    label 'collectCohortFusionNeopeptidesCohoHla'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/Cohortwide-Fusions-OUT/Cohortwide_Neopeptides-OUT", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
    path cohortLevelHlaPvacFuseOutTsvs
    
    output:
    path "Cohortwide_cohort_level_HLA_fusion_neopeptides.tsv", emit: cohortwideCohortHLAFusNeopeptides

    script:
    """
    concatenate-cohortwide-neopeptides-pvacfuse--nf.py Cohortwide_cohort_level_HLA_fusion_neopeptides.tsv ${cohortLevelHlaPvacFuseOutTsvs} && echo "Cohortwide HLA Fusion Neopeptides collection finished!"
    """
    stub:
    """
    touch Cohortwide_cohort_level_HLA_fusion_neopeptides.tsv
    echo "Stub run finished!"
    """
}
