process PREDICT_NEOPEPTIDES_COHORT_LEVEL_HLAS_PVACFUSE {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'predictCohortNeopeptides'
    
    container "${params.container__pvactools}"
    
    publishDir "${params.outputDir}/${sampleName}/PVACFUSE-COHORT-LEVEL-HLAS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(validatedAgfusionDir)
        path(cohortFivePercentFreqHLAs)
        path(metaDataDir) // Directory containing metadata files, including the reference proteome FASTA

    output:
        path("${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.fasta")
        path("${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"), emit: predictedCohortNeopeptides
        path("${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.tsv")
        path("${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.tsv"), emit: aggregatedEpitopes
        path("${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.reference_matches.tsv"), emit: referenceMatches
        path("${sampleName}-FI-validated-fusion-cohort-HLA-immunogenic-peptides-13aa.fasta"), emit: specializedFasta
        path("${sampleName}_cohort_HLA_pvacfuse_execution_report.txt"), emit: executionReport

    script:
    """
    # set flanking amino acid length for specialized FASTA generation
    FLANK_LENGTH=13
    # run pvacfuse with cohort-level HLAs
    bash predict-neopeptides-cohort-hla-pvacfuse--nf.sh ${sampleName} ${validatedAgfusionDir} ${cohortFivePercentFreqHLAs} ${metaDataDir} ${params.sharedHLANeoPred} "\$FLANK_LENGTH" ${params.numCores * 2}
    
    """

    stub:
    """
    # Create stub report
    REPORT_FILE="${sampleName}_cohort_HLA_pvacfuse_execution_report.txt"
    echo "=== PVACFUSE COHORT-LEVEL HLA NEOPEPTIDE PREDICTION REPORT ===" > "\$REPORT_FILE"

    # Create stub output file
    mkdir -p "${sampleName}_cohort-level-HLA-pred/MHC_Class_I"
    touch "${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"
    touch "${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.tsv"
    touch "${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.tsv"
    touch "${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.reference_matches.tsv"
    touch "${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.fasta"
    touch "${sampleName}-FI-validated-fusion-cohort-HLA-immunogenic-peptides-13aa.fasta"

    """
}
