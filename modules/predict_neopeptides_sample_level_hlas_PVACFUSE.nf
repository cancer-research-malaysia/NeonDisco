process PREDICT_NEOPEPTIDES_SAMPLE_LEVEL_HLAS_PVACFUSE {
    cpus params.numCores * 2
    
    label 'predictSampleNeopeptides'

    container "${params.container__pvactools}"
    
    publishDir "${params.outputDir}/${sampleName}/PVACFUSE-SAMPLE-LEVEL-HLAS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(validatedAgfusionDir)
        path(cohortWideHLAList)
        path(metaDataDir) // Directory containing metadata files, including the reference proteome FASTA

    output:
        path("${sampleName}_*-HLA-pred/MHC_Class_I/${sampleName}.fasta")
        path("${sampleName}_*-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"), emit: predictedSampleSpecificNeopeptides
        path("${sampleName}_*-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.tsv")
        path("${sampleName}_*-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.tsv"), emit: aggregatedEpitopes
        path("${sampleName}_*-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.reference_matches.tsv"), emit: referenceMatches
        path("${sampleName}-FI-validated-fusion-sample-HLA-immunogenic-peptides-13aa.fasta"), emit: specializedFasta
        path("${sampleName}_sample_HLA_pvacfuse_execution_report.txt"), emit: executionReport

    script:
    """
    # set flanking amino acid length for specialized FASTA generation
    FLANK_LENGTH=13
    # run pvacfuse with sample-level HLAs
    bash predict-neopeptides-sample-hla-pvacfuse--nf.sh ${sampleName} ${validatedAgfusionDir} ${cohortWideHLAList} ${metaDataDir} ${params.sampleHLANeoPred} \${FLANK_LENGTH} ${task.cpus}

    """

    stub:
    """
    # Create stub report
    REPORT_FILE="${sampleName}_sample_HLA_pvacfuse_execution_report.txt"
    echo "=== PVACFUSE SAMPLE-LEVEL HLA NEOPEPTIDE PREDICTION REPORT (STUB RUN) ===" > "\$REPORT_FILE"
    
    # Create stub output file - could be either sample-level or SEA-SET depending on HLA availability
    mkdir -p "${sampleName}_sample-level-HLA-pred/MHC_Class_I"
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.tsv"
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.tsv"
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.reference_matches.tsv"
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.fasta"
    touch "${sampleName}-FI-validated-fusion-sample-HLA-immunogenic-peptides-13aa.fasta"

    """
}
