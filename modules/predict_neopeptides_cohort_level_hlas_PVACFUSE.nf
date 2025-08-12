// 
process PREDICT_NEOPEPTIDES_COHORT_LEVEL_HLAS_PVACFUSE {
    
    label 'predictCohortNeopeptides'
    
    container "${params.container__pvactools}"
    
    publishDir "${params.outputDir}/${sampleName}/PVACFUSE-COHORT-LEVEL-HLAS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(validatedAgfusionDir)
        path(cohortFivePercentFreqHLAs)
        path(metaDataDir) // Directory containing metadata files, including the reference proteome FASTA

    output:
        path("${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"), emit: predictedCohortNeopeptides, optional: true

    script:
    """
    echo "Path to validated AGFusion directories: ${validatedAgfusionDir}"
    echo "Sample name: ${sampleName}"

    # Check if directory contains actual AGFusion results (ignore .empty files)
    agfusion_dirs=\$(find ${validatedAgfusionDir} -maxdepth 1 -type d ! -name "." ! -name ".." | wc -l)
    empty_marker=\$(find ${validatedAgfusionDir} -name ".empty" -type f | wc -l)
    
    if [ \$agfusion_dirs -eq 0 ] || [ \$empty_marker -gt 0 ]; then
        echo "No validated AGFusion directories found for ${sampleName}. Skipping neopeptide prediction process gracefully."
        exit 0
    fi

    echo "Running PVACFUSE to predict neoepitopes from validated AGFusion results..."
    echo "Path to cohort-wide 5% frequency HLAs: ${cohortFivePercentFreqHLAs}"
    echo "Prediction mode: Cohort-level --> ${params.cohortHLANeoPred}"
    echo "Number of cores to use: ${params.numCores * 2}"


    # Run pVacFuse with cohort-level HLA types
    echo "Extracting HLA types from cohort-wide 5% frequency HLAs file..."
    COHORT_HLAS=\$(awk '{print \$1}' ${cohortFivePercentFreqHLAs})
    echo "Cohort-wide HLA types: \${COHORT_HLAS}"


    echo "Running pVacfuse for cohort-wide HLA binding and immunogenicity prediction..."
    if pvacfuse run ${validatedAgfusionDir} ${sampleName} \${COHORT_HLAS} \
    BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL NetMHCpanEL NetMHCcons SMMPMBEC \
    "${sampleName}_cohort-level-HLA-pred" \
    --iedb-install-directory /opt/iedb \
    --allele-specific-binding-thresholds \
    --run-reference-proteome-similarity \
    --peptide-fasta ${metaDataDir}/Homo_sapiens.GRCh38.pep.all.fa.gz \
    --netmhc-stab \
    -t ${params.numCores * 2} \
    -a sample_name; then
        echo "pVacFuse run finished!"
    else
        echo "Something went wrong."
        exit 1
    fi

    """

    stub:
    """
    touch "${sampleName}_cohort-level-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"
    echo "Stub run finished!"
    """
}
