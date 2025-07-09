// 
process PREDICT_COHORTWIDE_NEOPEPTIDES_PVACFUSE {
    
    publishDir "${params.outputDir}/${sampleName}/PVACFUSE-COHORTWIDE-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pvactools}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name PREDICT_NEOPEPTIDES_COHORTWIDE -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts -v ${params.metaDataDir}:/home/app/metadata"
    
    input:
        tuple val(sampleName), path(validatedAgfusionDir)
        path(cohortFivePercentFreqHLAs)

    output:
        path("${sampleName}_cohortwide-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"), emit: predictedCohortNeopeptides, optional: true

    script:
    """
    echo "Path to validated AGFusion directories: ${validatedAgfusionDir}"
    echo "Sample name: ${sampleName}"

    # Check if the validated AGFusion directory is not empty; if empty, exit this process gracefully
    if [ -z "\$(ls -A ${validatedAgfusionDir})" ]; then
        echo "No validated AGFusion directories found for ${sampleName}. Skipping neopeptide prediction process gracefully."
        exit 0
    fi

    echo "Running PVACFUSE to predict neoepitopes from validated AGFusion results..."
    echo "Path to cohort-wide 5% frequency HLAs: ${cohortFivePercentFreqHLAs}"
    echo "Prediction parameters: Cohort-wide: ${params.cohortWideNeoPred}"
    echo "Number of cores to use: ${params.numCores * 3}"


    # Run pVacFuse with sample-specific HLA types
    echo "Extracting HLA types from cohort-wide 5% frequency HLAs file..."
    COHORT_HLAS=\$(awk '{print \$1}' ${cohortFivePercentFreqHLAs})
    echo "Cohort-wide HLA types: testcommand --cohort \${COHORT_HLAS}"


    echo "Running pVacfuse for cohort-wide HLA binding and immunogenicity prediction..."
    if pvacfuse run ${validatedAgfusionDir} ${sampleName} \${COHORT_HLAS} BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL NetMHCpanEL NetMHCcons SMMPMBEC "${sampleName}_cohortwide-HLA-pred" --iedb-install-directory /opt/iedb -t ${params.numCores * 3} --allele-specific-binding-thresholds --run-reference-proteome-similarity --peptide-fasta /home/app/metadata/Homo_sapiens.GRCh38.pep.all.fa.gz --netmhc-stab -a sample_name; then
        echo "pVacFuse run finished!"
    else
        echo "Something went wrong."
        exit 1
    fi

    """

    stub:
    """
    touch "${sampleName}_cohortwide-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"
    echo "Stub run finished!"
    """
}
