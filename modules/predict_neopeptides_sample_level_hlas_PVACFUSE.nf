// 
process PREDICT_NEOPEPTIDES_SAMPLE_LEVEL_HLAS_PVACFUSE {
    
    label 'predictSampleNeopeptides'

    container "${params.container__pvactools}"
    
    publishDir "${params.outputDir}/${sampleName}/PVACFUSE-SAMPLE-LEVEL-HLAS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(validatedAgfusionDir)
        path(cohortWideHLAList)

    output:
        path("${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"), emit: predictedSampleSpecificNeopeptides, optional: true

    script:
    """
    echo "Path to validated AGFusion directories: ${validatedAgfusionDir}"
    echo "Sample name: ${sampleName}"

    # Check if the validated AGFusion directory is not empty; if empty, exit this process gracefully
    if [ -z "\$(ls -A ${validatedAgfusionDir})" ]; then
        echo "No validated AGFusion directories found for ${sampleName}. Skipping neopeptide prediction process gracefully."
        exit 0
    fi

    echo "Running pVacfuse to predict neoepitopes from validated AGFusion results..."
    echo "Path to cohort-wide HLA allotype TSV: ${cohortWideHLAList}"
    echo "Prediction mode: Sample-level: ${params.sampleLevelHLANeoPred}"
    echo "Number of cores to use: ${params.numCores * 3}"


    # extract sample-specific HLA types from the cohort-wide HLA file
    echo "Extracting sample-specific HLA types from cohort-wide HLA file..."
    SSHLA=\$(grep "${sampleName}" ${cohortWideHLAList} | awk '{print \$2}')
    echo "Sample-specific HLA types: \${SSHLA}"

    # check if SSHLA is empty
    if [ -z "\${SSHLA}" ]; then
        echo "No sample-specific HLA types found for ${sampleName}. Exiting."
        exit 0
    fi

    # Run pVacFuse with sample-specific HLA types
    echo "Running pVacfuse for sample-specific prediction..."
    if pvacfuse run ${validatedAgfusionDir} ${sampleName} \${SSHLA} BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL NetMHCpanEL NetMHCcons SMMPMBEC "${sampleName}_sample-level-HLA-pred" --iedb-install-directory /opt/iedb -t ${params.numCores * 3} --allele-specific-binding-thresholds --run-reference-proteome-similarity --peptide-fasta /tmp/metadata/Homo_sapiens.GRCh38.pep.all.fa.gz --netmhc-stab -a sample_name; then
        echo "pVacFuse run finished!"
    else
        echo "Something went wrong."
        exit 1
    fi

    """

    stub:
    """
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"
    echo "Stub run finished!"
    """
}
