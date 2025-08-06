// 
process PREDICT_NEOPEPTIDES_SAMPLE_LEVEL_HLAS_PVACFUSE {
    
    label 'predictSampleNeopeptides'

    container "${params.container__pvactools}"
    
    publishDir "${params.outputDir}/${sampleName}/PVACFUSE-SAMPLE-LEVEL-HLAS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(validatedAgfusionDir)
        path(cohortWideHLAList)
        path(metaDataDir) // Directory containing metadata files, including the reference proteome FASTA

    output:
        path("${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"), emit: predictedSampleSpecificNeopeptides, optional: true

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

    echo "Running pVacfuse to predict neoepitopes from validated AGFusion results..."
    echo "Path to cohort-wide HLA allotype TSV: ${cohortWideHLAList}"
    echo "Prediction mode: Sample-level: ${params.sampleLevelHLANeoPred}"
    echo "Number of cores to use: ${params.numCores * 2}"


    # extract sample-specific HLA types from the cohort-wide HLA file
    echo "Extracting sample-specific HLA types from cohort-wide HLA file..."
    SSHLA=\$(grep "${sampleName}" ${cohortWideHLAList} | awk '{print \$2}')
    echo "Sample-specific HLA types: \${SSHLA}"

    # check if SSHLA is empty
    if [ -z "\${SSHLA}" ]; then
        echo "No sample-specific HLA types found for ${sampleName}. Assigning SEAsian-prevalent HLA types..."
        SSHLA="HLA-A*11:01 HLA-A*24:02 HLA-A*02:07 HLA-A*33:03 HLA-B*46:02 HLA-B*44:03 HLA-B*40:01 HLA-C*03:04 HLA-C*01:02"
    fi

    # Run pVacFuse with sample-specific HLA types
    echo "Running pVacfuse for sample-specific prediction..."
    if pvacfuse run ${validatedAgfusionDir} ${sampleName} \${SSHLA} \
    BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL NetMHCpanEL NetMHCcons SMMPMBEC \
    "${sampleName}_sample-level-HLA-pred" \
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
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"
    echo "Stub run finished!"
    """
}
