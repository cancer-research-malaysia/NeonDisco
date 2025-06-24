// 
process PREDICT_COHORTWIDE_NEOPEPTIDES_PVACFUSE {
    maxForks 1
    
    publishDir "${params.outputDir}/${sampleName}/PVACFUSE-COHORTWIDE-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pvactools}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name PREDICT_NEOPEPTIDES_COHORTWIDE -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts -v ${params.metaFilesLoc}:/home/app/metadata"
    
    input:
        tuple val(sampleName), path(validatedAgfusionDir)
        path(cohortFivePercentFreqHLAs)

    output:
        path("${sampleName}_cohortwide-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"), emit: predictedCohortNeopeptides

    script:
    """
    echo "Path to validated AGFusion directories: ${validatedAgfusionDir}"
    echo "Sample name: ${sampleName}"
    echo "Running PVACFUSE to predict neoepitopes from validated AGFusion results..."
    echo "Path to cohort-wide 5% frequency HLAs: ${cohortFivePercentFreqHLAs}"
    echo "Prediction parameters: Cohort-wide: ${params.cohortWideNeoPeptidePrediction}"
    echo "Number of cores to use: ${params.numCores}"


    # Run pVacFuse with sample-specific HLA types
    echo "Extracting HLA types from cohort-wide 5% frequency HLAs file..."
    COHORT_HLAS=\$(awk '{print \$1}' ${cohortFivePercentFreqHLAs})
    echo "Cohort-wide HLA types: testcommand --cohort \${COHORT_HLAS}"
    # echo "Running pVacfuse for cohort-wide HLA binding and immunogenicity prediction..."
    # if pvacfuse run ${validatedAgfusionDir} ${sampleName} \${SSHLA} all "${sampleName}_cohortwide-HLA-pred" --iedb-install-directory /opt/iedb -t ${params.numCores} --allele-specific-binding-thresholds --run-reference-proteome-similarity --peptide-fasta /home/app/metadata/Homo_sapiens.GRCh38.pep.all.fa.gz; then
    #    echo "pVacFuse run finished!"
    #else
    #    echo "Something went wrong."
    #    exit 1
    #fi

    """

    stub:
    """
    touch "${sampleName}_cohortwide-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"
    echo "Stub run finished!"
    """
}
