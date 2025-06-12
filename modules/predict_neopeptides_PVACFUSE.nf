// 
process PREDICT_NEOPEPTIDES_PVACFUSE {
    maxForks 1
    
    publishDir "${params.outputDir}/${sampleName}/PVACFUSE-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pvactools}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name PREDICT_NEOPEPTIDES -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts -v ${params.metaFilesLoc}:/home/app/metadata"
    
    input:
        tuple val(sampleName), path(validatedAgfusionDir)

    output:
        path("${sampleName}_sample-specific-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"), emit: predictedNeopeptides

    script:
    """
    echo "Path to validated AGFusion directories: ${validatedAgfusionDir}"
    echo "Sample name: ${sampleName}"
    echo "Running PVACFUSE to predict neoepitopes from validated AGFusion results..."
    echo "Using sample-specific HLA alleles: ${params.sampleSpecificHLA}"
    echo "Number of cores to use: ${params.numCores}"

    if pvacfuse run ${validatedAgfusionDir} ${sampleName} ${params.sampleSpecificHLA} all "${sampleName}_sample-specific-HLA-pred" --iedb-install-directory /opt/iedb -t ${params.numCores} --allele-specific-binding-thresholds --run-reference-proteome-similarity --peptide-fasta /home/app/metadata/Homo_sapiens.GRCh38.pep.all.fa.gz; then
        echo "pVacFuse run finished!"
    else
        echo "Something went wrong."
        exit 1
    fi
    """

    stub:
    """
    touch "${sampleName}_sample-specific-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"
    echo "Stub run finished!"
    """
}
