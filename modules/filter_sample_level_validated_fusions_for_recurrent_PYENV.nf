// Alternative approach: Create a separate recurrent-aware filtering step
process FILTER_SAMPLE_LEVEL_VALIDATED_FUSIONS_FOR_RECURRENT_PYENV {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'getSampleRecurrentFusions'

    container "${params.container__pyenv}"
    
    publishDir "${params.outputDir}/${sampleName}/RECURRENT-VALIDATED-FUSIONS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
    tuple val(sampleName), path(validatedFusions), path(validatedAgfusionDir)
    tuple val(cohortLabel), path(recurrentFusionsTsv)

    output:
    tuple val(sampleName), path("${sampleName}-FI-validated-recurrent-fusions-only.tsv"), emit: validatedRecurrentFusions
    tuple val(sampleName), path("FI-validated-recurrent-agfusion-dirs"), emit: validatedRecurrentAgfusionDir
    path("${sampleName}_recurrent_filter_report.txt"), emit: recurrentFilterReport

    script:
    """
    filter-copy-recurrent-fusion-agfusiondir--nf.sh ${sampleName} ${validatedFusions} ${validatedAgfusionDir} ${recurrentFusionsTsv} || echo "ERROR: Cannot run filtering for recurrent script for ${sampleName}. Check logs."

    """
    
    stub:
    """
    mkdir -p FI-validated-recurrent-agfusion-dirs
    echo -e "sample_name\\tgene_pair\\tfusion_name\\tgene1\\tgene2" > ${sampleName}-validated-recurrent-fusions-only.tsv
    echo "FILTER_VALIDATED_FOR_RECURRENT_PYENV: Stub run finished!" > ${sampleName}_recurrent_filter_report.txt
    echo "STATUS: Stub run" >> ${sampleName}_recurrent_filter_report.txt
    """
    
}
