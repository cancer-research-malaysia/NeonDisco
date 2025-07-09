// Alternative approach: Create a separate recurrent-aware filtering step
process FILTER_VALIDATED_FUSIONS_FOR_RECURRENT_PYENV {
    
    publishDir "${params.outputDir}/${sampleName}/RECURRENT-VALIDATED-FUSIONS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FILTER-FOR-RECURRENT -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
    tuple val(sampleName), path(validatedFusions), path(validatedAgfusionDir)
    tuple val(cohortLabel), path(recurrentFusionsTsv)

    output:
    tuple val(sampleName), path("${sampleName}-validated-recurrent-fusions-only.tsv"), emit: validatedRecurrentFusions
    tuple val(sampleName), path("validated-recurrent-agfusion-outdir"), emit: validatedRecurrentAgfusionDir
    path("${sampleName}_recurrent_filter_report.txt"), emit: recurrentFilterReport

    script:
    """
    bash /home/app/scripts/filter-copy-recurrent-fusion-agfusiondir--nf.sh ${sampleName} ${validatedFusions} ${recurrentFusionsTsv} || echo "ERROR: Cannot run filtering for recurrent script for ${sampleName}. Check logs."

    """
    
    stub:
    """
    mkdir -p validated-recurrent-agfusion-outdir
    echo -e "sample_name\\tgene_pair\\tfusion_name\\tgene1\\tgene2" > ${sampleName}-validated-recurrent-fusions-only.tsv
    echo "FILTER_VALIDATED_FOR_RECURRENT_PYENV: Stub run finished!" > ${sampleName}_recurrent_filter_report.txt
    echo "STATUS: Stub run" >> ${sampleName}_recurrent_filter_report.txt
    """
    
}
