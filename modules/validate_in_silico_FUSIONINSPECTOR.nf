// 
process VALIDATE_IN_SILICO_FUSIONINSPECTOR {
    
    publishDir "${params.outputDir}/${sampleName}/IN-SILICO-VALIDATION-FUSINS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__starfusion}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name PREP-FUSIONS-FOR-FUSINS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        path(agfusion_outdir)
        tuple val(sampleName), path(uniqueFiltFusionPairsForFusIns)
        

    // output:
    //     tuple val(sampleName), path("${sampleName}-collated-FT-normFiltered-unique-genePairs-for-FusIns.txt"), emit: uniqueFiltFusionPairsForFusIns

    script:
    """
    echo "Path to input file for FusionInspector gene pair input: ${uniqueFiltFusionPairsForFusIns}"
    echo "Path to filtered AGFusion output directory: ${agfusion_outdir}"
    echo "Sample name: ${sampleName}"

    echo "Running preprocessing script to filter for agfusion-compatible gene pairs..."
    if bash /home/app/scripts/fusins-preproc--nf.sh ${uniqueFiltFusionPairsForFusIns} ${agfusion_outdir} ${sampleName}-genePairs-for-FusIns-filtered.txt; then
        echo "Preprocessing script has finished running."
    fi

    """
    stub:
    """
    echo "stub run finished!" > ${sampleName}-genePairs-for-FusIns-filtered.txt
    """
}
