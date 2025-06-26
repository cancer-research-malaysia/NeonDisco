// 
process VALIDATE_IN_SILICO_FUSIONINSPECTOR {
    
    publishDir "${params.outputDir}/${sampleName}/IN-SILICO-VALIDATION-FUSINS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__starfusion}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FUSION_INSPECTOR_VALIDATION -v \$(pwd):/home/app/nf_work -v ${params.ctatDB}:/home/refs/ctat-db -v ${params.binDir}:/home/app/scripts"
    
    input:
        path(filtered_agfusion_outdir)
        tuple val(sampleName), path(uniqueFiltFusionPairsForFusIns)
        tuple val(dummyVar), path(trimmedReads)
        

    output:
        tuple val(sampleName), path("FI/${sampleName}.FusionInspector.fusions.abridged.tsv"), emit: fusInspectorTsv

    script:
    """
    echo "Path to input file for FusionInspector gene pair input: ${uniqueFiltFusionPairsForFusIns}"
    echo "Path to filtered AGFusion output directory: ${filtered_agfusion_outdir}"
    echo "Sample name: ${sampleName}"
    echo "Trimmed reads: ${trimmedReads}"
    # extract reads from the tuple
    READ1=${trimmedReads[0]}  # First file in the nested list will be read 1 file
    READ2=${trimmedReads[1]}  # Second file in the nested list will be read 2 file
    echo "Read 1: \${READ1}"
    echo "Read 2: \${READ2}"

    
    echo "Running preprocessing script to filter for agfusion-compatible gene pairs..."
    if bash /home/app/scripts/fusins-preproc--nf.sh ${uniqueFiltFusionPairsForFusIns} ${filtered_agfusion_outdir} ${sampleName}-genePairs-for-FusIns-filtered.txt; then
        echo "Preprocessing script has finished running."
        echo "Starting FusionInspector run with filtered gene pairs..."
        if bash /home/app/scripts/fusins--nf.sh ${sampleName}-genePairs-for-FusIns-filtered.txt /home/refs/ctat-db \$READ1 \$READ2 ${sampleName}; then
            echo "FusionInspector has finished running."
        else
            echo "FusionInspector run failed. Please check the logs for details."
            exit 1
        fi
    fi

    """
    stub:
    """
    echo "stub run finished!" > ${sampleName}-genePairs-for-FusIns-filtered.txt
    """
}
