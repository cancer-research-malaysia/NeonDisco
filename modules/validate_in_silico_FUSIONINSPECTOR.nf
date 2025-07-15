// 
process VALIDATE_IN_SILICO_FUSIONINSPECTOR {
    memory '100 GB'
    afterScript params.deleteIntMedFiles ? "find ./ -name \"${sampleName}*_trimmed.R?.f*q.*\" -type l -exec sh -c 'rm -f \$(readlink -f \"{}\")' \\; -delete" : "echo 'Skipping intermediate file cleanup...'"
    
    publishDir "${params.outputDir}/${sampleName}/IN-SILICO-VALIDATION-FUSINS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__starfusion}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name FUSION_INSPECTOR_VALIDATION -v ${params.ctatDB}:/home/refs/ctat-db"
    
    input:
        tuple val(sampleName), path(filtered_agfusion_outdir), path(uniqueFiltFusionPairsForFusIns), path(trimmedReads)
        
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
    if fusins-preproc--nf.sh ${uniqueFiltFusionPairsForFusIns} ${filtered_agfusion_outdir} ${sampleName}-genePairs-for-FusIns-filtered.txt; then
        echo "Preprocessing script has finished running."
        echo "Starting FusionInspector run with filtered gene pairs..."
        if fusins--nf.sh ${sampleName}-genePairs-for-FusIns-filtered.txt /home/refs/ctat-db \$READ1 \$READ2 ${sampleName}; then
            echo "FusionInspector has finished running."
        else
            echo "FusionInspector run failed. Please check the logs for details."
            exit 1
        fi
    fi

    """
    stub:
    """
    mkdir -p FI
    touch "FI/${sampleName}.FusionInspector.fusions.abridged.tsv"
    echo "stub run finished!" > "FI/${sampleName}.FusionInspector.fusions.abridged.tsv"
    """
}
