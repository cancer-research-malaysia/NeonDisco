// 
process VALIDATE_IN_SILICO_FUSIONINSPECTOR {
    maxForks 1

    label 'validateInSilico'

    container "${params.container__starfusion}"
    
    afterScript params.deleteIntMedFiles ? "find ./ -name \"${sampleName}*_trimmed.R?.f*q.*\" -type l -exec sh -c 'rm -f \$(readlink -f \"{}\")' \\; -delete" : "echo 'Skipping intermediate file cleanup...'"
    
    publishDir "${params.outputDir}/${sampleName}/IN-SILICO-VALIDATION-FUSINS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(filteredAgfusionOutdir), path(uniqueFiltFusionPairsForFusIns), path(trimmedReads)
        path ctatDB
    
    output:
        tuple val(sampleName), path("FI/${sampleName}.FusionInspector.fusions.abridged.tsv"), emit: fusInspectorTsv

    script:
    """
    echo "Path to input file for FusionInspector gene pair input: ${uniqueFiltFusionPairsForFusIns}"
    echo "Path to filtered AGFusion output directory: ${filteredAgfusionOutdir}"
    echo "Sample name: ${sampleName}"
    echo "Trimmed reads: ${trimmedReads}"
    # extract reads from the tuple
    READ1=${trimmedReads[0]}  # First file in the nested list will be read 1 file
    READ2=${trimmedReads[1]}  # Second file in the nested list will be read 2 file
    echo "Read 1: \${READ1}"
    echo "Read 2: \${READ2}"

    
    echo "Running preprocessing script to filter for agfusion-compatible gene pairs..."
    if fusins-preproc--nf.sh ${uniqueFiltFusionPairsForFusIns} ${filteredAgfusionOutdir} ${sampleName}-genePairs-for-FusIns-filtered.txt; then
        echo "Preprocessing script has finished running."
        echo "Starting FusionInspector run with filtered gene pairs..."
        if fusins--nf.sh ${sampleName}-genePairs-for-FusIns-filtered.txt ${ctatDB} \$READ1 \$READ2 ${sampleName} ${params.numCores}; then
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
