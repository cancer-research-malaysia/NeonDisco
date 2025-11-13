process VALIDATE_IN_SILICO_FUSIONINSPECTOR {
    //maxForks 1
    cpus params.numCores

    label 'validateInSilico'

    container "${params.container__starfusion}"
    
    afterScript params.deleteIntMedFiles ? "find ./ -name \"${sampleName}*_trimmed.R?.f*q.*\" -type l -exec sh -c 'rm -f \$(readlink -f \"{}\")' \\; -delete; rm -rf ${ctatDB}" : 
    "echo 'Skipping intermediate file cleanup...'; rm -rf ${ctatDB}"
    
    publishDir "${params.outputDir}/${sampleName}/IN-SILICO-VALIDATION-FUSINS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(filteredAgfusionTarball), path(uniqueFiltFusionPairsForFusIns), path(trimmedReads)
        path ctatDB
    
    output:
        tuple val(sampleName), path("FI/${sampleName}.FusionInspector.fusions.abridged.tsv"), optional: true, emit: fusInspectorTsv
        path("${sampleName}.FAILED"), optional: true

    script:
    """
    echo "=== VALIDATE_IN_SILICO_FUSIONINSPECTOR Process Log ==="
    echo "Input tarball of filtered AGFusion output: ${filteredAgfusionTarball}"
     # extract the filtered-agfusion-dirs.tar.gz to get the filtered-agfusion-dirs directory
    tar -xzf ${filteredAgfusionTarball}
    echo "Path to input file for FusionInspector gene pair input: ${uniqueFiltFusionPairsForFusIns}"
    echo "Path to filtered AGFusion output directory: filtered-agfusion-dirs"
    echo "Sample name: ${sampleName}"
    echo "Trimmed reads: ${trimmedReads}"
    # extract reads from the tuple
    READ1=${trimmedReads[0]}  # First file in the nested list will be read 1 file
    READ2=${trimmedReads[1]}  # Second file in the nested list will be read 2 file
    echo "Read 1: \${READ1}"
    echo "Read 2: \${READ2}"

    
    echo "Running preprocessing script to filter for agfusion-compatible gene pairs..."
    if fusins-preproc--nf.sh ${uniqueFiltFusionPairsForFusIns} filtered-agfusion-dirs ${sampleName}-genePairs-for-FusIns-filtered.txt; then
        echo "Preprocessing script has finished running."
        echo "Starting FusionInspector run with filtered gene pairs..."
        if fusins--nf.sh ${sampleName}-genePairs-for-FusIns-filtered.txt ${ctatDB} \$READ1 \$READ2 ${sampleName} ${task.cpus}; then
            echo "FusionInspector has finished running."
        else
            echo "FusionInspector run failed. Please check the logs for details."
            echo "FAILED: ${sampleName} at \$(date)" > ${sampleName}.FAILED
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
