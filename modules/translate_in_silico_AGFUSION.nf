// 
process TRANSLATE_IN_SILICO_AGFUSION {
    
    publishDir "${params.outputDir}/${sampleName}/AGFUSION-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__agfusion}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name TRANSLATE-IN-SILICO -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(filteredFusions)

    output:
        tuple val(sampleName), path("agfusion-dirs/"), emit: agfusion_outdir

    script:
    """
    echo "Path to filtered fusion transcripts: ${filteredFusions}"
    echo "Running parser to extract fusion transcripts from the filtered lists of fusion transcripts from the fusion calling tools..."
    if python /home/app/scripts/generate-agfusion-cmd--nf.py -i ${filteredFusions} -t -c; then
        echo "Parser has finished running the output of the selected FT calling tool of ${sampleName}."
    fi

    echo "Running AGFusion..."
    if bash agfusion-cmd.sh; then
        echo "AGFusion has finished."
    fi
    """
    stub:
    """
    mkdir -p agfusion-dirs
    echo "stub run finished!" > agfusion-dirs/stub.out
    """
}
