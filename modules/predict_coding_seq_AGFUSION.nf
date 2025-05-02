// 
process PREDICT_CODING_SEQ_AGFUSION {
    
    publishDir "${params.outputDir}/${sampleName}/AGFUSION-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__agfusion}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name AGFUSION-CODSEQ-PREDICTION -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(filteredFusionList)

    output:
        tuple val(sampleName), path("agfusion-OUT/**"), emit: agfusion_outdir

    script:
    """
    echo "Path to collated fusion transcripts: ${filteredFusionList}"
    echo "Running parser to extract fusion transcripts from the collated raw lists of fusion transcripts from Arriba and Fusioncatcher..."
    if python /home/app/scripts/parse-inp-and-generate-agf-cmd.py -i ${filteredFusionList} -t -c; then
        echo "Parser has finished running the output of the selected FT calling tool of ${sampleName}."
    fi

    echo "Running AGFusion to predict fusion protein sequences from the collated raw lists of fusion transcripts from Arriba and Fusioncatcher..."
    if bash agfusion-cmd.sh; then
        echo "AGFusion has finished running on the collated list of fusion transcripts of ${sampleName} from Arriba and Fusioncatcher."
    fi
    """
    stub:
    """
    mkdir -p agfusion-OUT
    echo "stub run finished!" > agfusion-OUT/stub.out
    """
}
