// Run HLA typing module
process PREDICT_CODING_SEQ_AGFUSION {
    
    publishDir "${params.outputDir}/${sampleName}/AGFUSION-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    container "${params.container__agfusion}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name AGFUSION-CODSEQ-PREDICTION -v ${params.agfusionDB}:/work/libs -v \$(pwd):/work/nf_work -v ${params.binDir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(collatedFTList)

    // output:
    //     tuple val(sampleName), path("${sampleName}_agfusion_*"), emit: agfusion_outdir

    script:
    """
    echo "Path to collated fusion transcripts: ${collatedFTList}"
    echo "Running parser to extract fusion transcripts from the collated raw lists of fusion transcripts from Arriba and Fusioncatcher..."
    if python /work/scripts/parse-collation-4-agf-nf.py -i ${collatedFTList} -t -c "./agfusion-cmd.sh" -o "./"; then
        echo "Parser has finished running the output of the selected FT calling tool of ${sampleName}."
    fi

    echo "Running AGFusion to predict fusion protein sequences from the collated raw lists of fusion transcripts from Arriba and Fusioncatcher..."
    if bash agfusion-cmd.sh; then
        echo "AGFusion has finished running on the collated list of fusion transcripts of ${sampleName} from Arriba and Fusioncatcher."
    fi
    """
    stub:
    """
    mkdir -p ${sampleName}_agfusion_test_stub
    echo "stub run finished!" > ${sampleName}_agfusion_test_stub/stub.out
    """
}
