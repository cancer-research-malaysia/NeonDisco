// Run HLA typing module
process PREDICT_CODING_SEQ_AGFUSION {
    publishDir "${params.output_dir}/${sampleName}/FT/AGFUSION-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__agfusion}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name predict-coding-seq-agfusion -v ${params.agfusion_db}:/work/libs -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        tuple val(sampleName), path(ftFile1), path(ftFile2)

    output:
        tuple val(sampleName), path("${sampleName}_agfusion_*"), emit: agfusion_outdir

    script:
    if (ftFile2 != "./assets/NO_FILE") {
        // This block runs when ftFile2 path leads to a valid path
        // Process both ftFile1 and ftFile2
        """
        echo "Running AGFusion to predict fusion protein sequences from both Arriba and FusionCatcher..."
        if bash /work/scripts/agfusion-nf.sh ${ftFile1} ${params.num_cores} && bash /work/scripts/agfusion-nf.sh ${ftFile2} ${params.num_cores}; then
            echo "AGFusion has finished running on the outputs of Arriba and FusionCatcher of ${sampleName}."
        fi
        """
    } else {
        """
        echo "Running AGFusion to predict fusion protein sequences from the output of the selected FT calling tool..."
        if bash /work/scripts/agfusion-nf.sh ${ftFile1} ${params.num_cores}; then
            echo "AGFusion has finished running the output of the selected FT calling tool of ${sampleName}."
        fi
        """
    }
    stub:
    """
    mkdir -p ${sampleName}_agfusion_test_stub
    echo "stub run finished!" > ${sampleName}_agfusion_test_stub/stub.out
    """
}
