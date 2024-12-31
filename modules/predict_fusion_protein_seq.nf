// Run HLA typing module
process PREDICT_FUSION_PROTEIN_SEQ
    publishDir "${params.output_dir}/${sampleName}/FT/AGFusion-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__agfusion}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name predict-protein-agfusion -v ${params.agfusion_db}:/work/libs -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
        
        val(numCores)

    output:
        

    script:
    """
    echo "Running AGFusion to predict fusion protein sequences..."
    """
    stub:
    """
    echo "stub run finished!"
    """
}
