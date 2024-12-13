// Run HLA typing module
process TYPE_HLA_ALLELES {
    //publishDir "${params.output_dir}/${sampleName}", mode: 'copy'
    container "${params.container__hlahd}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name hla-typing -v \$(pwd):/work/nf_work -v ${params.bin_dir}:/work/scripts"
    
    input:
    tuple val(sampleName), path(bamFile)

    script:
    """
    echo "Path to BAM input file (WES) for sample ${sampleName}: ${bamFile}"
    if bash /work/scripts/hlahd-nf.sh ${sampleName} ${bamFile} 12; then
        echo "HLA-HD has finished running on ${sampleName}."
    fi
    """
}
