// 
process REFORMAT_HLA_TYPES_PYENV {

    container "${params.container__pyenv}"
    containerOptions "--rm -e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name REFORMAT-HLAS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
    tuple val(sampleName), path(jsonFile)
  
    output:
    tuple val(sampleName), path("${sampleName}-HLA-types-reformatted.tsv"), emit: hlaTypingTsv
  
    script:
    """
    python /home/app/scripts/utils/reformat-hlas--nf.py ${sampleName} ${jsonFile} > ${sampleName}-HLA-types-reformatted.tsv
    """
    stub:
    """
    touch ${sampleName}-HLA-types-reformatted.tsv
    echo "Stub run finished!" > test_stub_reformat-hlas.log
    """
}
