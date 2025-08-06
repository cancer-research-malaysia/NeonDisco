// 
process REFORMAT_HLA_TYPES_PYENV {

    label 'reformatHLAs'

    container "${params.container__pyenv}"
    
    input:
    tuple val(sampleName), path(jsonFile)
  
    output:
    tuple val(sampleName), path("${sampleName}-HLA-types-reformatted.tsv"), emit: hlaTypingTsv
  
    script:
    """
    reformat-hlas--nf.py ${sampleName} ${jsonFile} > ${sampleName}-HLA-types-reformatted.tsv
    """
    stub:
    """
    touch ${sampleName}-HLA-types-reformatted.tsv
    echo "Stub run finished!" > test_stub_reformat-hlas.log
    """
}
