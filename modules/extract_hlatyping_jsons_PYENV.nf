// 
process EXTRACT_HLATYPING_JSONS_PYENV {
    publishDir "${params.outputDir}/${sampleName}/HLATYPING-COHORT-WIDE-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    container "${params.container__pyenv}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name PROCESS-JSONS -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
    tuple val(sampleName), path(jsonFile)
  
    output:
    tuple val(sampleName), path("${sampleName}-HLA-types-reformatted.tsv"), emit: hlaTypingTsv
  
    script:
    """
    python /home/app/scripts/utils/extract-hla-jsons.py ${sampleName} ${jsonFile} > ${sampleName}-HLA-types-reformatted.tsv
    """
    stub:
    """
    touch ${sampleName}-HLA-types-reformatted.tsv
    echo "Stub run finished!" > test_stub_extract-json.log
    """
}
