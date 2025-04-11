//
process ALIGN_READS_2PASS_STARSAM_S3LOCAL {
    maxForks 5
    afterScript "find ./ -name '${sampleName}*.f*.gz' -type l -exec sh -c 'rm -f \$(readlink -f \"{}\")' \\; -delete"
    //publishDir "${params.outputDir}/${sampleName}/STAR-out-2P", mode: 'copy',
    //   saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name ALIGNMENT-2P -v ${params.arribaDB}:/home/app/libs -v \$(pwd):/home/app/nf_work -v ${params.binDir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(readFile1), path(readFile2)

    output:
        tuple val(sampleName), path("*-STAR*Aligned.out.bam", arity: '1'), emit: aligned_bams

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}
    CORES=${params.numCores}
    STAR_INDEX="/home/app/libs/ref_genome.fa.star.idx"

    echo "Processing files of sample \${SAMPLE_ID}"
    echo "Number of cores to use: ${params.numCores}"
    echo "The index path: \${STAR_INDEX}"
    
    # STAR 2-pass alignment
    if bash /home/app/scripts/star-2pass-nf.sh ${readFile1} ${readFile2} "\${SAMPLE_ID}" ${params.numCores} "\${STAR_INDEX}"; then
        echo "STAR sample-level 2-pass alignment is complete!"
    else
        echo "STAR alignment failed. Check logs. Exiting..."
        exit 1
    fi

    """
    stub:
    """
    touch ${sampleName}-STAR_2pass_Aligned.out.bam
    touch ${sampleName}-STAR_2pass_Log.final.out
    touch ${sampleName}-STAR_2pass_SJ.out.tab
    echo "Stub run finished!" > test_stub_STAR-2pass-align.log
    """
}
