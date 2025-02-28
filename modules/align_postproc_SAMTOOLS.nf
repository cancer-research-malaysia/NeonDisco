// align trimmed reads
process ALIGN_READS_2PASS_STARSAM {
    errorStrategy 'finish'
    maxForks 2
    publishDir "${params.output_dir}/${sampleName}/STAR-out-2pass", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    container "${params.container__preproc}"
    containerOptions "-e \"MHF_HOST_UID=\$(id -u)\" -e \"MHF_HOST_GID=\$(id -g)\" --name star-align-reads -v ${params.arriba_db}:/home/app/libs -v \$(pwd):/home/app/nf_work -v ${params.bin_dir}:/home/app/scripts"
    
    input:
        tuple val(sampleName), path(readFiles)

    output:
        tuple val(sampleName), path("*-STAR*Aligned.out.bam", arity: '1'), emit: aligned_bams
        tuple path("*_Log.final.out"), path("*_SJ.out.tab")

    script:
    """
    # variables
    READ1=${readFiles[0]}
    READ2=${readFiles[1]}
    SAMPLE_ID=${sampleName}
    CORES=${params.num_cores}
    STAR_INDEX="/home/app/libs/ref_genome.fa.star.idx"

    echo "Processing files: \${READ1} & \${READ2} of sample \${SAMPLE_ID}"
    echo "Number of cores to use: ${params.num_cores}"
    echo "The index path: \${STAR_INDEX}"
	echo "Starting STAR sample-level 2-pass alignment..."

	if bash /home/app/scripts/star-2pass-nf.sh "\${READ1}" "\${READ2}" "\${SAMPLE_ID}" ${params.num_cores} "\${STAR_INDEX}"; then
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

// Fix mate information
// process FIX_MATE {
//     container "${params.container__samtools}"
    
//     input:
//         tuple val(sampleName), path(bam)
    
//     output:
//         tuple val(sampleName), path("*_fixmate.bam"), emit: fixmate_bam
    
//     script:
//     """
//     samtools view -h -F 4 --min-MQ ${params.min_mapq} ${bam} | \
//     samtools sort -n -@ ${params.num_cores} -m 4G -O SAM - | \
//     samtools fixmate -pcmu -O bam - ${sampleName}_fixmate.bam
//     """
// }

// // Mark duplicates
// process MARK_DUPLICATES {
//     container "${params.container__samtools}"
    
//     input:
//         tuple val(sampleName), path(bam)
    
//     output:
//         tuple val(sampleName), path("*_final.bam"), emit: final_bam
    
//     script:
//     """
//     samtools sort -@ ${params.num_cores} -m 4G -O BAM ${bam} | \
//     samtools markdup -r -@ ${params.num_cores} - ${sampleName}_final.bam
//     """
// }

// // Index final BAM
// process INDEX_BAM {
//     publishDir "${params.output_dir}/${sampleName}/STAR-alignment", mode: 'copy'
//     container "${params.container__samtools}"
    
//     input:
//         tuple val(sampleName), path(bam)
    
//     output:
//         tuple val(sampleName), path("*_final.bam"), path("*_final.bam.bai"), emit: indexed_bam
    
//     script:
//     """
//     samtools index ${bam}
//     """
// }