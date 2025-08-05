// type HLA allotypes using arcasHLA
// process TYPE_HLA_ALLELES_ARCASHLA {
    
//     label 'typeHLAs'

//     container "${params.container__arcashla}"
    
//     publishDir "${params.outputDir}/${sampleName}/HLA-TYPING-out", mode: 'copy',
//         saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }

    
//     input:
//         tuple val(sampleName), path(bam), path(bamIdx)

//     output:
//         tuple val(sampleName), path("*.genotype.json"), emit: allotype_json
    
//     script:
//     """
//     # Initialize variables
//     SAMPLE_ID=${sampleName}
//     BAM=${bam}
//     BAM_IDX=${bamIdx}
    
//     echo "Processing files: \${BAM} of sample \${SAMPLE_ID}"
//     echo "The index file is: \${BAM_IDX}"
//     echo "Number of cores to use: ${params.numCores}"
//     echo "Starting arcasHLA typing..."

//     # Running arcasHLA
//     if bash arcasHLA--nf.sh "\${SAMPLE_ID}" "\${BAM}" "${params.numCores}"; then
//         echo "arcasHLA typing finished successfully!"
//     else
//         echo "arcasHLA typing failed. Check logs. Exiting..."
//         exit 1
//     fi
//     """
//     stub:
//     """
//     touch ${sampleName}_test-stub.genotype.json

//     echo "Stub run finished!" > test_stub_type-hlas.log
//     """
// }

// Single process approach - arcasHLA with HLA-HD fallback in one instance

process TYPE_HLA_WITH_FALLBACK {
    
    label 'typeHLAs'
    
    container "${params.container__hlatyping}"  // Container with both tools
    
    publishDir "${params.outputDir}/${sampleName}/HLA-TYPING-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    input:
    tuple val(sampleName), path(bam), path(bamIdx)

    output:
    tuple val(sampleName), path("${sampleName}.genotype.json"), emit: hla_json
    path "${sampleName}_hla_typing.log", emit: logs
    
    script:
    """
    # Initialize variables and logging
    SAMPLE_ID=${sampleName}
    BAM=${bam}
    BAM_IDX=${bamIdx}
    LOG_FILE=${sampleName}_hla_typing.log
    
    echo "=== HLA Typing for \${SAMPLE_ID} ===" > \$LOG_FILE
    echo "Processing files: \${BAM}" >> \$LOG_FILE
    echo "Index file: \${BAM_IDX}" >> \$LOG_FILE
    echo "Number of cores: ${params.numCores}" >> \$LOG_FILE
    echo "Timestamp: \$(date)" >> \$LOG_FILE
    echo "" >> \$LOG_FILE
    
    # Try arcasHLA first
    echo "=== ATTEMPTING ARCASHLA ===" >> \$LOG_FILE
    echo "Starting arcasHLA typing..." >> \$LOG_FILE
    
    if bash arcasHLA--nf.sh "\${SAMPLE_ID}" "\${BAM}" "${params.numCores}" 2>&1 | tee -a \$LOG_FILE; then
        echo "arcasHLA command completed successfully" >> \$LOG_FILE
        
        # Check if we got meaningful results
        if [[ -f "\${SAMPLE_ID}.genotype.json" ]]; then
            # Count actual HLA alleles in the JSON
            allele_count=\$(grep -o '"[ABC][*][0-9]' "\${SAMPLE_ID}.genotype.json" | wc -l || echo "0")
            
            if [[ \$allele_count -gt 0 ]]; then
                echo "SUCCESS: arcasHLA found \${allele_count} HLA alleles for \${SAMPLE_ID}" >> \$LOG_FILE
                echo "Method used: arcasHLA" >> \$LOG_FILE
                echo "=== WORKFLOW COMPLETE ===" >> \$LOG_FILE
                exit 0  # Success - exit early
            else
                echo "WARNING: arcasHLA completed but found no HLA alleles" >> \$LOG_FILE
                echo "JSON content preview:" >> \$LOG_FILE
                head -5 "\${SAMPLE_ID}.genotype.json" >> \$LOG_FILE || echo "Could not read JSON" >> \$LOG_FILE
            fi
        else
            echo "WARNING: arcasHLA did not produce expected JSON output" >> \$LOG_FILE
        fi
    else
        echo "WARNING: arcasHLA command failed" >> \$LOG_FILE
    fi
    
    echo "" >> \$LOG_FILE
    echo "=== FALLING BACK TO HLA-HD ===" >> \$LOG_FILE
    echo "arcasHLA was unsuccessful, trying HLA-HD..." >> \$LOG_FILE
    
    # Convert BAM to FASTQ for HLA-HD
    echo "Converting BAM to FASTQ..." >> \$LOG_FILE
    samtools fastq -1 \${SAMPLE_ID}_R1.fastq -2 \${SAMPLE_ID}_R2.fastq \${BAM} 2>&1 | tee -a \$LOG_FILE
    
    if [[ ! -f "\${SAMPLE_ID}_R1.fastq" || ! -f "\${SAMPLE_ID}_R2.fastq" ]]; then
        echo "ERROR: Failed to convert BAM to FASTQ" >> \$LOG_FILE
        echo "Creating empty JSON output..." >> \$LOG_FILE
        echo '{"A": {}, "B": {}, "C": {}, "method": "failed", "sample": "'"\${SAMPLE_ID}"'"}' > \${SAMPLE_ID}.genotype.json
        exit 0
    fi
    
    # Run HLA-HD
    echo "Starting HLA-HD typing..." >> \$LOG_FILE
    if hlahd.sh -t \${task.cpus} -m 100 -f \${PWD}/freq_data \\
        \${SAMPLE_ID}_R1.fastq \${SAMPLE_ID}_R2.fastq \\
        gene_split_filt \${PWD}/HLA_gene.split.txt \${SAMPLE_ID} \\
        \${PWD}/result/\${SAMPLE_ID} 2>&1 | tee -a \$LOG_FILE; then
        
        echo "HLA-HD command completed" >> \$LOG_FILE
        
        # Convert HLA-HD output to arcasHLA-compatible JSON
        if [[ -f "\${PWD}/result/\${SAMPLE_ID}/\${SAMPLE_ID}_final.result.txt" ]]; then
            echo "Converting HLA-HD output to JSON format..." >> \$LOG_FILE
            
            # Convert HLA-HD TSV format to arcasHLA-compatible JSON
            convert_hlahd_to_json.py "\${PWD}/result/\${SAMPLE_ID}/\${SAMPLE_ID}_final.result.txt" "\${SAMPLE_ID}" > \${SAMPLE_ID}.genotype.json
            
            # Check if we got results from HLA-HD
            hd_allele_count=\$(grep -o '"[ABC][*][0-9]' "\${SAMPLE_ID}.genotype.json" | wc -l || echo "0")
            
            if [[ \$hd_allele_count -gt 0 ]]; then
                echo "SUCCESS: HLA-HD found \${hd_allele_count} HLA alleles for \${SAMPLE_ID}" >> \$LOG_FILE
                echo "Method used: HLA-HD (fallback)" >> \$LOG_FILE
            else
                echo "WARNING: HLA-HD completed but found no HLA alleles" >> \$LOG_FILE
                echo "Method used: HLA-HD (no results)" >> \$LOG_FILE
            fi
        else
            echo "ERROR: HLA-HD did not produce expected output file" >> \$LOG_FILE
            echo "Creating empty JSON output..." >> \$LOG_FILE
            echo '{"A": {}, "B": {}, "C": {}, "method": "HLA-HD_failed", "sample": "'"\${SAMPLE_ID}"'"}' > \${SAMPLE_ID}.genotype.json
        fi
    else
        echo "ERROR: HLA-HD command failed" >> \$LOG_FILE
        echo "Creating empty JSON output..." >> \$LOG_FILE
        echo '{"A": {}, "B": {}, "C": {}, "method": "both_failed", "sample": "'"\${SAMPLE_ID}"'"}' > \${SAMPLE_ID}.genotype.json
    fi
    
    # Cleanup temporary files to save space
    rm -f \${SAMPLE_ID}_R1.fastq \${SAMPLE_ID}_R2.fastq
    rm -rf result/
    
    echo "=== WORKFLOW COMPLETE ===" >> \$LOG_FILE
    echo "Final JSON preview:" >> \$LOG_FILE
    head -10 "\${SAMPLE_ID}.genotype.json" >> \$LOG_FILE || echo "Could not read final JSON" >> \$LOG_FILE
    """
    
    stub:
    """
    echo '{"A": {"A*01:01": 1.0}, "B": {"B*07:02": 1.0}, "C": {"C*07:02": 1.0}, "method": "stub", "sample": "${sampleName}"}' > ${sampleName}.genotype.json
    echo "Stub run finished!" > ${sampleName}_hla_typing.log
    """
}
