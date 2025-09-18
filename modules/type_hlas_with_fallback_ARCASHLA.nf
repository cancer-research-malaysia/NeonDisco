// Single process approach - arcasHLA with HLA-HD fallback in one instance
process TYPE_HLAS_WITH_FALLBACK_ARCASHLA {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'typeHLAs'
    
    container "${params.container__hlatyping}"  // Container with both tools
    
    publishDir "${params.outputDir}/${sampleName}/HLA-TYPING-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    input:
    tuple val(sampleName), path(bam), path(bamIdx), path(querynameBam)

    output:
    tuple val(sampleName), path("${sampleName}.genotype.json"), emit: hla_json
    path "${sampleName}_hla_typing.log"
    
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
                echo "SUCCESS: arcasHLA found \${allele_count} HLA alleles for \${SAMPLE_ID}" | tee -a \$LOG_FILE
                echo "Method used: arcasHLA" | tee -a \$LOG_FILE
                echo "=== WORKFLOW COMPLETE ===" | tee -a \$LOG_FILE
                exit 0  # Success - exit early
            else
                echo "WARNING: arcasHLA completed but found no HLA alleles" | tee -a \$LOG_FILE
                echo "JSON content preview:" | tee -a \$LOG_FILE
                head -5 "\${SAMPLE_ID}.genotype.json" | tee -a \$LOG_FILE || echo "Could not read JSON" | tee -a \$LOG_FILE
            fi
        else
            echo "WARNING: arcasHLA did not produce expected JSON output" >> \$LOG_FILE
            exit 1
        fi
    else
        echo "WARNING: arcasHLA command failed" >> \$LOG_FILE
        exit 1
    fi
    
    echo "" | tee -a \$LOG_FILE
    echo "=== FALLING BACK TO HLA-HD ===" | tee -a \$LOG_FILE
    echo "arcasHLA was unsuccessful, trying HLA-HD..." | tee -a \$LOG_FILE
    
    # Convert BAM to FASTQ for HLA-HD
    echo "Converting BAM to FASTQ..." | tee -a \$LOG_FILE
    if fish-sampicard--nf.sh \${SAMPLE_ID} \${BAM} 2>&1 | tee -a \$LOG_FILE; then
        echo "BAM to FASTQ conversion completed successfully" | tee -a \$LOG_FILE
        # Run HLA-HD with the generated FASTQ files
        echo "Starting HLA-HD typing..." | tee -a \$LOG_FILE
        if [[ -f "\${SAMPLE_ID}_Bam2Fq_R1.fastq" && -f "\${SAMPLE_ID}_Bam2Fq_R2.fastq" ]]; then
            echo "FASTQ files generated: \${SAMPLE_ID}_Bam2Fq_R1.fastq and \${SAMPLE_ID}_Bam2Fq_R2.fastq" >> \$LOG_FILE
            if hlahd--nf.sh \${SAMPLE_ID} \${SAMPLE_ID}_Bam2Fq_R1.fastq \${SAMPLE_ID}_Bam2Fq_R2.fastq ${params.numCores} 2>&1 | tee -a \$LOG_FILE; then
                echo "HLA-HD command completed successfully" | tee -a \$LOG_FILE
            else
                echo "ERROR: HLA-HD command failed" | tee -a \$LOG_FILE
                exit 1
            fi
        else
            echo "ERROR: FASTQ files not found after conversion" | tee -a \$LOG_FILE
            exit 1
        fi
    else
        echo "ERROR: BAM to FASTQ conversion failed" | tee -a \$LOG_FILE
        exit 1
    fi
    
        
    # Convert HLA-HD output to arcasHLA-compatible JSON
    if [[ -f "HLAHD/\${SAMPLE_ID}/result/\${SAMPLE_ID}_final.result.txt" ]]; then
            echo "Converting HLA-HD output to JSON format..." | tee -a \$LOG_FILE
            
            # Convert HLA-HD TSV format to arcasHLA-compatible JSON
            convert-hlahd-to-json--nf.py "HLAHD/\${SAMPLE_ID}/result/\${SAMPLE_ID}_final.result.txt" "\${SAMPLE_ID}" > \${SAMPLE_ID}.genotype.json
            
            # Check if we got results from HLA-HD
            hd_allele_count=\$(grep -o '"[ABC][*][0-9]' "\${SAMPLE_ID}.genotype.json" | wc -l || echo "0")
            
            if [[ \$hd_allele_count -gt 0 ]]; then
                echo "SUCCESS: HLA-HD found \${hd_allele_count} HLA alleles for \${SAMPLE_ID}" | tee -a \$LOG_FILE
                echo "Method used: HLA-HD (fallback)" | tee -a \$LOG_FILE
            else
                echo "WARNING: HLA-HD completed but found no HLA alleles" | tee -a \$LOG_FILE
                echo "Method used: HLA-HD (no results)" | tee -a \$LOG_FILE
            fi
    else
        echo "ERROR: HLA-HD did not produce expected output file" | tee -a \$LOG_FILE
        echo "Creating empty JSON output..." | tee -a \$LOG_FILE
        echo '{}' > \${SAMPLE_ID}.genotype.json
    fi

    echo "=== WORKFLOW COMPLETE ===" | tee -a \$LOG_FILE
    echo "Final JSON preview:" | tee -a \$LOG_FILE
    head "\${SAMPLE_ID}.genotype.json" >> \$LOG_FILE || echo "Could not read final JSON" | tee -a \$LOG_FILE
    """
    
    stub:
    """
    echo '{"A": {"A*01:01": 1.0}, "B": {"B*07:02": 1.0}, "C": {"C*07:02": 1.0}, "method": "stub", "sample": "${sampleName}"}' > ${sampleName}.genotype.json
    echo "Stub run finished!" > ${sampleName}_hla_typing.log
    """
}
