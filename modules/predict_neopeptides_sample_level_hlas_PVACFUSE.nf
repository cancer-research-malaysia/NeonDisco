process PREDICT_NEOPEPTIDES_SAMPLE_LEVEL_HLAS_PVACFUSE {
    errorStrategy 'retry'
    maxRetries 3
    
    label 'predictSampleNeopeptides'

    container "${params.container__pvactools}"
    
    publishDir "${params.outputDir}/${sampleName}/PVACFUSE-SAMPLE-LEVEL-HLAS-out", mode: 'copy',
        saveAs: { filename -> workflow.stubRun ? filename + ".stub" : filename }
    
    
    input:
        tuple val(sampleName), path(validatedAgfusionDir)
        path(cohortWideHLAList)
        path(metaDataDir) // Directory containing metadata files, including the reference proteome FASTA

    output:
        path("${sampleName}_*-HLA-pred/MHC_Class_I/${sampleName}.fasta")
        path("${sampleName}_*-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"), emit: predictedSampleSpecificNeopeptides
        path("${sampleName}_*-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.tsv")
        path("${sampleName}_*-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.tsv"), emit: aggregatedEpitopes
        path("${sampleName}_*-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.reference_matches.tsv"), emit: referenceMatches
        path("${sampleName}_sample_pvacfuse_execution_report.txt"), emit: executionReport
        path("${sampleName}-FI-validated-fusion-sample-hla-immunogenic-peptides-13aa.fasta"), emit: specializedFasta

    script:
    """
    # Initialize the execution report
    REPORT_FILE="${sampleName}_sample_pvacfuse_execution_report.txt"

    # Initialize parameter for fasta output amino acid flanking length
    FLANK_LENGTH=13
    
    # Function to log messages to both stdout and report file
    log_message() {
        local message="\$1"
        local timestamp=\$(date '+%Y-%m-%d %H:%M:%S')
        echo "[\$timestamp] \$message" | tee -a "\$REPORT_FILE"
    }
    
    # Start logging
    log_message "=== PVACFUSE SAMPLE-LEVEL HLA NEOPEPTIDE PREDICTION REPORT ==="
    log_message "Sample Name: ${sampleName}"
    log_message "Process Started: \$(date)"
    log_message ""
    
    log_message "INPUT VALIDATION:"
    log_message "Path to validated AGFusion directories: ${validatedAgfusionDir}"
    
    # Check if directory contains actual AGFusion results (ignore .empty files)
    agfusion_dirs=\$(find ${validatedAgfusionDir} -maxdepth 1 -type d ! -name "." ! -name ".." | wc -l)
    empty_marker=\$(find ${validatedAgfusionDir} -name ".empty" -type f | wc -l)
    
    log_message "Number of AGFusion directories found: \$agfusion_dirs"
    log_message "Number of .empty marker files: \$empty_marker"
    
    if [ \$agfusion_dirs -eq 0 ] || [ \$empty_marker -gt 0 ]; then
        log_message "WARNING: No validated AGFusion directories found for ${sampleName}"
        log_message "RESULT: Skipping neopeptide prediction process gracefully"
        log_message "Process Completed: \$(date)"
        log_message "STATUS: SKIPPED - No valid input data"
        exit 0
    fi
    
    log_message "INPUT VALIDATION: PASSED"
    log_message ""
    
    log_message "CONFIGURATION:"
    log_message "Path to cohort-wide HLA allotype TSV: ${cohortWideHLAList}"
    log_message "Prediction mode: Sample-level --> ${params.sampleHLANeoPred}"
    log_message "Number of cores to use: ${params.numCores * 2}"
    log_message "Reference proteome: ${metaDataDir}/Homo_sapiens.GRCh38.pep.all.fa.gz"
    log_message ""
    
    log_message "SAMPLE-SPECIFIC HLA EXTRACTION:"
    log_message "Extracting sample-specific HLA types from cohort-wide HLA file..."
    SSHLA=\$(awk -v sample="${sampleName}" 'NR > 1 && \$1 == sample {print \$2}' ${cohortWideHLAList})
    
    # Initialize variables for execution tracking
    PREDICTION_TYPE=""
    OUTPUT_DIR=""
    
    # Check if SSHLA is empty and handle accordingly
    if [ -z "\${SSHLA}" ]; then
        log_message "WARNING: No sample-specific HLA types found for ${sampleName}"
        log_message "FALLBACK: Assigning Southeast Asian-prevalent HLA types"
        SSHLA="HLA-A*11:01,HLA-A*24:02,HLA-A*02:07,HLA-A*33:03,HLA-B*46:02,HLA-B*44:03,HLA-B*40:01,HLA-C*03:04,HLA-C*01:02"
        PREDICTION_TYPE="SEA-SET"
        OUTPUT_DIR="${sampleName}_SEA-SET-HLA-pred"
        log_message "Using Southeast Asian HLA set: \${SSHLA}"
    else
        log_message "Sample-specific HLA types found: \${SSHLA}"
        PREDICTION_TYPE="SAMPLE-SPECIFIC"
        OUTPUT_DIR="${sampleName}_sample-level-HLA-pred"
    fi
    
    # Count HLA alleles
    HLA_COUNT=\$(echo "\${SSHLA}" | tr ',' '\n' | wc -l)
    log_message "Total number of HLA alleles: \${HLA_COUNT}"
    log_message "Prediction type: \${PREDICTION_TYPE}"
    log_message "Output directory: \${OUTPUT_DIR}"
    log_message ""
    
    log_message "PVACFUSE EXECUTION:"
    log_message "Starting pVacfuse for sample-specific HLA binding and immunogenicity prediction..."
    log_message "Algorithms: BigMHC_EL, BigMHC_IM, DeepImmuno, MHCflurry, MHCflurryEL, NetMHCpanEL, NetMHCcons, SMMPMBEC"
    log_message "Additional options:"
    log_message "  - IEDB install directory: /opt/iedb"
    log_message "  - Allele-specific binding thresholds: enabled"
    log_message "  - Reference proteome similarity: enabled"
    log_message "  - NetMHC stability prediction: enabled"
    log_message ""
    
    # Capture start time for duration calculation
    START_TIME=\$(date +%s)
    
    # Run pVacFuse with the determined HLA types
    log_message "Executing pVacfuse run command..."
    if pvacfuse run ${validatedAgfusionDir} ${sampleName} \${SSHLA} \\
    BigMHC_EL BigMHC_IM DeepImmuno MHCflurry MHCflurryEL NetMHCpanEL NetMHCcons SMMPMBEC \\
    "\${OUTPUT_DIR}" \\
    --iedb-install-directory /opt/iedb \\
    --allele-specific-binding-thresholds \\
    --run-reference-proteome-similarity \\
    --peptide-fasta ${metaDataDir}/Homo_sapiens.GRCh38.pep.all.fa.gz \\
    --netmhc-stab \\
    -t ${params.numCores * 2} \\
    -a sample_name 2>&1 | tee -a "\$REPORT_FILE"; then
        
        # Calculate execution time
        END_TIME=\$(date +%s)
        DURATION=\$((END_TIME - START_TIME))
        DURATION_MIN=\$((DURATION / 60))
        DURATION_SEC=\$((DURATION % 60))
        
        log_message ""
        log_message "EXECUTION COMPLETED SUCCESSFULLY"
        log_message "Execution time: \${DURATION_MIN}m \${DURATION_SEC}s"
        log_message "Prediction type used: \${PREDICTION_TYPE}"
        
        # Check if output file was created and get basic stats
        OUTPUT_FILE="="\${OUTPUT_DIR}/MHC_Class_I/${sampleName}.filtered.tsv"
        if [ -f "\$OUTPUT_FILE" ]; then
            RESULT_COUNT=\$(tail -n +2 "\$OUTPUT_FILE" | wc -l)
            log_message "Output file created: \$OUTPUT_FILE"
            log_message "Number of predicted neopeptides: \$RESULT_COUNT"

            # copy and rename reference matches file if it exists
            REF_MATCHES_FILE="\${OUTPUT_DIR}/MHC_Class_I/${sampleName}.all_epitopes.aggregated.tsv.reference_matches"
            if [ -f "\$REF_MATCHES_FILE" ]; then
                cp "\$REF_MATCHES_FILE" "\${OUTPUT_DIR}/MHC_Class_I/${sampleName}.all_epitopes.aggregated.reference_matches.tsv"
                log_message "Reference matches file renamed to: ${sampleName}.all_epitopes.aggregated.reference_matches.tsv"
            else
                log_message "No reference matches file generated."
            fi
            
            AGGREGATED_FILE="\${OUTPUT_DIR}/MHC_Class_I/${sampleName}.all_epitopes.aggregated.tsv"
            if [ -f "\$AGGREGATED_FILE" ]; then
                AGG_COUNT=\$(tail -n +2 "\$AGGREGATED_FILE" | wc -l)
                log_message "Aggregated epitopes file found: \$AGGREGATED_FILE"

                # now run pvacfuse generate_protein_fasta to create a specialized FASTA file
                log_message ""
                log_message "Generating specialized FASTA file with pVacfuse generate_protein_fasta to get 13 aa upstream and downstream of fusion junctions..."
                if pvacfuse generate_protein_fasta --input-tsv "\$AGGREGATED_FILE" --aggregate-report-evaluation Pending ${validatedAgfusionDir} "\$FLANK_LENGTH" ${sampleName}-FI-validated-fusion-sample-hla-immunogenic-peptides-13aa.fasta 2>&1 | tee -a "\$REPORT_FILE"; then
                    log_message "${sampleName}-FI-validated-fusion-sample-hla-immunogenic-peptides-13aa.fasta created"
                else
                    log_message "WARNING: pVacfuse generate_protein_fasta execution failed"
                    exit 1
            fi

            # Additional statistics if results exist
            if [ \$RESULT_COUNT -gt 0 ]; then
                # Get unique fusion count if available
                UNIQUE_FUSIONS=\$(tail -n +2 "\$OUTPUT_FILE" | cut -f1 | sort -u | wc -l)
                log_message "Number of unique fusions with neopeptides: \$UNIQUE_FUSIONS"
            fi
        else
            log_message "WARNING: Expected output file not found: \$OUTPUT_FILE"
            exit 1
        fi

        log_message ""
        log_message "HLA TYPE SUMMARY:"
        log_message "HLA types used for prediction: \${SSHLA}"
        log_message "Source: \${PREDICTION_TYPE}"
        if [ "\${PREDICTION_TYPE}" = "SEA-SET" ]; then
            log_message "Note: Sample-specific HLAs were not available, so Southeast Asian prevalent HLAs were used as fallback"
        fi
        
        log_message ""
        log_message "Process Completed: \$(date)"
        log_message "STATUS: SUCCESS"
        
    else
        log_message ""
        log_message "ERROR: pVacFuse execution failed"
        log_message "HLA types attempted: \${SSHLA}"
        log_message "Prediction type: \${PREDICTION_TYPE}"
        log_message "Output directory: \${OUTPUT_DIR}"
        log_message "Process Completed: \$(date)"
        log_message "STATUS: FAILED"
        log_message ""
        log_message "Please check the error messages above for troubleshooting information."
        exit 1
    fi
    
    log_message ""
    log_message "=== END OF REPORT ==="
    """

    stub:
    """
    # Create stub report
    REPORT_FILE="${sampleName}_sample_pvacfuse_execution_report.txt"
    cat > "\$REPORT_FILE" << EOF
=== PVACFUSE SAMPLE-LEVEL HLA NEOPEPTIDE PREDICTION REPORT (STUB RUN) ===
Sample Name: ${sampleName}
Process Started: \$(date)

This is a stub run - no actual processing was performed.
All outputs are placeholder files for workflow testing purposes.

STUB CONFIGURATION:
- Input: ${validatedAgfusionDir}
- HLA List: ${cohortWideHLAList} 
- Metadata: ${metaDataDir}
- Prediction Mode: Sample-level
- Cores: ${params.numCores * 2}

Process Completed: \$(date)
STATUS: STUB RUN COMPLETE
=== END OF REPORT ===
EOF
    
    # Create stub output file - could be either sample-level or SEA-SET depending on HLA availability
    mkdir -p "${sampleName}_sample-level-HLA-pred/MHC_Class_I"
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.filtered.tsv"
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.tsv"
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.tsv"
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.all_epitopes.aggregated.reference_matches.tsv
    touch "${sampleName}_sample-level-HLA-pred/MHC_Class_I/${sampleName}.fasta"
    touch "${sampleName}-FI-validated-fusion-sample-hla-immunogenic-peptides-13aa.fasta"
    echo "Stub run finished!" | tee -a "\$REPORT_FILE"
    """
}
