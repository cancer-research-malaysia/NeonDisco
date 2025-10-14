#!/usr/bin/env bash

# Script to run pVacfuse for sample-level HLA neopeptide prediction using validated AGFusion results
# Handles cases with missing sample-specific HLA types by falling back to a predefined Southeast Asian HLA set
# Logs detailed execution report including input validation, configuration, execution status, and summary statistics
# Usage: bash predict-neopeptides-sample-hla-pvacfuse--nf.sh <sample_name> <validated_agfusion_dir> <cohort_wide_hla_list> <meta_data_dir> <prediction_mode> <flank_length> <num_cores>

set -euo pipefail

# Function to log messages to both stdout and report file
log_message() {
	local message="$1"
	local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
	echo "[$timestamp] $message" | tee -a "$REPORT_FILE"
}

# Input arguments
SAMPLE_NAME="$1"  # Sample name identifier
VALIDATED_AGF_DIR="$2"  # Directory containing validated AGFusion results for the sample
COHORT_HLA_LIST="$3"  # Path to cohort-wide HLA allotype TSV file
METADATA_DIR="$4"  # Directory containing metadata files (e.g., reference proteome)
PREDICTION_MODE="$5"  # Prediction mode (should be params.sampleHLANeoPred for this script)
FLANK_LENGTH="$6"  # Length of amino acid flanking region for FASTA output
NUM_CORES="$7"  # Number of CPU cores to use for parallel processing

# Initialize the execution report
REPORT_FILE="${SAMPLE_NAME}_sample_HLA_pvacfuse_execution_report.txt"

# Start logging
log_message "=== PVACFUSE SAMPLE-LEVEL HLA NEOPEPTIDE PREDICTION REPORT ==="
log_message "Sample Name: ${SAMPLE_NAME}"
log_message "Process Started: $(date)"
log_message ""

log_message "INPUT VALIDATION:"
log_message "Path to validated AGFusion directories: ${VALIDATED_AGF_DIR}"

# Check if directory contains actual AGFusion results (ignore .empty files)
agfusion_dirs=$(find ${VALIDATED_AGF_DIR} -maxdepth 1 -type d ! -name "." ! -name ".." | wc -l)
empty_marker=$(find ${VALIDATED_AGF_DIR} -name ".empty" -type f | wc -l)

log_message "Number of AGFusion directories found: $agfusion_dirs"
log_message "Number of .empty marker files: $empty_marker"

if [ $agfusion_dirs -eq 0 ] || [ $empty_marker -gt 0 ]; then
	log_message "WARNING: No validated AGFusion directories found for ${SAMPLE_NAME}"
	log_message "RESULT: Skipping neopeptide prediction process gracefully"
	log_message "Process Completed: $(date)"
	log_message "STATUS: SKIPPED - No valid input data"
	exit 0
fi

log_message "INPUT VALIDATION: PASSED"
log_message ""

log_message "CONFIGURATION:"
log_message "Path to cohort-wide HLA allotype TSV: ${COHORT_HLA_LIST}"
log_message "Prediction mode: Sample-level --> ${PREDICTION_MODE}"
log_message "Number of cores to use: ${NUM_CORES}"
log_message "Reference proteome: ${METADATA_DIR}/Homo_sapiens.GRCh38.pep.all.fa.gz"
log_message ""

log_message "SAMPLE-SPECIFIC HLA EXTRACTION:"
log_message "Extracting sample-specific HLA types from cohort-wide HLA file..."
SSHLA=$(awk -v sample="${SAMPLE_NAME}" 'NR > 1 && $1 == sample {print $2}' ${COHORT_HLA_LIST})

# Initialize variables for execution tracking
PREDICTION_TYPE=""
OUTPUT_DIR=""

# Check if SSHLA is empty and handle accordingly
if [ -z "${SSHLA}" ]; then
	log_message "WARNING: No sample-specific HLA types found for ${SAMPLE_NAME}"
	log_message "FALLBACK: Assigning Southeast Asian-prevalent HLA types"
	SSHLA="HLA-A*11:01,HLA-A*24:02,HLA-A*02:07,HLA-A*33:03,HLA-B*46:02,HLA-B*44:03,HLA-B*40:01,HLA-C*03:04,HLA-C*01:02"
	PREDICTION_TYPE="SEA-SET"
	OUTPUT_DIR="${SAMPLE_NAME}_SEA-SET-HLA-pred"
	log_message "Using Southeast Asian HLA set: ${SSHLA}"
else
	log_message "Sample-specific HLA types found: ${SSHLA}"
	PREDICTION_TYPE="SAMPLE-SPECIFIC"
	OUTPUT_DIR="${SAMPLE_NAME}_sample-level-HLA-pred"
fi

# Count HLA alleles
HLA_COUNT=$(echo "${SSHLA}" | tr ',' '\n' | wc -l)
log_message "Total number of HLA alleles: ${HLA_COUNT}"
log_message "Prediction type: ${PREDICTION_TYPE}"
log_message "Output directory: ${OUTPUT_DIR}"
log_message ""

log_message "PVACFUSE EXECUTION:"
log_message "Starting pVacfuse for sample-specific HLA binding and immunogenicity prediction..."
log_message "Algorithms: BigMHC_IM, DeepImmuno, MHCflurry, MHCflurryEL, NetMHCpanEL, NetMHCcons, SMMPMBEC"
log_message "Additional options:"
log_message "  - IEDB install directory: /opt/iedb"
log_message "  - Allele-specific binding thresholds: enabled"
log_message "  - Reference proteome similarity: enabled"
log_message "  - NetMHC stability prediction: enabled"
log_message ""

# Capture start time for duration calculation
START_TIME=$(date +%s)

# Run pVacFuse with the determined HLA types
log_message "Executing pVacfuse run command..."
if pvacfuse run ${VALIDATED_AGF_DIR} ${SAMPLE_NAME} ${SSHLA} BigMHC_IM DeepImmuno MHCflurry MHCflurryEL NetMHCpanEL NetMHCcons SMMPMBEC "${OUTPUT_DIR}" --iedb-install-directory /opt/iedb --allele-specific-binding-thresholds --downstream-sequence-length full --run-reference-proteome-similarity --peptide-fasta ${METADATA_DIR}/Homo_sapiens.GRCh38.pep.all.fa.gz --netmhc-stab -t ${NUM_CORES} -a sample_name 2>&1 | tee -a "$REPORT_FILE"; then
	# Calculate execution time
	END_TIME=$(date +%s)
	DURATION=$((END_TIME - START_TIME))
	DURATION_MIN=$((DURATION / 60))
	DURATION_SEC=$((DURATION % 60))
	
	log_message ""
	log_message "EXECUTION COMPLETED SUCCESSFULLY"
	log_message "Execution time: ${DURATION_MIN}m ${DURATION_SEC}s"
	log_message "Prediction type used: ${PREDICTION_TYPE}"
	
	# Check if output file was created and get basic stats
	OUTPUT_FILE="${OUTPUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.filtered.tsv"

	if [ -f "$OUTPUT_FILE" ]; then
		RESULT_COUNT=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
		log_message "Output file created: $OUTPUT_FILE"
		log_message "Number of predicted neopeptides: $RESULT_COUNT"

		# copy and rename reference matches file if it exists
		REF_MATCHES_FILE="${OUTPUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.all_epitopes.aggregated.tsv.reference_matches"
		if [ -f "$REF_MATCHES_FILE" ]; then
			cp "$REF_MATCHES_FILE" "${OUTPUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.all_epitopes.aggregated.reference_matches.tsv"
			log_message "Reference matches file renamed to: ${SAMPLE_NAME}.all_epitopes.aggregated.reference_matches.tsv"
		else
			log_message "No reference matches file generated."
		fi
		
		AGGREGATED_FILE="${OUTPUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.all_epitopes.aggregated.tsv"
		if [ -f "$AGGREGATED_FILE" ]; then
			AGG_COUNT=$(tail -n +2 "$AGGREGATED_FILE" | wc -l)
			log_message "Aggregated epitopes file found: $AGGREGATED_FILE"

			# now run pvacfuse generate_protein_fasta to create a specialized FASTA file
			log_message ""
			log_message "Generating specialized FASTA file with pVacfuse generate_protein_fasta to get 13 aa upstream and downstream of fusion junctions..."
			if pvacfuse generate_protein_fasta --input-tsv "$AGGREGATED_FILE" --aggregate-report-evaluation Pending ${VALIDATED_AGF_DIR} ${FLANK_LENGTH} ${SAMPLE_NAME}-FI-validated-fusion-sample-HLA-immunogenic-peptides-13aa.fasta 2>&1 | tee -a "$REPORT_FILE"; then
				log_message "${SAMPLE_NAME}-FI-validated-fusion-sample-HLA-immunogenic-peptides-13aa.fasta created"
			else
				log_message "WARNING: pVacfuse generate_protein_fasta execution failed"
				exit 1
			fi
		fi

		# Additional statistics if results exist
		if [ $RESULT_COUNT -gt 0 ]; then
			# Get unique fusion count if available
			UNIQUE_FUSIONS=$(tail -n +2 "$OUTPUT_FILE" | cut -f1 | sort -u | wc -l)
			log_message "Number of unique fusions with neopeptides: $UNIQUE_FUSIONS"
		fi
	else
		log_message "WARNING: Expected output file not found: $OUTPUT_FILE"
		touch "${OUTPUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.filtered.tsv"
		log_message "Created empty output file to maintain workflow consistency."
		RESULT_COUNT=0
		UNIQUE_FUSIONS=0
		log_message "Number of predicted neopeptides: $RESULT_COUNT"
		log_message "Creating other expected output files as empty to maintain workflow consistency."
		touch "${OUTPUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.all_epitopes.tsv"
		touch "${OUTPUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.all_epitopes.aggregated.tsv"
		touch "${OUTPUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.all_epitopes.aggregated.reference_matches.tsv"
		touch "${OUTPUT_DIR}/MHC_Class_I/${SAMPLE_NAME}.fasta"
		touch "${SAMPLE_NAME}-FI-validated-fusion-sample-HLA-immunogenic-peptides-13aa.fasta"
		
		log_message "Created empty aggregated epitopes and reference matches files."
	fi

	log_message ""
	log_message "HLA TYPE SUMMARY:"
	log_message "HLA types used for prediction: ${SSHLA}"
	log_message "Source: ${PREDICTION_TYPE}"
	if [ "${PREDICTION_TYPE}" = "SEA-SET" ]; then
		log_message "Note: Sample-specific HLAs were not available, so Southeast Asian prevalent HLAs were used as fallback"
	fi
	
	log_message ""
	log_message "Process Completed: $(date)"
	log_message "STATUS: SUCCESS"
	
else
	log_message ""
	log_message "ERROR: pVacFuse execution failed"
	log_message "HLA types attempted: ${SSHLA}"
	log_message "Prediction type: ${PREDICTION_TYPE}"
	log_message "Output directory: ${OUTPUT_DIR}"
	log_message "Process Completed: $(date)"
	log_message "STATUS: FAILED"
	log_message ""
	log_message "Please check the error messages above for troubleshooting information."
	exit 1
fi

log_message ""
log_message "=== END OF REPORT ==="
