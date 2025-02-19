#!/bin/bash

# Default values
CORES=2
MODE=""
INPUT="/mnt"  # Default to Docker mount point
APP_PATH="/home/app"  # Path to application in Docker container

# Helper function to print usage
print_usage() {
    cat << EOF
Usage: $(basename $0) [options]

Options:
    -m, --mode MODE       	Operation mode (Required)
                         		Valid modes: BAM2BED, BED2JUNC, CALLJUNC, DE, GO, DAS
    -i, --input PATH     	Input file or directory path
    -g, --group FILE     	Group file for DE/DAS analysis
    -c, --cores N        	Number of CPU cores to use (Default: 2)
    -h, --help          	Show this help message

Modes:
    BAM2BED      	    - Convert single BAM file to BED format
    BED2JUNC 		    - Run BED files to junction analysis
    CALLJUNC            - Process multiple BAM files and identify junctions
    DE             	    - Differential expression analysis
    GO           	    - Gene Ontology enrichment analysis
    DAS            	    - Differential alternative splicing analysis
EOF
}

# Function to validate mode
validate_mode() {
    local mode=$1
    local valid_modes=("BAM2BED" "BED2JUNC" "CALLJUNC" "DE" "GO" "DAS")
    for valid_mode in "${valid_modes[@]}"; do
        if [[ "$mode" == "$valid_mode" ]]; then
            return 0
        fi
    done
    echo "Error: Invalid mode '$mode'" >&2
    return 1
}

# Parse command line arguments
while getopts ":m:i:g:c:h-:" opt; do
    case $opt in
        -)
            case "${OPTARG}" in
                mode)        	MODE="${!OPTIND}"; OPTIND=$((OPTIND + 1));;
                input)       	INPUT="${!OPTIND}"; OPTIND=$((OPTIND + 1));;
                group)       	GROUP_FILE="${!OPTIND}"; OPTIND=$((OPTIND + 1));;
                cores)       	CORES="${!OPTIND}"; OPTIND=$((OPTIND + 1));;
                help)        	print_usage; exit 0;;
                *)          	echo "Invalid option: --${OPTARG}" >&2; exit 1;;
            esac;;
        m)  MODE=$OPTARG;;
        i)  INPUT=$OPTARG;;
        g)  GROUP_FILE=$OPTARG;;
        c)  CORES=$OPTARG;;
        h)  print_usage; exit 0;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
        :)  echo "Option -$OPTARG requires an argument" >&2; exit 1;;
    esac
done

# Validate required arguments
if [[ -z "$MODE" ]]; then
    echo "Error: Mode (-m, --mode) is required" >&2
    print_usage
    exit 1
fi

# Validate mode
if ! validate_mode "$MODE"; then
    print_usage
    exit 1
fi

# Change to working directory if INPUT is dir
if [[ -d "${INPUT}" ]]; then
	cd "$INPUT" 
fi

echo "Current folder is $PWD..."

# Function for BAM to BED conversion
run_BAMtoBED() {
    local bam_file=$1
    echo "Processing BAM file: $bam_file"
    
    echo "Start to get junction bed"
    python "$APP_PATH/altanalyze/import_scripts/BAMtoJunctionBED.py" \
        --i "$bam_file" \
        --species Hs \
        --r "$APP_PATH/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt"

    echo "Start to get exon bed"
    python "$APP_PATH/altanalyze/import_scripts/BAMtoExonBED.py" \
        --i "$bam_file" \
        --r "$APP_PATH/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs.bed" \
        --s Hs
}

# Function to setup MultiPath-PSI analysis
setup_multipath_PSI() {
    local bed_dir=$1
    local input=$2
	local task="original"

    mkdir -p "${input}/altanalyze_output/ExpressionInput"
    
    # Create groups file
    local groups_file="${input}/altanalyze_output/ExpressionInput/groups.${task}.txt"
    touch "$groups_file"
    
    cd "$bed_dir" || exit 1
    local count=0
    for file in *__junction.bed; do
        [[ -f "$file" ]] || continue
        local stream
        stream=$(echo "$file" | sed 's/__junction.bed/.bed/g')
        if ((count % 2 == 0)); then
            stream+=$'\t1\texp'
        else
            stream+=$'\t2\tctl'
        fi
        echo -e "$stream" >> "../$groups_file"
        ((count++))
    done
    cd - || exit 1
    
    # Create comps file
    echo -e '1\t2' > "${input}/altanalyze_output/ExpressionInput/comps.${task}.txt"
    
    # Run AltAnalyze
    python "$APP_PATH/altanalyze/AltAnalyze.py" \
        --species Hs \
        --platform RNASeq \
        --version EnsMart91 \
        --bedDir "$bed_dir" \
        --output "$input/altanalyze_output" \
        --groupdir "$input/altanalyze_output/ExpressionInput/groups.${task}.txt" \
        --compdir "$input/altanalyze_output/ExpressionInput/comps.${task}.txt" \
        --expname "$task" \
        --runGOElite no
        
    echo "Now pruning the raw junction count matrix..."
    python "$APP_PATH/prune.py"
}

# Main execution based on mode
case "$MODE" in
    BAM2BED)
        [[ -f "$INPUT" ]] || { echo "[ERROR] BAM file not found or not a file: $INPUT" >&2; exit 1; }
        echo "Running bam to bed workflow; bam file is $INPUT"
        run_BAMtoBED "$INPUT"
        ;;
        
    BED2JUNC)
        [[ -d "${INPUT}/bed" ]] || { echo "[ERROR] BED directory not found or input is a file: $INPUT/bed" >&2; exit 1; }
        echo "Running bed to junction workflow; BED folder is $INPUT/bed"
        setup_multipath_PSI "${INPUT}/bed" "${INPUT}"
        ;;
        
    CALLJUNC)
        [[ -d "${INPUT}/bam" ]] || { echo "[ERROR] BAM directory not found or input is a file: $INPUT/bam" >&2; exit 1; }
        echo "Identify splicing junction; bam folder is $INPUT/bam, using $CORES cores"
        
        # Collect BAM files for parallelization
        find "${INPUT}/bam" -name "*.bam" > "${INPUT}/samples.txt"
        
        # Export necessary variables for parallel execution
        export -f run_BAMtoBED
        export TMPDIR=/tmp
        export g_bam_folder="${INPUT}/bam"
        
        # Process BAM files in parallel
        parallel -P "$CORES" run_BAMtoBED {} < "${INPUT}/samples.txt"
        
        # Move BED files to bed folder
        mkdir -p "${INPUT}/bed"
        find "$INPUT" -name "*.bed" -exec mv {} "${INPUT}/bed/" \;
        
        # Run MultiPath-PSI analysis
        setup_multipath_PSI "${INPUT}/bed" "${INPUT}"
        ;;
        
    DE)
        [[ -f "$GROUP_FILE" ]] || { echo "[ERROR] Group file not found: $GROUP_FILE" >&2; exit 1; }
        [[ -d "$INPUT" ]] || { echo "[ERROR] Output directory not found or is a file: $INPUT" >&2; exit 1; }
        
        echo "Identify differentially expressed genes; AltAnalyze output folder is $INPUT, group file is $GROUP_FILE"
        
		mkdir -p "${INPUT}/ExpressionInput"
		
		python "$APP_PATH/altanalyze/stats_scripts/metaDataAnalysis.py" \
            --p RNASeq \
            --s Hs \
            --adjp yes \
            --pval 1 \
            --f 1 \
            --i "${INPUT}/ExpressionInput/exp.original-steady-state.txt" \
            --m "$GROUP_FILE"
        ;;
        
    GO)
        [[ -f "$INPUT" ]] || { echo "[ERROR] Gene list file not found or is a directory: $INPUT" >&2; exit 1; }
        echo "Gene enrichment analysis using GO-Elite, gene list file is $INPUT"
        
        # BioMarkers analysis
        mkdir -p "$INPUT/GO_Elite_result_BioMarkers"
        
		python "$APP_PATH/altanalyze/GO_Elite.py" \
            --species Hs \
            --mod Ensembl \
            --pval 0.05 \
            --num 3 \
            --input "$INPUT" \
            --output "$INPUT/GO_Elite_result_BioMarkers" \
            --dataToAnalyze BioMarkers
            
        # Gene Ontology analysis
        mkdir -p "$INPUT/GO_Elite_result_GeneOntology"

        python "$APP_PATH/altanalyze/GO_Elite.py" \
            --species Hs \
            --mod Ensembl \
            --pval 0.05 \
            --num 3 \
            --input "$INPUT" \
            --output "$INPUT/GO_Elite_result_GeneOntology" \
            --dataToAnalyze GeneOntology
        ;;
        
    DAS)
        [[ -f "$GROUP_FILE" ]] || { echo "[ERROR] Group file not found: $GROUP_FILE" >&2; exit 1; }
        [[ -d "$INPUT" ]] || { echo "Error: Output directory not found or is not a directory: $INPUT" >&2; exit 1; }
        
        echo "Identify differentially spliced event; AltAnalyze output folder is $INPUT, group file is $GROUP_FILE"

		mkdir -p "$INPUT/AltResults/AlternativeOutput"

        python "$APP_PATH/altanalyze/stats_scripts/metaDataAnalysis.py" \
            --p PSI \
            --dPSI 0 \
            --pval 1 \
            --adjp no \
            --i "$INPUT/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt" \
            --m "$GROUP_FILE"
        ;;
esac
