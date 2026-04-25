#!/bin/bash
set -euo pipefail

# Base directories
readonly BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Logging configuration
readonly LOG_DIR="$BASE_DIR/logs"
readonly LOG_FILE="$LOG_DIR/blast_run_${TIMESTAMP}.log"
readonly CLEAR_LOGS=true  # Set to true to clear logs folder before each run

# Database and BLAST configuration
readonly RECREATE_DATABASE=false  # Set to true to recreate BLAST databases even if they exist
readonly OVERWRITE_BLAST_RESULTS=true  # Set to true to overwrite existing BLAST results, false to skip existing runs

# Input FASTA files
readonly H_FINAL_PHYLO_GRF="$BASE_DIR/1_INPUTS/h_Final_Phylo_Datasets/Final_Phylo_GRF_ProtSeq_with_SmelGRF_e-value_1e-5.fasta"
readonly H_FINAL_PHYLO_GIF="$BASE_DIR/1_INPUTS/h_Final_Phylo_Datasets/Final_Phylo_GIF_ProtSeq_with_SmelGIF_e-value_1e-5.fasta"

# Database and Output directories
readonly DB_DIR="$BASE_DIR/0_BLAST_DB"
readonly OUT_FINAL_GRF="$BASE_DIR/2_OUTPUTS/Final_GRF"
readonly OUT_FINAL_GIF="$BASE_DIR/2_OUTPUTS/Final_GIF"
readonly EGGPLANT_v4_1_PROTEINS="$BASE_DIR/1_INPUTS/Eggplant_V4_1_protein_function.fa"

# Available genomes (uncomment as needed)
readonly GENOME_LIST=(
    #"$EGGPLANT_V4_1_TRANSCRIPTS"
    "$EGGPLANT_v4_1_PROTEINS"
    # Add other genomes here when needed
)

# BLAST parameters - Recommended for gene family identification
readonly E_VALUE_LIST=( 
    #1e-10 
    1e-15
    #1e-20   # stringent, ensures only strong homology is retained.
)
readonly WORD_SIZE_LIST=(
    #3  # Default word size for BLASTP (valid: 2-7)
    2  # More sensitive but slower
)
readonly NUM_THREADS=4  # Number of CPU threads to use (adjust based on your system)

# ========================================
# Functions
# ========================================

# Setup logging system
setup_logging() {
    mkdir -p "$LOG_DIR"
    
    # Clear logs if toggle is enabled
    if [[ "$CLEAR_LOGS" == true ]]; then
        echo "Clearing logs folder..."
        rm -f "$LOG_DIR"/*.log 2>/dev/null
        echo "Logs cleared."
    fi
    
    # Start logging
    exec > >(tee -a "$LOG_FILE") 2>&1
    echo "Log file: $LOG_FILE"
    
    # Suppress numpy warnings in Python scripts
    export PYTHONWARNINGS="ignore::UserWarning"
}

# Discover and categorize FASTA files (handles single file or directory)
discover_fasta_files() {
    local query_path=$1
    local prefix=$2
    local -n all_ref="${prefix}_GRF_GIF_FASTA_LIST"
    local -n grf_ref="${prefix}_GRF_FASTA_LIST"
    local -n gif_ref="${prefix}_GIF_FASTA_LIST"
    
    # Clear arrays
    all_ref=() grf_ref=() gif_ref=()
    
    [[ ! -e "$query_path" ]] && { echo "Path not found: $query_path"; return 1; }
    
    # If it's a file, process it directly
    if [[ -f "$query_path" ]]; then
        filepath=$(realpath "$query_path")
        filename=$(basename "$query_path")
        case "$filename" in
            *GRF*) all_ref+=("$filepath"); grf_ref+=("$filepath") ;;
            *GIF*) all_ref+=("$filepath"); gif_ref+=("$filepath") ;;
        esac
        printf "  Found: %d total | %d GRF | %d GIF\n" "${#all_ref[@]}" "${#grf_ref[@]}" "${#gif_ref[@]}"
        return 0
    fi
    
    # Otherwise treat as directory
    [[ ! -d "$query_path" ]] && { echo "Not a file or directory: $query_path"; return 1; }
    
    echo "Scanning: $query_path"
    while IFS= read -r -d '' file; do
        filepath=$(realpath "$file")
        filename=$(basename "$file")
        
        case "$filename" in
            *GRF*) 
                all_ref+=("$filepath")
                grf_ref+=("$filepath")
                ;;
            *GIF*)
                all_ref+=("$filepath")
                gif_ref+=("$filepath")
                ;;
        esac
    done < <(find "$query_path" -type f \( -iname "*.fa" -o -iname "*.fasta" \) -print0)
    
    printf "  Found: %d total | %d GRF | %d GIF\n" "${#all_ref[@]}" "${#grf_ref[@]}" "${#gif_ref[@]}"
}

# Execute BLASTP for array of query files against genome databases
run_blast() {
    local -n genomes=$1    # Array of genome databases
    local -n queries=$2    # Array of query FASTA files
    local gene_type=$3     # Gene type label (e.g., "GRF", "GIF", "GRF-GIF")
    local output_dir=$4    # Base output directory

    echo "Running BLASTP for $gene_type (${#queries[@]} queries vs ${#genomes[@]} genomes)"
    
    for genome in "${genomes[@]}"; do
        # Validate genome file exists
        if [[ ! -f "$genome" ]]; then
            echo "  ERROR: Genome file not found: $genome"
            continue
        fi
        
        # Extract basename and remove extension
        local genome_base=$(basename "$genome")
        genome_base="${genome_base%.*}"  # Remove extension after last dot
        
        # Validate genome_base is not empty
        if [[ -z "$genome_base" ]]; then
            echo "  ERROR: Failed to extract basename from: $genome"
            continue
        fi
        
        local db_path="$DB_DIR/${genome_base}_DB"
        
        # Validate database exists
        if [[ ! -f "${db_path}.psq" ]]; then
            echo "  ERROR: Database not found for $genome_base at $db_path"
            continue
        fi
        
        echo "  Processing genome: $genome_base"
        
        for query in "${queries[@]}"; do
            [[ ! -s "$query" ]] && { echo "    Skipping missing FASTA: $query"; continue; }
            
            local query_base=$(basename "$query")
            echo "    Query: $query_base"
            
            for e_val in "${E_VALUE_LIST[@]}"; do
                for word_size in "${WORD_SIZE_LIST[@]}"; do
                    local result_dir="$output_dir/${gene_type}_blastp_results/ev_${e_val}_ws_${word_size}_${genome_base}_DB_${TIMESTAMP}"
                    mkdir -p "$result_dir"
                    
                    local csv_out="$result_dir/blastp_${genome_base}_VS_${query_base}_${e_val}_${word_size}.csv"
                    local txt_out="$result_dir/blastp_${genome_base}_VS_${query_base}_${e_val}_${word_size}.txt"
                    
                    # Check if results already exist and skip if OVERWRITE is false
                    if [[ -f "$csv_out" ]] && [[ -f "$txt_out" ]] && [[ "$OVERWRITE_BLAST_RESULTS" != true ]]; then
                        echo "      Skipping existing result: e-value=$e_val, word_size=$word_size"
                        continue
                    fi
                    
                    # CSV output with header
                    echo "Subject ID,Query ID,E-value,Percent Identity,Alignment Length,Query Length,Target Length,Query Coverage,Mismatches,Gaps,Bit Score,Max Score,Total Score" > "$csv_out"
                    blastp -query "$query" -db "$db_path" -evalue "$e_val" -word_size "$word_size" -num_threads "$NUM_THREADS" \
                           -outfmt "10 sseqid qseqid evalue pident length qlen slen qcovs mismatch gaps bitscore score score" >> "$csv_out"
                    
                    # Standard text output
                    blastp -query "$query" -db "$db_path" -evalue "$e_val" -word_size "$word_size" -num_threads "$NUM_THREADS" > "$txt_out"
                done
            done
        done
    done
}

# Organize CSV file by sorting on Subject ID, then E-value (ascending)
organize_csv() {
    local input_csv=$1
    
    # Check if input file exists and is not empty
    [[ ! -s "$input_csv" ]] && { echo "    Warning: Input CSV not found or empty: $input_csv"; return 1; }
    
    echo "    Organizing CSV: $(basename "$input_csv")"
    
    # Use the external csv_organizer.py script
    local organizer_script="$BASE_DIR/modules/csv_organizer.py"
    
    if [[ ! -f "$organizer_script" ]]; then
        echo "    Warning: csv_organizer.py not found at $organizer_script"
        return 1
    fi
    
    # Run the external Python organizer with automatic naming
    # This will create both {filename}_raw.csv and {filename}_organized.csv
    if python3 "$organizer_script" "$input_csv"; then
        echo "    CSV organized successfully with automatic naming"
    else
        echo "    Error: Failed to organize CSV"
        return 1
    fi
}

# Check if required Python dependencies are available
check_python_dependencies() {
    echo "Checking Python dependencies..."
    
    if ! command -v python3 &> /dev/null; then
        echo "  Warning: python3 not found. CSV organization will be skipped."
        return 1
    fi
    
    python3 -c "import pandas" 2>/dev/null || {
        echo "  Warning: pandas not available. CSV organization will be skipped."
        echo "  Install with: pip install pandas"
        return 1
    }
    
    local organizer_script="$BASE_DIR/modules/csv_organizer.py"
    if [[ ! -f "$organizer_script" ]]; then
        echo "  Warning: csv_organizer.py not found at $organizer_script"
        return 1
    fi
    
    echo "  Python dependencies and csv_organizer.py OK"
    return 0
}

# Merge all CSV files from BLAST results into a single file
merge_csv_results() {
    local gene_type=$1
    local output_base=$2
    local input_dir="$output_base/${gene_type}_blastp_results"
    local output_dir="$output_base/${gene_type}_curated_blastp_results"
    local raw_file="$output_dir/merged_output_${TIMESTAMP}.csv"

    mkdir -p "$output_dir"

    mapfile -t csv_files < <(find "$input_dir" -type f -name "*.csv" 2>/dev/null)

    if [[ ${#csv_files[@]} -eq 0 ]]; then
        echo "  No CSV files found for $gene_type"
        return
    fi

    # Merge all CSV files (keep header from first file, skip headers from others)
    awk 'NR==1 || FNR>1' "${csv_files[@]}" > "$raw_file"
    echo "  Merged ${#csv_files[@]} files"
    
    # Organize the merged CSV file if Python dependencies are available
    if check_python_dependencies >/dev/null 2>&1; then
        organize_csv "$raw_file"
        echo "  Organized and saved both raw and organized versions"
        
        # Delete the original merged file (without _raw or _organized suffix)
        if [[ -f "$raw_file" ]]; then
            rm "$raw_file"
            echo "  Deleted original merged file: $(basename "$raw_file")"
        fi
    else
        echo "  Saved (unsorted) -> $raw_file"
        echo "  Note: Install Python + pandas for automatic CSV sorting"
    fi
}

# Create BLAST databases for all genomes
create_blast_databases() {
    echo "Creating BLAST databases..."
    mkdir -p "$DB_DIR"
    
    for genome in "${GENOME_LIST[@]}"; do
        # Validate genome file exists
        if [[ ! -f "$genome" ]]; then
            echo "  ERROR: Genome file not found: $genome"
            continue
        fi
        
        # Extract basename and remove extension
        local genome_base=$(basename "$genome")
        genome_base="${genome_base%.*}"  # Remove extension after last dot
        
        # Validate genome_base is not empty
        if [[ -z "$genome_base" ]]; then
            echo "  ERROR: Failed to extract basename from: $genome"
            continue
        fi
        
        local db_path="$DB_DIR/${genome_base}_DB"

        if [[ -f "${db_path}.psq" ]] && [[ "$RECREATE_DATABASE" != true ]]; then
            echo "  Database exists for $genome_base"
        else
            if [[ -f "${db_path}.psq" ]] && [[ "$RECREATE_DATABASE" == true ]]; then
                echo "  Recreating database for $genome_base..."
                rm -f "${db_path}."* 2>/dev/null
            else
                echo "  Creating database for $genome_base..."
            fi
            makeblastdb -in "$genome" -dbtype prot -out "$db_path"
        fi
    done
}

# Execute BLAST analysis for a dataset type (Full or WS)
run_blast_dataset_analysis() {
    local dataset_type=$1  # "Final_GRF" or "Final_GIF"
    local output_dir=$2
    local -n all_files="${dataset_type}_GRF_GIF_FASTA_LIST"
    local -n grf_files="${dataset_type}_GRF_FASTA_LIST"
    local -n gif_files="${dataset_type}_GIF_FASTA_LIST"
    
    echo "Running $dataset_type dataset analysis..."
    
    # Run BLAST for each gene type if files exist
    [[ ${#all_files[@]} -gt 0 ]] && {
        run_blast GENOME_LIST all_files "${dataset_type}_Both" "$output_dir"
        merge_csv_results "${dataset_type}_Both" "$output_dir"
    }
    
    [[ ${#grf_files[@]} -gt 0 ]] && {
        run_blast GENOME_LIST grf_files "${dataset_type}" "$output_dir"
        merge_csv_results "${dataset_type}" "$output_dir"
    }
    
    [[ ${#gif_files[@]} -gt 0 ]] && {
        run_blast GENOME_LIST gif_files "${dataset_type}" "$output_dir"
        merge_csv_results "${dataset_type}" "$output_dir"
    }
}

# ========================================
# Main Execution
# ========================================

main() {
    # Setup logging
    setup_logging
    
    echo "=============================== GRF-GIF BLASTp Pipeline ==============================="
    echo "Timestamp: $TIMESTAMP"
    echo ""
    
    # Check Python dependencies for CSV organization
    check_python_dependencies
    echo ""
    
    # Discover FASTA files in query directories
    echo "Discovering FASTA files..."
    discover_fasta_files "$H_FINAL_PHYLO_GRF" "Final_GRF"
    discover_fasta_files "$H_FINAL_PHYLO_GIF" "Final_GIF"
    echo ""
    
    # Create BLAST databases
    create_blast_databases
    echo ""
    
    # Run analyses
    run_blast_dataset_analysis "Final_GRF" "$OUT_FINAL_GRF"
    run_blast_dataset_analysis "Final_GIF" "$OUT_FINAL_GIF"
    
    echo ""
    echo "=============================== Pipeline Completed ==============================="
    echo "Log saved to: $LOG_FILE"
}

# Execute main function
main "$@"
