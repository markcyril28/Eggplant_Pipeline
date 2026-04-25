#!/bin/bash
set -euo pipefail

# ========================================
# GRF-GIF BLASTN Analysis Pipeline
# ========================================
# This script performs BLASTN analysis of GRF and GIF gene families
# against eggplant (Solanum melongena) reference genomes.
#
# Features:
# - Automatic discovery of FASTA query files
# - Categorization by gene type (GRF, GIF, or combined)
# - Parallel processing of multiple genomes and parameters
# - CSV and text output formats
# - Result merging and organization
# - CSV sorting by Subject ID and E-value (requires Python + pandas)
#
# Dependencies:
# - BLAST+ (makeblastdb, blastn)
# - Python 3 + pandas (optional, for CSV organization)
#
# Input:  FASTA files in designated query directories
# Output: BLAST results in CSV and TXT formats, organized CSV files
# ========================================

# NOTE: This is a legacy script. The active pipeline uses 2_BLAST_Alignment.sh + modules/01_identification/blastp.sh

# ===================== IMPORTANT VARIABLES =====================
readonly PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
readonly REFSEQS_DIR="$PIPELINE_DIR/1_RefSeqs/a_Smel_RefSeqs"
readonly BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly BASE_IDENTIFICATION_DIR="$BASE_DIR/2_Identification/z_BLAST_Gene_Identification"
readonly BASE_RAW_DATA_DIR="$BASE_IDENTIFICATION_DIR/1_INPUTS"
readonly TIMESTAMP=$(date +%F_%H-%M)

# Query directories
readonly FULL_QUERY_DIR="$BASE_IDENTIFICATION_DIR/1_INPUTS/a_Full_Datasets_GRF-GIF_curated_renamed_QUERY"
readonly WS_QUERY_DIR="$BASE_IDENTIFICATION_DIR/1_INPUTS/b_Well_Studied_Datasets_GRF-GIF_curated_QUERY"
H_FINAL_PHYLO_GRF="1_INPUTS/h_Final_Phylo_Datasets/Final_Phylo_GRF_ProtSeq_with_SmelGRF_e-value_1e-5.fasta"
H_FINAL_PHYLO_GIF="1_INPUTS/h_Final_Phylo_Datasets/Final_Phylo_GIF_ProtSeq_with_SmelGIF_e-value_1e-5.fasta"

# Database and Output directories
readonly DB_DIR="$BASE_DIR/2_Identification/0_BLAST_DB"
readonly OUT_FULL="$BASE_IDENTIFICATION_DIR/3_RESULTS/a_Full_Datasets_GRF-GIF_BLAST_RESULTS"
readonly OUT_WS="$BASE_IDENTIFICATION_DIR/3_RESULTS/b_Well_Studied_Datasets_GRF-GIF_BLAST_RESULTS"
#mkdir -p $OUT_WS

# Genome references — reference from 1_RefSeqs, do not copy
readonly EGGPLANT_V3_CHR="$REFSEQS_DIR/Solanum_melongena_consortium/assembly/V3/Eggplant_V3_Chromosomes.fa"
readonly EGGPLANT_V3_CDS="$REFSEQS_DIR/Solanum_melongena_consortium/assembly/V3/Eggplant_V3.CDS.putative_function.fa"
readonly EGGPLANT_V4_1_GENOME="$REFSEQS_DIR/Solanum_melongena_v4.1/Eggplant_V4.1.fa"
readonly EGGPLANT_V4_1_TRANSCRIPTS="$REFSEQS_DIR/Solanum_melongena_v4.1/Eggplant_V4.1_transcripts.function.fa"
readonly EGGPLANT_v4_1_PROTEINS="$REFSEQS_DIR/Solanum_melongena_v4.1/Eggplant_V4.1_protein.function.fa"
# ===============================================================

# Available genomes (uncomment as needed)
GENOME_LIST=(
    #"$EGGPLANT_V4_1_TRANSCRIPTS"
    "$EGGPLANT_V4_1_PROTEINS"
    # Add other genomes here when needed
)

# BLAST parameters - Recommended for gene family identification

readonly E_VALUE_LIST=( 
    #1e-10 
    1e-15
    #1e-20   # stringent, ensures only strong homology is retained.
)

readonly WORD_SIZE_LIST=(
    11  # Smaller word size for sensitivity if looking for divergent homologs.
    #15
)

# Initialize arrays for FASTA files
declare -a FULL_GRF_GIF_FASTA_LIST FULL_GRF_FASTA_LIST FULL_GIF_FASTA_LIST
declare -a WS_GRF_GIF_FASTA_LIST WS_GRF_FASTA_LIST WS_GIF_FASTA_LIST

# ========================================
# Functions
# ========================================

# Discover and categorize FASTA files in a directory
discover_fasta_files() {
    local query_dir=$1
    local prefix=$2
    local -n all_ref="${prefix}_GRF_GIF_FASTA_LIST"
    local -n grf_ref="${prefix}_GRF_FASTA_LIST"
    local -n gif_ref="${prefix}_GIF_FASTA_LIST"
    
    # Clear arrays
    all_ref=() grf_ref=() gif_ref=()
    
    [[ ! -d "$query_dir" ]] && { echo "Directory not found: $query_dir"; return 1; }
    
    echo "Scanning: $query_dir"
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
    done < <(find "$query_dir" -type f \( -iname "*.fa" -o -iname "*.fasta" \) -print0)
    
    printf "  Found: %d total | %d GRF | %d GIF\n" "${#all_ref[@]}" "${#grf_ref[@]}" "${#gif_ref[@]}"
}

# Execute BLASTN for array of query files against genome databases
run_blast() {
    local -n genomes=$1    # Array of genome databases
    local -n queries=$2    # Array of query FASTA files
    local gene_type=$3     # Gene type label (e.g., "GRF", "GIF", "GRF-GIF")
    local output_dir=$4    # Base output directory

    echo "Running BLAST for $gene_type (${#queries[@]} queries vs ${#genomes[@]} genomes)"
    
    for genome in "${genomes[@]}"; do
        local genome_base=$(basename "$genome" | cut -d. -f1)
        local db_path="$DB_DIR/${genome_base}_DB"
        
        echo "  Processing genome: $genome_base"
        
        for query in "${queries[@]}"; do
            [[ ! -s "$query" ]] && { echo "    Skipping missing FASTA: $query"; continue; }
            
            local query_base=$(basename "$query")
            echo "    Query: $query_base"
            
            for e_val in "${E_VALUE_LIST[@]}"; do
                for word_size in "${WORD_SIZE_LIST[@]}"; do
                    local result_dir="$output_dir/${gene_type}_blastn_results/ev_${e_val}_ws_${word_size}_${genome_base}_DB_${TIMESTAMP}"
                    mkdir -p "$result_dir"
                    
                    local csv_out="$result_dir/blastn_${genome_base}_VS_${query_base}_${e_val}_${word_size}.csv"
                    local txt_out="$result_dir/blastn_${genome_base}_VS_${query_base}_${e_val}_${word_size}.txt"
                    
                    # CSV output with header
                    echo "Subject ID,Query ID,E-value,Percent Identity,Alignment Length,Query Length,Target Length,Query Coverage,Mismatches,Gaps,Bit Score,Max Score,Total Score" > "$csv_out"
                    blastn -query "$query" -db "$db_path" -evalue "$e_val" -word_size "$word_size" \
                           -outfmt "10 sseqid qseqid evalue pident length qlen slen qcovs mismatch gaps bitscore score score" >> "$csv_out"
                    
                    # Standard text output
                    blastn -query "$query" -db "$db_path" -evalue "$e_val" -word_size "$word_size" > "$txt_out"
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
    local organizer_script="$BASE_DIR/2_Identification/z_BLAST_Gene_Identification/csv_organizer.py"
    
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
    
    local organizer_script="$BASE_IDENTIFICATION_DIR/csv_organizer.py"
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
    local input_dir="$output_base/${gene_type}_blastn_results"
    local output_dir="$output_base/${gene_type}_curated_blastn_results"
    local merged_file="$output_dir/merged_output_${TIMESTAMP}.csv"
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
        local genome_base=$(basename "$genome" | cut -d. -f1)
        local db_path="$DB_DIR/${genome_base}_DB"

        if [[ -f "${db_path}.nsq" ]]; then
            echo "  Database exists for $genome_base"
        else
            echo "  Creating database for $genome_base..."
            makeblastdb -in "$genome" -dbtype nucl -out "$db_path"
        fi
    done
}

# Execute BLAST analysis for a dataset type (Full or WS)
run_blast_dataset_analysis() {
    local dataset_type=$1  # "Full" or "WS"
    local output_dir=$2
    local -n all_files="${dataset_type}_GRF_GIF_FASTA_LIST"
    local -n grf_files="${dataset_type}_GRF_FASTA_LIST"
    local -n gif_files="${dataset_type}_GIF_FASTA_LIST"
    
    echo "Running $dataset_type dataset analysis..."
    
    # Run BLAST for each gene type if files exist
    [[ ${#all_files[@]} -gt 0 ]] && {
        run_blast GENOME_LIST all_files "${dataset_type}_GRF-GIF" "$output_dir"
        merge_csv_results "${dataset_type}_GRF-GIF" "$output_dir"
    }
    
    [[ ${#grf_files[@]} -gt 0 ]] && {
        run_blast GENOME_LIST grf_files "${dataset_type}_GRF" "$output_dir"
        merge_csv_results "${dataset_type}_GRF" "$output_dir"
    }
    
    [[ ${#gif_files[@]} -gt 0 ]] && {
        run_blast GENOME_LIST gif_files "${dataset_type}_GIF" "$output_dir"
        merge_csv_results "${dataset_type}_GIF" "$output_dir"
    }
}

# ========================================
# Main Execution
# ========================================

main() {
    echo "=============================== GRF-GIF BLASTn Pipeline ==============================="
    echo "Timestamp: $TIMESTAMP"
    echo ""
    
    # Check Python dependencies for CSV organization
    check_python_dependencies
    echo ""
    
    # Discover FASTA files in query directories
    echo "Discovering FASTA files..."
    #discover_fasta_files "$FULL_QUERY_DIR" "FULL"
    discover_fasta_files "$WS_QUERY_DIR" "WS"
    echo ""
    
    # Create BLAST databases
    create_blast_databases
    echo ""
    
    # Run analyses: $1 "Full" or "WS"; $2 output directory
    #run_blast_dataset_analysis "FULL" "$OUT_FULL"
    run_blast_dataset_analysis "WS" "$OUT_WS"
    
    echo ""
    echo "=============================== Pipeline Completed ==============================="
}

# Execute main function
main "$@"
