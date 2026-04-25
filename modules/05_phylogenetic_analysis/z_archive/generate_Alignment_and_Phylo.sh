#!/bin/bash
set -euo pipefail

# ========================================
# Phylogenetic Analysis Pipeline
# ========================================
# Performs phylogenetic analysis.
# 1. Discovers and merges FASTA files from query directories
# 2. Performs multiple sequence alignment (MUSCLE, CLUSTAL, MAFFT, PROBCONS)
# 3. Constructs phylogenetic trees using MEGACC/IQTREE2
#
# Dependencies: MUSCLE, ClustalW, MAFFT, PROBCONS, MEGACC, IQ-TREE2
# ========================================

# ---------------- INPUTS ----------------
readonly INPUT_BASE_DIR="$PWD"
readonly INPUT_DIR="0_INPUT_RAW_FASTA_and_ALIGNMENT"
readonly CONFIG_DIR="2_CONFIG_FILES"

mkdir -p "$INPUT_DIR" "$CONFIG_DIR" 

INPUT_GROUP_SCRIPT=(
    #"f_Curated"
    #"g_Renamed"
    "h_Final"
)

# Accept from CLI or auto-discover later if empty
#INPUT_GROUP=()

# Alignment methods to use
readonly ALIGNMENT_METHODS=(
    "CLUSTALO"
    #"CLUSTALW"
    #"MAFFT"
    #"PROBCONS"
    #"MUSCLE"
)

# Phylogenetic software to use
readonly PHYLO_SOFTWARE=(
    "MEGA_CC_12_Ubuntu"
    #"IQTREE2"
)

readonly CONFIG_FILE=(
	#"$CONFIG_DIR/infer_ML_nucleotide_18s.mao"
    #"$CONFIG_DIR/infer_ML_nucleotide_matK_and_concat.mao"
    "$CONFIG_DIR/infer_ML_amino_acid.mao"
)

CPU=4               # Optimal Number of CPU cores to use for Phylo is 8  
#RUN_ALIGNMENT=TRUE
#RUN_PHYLO=FALSE

# ---------------- OUTPUTS ----------------
readonly OUTPUT_DIR="3_PHYLOGENETIC_TREE_RESULTS"

# ========================================================================
# LOGGING
# ========================================================================
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="${LOG_DIR:-logs}"
LOG_FILE="${LOG_FILE:-$LOG_DIR/phylo_pipeline_${RUN_ID}_full_log.log}"
#rm -rf "$LOG_DIR"/*.log

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { local level="$1"; shift; printf '[%s] [%s] %s\n' "$(timestamp)" "$level" "$*"; }
log_info() { log INFO "$@"; }
log_warn() { log WARN "$@"; }
log_error() { log ERROR "$@"; }
log_step() { log INFO "============================== $* =============================="; }

setup_logging() {
	mkdir -p "$LOG_DIR"
    
	log_choice="${log_choice:-1}"
	if [[ "$log_choice" == "2" ]]; then
		exec >"$LOG_FILE" 2>&1
	else
		exec > >(tee -a "$LOG_FILE") 2>&1
	fi
	log_info "Log: $LOG_FILE"
}

trap 'log_error "Command failed (rc=$?) at line $LINENO: ${BASH_COMMAND:-unknown}"; exit 1' ERR
trap 'log_info "Finished"' EXIT

run_with_time_to_log() {
    if [[ $# -eq 0 ]]; then
        log_error "run_with_time_to_log called with no command"
        return 1
    fi
    /usr/bin/time -v "$@" >> "$LOG_FILE" 2>&1
}

# ========================================================================
# FUNCTIONS
# ========================================================================

print_usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  -g, --group NAME         Add a single input group (can be repeated)
  -G, --groups LIST        Comma-separated list of groups
  -h, --help               Show this help

Examples:
  $(basename "$0") --group curated_21_genes_version
  $(basename "$0") --groups curated_21_genes_version,curated_64_genes_version
EOF
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -g|--group)
                [[ -n "${2:-}" ]] || { echo "Missing value for $1"; exit 2; }
                INPUT_GROUP+=("$2"); shift 2 ;;
            -G|--groups)
                [[ -n "${2:-}" ]] || { echo "Missing value for $1"; exit 2; }
                IFS=',' read -r -a __groups <<< "$2"
                for g in "${__groups[@]}"; do
                    [[ -n "$g" ]] && INPUT_GROUP+=("$g")
                done
                shift 2 ;;
            --run-alignment)
                RUN_ALIGNMENT=TRUE; shift ;;
            --skip-alignment)
                RUN_ALIGNMENT=FALSE; shift ;;
            --alignment)
                [[ -n "${2:-}" ]] || { echo "Missing value for $1 (true/false)"; exit 2; }
                case "${2,,}" in
                    TRUE|true|1|yes|on) RUN_ALIGNMENT=TRUE ;;
                    FALSE|false|0|no|off) RUN_ALIGNMENT=FALSE ;;
                    *) echo "Invalid value for --alignment: $2 (use true/false)"; exit 2 ;;
                esac
                shift 2 ;;
            --run-phylo)
                RUN_PHYLO=TRUE; shift ;;
            --skip-phylo)
                RUN_PHYLO=FALSE; shift ;;
            --phylo)
                [[ -n "${2:-}" ]] || { echo "Missing value for $1 (true/false)"; exit 2; }
                case "${2,,}" in
                    TRUE|true|1|yes|on) RUN_PHYLO=TRUE ;;
                    FALSE|false|0|no|off) RUN_PHYLO=FALSE ;;
                    *) echo "Invalid value for --phylo: $2 (use true/false)"; exit 2 ;;
                esac
                shift 2 ;;
            -h|--help)
                print_usage; exit 0 ;;
            --)
                shift; break ;;
            *)
                echo "Unknown argument: $1"
                print_usage; exit 2 ;;
        esac
    done
}

format_fasta_fold_80() {
    # Format a FASTA sequence in standard 80-character lines
    local input_file=$1
    [[ ! -f "$input_file" ]] && { log_error "File not found: $input_file"; return 1; }
    
    local temp_file="${input_file}.fmt_tmp"
    > "$temp_file"
    
    local first_sequence=true
    while IFS= read -r line; do
        if [[ "$line" =~ ^">" ]]; then
            if [[ "$first_sequence" == true ]]; then
                first_sequence=false
            else
                echo "" >> "$temp_file"
            fi
            echo "$line" >> "$temp_file"
        else
            echo "$line" | fold -w 80 >> "$temp_file"
        fi
    done < "$input_file"
    
    mv "$temp_file" "$input_file"
    log_info "Formatted: $input_file"
}

validate_fasta_sequences() {
    local file=$1
    local has_valid_sequence=false
    local current_header=""
    local current_sequence=""
    local line_count=0
    
    while IFS= read -r line; do
        ((line_count++))
        [[ -z "$line" ]] && continue
        if [[ "$line" =~ ^">" ]]; then
            if [[ -n "$current_header" ]]; then
                local clean_sequence=$(echo "$current_sequence" | tr -d '[:space:]')
                if [[ -n "$clean_sequence" ]]; then
                    has_valid_sequence=true; break
                fi
            fi
            current_header="$line"; current_sequence=""
        else
            current_sequence+="$line"
        fi
    done < "$file"
    
    if [[ -n "$current_header" ]]; then
        local clean_sequence=$(echo "$current_sequence" | tr -d '[:space:]')
        if [[ -n "$clean_sequence" ]]; then
            has_valid_sequence=true
        fi
    fi
    
    if [[ $line_count -eq 0 ]]; then
        log_warn "Empty: $file"
        return 1
    fi
    $has_valid_sequence
}

clean_merged_fasta() {
    local input_file=$1
    local temp_file="${input_file}.tmp"
    local current_header=""
    local current_sequence=""
    local entries_removed=0
    
    > "$temp_file"
    while IFS= read -r line; do
        if [[ "$line" =~ ^">" ]]; then
            if [[ -n "$current_header" ]]; then
                local clean_sequence=$(echo "$current_sequence" | tr -d '[:space:]')
                if [[ -n "$clean_sequence" ]]; then
                    echo "$current_header" >> "$temp_file"
                    echo "$current_sequence" >> "$temp_file"
                else
                    ((entries_removed++))
                fi
            fi
            current_header="$line"; current_sequence=""
        else
            current_sequence+="$line"$'\n'
        fi
    done < "$input_file"
    
    if [[ -n "$current_header" ]]; then
        local clean_sequence=$(echo "$current_sequence" | tr -d '[:space:]')
        if [[ -n "$clean_sequence" ]]; then
            echo "$current_header" >> "$temp_file"
            echo "$current_sequence" >> "$temp_file"
        else
            ((entries_removed++))
        fi
    fi
    
    mv "$temp_file" "$input_file"
    [[ $entries_removed -gt 0 ]] && log_info "Cleaned: removed $entries_removed empty sequences"
}

merge_fasta_by_gene() {
    local query_dir=$1; local prefix=$2; local gene_type=$3; local output_dir="$4"
    local output_file="$output_dir/${prefix}_Smel_${gene_type}_merged.fasta"
    [[ ! -d "$query_dir" ]] && { log_error "Directory not found: $query_dir"; return 1; }
    
    if [[ -s "$output_file" ]]; then
        log_info "$gene_type merge: SKIPPED (exists)"
        return 0
    fi
    log_step "Merging $gene_type"
    
    > "$output_file"
    local count=0
    while IFS= read -r -d '' file; do
        local filename=$(basename "$file")
        if [[ "$filename" == *"$gene_type"* ]]; then
            if [[ -s "$file" ]]; then
                if validate_fasta_sequences "$file"; then
                    [[ $count -gt 0 ]] && echo "" >> "$output_file"
                    cat "$file" >> "$output_file" && ((count++))
                fi
            fi
        fi
    done < <(find "$query_dir" -type f \( -iname "*.fa" -o -iname "*.fasta" \) -print0)
    
    log_info "Merged $count files"
    [[ -s "$output_file" ]] && clean_merged_fasta "$output_file"
}

align_sequences() {
    local input_file=$1; local method=$2; local output_dir=$3
    local basename=$(basename "$input_file" .fasta)
    local output_file="$output_dir/${basename}.fas"
    
    [[ ! -s "$input_file" ]] && { log_warn "Empty input: $input_file"; return 1; }
    if [[ -s "$output_file" ]]; then
        log_info "Align $method: SKIPPED (exists)"
        return 0
    fi

    log_step "Aligning $basename with $method"
    case "$method" in
        "MUSCLE") 
            muscle -in "$input_file" -out "$output_file" -maxiters 1000 -diags0 -threads $CPU ;;
        
        "CLUSTALO")
            run_with_time_to_log \
                clustalo -i "$input_file" -o "$output_file" --outfmt=fasta ;;

        "CLUSTALW")
            run_with_time_to_log \
                clustalw -INFILE="$input_file" -OUTFILE="$output_file" -OUTPUT=FASTA ;;

        #"CLUSTAL") 
        #    clustalo -i "$input_file" -o "$output_file" --outfmt=fasta \
        #        --full --full-iter --iter=10 \
        #        --max-guidetree-iterations=10 --max-hmm-iterations=10 \
        #        --threads $CPU ;;
        
        "MAFFT") 
            mafft --thread $CPU --localpair --maxiterate 1000 "$input_file" > "$output_file" ;;
        
        "PROBCONS") 
            probcons -c 5 -ir 1000 -pre 20 "$input_file" > "$output_file" ;; 
        *) 
            log_error "Unknown alignment method: $method"; return 1 ;;
    esac

    log_info "Output alignment: $output_file"
}

generate_MEGA_CC_12_HPC_Docker_tree() {
    local aligned_file=$1
    local method=$2
    local config_file=$3
    local output_dir=$4

    local basename=$(basename "$aligned_file" .fas)
    local config_base=$(basename "$config_file" .mao)
    local tree_dir="$output_dir/${method}_aligned/MEGA12_HPC_Docker"
    local output_file="$tree_dir/${basename}_${config_base}.nwk"
    local mega_log="$tree_dir/${basename}_MEGA.log"

    mkdir -p "$tree_dir"

    # Pre-checks
    if [[ ! -s "$aligned_file" ]]; then
        log_warn "Aligned file empty: $aligned_file"
        return 1
    fi
    if [[ ! -f "$config_file" ]]; then
        log_error "Config file not found: $config_file"
        return 1
    fi

    # Skip if already generated
    if [[ -s "$output_file" ]]; then
        log_info "Tree already exists: $output_file (skipped)"
        return 0
    fi

    log_info "Generating MEGA tree for $(basename "$aligned_file") | Aligned with $method | Config: $(basename "$config_file")"

    # Convert paths to absolute paths for Docker mounting
    local abs_aligned_file=$(realpath "$aligned_file")
    local abs_config_file=$(realpath "$config_file")
    local abs_tree_dir=$(realpath "$tree_dir")

    # Extract directory paths for Docker volume mounting
    local data_dir=$(dirname "$abs_aligned_file")
    local config_dir=$(dirname "$abs_config_file")

    # Run MEGA via Docker with volume mounts
    sudo docker run --rm \
        -v "$data_dir:/data" \
        -v "$config_dir:/config" \
        -v "$abs_tree_dir:/output" \
        megacc:latest \
        megacc \
        -d "/data/$(basename "$abs_aligned_file")" \
        -a "/config/$(basename "$abs_config_file")" \
        -o "/output/$(basename "$output_file")" \
        > "$mega_log" 2>&1
        #--cpu $CPU \

    if [[ -s "$output_file" ]]; then
        log_info "✅ Tree: $output_file"
    else
        log_error "MEGA failed (see $mega_log)"
        return 1
    fi
}

generate_MEGA_CC_12_Ubuntu_tree() {
    local aligned_file=$1
    local method=$2
    local config_file=$3
    local output_dir=$4

    local basename=$(basename "$aligned_file" .fas)
    local config_base=$(basename "$config_file" .mao)
    local tree_dir="$output_dir/${method}_aligned/MEGA12_Ubuntu"
    local output_file="$tree_dir/${basename}_${config_base}.nwk"
    local mega_log="$tree_dir/${basename}_MEGA.log"

    mkdir -p "$tree_dir"
    touch "$mega_log"

    # Pre-checks
    if [[ ! -s "$aligned_file" ]]; then
        log_warn "Aligned file empty: $aligned_file"
        return 1
    fi
    if [[ ! -f "$config_file" ]]; then
        log_error "Config file not found: $config_file"
        return 1
    fi

    # Skip if already generated
    if [[ -s "$output_file" ]]; then
        log_info "Tree already exists: $output_file (skipped)"
        return 0
    fi

    log_info "Generating MEGA tree for $(basename "$aligned_file") | Aligned with $method | Config: $(basename "$config_file")"

    # Run MEGA with timing and log output
    run_with_time_to_log \
        megacc \
            -d "$aligned_file" \
            -a "$config_file" \
            -o "$output_file" \
            --cpu $CPU \
            > "$mega_log" 2>&1

    if [[ -s "$output_file" ]]; then
        log_info "✅ Tree: $output_file"
    else
        log_error "MEGA12 failed (see $mega_log)"
        return 1
    fi
}

generate_IQTREE2_tree() {
    local aligned_file=$1
    local method=$2
    local output_dir=$3

    local basename=$(basename "$aligned_file" .fas)
    local tree_dir="$output_dir/${method}_aligned/IQTREE2"
    local output_prefix="$tree_dir/${basename}_IQTREE2"
    local tree_file="${output_prefix}.treefile"
    local log_file="${output_prefix}.log"

    mkdir -p "$tree_dir"
    touch "$log_file"

    # Pre-checks
    if [[ ! -s "$aligned_file" ]]; then
        log_warn "Empty: $aligned_file"
        return 1
    fi
    if ! command -v iqtree &>/dev/null; then
        log_error "iqtree not in PATH"
        return 1
    fi

    # Skip if already generated
    if [[ -s "$tree_file" ]]; then
        log_info "IQ-TREE2: SKIPPED (exists)"
        return 0
    fi

    log_step "IQ-TREE2: $basename | $method"

    # Run IQ-TREE2 with timing and bootstrap support
    run_with_time_to_log \
        iqtree2 \
            -s "$aligned_file" \
            -m MFP+MERGE \
            -nt AUTO \
            --nstop 500 --nbest 20 \
            --allnni \
            --polytomy --safe \
            --boot-trees 2000 \
            -bb 5000 -alrt 2000 \
            --redo \
            -pre "$output_prefix"
            > "$log_file" 2>&1

    if [[ -s "$tree_file" ]]; then
        log_info "✅ Tree: $tree_file"
    else
        log_error "IQ-TREE2 failed (see $log_file)"
        return 1
    fi
}

generate_IQTREE2_tree_Full_Accuracy() {
    local aligned_file=$1
    local method=$2
    local output_dir=$3

    local basename=$(basename "$aligned_file" .fas)
    local tree_dir="$output_dir/${method}_aligned/IQTREE2"
    local output_prefix="$tree_dir/${basename}_IQTREE2"
    local tree_file="${output_prefix}.treefile"
    local log_file="${output_prefix}.log"

    mkdir -p "$tree_dir"
    touch "$log_file"

    # Pre-checks
    if [[ ! -s "$aligned_file" ]]; then
        log_warn "Empty: $aligned_file"
        return 1
    fi
    if ! command -v iqtree &>/dev/null; then
        log_error "iqtree not in PATH"
        return 1
    fi

    # Skip if already generated
    if [[ -s "$tree_file" ]]; then
        log_info "IQ-TREE2: SKIPPED (exists)"
        return 0
    fi

    log_step "IQ-TREE2: $basename | $method"
    log_info "IQ-TREE2 version: $(iqtree --version)"

    # Run IQ-TREE2 with timing and bootstrap support
    run_with_time_to_log \
        iqtree \
            -s "$aligned_file" \
            -m MFP+MERGE -st AA \
            -nt AUTO \
            --allnni \
            --polytomy --safe \
            -bb 40000 -alrt 10000 --lbp 10000 \
            --bnni \
            --sampling GENESITE \
            --redo \
            -pre "$output_prefix"
            > "$log_file" 2>&1

    if [[ -s "$tree_file" ]]; then
        log_info "✅ Tree: $tree_file"
    else
        log_error "IQ-TREE2 failed (see $log_file)"
        return 1
    fi
}

generate_IQTREE2_tree_Standard_Bootstrap() {
    local aligned_file=$1
    local method=$2
    local output_dir=$3

    local basename=$(basename "$aligned_file" .fas)
    local tree_dir="$output_dir/${method}_aligned/IQTREE2"
    local output_prefix="$tree_dir/${basename}_IQTREE2"
    local tree_file="${output_prefix}.treefile"
    local log_file="${output_prefix}.log"

    mkdir -p "$tree_dir"
    touch "$log_file"

    # Pre-checks
    if [[ ! -s "$aligned_file" ]]; then
        log_warn "Empty: $aligned_file"
        return 1
    fi
    if ! command -v iqtree &>/dev/null; then
        log_error "iqtree not in PATH"
        return 1
    fi

    # Skip if already generated
    if [[ -s "$tree_file" ]]; then
        log_info "IQ-TREE2: SKIPPED (exists)"
        return 0
    fi

    log_step "IQ-TREE2: $basename | $method"
    log_info "IQ-TREE2 version: $(iqtree --version)"

    # Run IQ-TREE2 with timing and bootstrap support
    run_with_time_to_log \
        iqtree \
            -s "$aligned_file" \
            -m MFP+MERGE -st AA \
            -nt AUTO \
            --nstop 1000 --nbest 50 --allnni \
            --polytomy --safe \
            -b 1000 -alrt 1000 \
            --redo \
            -pre "$output_prefix"
            > "$log_file" 2>&1

    if [[ -s "$tree_file" ]]; then
        log_info "✅ Tree: $tree_file"
    else
        log_error "IQ-TREE2 failed (see $log_file)"
        return 1
    fi
}

# ========================================================================
# MAIN
# ========================================================================
main() {
    setup_logging
    parse_args "$@"

    # Use the INPUT_GROUP_SCRIPT if none were provided: choose subfolders in INPUT_DIR that contain a b_RAW directory
    if [[ ${#INPUT_GROUP[@]} -eq 0 ]]; then
        INPUT_GROUP=$INPUT_GROUP_SCRIPT
    fi

    if [[ ${#INPUT_GROUP[@]} -eq 0 ]]; then
        log_error "No input groups provided and none found in '$INPUT_DIR'."
        log_info "Use --group NAME or --groups a,b,c"
        exit 1
    fi

    log_step "Starting Phylogenetic Analysis Pipeline"
    log_info "Input groups: ${INPUT_GROUP[*]}"

    for group in "${INPUT_GROUP[@]}"; do
        local query_dir="$INPUT_DIR/$group"
        local output_subdir="$OUTPUT_DIR/$group"
        mkdir -p "$query_dir/b_RAW" "$output_subdir"

        if [ "$RUN_ALIGNMENT" = TRUE ]; then
            log_step "Step 2: Sequence Alignments for $group"
            for b_RAW_file in "$query_dir/b_RAW/"*.fasta; do
                [[ ! -f "$b_RAW_file" ]] && continue
                format_fasta_fold_80 "$b_RAW_file"
                for align_method in "${ALIGNMENT_METHODS[@]}"; do
                    mkdir -p "$query_dir/c_ALIGNMENT/${align_method}_aligned"
                    align_sequences "$b_RAW_file" "$align_method" "$query_dir/c_ALIGNMENT/${align_method}_aligned"
                done
            done
        else
            log_warn "Skipping alignment (RUN_ALIGNMENT=FALSE)"
        fi

        if [ "$RUN_PHYLO" = TRUE ]; then
    log_step "Step 3: Phylogenetic Trees for $group"

    for align_method in "${ALIGNMENT_METHODS[@]}"; do
        aligned_files=("$query_dir/c_ALIGNMENT/${align_method}_aligned/"*.fas)
        #aligned_files=("$query_dir/c_ALIGNMENT/${align_method}_aligned/concatenated_sequences.fas")
        for aligned_file in "${aligned_files[@]}"; do
            [[ ! -f "$aligned_file" ]] && continue
                
            for software in "${PHYLO_SOFTWARE[@]}"; do
                case "$software" in
                    "MEGA_CC_12_Ubuntu")
                        log_step "$software"
                        generate_MEGA_CC_12_Ubuntu_tree "$aligned_file" "$align_method" "$config_file" "$output_subdir"
                        ;;
                    "IQTREE2")
                        log_step "$software"
                        #generate_IQTREE2_tree "$aligned_file" "$align_method" "$output_subdir"
                        generate_IQTREE2_tree_Full_Accuracy "$aligned_file" "$align_method" "$output_subdir"
                        #generate_IQTREE2_tree_Standard_Bootstrap "$aligned_file" "$align_method" "$output_subdir"
                        ;;
                    *)
                        log_error "Unknown software: $software"
                        ;;
                esac
            done
        done
    done

else
    log_warn "Skipping phylogenetic tree generation (RUN_PHYLO=FALSE)"
fi

    done

    log_step "Pipeline Completed"
}

main "$@"
