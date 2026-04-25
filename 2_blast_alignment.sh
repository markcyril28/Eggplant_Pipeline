#!/bin/bash
# ============================================================================
# Program 2: BLAST Identification & Ortholog Alignment
# ============================================================================
# Edit gene_groups in 2_blast_alignmentCONFIG.toml, then run:
#   bash b_blast_alignment.sh
# ============================================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"

# ===================== IMPORTANT VARIABLES =====================
# Gene groups: edit [pipeline].gene_groups in 2_blast_alignmentCONFIG.toml
# ===============================================================
SHARED_CONFIG="$PIPELINE_DIR/2_blast_alignmentCONFIG.toml"
mapfile -t GENE_GROUPS < <(
    python3 "$MODULES/utils/parse_toml.py" "$SHARED_CONFIG" pipeline gene_groups 2>/dev/null
)
[[ ${#GENE_GROUPS[@]} -eq 0 ]] && { echo "ERROR: No gene_groups in 2_blast_alignmentCONFIG.toml [pipeline].gene_groups" >&2; exit 1; }

# Set PROJECT_ROOT before sourcing logging utility so log paths are absolute
PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"

# Source the logging utility — setup_logging() creates all subdirectories
source "$MODULES/logging/logging_utils.sh"


TOML_PARSER="$MODULES/utils/parse_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

TEMP_FILES=()
cleanup_all() {
    teardown_logging 2>/dev/null
    rm -f "${TEMP_FILES[@]}"
}
trap cleanup_all EXIT

# parse_toml emits TOML arrays one item per line.
# Output directly for mapfile consumers.
read_toml_list() {
    local section="$1"
    local key="$2"
    local raw

    raw=$(get_toml "$section" "$key" 2>/dev/null || true)
    [[ -z "${raw// }" ]] && return 0

    local items=()
    mapfile -t items <<< "$raw"
    printf '%s\n' "${items[@]}"
}

# Wait until a parallel slot opens before launching the next job
# Parameterized form (preferred) — pass the concurrency limit explicitly.
wait_for_slot() { local limit="$1"; while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done; }

# ---------------------------------------------------------------------------
# resolve_top_folder <target_path>
#   Route output to a genome-specific first-level folder when possible.
# ---------------------------------------------------------------------------
resolve_top_folder() {
    local target_path="$1"
    case "$target_path" in
        *Solanum_melongena_v4.1*|*smel_v4_1*|*Eggplant*)
            echo "Solanum_melongena_v4.1" ;;
        *GPE001970*|*unito*)
            echo "GPE001970_SMEL5" ;;
        *)
            echo "shared" ;;
    esac
}

# Resolve a configured path against the pipeline root and return it only when it exists.
resolve_configured_path() {
    local configured_path="$1"

    if [[ -f "$configured_path" ]]; then
        echo "$configured_path"
        return 0
    fi

    if [[ "$configured_path" != /* && ! "$configured_path" =~ ^[A-Za-z]:[\\/] ]]; then
        local pipeline_path="$PIPELINE_DIR/$configured_path"
        if [[ -f "$pipeline_path" ]]; then
            echo "$pipeline_path"
            return 0
        fi
    fi

    # Fallback for relocated DMP query FASTAs.
    if [[ "$configured_path" == II_INPUTS/DMP/query_fasta/* ]]; then
        local rel_path="${configured_path#II_INPUTS/DMP/query_fasta/}"
        local copied_input_path="$PIPELINE_DIR/II_INPUTS/DMP_query_fasta_file/$rel_path"
        if [[ -f "$copied_input_path" ]]; then
            log_info "Using copied DMP query FASTA source: $configured_path -> ${copied_input_path#$PIPELINE_DIR/}" >&2
            echo "$copied_input_path"
            return 0
        fi

        local archived_path="$PIPELINE_DIR/z_archive/3_RESULT_v1/DMP/02_BLAST_Alignment/BLAST_Result/DMP_query_fasta_file/$rel_path"
        if [[ -f "$archived_path" ]]; then
            log_warn "Using archived DMP query FASTA fallback: $configured_path -> ${archived_path#$PIPELINE_DIR/}" >&2
            echo "$archived_path"
            return 0
        fi
    fi

    return 1
}

# Resolve configured path to absolute/canonical pipeline path without existence checks.
resolve_pipeline_path() {
    local configured_path="$1"

    if [[ "$configured_path" == /* || "$configured_path" =~ ^[A-Za-z]:[\\/] ]]; then
        echo "$configured_path"
    else
        echo "$PIPELINE_DIR/$configured_path"
    fi
}

# Auto-generate configured query_protein_fastas if any are missing.
ensure_query_protein_fastas() {
    local converter_script="$MODULES/02_blast_ortholog_alignment/convert_query_fastas_to_proteins.sh"

    mapfile -t CONFIGURED_PROTEIN_QUERIES < <(read_toml_list ortholog_blast query_protein_fastas)
    if [[ ${#CONFIGURED_PROTEIN_QUERIES[@]} -eq 0 ]]; then
        log_info "No query_protein_fastas configured; skipping protein query FASTA generation."
        return 0
    fi

    local missing_count=0
    local protein_query_abs
    for protein_query in "${CONFIGURED_PROTEIN_QUERIES[@]}"; do
        protein_query_abs="$(resolve_pipeline_path "$protein_query")"
        if [[ ! -f "$protein_query_abs" ]]; then
            ((++missing_count))
        fi
    done

    if (( missing_count == 0 )); then
        log_info "All configured protein query FASTAs already exist (${#CONFIGURED_PROTEIN_QUERIES[@]} files)."
        return 0
    fi

    if [[ ! -f "$converter_script" ]]; then
        log_error "Missing converter script: $converter_script"
        return 1
    fi

    log_warn "$missing_count configured protein query FASTA file(s) are missing; generating from query_fastas..."
    if ! bash "$converter_script" --config "$CONFIG_FILE" --threads "$THREADS_PER_JOB"; then
        log_error "Protein query FASTA generation failed."
        return 1
    fi

    # Verify outputs after generation.
    local still_missing=0
    for protein_query in "${CONFIGURED_PROTEIN_QUERIES[@]}"; do
        protein_query_abs="$(resolve_pipeline_path "$protein_query")"
        if [[ ! -f "$protein_query_abs" ]]; then
            ((++still_missing))
        fi
    done

    if (( still_missing > 0 )); then
        log_error "$still_missing configured protein query FASTA file(s) are still missing after conversion."
        return 1
    fi

    log_info "Protein query FASTA auto-generation completed successfully."
}

for GENE_GROUP in "${GENE_GROUPS[@]}"; do

    # Resolve config: split directory or monolithic file
    CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
    if [[ -d "$CONFIG_DIR" ]]; then
        CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_blast_cfg_XXXXXX.toml")
        TEMP_FILES+=("$CONFIG_FILE")
        python3 "$MODULES/utils/merge_toml.py" \
            "$PIPELINE_DIR/2_blast_alignmentCONFIG.toml" \
            "$CONFIG_DIR/00_common.toml" \
            "$CONFIG_DIR/01_hmmer_gene_identification.toml" \
            "$CONFIG_DIR/02_blast_ortholog_alignment.toml" > "$CONFIG_FILE"
    else
        CONFIG_FILE="$PIPELINE_DIR/config/${GENE_GROUP}.toml"
    fi

    setup_logging
    log_step "Processing gene group: $GENE_GROUP"
    log_info "Pipeline started at: $(date)"
    log_info "Processing gene group: $GENE_GROUP with ${#GENE_GROUPS[@]} total gene groups"

BASE_DIR="$PIPELINE_DIR/$(get_toml general base_dir)"
MACHINE=$(get_toml pipeline machine 2>/dev/null || echo "Local")
CPU=$(get_toml pipeline compute "$MACHINE" threads 2>/dev/null || nproc)
BLAST_OPTIMAL_THREADS=$(get_toml pipeline compute "$MACHINE" blast_optimal_threads 2>/dev/null || echo "4")
MAX_PARALLEL=$(get_toml pipeline compute "$MACHINE" max_parallel 2>/dev/null || echo "$(( CPU / BLAST_OPTIMAL_THREADS ))")
(( MAX_PARALLEL < 1 )) && MAX_PARALLEL=1
THREADS_PER_JOB=$BLAST_OPTIMAL_THREADS
OVERWRITE=$(get_toml pipeline overwrite 2>/dev/null || echo "true")
CLEAR_OUTPUT=$(get_toml pipeline clear_output 2>/dev/null || echo "false")
BLASTN_EVALUE=$(get_toml ortholog_blast blastn_e_value 2>/dev/null || echo "1e-10")
BLASTN_WORD_SIZE=$(get_toml ortholog_blast blastn_word_size 2>/dev/null || echo "11")
BLASTN_MAX_TARGET_SEQS=$(get_toml ortholog_blast blastn_max_target_seqs 2>/dev/null || echo "500")
BLASTN_MAX_HSPS=$(get_toml ortholog_blast blastn_max_hsps 2>/dev/null || echo "1")
BLASTN_PERC_IDENTITY=$(get_toml ortholog_blast blastn_perc_identity 2>/dev/null || echo "0")
BLASTN_QCOV_HSP_PERC=$(get_toml ortholog_blast blastn_qcov_hsp_perc 2>/dev/null || echo "0")
BLASTN_TASK=$(get_toml ortholog_blast blastn_task 2>/dev/null || echo "blastn")
BLASTN_DUST=$(get_toml ortholog_blast blastn_dust 2>/dev/null || echo "yes")
BLASTX_EVALUE=$(get_toml ortholog_blast blastx_e_value 2>/dev/null || echo "1e-15")
BLASTX_WORD_SIZE=$(get_toml ortholog_blast blastx_word_size 2>/dev/null || echo "3")
BLASTX_MAX_TARGET_SEQS=$(get_toml ortholog_blast blastx_max_target_seqs 2>/dev/null || echo "500")
BLASTX_MAX_HSPS=$(get_toml ortholog_blast blastx_max_hsps 2>/dev/null || echo "1")
BLASTX_QCOV_HSP_PERC=$(get_toml ortholog_blast blastx_qcov_hsp_perc 2>/dev/null || echo "0")
BLASTX_MATRIX=$(get_toml ortholog_blast blastx_matrix 2>/dev/null || echo "BLOSUM62")
BLASTX_SEG=$(get_toml ortholog_blast blastx_seg 2>/dev/null || echo "yes")
BLASTP_EVALUE=$(get_toml ortholog_blast blastp_e_value 2>/dev/null || echo "1e-15")
BLASTP_WORD_SIZE=$(get_toml ortholog_blast blastp_word_size 2>/dev/null || echo "2")
BLASTP_MAX_TARGET_SEQS=$(get_toml ortholog_blast blastp_max_target_seqs 2>/dev/null || echo "500")
BLASTP_MAX_HSPS=$(get_toml ortholog_blast blastp_max_hsps 2>/dev/null || echo "1")
BLASTP_QCOV_HSP_PERC=$(get_toml ortholog_blast blastp_qcov_hsp_perc 2>/dev/null || echo "0")
BLASTP_MATRIX=$(get_toml ortholog_blast blastp_matrix 2>/dev/null || echo "BLOSUM62")
BLASTP_SEG=$(get_toml ortholog_blast blastp_seg 2>/dev/null || echo "yes")
FORCE_REINDEX=$(get_toml ortholog_blast force_reindex 2>/dev/null || echo "false")

BLAST_DIR="$BASE_DIR/02_BLAST_Alignment"

# ── Clear output directories ──────────────────────────────────────────────────
if [[ "$CLEAR_OUTPUT" == "true" || "$CLEAR_OUTPUT" == "True" ]]; then
    log_warn "clear_output=true — deleting all existing BLAST output for $GENE_GROUP"
    if [[ -d "$BLAST_DIR" ]]; then
        rm -rf "$BLAST_DIR"
        log_info "  Cleared: $BLAST_DIR"
    fi
fi

log_step "BLAST Alignment: $GENE_GROUP"
log_info "Base directory: $BASE_DIR"
log_info "BLAST parameters:"
log_info "  BLASTn: e=$BLASTN_EVALUE ws=$BLASTN_WORD_SIZE task=$BLASTN_TASK dust=$BLASTN_DUST"
log_info "  BLASTn: max_target_seqs=$BLASTN_MAX_TARGET_SEQS max_hsps=$BLASTN_MAX_HSPS perc_identity=$BLASTN_PERC_IDENTITY qcov=$BLASTN_QCOV_HSP_PERC"
log_info "  BLASTx: e=$BLASTX_EVALUE ws=$BLASTX_WORD_SIZE matrix=$BLASTX_MATRIX seg=$BLASTX_SEG"
log_info "  BLASTx: max_target_seqs=$BLASTX_MAX_TARGET_SEQS max_hsps=$BLASTX_MAX_HSPS qcov=$BLASTX_QCOV_HSP_PERC"
log_info "  BLASTp: e=$BLASTP_EVALUE ws=$BLASTP_WORD_SIZE matrix=$BLASTP_MATRIX seg=$BLASTP_SEG"
log_info "  BLASTp: max_target_seqs=$BLASTP_MAX_TARGET_SEQS max_hsps=$BLASTP_MAX_HSPS qcov=$BLASTP_QCOV_HSP_PERC"
log_info "  Force Reindex: $FORCE_REINDEX"
log_info "Threading: CPU=$CPU  THREADS_PER_JOB=$THREADS_PER_JOB  MAX_PARALLEL=$MAX_PARALLEL  (machine=$MACHINE)"

# Build makeblastdb flags
MAKEDB_FLAGS=""
if [[ "$FORCE_REINDEX" == "true" || "$FORCE_REINDEX" == "True" ]]; then
    MAKEDB_FLAGS="--force"
    log_info "makeblastdb will use --force flag"
fi

log_configuration  # Log pipeline configuration settings
catalog_all_software \
    "blastn:blastn -version" \
    "blastp:blastp -version" \
    "blastx:blastx -version" \
    "makeblastdb:makeblastdb -version" \
    "python:python3 --version"

# ── Operations ────────────────────────────────────────────────────────────────
# Read the operations list from TOML.  Default: run everything when absent.
# Each BLAST type (blastn, blastx, blastp) is now a first-class operation.
# Legacy "blast_alignment" enables all three BLAST types for backward compat.
_ops=$(get_toml ortholog_blast operations 2>/dev/null | tr '\n' ' ' || true)
if [[ -z "${_ops// }" ]]; then
    DO_BLASTN=true; DO_BLASTX=true; DO_BLASTP=true
    DO_VIZ_HEATMAP=true; DO_VIZ_HEATMAP_EVALUE=true; DO_VIZ_LOLLIPOP=true
else
    # Legacy: "blast_alignment" enables all three BLAST programs
    if [[ "$_ops" == *"blast_alignment"* ]]; then
        DO_BLASTN=true; DO_BLASTX=true; DO_BLASTP=true
    else
        [[ " $_ops " == *" blastn "* ]]  && DO_BLASTN=true || DO_BLASTN=false
        [[ " $_ops " == *" blastx "* ]]  && DO_BLASTX=true || DO_BLASTX=false
        [[ " $_ops " == *" blastp "* ]]  && DO_BLASTP=true || DO_BLASTP=false
    fi
    # visualize_blastn (legacy) enables all three figure sub-steps
    if [[ "$_ops" == *"visualize_blastn"* ]]; then
        DO_VIZ_HEATMAP=true; DO_VIZ_HEATMAP_EVALUE=true; DO_VIZ_LOLLIPOP=true
    else
        [[ " $_ops " == *" visualize_heatmap_evalue "* ]] && DO_VIZ_HEATMAP_EVALUE=true || DO_VIZ_HEATMAP_EVALUE=false
        [[ " $_ops " == *" visualize_heatmap "*        ]] && DO_VIZ_HEATMAP=true        || DO_VIZ_HEATMAP=false
        [[ " $_ops " == *" visualize_lollipop "*       ]] && DO_VIZ_LOLLIPOP=true       || DO_VIZ_LOLLIPOP=false
    fi
fi
DO_VISUALIZE=false
[[ "$DO_VIZ_HEATMAP" == "true" || "$DO_VIZ_HEATMAP_EVALUE" == "true" || "$DO_VIZ_LOLLIPOP" == "true" ]] && DO_VISUALIZE=true
log_info "Operations — BLASTn: $DO_BLASTN | BLASTx: $DO_BLASTX | BLASTp: $DO_BLASTP | visualize_heatmap: $DO_VIZ_HEATMAP | visualize_heatmap_evalue: $DO_VIZ_HEATMAP_EVALUE | visualize_lollipop: $DO_VIZ_LOLLIPOP"

# ---- BLASTn Identification ----
METHOD=$(get_toml identification method 2>/dev/null || echo "")
if [[ "$METHOD" == "blastn" && "$DO_BLASTN" == "true" ]]; then
    log_step "BLASTn Identification: $GENE_GROUP"
    log_info "Starting BLASTn identification for $GENE_GROUP"

    IDENT_DIR="$BASE_DIR/01_Identification"
    mkdir -p "$IDENT_DIR"
    log_info "Created identification directory: $IDENT_DIR"

    # Clean stale BLASTn identification results to prevent ev_* folder accumulation
    for _gdir in "$IDENT_DIR"/*/; do
        [[ -d "${_gdir}blastn_results" ]] && { log_info "Cleaning stale BLASTn ident results in $(basename "${_gdir%/}")"; rm -rf "${_gdir}blastn_results"; }
    done

    # Build nucl BLAST databases in parallel and collect DB paths
    DB_PATHS=()
    DB_PIDS=()
    log_info "Building nucleotide BLAST databases for target genomes..."
    mapfile -t _target_genomes < <(get_toml ortholog_blast target_genomes)
    for genome_path in "${_target_genomes[@]}"; do
        GENOME_FULL="$PIPELINE_DIR/$genome_path"
        GENOME_BASENAME=$(basename "$GENOME_FULL")
        log_info "Building database for genome: $GENOME_BASENAME"
        wait_for_slot "$MAX_PARALLEL"
        bash "$MODULES/02_blast_ortholog_alignment/makeblastdb.sh" \
            --input "${GENOME_FULL}.fa" --dbtype nucl --outdir "$(dirname "$GENOME_FULL")" $MAKEDB_FLAGS &
        DB_PIDS+=("$!")
        DB_PATHS+=("$(dirname "$GENOME_FULL")/${GENOME_BASENAME}_db/$GENOME_BASENAME")
    done
    FAILED=0
    for pid in "${DB_PIDS[@]}"; do wait "$pid" || ((++FAILED)); done
    if (( FAILED > 0 )); then
        log_error "$FAILED database build(s) failed"
        exit 1
    else
        log_info "Successfully built ${#DB_PIDS[@]} nucleotide databases"
    fi

    # Run BLASTn queries against ALL databases in parallel
    BLAST_PIDS=()
    log_info "Running BLASTn queries against ${#DB_PATHS[@]} databases..."
    mapfile -t _query_fastas < <(get_toml ortholog_blast query_fastas)
    for DB_PATH in "${DB_PATHS[@]}"; do
        for query in "${_query_fastas[@]}"; do
            QUERY_FULL="$PIPELINE_DIR/$query"
            [[ -f "$QUERY_FULL" ]] || { log_warn "Query not found: $QUERY_FULL"; continue; }
            
            # Determine genome-specific output directory
            DB_DIR=$(dirname "$DB_PATH")
            DB_NAME=$(basename "$DB_DIR")
            TOP_FOLDER=$(resolve_top_folder "$DB_NAME")
            GENOME_BLAST_DIR="$IDENT_DIR/$TOP_FOLDER"
            mkdir -p "$GENOME_BLAST_DIR"
            
            log_info "Running BLASTn: query=$(basename "$QUERY_FULL") vs db=$(basename "$DB_DIR")"
            wait_for_slot "$MAX_PARALLEL"
            bash "$MODULES/02_blast_ortholog_alignment/blastn.sh" \
                --query "$QUERY_FULL" \
                --db "$DB_PATH" \
                --evalue "$BLASTN_EVALUE" \
                --word-size "$BLASTN_WORD_SIZE" \
                --max-target-seqs "$BLASTN_MAX_TARGET_SEQS" \
                --max-hsps "$BLASTN_MAX_HSPS" \
                --perc-identity "$BLASTN_PERC_IDENTITY" \
                --qcov-hsp-perc "$BLASTN_QCOV_HSP_PERC" \
                --task "$BLASTN_TASK" \
                --dust "$BLASTN_DUST" \
                --outdir "$GENOME_BLAST_DIR" \
                --threads "$THREADS_PER_JOB" &
            BLAST_PIDS+=("$!")
        done
    done

    # Wait and check for failures
    log_info "Waiting for ${#BLAST_PIDS[@]} BLASTn jobs to complete..."
    FAILED=0
    for pid in "${BLAST_PIDS[@]}"; do
        wait "$pid" || ((++FAILED))
    done
    if (( FAILED > 0 )); then
        log_error "$FAILED BLASTn job(s) failed"
        exit 1
    else
        log_info "All BLASTn jobs completed successfully"
    fi

    # Create curated results per genome directory (3 versions: hmmer_only, combined, plant_only)
    for genome_dir in "$IDENT_DIR"/*/; do
        genome_dir="${genome_dir%/}"
        if [[ -d "$genome_dir" && "$genome_dir" != "$IDENT_DIR/BLAST_DB" && "$genome_dir" != "$IDENT_DIR/curated_results" ]]; then
            genome_name=$(basename "$genome_dir")
            mkdir -p "$genome_dir/curated_results"
            # Clean old merged BLASTn CSVs to prevent stale file accumulation
            rm -f "$genome_dir/curated_results"/merged_blastn_*.csv
            if [[ -d "$genome_dir/blastn_results" ]]; then
                shopt -s nullglob; _csvs=("$genome_dir/blastn_results"/ev_*/*.csv); shopt -u nullglob
                if [[ ${#_csvs[@]} -eq 0 ]]; then
                    log_warn "No BLASTn CSVs found for genome: $genome_name"
                else
                    DATE_TAG=$(date +%F)
                    log_info "Merging BLASTn results (3 versions) for genome: $genome_name"
                    bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                        --input-dir "$genome_dir/blastn_results" \
                        --include-pattern "hmmer_identified" \
                        --output "$genome_dir/curated_results/merged_blastn_${DATE_TAG}_${genome_name}_hmmer_only.csv"
                    bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                        --input-dir "$genome_dir/blastn_results" \
                        --output "$genome_dir/curated_results/merged_blastn_${DATE_TAG}_${genome_name}_combined.csv"
                    bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                        --input-dir "$genome_dir/blastn_results" \
                        --exclude-pattern "hmmer_identified" \
                        --output "$genome_dir/curated_results/merged_blastn_${DATE_TAG}_${genome_name}_plant_only.csv"
                fi
            fi
        fi
    done

    log_step "BLASTn Identification complete: $GENE_GROUP"
fi

# ---- BLASTx Ortholog Alignment (translated nucleotide query vs protein DB) ----
ORTHO_ENABLED=$(get_toml ortholog_blast enabled 2>/dev/null || echo "false")
USE_HMMER=$(get_toml ortholog_blast use_hmmer_output 2>/dev/null || echo "false")
PROT_DB_PATHS=()   # initialised here so BLASTp can reuse DBs built by BLASTx
if [[ "$DO_BLASTX" != "true" ]] || [[ "$ORTHO_ENABLED" != "True" && "$ORTHO_ENABLED" != "true" ]]; then
    log_info "BLASTx Ortholog Alignment is not enabled for $GENE_GROUP. Skipping BLASTx section."
else

log_info "BLASTx Ortholog Alignment is enabled for $GENE_GROUP"

BLAST_DIR="$BASE_DIR/02_BLAST_Alignment"
mkdir -p "$BLAST_DIR"
log_info "Created BLAST alignment directory: $BLAST_DIR"

# Check for existing BLASTx results
BLASTX_SKIP=false
if [[ "$OVERWRITE" != "true" ]]; then
    shopt -s nullglob
    EXISTING_BLASTX=("$BLAST_DIR"/*/blastx_results/ev_*/*.csv)
    shopt -u nullglob
    if [[ ${#EXISTING_BLASTX[@]} -gt 0 ]]; then
        log_info "BLASTx results exist (${#EXISTING_BLASTX[@]} CSVs). Skipping (OVERWRITE=false)."
        BLASTX_SKIP=true
    else
        log_info "No existing BLASTx results found, proceeding with analysis."
    fi
else
    log_info "OVERWRITE=true, will re-run BLASTx analysis."
fi

if ! $BLASTX_SKIP; then
log_step "BLASTx Ortholog Alignment: $GENE_GROUP"
log_info "Starting BLASTx ortholog alignment for $GENE_GROUP"

# Clean stale BLASTx results to prevent ev_* folder accumulation
for _gdir in "$BLAST_DIR"/*/; do
    [[ -d "${_gdir}blastx_results" ]] && { log_info "Cleaning stale BLASTx results in $(basename "${_gdir%/}")"; rm -rf "${_gdir}blastx_results"; }
done

if ! ensure_query_protein_fastas; then
    log_error "Unable to prepare query_protein_fastas for BLASTx stage."
    exit 1
fi

# ---- Build protein subject database ----
if [[ "$USE_HMMER" == "true" || "$USE_HMMER" == "True" ]]; then
    # Merge all HMMER-identified CD-HIT-reduced proteins into a single FASTA
    HMMER_MERGED="$BLAST_DIR/hmmer_identified_proteins.fa"
    log_info "Merging HMMER-identified proteins from 01_Identification..."

    # Collect matching files safely (glob may not match)
    shopt -s nullglob
    HMMER_FILES=("$BASE_DIR/01_Identification"/*/*/c_CD_HIT_Reduced/*/*_proteins_cdhit.fa)
    shopt -u nullglob

    if [[ ${#HMMER_FILES[@]} -eq 0 ]]; then
        log_error "No HMMER-identified proteins found. Run a_hmmer_identify.sh first."
        exit 1
    fi

    log_info "Found ${#HMMER_FILES[@]} HMMER output files to merge"
    cat "${HMMER_FILES[@]}" > "$HMMER_MERGED"
    MERGED_COUNT=$(grep -c '^>' "$HMMER_MERGED" || echo 0)
    log_info "Merged $MERGED_COUNT proteins from ${#HMMER_FILES[@]} HMMER output files"

    log_info "Building protein BLAST database from merged HMMER results..."
    bash "$MODULES/02_blast_ortholog_alignment/makeblastdb.sh" \
        --input "$HMMER_MERGED" --dbtype prot --outdir "$BLAST_DIR/BLAST_DB" $MAKEDB_FLAGS

    DB_BASE=$(basename "$HMMER_MERGED" .fa)
    DB_PATH="$BLAST_DIR/BLAST_DB/${DB_BASE}_db/$DB_BASE"
    log_info "Protein database created at: $DB_PATH"
else
    # Fallback: build protein DB from each target_proteins entry (parallel)
    log_info "Using fallback method: building protein databases from target_proteins..."
    PROT_DB_PATHS=()
    DB_PIDS=()
    mapfile -t _prot_paths < <(get_toml ortholog_blast target_proteins)
    for prot_path in "${_prot_paths[@]}"; do
        PROT_FULL="$PIPELINE_DIR/$prot_path"
        [[ -f "$PROT_FULL" ]] || { log_warn "Protein file not found: $PROT_FULL"; continue; }
        log_info "Building protein database for: $(basename "$PROT_FULL")"
        wait_for_slot "$MAX_PARALLEL"
        bash "$MODULES/02_blast_ortholog_alignment/makeblastdb.sh" \
            --input "$PROT_FULL" --dbtype prot --outdir "$BLAST_DIR/BLAST_DB" $MAKEDB_FLAGS &
        DB_PIDS+=("$!")
        PROT_BASE=$(basename "$PROT_FULL" .fa)
        PROT_DB_PATHS+=("$BLAST_DIR/BLAST_DB/${PROT_BASE}_db/$PROT_BASE")
    done
    FAILED=0
    for pid in "${DB_PIDS[@]}"; do wait "$pid" || ((++FAILED)); done
    if (( FAILED > 0 )); then
        log_error "$FAILED protein database build(s) failed"
        exit 1
    else
        log_info "Successfully built ${#DB_PIDS[@]} protein databases"
    fi
fi

# Initialize PROT_DB_PATHS based on USE_HMMER setting
    if [[ "$USE_HMMER" == "true" || "$USE_HMMER" == "True" ]]; then
        # When using HMMER, we have a single merged database
        PROT_DB_PATHS=("$DB_PATH")
        log_info "Using HMMER-merged protein database: $(basename "$DB_PATH")"
    fi
    
# Run BLASTx queries in parallel (MAX_PARALLEL jobs × THREADS_PER_JOB threads)
    log_info "Running BLASTx queries against ${#PROT_DB_PATHS[@]} protein database(s)..."
    BLAST_PIDS=()
    mapfile -t CONFIGURED_BLASTX_QUERIES < <(read_toml_list ortholog_blast query_fastas)
    RESOLVED_BLASTX_QUERIES=()
    MISSING_BLASTX_QUERIES=0
    for query in "${CONFIGURED_BLASTX_QUERIES[@]}"; do
        if QUERY_FULL=$(resolve_configured_path "$query"); then
            RESOLVED_BLASTX_QUERIES+=("$QUERY_FULL")
        else
            ((++MISSING_BLASTX_QUERIES))
        fi
    done

    if (( MISSING_BLASTX_QUERIES > 0 && ${#RESOLVED_BLASTX_QUERIES[@]} > 0 )); then
        log_warn "Skipped $MISSING_BLASTX_QUERIES configured BLASTx query FASTA(s) that were not found; check config/${GENE_GROUP}/02_blast_alignment.toml"
    fi

    if [[ ${#RESOLVED_BLASTX_QUERIES[@]} -eq 0 ]]; then
        log_warn "No BLASTx query FASTAs are available for $GENE_GROUP; skipping BLASTx ortholog alignment."
        BLAST_PIDS=()
    fi

    for DB_PATH in "${PROT_DB_PATHS[@]}"; do
        for QUERY_FULL in "${RESOLVED_BLASTX_QUERIES[@]}"; do
            
            # Determine genome-specific output directory
            DB_DIR=$(dirname "$DB_PATH")
            DB_NAME=$(basename "$DB_DIR")
            TOP_FOLDER=$(resolve_top_folder "$DB_NAME")
            GENOME_BLAST_DIR="$BLAST_DIR/$TOP_FOLDER"
            mkdir -p "$GENOME_BLAST_DIR"
            
            log_info "Running BLASTx: query=$(basename "$QUERY_FULL") vs db=$(basename "$DB_PATH")"
            wait_for_slot "$MAX_PARALLEL"
            bash "$MODULES/02_blast_ortholog_alignment/blastx.sh" \
                --query "$QUERY_FULL" \
                --db "$DB_PATH" \
                --evalue "$BLASTX_EVALUE" \
                --word-size "$BLASTX_WORD_SIZE" \
                --max-target-seqs "$BLASTX_MAX_TARGET_SEQS" \
                --max-hsps "$BLASTX_MAX_HSPS" \
                --qcov-hsp-perc "$BLASTX_QCOV_HSP_PERC" \
                --matrix "$BLASTX_MATRIX" \
                --seg "$BLASTX_SEG" \
                --outdir "$GENOME_BLAST_DIR" \
                --threads "$THREADS_PER_JOB" &
            BLAST_PIDS+=("$!")
        done
    done

    # Wait and check for failures
    log_info "Waiting for ${#BLAST_PIDS[@]} BLASTx jobs to complete..."
    FAILED=0
    for pid in "${BLAST_PIDS[@]}"; do
        wait "$pid" || ((++FAILED))
    done
    if (( FAILED > 0 )); then
        log_error "$FAILED BLASTx job(s) failed"
        exit 1
    else
        log_info "All BLASTx jobs completed successfully"

        # Create curated results per genome directory (3 versions: hmmer_only, combined, plant_only)
        for genome_dir in "$BLAST_DIR"/*/; do
            genome_dir="${genome_dir%/}"
            if [[ -d "$genome_dir" && "$genome_dir" != "$BLAST_DIR/BLAST_DB" && "$genome_dir" != "$BLAST_DIR/curated_results" ]]; then
                genome_name=$(basename "$genome_dir")
                mkdir -p "$genome_dir/curated_results"
                # Clean old merged BLASTx CSVs to prevent stale file accumulation
                rm -f "$genome_dir/curated_results"/merged_blastx_*.csv
                if [[ -d "$genome_dir/blastx_results" ]]; then
                    shopt -s nullglob; _csvs=("$genome_dir/blastx_results"/ev_*/*.csv); shopt -u nullglob
                    if [[ ${#_csvs[@]} -eq 0 ]]; then
                        log_warn "No BLASTx CSVs found for genome: $genome_name"
                    else
                        DATE_TAG=$(date +%F)
                        log_info "Merging BLASTx results (3 versions) for genome: $genome_name"
                        bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                            --input-dir "$genome_dir/blastx_results" \
                            --include-pattern "hmmer_identified" \
                            --output "$genome_dir/curated_results/merged_blastx_${DATE_TAG}_${genome_name}_hmmer_only.csv"
                        bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                            --input-dir "$genome_dir/blastx_results" \
                            --output "$genome_dir/curated_results/merged_blastx_${DATE_TAG}_${genome_name}_combined.csv"
                        bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                            --input-dir "$genome_dir/blastx_results" \
                            --exclude-pattern "hmmer_identified" \
                            --output "$genome_dir/curated_results/merged_blastx_${DATE_TAG}_${genome_name}_plant_only.csv"
                    fi
                fi
            fi
        done
    fi  # end if (( FAILED > 0 ))

log_step "BLASTx Ortholog Alignment complete: $GENE_GROUP"
fi  # end if ORTHO_ENABLED
fi  # end BLASTX_SKIP

# ---- BLASTn Ortholog Alignment (nucleotide query vs genome/transcript DB) ----
if [[ "$DO_BLASTN" == "true" ]]; then

# Check for existing BLASTn results
BLASTN_SKIP=false
if [[ "$OVERWRITE" != "true" ]]; then
    shopt -s nullglob
    EXISTING_BLASTN=("$BLAST_DIR"/*/blastn_results/ev_*/*.csv)
    shopt -u nullglob
    if [[ ${#EXISTING_BLASTN[@]} -gt 0 ]]; then
        log_info "BLASTn results exist (${#EXISTING_BLASTN[@]} CSVs). Skipping (OVERWRITE=false)."
        BLASTN_SKIP=true
    else
        log_info "No existing BLASTn results found, proceeding with analysis."
    fi
else
    log_info "OVERWRITE=true, will re-run BLASTn analysis."
fi

if ! $BLASTN_SKIP; then
log_step "BLASTn Ortholog Alignment: $GENE_GROUP"
log_info "Starting BLASTn ortholog alignment for $GENE_GROUP"

# Clean stale BLASTn ortholog results to prevent ev_* folder accumulation
for _gdir in "$BLAST_DIR"/*/; do
    [[ -d "${_gdir}blastn_results" ]] && { log_info "Cleaning stale BLASTn results in $(basename "${_gdir%/}")"; rm -rf "${_gdir}blastn_results"; }
done

# Build all nucl BLAST databases in parallel (genomes + transcripts)
    log_info "Building nucleotide BLAST databases for genomes and transcripts..."
    NUCL_DB_PATHS=()
    DB_PIDS=()

    mapfile -t _orth_genomes < <(get_toml ortholog_blast target_genomes)
    for genome_path in "${_orth_genomes[@]}"; do
        GENOME_FULL="$PIPELINE_DIR/$genome_path"
        [[ -f "${GENOME_FULL}.fa" ]] || { log_warn "Genome not found: ${GENOME_FULL}.fa"; continue; }
        log_info "Building nucleotide database for genome: $(basename "$GENOME_FULL")"
        wait_for_slot "$MAX_PARALLEL"
        bash "$MODULES/02_blast_ortholog_alignment/makeblastdb.sh" \
            --input "${GENOME_FULL}.fa" --dbtype nucl --outdir "$(dirname "$GENOME_FULL")" $MAKEDB_FLAGS &
        DB_PIDS+=("$!")
        GENOME_BASENAME=$(basename "$GENOME_FULL")
        NUCL_DB_PATHS+=("$(dirname "$GENOME_FULL")/${GENOME_BASENAME}_db/$GENOME_BASENAME")
    done

    mapfile -t _orth_transcripts < <(get_toml ortholog_blast target_transcripts)
    for transcript_path in "${_orth_transcripts[@]}"; do
        TRANS_FULL="$PIPELINE_DIR/$transcript_path"
        [[ -f "$TRANS_FULL" ]] || { log_warn "Transcript not found: $TRANS_FULL"; continue; }
        log_info "Building nucleotide database for transcript: $(basename "$TRANS_FULL")"
        wait_for_slot "$MAX_PARALLEL"
        bash "$MODULES/02_blast_ortholog_alignment/makeblastdb.sh" \
            --input "$TRANS_FULL" --dbtype nucl --outdir "$(dirname "$TRANS_FULL")" $MAKEDB_FLAGS &
        DB_PIDS+=("$!")
        TRANS_BASE=$(basename "$TRANS_FULL" .fa)
        NUCL_DB_PATHS+=("$(dirname "$TRANS_FULL")/${TRANS_BASE}_db/$TRANS_BASE")
    done

    FAILED=0
    for pid in "${DB_PIDS[@]}"; do wait "$pid" || ((++FAILED)); done
    if (( FAILED > 0 )); then
        log_error "$FAILED nucleotide database build(s) failed"
        exit 1
    else
        log_info "Successfully built ${#DB_PIDS[@]} nucleotide databases"
    fi
if [[ ${#NUCL_DB_PATHS[@]} -eq 0 ]]; then
    log_warn "No nucleotide databases built. Skipping BLASTn ortholog alignment."
else
    # Collect BLASTn query files: HMMER transcripts + external query_fastas
    log_info "Collecting BLASTn query files..."
    BLASTN_QUERIES=()

    if [[ "$USE_HMMER" == "true" || "$USE_HMMER" == "True" ]]; then
        HMMER_MERGED_TRANS="$BLAST_DIR/hmmer_identified_transcripts.fa"
        shopt -s nullglob
        HMMER_TRANS_FILES=("$BASE_DIR/01_Identification"/*/*/c_CD_HIT_Reduced/*/*_transcripts_cdhit.fa)
        shopt -u nullglob

        if [[ ${#HMMER_TRANS_FILES[@]} -gt 0 ]]; then
            log_info "Found ${#HMMER_TRANS_FILES[@]} HMMER transcript files to merge"
            cat "${HMMER_TRANS_FILES[@]}" > "$HMMER_MERGED_TRANS"
            TRANS_COUNT=$(grep -c '^>' "$HMMER_MERGED_TRANS" || echo 0)
            log_info "Merged $TRANS_COUNT HMMER-identified transcripts for BLASTn queries"
            BLASTN_QUERIES+=("$HMMER_MERGED_TRANS")
        else
            log_warn "No HMMER transcript files found for BLASTn queries"
        fi
    fi

    mapfile -t CONFIGURED_BLASTN_QUERIES < <(read_toml_list ortholog_blast query_fastas)
    MISSING_BLASTN_QUERIES=0
    for query in "${CONFIGURED_BLASTN_QUERIES[@]}"; do
        if QUERY_FULL=$(resolve_configured_path "$query"); then
            BLASTN_QUERIES+=("$QUERY_FULL")
        else
            ((++MISSING_BLASTN_QUERIES))
        fi
    done

    if (( MISSING_BLASTN_QUERIES > 0 && ${#BLASTN_QUERIES[@]} > 0 )); then
        log_warn "Skipped $MISSING_BLASTN_QUERIES configured BLASTn query FASTA(s) that were not found; check config/${GENE_GROUP}/02_blast_alignment.toml"
    fi

    if [[ ${#BLASTN_QUERIES[@]} -eq 0 ]]; then
        log_warn "No BLASTn query FASTAs are available for $GENE_GROUP; skipping BLASTn ortholog alignment."
    else
        # Run BLASTn: all queries × all genome/transcript DBs
        log_info "Running BLASTn queries against ${#NUCL_DB_PATHS[@]} nucleotide database(s)..."
        BLAST_PIDS=()
        for DB_PATH in "${NUCL_DB_PATHS[@]}"; do
            for QUERY_FULL in "${BLASTN_QUERIES[@]}"; do

                # Determine genome-specific output directory
                DB_DIR=$(dirname "$DB_PATH")
                DB_NAME=$(basename "$DB_DIR")
                TOP_FOLDER=$(resolve_top_folder "$DB_NAME")
                GENOME_BLAST_DIR="$BLAST_DIR/$TOP_FOLDER"
                mkdir -p "$GENOME_BLAST_DIR"

                log_info "Running BLASTn: query=$(basename "$QUERY_FULL") vs db=$(basename "$DB_DIR")"
                wait_for_slot "$MAX_PARALLEL"
                bash "$MODULES/02_blast_ortholog_alignment/blastn.sh" \
                    --query "$QUERY_FULL" \
                    --db "$DB_PATH" \
                    --evalue "$BLASTN_EVALUE" \
                    --word-size "$BLASTN_WORD_SIZE" \
                    --max-target-seqs "$BLASTN_MAX_TARGET_SEQS" \
                    --max-hsps "$BLASTN_MAX_HSPS" \
                    --perc-identity "$BLASTN_PERC_IDENTITY" \
                    --qcov-hsp-perc "$BLASTN_QCOV_HSP_PERC" \
                    --task "$BLASTN_TASK" \
                    --dust "$BLASTN_DUST" \
                    --outdir "$GENOME_BLAST_DIR" \
                    --threads "$THREADS_PER_JOB" &
                BLAST_PIDS+=("$!")
            done
        done

        log_info "Waiting for ${#BLAST_PIDS[@]} BLASTn jobs to complete..."
        FAILED=0
        for pid in "${BLAST_PIDS[@]}"; do
            wait "$pid" || ((++FAILED))
        done
        if (( FAILED > 0 )); then
            log_error "$FAILED BLASTn ortholog job(s) failed"
            exit 1
        else
            log_info "All BLASTn jobs completed successfully"
        fi

        # Create curated results per genome directory (3 versions: hmmer_only, combined, plant_only)
        for genome_dir in "$BLAST_DIR"/*/; do
            genome_dir="${genome_dir%/}"
            if [[ -d "$genome_dir" && "$genome_dir" != "$BLAST_DIR/BLAST_DB" && "$genome_dir" != "$BLAST_DIR/curated_results" ]]; then
                genome_name=$(basename "$genome_dir")
                mkdir -p "$genome_dir/curated_results"
                # Clean old merged BLASTn CSVs to prevent stale file accumulation
                rm -f "$genome_dir/curated_results"/merged_blastn_*.csv
                if [[ -d "$genome_dir/blastn_results" ]]; then
                    shopt -s nullglob; _csvs=("$genome_dir/blastn_results"/ev_*/*.csv); shopt -u nullglob
                    if [[ ${#_csvs[@]} -eq 0 ]]; then
                        log_warn "No BLASTn CSVs found for genome: $genome_name"
                    else
                        DATE_TAG=$(date +%F)
                        log_info "Merging BLASTn results (3 versions) for genome: $genome_name"
                        bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                            --input-dir "$genome_dir/blastn_results" \
                            --include-pattern "hmmer_identified" \
                            --output "$genome_dir/curated_results/merged_blastn_${DATE_TAG}_${genome_name}_hmmer_only.csv"
                        bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                            --input-dir "$genome_dir/blastn_results" \
                            --output "$genome_dir/curated_results/merged_blastn_${DATE_TAG}_${genome_name}_combined.csv"
                        bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                            --input-dir "$genome_dir/blastn_results" \
                            --exclude-pattern "hmmer_identified" \
                            --output "$genome_dir/curated_results/merged_blastn_${DATE_TAG}_${genome_name}_plant_only.csv"
                    fi
                fi
            fi
        done
    fi

fi  # end if [[ ${#NUCL_DB_PATHS[@]} -eq 0 ]]
log_step "BLASTn Ortholog Alignment complete: $GENE_GROUP"
fi  # end BLASTN_SKIP
else
    log_info "BLASTn Ortholog Alignment is not enabled. Skipping."
fi  # end DO_BLASTN

# ---- BLASTp Ortholog Alignment (protein query vs protein DB) ----
if [[ "$DO_BLASTP" == "true" ]]; then

mkdir -p "$BLAST_DIR"

BLASTP_SKIP=false
if [[ "$OVERWRITE" != "true" ]]; then
    shopt -s nullglob
    EXISTING_BLASTP=("$BLAST_DIR"/*/blastp_results/ev_*/*.csv)
    shopt -u nullglob
    if [[ ${#EXISTING_BLASTP[@]} -gt 0 ]]; then
        log_info "BLASTp results exist (${#EXISTING_BLASTP[@]} CSVs). Skipping (OVERWRITE=false)."
        BLASTP_SKIP=true
    else
        log_info "No existing BLASTp results found, proceeding with analysis."
    fi
else
    log_info "OVERWRITE=true, will re-run BLASTp analysis."
fi

if ! $BLASTP_SKIP; then
log_step "BLASTp Ortholog Alignment: $GENE_GROUP"
log_info "Starting BLASTp ortholog alignment for $GENE_GROUP"

# Clean stale BLASTp results to prevent ev_* folder accumulation
for _gdir in "$BLAST_DIR"/*/; do
    [[ -d "${_gdir}blastp_results" ]] && { log_info "Cleaning stale BLASTp results in $(basename "${_gdir%/}")"; rm -rf "${_gdir}blastp_results"; }
done

if ! ensure_query_protein_fastas; then
    log_error "Unable to prepare query_protein_fastas for BLASTp stage."
    exit 1
fi

# Reuse protein DBs built by BLASTx if available; otherwise build them now.
if [[ ${#PROT_DB_PATHS[@]} -eq 0 ]]; then
    log_info "Protein databases not yet built (BLASTx was skipped); building for BLASTp..."
    DB_PIDS=()
    if [[ "$USE_HMMER" == "true" || "$USE_HMMER" == "True" ]]; then
        HMMER_MERGED="$BLAST_DIR/hmmer_identified_proteins.fa"
        if [[ ! -f "$HMMER_MERGED" ]]; then
            shopt -s nullglob
            HMMER_FILES=("$BASE_DIR/01_Identification"/*/*/c_CD_HIT_Reduced/*/*_proteins_cdhit.fa)
            shopt -u nullglob
            if [[ ${#HMMER_FILES[@]} -eq 0 ]]; then
                log_error "No HMMER-identified proteins found. Run a_hmmer_identify.sh first."
                exit 1
            fi
            cat "${HMMER_FILES[@]}" > "$HMMER_MERGED"
        fi
        bash "$MODULES/02_blast_ortholog_alignment/makeblastdb.sh" \
            --input "$HMMER_MERGED" --dbtype prot --outdir "$BLAST_DIR/BLAST_DB" $MAKEDB_FLAGS
        DB_BASE=$(basename "$HMMER_MERGED" .fa)
        PROT_DB_PATHS=("$BLAST_DIR/BLAST_DB/${DB_BASE}_db/$DB_BASE")
    else
        mapfile -t _prot_paths < <(get_toml ortholog_blast target_proteins)
        for prot_path in "${_prot_paths[@]}"; do
            PROT_FULL="$PIPELINE_DIR/$prot_path"
            [[ -f "$PROT_FULL" ]] || { log_warn "Protein file not found: $PROT_FULL"; continue; }
            wait_for_slot "$MAX_PARALLEL"
            bash "$MODULES/02_blast_ortholog_alignment/makeblastdb.sh" \
                --input "$PROT_FULL" --dbtype prot --outdir "$BLAST_DIR/BLAST_DB" $MAKEDB_FLAGS &
            DB_PIDS+=("$!")
            PROT_BASE=$(basename "$PROT_FULL" .fa)
            PROT_DB_PATHS+=("$BLAST_DIR/BLAST_DB/${PROT_BASE}_db/$PROT_BASE")
        done
        FAILED=0
        for pid in "${DB_PIDS[@]}"; do wait "$pid" || ((++FAILED)); done
        if (( FAILED > 0 )); then
            log_error "$FAILED protein database build(s) failed"
            exit 1
        fi
    fi
fi

# Run BLASTp: all protein queries × all protein DBs
log_info "Running BLASTp queries against ${#PROT_DB_PATHS[@]} protein database(s)..."
BLAST_PIDS=()
mapfile -t CONFIGURED_BLASTP_QUERIES < <(read_toml_list ortholog_blast query_protein_fastas)
RESOLVED_BLASTP_QUERIES=()
MISSING_BLASTP_QUERIES=0
for query in "${CONFIGURED_BLASTP_QUERIES[@]}"; do
    if QUERY_FULL=$(resolve_configured_path "$query"); then
        RESOLVED_BLASTP_QUERIES+=("$QUERY_FULL")
    else
        ((++MISSING_BLASTP_QUERIES))
    fi
done

if (( MISSING_BLASTP_QUERIES > 0 && ${#RESOLVED_BLASTP_QUERIES[@]} > 0 )); then
    log_warn "Skipped $MISSING_BLASTP_QUERIES configured BLASTp query FASTA(s) not found; check config/${GENE_GROUP}/02_blast_ortholog_alignment.toml"
fi

if [[ ${#RESOLVED_BLASTP_QUERIES[@]} -eq 0 ]]; then
    log_warn "No BLASTp protein query FASTAs available for $GENE_GROUP; skipping BLASTp."
else
    for DB_PATH in "${PROT_DB_PATHS[@]}"; do
        for QUERY_FULL in "${RESOLVED_BLASTP_QUERIES[@]}"; do
            DB_DIR=$(dirname "$DB_PATH")
            DB_NAME=$(basename "$DB_DIR")
            TOP_FOLDER=$(resolve_top_folder "$DB_NAME")
            GENOME_BLAST_DIR="$BLAST_DIR/$TOP_FOLDER"
            mkdir -p "$GENOME_BLAST_DIR"
            log_info "Running BLASTp: query=$(basename "$QUERY_FULL") vs db=$(basename "$DB_PATH")"
            wait_for_slot "$MAX_PARALLEL"
            bash "$MODULES/02_blast_ortholog_alignment/blastp.sh" \
                --query "$QUERY_FULL" \
                --db "$DB_PATH" \
                --evalue "$BLASTP_EVALUE" \
                --word-size "$BLASTP_WORD_SIZE" \
                --max-target-seqs "$BLASTP_MAX_TARGET_SEQS" \
                --max-hsps "$BLASTP_MAX_HSPS" \
                --qcov-hsp-perc "$BLASTP_QCOV_HSP_PERC" \
                --matrix "$BLASTP_MATRIX" \
                --seg "$BLASTP_SEG" \
                --outdir "$GENOME_BLAST_DIR" \
                --threads "$THREADS_PER_JOB" &
            BLAST_PIDS+=("$!")
        done
    done

    log_info "Waiting for ${#BLAST_PIDS[@]} BLASTp jobs to complete..."
    FAILED=0
    for pid in "${BLAST_PIDS[@]}"; do
        wait "$pid" || ((++FAILED))
    done
    if (( FAILED > 0 )); then
        log_error "$FAILED BLASTp job(s) failed"
        exit 1
    else
        log_info "All BLASTp jobs completed successfully"

        # Merge results per genome (3 versions: hmmer_only, combined, plant_only)
        for genome_dir in "$BLAST_DIR"/*/; do
            genome_dir="${genome_dir%/}"
            if [[ -d "$genome_dir" && "$genome_dir" != "$BLAST_DIR/BLAST_DB" && "$genome_dir" != "$BLAST_DIR/curated_results" ]]; then
                genome_name=$(basename "$genome_dir")
                mkdir -p "$genome_dir/curated_results"
                # Clean old merged BLASTp CSVs to prevent stale file accumulation
                rm -f "$genome_dir/curated_results"/merged_blastp_*.csv
                if [[ -d "$genome_dir/blastp_results" ]]; then
                    shopt -s nullglob; _csvs=("$genome_dir/blastp_results"/ev_*/*.csv); shopt -u nullglob
                    if [[ ${#_csvs[@]} -eq 0 ]]; then
                        log_warn "No BLASTp CSVs found for genome: $genome_name"
                    else
                        DATE_TAG=$(date +%F)
                        log_info "Merging BLASTp results (3 versions) for genome: $genome_name"
                        bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                            --input-dir "$genome_dir/blastp_results" \
                            --include-pattern "hmmer_identified" \
                            --output "$genome_dir/curated_results/merged_blastp_${DATE_TAG}_${genome_name}_hmmer_only.csv"
                        bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                            --input-dir "$genome_dir/blastp_results" \
                            --output "$genome_dir/curated_results/merged_blastp_${DATE_TAG}_${genome_name}_combined.csv"
                        bash "$MODULES/02_blast_ortholog_alignment/merge_blast_csv.sh" \
                            --input-dir "$genome_dir/blastp_results" \
                            --exclude-pattern "hmmer_identified" \
                            --output "$genome_dir/curated_results/merged_blastp_${DATE_TAG}_${genome_name}_plant_only.csv"
                    fi
                fi
            fi
        done
    fi
fi

log_step "BLASTp Ortholog Alignment complete: $GENE_GROUP"
fi  # end BLASTP_SKIP
else
    log_info "BLASTp Ortholog Alignment is not enabled. Skipping."
fi  # end DO_BLASTP

# ── BLASTn Visualisation ────────────────────────────────────────────────────
if [[ "$DO_VISUALIZE" == "true" ]]; then
    LOLLIPOP_TOP_N=$(get_toml blast_visualize lollipop_top_n 2>/dev/null || echo "10")
    VIZ_COLORMAP=$(get_toml blast_visualize colormap 2>/dev/null || echo "RdYlGn")
    VIZ_FIGURE_DPI=$(get_toml blast_visualize figure_dpi 2>/dev/null || echo "150")
    VIZ_SAVE_DPI=$(get_toml blast_visualize save_dpi 2>/dev/null || echo "300")
    VIZ_HEATMAP_VMIN=$(get_toml blast_visualize heatmap_vmin 2>/dev/null || echo "65.0")
    VIZ_HEATMAP_VMAX=$(get_toml blast_visualize heatmap_vmax 2>/dev/null || echo "100.0")
    VIZ_HEATMAP_W_SCALE=$(get_toml blast_visualize heatmap_w_scale 2>/dev/null || echo "1.20")
    VIZ_HEATMAP_H_SCALE=$(get_toml blast_visualize heatmap_h_scale 2>/dev/null || echo "0.92")
    VIZ_LOLLIPOP_NCOLS=$(get_toml blast_visualize lollipop_ncols 2>/dev/null || echo "2")
    VIZ_LOLLIPOP_X_PAD=$(get_toml blast_visualize lollipop_x_pad 2>/dev/null || echo "1.60")
    VIZ_LOLLIPOP_DOT_SIZE=$(get_toml blast_visualize lollipop_dot_size 2>/dev/null || echo "100")
    VIZ_LOLLIPOP_DOT_SIZE_HI=$(get_toml blast_visualize lollipop_dot_size_hi 2>/dev/null || echo "150")
    VIZ_HI_STEM_COLOR=$(get_toml blast_visualize hi_stem_color 2>/dev/null || echo "#c7920a")
    VIZ_STEM_COLOR=$(get_toml blast_visualize stem_color 2>/dev/null || echo "#d1d5db")
    VIZ_SCRIPT="$MODULES/02_blast_ortholog_alignment/visualize_blast.sh"

    # Build comma-separated figures list from per-figure operation flags
    VIZ_FIGURES=""
    [[ "$DO_VIZ_HEATMAP" == "true" ]]        && VIZ_FIGURES+="heatmap,"
    [[ "$DO_VIZ_HEATMAP_EVALUE" == "true" ]] && VIZ_FIGURES+="heatmap_evalue,"
    [[ "$DO_VIZ_LOLLIPOP" == "true" ]]       && VIZ_FIGURES+="lollipop,"
    VIZ_FIGURES="${VIZ_FIGURES%,}"   # trim trailing comma

    for genome_dir in "$BLAST_DIR"/*/; do
        if [[ -d "$genome_dir" && \
              "$genome_dir" != "$BLAST_DIR/BLAST_DB/" && \
              "$genome_dir" != "$BLAST_DIR/curated_results/" && \
              "$genome_dir" != "$BLAST_DIR/shared/" ]]; then
            curated="$genome_dir/curated_results"
            if [[ -d "$curated" ]] && compgen -G "$curated/*_plant_only.csv" >/dev/null 2>&1; then
                log_info "  Generating visualisations for: $(basename "$genome_dir")"
                bash "$VIZ_SCRIPT" \
                    --results-dir "$curated" \
                    --gene-group "$GENE_GROUP" \
                    --top-n "$LOLLIPOP_TOP_N" \
                    --figures "$VIZ_FIGURES" \
                    --colormap "$VIZ_COLORMAP" \
                    --figure-dpi "$VIZ_FIGURE_DPI" \
                    --save-dpi "$VIZ_SAVE_DPI" \
                    --heatmap-vmin "$VIZ_HEATMAP_VMIN" \
                    --heatmap-vmax "$VIZ_HEATMAP_VMAX" \
                    --heatmap-w-scale "$VIZ_HEATMAP_W_SCALE" \
                    --heatmap-h-scale "$VIZ_HEATMAP_H_SCALE" \
                    --lollipop-ncols "$VIZ_LOLLIPOP_NCOLS" \
                    --lollipop-x-pad "$VIZ_LOLLIPOP_X_PAD" \
                    --lollipop-dot-size "$VIZ_LOLLIPOP_DOT_SIZE" \
                    --lollipop-dot-size-hi "$VIZ_LOLLIPOP_DOT_SIZE_HI" \
                    --hi-stem-color "$VIZ_HI_STEM_COLOR" \
                    --stem-color "$VIZ_STEM_COLOR"
            fi
        fi
    done
fi

    log_info "Completed processing for gene group: $GENE_GROUP"
done

log_step "All BLAST alignment processes completed"
log_info "Pipeline execution finished successfully"
log_info "Pipeline ended at: $(date)"
wait
teardown_logging
