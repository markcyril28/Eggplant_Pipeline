#!/bin/bash
# ============================================================================
# Stage 0: Reference Sequence Setup and Download
# ============================================================================
# Bootstrap stage for a fresh machine. Creates the I_RefSeqs/ output tree and
# downloads:
#   - DMP query FASTAs from NCBI (At, Zm, Ip + 12 species from the DMP HIR table)
#   - Solanum melongena Unito Genomics genomes + annotations
#   - Various crop species genomes
#
# All behaviour driven by 0_RefSeq_setup_and_downloadCONFIG.toml
#
# Usage:
#   bash 0_RefSeq_setup_and_download.sh                  # run all enabled ops
#   bash 0_RefSeq_setup_and_download.sh --list           # list operations
#   bash 0_RefSeq_setup_and_download.sh --dry-run        # show planned actions
#   bash 0_RefSeq_setup_and_download.sh --only OP[,OP]   # run subset
#   bash 0_RefSeq_setup_and_download.sh --overwrite      # re-download
# ============================================================================

set -euo pipefail

# ===================== IMPORTANT VARIABLES =====================
# Everything is loaded from 0_RefSeq_setup_and_downloadCONFIG.toml.
# Edit the TOML to change operations, parallelism, species coverage, etc.
# ===============================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"
FASTA_DOWNLOAD_MODULES="$MODULES/00_refseq/fasta_download"
CONFIG_FILE="$PIPELINE_DIR/0_RefSeq_setup_and_downloadCONFIG.toml"
TOML_PARSER="$MODULES/utils/parse_toml.py"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE" >&2
    exit 1
fi

PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"
source "$MODULES/logging/logging_utils.sh"

get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }
should_run() { [[ " ${OPERATIONS[*]} " =~ " $1 " ]]; }

# ----- Load config ----------------------------------------------------------
mapfile -t OPERATIONS    < <(get_toml pipeline operations)
MAX_PARALLEL=$(get_toml pipeline max_parallel)
OVERWRITE=$(get_toml pipeline overwrite)
NCBI_API_KEY_CFG=$(get_toml pipeline ncbi_api_key 2>/dev/null || echo "")

REFSEQ_ROOT_REL=$(get_toml output refseq_root)
SMEL_REL=$(get_toml output smel_subdir)
CROP_REL=$(get_toml output crop_genomes_subdir)
EXTRA_REL=$(get_toml output extra_genomes_subdir)
DMP_REL=$(get_toml output dmp_query_subdir)

REFSEQ_OUT_DIR="$PIPELINE_DIR/$REFSEQ_ROOT_REL"
SMEL_DIR="$REFSEQ_OUT_DIR/$SMEL_REL"
CROP_DIR="$REFSEQ_OUT_DIR/$CROP_REL"
EXTRA_DIR="$REFSEQ_OUT_DIR/$EXTRA_REL"
DMP_QUERY_DIR="$REFSEQ_OUT_DIR/$DMP_REL"

mapfile -t DMP_DIRS      < <(get_toml dmp_queries dirs)
mapfile -t DMP_DOWNLOADS < <(get_toml dmp_queries downloads)
mapfile -t DMP_MERGES    < <(get_toml dmp_queries merge_outputs)

UNITO_REL=$(get_toml smel_unito downloader_path)
UNITO_DEST_REL=$(get_toml smel_unito dest_subdir)
EXTRACT_REL=$(get_toml smel_unito extract_transcripts_path)
CROP_DOWNLOADER_REL=$(get_toml crop_genomes downloader_path)

# Export NCBI API key so sourced library scripts pick it up
if [[ -n "$NCBI_API_KEY_CFG" ]]; then
    export NCBI_API_KEY="$NCBI_API_KEY_CFG"
fi
export OVERWRITE

DRY_RUN=false
ONLY_OPS=""

# ----- CLI overrides --------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --list)
            echo "Available operations (from $CONFIG_FILE):"
            for op in "${OPERATIONS[@]}"; do echo "  - $op"; done
            exit 0
            ;;
        --dry-run|-n)  DRY_RUN=true; shift ;;
        --only)        ONLY_OPS="$2"; shift 2 ;;
        --overwrite)   OVERWRITE=true; export OVERWRITE; shift ;;
        --parallel)    MAX_PARALLEL="$2"; shift 2 ;;
        -h|--help)     sed -n '2,21p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -n "$ONLY_OPS" ]]; then
    IFS=',' read -r -a OPERATIONS <<< "$ONLY_OPS"
fi

# ----- Helpers --------------------------------------------------------------
wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done
}

# Skip vs. overwrite policy
# OVERWRITE=true  -> always re-run (re-download / re-merge / re-extract)
# OVERWRITE=false -> skip when expected outputs already exist
# Returns 0 (skip) or 1 (proceed). Args: <work_dir> [pattern]
should_skip_dir() {
    local dir="$1" pattern="${2:-*.fasta}"
    [[ "$OVERWRITE" == "true" ]] && return 1
    [[ -d "$dir" ]] || return 1
    # Find any non-empty file matching pattern
    local existing
    existing=$(find "$dir" -maxdepth 1 -type f -name "$pattern" -size +0c 2>/dev/null | head -1)
    [[ -n "$existing" ]]
}

should_skip_file() {
    local path="$1"
    [[ "$OVERWRITE" == "true" ]] && return 1
    [[ -s "$path" ]]
}

run_module() {
    local module="$1" work_dir="$2"
    if [[ ! -f "$module" ]]; then
        log_warn "  [SKIP-MISSING] Module not found: $(basename "$module")"
        return 0
    fi
    if should_skip_dir "$work_dir" "*.fasta"; then
        log_info "  [SKIP-EXISTS] $(basename "$work_dir")/ already has FASTA(s) — set overwrite=true to redownload"
        return 0
    fi
    mkdir -p "$work_dir"
    if $DRY_RUN; then
        echo "[DRY-RUN] (cd $work_dir && bash $module)"
        return 0
    fi
    ( cd "$work_dir" && bash "$module" )
}

run_merge() {
    local out_filename="$1" work_dir="$2"
    local merger="$FASTA_DOWNLOAD_MODULES/merge_fasta_pwd.sh"
    if [[ ! -d "$work_dir" ]]; then
        log_warn "  [SKIP-MISSING] merge target dir absent: $(basename "$work_dir")"
        return 0
    fi
    if should_skip_file "$work_dir/$out_filename"; then
        log_info "  [SKIP-EXISTS] $(basename "$work_dir")/$out_filename — set overwrite=true to regenerate"
        return 0
    fi
    if $DRY_RUN; then
        echo "[DRY-RUN] (cd $work_dir && bash $merger $out_filename)"
        return 0
    fi
    # Reuse legacy merge scripts when present (Arabidopsis/Zm/Ip already have them)
    local legacy_lookup
    case "$(basename "$work_dir")" in
        Arabidopsis_thaliana)  legacy_lookup="$FASTA_DOWNLOAD_MODULES/merged_fasta_At_v1.sh" ;;
        Zea_mays)              legacy_lookup="$FASTA_DOWNLOAD_MODULES/merged_ZmDMPs_fasta.sh" ;;
        Ipomoea_batatas)       legacy_lookup="$FASTA_DOWNLOAD_MODULES/merged_IpDMPs_fasta.sh" ;;
        Other_DMPs)            legacy_lookup="$FASTA_DOWNLOAD_MODULES/merged_OtherDMPs_fasta.sh" ;;
        *)                     legacy_lookup="" ;;
    esac
    if [[ -n "$legacy_lookup" && -f "$legacy_lookup" ]]; then
        ( cd "$work_dir" && bash "$legacy_lookup" )
    else
        ( cd "$work_dir" && bash "$merger" "$out_filename" )
    fi
}

setup_logging
log_step "Stage 0: Reference Sequence Setup and Download"
log_info "PIPELINE_DIR : $PIPELINE_DIR"
log_info "CONFIG_FILE  : $CONFIG_FILE"
log_info "REFSEQ_OUT   : $REFSEQ_OUT_DIR"
log_info "OPERATIONS   : ${OPERATIONS[*]}"
log_info "MAX_PARALLEL : $MAX_PARALLEL"
log_info "OVERWRITE    : $OVERWRITE"
log_info "DRY_RUN      : $DRY_RUN"
log_info "DMP species  : ${#DMP_DIRS[@]}"

trap 'teardown_logging 2>/dev/null || true' EXIT

# ============================================================================
# OP: SETUP_DIRS
# ============================================================================
if should_run "SETUP_DIRS"; then
    log_step "[SETUP_DIRS] Creating I_RefSeqs/ directory tree"
    target_dirs=(
        "$REFSEQ_OUT_DIR"
        "$SMEL_DIR"
        "$SMEL_DIR/Solanum_melongena_v4.1"
        "$SMEL_DIR/Solanum_melongena_consortium"
        "$SMEL_DIR/ncbi"
        "$SMEL_DIR/unito_genomics"
        "$SMEL_DIR/unito_genomics/unito_genomics_data"
        "$CROP_DIR"
        "$EXTRA_DIR"
        "$DMP_QUERY_DIR"
    )
    for d in "${DMP_DIRS[@]}"; do target_dirs+=("$DMP_QUERY_DIR/$d"); done
    if $DRY_RUN; then
        printf '[DRY-RUN] mkdir -p %s\n' "${target_dirs[@]}"
    else
        mkdir -p "${target_dirs[@]}"
    fi
fi

# ============================================================================
# OP: DOWNLOAD_DMP_QUERIES
# ============================================================================
if should_run "DOWNLOAD_DMP_QUERIES"; then
    log_step "[DOWNLOAD_DMP_QUERIES] Downloading DMP query FASTAs from NCBI"
    for ((i=0; i<${#DMP_DIRS[@]}; i++)); do
        species_dir="${DMP_DIRS[i]}"
        download_module="${DMP_DOWNLOADS[i]}"
        [[ -z "$download_module" ]] && { log_info "  -> $species_dir: no download module (skip)"; continue; }
        out_dir="$DMP_QUERY_DIR/$species_dir"
        module_path="$FASTA_DOWNLOAD_MODULES/$download_module"
        wait_for_slot "$MAX_PARALLEL"
        log_info "  -> $species_dir via $download_module"
        run_module "$module_path" "$out_dir" &
    done
    wait
    log_info "[DOWNLOAD_DMP_QUERIES] All downloads complete"
fi

# ============================================================================
# OP: MERGE_DMP_QUERIES
# ============================================================================
if should_run "MERGE_DMP_QUERIES"; then
    log_step "[MERGE_DMP_QUERIES] Merging per-species FASTAs"
    for ((i=0; i<${#DMP_DIRS[@]}; i++)); do
        species_dir="${DMP_DIRS[i]}"
        merge_out="${DMP_MERGES[i]}"
        [[ -z "$merge_out" ]] && continue
        out_dir="$DMP_QUERY_DIR/$species_dir"
        log_info "  -> $species_dir -> $merge_out"
        run_merge "$merge_out" "$out_dir"
    done
fi

# ============================================================================
# OP: DOWNLOAD_SMEL_UNITO
# ============================================================================
if should_run "DOWNLOAD_SMEL_UNITO"; then
    log_step "[DOWNLOAD_SMEL_UNITO] Mirroring Solanum melongena Unito Genomics"
    UNITO_MODULE="$REFSEQ_OUT_DIR/$UNITO_REL"
    UNITO_DEST="$REFSEQ_OUT_DIR/$UNITO_DEST_REL"
    # The downloader uses wget --continue --timestamping which is per-file safe.
    # At the orchestrator level we still skip the whole stage when target subdirs
    # already contain genome FASTAs and OVERWRITE=false.
    if should_skip_dir "$UNITO_DEST/genomes" "*.fa.gz" \
       || should_skip_dir "$UNITO_DEST/genomes" "*.fa"; then
        log_info "  [SKIP-EXISTS] $UNITO_DEST/genomes already populated — set overwrite=true to refresh"
    elif [[ -f "$UNITO_MODULE" ]]; then
        if $DRY_RUN; then
            echo "[DRY-RUN] bash $UNITO_MODULE $UNITO_DEST"
        else
            mkdir -p "$UNITO_DEST"
            bash "$UNITO_MODULE" "$UNITO_DEST"
        fi
    else
        log_warn "  [SKIP-MISSING] Unito Genomics downloader not found at $UNITO_MODULE"
    fi
fi

# ============================================================================
# OP: EXTRACT_TRANSCRIPTS
# ============================================================================
if should_run "EXTRACT_TRANSCRIPTS"; then
    log_step "[EXTRACT_TRANSCRIPTS] Extracting transcripts/CDS/proteins via gffread"
    EXTRACT_MODULE="$REFSEQ_OUT_DIR/$EXTRACT_REL"
    # extract_transcripts.sh has its own per-genome OVERWRITE check; we export
    # OVERWRITE so it honors the same toggle, and short-circuit the whole stage
    # if the transcripts/ output dir is already populated.
    EXTRACT_OUT_DIR="$(dirname "$EXTRACT_MODULE")/transcripts"
    if should_skip_dir "$EXTRACT_OUT_DIR" "*_transcripts.fa"; then
        log_info "  [SKIP-EXISTS] transcripts/ already populated — set overwrite=true to re-extract"
    elif [[ -f "$EXTRACT_MODULE" ]]; then
        if $DRY_RUN; then
            echo "[DRY-RUN] OVERWRITE=$OVERWRITE bash $EXTRACT_MODULE --parallel $MAX_PARALLEL"
        else
            OVERWRITE="$OVERWRITE" bash "$EXTRACT_MODULE" --parallel "$MAX_PARALLEL"
        fi
    else
        log_warn "  [SKIP-MISSING] extract_transcripts.sh not found at $EXTRACT_MODULE"
    fi
fi

# ============================================================================
# OP: DOWNLOAD_CROP_GENOMES
# ============================================================================
if should_run "DOWNLOAD_CROP_GENOMES"; then
    log_step "[DOWNLOAD_CROP_GENOMES] Downloading crop species genomes"
    CROP_MODULE="$REFSEQ_OUT_DIR/$CROP_DOWNLOADER_REL"
    # RefSeq_Downloader.sh already does per-file existence checks. We only
    # short-circuit the entire stage when the crop output dir is fully empty
    # of subdirs is unsafe (each species runs independently). If overwrite=true,
    # we delete all *.fa/*.fa.gz in the crop dir so the downloader re-fetches.
    if [[ -f "$CROP_MODULE" ]]; then
        if [[ "$OVERWRITE" == "true" ]]; then
        if $DRY_RUN; then
                echo "[DRY-RUN] find $CROP_DIR -type f \( -name '*.fa' -o -name '*.fa.gz' \) -delete"
            echo "[DRY-RUN] bash $CROP_MODULE"
        else
                log_info "  overwrite=true: clearing existing crop FASTA files"
                find "$CROP_DIR" -type f \( -name '*.fa' -o -name '*.fa.gz' -o -name '*.fna' \
                    -o -name '*.fna.gz' -o -name '*.faa' -o -name '*.faa.gz' \) -delete 2>/dev/null || true
                bash "$CROP_MODULE"
            fi
        else
            if $DRY_RUN; then
                echo "[DRY-RUN] bash $CROP_MODULE  # per-file skip honored by module"
            else
                log_info "  overwrite=false: per-file skip is honored by RefSeq_Downloader.sh"
            bash "$CROP_MODULE"
        fi
        fi
    else
        log_warn "  [SKIP-MISSING] Crop genomes downloader not found at $CROP_MODULE"
    fi
fi

log_step "Stage 0 complete"
log_info "Output tree at: $REFSEQ_OUT_DIR"
