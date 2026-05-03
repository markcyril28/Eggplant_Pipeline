#!/bin/bash
# ============================================================================
# Program 12: In Silico PCR (MFEprimer + UCSC isPcr)
# ============================================================================
# Two-engine in silico PCR validation against eggplant reference genomes.
#
#   MFEprimer-3.0  (https://github.com/quwubin/MFEprimer-3.0)
#       - conda-installable (bioconda::mfeprimer)
#       - k-mer indexed, fast, supports degenerate primers
#       - thermodynamics: Tm, dG, dimers, hairpins
#       - JSON + TSV output (machine readable)
#
#   isPcr (UCSC Jim Kent)  (https://hgdownload.soe.ucsc.edu/admin/exe/)
#       - gold-standard, blazing fast
#       - .ooc index (built once per genome) + amplicon FASTA out
#       - no conda; binary downloaded by setup helper
#
# Operations (configured via [in_silico_pcr].operations in TOML):
#   index_genomes  — build MFEprimer index (.ufm) + isPcr .ooc per genome
#   run_engines    — run all engines listed in `engines = [...]` (mfeprimer, ispcr)
#   summarize      — merge per-genome results into one TSV per primer set
#
# Output layout under III_RESULT/{GROUP}/12_In_Silico_PCR/{genome}/:
#   01_Indices/                 (MFEprimer .ufm, isPcr .ooc)
#   02_MFEprimer/{set}.{json,tsv}
#   03_isPcr/{set}.fa           (amplicon FASTA)
#   04_Summary/{set}_amplicons.tsv  (combined per-primer results)
#
# Edit 12_in_silico_pcrCONFIG.toml to choose gene groups, primer sets, and
# operations, then run:
#   bash 12_In_Silico_PCR.sh
# ============================================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"

PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"

source "$MODULES/logging/logging_utils.sh"

TOML_PARSER="$MODULES/utils/parse_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

SHARED_CONFIG="$PIPELINE_DIR/12_in_silico_pcrCONFIG.toml"
mapfile -t GENE_GROUPS < <(python3 "$TOML_PARSER" "$SHARED_CONFIG" pipeline gene_groups 2>/dev/null)
if [[ ${#GENE_GROUPS[@]} -eq 0 ]]; then
    echo "ERROR: pipeline.gene_groups is empty in 12_in_silico_pcrCONFIG.toml" >&2
    exit 1
fi

wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done
}

op_enabled() {
    local op="$1"
    for o in "${OPERATIONS[@]}"; do [[ "$o" == "$op" ]] && return 0; done
    return 1
}

TMP_CONFIG_FILES=()
cleanup_tmp_configs() {
    local cfg
    for cfg in "${TMP_CONFIG_FILES[@]:-}"; do
        [[ -n "$cfg" && -f "$cfg" ]] && rm -f "$cfg"
    done
}
trap 'cleanup_tmp_configs; safe_teardown_logging' EXIT

for GENE_GROUP in "${GENE_GROUPS[@]}"; do

# ── Config resolution ─────────────────────────────────────────────────────
CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
PER_GROUP_PCR="$CONFIG_DIR/12_in_silico_pcr_${GENE_GROUP}.toml"
if [[ -f "$PER_GROUP_PCR" ]]; then
    CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_pcr_cfg_XXXXXX.toml")
    cat "$SHARED_CONFIG" "$PER_GROUP_PCR" > "$CONFIG_FILE"
    TMP_CONFIG_FILES+=("$CONFIG_FILE")
else
    CONFIG_FILE="$SHARED_CONFIG"
fi

# ── Compute profile ───────────────────────────────────────────────────────
MACHINE=$(get_toml pipeline machine 2>/dev/null || echo "Local")
CPU=$(get_toml pipeline compute "$MACHINE" threads 2>/dev/null || nproc)
MAX_PARALLEL=$(get_toml pipeline compute "$MACHINE" max_parallel 2>/dev/null || echo "$CPU")
(( MAX_PARALLEL < 1 )) && MAX_PARALLEL=1
OVERWRITE=$(get_toml pipeline overwrite 2>/dev/null || echo "true")

BASE_DIR="$PIPELINE_DIR/$(get_toml general base_dir 2>/dev/null || echo "III_RESULT/${GENE_GROUP}")"
PCR_DIR="$BASE_DIR/12_In_Silico_PCR"
mkdir -p "$PCR_DIR"

setup_logging

PCR_ENABLED=$(get_toml in_silico_pcr enabled 2>/dev/null || echo "true")
if [[ "$PCR_ENABLED" != "true" && "$PCR_ENABLED" != "True" ]]; then
    log_info "In silico PCR not enabled for $GENE_GROUP. Skipping."
    teardown_logging
    continue
fi

mapfile -t OPERATIONS < <(get_toml in_silico_pcr operations 2>/dev/null \
    || printf '%s\n' "index_genomes" "run_engines" "summarize")

# Engine selection — list of programs commented in/out in the TOML
mapfile -t ENGINES < <(get_toml in_silico_pcr engines 2>/dev/null \
    || printf '%s\n' "mfeprimer" "ispcr")
engine_enabled() {
    local eng="$1"
    for e in "${ENGINES[@]}"; do [[ "$e" == "$eng" ]] && return 0; done
    return 1
}

# Genome list — names parallel to fasta paths
mapfile -t GENOME_NAMES < <(get_toml in_silico_pcr genomes names 2>/dev/null)
mapfile -t GENOME_FASTAS < <(get_toml in_silico_pcr genomes fastas 2>/dev/null)
if [[ ${#GENOME_NAMES[@]} -eq 0 || ${#GENOME_FASTAS[@]} -eq 0 ]]; then
    log_error "in_silico_pcr.genomes.names / fastas missing — nothing to do for $GENE_GROUP"
    teardown_logging
    continue
fi
if [[ ${#GENOME_NAMES[@]} -ne ${#GENOME_FASTAS[@]} ]]; then
    log_error "names / fastas array length mismatch in $GENE_GROUP config"
    teardown_logging
    continue
fi

# Primer sets — each row: "set_name<TAB>tsv_path"
mapfile -t PRIMER_SET_NAMES < <(get_toml in_silico_pcr primers set_names 2>/dev/null)
mapfile -t PRIMER_SET_FILES < <(get_toml in_silico_pcr primers set_files 2>/dev/null)
if [[ ${#PRIMER_SET_NAMES[@]} -ne ${#PRIMER_SET_FILES[@]} ]]; then
    log_error "primers.set_names / set_files length mismatch in $GENE_GROUP config"
    teardown_logging
    continue
fi

ISPCR_BIN=$(get_toml in_silico_pcr ispcr binary 2>/dev/null \
    || echo "$PIPELINE_DIR/modules/12_in_silico_pcr/bin/isPcr")
# Resolve to absolute path if TOML returned a relative value
[[ "$ISPCR_BIN" != /* ]] && ISPCR_BIN="$PIPELINE_DIR/$ISPCR_BIN"

log_step "In Silico PCR Pipeline: $GENE_GROUP"
log_info "Operations: ${OPERATIONS[*]}"
log_info "Engines:    ${ENGINES[*]}"
log_info "Genomes:    ${GENOME_NAMES[*]}"
log_info "Primer sets: ${PRIMER_SET_NAMES[*]:-<none>}"

# ============================================================================
# Operation 1: Build per-genome indices (MFEprimer .ufm, isPcr .ooc)
# ============================================================================
if op_enabled "index_genomes"; then
    log_step "Operation 1/4: Building per-genome indices"

    INDEX_K=$(get_toml in_silico_pcr mfeprimer index_k 2>/dev/null || echo "9")
    OOC_TILE=$(get_toml in_silico_pcr ispcr ooc_tile_size 2>/dev/null || echo "11")
    OOC_REPEAT=$(get_toml in_silico_pcr ispcr ooc_repeat 2>/dev/null || echo "1024")

    for i in "${!GENOME_NAMES[@]}"; do
        gname="${GENOME_NAMES[$i]}"
        gfa="$PIPELINE_DIR/${GENOME_FASTAS[$i]}"
        idx_dir="$PCR_DIR/${gname}/01_Indices"
        mkdir -p "$idx_dir"
        wait_for_slot "$MAX_PARALLEL"
        bash "$MODULES/12_in_silico_pcr/build_indices.sh" \
            --genome      "$gfa" \
            --genome-name "$gname" \
            --outdir      "$idx_dir" \
            --kmer        "$INDEX_K" \
            --ooc-tile    "$OOC_TILE" \
            --ooc-repeat  "$OOC_REPEAT" \
            --ispcr-bin   "$ISPCR_BIN" \
            --overwrite   "$OVERWRITE" &
    done
    wait
fi

# ============================================================================
# Operation 2: MFEprimer search (gated by `engines = [..., "mfeprimer", ...]`)
# ============================================================================
if op_enabled "run_engines" && engine_enabled "mfeprimer"; then
    log_step "Operation 2/4: MFEprimer in silico PCR"

    MFE_MIN_AMPLICON=$(get_toml in_silico_pcr mfeprimer min_amplicon 2>/dev/null || echo "75")
    MFE_MAX_AMPLICON=$(get_toml in_silico_pcr mfeprimer max_amplicon 2>/dev/null || echo "1000")
    MFE_MIS_END=$(get_toml in_silico_pcr mfeprimer misEnd 2>/dev/null || echo "3")
    MFE_INDEX_K=$(get_toml in_silico_pcr mfeprimer index_k 2>/dev/null || echo "9")
    MFE_DIVALENT=$(get_toml in_silico_pcr mfeprimer divalent 2>/dev/null || echo "1.5")
    MFE_MONOVALENT=$(get_toml in_silico_pcr mfeprimer monovalent 2>/dev/null || echo "50")
    MFE_DNTP=$(get_toml in_silico_pcr mfeprimer dntp 2>/dev/null || echo "0.25")
    MFE_OLIGO=$(get_toml in_silico_pcr mfeprimer oligo 2>/dev/null || echo "50")

    for i in "${!GENOME_NAMES[@]}"; do
        gname="${GENOME_NAMES[$i]}"
        idx_dir="$PCR_DIR/${gname}/01_Indices"
        mfe_dir="$PCR_DIR/${gname}/02_MFEprimer"
        mkdir -p "$mfe_dir"
        for j in "${!PRIMER_SET_NAMES[@]}"; do
            set_name="${PRIMER_SET_NAMES[$j]}"
            primers_tsv="$PIPELINE_DIR/${PRIMER_SET_FILES[$j]}"
            [[ -f "$primers_tsv" ]] || { log_warn "primers TSV not found: $primers_tsv"; continue; }
            wait_for_slot "$MAX_PARALLEL"
            bash "$MODULES/12_in_silico_pcr/mfeprimer_run.sh" \
                --primers-tsv "$primers_tsv" \
                --set-name    "$set_name" \
                --index-dir   "$idx_dir" \
                --genome-name "$gname" \
                --outdir      "$mfe_dir" \
                --min-amplicon "$MFE_MIN_AMPLICON" \
                --max-amplicon "$MFE_MAX_AMPLICON" \
                --mis-end      "$MFE_MIS_END" \
                --index-k      "$MFE_INDEX_K" \
                --divalent     "$MFE_DIVALENT" \
                --monovalent   "$MFE_MONOVALENT" \
                --dntp         "$MFE_DNTP" \
                --oligo        "$MFE_OLIGO" \
                --threads      "$CPU" \
                --overwrite    "$OVERWRITE" &
        done
    done
    wait
fi

# ============================================================================
# Operation 3: UCSC isPcr search (gated by `engines = [..., "ispcr", ...]`)
# ============================================================================
if op_enabled "run_engines" && engine_enabled "ispcr"; then
    log_step "Operation 3/4: UCSC isPcr in silico PCR"

    if [[ ! -x "$ISPCR_BIN" ]]; then
        log_warn "isPcr binary not found at $ISPCR_BIN"
        log_warn "  Run: bash $MODULES/12_in_silico_pcr/download_ispcr.sh"
        log_warn "  Skipping ispcr operation."
    else
        ISPCR_MIN=$(get_toml in_silico_pcr ispcr min_amplicon 2>/dev/null || echo "75")
        ISPCR_MAX=$(get_toml in_silico_pcr ispcr max_amplicon 2>/dev/null || echo "4000")
        ISPCR_PERFECT=$(get_toml in_silico_pcr ispcr min_perfect 2>/dev/null || echo "15")
        ISPCR_GOOD=$(get_toml in_silico_pcr ispcr min_good 2>/dev/null || echo "15")
        ISPCR_TILE=$(get_toml in_silico_pcr ispcr tile_size 2>/dev/null || echo "11")
        ISPCR_FLIP=$(get_toml in_silico_pcr ispcr flip_reverse 2>/dev/null || echo "true")

        for i in "${!GENOME_NAMES[@]}"; do
            gname="${GENOME_NAMES[$i]}"
            gfa="$PIPELINE_DIR/${GENOME_FASTAS[$i]}"
            idx_dir="$PCR_DIR/${gname}/01_Indices"
            isp_dir="$PCR_DIR/${gname}/03_isPcr"
            mkdir -p "$isp_dir"
            for j in "${!PRIMER_SET_NAMES[@]}"; do
                set_name="${PRIMER_SET_NAMES[$j]}"
                primers_tsv="$PIPELINE_DIR/${PRIMER_SET_FILES[$j]}"
                [[ -f "$primers_tsv" ]] || { log_warn "primers TSV not found: $primers_tsv"; continue; }
                wait_for_slot "$MAX_PARALLEL"
                bash "$MODULES/12_in_silico_pcr/ispcr_run.sh" \
                    --primers-tsv "$primers_tsv" \
                    --set-name    "$set_name" \
                    --genome      "$gfa" \
                    --genome-name "$gname" \
                    --index-dir   "$idx_dir" \
                    --outdir      "$isp_dir" \
                    --ispcr-bin   "$ISPCR_BIN" \
                    --min-size    "$ISPCR_MIN" \
                    --max-size    "$ISPCR_MAX" \
                    --min-perfect "$ISPCR_PERFECT" \
                    --min-good    "$ISPCR_GOOD" \
                    --tile-size   "$ISPCR_TILE" \
                    --flip-reverse "$ISPCR_FLIP" \
                    --overwrite   "$OVERWRITE" &
            done
        done
        wait
    fi
fi

# ============================================================================
# Operation 4: Summarize amplicons across genomes per primer set
# ============================================================================
if op_enabled "summarize"; then
    log_step "Operation 4/4: Summarizing amplicons per primer set"

    SUM_DIR="$PCR_DIR/04_Summary"
    mkdir -p "$SUM_DIR"

    for j in "${!PRIMER_SET_NAMES[@]}"; do
        set_name="${PRIMER_SET_NAMES[$j]}"
        python3 "$MODULES/12_in_silico_pcr/summarize_amplicons.py" \
            --pcr-dir   "$PCR_DIR" \
            --set-name  "$set_name" \
            --genomes   "${GENOME_NAMES[@]}" \
            --output    "$SUM_DIR/${set_name}_amplicons.tsv" \
            | while IFS= read -r line; do log_info "$line"; done || true
    done
fi

log_step "In Silico PCR Pipeline complete: $GENE_GROUP"
teardown_logging

done
