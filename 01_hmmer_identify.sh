#!/bin/bash
# ============================================================================
# Program 1: HMMER Gene Identification
# ============================================================================
# Comment in/out the gene groups below, then run:
#   bash a_hmmer_identify.sh
# ============================================================================

set -euo pipefail

# ===================== IMPORTANT VARIABLES =====================
# GENE_GROUPS, CPU, MAX_PARALLEL, OVERWRITE, and OPERATIONS are all loaded
# from 01_hmmer_identifyCONFIG.toml [pipeline]  — edit gene_groups there.
# ===============================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"

PROJECT_ROOT="$PIPELINE_DIR"
mkdir -p "$PROJECT_ROOT/logs"

source "$MODULES/logging/logging_utils.sh"

TOML_PARSER="$MODULES/utils/parse_toml.py"
get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }
should_run() { [[ " ${OPERATIONS[*]} " =~ " $1 " ]]; }

# Load GENE_GROUPS from shared config (read before the per-group loop)
SHARED_CONFIG="$PIPELINE_DIR/01_hmmer_identifyCONFIG.toml"
mapfile -t GENE_GROUPS < <(python3 "$TOML_PARSER" "$SHARED_CONFIG" pipeline gene_groups 2>/dev/null)
if [[ ${#GENE_GROUPS[@]} -eq 0 ]]; then
    echo "ERROR: pipeline.gene_groups is empty in 01_hmmer_identifyCONFIG.toml" >&2
    exit 1
fi

TEMP_FILES=()
cleanup_all() {
    rm -f "${TEMP_FILES[@]+"${TEMP_FILES[@]}"}" 2>/dev/null || true
    safe_teardown_logging
}
trap cleanup_all EXIT

for GENE_GROUP in "${GENE_GROUPS[@]}"; do

# Resolve config: deep-merge shared + group configs
CONFIG_DIR="$PIPELINE_DIR/config/${GENE_GROUP}"
MERGE_TOML="$MODULES/utils/merge_toml.py"
if [[ -d "$CONFIG_DIR" ]]; then
    CONFIG_FILE=$(mktemp "${TMPDIR:-/tmp}/${GENE_GROUP}_hmmer_cfg_XXXXXX.toml")
    TEMP_FILES+=("$CONFIG_FILE")
    python3 "$MERGE_TOML" \
        "$PIPELINE_DIR/01_hmmer_identifyCONFIG.toml" \
        "$CONFIG_DIR/00_common.toml" \
        "$CONFIG_DIR/01_hmmer_gene_identification.toml" > "$CONFIG_FILE"
else
    CONFIG_FILE="$PIPELINE_DIR/config/${GENE_GROUP}.toml"
fi

BASE_DIR="$PIPELINE_DIR/$(get_toml general base_dir)"
E_VALUE=$(get_toml identification hmmer_params e_value)
CDHIT_THRESHOLD=$(get_toml identification cdhit_params identity)
DB_TYPE=$(get_toml identification db_type)
MACHINE=$(get_toml pipeline machine)
CPU=$(get_toml pipeline compute "$MACHINE" threads)
MAX_PARALLEL=$(get_toml pipeline compute "$MACHINE" max_parallel)
OVERWRITE=$(get_toml pipeline overwrite)

ops_str=$(get_toml pipeline operations 2>/dev/null || true)
if [[ -n "$ops_str" ]]; then
    mapfile -t OPERATIONS <<< "$ops_str"
else
    OPERATIONS=("BUILD_PROFILES" "SEARCH_GENOMES" "GENERATE_REPORT")
fi

setup_logging
log_step "HMMER Identification: $GENE_GROUP"

IDENT_DIR="$BASE_DIR/01_Identification"
mkdir -p "$IDENT_DIR"

# Discover genes within this gene group and the HMM profiles each gene requires.
# The .genes table maps gene name -> list of profile paths (logical AND across profiles).
# parse_toml.py flattens dict subsections to "GENES__<NAME>=(...)" lines; we extract names.
mapfile -t GENE_NAMES < <(
    python3 "$TOML_PARSER" "$CONFIG_FILE" identification gene_groups "$GENE_GROUP" genes 2>/dev/null \
        | awk -F'=' '/^GENES__/{sub(/^GENES__/, "", $1); print $1}'
)
if [[ ${#GENE_NAMES[@]} -eq 0 ]]; then
    log_error "No genes defined under [identification.gene_groups.$GENE_GROUP.genes]"
    exit 1
fi

# Build the flat (deduped) profile list as the union of every gene's profile entries,
# and remember per-gene profile basenames for the later intersection step.
declare -A _SEEN_PROFILE=()
declare -A GENE_PROFILE_BASES=()   # gene_name -> space-separated profile basenames
PROFILES=()
for _GENE in "${GENE_NAMES[@]}"; do
    mapfile -t _GENE_HMMS < <(get_toml identification gene_groups "$GENE_GROUP" genes "$_GENE")
    if (( ${#_GENE_HMMS[@]} == 0 )); then
        log_error "Gene '$_GENE' has no HMM profiles assigned"
        exit 1
    fi
    _BASES=()
    for _HMM in "${_GENE_HMMS[@]}"; do
        if [[ -z "${_SEEN_PROFILE[$_HMM]:-}" ]]; then
            PROFILES+=("$_HMM")
            _SEEN_PROFILE[$_HMM]=1
        fi
        _BASES+=("$(basename "$_HMM" .hmm)")
    done
    GENE_PROFILE_BASES[$_GENE]="${_BASES[*]}"
done
unset _SEEN_PROFILE _GENE _GENE_HMMS _BASES _HMM
log_info "Gene group $GENE_GROUP: ${#GENE_NAMES[@]} gene(s) [${GENE_NAMES[*]}], ${#PROFILES[@]} unique profile(s)"

ALIGNMENTS=()  # seed alignments deprecated; HMM profiles are pre-built (Pfam)
mapfile -t TARGET_LABELS < <(get_toml identification targets labels)
mapfile -t TARGET_PROTEINS < <(get_toml identification targets proteins)
mapfile -t TARGET_TRANSCRIPTS < <(get_toml identification targets transcripts)
mapfile -t TARGET_CDS < <(get_toml identification targets cds 2>/dev/null || true)
mapfile -t TARGET_ANNOTATIONS < <(get_toml identification targets annotations 2>/dev/null || true)
mapfile -t TARGET_GENOMES < <(get_toml identification targets genomes 2>/dev/null || true)

# Gene-structure / GenBank settings (shared by EXTRACT_GENE_STRUCTURES + GENERATE_GENBANK)
GS_FLANK_BP=$(get_toml identification gene_structures flank_bp 2>/dev/null || echo 1000)
GS_ORGANISM=$(get_toml identification gene_structures organism 2>/dev/null || echo "Solanum melongena")

# --- Phase 1: Build & press HMM profiles (once, before genome loop) ---
if should_run "BUILD_PROFILES"; then
log_step "Building HMM profiles"
for i in "${!PROFILES[@]}"; do
    HMM="$PIPELINE_DIR/${PROFILES[$i]}"
    PROFILE_NAME=$(basename "$HMM" .hmm)

    if [[ $i -lt ${#ALIGNMENTS[@]} ]]; then
        SEED="$PIPELINE_DIR/${ALIGNMENTS[$i]}"
        if [[ -f "$SEED" ]]; then
            log_info "hmmbuild: $PROFILE_NAME"
            hmmbuild --cpu "$CPU" "$HMM" "$SEED"
        fi
    fi
    log_info "hmmpress: $PROFILE_NAME"
    hmmpress -f "$HMM"
done
fi # BUILD_PROFILES

# ---------------------------------------------------------------------------
# resolve_top_folder <genome_label>
#   Route output to a genome-specific first-level folder when possible.
# ---------------------------------------------------------------------------
resolve_top_folder() {
    local genome_label="$1"
    case "$genome_label" in
        *Solanum_melongena_v4.1*|*smel_v4_1*|*Eggplant*)
            echo "Solanum_melongena_v4.1" ;;
        *GPE001970*|*unito*)
            echo "GPE001970_SMEL5" ;;
        *)
            echo "shared" ;;
    esac
}

# --- Phase 2: Search genomes ---
if should_run "SEARCH_GENOMES"; then
NUM_GENOMES=${#TARGET_PROTEINS[@]}
if (( MAX_PARALLEL > 1 )); then
    # Divide CPU among the number of concurrently running genomes (not total genomes)
    _concurrent=$(( NUM_GENOMES < MAX_PARALLEL ? NUM_GENOMES : MAX_PARALLEL ))
    THREADS_PER_GENOME=$(( CPU / _concurrent ))
    (( THREADS_PER_GENOME < 1 )) && THREADS_PER_GENOME=1
    unset _concurrent
    log_info "Parallel genomes: $NUM_GENOMES  |  Threads per genome: $THREADS_PER_GENOME"
else
    THREADS_PER_GENOME=$CPU
    log_info "Sequential genomes: $NUM_GENOMES  |  Threads per genome: $THREADS_PER_GENOME"
fi

wait_for_slot() { while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do sleep 0.5; done; }
GENOME_PIDS=()
for g in "${!TARGET_PROTEINS[@]}"; do
    wait_for_slot
    (
        PROTEINS="$PIPELINE_DIR/${TARGET_PROTEINS[$g]}"
        TRANSCRIPTS="$PIPELINE_DIR/${TARGET_TRANSCRIPTS[$g]}"
        _CDS_REL="${TARGET_CDS[$g]:-}"
        CDS="$( [[ -n "$_CDS_REL" ]] && echo "$PIPELINE_DIR/$_CDS_REL" || echo "" )"
        GENOME_LABEL="${TARGET_LABELS[$g]}"

        # Determine genome-specific top-level folder.
        # Primary genomes get their own canonical folder directly under IDENT_DIR.
        # All other genomes map to "shared" but must each get a per-genome
        # subdirectory to prevent parallel jobs from racing on the same
        # intermediate files (seed_aln.fa, nucl.hmm, nhmmer_hits.tbl, etc.).
        TOP_FOLDER=$(resolve_top_folder "$GENOME_LABEL")
        if [[ "$TOP_FOLDER" == "shared" ]]; then
            GENOME_DIR="$IDENT_DIR/$TOP_FOLDER/$GENOME_LABEL/e-value_$E_VALUE"
        else
            GENOME_DIR="$IDENT_DIR/$TOP_FOLDER/e-value_$E_VALUE"
        fi

        # Skip if overwrite=false and genome results already exist
        if [[ "$OVERWRITE" == "false" && -d "$GENOME_DIR" && -n "$(ls -A "$GENOME_DIR" 2>/dev/null)" ]]; then
            log_info "Skipping $GENOME_LABEL (results exist, overwrite=false)"
            exit 0
        fi
        mkdir -p "$GENOME_DIR"

        log_step "Target genome: $GENOME_LABEL"

        for i in "${!PROFILES[@]}"; do
            HMM="$PIPELINE_DIR/${PROFILES[$i]}"

            bash "$MODULES/01_hmmer_gene_identification/hmmer.sh" \
                --hmm "$HMM" \
                --proteins "$PROTEINS" \
                --transcripts "$TRANSCRIPTS" \
                --evalue "$E_VALUE" \
                --outdir "$GENOME_DIR" \
                --threads "$THREADS_PER_GENOME" \
                --cdhit-threshold "$CDHIT_THRESHOLD" \
                --search-mode "$DB_TYPE" \
                --config "$CONFIG_FILE" \
                --skip-build
        done

        # ---- Per-gene domain intersection ------------------------------------
        # Each gene's hit set = intersection of hit IDs from its required profiles
        # (logical AND across the profiles listed under [.genes.<gene>]).
        # Single-profile genes pass through unchanged.
        GENES_DIR="$GENOME_DIR/d_GENES"
        mkdir -p "$GENES_DIR"
        for GENE in "${GENE_NAMES[@]}"; do
            read -r -a SG_BASES <<< "${GENE_PROFILE_BASES[$GENE]}"
            SG_OUT="$GENES_DIR/$GENE"
            mkdir -p "$SG_OUT"

            HIT_FILES=()
            MISSING=false
            for PB in "${SG_BASES[@]}"; do
                HF="$GENOME_DIR/a_HMMER_RESULTS/$PB/${PB}_hit_ids.txt"
                if [[ ! -s "$HF" ]]; then
                    log_warn "[$GENOME_LABEL] gene $GENE: profile $PB has no hits — gene set empty"
                    MISSING=true
                    break
                fi
                HIT_FILES+=("$HF")
            done

            SG_IDS="$SG_OUT/${GENE}_hit_ids.txt"
            if [[ "$MISSING" == "true" ]]; then
                : > "$SG_IDS"
                continue
            fi

            # Intersect: each hit-IDs file is already sort -u; an ID present in all
            # N files appears N times across cat input -> awk count == N.
            n_files=${#HIT_FILES[@]}
            awk -v n="$n_files" '{c[$1]++} END{for(k in c) if(c[k]==n) print k}' \
                "${HIT_FILES[@]}" | sort -u > "$SG_IDS"
            HIT_N=$(wc -l < "$SG_IDS")
            log_info "[$GENOME_LABEL] gene $GENE: $HIT_N hit(s) from intersection of ${n_files} profile(s) [${SG_BASES[*]}]"
            (( HIT_N == 0 )) && continue

            # Extract proteins and transcripts for the intersected ID set.
            awk 'NR==FNR{ids[$1]=1; next}
                 /^>/{ if(p) print ""; p=0; split($0,a," "); s=substr(a[1],2);
                       if(s in ids) p=1; if(p) print; next }
                 p{ gsub(/\.$/,""); print } END{ if(p) print "" }' \
                "$SG_IDS" "$PROTEINS" > "$SG_OUT/${GENE}_proteins.fa"

            if [[ -f "$TRANSCRIPTS" ]]; then
                awk 'NR==FNR{ids[$1]=1; next}
                     /^>/{ if(p) print ""; p=0; split($0,a," "); s=substr(a[1],2);
                           if(s in ids) p=1; if(p) print; next }
                     p{ gsub(/\.$/,""); print } END{ if(p) print "" }' \
                    "$SG_IDS" "$TRANSCRIPTS" > "$SG_OUT/${GENE}_transcripts.fa"
            fi

            if [[ -n "$CDS" && -f "$CDS" ]]; then
                awk 'NR==FNR{ids[$1]=1; next}
                     /^>/{ if(p) print ""; p=0; split($0,a," "); s=substr(a[1],2);
                           if(s in ids) p=1; if(p) print; next }
                     p{ gsub(/\.$/,""); print } END{ if(p) print "" }' \
                    "$SG_IDS" "$CDS" > "$SG_OUT/${GENE}_cds.fa"
            fi

            # CD-HIT cluster the gene-level proteins/transcripts.
            # Guard with || log_warn: a cd-hit failure on a tiny intersected set
            # should not kill the genome subshell (set -e propagates here).
            if [[ -s "$SG_OUT/${GENE}_proteins.fa" ]]; then
                cd-hit -i "$SG_OUT/${GENE}_proteins.fa" \
                       -o "$SG_OUT/${GENE}_proteins_cdhit.fa" \
                       -c "$CDHIT_THRESHOLD" -n 5 -T "$THREADS_PER_GENOME" \
                       -M 16000 -g 1 -aL 0.8 -aS 0.8 \
                       > "$SG_OUT/cdhit_${GENE}.log" 2>&1 \
                || log_warn "[$GENOME_LABEL] gene $GENE: cd-hit (proteins) failed; see $SG_OUT/cdhit_${GENE}.log"
            fi
            if [[ -s "$SG_OUT/${GENE}_transcripts.fa" ]]; then
                cd-hit-est -i "$SG_OUT/${GENE}_transcripts.fa" \
                           -o "$SG_OUT/${GENE}_transcripts_cdhit.fa" \
                           -c "$CDHIT_THRESHOLD" -n 8 -T "$THREADS_PER_GENOME" \
                           -M 16000 -g 1 -l 30 -r 0 -aL 0.8 -aS 0.8 \
                           > "$SG_OUT/cdhit_est_${GENE}.log" 2>&1 \
                || log_warn "[$GENOME_LABEL] gene $GENE: cd-hit-est (transcripts) failed; see $SG_OUT/cdhit_est_${GENE}.log"
            fi
            if [[ -s "$SG_OUT/${GENE}_cds.fa" ]]; then
                cd-hit-est -i "$SG_OUT/${GENE}_cds.fa" \
                           -o "$SG_OUT/${GENE}_cds_cdhit.fa" \
                           -c "$CDHIT_THRESHOLD" -n 8 -T "$THREADS_PER_GENOME" \
                           -M 16000 -g 1 -l 30 -r 0 -aL 0.8 -aS 0.8 \
                           > "$SG_OUT/cdhit_est_cds_${GENE}.log" 2>&1 \
                || log_warn "[$GENOME_LABEL] gene $GENE: cd-hit-est (CDS) failed; see $SG_OUT/cdhit_est_cds_${GENE}.log"
            fi
        done
    ) &
    GENOME_PIDS+=($!)
    log_info "Launched genome ${TARGET_LABELS[$g]} (PID $!)"
done

# Wait for all parallel jobs and collect failures
FAILED=0
for pid in "${GENOME_PIDS[@]}"; do
    wait "$pid" || FAILED=$((FAILED + 1))
done
if (( FAILED > 0 )); then
    log_error "$FAILED genome job(s) failed"
    exit 1
fi
fi # SEARCH_GENOMES

# --- Phase 3: Extract gene structures from GFF annotations (parallel per genome) ---
if should_run "EXTRACT_GENE_STRUCTURES"; then
log_step "Extracting gene structures for $GENE_GROUP (MAX_PARALLEL=$MAX_PARALLEL)"
STRUCT_SCRIPT="$MODULES/01_hmmer_gene_identification/extract_gene_structures.py"
struct_wait_for_slot() { while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do sleep 0.5; done; }
STRUCT_PIDS=()

for g in "${!TARGET_LABELS[@]}"; do
    GENOME_LABEL="${TARGET_LABELS[$g]}"
    _ANN_REL="${TARGET_ANNOTATIONS[$g]:-}"
    [[ -z "$_ANN_REL" ]] && continue
    GFF3="$PIPELINE_DIR/$_ANN_REL"
    [[ ! -f "$GFF3" ]] && { log_warn "[$GENOME_LABEL] annotation not found: $GFF3 — skipping"; continue; }

    _GEN_REL="${TARGET_GENOMES[$g]:-}"
    GENOME_FASTA=""
    if [[ -n "$_GEN_REL" ]]; then
        _CAND="$PIPELINE_DIR/$_GEN_REL"
        if [[ -f "$_CAND" && -f "$_CAND.fai" ]]; then
            GENOME_FASTA="$_CAND"
        elif [[ -f "$_CAND" ]]; then
            log_warn "[$GENOME_LABEL] genome FASTA $_CAND has no .fai index — structure.fa output disabled"
        fi
    fi

    TOP_FOLDER=$(resolve_top_folder "$GENOME_LABEL")
    if [[ "$TOP_FOLDER" == "shared" ]]; then
        GENOME_DIR="$IDENT_DIR/$TOP_FOLDER/$GENOME_LABEL/e-value_$E_VALUE"
    else
        GENOME_DIR="$IDENT_DIR/$TOP_FOLDER/e-value_$E_VALUE"
    fi
    STRUCT_DIR="$GENOME_DIR/e_GENE_Structures"

    # Collect all (gene, hit_ids) pairs with non-empty hits — parse GFF once per genome.
    GENE_SPEC_ARGS=()
    for GENE in "${GENE_NAMES[@]}"; do
        SG_IDS="$GENOME_DIR/d_GENES/$GENE/${GENE}_hit_ids.txt"
        [[ ! -s "$SG_IDS" ]] && continue
        GENE_SPEC_ARGS+=(--gene-spec "${GENE}:${SG_IDS}")
    done
    (( ${#GENE_SPEC_ARGS[@]} == 0 )) && continue

    struct_wait_for_slot
    (
        python3 "$STRUCT_SCRIPT" \
            --gff3 "$GFF3" \
            --output "$STRUCT_DIR" \
            --overwrite "$OVERWRITE" \
            --genome-fasta "$GENOME_FASTA" \
            --flank-bp "$GS_FLANK_BP" \
            --organism "$GS_ORGANISM" \
            "${GENE_SPEC_ARGS[@]}" \
        && log_info "[$GENOME_LABEL] structures → $STRUCT_DIR/" \
        || log_warn "[$GENOME_LABEL] structure extraction failed (non-fatal)"
    ) &
    STRUCT_PIDS+=($!)
done

# Wait for all background structure jobs
STRUCT_FAILED=0
for pid in "${STRUCT_PIDS[@]}"; do
    wait "$pid" || STRUCT_FAILED=$((STRUCT_FAILED + 1))
done
if (( STRUCT_FAILED > 0 )); then
    log_warn "$STRUCT_FAILED genome(s) failed during gene-structure extraction (non-fatal)"
fi
fi # EXTRACT_GENE_STRUCTURES

# --- Phase 3b: Rebuild GenBank files from existing extraction outputs (parallel) ---
# Runs generate_genbank.py across all genomes; each reads its own e_GENE_Structures/
# tree and writes structure.gb for every {GENE}/{mRNA_ID}/ folder. Useful for
# regenerating .gb without re-parsing the source GFF (EXTRACT_GENE_STRUCTURES
# already emits .gb inline; this operation is a standalone re-run path).
if should_run "GENERATE_GENBANK"; then
log_step "Generating GenBank files for $GENE_GROUP (MAX_PARALLEL=$MAX_PARALLEL)"
GB_SCRIPT="$MODULES/01_hmmer_gene_identification/generate_genbank.py"
gb_wait_for_slot() { while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do sleep 0.5; done; }
GB_PIDS=()

for g in "${!TARGET_LABELS[@]}"; do
    GENOME_LABEL="${TARGET_LABELS[$g]}"
    TOP_FOLDER=$(resolve_top_folder "$GENOME_LABEL")
    if [[ "$TOP_FOLDER" == "shared" ]]; then
        GENOME_DIR="$IDENT_DIR/$TOP_FOLDER/$GENOME_LABEL/e-value_$E_VALUE"
    else
        GENOME_DIR="$IDENT_DIR/$TOP_FOLDER/e-value_$E_VALUE"
    fi
    STRUCT_DIR="$GENOME_DIR/e_GENE_Structures"
    [[ ! -d "$STRUCT_DIR" ]] && continue

    _GB_GEN_REL="${TARGET_GENOMES[$g]:-}"
    GB_GENOME_FASTA=""
    if [[ -n "$_GB_GEN_REL" ]]; then
        _GB_CAND="$PIPELINE_DIR/$_GB_GEN_REL"
        [[ -f "$_GB_CAND" && -f "$_GB_CAND.fai" ]] && GB_GENOME_FASTA="$_GB_CAND"
    fi

    gb_wait_for_slot
    (
        python3 "$GB_SCRIPT" \
            --root "$STRUCT_DIR" \
            --overwrite "$OVERWRITE" \
            --organism "$GS_ORGANISM" \
            --flank-bp "$GS_FLANK_BP" \
            --genome-fasta "$GB_GENOME_FASTA" \
        && log_info "[$GENOME_LABEL] GenBank → $STRUCT_DIR/" \
        || log_warn "[$GENOME_LABEL] GenBank generation failed (non-fatal)"
    ) &
    GB_PIDS+=($!)
done

GB_FAILED=0
for pid in "${GB_PIDS[@]}"; do
    wait "$pid" || GB_FAILED=$((GB_FAILED + 1))
done
if (( GB_FAILED > 0 )); then
    log_warn "$GB_FAILED genome(s) failed during GenBank generation (non-fatal)"
fi
fi # GENERATE_GENBANK

# --- Phase 4: Generate report (tables, Markdown) ---
# Also builds the BLASTn-based ortholog table (section 2 in report.md)
# rooted on Solanum_melongena_V3; requires blastn + makeblastdb on PATH
# (provided by the 'egg' conda env). If BLAST is unavailable the report
# still generates, just without the ortholog section.
if should_run "GENERATE_REPORT"; then
log_step "Generating HMMER report for $GENE_GROUP"
REPORT_SCRIPT="$MODULES/01_hmmer_gene_identification/hmmer_report.py"
if [[ -f "$REPORT_SCRIPT" ]]; then
    python3 "$REPORT_SCRIPT" "$IDENT_DIR" \
        --gene-group "$GENE_GROUP" \
        --evalue "$E_VALUE" \
        --threads "$CPU" \
        && log_info "Report saved to $IDENT_DIR/d_REPORT/" \
        || log_warn "Report generation failed (non-fatal)"
fi
fi # GENERATE_REPORT

log_step "HMMER Identification complete: $GENE_GROUP"
teardown_logging

done

# --- Phase 5: Cross-gene-group summary report ---
REPORT_SCRIPT="$MODULES/01_hmmer_gene_identification/hmmer_report.py"
if (( ${#GENE_GROUPS[@]} >= 1 )) && [[ -f "$REPORT_SCRIPT" ]]; then
    setup_logging
    log_step "Cross-gene-group HMMER summary"
    RESULT_BASE="$PIPELINE_DIR/III_RESULT"
    CROSS_REPORT_DIR="$RESULT_BASE/cross_group_hmmer_report"

    # Collect target labels from last merged config (targets are shared across groups)
    mapfile -t ALL_TARGET_LABELS < <(get_toml identification targets labels 2>/dev/null || true)

    python3 "$REPORT_SCRIPT" --cross-group \
        --result-base "$RESULT_BASE" \
        --gene-groups "${GENE_GROUPS[@]}" \
        --output "$CROSS_REPORT_DIR" \
        --evalue "$E_VALUE" \
        --target-labels "${ALL_TARGET_LABELS[@]}" \
        && log_info "Cross-group report saved to $CROSS_REPORT_DIR/" \
        || log_warn "Cross-group report generation failed (non-fatal)"

    teardown_logging
fi
