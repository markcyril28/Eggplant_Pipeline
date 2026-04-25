#!/bin/bash
# ============================================================================
# Module: MEME Suite Promoter Motif Analysis
# ============================================================================
# Three-stage motif analysis pipeline for promoter sequences:
#   1. MEME  — de novo motif discovery in upstream/downstream flanking regions
#   2. TOMTOM — compare discovered motifs to plant TF databases
#   3. FIMO  — scan promoter sequences for known plant TF binding sites
#
# Input:  Per-gene promoter FASTA files from 02_FASTA_with_upstream_and_downstream/
# Output: 04_MEME_Analysis/ per genome label
#
# Standalone (development):
#   bash run_meme_pipeline.sh \
#       --fasta-dir /path/to/02_FASTA_* \
#       --outdir /path/to/04_MEME_Analysis \
#       --databases-dir 2_INPUTS/meme_motif_databases
#
# Orchestrated (production — called from f_motif_analysis.sh):
#   bash run_meme_pipeline.sh \
#       --fasta-dir <dir> --outdir <dir> \
#       --databases-dir <dir> \
#       --label <genome_label> \
#       [--fimo-dbs "JASPAR/...,ARABD/..."] \
#       [--tomtom-dbs "JASPAR/...,ARABD/..."] \
#       [--nmotifs N] [--minw N] [--maxw N] [--mod anr|oops|zoops] \
#       [--threads N] [--overwrite]
#
# Steps controlled by --steps flag (comma-separated):
#   meme, tomtom, fimo   (default: all three)
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../../logging/logging_utils.sh"

# ===================== IMPORTANT VARIABLES =====================
FASTA_DIR=""      # directory of per-gene FASTA files (merged on the fly)
FASTA_FILE=""     # single pre-merged FASTA (skips merge step; takes priority)
OUTDIR=""
DB_DIR=""
LABEL="sequences"
STEPS="meme,tomtom,fimo"
OVERWRITE=false
THREADS=4

# MEME de novo parameters
NMOTIFS=10
MINW=6
MAXW=50
MOD="anr"           # any number of repetitions (best for TF binding)
MARKOV_ORDER=0
TIME_LIMIT=300      # seconds per MEME run

# Plant-relevant databases (relative to DB_DIR, comma-separated)
FIMO_DBS="JASPAR/JASPAR2024_CORE_plants_non-redundant_v2.meme,ARABD/ArabidopsisDAPv1.meme"
TOMTOM_DBS="JASPAR/JASPAR2024_CORE_plants_non-redundant_v2.meme,ARABD/ArabidopsisDAPv1.meme"

# Parallelization: TOMTOM/FIMO are single-threaded; run DB jobs concurrently
MAX_PARALLEL=4      # max concurrent TOMTOM/FIMO database comparison jobs
# ===============================================================

should_run() { [[ ",$STEPS," == *",$1,"* ]]; }

usage() {
    cat <<EOF
Usage: $(basename "$0") --fasta-dir <dir> --outdir <dir> --databases-dir <dir> [options]

Required (one of):
  --fasta-file      Pre-merged multi-sequence FASTA (used directly — skip merge)
  --fasta-dir       Directory of per-gene FASTA files (merged on the fly)
  --outdir          Output directory for MEME analysis results
  --databases-dir   Path to extracted motif databases (2_INPUTS/meme_motif_databases)

Options:
  --label           Genome label (default: sequences)
  --steps           Comma-separated steps: meme,tomtom,fimo (default: all)
  --threads         CPU threads (default: 4)
  --nmotifs         Number of motifs for MEME (default: 10)
  --minw            Minimum motif width (default: 6)
  --maxw            Maximum motif width (default: 50)
  --mod             MEME model: anr|oops|zoops (default: anr)
  --fimo-dbs        FIMO databases, comma-separated relative paths (default: JASPAR2024+ARABD)
  --tomtom-dbs      TOMTOM databases, comma-separated relative paths (default: JASPAR2024+ARABD)
  --max-parallel    Max concurrent TOMTOM/FIMO background jobs (default: 4)
  --time-limit      MEME time limit in seconds per run (default: 300)
  --markov-order    Markov background order for MEME (default: 0)
  --overwrite       Overwrite existing outputs
  -h, --help        Show this help
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta-file)   FASTA_FILE="$2";  shift 2 ;;
        --fasta-dir)    FASTA_DIR="$2";   shift 2 ;;
        --outdir)       OUTDIR="$2";      shift 2 ;;
        --databases-dir) DB_DIR="$2";     shift 2 ;;
        --label)        LABEL="$2";       shift 2 ;;
        --steps)        STEPS="$2";       shift 2 ;;
        --threads)      THREADS="$2";     shift 2 ;;
        --nmotifs)      NMOTIFS="$2";     shift 2 ;;
        --minw)         MINW="$2";        shift 2 ;;
        --maxw)         MAXW="$2";        shift 2 ;;
        --mod)          MOD="$2";         shift 2 ;;
        --fimo-dbs)     FIMO_DBS="$2";      shift 2 ;;
        --tomtom-dbs)   TOMTOM_DBS="$2";    shift 2 ;;
        --max-parallel) MAX_PARALLEL="$2";  shift 2 ;;
        --time-limit)   TIME_LIMIT="$2";    shift 2 ;;
        --markov-order) MARKOV_ORDER="$2";  shift 2 ;;
        --overwrite)    OVERWRITE=true;     shift ;;
        -h|--help)      usage ;;
        *) log_error "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
[[ -z "$FASTA_FILE" && -z "$FASTA_DIR" ]] && { log_error "--fasta-file or --fasta-dir required"; exit 1; }
[[ -z "$OUTDIR" ]]     && { log_error "--outdir required"; exit 1; }
[[ -z "$DB_DIR" ]]     && { log_error "--databases-dir required"; exit 1; }
[[ -n "$FASTA_FILE" && ! -f "$FASTA_FILE" ]] && { log_error "FASTA file not found: $FASTA_FILE"; exit 1; }
[[ -n "$FASTA_DIR"  && ! -d "$FASTA_DIR"  ]] && { log_error "FASTA directory not found: $FASTA_DIR"; exit 1; }
[[ ! -d "$DB_DIR" ]] && { log_error "Motif databases directory not found: $DB_DIR"; exit 1; }

# Check that meme tools are available
for tool in meme fimo tomtom; do
    if ! command -v "$tool" &>/dev/null; then
        log_error "$tool not found in PATH. Install with:"
        log_error "  mamba install -n eggplant -c bioconda meme"
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# Resolve effective alphabet (used by MEME, TOMTOM, and FIMO alphabet checks)
# ---------------------------------------------------------------------------
EFFECTIVE_ALPH="dna"
if [[ "$LABEL" =~ (amino_acid|protein) ]]; then
    EFFECTIVE_ALPH="protein"
fi

# ---------------------------------------------------------------------------
# Output directories
# ---------------------------------------------------------------------------
MERGED_DIR="$OUTDIR/01_merged_promoters"
MEME_DIR="$OUTDIR/02_MEME"
TOMTOM_DIR="$OUTDIR/03_TOMTOM"
FIMO_DIR="$OUTDIR/04_FIMO"

mkdir -p "$MERGED_DIR" "$MEME_DIR" "$TOMTOM_DIR" "$FIMO_DIR"

log_step "MEME Suite Promoter Motif Analysis: $LABEL"
log_info "FASTA input:  ${FASTA_FILE:-$FASTA_DIR}"
log_info "Output:       $OUTDIR"
log_info "Databases:    $DB_DIR"
log_info "Steps:        $STEPS"
log_info "Threads:      $THREADS (MEME) | max parallel (TOMTOM/FIMO): $MAX_PARALLEL"

wait_for_slot() { while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do sleep 0.5; done; }

# ===========================================================================
# Step 1: Merge per-gene FASTA files into a single multi-sequence FASTA
#         (skipped when --fasta-file is given — file is used directly)
# ===========================================================================
MERGED_FASTA="$MERGED_DIR/${LABEL}_promoters.fa"

if [[ -n "$FASTA_FILE" ]]; then
    # Pre-merged input: use as-is; symlink into MERGED_DIR for traceability
    MERGED_FASTA="$FASTA_FILE"
    SEQ_COUNT=$(grep -c '^>' "$MERGED_FASTA" 2>/dev/null || echo 0)
    log_info "Pre-merged FASTA: $MERGED_FASTA ($SEQ_COUNT sequences)"
elif [[ "$OVERWRITE" == true || ! -s "$MERGED_FASTA" ]]; then
    log_step "Merging promoter FASTA files"

    # Collect all .fa / .fasta files in the input directory
    FASTA_COUNT=0
    : > "$MERGED_FASTA"                         # truncate/create
    while IFS= read -r fa; do
        # Sanitise FASTA header: keep only the gene ID (first field before space/colon)
        # Input: >GeneID:chr:start-end:strand region=... upstream=...
        # Output: >GeneID  (sequence lines wrapped at 80 characters)
        awk 'BEGIN { hdr = ""; seq = "" }
/^>/ {
    if (hdr != "") {
        print hdr
        n = length(seq)
        for (i = 1; i <= n; i += 80) print substr(seq, i, 80)
        seq = ""
    }
    split($1, a, ":")
    sub(/^>/, "", a[1])
    hdr = ">" a[1]
    next
}
{ seq = seq $0 }
END {
    if (hdr != "") {
        print hdr
        n = length(seq)
        for (i = 1; i <= n; i += 80) print substr(seq, i, 80)
    }
}' "$fa" >> "$MERGED_FASTA"
        (( ++FASTA_COUNT ))
    done < <(find "$FASTA_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) | sort)

    if [[ "$FASTA_COUNT" -eq 0 ]]; then
        log_error "No FASTA files found in: $FASTA_DIR"
        exit 1
    fi
    log_info "Merged $FASTA_COUNT FASTA files -> $MERGED_FASTA"
else
    log_info "Merged FASTA exists, skipping merge (use --overwrite to redo)"
fi

if [[ -z "$FASTA_FILE" ]]; then
    SEQ_COUNT=$(grep -c '^>' "$MERGED_FASTA" 2>/dev/null || echo 0)
    log_info "Total promoter sequences: $SEQ_COUNT"
fi

if [[ "$SEQ_COUNT" -lt 2 ]]; then
    log_warn "MEME requires at least 2 sequences. Found $SEQ_COUNT — skipping MEME/TOMTOM/FIMO."
    exit 0
fi

# ===========================================================================
# Step 2: MEME — de novo motif discovery
# ===========================================================================
if should_run "meme"; then
    MEME_OUT="$MEME_DIR/$LABEL"

    if [[ "$OVERWRITE" == true || ! -f "$MEME_OUT/meme.txt" ]]; then
        log_step "MEME de novo motif discovery"
        log_info "Model: $MOD  |  nmotifs: $NMOTIFS  |  width: $MINW-$MAXW  |  threads: $THREADS"

        # Alphabet flag (uses EFFECTIVE_ALPH resolved after arg parse)
        meme_alph_flag="-dna"
        [[ "$EFFECTIVE_ALPH" == "protein" ]] && meme_alph_flag="-protein"

        # For protein sequences, strip alignment gap characters ('-') that MEME rejects
        MEME_INPUT="$MERGED_FASTA"
        if [[ "$EFFECTIVE_ALPH" == "protein" ]]; then
            DEGAPPED_FASTA="$MEME_DIR/${LABEL}_degapped.fa"
            sed '/^>/!s/-//g' "$MERGED_FASTA" > "$DEGAPPED_FASTA"
            MEME_INPUT="$DEGAPPED_FASTA"
            log_info "Degapped protein FASTA -> $DEGAPPED_FASTA"
        fi

        # Build MEME command
        meme_cmd=(meme "$MEME_INPUT" \
            "$meme_alph_flag" \
            -oc "$MEME_OUT" \
            -time "$TIME_LIMIT" \
            -mod "$MOD" \
            -nmotifs "$NMOTIFS" \
            -minw "$MINW" \
            -maxw "$MAXW")

        # -objfun and -markov_order require MEME >= 5.x
        _meme_ver=$(meme -version 2>/dev/null || echo "0")
        if [[ "$_meme_ver" == 5.* ]]; then
            meme_cmd+=(-objfun classic -markov_order "$MARKOV_ORDER")
        fi

        # Only add -p if parallel MEME (MPI build) is available
        if meme -p 2 -version &>/dev/null; then
            # Cap threads to available CPU cores to avoid Open MPI slot errors
            _avail_cores=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo "$THREADS")
            if (( THREADS > _avail_cores )); then
                log_warn "Requested $THREADS MEME threads but only $_avail_cores cores available — capping."
                THREADS=$_avail_cores
            fi
            # Allow oversubscription so MEME never fails on slot limits
            # (Open MPI sees physical cores, not hyperthreads, on WSL2)
            export OMPI_MCA_rmaps_base_oversubscribe=1
            # Suppress harmless MPI warnings: InfiniBand BTL + CMA/ptrace on WSL2
            export OMPI_MCA_btl=^openib
            export OMPI_MCA_btl_vader_single_copy_mechanism=none
            meme_cmd+=(-p "$THREADS")
        elif (( THREADS > 1 )); then
            log_warn "Parallel MEME not configured (no MPI build). Running single-threaded."
        fi

        "${meme_cmd[@]}" 2>&1 | tee "$MEME_DIR/${LABEL}_meme.log"
        meme_exit=${PIPESTATUS[0]}

        if [[ $meme_exit -ne 0 ]]; then
            log_error "MEME failed (exit=$meme_exit) — see $MEME_DIR/${LABEL}_meme.log"
            exit $meme_exit
        fi

        log_info "MEME complete -> $MEME_OUT/meme.html"
    else
        log_info "MEME results exist, skipping (use --overwrite to redo)"
    fi
fi

MEME_XML="$MEME_DIR/$LABEL/meme.xml"

# ===========================================================================
# Step 3: TOMTOM — compare discovered motifs to plant TF databases
# ===========================================================================
if should_run "tomtom"; then
    if [[ ! -f "$MEME_XML" ]]; then
        log_warn "MEME XML not found, skipping TOMTOM: $MEME_XML"
    else
        log_step "TOMTOM motif comparison"

        IFS=',' read -ra TOMTOM_DB_LIST <<< "$TOMTOM_DBS"
        tomtom_pids=()
        for db_rel in "${TOMTOM_DB_LIST[@]}"; do
            db_rel="${db_rel// /}"        # strip spaces
            DB_PATH="$DB_DIR/$db_rel"

            if [[ ! -f "$DB_PATH" ]]; then
                log_warn "Database not found, skipping TOMTOM: $DB_PATH"
                continue
            fi

            # Skip database if its alphabet mismatches the query alphabet
            db_alph_line=$(grep -m1 '^ALPHABET=' "$DB_PATH" 2>/dev/null || echo "")
            db_is_dna=false
            [[ "$db_alph_line" == "ALPHABET= ACGT" || "$db_alph_line" == "ALPHABET= ACGTU" ]] && db_is_dna=true
            if [[ "$EFFECTIVE_ALPH" == "protein" && "$db_is_dna" == true ]]; then
                log_warn "Skipping TOMTOM vs $db_rel: DNA database incompatible with protein query"
                continue
            fi
            if [[ "$EFFECTIVE_ALPH" == "dna" && "$db_is_dna" == false && -n "$db_alph_line" ]]; then
                log_warn "Skipping TOMTOM vs $db_rel: protein database incompatible with DNA query"
                continue
            fi

            # Derive a safe directory name from the database basename
            db_name="$(basename "${db_rel%.meme}")"
            TOMTOM_OUT="$TOMTOM_DIR/${LABEL}_vs_${db_name}"

            if [[ "$OVERWRITE" == true || ! -f "$TOMTOM_OUT/tomtom.tsv" ]]; then
                log_info "Comparing vs: $db_rel"
                wait_for_slot
                (
                    tomtom \
                        -oc "$TOMTOM_OUT" \
                        -xalph \
                        -no-ssc \
                        -thresh 0.05 \
                        -min-overlap 5 \
                        "$MEME_XML" "$DB_PATH" \
                        2>&1 | tee "$TOMTOM_DIR/${LABEL}_vs_${db_name}.log"
                    log_info "TOMTOM done -> $TOMTOM_OUT/tomtom.tsv"
                ) &
                tomtom_pids+=("$!:$db_name")
            else
                log_info "TOMTOM results exist for $db_name, skipping"
            fi
        done
        tomtom_failed=0
        for entry in "${tomtom_pids[@]+"${tomtom_pids[@]}"}" ; do
            pid="${entry%%:*}" db="${entry#*:}"
            if ! wait "$pid"; then
                log_error "TOMTOM failed for $db (PID $pid)"
                ((tomtom_failed++)) || true
            fi
        done
        (( tomtom_failed > 0 )) && log_warn "$tomtom_failed TOMTOM comparison(s) failed"
    fi
fi

# ===========================================================================
# Step 4: FIMO — scan sequences for known plant TF binding sites
# ===========================================================================
if should_run "fimo"; then
    log_step "FIMO known-motif scanning"

    IFS=',' read -ra FIMO_DB_LIST <<< "$FIMO_DBS"
    fimo_pids=()
    for db_rel in "${FIMO_DB_LIST[@]}"; do
        db_rel="${db_rel// /}"
        DB_PATH="$DB_DIR/$db_rel"

        if [[ ! -f "$DB_PATH" ]]; then
            log_warn "Database not found, skipping FIMO: $DB_PATH"
            continue
        fi

        # Skip database if its alphabet mismatches the query alphabet
        db_alph_line=$(grep -m1 '^ALPHABET=' "$DB_PATH" 2>/dev/null || echo "")
        db_is_dna=false
        [[ "$db_alph_line" == "ALPHABET= ACGT" || "$db_alph_line" == "ALPHABET= ACGTU" ]] && db_is_dna=true
        if [[ "$EFFECTIVE_ALPH" == "protein" && "$db_is_dna" == true ]]; then
            log_warn "Skipping FIMO vs $db_rel: DNA database incompatible with protein sequences"
            continue
        fi
        if [[ "$EFFECTIVE_ALPH" == "dna" && "$db_is_dna" == false && -n "$db_alph_line" ]]; then
            log_warn "Skipping FIMO vs $db_rel: protein database incompatible with DNA sequences"
            continue
        fi

        db_name="$(basename "${db_rel%.meme}")"
        FIMO_OUT="$FIMO_DIR/${LABEL}_${db_name}"

        if [[ "$OVERWRITE" == true || ! -f "$FIMO_OUT/fimo.tsv" ]]; then
            log_info "Scanning with: $db_rel"
            wait_for_slot
            (
                fimo \
                    --oc "$FIMO_OUT" \
                    --thresh 1e-4 \
                    "$DB_PATH" "$MERGED_FASTA" \
                    2>&1 | tee "$FIMO_DIR/${LABEL}_${db_name}.log"
                log_info "FIMO done -> $FIMO_OUT/fimo.tsv"
            ) &
            fimo_pids+=("$!:$db_name")
        else
            log_info "FIMO results exist for $db_name, skipping"
        fi
    done
    fimo_failed=0
    for entry in "${fimo_pids[@]+"${fimo_pids[@]}"}" ; do
        pid="${entry%%:*}" db="${entry#*:}"
        if ! wait "$pid"; then
            log_error "FIMO failed for $db (PID $pid)"
            ((fimo_failed++)) || true
        fi
    done
    (( fimo_failed > 0 )) && log_warn "$fimo_failed FIMO scan(s) failed"
fi

log_step "MEME Suite analysis complete: $LABEL"
log_info "Results:"
should_run "meme"   && log_info "  MEME:   $MEME_DIR/$LABEL/meme.html"
should_run "tomtom" && log_info "  TOMTOM: $TOMTOM_DIR/"
should_run "fimo"   && log_info "  FIMO:   $FIMO_DIR/"
