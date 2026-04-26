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
NMOTIFS=15
MINW=6
MAXW=50
MOD="zoops"         # zero or one occurrence per sequence (best for TF binding)
MARKOV_ORDER=3
TIME_LIMIT=18000    # seconds per MEME run (5 h)

# Plant-relevant databases (relative to DB_DIR, comma-separated)
FIMO_DBS="JASPAR/JASPAR2024_CORE_plants_non-redundant_v2.meme,ARABD/ArabidopsisDAPv1.meme"
TOMTOM_DBS="JASPAR/JASPAR2024_CORE_plants_non-redundant_v2.meme,ARABD/ArabidopsisDAPv1.meme"

# Parallelization: TOMTOM/FIMO are single-threaded; run DB jobs concurrently
MAX_PARALLEL=4      # max concurrent TOMTOM/FIMO database comparison jobs

# Sequence alphabet: dna | protein | auto (auto-detects from LABEL name)
ALPHABET="auto"

# JPEG export: assemble all motif logos into a single grid image after MEME
# color_scheme: preset name (default|nature|colorblind|dmp) or path to .alph file.
# Colors are baked into MEME XML at discovery time via -alph flag.
COLOR_SCHEME="default"
JPEG_DPI=300
JPEG_COLUMNS=5
JPEG_QUALITY=92
JPEG_LOGO_PADDING=10

# Motif location diagram palette: version name from meme_motif_colors.toml
# (wong|warm|cool|pastel|high_contrast|dark|dmp_interaction) or comma-separated hex codes.
MOTIF_LOCATION_PALETTE="wong"
# Background color for motif locations diagram: dark|light or hex (empty = white)
MOTIF_LOCATION_BG=""
# Phylo-ordered alignment file: CLUSTAL format.  When set, sequences in the
# motif location diagram are sorted top-to-bottom by the order they appear
# in this file (typically a phylogenetically sorted alignment).
PHYLO_ORDER_FILE=""
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
  --threads         CPU threads for MEME (default: 4)
  --nmotifs         Number of motifs for MEME (default: 15)
  --minw            Minimum motif width (default: 6)
  --maxw            Maximum motif width (default: 50)
  --mod             MEME model: anr|oops|zoops (default: zoops)
  --fimo-dbs        FIMO databases, comma-separated relative paths (default: JASPAR2024+ARABD)
  --tomtom-dbs      TOMTOM databases, comma-separated relative paths (default: JASPAR2024+ARABD)
  --max-parallel    Max concurrent TOMTOM/FIMO background jobs (default: 4)
  --time-limit      MEME time limit in seconds per run (default: 18000)
  --markov-order    Markov background order for MEME (default: 3)
  --color-scheme      Color palette: default|nature|colorblind|dmp or path to .alph file (default: default)
  --jpeg-dpi          DPI for JPEG output (default: 300)
  --jpeg-columns      Motif logos per row in JPEG grid (default: 5)
  --jpeg-quality      JPEG compression quality 0–100 (default: 92)
  --jpeg-logo-padding Pixel gap between logo tiles in montage grid (default: 10)
  --motif-palette       Motif location palette: version name (wong|warm|cool|pastel|high_contrast|dark|dmp_interaction)
                        or comma-separated hex codes (default: wong)
  --motif-bg            Background for motif locations plot: dark|light or hex code (default: empty = white)
  --phylo-order-file    CLUSTAL alignment file whose sequence order defines row order (default: empty)
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
        --alphabet)       ALPHABET="$2";       shift 2 ;;
        --color-scheme)      COLOR_SCHEME="$2";      shift 2 ;;
        --jpeg-dpi)          JPEG_DPI="$2";          shift 2 ;;
        --jpeg-columns)      JPEG_COLUMNS="$2";      shift 2 ;;
        --jpeg-quality)      JPEG_QUALITY="$2";      shift 2 ;;
        --jpeg-logo-padding) JPEG_LOGO_PADDING="$2"; shift 2 ;;
        --motif-palette)     MOTIF_LOCATION_PALETTE="$2"; shift 2 ;;
        --motif-bg)          MOTIF_LOCATION_BG="$2";      shift 2 ;;
        --phylo-order-file)  PHYLO_ORDER_FILE="$2";       shift 2 ;;
        --overwrite)      OVERWRITE=true;       shift ;;
        -h|--help)        usage ;;
        *) log_error "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Resolve effective alphabet (used by MEME, TOMTOM, and FIMO alphabet checks)
# ---------------------------------------------------------------------------
EFFECTIVE_ALPH="dna"
if [[ "$ALPHABET" == "auto" ]]; then
    if [[ "$LABEL" =~ (amino_acid|protein) ]]; then
        EFFECTIVE_ALPH="protein"
    fi
elif [[ "$ALPHABET" == "protein" ]]; then
    EFFECTIVE_ALPH="protein"
fi

# ---------------------------------------------------------------------------
# Resolve color scheme → MEME custom alphabet (.alph) file
# ---------------------------------------------------------------------------
# MEME custom alphabet files embed letter colors at discovery time via -alph.
# All downstream outputs (HTML report, logo PNGs, meme2images) automatically
# use the embedded colors — no post-hoc re-rendering needed.
COLORS_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)/config/colors_config/meme_colors"
MOTIF_COLORS_TOML="$(cd "$SCRIPT_DIR/../../.." && pwd)/config/colors_config/meme_motif_colors.toml"
TOML_PARSER="$(cd "$SCRIPT_DIR/../../.." && pwd)/modules/utils/parse_toml.py"
ALPH_FILE=""

if [[ -n "$COLOR_SCHEME" ]]; then
    case "$COLOR_SCHEME" in
        default|nature|colorblind|dmp)
            _alph_candidate="$COLORS_DIR/${COLOR_SCHEME}_${EFFECTIVE_ALPH}.alph"
            ;;
        "") ;;
        *)  # Treat as direct path to .alph file
            _alph_candidate="$COLOR_SCHEME"
            ;;
    esac
    if [[ -n "${_alph_candidate:-}" && -f "$_alph_candidate" ]]; then
        ALPH_FILE="$_alph_candidate"
    elif [[ -n "${_alph_candidate:-}" ]]; then
        log_warn "Color scheme .alph file not found: $_alph_candidate — using default MEME colors"
    fi
    unset _alph_candidate
fi

# ---------------------------------------------------------------------------
# Resolve motif-location palette → comma-separated hex string for Python
# ---------------------------------------------------------------------------
# Accepts: a named version (wong|warm|cool|pastel|high_contrast|dark|dmp_interaction)
# which is looked up from meme_motif_colors.toml, or raw comma-separated hex codes.
RESOLVED_MOTIF_PALETTE=""
if [[ -n "$MOTIF_LOCATION_PALETTE" ]]; then
    case "$MOTIF_LOCATION_PALETTE" in
        wong|warm|cool|pastel|high_contrast|dark|dmp_interaction)
            if [[ -f "$MOTIF_COLORS_TOML" ]]; then
                # parse_toml.py returns space-separated values; convert to comma-separated
                _raw=$(python3 "$TOML_PARSER" "$MOTIF_COLORS_TOML" "$MOTIF_LOCATION_PALETTE" colors 2>/dev/null) || true
                if [[ -n "$_raw" ]]; then
                    # parse_toml.py prints list items one-per-line; join with commas
                    RESOLVED_MOTIF_PALETTE="${_raw//$'\n'/,}"
                    log_info "Motif location palette: $MOTIF_LOCATION_PALETTE (from meme_motif_colors.toml)"
                else
                    log_warn "Palette '$MOTIF_LOCATION_PALETTE' not found in $MOTIF_COLORS_TOML — using default"
                fi
                unset _raw
            else
                log_warn "meme_motif_colors.toml not found at $MOTIF_COLORS_TOML — using default palette"
            fi
            ;;
        *)
            # Treat as raw comma-separated hex codes passed through directly
            RESOLVED_MOTIF_PALETTE="$MOTIF_LOCATION_PALETTE"
            log_info "Motif location palette: custom hex codes"
            ;;
    esac
fi

# ---------------------------------------------------------------------------
# Resolve per-motif domain colors and grouped legend labels
# ---------------------------------------------------------------------------
# These are read from meme_motif_colors.toml alongside the palette colors.
# Applied only for protein (amino-acid) alphabets where structural mapping
# is meaningful.  For nucleotide datasets, domain colors are skipped and
# the plain palette (--palette) is used instead.
DOMAIN_COLORS=""
DOMAIN_LABELS=""
if [[ "$MOTIF_LOCATION_PALETTE" == "dmp_interaction" ]] && \
   [[ "$EFFECTIVE_ALPH" == "protein" ]] && \
   [[ -f "$MOTIF_COLORS_TOML" ]]; then
    _dc=$(python3 "$TOML_PARSER" "$MOTIF_COLORS_TOML" dmp_interaction domain_colors_json 2>/dev/null) || true
    _dl=$(python3 "$TOML_PARSER" "$MOTIF_COLORS_TOML" dmp_interaction domain_labels_json 2>/dev/null) || true
    [[ -n "$_dc" ]] && DOMAIN_COLORS="$_dc" && log_info "Domain colors: structural mapping loaded (dmp_interaction)"
    [[ -n "$_dl" ]] && DOMAIN_LABELS="$_dl"
    unset _dc _dl
fi

# ---------------------------------------------------------------------------
# Resolve motif-location background colour
# ---------------------------------------------------------------------------
RESOLVED_MOTIF_BG=""
if [[ -n "$MOTIF_LOCATION_BG" ]]; then
    case "$MOTIF_LOCATION_BG" in
        dark)  RESOLVED_MOTIF_BG="#111111" ;;
        light) RESOLVED_MOTIF_BG="white"   ;;
        *)     RESOLVED_MOTIF_BG="$MOTIF_LOCATION_BG" ;;
    esac
fi

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
        log_error "  mamba install -n egg -c bioconda meme"
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# Output directories
# ---------------------------------------------------------------------------
# Per-alphabet subfolder so AA and NT runs sharing the same OUTDIR (the
# orchestrator strips alphabet tokens from the label to share the parent
# folder) keep their merged FASTA, MEME, TOMTOM, FIMO, and JPEG outputs
# cleanly separated.
case "$EFFECTIVE_ALPH" in
    protein) ALPH_SUBDIR="amino_acid" ;;
    *)       ALPH_SUBDIR="nucleotide" ;;
esac
MERGED_DIR="$OUTDIR/01_merged_promoters/$ALPH_SUBDIR"
MEME_DIR="$OUTDIR/02_MEME/$ALPH_SUBDIR"
TOMTOM_DIR="$OUTDIR/03_TOMTOM/$ALPH_SUBDIR"
FIMO_DIR="$OUTDIR/04_FIMO/$ALPH_SUBDIR"

# Only create the dirs whose steps will actually run, so disabled steps
# don't leave empty 03_TOMTOM/ or 04_FIMO/ trees behind.
mkdir -p "$MERGED_DIR" "$MEME_DIR"
should_run "tomtom" && mkdir -p "$TOMTOM_DIR"
should_run "fimo"   && mkdir -p "$FIMO_DIR"

log_step "MEME Suite Promoter Motif Analysis: $LABEL"
log_info "FASTA input:  ${FASTA_FILE:-$FASTA_DIR}"
log_info "Output:       $OUTDIR"
log_info "Databases:    $DB_DIR"
log_info "Steps:        $STEPS"
log_info "Threads:      $THREADS (MEME) | max parallel (TOMTOM/FIMO): $MAX_PARALLEL"
log_info "MEME params:  nmotifs=$NMOTIFS  minw=$MINW  maxw=$MAXW  mod=$MOD  markov=$MARKOV_ORDER  time_limit=${TIME_LIMIT}s"
log_info "TOMTOM dbs:   $TOMTOM_DBS"
log_info "FIMO dbs:     $FIMO_DBS"
log_info "Overwrite:    $OVERWRITE"
if [[ -n "$ALPH_FILE" ]]; then
    log_info "Color scheme: $COLOR_SCHEME -> $ALPH_FILE (baked into MEME XML via -alph)"
else
    log_info "Color scheme: (default MEME colors, no .alph file)"
fi

wait_for_slot() { local limit="${1:-$MAX_PARALLEL}"; while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done; }

# ===========================================================================
# Step 1: Merge per-gene FASTA files into a single multi-sequence FASTA
#         (skipped when --fasta-file is given — file is used directly)
# ===========================================================================
MERGED_FASTA="$MERGED_DIR/${LABEL}_promoters.fa"

if [[ -n "$FASTA_FILE" ]]; then
    # Pre-merged input: use as-is; symlink into MERGED_DIR for traceability
    MERGED_FASTA="$FASTA_FILE"
    SEQ_COUNT=$(grep -c '^>' "$MERGED_FASTA" 2>/dev/null || echo 0)
    fasta_size_kb=$(du -k "$MERGED_FASTA" 2>/dev/null | cut -f1)
    log_info "Pre-merged FASTA: $MERGED_FASTA ($SEQ_COUNT sequences, ${fasta_size_kb:-0} KB)"
elif [[ "$OVERWRITE" == true || ! -s "$MERGED_FASTA" ]]; then
    log_step "Merging promoter FASTA files"

    # If FASTA_DIR contains no .fa files directly (new param-subdir layout: e.g. 2000up_0down/),
    # auto-resolve to the single param subdirectory.  Multiple subdirs → error (ambiguous).
    if ! find "$FASTA_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) -quit 2>/dev/null | grep -q .; then
        mapfile -t _param_subdirs < <(find "$FASTA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
        if [[ "${#_param_subdirs[@]}" -eq 1 ]]; then
            log_info "Auto-resolving param subdir: $(basename "${_param_subdirs[0]}")"
            FASTA_DIR="${_param_subdirs[0]}"
        elif [[ "${#_param_subdirs[@]}" -gt 1 ]]; then
            log_error "Multiple param subdirs found in $FASTA_DIR — pass the specific subdir as --fasta-dir"
            exit 1
        fi
    fi

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
    fasta_size_kb=$(du -k "$MERGED_FASTA" 2>/dev/null | cut -f1)
    log_info "Total promoter sequences: $SEQ_COUNT (${fasta_size_kb:-0} KB)"
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
        log_info "Time limit: ${TIME_LIMIT}s  |  Markov order: $MARKOV_ORDER"

        meme_start=$SECONDS

        # Alphabet: prefer custom .alph file (embeds letter colors) over -dna/-protein
        if [[ -n "$ALPH_FILE" ]]; then
            meme_alph_flags=(-alph "$ALPH_FILE")
        elif [[ "$EFFECTIVE_ALPH" == "protein" ]]; then
            meme_alph_flags=(-protein)
        else
            meme_alph_flags=(-dna)
        fi

        # Strip characters MEME rejects from MSA-derived FASTA inputs:
        #   '-' alignment gaps      (both alphabets)
        #   '*' stop codons         (protein only)
        # MEME chokes on gap characters in either alphabet, so this runs for
        # both protein and DNA so that aligned .fas inputs (from 04_MSA) work.
        MEME_INPUT="$MERGED_FASTA"
        DEGAPPED_FASTA="$MEME_DIR/${LABEL}_degapped.fa"
        if [[ "$EFFECTIVE_ALPH" == "protein" ]]; then
            sed '/^>/!{s/-//g; s/\*//g;}' "$MERGED_FASTA" > "$DEGAPPED_FASTA"
            log_info "Sanitized protein FASTA (stripped gaps and stop codons) -> $DEGAPPED_FASTA"
        else
            sed '/^>/!{s/-//g;}' "$MERGED_FASTA" > "$DEGAPPED_FASTA"
            log_info "Sanitized DNA FASTA (stripped gaps) -> $DEGAPPED_FASTA"
        fi
        MEME_INPUT="$DEGAPPED_FASTA"

        # Build MEME command: -objfun and -markov_order require MEME >= 5.x
        meme_cmd=(meme "$MEME_INPUT" \
            "${meme_alph_flags[@]}" \
            -oc "$MEME_OUT" \
            -time "$TIME_LIMIT" \
            -mod "$MOD" \
            -nmotifs "$NMOTIFS" \
            -minw "$MINW" \
            -maxw "$MAXW")

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
            # Suppress Open MPI warnings and avoid slot errors
            export OMPI_MCA_rmaps_base_oversubscribe=1
            export OMPI_MCA_btl="^openib"                        # no InfiniBand
            export OMPI_MCA_btl_vader_single_copy_mechanism=none   # no CMA in WSL
            meme_cmd+=(-p "$THREADS")
        elif (( THREADS > 1 )); then
            log_warn "Parallel MEME not configured (no MPI build). Running single-threaded."
        fi

        # Temporarily disable errexit so PIPESTATUS can be read even if
        # MEME returns non-zero (pipefail would otherwise abort the script
        # before the error handler below can run).
        set +e
        "${meme_cmd[@]}" 2>&1 | tee "$MEME_DIR/${LABEL}_meme.log"
        meme_exit=${PIPESTATUS[0]}
        set -e
        meme_elapsed=$(( SECONDS - meme_start ))

        if [[ $meme_exit -ne 0 ]]; then
            log_error "MEME failed (exit=$meme_exit) after ${meme_elapsed}s — see $MEME_DIR/${LABEL}_meme.log"
            exit $meme_exit
        fi

        motifs_found=0
        [[ -f "$MEME_OUT/meme.xml" ]] && motifs_found=$(grep -c '<motif ' "$MEME_OUT/meme.xml" 2>/dev/null || echo 0)
        log_info "MEME complete in ${meme_elapsed}s -> $MEME_OUT/meme.html ($motifs_found motifs discovered)"
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
        TOMTOM_PIDS=()

        IFS=',' read -ra TOMTOM_DB_LIST <<< "$TOMTOM_DBS"
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
                wait_for_slot "$MAX_PARALLEL"
                (
                    t_start=$SECONDS
                    set +e
                    tomtom \
                        -oc "$TOMTOM_OUT" \
                        -xalph \
                        -no-ssc \
                        -thresh 0.05 \
                        -min-overlap 5 \
                        "$MEME_XML" "$DB_PATH" \
                        2>&1 | tee "$TOMTOM_DIR/${LABEL}_vs_${db_name}.log"
                    t_exit=${PIPESTATUS[0]}
                    set -e
                    t_elapsed=$(( SECONDS - t_start ))
                    if [[ $t_exit -ne 0 ]]; then
                        log_error "TOMTOM failed for $db_rel (exit=$t_exit) after ${t_elapsed}s"
                    else
                        matches=0
                        [[ -f "$TOMTOM_OUT/tomtom.tsv" ]] && matches=$(( $(wc -l < "$TOMTOM_OUT/tomtom.tsv") - 1 ))
                        (( matches < 0 )) && matches=0
                        log_info "TOMTOM done in ${t_elapsed}s -> $TOMTOM_OUT/tomtom.tsv ($matches matches)"
                    fi
                    exit $t_exit
                ) &
                TOMTOM_PIDS+=("$!:$db_name")
            else
                log_info "TOMTOM results exist for $db_name, skipping"
            fi
        done
        tomtom_failed=0
        for entry in "${TOMTOM_PIDS[@]}"; do
            pid="${entry%%:*}"; db="${entry#*:}"
            if ! wait "$pid"; then
                ((tomtom_failed++)) || true
            fi
        done
        if (( tomtom_failed > 0 )); then
            log_warn "$tomtom_failed of ${#TOMTOM_PIDS[@]} TOMTOM job(s) failed — see errors above"
        fi
    fi
fi

# ===========================================================================
# Step 4: FIMO — scan sequences for known plant TF binding sites
# ===========================================================================
if should_run "fimo"; then
    log_step "FIMO known-motif scanning"
    FIMO_PIDS=()

    IFS=',' read -ra FIMO_DB_LIST <<< "$FIMO_DBS"
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
            log_warn "Skipping FIMO vs $db_rel: DNA database incompatible with protein query"
            continue
        fi
        if [[ "$EFFECTIVE_ALPH" == "dna" && "$db_is_dna" == false && -n "$db_alph_line" ]]; then
            log_warn "Skipping FIMO vs $db_rel: protein database incompatible with DNA query"
            continue
        fi

        db_name="$(basename "${db_rel%.meme}")"
        FIMO_OUT="$FIMO_DIR/${LABEL}_${db_name}"

        if [[ "$OVERWRITE" == true || ! -f "$FIMO_OUT/fimo.tsv" ]]; then
            log_info "Scanning with: $db_rel"
            wait_for_slot "$MAX_PARALLEL"
            (
                f_start=$SECONDS
                set +e
                fimo \
                    --oc "$FIMO_OUT" \
                    --thresh 1e-4 \
                    "$DB_PATH" "$MERGED_FASTA" \
                    2>&1 | tee "$FIMO_DIR/${LABEL}_${db_name}.log"
                f_exit=${PIPESTATUS[0]}
                set -e
                f_elapsed=$(( SECONDS - f_start ))
                if [[ $f_exit -ne 0 ]]; then
                    log_error "FIMO failed for $db_rel (exit=$f_exit) after ${f_elapsed}s"
                else
                    sites=0
                    [[ -f "$FIMO_OUT/fimo.tsv" ]] && sites=$(( $(wc -l < "$FIMO_OUT/fimo.tsv") - 1 ))
                    (( sites < 0 )) && sites=0
                    log_info "FIMO done in ${f_elapsed}s -> $FIMO_OUT/fimo.tsv ($sites binding sites)"
                fi
                exit $f_exit
            ) &
            FIMO_PIDS+=("$!:$db_name")
        else
            log_info "FIMO results exist for $db_name, skipping"
        fi
    done
    fimo_failed=0
    for entry in "${FIMO_PIDS[@]}"; do
        pid="${entry%%:*}"; db="${entry#*:}"
        if ! wait "$pid"; then
            ((fimo_failed++)) || true
        fi
    done
    if (( fimo_failed > 0 )); then
        log_warn "$fimo_failed of ${#FIMO_PIDS[@]} FIMO job(s) failed — see errors above"
    fi
fi

# ===========================================================================
# Step 5: JPEG — assemble motif logo grid image
# Requires: MEME output (logo*.png files); ImageMagick must be in PATH.
# Color re-rendering: if COLOR_SCHEME != "default", meme2images is used to
#   regenerate PNGs with the chosen palette before montage assembly.
# ===========================================================================
if should_run "jpeg"; then
    JPEG_DIR="$OUTDIR/05_JPEG/$ALPH_SUBDIR"
    JPEG_OUT="$JPEG_DIR/${LABEL}_motifs.jpg"

    if [[ "$OVERWRITE" == true || ! -f "$JPEG_OUT" ]]; then
        if [[ ! -f "$MEME_XML" ]]; then
            log_warn "MEME XML not found — skipping JPEG export: $MEME_XML"
        else
            mkdir -p "$JPEG_DIR"
            log_step "JPEG motif logo assembly: $LABEL"

            # Logo source: MEME output dir. When -alph was used at discovery
            # time, logos already embed the custom palette — no re-rendering.
            LOGO_SRC_DIR="$MEME_DIR/$LABEL"

            # ----------------------------------------------------------------
            # Collect forward-strand logos (logo1.png, logo2.png, ...)
            # ----------------------------------------------------------------
            mapfile -t logo_files < <(find "$LOGO_SRC_DIR" -maxdepth 1 -name 'logo[0-9]*.png' | sort -V)

            if [[ ${#logo_files[@]} -eq 0 ]]; then
                log_warn "No logo PNG files found in $LOGO_SRC_DIR — skipping JPEG assembly"
            elif ! command -v magick &>/dev/null && ! command -v montage &>/dev/null && ! command -v convert &>/dev/null; then
                log_warn "ImageMagick not found — skipping JPEG assembly"
                log_warn "Install: mamba install -n egg -c conda-forge imagemagick"
            else
                log_info "Assembling ${#logo_files[@]} motif logo(s) → $JPEG_OUT"
                log_info "Grid: ${JPEG_COLUMNS} columns | DPI: ${JPEG_DPI}"

                # Prefer IMv7 'magick' over deprecated 'convert'
                _im_cmd=""
                if command -v magick &>/dev/null; then _im_cmd="magick"
                elif command -v convert &>/dev/null; then _im_cmd="convert"
                fi

                # Estimate thumbnail size: MEME logos are ~700×200 px; scale to DPI
                thumb_w=$(( JPEG_DPI * 700 / 72 ))
                thumb_h=$(( JPEG_DPI * 200 / 72 ))

                if command -v montage &>/dev/null; then
                    montage "${logo_files[@]}" \
                        -geometry "${thumb_w}x${thumb_h}+${JPEG_LOGO_PADDING}+${JPEG_LOGO_PADDING}" \
                        -tile "${JPEG_COLUMNS}x" \
                        -background white \
                        -label '' \
                        -density "$JPEG_DPI" \
                        -quality "$JPEG_QUALITY" \
                        "$JPEG_OUT" 2>&1 | tee "$JPEG_DIR/${LABEL}_montage.log" || true
                else
                    # Fallback: vertical strip
                    $_im_cmd "${logo_files[@]}" \
                        -background white \
                        -gravity Center \
                        -append \
                        -density "$JPEG_DPI" \
                        "$JPEG_OUT" 2>&1 | tee "$JPEG_DIR/${LABEL}_montage.log" || true
                fi

                if [[ -f "$JPEG_OUT" ]]; then
                    size_kb=$(du -k "$JPEG_OUT" 2>/dev/null | cut -f1)
                    log_info "JPEG written: $JPEG_OUT (${size_kb} KB)"
                else
                    log_warn "JPEG not created — check $JPEG_DIR/${LABEL}_montage.log"
                fi
            fi

            # ----------------------------------------------------------------
            # Motif locations diagram: convert meme-motif-locations.svg → JPEG
            # Shows all sequences with coloured motif blocks (full height).
            # Prefers rsvg-convert (librsvg) for accurate SVG rendering;
            # falls back to ImageMagick convert.
            # ----------------------------------------------------------------
            SVG_SRC="$MEME_DIR/$LABEL/meme-motif-locations.svg"
            LOCATIONS_JPEG="$JPEG_DIR/${LABEL}_motif_locations.jpg"

            if [[ ! -f "$SVG_SRC" ]]; then
                log_info "meme-motif-locations.svg not generated by MEME — skipping locations JPEG (full-view PNG fallback used)"
            elif [[ "$OVERWRITE" == false && -f "$LOCATIONS_JPEG" ]]; then
                log_info "Motif locations JPEG exists, skipping: $LOCATIONS_JPEG"
            else
                # Determine pixel height proportional to sequence count
                # (SVG viewport height ≈ 12 px per sequence + header; min 400)
                seq_h=12
                # Use || true to avoid grep exit-1 (no matches) mixing with fallback
                # echo inside $() would combine both outputs into an invalid multi-line value
                svg_seq_count=$(grep -c 'class="sequence"' "$SVG_SRC" 2>/dev/null || true)
                [[ -z "$svg_seq_count" || "$svg_seq_count" -eq 0 ]] && svg_seq_count="$SEQ_COUNT"
                svg_height=$(( svg_seq_count * seq_h + 150 ))
                (( svg_height < 400 )) && svg_height=400

                # Prefer IMv7 'magick' over deprecated 'convert'
                _svg_im_cmd=""
                if command -v magick &>/dev/null; then _svg_im_cmd="magick"
                elif command -v convert &>/dev/null; then _svg_im_cmd="convert"
                fi

                if command -v rsvg-convert &>/dev/null && [[ -n "$_svg_im_cmd" ]]; then
                    log_info "Converting motif locations SVG → JPEG (rsvg-convert, DPI=$JPEG_DPI)"
                    rsvg-convert \
                        --dpi-x "$JPEG_DPI" --dpi-y "$JPEG_DPI" \
                        --format png \
                        "$SVG_SRC" \
                    | $_svg_im_cmd - \
                        -background white \
                        -flatten \
                        -quality "$JPEG_QUALITY" \
                        "$LOCATIONS_JPEG" \
                        2>&1 | tee "$JPEG_DIR/${LABEL}_locations.log" || true
                elif [[ -n "$_svg_im_cmd" ]]; then
                    log_info "Converting motif locations SVG → JPEG ($_svg_im_cmd, DPI=$JPEG_DPI)"
                    $_svg_im_cmd \
                        -density "$JPEG_DPI" \
                        -background white \
                        "$SVG_SRC" \
                        -flatten \
                        -quality "$JPEG_QUALITY" \
                        "$LOCATIONS_JPEG" \
                        2>&1 | tee "$JPEG_DIR/${LABEL}_locations.log" || true
                else
                    log_warn "Neither rsvg-convert nor magick/convert found — skipping locations JPEG"
                    log_warn "Install: mamba install -n egg -c conda-forge librsvg imagemagick"
                fi

                if [[ -f "$LOCATIONS_JPEG" ]]; then
                    sz=$(du -k "$LOCATIONS_JPEG" 2>/dev/null | cut -f1)
                    log_info "Motif locations JPEG: $LOCATIONS_JPEG (${sz} KB)"
                else
                    log_warn "Motif locations JPEG not created — check $JPEG_DIR/${LABEL}_locations.log"
                fi
            fi

        fi
    else
        log_info "JPEG exists, skipping (use --overwrite to redo): $JPEG_OUT"
    fi

    # ----------------------------------------------------------------
    # Full-view motif locations PNG: Python-rendered, shows ALL
    # sequences regardless of dataset size (unlike the MEME SVG which
    # truncates long sequence lists).
    # Single render driven by [meme].motif_location_bg (resolved into
    # RESOLVED_MOTIF_BG above): empty/white = light, non-empty = that
    # colour (e.g. dark = #111111).
    # Runs independently — not gated by the logo-grid JPEG above.
    # ----------------------------------------------------------------
    FULLVIEW_PNG="$JPEG_DIR/${LABEL}_motif_locations_full.png"
    PLOT_SCRIPT="$SCRIPT_DIR/plot_motif_locations.py"

    if [[ ! -f "$MEME_XML" ]]; then
        log_warn "meme.xml not found — skipping full-view motif locations PNG"
    elif [[ "$OVERWRITE" == false && -f "$FULLVIEW_PNG" ]]; then
        log_info "Full-view motif locations PNG exists, skipping: $FULLVIEW_PNG"
    elif ! python3 -c "import matplotlib" &>/dev/null; then
        log_warn "matplotlib not available — skipping full-view motif locations PNG"
        log_warn "Install: pip install matplotlib"
    else
        mkdir -p "$JPEG_DIR"

        _plot_cmd=(
            python3 "$PLOT_SCRIPT"
            --meme-xml "$MEME_XML"
            --outdir   "$JPEG_DIR"
            --label    "$LABEL"
            --dpi      "$JPEG_DPI"
        )
        [[ -n "$RESOLVED_MOTIF_PALETTE" ]] && _plot_cmd+=(--palette       "$RESOLVED_MOTIF_PALETTE")
        [[ -n "$RESOLVED_MOTIF_BG"      ]] && _plot_cmd+=(--bg-color      "$RESOLVED_MOTIF_BG")
        [[ -n "$PHYLO_ORDER_FILE"       ]] && _plot_cmd+=(--phylo-order   "$PHYLO_ORDER_FILE")
        [[ -n "$DOMAIN_COLORS"          ]] && _plot_cmd+=(--domain-colors "$DOMAIN_COLORS")
        [[ -n "$DOMAIN_LABELS"          ]] && _plot_cmd+=(--domain-labels "$DOMAIN_LABELS")

        log_info "Generating motif locations PNG → $FULLVIEW_PNG (bg=${RESOLVED_MOTIF_BG:-white})"
        "${_plot_cmd[@]}" \
            2>&1 | tee "$JPEG_DIR/${LABEL}_fullview.log" || \
            log_warn "Full-view plot failed — check $JPEG_DIR/${LABEL}_fullview.log"

        if [[ -f "$FULLVIEW_PNG" ]]; then
            sz=$(du -k "$FULLVIEW_PNG" 2>/dev/null | awk '{print $1}')
            log_info "Full-view PNG: $FULLVIEW_PNG (${sz} KB)"
        fi
    fi
fi

# ---------------------------------------------------------------------------
# Results summary with validation
# ---------------------------------------------------------------------------
log_step "MEME Suite analysis complete: $LABEL"
log_info "Results:"
if should_run "meme"; then
    if [[ -f "$MEME_DIR/$LABEL/meme.xml" ]]; then
        n_motifs=$(grep -c '<motif ' "$MEME_DIR/$LABEL/meme.xml" 2>/dev/null || echo 0)
        log_info "  MEME:   $MEME_DIR/$LABEL/meme.html ($n_motifs motifs)"
    else
        log_warn "  MEME:   $MEME_DIR/$LABEL/ — meme.xml NOT FOUND"
    fi
fi
if should_run "tomtom"; then
    tomtom_count=$(find "$TOMTOM_DIR" -name 'tomtom.tsv' 2>/dev/null | wc -l)
    log_info "  TOMTOM: $TOMTOM_DIR/ ($tomtom_count result file(s))"
fi
if should_run "fimo"; then
    fimo_count=$(find "$FIMO_DIR" -name 'fimo.tsv' 2>/dev/null | wc -l)
    log_info "  FIMO:   $FIMO_DIR/ ($fimo_count result file(s))"
fi
if should_run "jpeg"; then
    if [[ -f "$JPEG_DIR/${LABEL}_motifs.jpg" ]]; then
        log_info "  JPEG logos:     $JPEG_DIR/${LABEL}_motifs.jpg"
    else
        log_warn "  JPEG logos:     $JPEG_DIR/ — ${LABEL}_motifs.jpg NOT FOUND"
    fi
    if [[ -f "$JPEG_DIR/${LABEL}_motif_locations.jpg" ]]; then
        log_info "  JPEG locations: $JPEG_DIR/${LABEL}_motif_locations.jpg"
    elif [[ -f "$JPEG_DIR/${LABEL}_motif_locations_full.png" ]]; then
        log_info "  JPEG locations: N/A (full-view PNG used as fallback — see above)"
    else
        log_warn "  JPEG locations: $JPEG_DIR/ — ${LABEL}_motif_locations.jpg NOT FOUND"
    fi
    if [[ -f "$JPEG_DIR/${LABEL}_motif_locations_full.png" ]]; then
        log_info "  Full-view PNG:  $JPEG_DIR/${LABEL}_motif_locations_full.png"
    else
        log_warn "  Full-view PNG:  $JPEG_DIR/ — ${LABEL}_motif_locations_full.png NOT FOUND"
    fi
fi

# Log total output size
total_size_mb=$(du -sk "$OUTDIR" 2>/dev/null | awk '{printf "%.1f", $1/1024}' || echo "?")
log_info "Total output size: ${total_size_mb} MB"
