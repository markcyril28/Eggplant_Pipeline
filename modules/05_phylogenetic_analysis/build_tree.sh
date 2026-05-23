#!/bin/bash
# Module: Phylogenetic Tree Construction
# Usage: bash build_tree.sh --input <aligned.fas> --software <MEGA_CC|IQTREE2|RAXML> --outdir <dir> [--config <file.mao>] [--threads N]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../utils/logging.sh"

CPU=4
SOFTWARE="MEGA_CC"
INPUT_FILE=""
OUTPUT_DIR="."
CONFIG_FILE=""
MEGACC_CONFIG_NUC=""       # nucleotide-specific MEGA_CC .mao config
MEGACC_CONFIG_AA=""        # protein-specific MEGA_CC .mao config
MEGACC_MODELSEL_NUC=""     # nucleotide model-selection .mao (optional; enables auto model fit)
MEGACC_MODELSEL_AA=""      # protein model-selection .mao (optional; enables auto model fit)
BOOTSTRAP=5000
ALRT=5000
IQTREE2_MODEL="MFP+MERGE"
IQTREE2_MODEL_NUC=""       # nucleotide-specific IQ-TREE2 model
IQTREE2_MODEL_AA=""        # protein-specific IQ-TREE2 model
IQTREE2_ALLNNI="true"
IQTREE2_POLYTOMY="false"
IQTREE2_SAFE="true"
IQTREE2_BNNI="true"
IQTREE2_FAST="false"
IQTREE2_PERS="0.05"
IQTREE2_REDO="true"
RAXML_BINARY="raxml-ng"
RAXML_MODEL="GTR+FC+R4"
RAXML_MODEL_NUC=""          # nucleotide-specific RAxML model
RAXML_MODEL_AA=""           # protein-specific RAxML model
RAXML_SEED="12345"
RAXML_MODE="all"
RAXML_BS_TREES="5000"
RAXML_SEARCH_REPLICATES="50"
RAXML_REDO="true"
GENOME_NAME=""
SUBPATH=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)     INPUT_FILE="$2"; shift 2 ;;
        --software)  SOFTWARE="$2"; shift 2 ;;
        --outdir)    OUTPUT_DIR="$2"; shift 2 ;;
        --config)    CONFIG_FILE="$2"; shift 2 ;;
        --megacc-config-nuc)   MEGACC_CONFIG_NUC="$2"; shift 2 ;;
        --megacc-config-aa)    MEGACC_CONFIG_AA="$2"; shift 2 ;;
        --megacc-modelsel-nuc) MEGACC_MODELSEL_NUC="$2"; shift 2 ;;
        --megacc-modelsel-aa)  MEGACC_MODELSEL_AA="$2"; shift 2 ;;
        --threads)   CPU="$2"; shift 2 ;;
        --bootstrap) BOOTSTRAP="$2"; shift 2 ;;
        --alrt)      ALRT="$2"; shift 2 ;;
        --model)     IQTREE2_MODEL="$2"; shift 2 ;;
        --iqtree2-model-nuc)   IQTREE2_MODEL_NUC="$2"; shift 2 ;;
        --iqtree2-model-aa)    IQTREE2_MODEL_AA="$2"; shift 2 ;;
        --allnni)    IQTREE2_ALLNNI="$2"; shift 2 ;;
        --polytomy)  IQTREE2_POLYTOMY="$2"; shift 2 ;;
        --safe)      IQTREE2_SAFE="$2"; shift 2 ;;
        --bnni)      IQTREE2_BNNI="$2"; shift 2 ;;
        --fast)      IQTREE2_FAST="$2"; shift 2 ;;
        --pers)      IQTREE2_PERS="$2"; shift 2 ;;
        --redo)      IQTREE2_REDO="$2"; shift 2 ;;
        --raxml-binary)            RAXML_BINARY="$2"; shift 2 ;;
        --raxml-model)             RAXML_MODEL="$2"; shift 2 ;;
        --raxml-model-nuc)         RAXML_MODEL_NUC="$2"; shift 2 ;;
        --raxml-model-aa)          RAXML_MODEL_AA="$2"; shift 2 ;;
        --raxml-seed)              RAXML_SEED="$2"; shift 2 ;;
        --raxml-mode)              RAXML_MODE="$2"; shift 2 ;;
        --raxml-bs-trees)          RAXML_BS_TREES="$2"; shift 2 ;;
        --raxml-search-replicates) RAXML_SEARCH_REPLICATES="$2"; shift 2 ;;
        --raxml-redo)              RAXML_REDO="$2"; shift 2 ;;
        --genome-name)             GENOME_NAME="$2"; shift 2 ;;
        --subpath)                 SUBPATH="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$INPUT_FILE" ]] && { log_error "Missing --input"; exit 1; }
[[ ! -s "$INPUT_FILE" ]] && { log_error "Input file empty: $INPUT_FILE"; exit 1; }

# Verify required tools are available (after arg parsing so $SOFTWARE is set)
case "$SOFTWARE" in
    IQTREE2)
        if ! command -v iqtree &> /dev/null; then
            log_error "IQTREE (iqtree) is not available in PATH"
            exit 1
        fi
        ;;
    RAXML)
        if ! command -v "$RAXML_BINARY" &> /dev/null; then
            log_error "RAxML-NG ($RAXML_BINARY) is not available in PATH"
            exit 1
        fi
        ;;
    MEGA_CC)
        # Need at least one valid .mao path; alphabet selection in the body picks
        # the right one. CONFIG_FILE is just the alphabet-unknown fallback.
        _have_megacc_config=0
        for _f in "$CONFIG_FILE" "$MEGACC_CONFIG_NUC" "$MEGACC_CONFIG_AA"; do
            [[ -n "$_f" && -f "$_f" ]] && { _have_megacc_config=1; break; }
        done
        if (( _have_megacc_config == 0 )); then
            log_error "MEGA_CC requires at least one valid .mao (--config / --megacc-config-nuc / --megacc-config-aa)"
            exit 1
        fi
        unset _have_megacc_config _f
        ;;
esac

BASENAME=$(basename "$INPUT_FILE")
BASENAME="${BASENAME%.*}"
# Output layout (priority: --subpath  →  --genome-name  →  flat):
#   $OUTPUT_DIR/<subpath>/<software>/   ← mirrors MSA folder structure
#   $OUTPUT_DIR/<genome>/<software>/    ← legacy, single-genome layout
#   $OUTPUT_DIR/<software>/             ← fallback
if [[ -n "$SUBPATH" ]]; then
    TREE_DIR="$OUTPUT_DIR/$SUBPATH/$SOFTWARE"
elif [[ -n "$GENOME_NAME" ]]; then
    TREE_DIR="$OUTPUT_DIR/$GENOME_NAME/$SOFTWARE"
else
    TREE_DIR="$OUTPUT_DIR/$SOFTWARE"
fi
mkdir -p "$TREE_DIR"

# Detect whether alignment content looks nucleotide-like or protein-like.
# Strategy:
#   1. Filename keywords — fast and reliable for pipeline-generated files.
#   2. Frequency-based content analysis — if >85% of non-gap residues are
#      ACGTUN, call it nucleotide. Handles alignments that contain MAFFT
#      terminal-gap chars (~), ambiguity codes, or other edge cases that
#      break a strict all-or-nothing IUPAC regex.
detect_sequence_alphabet() {
    local fasta="$1"
    local fname
    fname=$(basename "$fasta" | tr '[:upper:]' '[:lower:]')

    # --- Primary: filename keywords ---
    if [[ "$fname" == *nucleotide* || "$fname" == *_nuc_* || "$fname" == *_nt_* || "$fname" == *_dna_* ]]; then
        echo "nucleotide"; return
    fi
    if [[ "$fname" == *amino_acid* || "$fname" == *protein* || "$fname" == *_aa_* || "$fname" == *polypeptide* ]]; then
        echo "protein"; return
    fi

    # --- Fallback: frequency-based content analysis ---
    local seq_chars nuc_only total nuc_count
    # Strip headers, gaps, and ambiguity-only chars before counting
    seq_chars=$(grep -v '^>' "$fasta" | tr -d '\r\n[:space:]\-\.?\~\*' | tr '[:lower:]' '[:upper:]')
    [[ -z "$seq_chars" ]] && { echo "unknown"; return; }

    total=${#seq_chars}
    nuc_only=$(echo "$seq_chars" | tr -cd 'ACGTUN')
    nuc_count=${#nuc_only}

    # >85% ACGTUN residues → nucleotide; otherwise protein
    if (( nuc_count * 100 / total >= 85 )); then
        echo "nucleotide"
    else
        echo "protein"
    fi
}

# Deduplicate sequence names that collide after phylo tools sanitize headers
# (spaces -> _, brackets [] -> _). Appends _dup2, _dup3, etc. to collisions.
# Kept in TMPDIR so it doesn't pollute the per-software output folder.
DEDUP_FASTA=$(mktemp "${TMPDIR:-/tmp}/${BASENAME}_dedup_XXXXXX.fas")
awk '
/^>/ {
    name = substr($0, 2)
    gsub(/[^[:alnum:]_.-]+/, "_", name)
    gsub(/_+/, "_", name)
    sub(/^_+/, "", name)
    sub(/_+$/, "", name)
    if (name == "")
        name = "taxon"
    # Keep IDs in a conservative length range for downstream tool compatibility.
    if (length(name) > 200)
        name = substr(name, 1, 200)
    count[name]++
    if (count[name] > 1)
        print ">" name "_dup" count[name]
    else
        print ">" name
    next
}
{ print }
' "$INPUT_FILE" > "$DEDUP_FASTA"

cleanup_dedup() { rm -f "$DEDUP_FASTA"; }
trap cleanup_dedup EXIT

# Detect alphabet from the original input file so the filename hint is available.
SEQ_ALPHABET=$(detect_sequence_alphabet "$INPUT_FILE")

# Normalize alignment residues for downstream tools.
# - Protein mode: replace non-standard symbols (e.g., U/O) with X.
# - Nucleotide mode: convert U->T and replace invalid symbols with N.
normalize_sequence_chars() {
    local fasta="$1"
    local alphabet="$2"
    local tmp
    tmp=$(mktemp "${TMPDIR:-/tmp}/${BASENAME}_norm_XXXXXX.fas")

    awk -v mode="$alphabet" '
    /^>/ { print; next }
    {
        line = toupper($0)
        if (mode == "protein") {
            gsub(/[^ACDEFGHIKLMNPQRSTVWYBXZJ\-\?\.\*]/, "X", line)
        } else if (mode == "nucleotide") {
            gsub(/U/, "T", line)
            gsub(/[^ACGTRYSWKMBDHVN\-\?\.]/, "N", line)
        }
        print line
    }
    ' "$fasta" > "$tmp"

    mv "$tmp" "$fasta"
}

if [[ "$SEQ_ALPHABET" == "protein" || "$SEQ_ALPHABET" == "nucleotide" ]]; then
    normalize_sequence_chars "$DEDUP_FASTA" "$SEQ_ALPHABET"
fi

case "$SOFTWARE" in
    MEGA_CC)
        # Select config based on detected alphabet; fall back to generic CONFIG_FILE
        EFFECTIVE_MEGA_CONFIG="$CONFIG_FILE"
        if [[ "$SEQ_ALPHABET" == "nucleotide" && -n "$MEGACC_CONFIG_NUC" ]]; then
            EFFECTIVE_MEGA_CONFIG="$MEGACC_CONFIG_NUC"
            log_info "MEGA_CC: nucleotide input → $(basename "$EFFECTIVE_MEGA_CONFIG")"
        elif [[ "$SEQ_ALPHABET" == "protein" && -n "$MEGACC_CONFIG_AA" ]]; then
            EFFECTIVE_MEGA_CONFIG="$MEGACC_CONFIG_AA"
            log_info "MEGA_CC: protein input → $(basename "$EFFECTIVE_MEGA_CONFIG")"
        else
            log_info "MEGA_CC: $SEQ_ALPHABET input → $(basename "$EFFECTIVE_MEGA_CONFIG") (fallback)"
        fi
        [[ -z "$EFFECTIVE_MEGA_CONFIG" || ! -f "$EFFECTIVE_MEGA_CONFIG" ]] && { log_error "MEGA_CC requires a valid config file (alphabet: $SEQ_ALPHABET)"; exit 1; }

        # Optional model-selection sub-step: run "Find Best Models" first, then
        # template the chosen model into a derived inference .mao. Caches on
        # the derived .mao so reruns skip the (slow) selection scan.
        EFFECTIVE_MODELSEL_CONFIG=""
        if [[ "$SEQ_ALPHABET" == "nucleotide" && -n "$MEGACC_MODELSEL_NUC" ]]; then
            EFFECTIVE_MODELSEL_CONFIG="$MEGACC_MODELSEL_NUC"
        elif [[ "$SEQ_ALPHABET" == "protein" && -n "$MEGACC_MODELSEL_AA" ]]; then
            EFFECTIVE_MODELSEL_CONFIG="$MEGACC_MODELSEL_AA"
        elif [[ -n "$MEGACC_MODELSEL_NUC" || -n "$MEGACC_MODELSEL_AA" ]]; then
            log_info "MEGA_CC model selection skipped: no template configured for alphabet '$SEQ_ALPHABET'"
        fi
        if [[ -n "$EFFECTIVE_MODELSEL_CONFIG" && -f "$EFFECTIVE_MODELSEL_CONFIG" ]]; then
            DERIVED_MAO="$TREE_DIR/${BASENAME}_${SEQ_ALPHABET}_modelfit.mao"
            MODELSEL_BASE="$TREE_DIR/${BASENAME}_modelfit"
            MODELSEL_CSV="${MODELSEL_BASE}.csv"
            MODELSEL_LOG="${MODELSEL_BASE}.log"
            if [[ ! -s "$DERIVED_MAO" ]]; then
                log_step "MEGA_CC model selection: $BASENAME ($SEQ_ALPHABET)"
                MODELSEL_SNAPSHOT="${MODELSEL_BASE}_template.mao"
                cp -f "$EFFECTIVE_MODELSEL_CONFIG" "$MODELSEL_SNAPSHOT"
                # Inject the orchestrator's thread budget into the snapshot so
                # megacc respects it (the template's Number of Threads is fixed).
                sed -i "s|^\(Number of Threads[[:space:]]*=\).*|\1 $CPU|" "$MODELSEL_SNAPSHOT"
                rc=0
                megacc \
                    -d "$DEDUP_FASTA" \
                    -a "$MODELSEL_SNAPSHOT" \
                    -o "$MODELSEL_BASE" \
                    > "$MODELSEL_LOG" 2>&1 || rc=$?
                if (( rc != 0 )); then
                    log_error "MEGA_CC model selection failed with exit code $rc (see $MODELSEL_LOG)"
                    exit 1
                fi
                if [[ ! -s "$MODELSEL_CSV" ]]; then
                    log_error "MEGA_CC model selection produced no CSV at $MODELSEL_CSV (see $MODELSEL_LOG)"
                    exit 1
                fi
                if ! python3 "$SCRIPT_DIR/megacc_apply_model.py" \
                        --csv      "$MODELSEL_CSV" \
                        --template "$EFFECTIVE_MEGA_CONFIG" \
                        --output   "$DERIVED_MAO" \
                        --threads  "$CPU" 2>>"$MODELSEL_LOG"; then
                    log_error "megacc_apply_model.py failed (see $MODELSEL_LOG)"
                    exit 1
                fi
                log_info "Derived inference .mao: $DERIVED_MAO"
            else
                log_info "Derived .mao exists: $DERIVED_MAO (skipping model selection)"
            fi
            EFFECTIVE_MEGA_CONFIG="$DERIVED_MAO"
        fi

        # Always snapshot the inference template into the output folder so the
        # run is self-documenting whether or not model selection ran.
        cp -f "$EFFECTIVE_MEGA_CONFIG" "$TREE_DIR/${BASENAME}_${SEQ_ALPHABET}_inference_template.mao" 2>/dev/null || true

        CONFIG_BASE=$(basename "$EFFECTIVE_MEGA_CONFIG" .mao)
        # megacc -o is a base prefix, not a filename; appending .nwk to the prefix
        # would yield <base>.nwk.nwk. DERIVED_MAO already encodes BASENAME, so reuse it.
        if [[ -n "${DERIVED_MAO:-}" && "$EFFECTIVE_MEGA_CONFIG" == "$DERIVED_MAO" ]]; then
            OUTPUT_BASE="${DERIVED_MAO%.mao}"
        else
            OUTPUT_BASE="$TREE_DIR/${BASENAME}_${CONFIG_BASE}"
        fi
        OUTPUT_FILE="${OUTPUT_BASE}.nwk"
        MEGA_LOG="$TREE_DIR/${BASENAME}_MEGA.log"

        # Cache hit on either the ML best tree or the bootstrap consensus tree.
        EXISTING_TREE=""
        if [[ -s "$OUTPUT_FILE" ]]; then
            EXISTING_TREE="$OUTPUT_FILE"
        elif [[ -s "${OUTPUT_BASE}_consensus.nwk" ]]; then
            EXISTING_TREE="${OUTPUT_BASE}_consensus.nwk"
        fi
        if [[ -n "$EXISTING_TREE" ]]; then
            log_info "Tree exists: $EXISTING_TREE (skipped)"
            echo "$EXISTING_TREE"
            exit 0
        fi

        log_step "MEGA_CC: $BASENAME"
        rc=0
        # No --cpu flag: megacc reads thread count from the .mao "Number of Threads" line.
        megacc \
            -d "$DEDUP_FASTA" \
            -a "$EFFECTIVE_MEGA_CONFIG" \
            -o "$OUTPUT_BASE" \
            > "$MEGA_LOG" 2>&1 || rc=$?
        if (( rc != 0 )); then
            log_error "MEGA_CC failed with exit code $rc (see $MEGA_LOG)"
            exit 1
        fi

        # Prefer the ML best tree; fall back to the bootstrap consensus tree
        # if megacc only emitted that variant for this analysis.
        if [[ ! -s "$OUTPUT_FILE" && -s "${OUTPUT_BASE}_consensus.nwk" ]]; then
            OUTPUT_FILE="${OUTPUT_BASE}_consensus.nwk"
        fi

        if [[ -s "$OUTPUT_FILE" ]]; then
            log_info "Tree: $OUTPUT_FILE"
        else
            log_error "MEGA_CC produced no output (see $MEGA_LOG)"
            exit 1
        fi
        ;;

    IQTREE2)
        OUTPUT_PREFIX="$TREE_DIR/${BASENAME}_IQTREE2"
        TREE_FILE="${OUTPUT_PREFIX}.treefile"
        IQ_LOG="${OUTPUT_PREFIX}.log"
        IQTREE_FLAGS=()

        # Select model based on detected alphabet; fall back to generic IQTREE2_MODEL
        EFFECTIVE_IQTREE2_MODEL="$IQTREE2_MODEL"
        SEQTYPE_FLAG=()
        if [[ "$SEQ_ALPHABET" == "nucleotide" && -n "$IQTREE2_MODEL_NUC" ]]; then
            EFFECTIVE_IQTREE2_MODEL="$IQTREE2_MODEL_NUC"
            SEQTYPE_FLAG=(--seqtype DNA)
            log_info "IQ-TREE2: nucleotide input → model: $EFFECTIVE_IQTREE2_MODEL"
        elif [[ "$SEQ_ALPHABET" == "protein" && -n "$IQTREE2_MODEL_AA" ]]; then
            EFFECTIVE_IQTREE2_MODEL="$IQTREE2_MODEL_AA"
            SEQTYPE_FLAG=(--seqtype AA)
            log_info "IQ-TREE2: protein input → model: $EFFECTIVE_IQTREE2_MODEL"
        else
            log_info "IQ-TREE2: $SEQ_ALPHABET input → model: $EFFECTIVE_IQTREE2_MODEL (fallback)"
        fi

        [[ "$IQTREE2_ALLNNI" == "true" ]] && IQTREE_FLAGS+=(--allnni)
        [[ "$IQTREE2_POLYTOMY" == "true" ]] && IQTREE_FLAGS+=(--polytomy)
        [[ "$IQTREE2_SAFE" == "true" ]] && IQTREE_FLAGS+=(--safe)
        [[ "$IQTREE2_BNNI" == "true" ]] && IQTREE_FLAGS+=(--bnni)
        [[ "$IQTREE2_FAST" == "true" ]] && IQTREE_FLAGS+=(--fast)
        [[ -n "$IQTREE2_PERS" ]] && IQTREE_FLAGS+=(-pers "$IQTREE2_PERS")
        [[ "$IQTREE2_REDO" == "true" ]] && IQTREE_FLAGS+=(--redo)

        if [[ -s "$TREE_FILE" ]] && [[ "$IQTREE2_REDO" != "true" ]]; then
            log_info "Tree exists: $TREE_FILE (skipped)"
            echo "$TREE_FILE"
            exit 0
        fi

        log_step "IQ-TREE2: $BASENAME"
        rc=0
        iqtree \
            -s "$DEDUP_FASTA" \
            -m "$EFFECTIVE_IQTREE2_MODEL" \
            "${SEQTYPE_FLAG[@]:-}" \
            -T "$CPU" \
            -bb "$BOOTSTRAP" -alrt "$ALRT" \
            "${IQTREE_FLAGS[@]}" \
            -pre "$OUTPUT_PREFIX" \
            > "$IQ_LOG" 2>&1 || rc=$?
        if (( rc != 0 )); then
            log_error "IQ-TREE2 failed with exit code $rc (see $IQ_LOG)"
            exit 1
        fi

        if [[ -s "$TREE_FILE" ]]; then
            log_info "Tree: $TREE_FILE"
        else
            log_error "IQ-TREE2 produced no output (see $IQ_LOG)"
            exit 1
        fi
        ;;

    RAXML)
        OUTPUT_PREFIX="$TREE_DIR/${BASENAME}_RAXML"
        TREE_FILE="${OUTPUT_PREFIX}.raxml.bestTree"
        RAXML_LOG="${OUTPUT_PREFIX}.log"
        RAXML_FLAGS=()

        # Select model proactively based on detected alphabet
        if [[ "$SEQ_ALPHABET" == "nucleotide" ]]; then
            EFFECTIVE_RAXML_MODEL="${RAXML_MODEL_NUC:-$RAXML_MODEL}"
            log_info "RAxML: nucleotide input → model: $EFFECTIVE_RAXML_MODEL"
        elif [[ "$SEQ_ALPHABET" == "protein" ]]; then
            EFFECTIVE_RAXML_MODEL="${RAXML_MODEL_AA:-$RAXML_MODEL}"
            log_info "RAxML: protein input → model: $EFFECTIVE_RAXML_MODEL"
        else
            EFFECTIVE_RAXML_MODEL="$RAXML_MODEL"
            log_warn "RAxML: unknown alphabet, using default model: $EFFECTIVE_RAXML_MODEL"
        fi

        case "$RAXML_MODE" in
            all)
                RAXML_FLAGS+=(--all)
                RAXML_FLAGS+=(--bs-trees "$RAXML_BS_TREES")
                ;;
            search)
                RAXML_FLAGS+=(--search)
                RAXML_FLAGS+=(--tree "pars{${RAXML_SEARCH_REPLICATES}}")
                ;;
            *)
                log_error "Invalid RAxML mode: $RAXML_MODE (use: all or search)"
                exit 1
                ;;
        esac

        [[ "$RAXML_REDO" == "true" ]] && RAXML_FLAGS+=(--redo)

        if [[ -s "$TREE_FILE" ]] && [[ "$RAXML_REDO" != "true" ]]; then
            log_info "Tree exists: $TREE_FILE (skipped)"
            echo "$TREE_FILE"
            exit 0
        fi

        log_step "RAxML: $BASENAME"
        rc=0
        "$RAXML_BINARY" \
            --msa "$DEDUP_FASTA" \
            --model "$EFFECTIVE_RAXML_MODEL" \
            --threads "$CPU" \
            --force perf_threads \
            --seed "$RAXML_SEED" \
            --prefix "$OUTPUT_PREFIX" \
            "${RAXML_FLAGS[@]}" \
            > "$RAXML_LOG" 2>&1 || rc=$?
        if (( rc != 0 )); then
            log_error "RAxML failed with exit code $rc (see $RAXML_LOG)"
            exit 1
        fi

        if [[ -s "$TREE_FILE" ]]; then
            log_info "Tree: $TREE_FILE"
        else
            log_error "RAxML produced no output (see $RAXML_LOG)"
            exit 1
        fi
        ;;

    *)
        log_error "Unknown software: $SOFTWARE"; exit 1 ;;
esac
