#!/bin/bash
# Module: Multiple Sequence Alignment
# Usage: bash align.sh --input <fasta> --method <CLUSTALO|MAFFT|MUSCLE|PROBCONS|CLUSTALW> \
#                       --outdir <dir> [--threads N] [--config <toml>]
#
# When --config is supplied, per-method parameters are read from the TOML
# (e.g. [alignment.clustalo], [alignment.mafft], …).  Without --config the
# script falls back to sensible hardcoded defaults, so it remains usable
# standalone.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../utils/logging.sh"

TOML_PARSER="$SCRIPT_DIR/../utils/parse_toml.py"

CPU=4
METHOD="CLUSTALO"
INPUT_FASTA=""
OUTPUT_DIR="."
CONFIG_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)   INPUT_FASTA="$2"; shift 2 ;;
        --method)  METHOD="$2"; shift 2 ;;
        --outdir)  OUTPUT_DIR="$2"; shift 2 ;;
        --threads) CPU="$2"; shift 2 ;;
        --config)  CONFIG_FILE="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$INPUT_FASTA" ]] && { log_error "Missing --input"; exit 1; }
[[ ! -s "$INPUT_FASTA" ]] && { log_error "Input file empty or missing: $INPUT_FASTA"; exit 1; }

# Strip any common FASTA extension (.fasta, .fa, .fas)
BASENAME=$(basename "$INPUT_FASTA" .fasta)
BASENAME=$(basename "$BASENAME" .fa)
BASENAME=$(basename "$BASENAME" .fas)

ALIGN_DIR="$OUTPUT_DIR/${METHOD}_aligned"
mkdir -p "$ALIGN_DIR"
OUTPUT_FILE="$ALIGN_DIR/${BASENAME}.fas"
ALN_FILE="$ALIGN_DIR/${BASENAME}.aln"
INPUT_COPY="$ALIGN_DIR/input_fasta.fa"
FASTA_TO_CLUSTAL="$SCRIPT_DIR/../utils/fasta_to_clustal.py"

# Keep a method-local copy of the input for traceability and downstream checks.
cp -f "$INPUT_FASTA" "$INPUT_COPY"

if [[ -s "$OUTPUT_FILE" ]]; then
    # Generate .aln if missing even when .fas is cached
    if [[ ! -s "$ALN_FILE" ]]; then
        if python3 "$FASTA_TO_CLUSTAL" "$OUTPUT_FILE" "$ALN_FILE"; then
            log_info "Generated .aln: $ALN_FILE"
        else
            log_warn ".aln conversion failed (non-fatal): $ALN_FILE"
        fi
    fi
    log_info "Alignment exists, skipping: $OUTPUT_FILE"
    echo "$OUTPUT_FILE"
    exit 0
fi

# ---------------------------------------------------------------------------
# cfg <method_key> <param_key> <default>
#   Read a value from [alignment.<method_key>].<param_key> in the TOML config.
#   Falls back to <default> when no config is provided or key is missing.
# ---------------------------------------------------------------------------
cfg() {
    local section="$1" key="$2" default="$3"
    if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
        local val
        val=$(python3 "$TOML_PARSER" "$CONFIG_FILE" alignment "$section" "$key" 2>/dev/null) || true
        [[ -n "$val" ]] && echo "$val" || echo "$default"
    else
        echo "$default"
    fi
}

log_step "Aligning $BASENAME with $METHOD"

case "$METHOD" in
    CLUSTALO)
        cmd=(clustalo -i "$INPUT_FASTA" -o "$OUTPUT_FILE"
             --outfmt="$(cfg clustalo outfmt fasta)")
        [[ "$(cfg clustalo full false)" == "true" ]]      && cmd+=(--full)
        [[ "$(cfg clustalo full_iter false)" == "true" ]]  && cmd+=(--full-iter)
        iter=$(cfg clustalo iterations 0)
        (( iter > 0 )) 2>/dev/null && cmd+=(--iter="$iter")
        gt=$(cfg clustalo max_guidetree_iterations 0)
        (( gt > 0 )) 2>/dev/null && cmd+=(--max-guidetree-iterations="$gt")
        hmm=$(cfg clustalo max_hmm_iterations 0)
        (( hmm > 0 )) 2>/dev/null && cmd+=(--max-hmm-iterations="$hmm")
        (( CPU > 1 )) && cmd+=(--threads="$CPU")
        "${cmd[@]}"
        ;;

    MAFFT)
        # Sanitize input: strip illegal chars (= < >) from sequence lines only
        cleaned_fasta="$(mktemp --suffix=.fa)"
        awk '/^>/{print; next} {gsub(/[^A-Za-z*\-]/, ""); print}' "$INPUT_FASTA" > "$cleaned_fasta"

        cmd=(mafft --thread "$CPU")
        algo=$(cfg mafft algorithm "--localpair")
        cmd+=($algo)
        cmd+=(--maxiterate "$(cfg mafft maxiterate 1000)")
        ep=$(cfg mafft ep "")
        [[ -n "$ep" ]] && cmd+=(--ep "$ep")
        [[ "$(cfg mafft anysymbol false)" == "true" ]] && cmd+=(--anysymbol)
        cmd+=("$cleaned_fasta")
        "${cmd[@]}" > "$OUTPUT_FILE"
        rm -f "$cleaned_fasta"
        ;;

    MUSCLE)
        # Detect MUSCLE version: v5 uses -align/-output/-threads; v3 uses -in/-out/-maxiters
        # Pattern handles both: "MUSCLE v3.8.1551" (v3) and "muscle 5.1.0" (v5) formats
        muscle_ver=$(muscle -version 2>&1 | grep -oiP 'muscle\s*v?\K\d+' | head -1) || muscle_ver=""
        if [[ "$muscle_ver" == "5" ]]; then
            cmd=(muscle -align "$INPUT_FASTA" -output "$OUTPUT_FILE")
            (( CPU > 1 )) && cmd+=(-threads "$CPU")
        else
            # v3 (or unknown version) flags
            cmd=(muscle -in "$INPUT_FASTA" -out "$OUTPUT_FILE")
            maxiters=$(cfg muscle maxiters 100)
            (( maxiters > 0 )) 2>/dev/null && cmd+=(-maxiters "$maxiters")
        fi
        "${cmd[@]}"
        ;;

    PROBCONS)
        probcons \
            -c  "$(cfg probcons consistency   5)" \
            -ir "$(cfg probcons iterative_ref 1000)" \
            -pre "$(cfg probcons pre_training  20)" \
            "$INPUT_FASTA" > "$OUTPUT_FILE"
        ;;

    CLUSTALW)
        # ClustalW ignores -OUTFILE for FASTA output; it writes
        # {input_basename}.fasta next to the input file.  We also
        # capture the native .aln via -OUTFILE.
        cw_bin=$(command -v clustalw2 2>/dev/null || command -v clustalw 2>/dev/null) \
            || { log_error "Neither clustalw2 nor clustalw found"; exit 1; }
        input_dir="$(dirname "$INPUT_FASTA")"
        input_base="$(basename "$INPUT_FASTA" .fasta)"
        input_base="$(basename "$input_base" .fa)"
        input_base="$(basename "$input_base" .fas)"
        cw_fasta="$input_dir/${input_base}.fasta"
        cw_aln="$input_dir/${input_base}.aln"

        cmd=("$cw_bin"
             -INFILE="$INPUT_FASTA"
             -OUTPUT="$(cfg clustalw output FASTA)")
        iter_type=$(cfg clustalw iteration "")
        [[ -n "$iter_type" ]] && cmd+=(-ITERATION="$iter_type")
        numiter=$(cfg clustalw numiter 0)
        (( numiter > 0 )) 2>/dev/null && cmd+=(-NUMITER="$numiter")
        "${cmd[@]}"

        # Move ClustalW's FASTA output to expected location
        [[ -s "$cw_fasta" ]] && mv "$cw_fasta" "$OUTPUT_FILE"
        # Move native .aln if produced (avoid overwriting our converted one later)
        [[ -s "$cw_aln" ]] && mv "$cw_aln" "$ALN_FILE"
        ;;

    *)
        log_error "Unknown method: $METHOD"; exit 1 ;;
esac

# Convert aligned FASTA to Clustal (.aln) format
if python3 "$FASTA_TO_CLUSTAL" "$OUTPUT_FILE" "$ALN_FILE"; then
    log_info "Generated .aln: $ALN_FILE"
else
    log_warn ".aln conversion failed (non-fatal): $ALN_FILE"
fi

log_info "Alignment complete -> $OUTPUT_FILE"
echo "$OUTPUT_FILE"
