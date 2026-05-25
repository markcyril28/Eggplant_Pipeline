#!/bin/bash
# ============================================================================
# Module: Gene-Specific Structure Processing (full 5-operation pipeline)
# ============================================================================
# Operations (comment in/out the DO_* flags below to toggle):
#   1. Extract AlphaFold3 zip files
#   2. Copy model 0 to AlphaFold3_Results/
#   3. Convert .cif → .pdb  (gemmi)
#   4. Update PDB HEADER with gene name
#   5. Generate publication-quality images (PyMOL)
#   6. Extract quality metrics (AF3 pLDDT/ptm vs SWISS GMQE/QMEAN)
#   7. Structural alignment AF3 vs SWISS (PyMOL super)
#   8. Comparative renders (overlay, confidence, deviation)
#   9. Summary comparison report + matplotlib figures
#
# Usage (orchestrated):
#   bash process_structures.sh \
#       --input-dir /path/to/run_dir \
#       --color-config /path/to/protein_structure_colors.toml \
#       --gene-group DMP-HAP2 \
#       --overwrite true
#
# Usage (standalone):
#   bash process_structures.sh --input-dir /path/to/run_dir
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Operations: defaults allow orchestrator env-var override ─────────────────
# Ops 1-5 default true (standalone); orchestrator sets them via exported env vars.
# Ops 6-9 default false; orchestrator exports them when the TOML operation list enables them.
DO_EXTRACT=${DO_EXTRACT:-true}              # 1. Extract AlphaFold3 zip files
DO_COPY_MODEL0=${DO_COPY_MODEL0:-true}     # 2. Copy model 0 to AlphaFold3_Results/
DO_CIF_TO_PDB=${DO_CIF_TO_PDB:-true}      # 3. Convert .cif → .pdb (gemmi, AlphaFold3_Results/ only)
DO_CIF_TO_PDB_ALL=${DO_CIF_TO_PDB_ALL:-false}  # 3b. Convert ALL .cif → .pdb recursively under run dir
DO_UPDATE_HEADER=${DO_UPDATE_HEADER:-true} # 4. Update PDB HEADER with gene name
DO_RENDER=${DO_RENDER:-true}               # 5. Generate publication-quality images
DO_EXTRACT_METRICS=${DO_EXTRACT_METRICS:-false}  # 6. Extract quality metrics (AF3 vs SWISS)
DO_STRUCT_ALIGN=${DO_STRUCT_ALIGN:-false}        # 7. Structural alignment (AF3 vs SWISS, PyMOL)
DO_COMPARE_RENDER=${DO_COMPARE_RENDER:-false}    # 8. Comparative visualisation (PyMOL)
DO_COMPARE_REPORT=${DO_COMPARE_REPORT:-false}    # 9. Summary comparison report + figures

# ── Defaults (standalone mode) ──────────────────────────────────────────────
INPUT_DIR=""
COLOR_CONFIG=""
GENE_GROUP="DMP"
OVERWRITE="false"
THREADS=4                # parallel workers for CIF→PDB, PDB header, and render
BACKGROUNDS="black white"       # space-separated list; passed to comparative_render.py
COLOR_VERSIONS=""               # space-separated color versions; empty = all versions from color config
COMPARE_RENDER_TYPES=""         # space-separated Op 8 sub-types (overlay confidence deviation); empty → default "overlay"
MULTI_MODEL_GENES=""            # space-separated substrings; matching genes also get extra models
MULTI_MODEL_NUMBERS="1"         # space-separated model numbers to copy for multi-model genes

wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do sleep 0.5; done
}

# ── Parse flags ─────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)     INPUT_DIR="$2";     shift 2 ;;
        --color-config)  COLOR_CONFIG="$2";  shift 2 ;;
        --gene-group)    GENE_GROUP="$2";    shift 2 ;;
        --overwrite)     OVERWRITE="$2";     shift 2 ;;
        --threads)       THREADS="$2";       shift 2 ;;
        --backgrounds)          BACKGROUNDS="$2";           shift 2 ;;
        --color-versions)       COLOR_VERSIONS="$2";        shift 2 ;;
        --compare-render-types) COMPARE_RENDER_TYPES="$2";  shift 2 ;;
        --multi-model-genes)    MULTI_MODEL_GENES="$2";     shift 2 ;;
        --multi-model-numbers)  MULTI_MODEL_NUMBERS="$2";   shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$INPUT_DIR" ]] && { echo "ERROR: --input-dir is required"; exit 1; }
[[ ! -d "$INPUT_DIR" ]] && { echo "ERROR: directory does not exist: $INPUT_DIR"; exit 1; }

# ── Locate / create the AlphaFold3_Results output directory ─────────────────
# If INPUT_DIR already IS an AlphaFold3_Results (or legacy Protein_Structures)
# directory, treat its parent as the run directory.
RUN_DIR="$INPUT_DIR"
STRUCT_DIR=""
basename_input="$(basename "$INPUT_DIR")"
if [[ "$basename_input" == AlphaFold3_Results || "$basename_input" == Protein_Structures || "$basename_input" == Protein_Structure_* ]]; then
    STRUCT_DIR="$INPUT_DIR"
    RUN_DIR="$(dirname "$INPUT_DIR")"
else
    # Search for existing AlphaFold3_Results (or legacy names) inside
    for d in "$INPUT_DIR"/AlphaFold3_Results "$INPUT_DIR"/Protein_Structure*; do
        if [[ -d "$d" ]]; then
            STRUCT_DIR="$d"
            break
        fi
    done
fi
# Default: create AlphaFold3_Results/ inside the run directory
[[ -z "$STRUCT_DIR" ]] && STRUCT_DIR="$RUN_DIR/AlphaFold3_Results"

# ── Helper: strip AlphaFold3 timestamp suffix ────────────────────────────────
# Removes trailing _YYYY_MM_DD_HH_MM pattern from a name
strip_timestamp() {
    echo "$1" | sed -E 's/_[0-9]{4}_[0-9]{2}_[0-9]{2}_[0-9]{2}_[0-9]{2}$//'
}

# ════════════════════════════════════════════════════════════════════════════
# Op 1: Extract AlphaFold3 zip files
# ════════════════════════════════════════════════════════════════════════════
if [[ "${DO_EXTRACT:-}" == "true" ]]; then
    echo "=== Op 1: Extracting zip files ==="
    found_zip=false
    for zipfile in "$RUN_DIR"/*.zip; do
        [[ -f "$zipfile" ]] || continue
        found_zip=true
        dirname_only="$(basename "${zipfile%.zip}")"
        if [[ -d "$RUN_DIR/$dirname_only" && "$OVERWRITE" != "true" ]]; then
            echo "  Already extracted: $dirname_only — skipping"
            continue
        fi
        echo "  Extracting: $(basename "$zipfile") → $dirname_only/"
        mkdir -p "$RUN_DIR/$dirname_only"
        unzip -q -o "$zipfile" -d "$RUN_DIR/$dirname_only"
    done
    [[ "$found_zip" == "false" ]] && echo "  No zip files found — nothing to extract"
fi

# ════════════════════════════════════════════════════════════════════════════
# Op 2: Copy model 0 to Protein_Structures/
# ════════════════════════════════════════════════════════════════════════════
if [[ "${DO_COPY_MODEL0:-}" == "true" ]]; then
    echo ""
    echo "=== Op 2: Copying model 0 to AlphaFold3_Results/ ==="
    mkdir -p "$STRUCT_DIR"
    found_models=false
    # Recursively find all model_0.cif files under RUN_DIR, excluding the
    # output Protein_Structures directory.  Derive gene name from the
    # immediate parent folder of each CIF (strip AlphaFold3 timestamp).
    while IFS= read -r model_cif; do
        found_models=true
        parent_name="$(basename "$(dirname "$model_cif")")"
        gene_name=$(strip_timestamp "$parent_name")
        dest="$STRUCT_DIR/$gene_name"
        mkdir -p "$dest"
        if [[ -f "$dest/$(basename "$model_cif")" && "$OVERWRITE" != "true" ]]; then
            echo "  Already copied: $gene_name — skipping"
            continue
        fi
        cp "$model_cif" "$dest/"
        echo "  Copied: $gene_name/$(basename "$model_cif")"
    done < <(find "$RUN_DIR" -path "$STRUCT_DIR" -prune -o -name "*_model_0.cif" -print 2>/dev/null)
    [[ "$found_models" == "false" ]] && echo "  No model_0.cif files found — nothing to copy"

    # ── Also copy extra models for genes listed in MULTI_MODEL_GENES ─────
    if [[ -n "$MULTI_MODEL_GENES" ]]; then
        for model_num in $MULTI_MODEL_NUMBERS; do
            while IFS= read -r model_cif; do
                parent_name="$(basename "$(dirname "$model_cif")")"
                gene_name=$(strip_timestamp "$parent_name")
                # Check if this gene matches any multi-model gene substring
                match=false
                for mm_gene in $MULTI_MODEL_GENES; do
                    if [[ "$gene_name" == *"$mm_gene"* ]]; then
                        match=true
                        break
                    fi
                done
                [[ "$match" == "false" ]] && continue
                dest="${STRUCT_DIR}/${gene_name}_model_${model_num}"
                mkdir -p "$dest"
                if [[ -f "$dest/$(basename "$model_cif")" && "$OVERWRITE" != "true" ]]; then
                    echo "  Already copied: ${gene_name}_model_${model_num} — skipping"
                    continue
                fi
                cp "$model_cif" "$dest/"
                echo "  Copied: ${gene_name}_model_${model_num}/$(basename "$model_cif")"
            done < <(find "$RUN_DIR" -path "$STRUCT_DIR" -prune -o -name "*_model_${model_num}.cif" -print 2>/dev/null)
        done
    fi
fi

# ════════════════════════════════════════════════════════════════════════════
# Op 3: CIF → PDB conversion
# ════════════════════════════════════════════════════════════════════════════
if [[ "${DO_CIF_TO_PDB:-}" == "true" ]]; then
    echo ""
    echo "=== Op 3: CIF → PDB conversion ==="
    if [[ -d "$STRUCT_DIR" ]]; then
        python3 "$SCRIPT_DIR/cif_to_pdb.py" --input-dir "$STRUCT_DIR" --workers "$THREADS"
    else
        echo "  $(basename "$STRUCT_DIR")/ not found — run copy-model-0 first"
    fi
fi

# ════════════════════════════════════════════════════════════════════════════
# Op 3b: CIF → PDB conversion (recursive, entire run directory)
# Converts every .cif under RUN_DIR (including raw AlphaFold3 outputs and
# templates/ subfolders). Skips reference.cif. Useful when you need PDBs of
# template hits or extra models, not just the curated model_0 copies.
# ════════════════════════════════════════════════════════════════════════════
if [[ "${DO_CIF_TO_PDB_ALL:-}" == "true" ]]; then
    echo ""
    echo "=== Op 3b: CIF → PDB conversion (recursive, run-wide) ==="
    python3 "$SCRIPT_DIR/cif_to_pdb.py" --input-dir "$RUN_DIR" --workers "$THREADS"
fi

# ════════════════════════════════════════════════════════════════════════════
# Op 4: Update PDB HEADER with gene name
# ════════════════════════════════════════════════════════════════════════════
if [[ "${DO_UPDATE_HEADER:-}" == "true" ]]; then
    echo ""
    echo "=== Op 4: Updating PDB HEADER ==="
    if [[ ! -d "$STRUCT_DIR" ]]; then
        echo "  $(basename "$STRUCT_DIR")/ not found — skipping"
    else
        for gene_dir in "$STRUCT_DIR"/*/; do
            gene_dir="${gene_dir%/}"
            [[ -d "$gene_dir" ]] || continue
            gene_name="$(basename "$gene_dir")"
            [[ "$gene_name" == __pycache__ ]] && continue
            target_pdb="$gene_dir/${gene_name}.pdb"
            # Find the converted model PDB as source (prefer model_0, fall back to any model_N)
            src_pdb=$(find "$gene_dir" -maxdepth 1 -name "*_model_0.pdb" -print -quit 2>/dev/null)
            [[ -z "$src_pdb" ]] && src_pdb=$(find "$gene_dir" -maxdepth 1 -name "*_model_[0-9].pdb" -print -quit 2>/dev/null)
            [[ -z "$src_pdb" ]] && continue
            if [[ -f "$target_pdb" && "$OVERWRITE" != "true" ]]; then
                echo "  Already exists: ${gene_name}.pdb — skipping"
                continue
            fi
            wait_for_slot "$THREADS"
            (
                # Uppercase gene name for HEADER classification field
                header_name=$(echo "$gene_name" | tr '[:lower:]' '[:upper:]')
                # Extract hash ID from the original HEADER (position 63+)
                hash_id=""
                if head -1 "$src_pdb" | grep -q "^HEADER"; then
                    hash_id=$(head -1 "$src_pdb" | cut -c63- | tr -d '[:space:]')
                fi
                {
                    printf "HEADER    %-50s%s\n" "$header_name" "$hash_id"
                    tail -n +2 "$src_pdb"
                } > "$target_pdb"
                echo "  Created: ${gene_name}.pdb"
            ) &
        done
        wait
    fi
fi

# ════════════════════════════════════════════════════════════════════════════
# Op 5: Generate publication-quality images
# ════════════════════════════════════════════════════════════════════════════
if [[ "${DO_RENDER:-}" == "true" ]]; then
    echo ""
    echo "=== Op 5: Rendering structures ==="
    if [[ ! -d "$STRUCT_DIR" ]]; then
        echo "  $(basename "$STRUCT_DIR")/ not found — skipping"
    else
        if [[ "$OVERWRITE" == "true" ]] || \
           ! find "$STRUCT_DIR" -name "*.jpg" -type f -print -quit 2>/dev/null | grep -q .; then
            # Select renderer based on gene group
            case "$GENE_GROUP" in
                DMP-HAP2)   RENDERER="$SCRIPT_DIR/render_structures_dmp_hap2.py" ;;
                HAP2)       RENDERER="$SCRIPT_DIR/render_structures_hap2.py" ;;
                *)          RENDERER="$SCRIPT_DIR/render_structures_dmp.py" ;;
            esac

            RENDER_ARGS=(--input-dir "$STRUCT_DIR" --workers "$THREADS")
            if [[ -n "$COLOR_CONFIG" && -f "$COLOR_CONFIG" ]]; then
                RENDER_ARGS+=(--color-config "$COLOR_CONFIG")
            fi
            # HAP2 renderer accepts --gene-group
            case "$GENE_GROUP" in
                HAP2) RENDER_ARGS+=(--gene-group "HAP2") ;;
            esac
            # Subset color versions when requested via --color-versions
            [[ -n "$COLOR_VERSIONS" ]] && RENDER_ARGS+=(--color-versions "$COLOR_VERSIONS")
            # Pass backgrounds from TOML (overrides color config default)
            [[ -n "$BACKGROUNDS" ]] && RENDER_ARGS+=(--backgrounds "$BACKGROUNDS")

            echo "  Using renderer: $(basename "$RENDERER")"
            python3 "$RENDERER" "${RENDER_ARGS[@]}"
        else
            echo "  Images already exist — skipping (OVERWRITE=$OVERWRITE)"
        fi  # end overwrite / no-existing-images guard
    fi
fi

# ════════════════════════════════════════════════════════════════════════════
# Op 6: Extract quality metrics (AF3 vs SWISS)
# ════════════════════════════════════════════════════════════════════════════
if [[ "${DO_EXTRACT_METRICS:-}" == "true" ]]; then
    echo ""
    echo "=== Op 6: Extracting quality metrics ==="
    if [[ ! -d "$RUN_DIR/AlphaFold3_Results" || ! -d "$RUN_DIR/SWISS_Results" ]]; then
        echo "  Missing AlphaFold3_Results/ or SWISS_Results/ in $(basename "$RUN_DIR") — skipping"
    else
        python3 "$SCRIPT_DIR/extract_quality_metrics.py" \
            --run-dir "$RUN_DIR" \
            --overwrite "$OVERWRITE"
    fi
fi

# ════════════════════════════════════════════════════════════════════════════
# Op 7: Structural alignment AF3 vs SWISS (PyMOL super)
# ════════════════════════════════════════════════════════════════════════════
if [[ "${DO_STRUCT_ALIGN:-}" == "true" ]]; then
    echo ""
    echo "=== Op 7: Structural alignment (AF3 vs SWISS) ==="
    if [[ ! -d "$RUN_DIR/AlphaFold3_Results" || ! -d "$RUN_DIR/SWISS_Results" ]]; then
        echo "  Missing AlphaFold3_Results/ or SWISS_Results/ in $(basename "$RUN_DIR") — skipping"
    else
        python3 "$SCRIPT_DIR/structural_alignment.py" \
            --run-dir "$RUN_DIR" \
            --overwrite "$OVERWRITE"
    fi
fi

# ════════════════════════════════════════════════════════════════════════════
# Op 8: Comparative renders (overlay, confidence, deviation)
# ════════════════════════════════════════════════════════════════════════════
if [[ "${DO_COMPARE_RENDER:-}" == "true" ]]; then
    echo ""
    echo "=== Op 8: Comparative renders ==="
    if [[ ! -d "$RUN_DIR/AlphaFold3_Results" || ! -d "$RUN_DIR/SWISS_Results" ]]; then
        echo "  Missing AlphaFold3_Results/ or SWISS_Results/ in $(basename "$RUN_DIR") — skipping"
    else
        COMPARE_ARGS=(--run-dir "$RUN_DIR" --overwrite "$OVERWRITE")
        [[ -n "${COLOR_CONFIG:-}" && -f "${COLOR_CONFIG:-}" ]] && COMPARE_ARGS+=(--color-config "$COLOR_CONFIG")
        [[ -n "${BACKGROUNDS:-}" ]] && COMPARE_ARGS+=(--backgrounds "$BACKGROUNDS")
        [[ -n "${COMPARE_RENDER_TYPES:-}" ]] && COMPARE_ARGS+=(--render-types "$COMPARE_RENDER_TYPES")
        python3 "$SCRIPT_DIR/comparative_render.py" "${COMPARE_ARGS[@]}"
    fi
fi

# ════════════════════════════════════════════════════════════════════════════
# Op 9: Summary comparison report + figures
# ════════════════════════════════════════════════════════════════════════════
if [[ "${DO_COMPARE_REPORT:-}" == "true" ]]; then
    echo ""
    echo "=== Op 9: Comparison report ==="
    if [[ ! -d "$RUN_DIR/Comparison_Results" ]]; then
        echo "  Comparison_Results/ not found — run Ops 6-8 first"
    else
        python3 "$SCRIPT_DIR/comparison_report.py" \
            --run-dir "$RUN_DIR" \
            --overwrite "$OVERWRITE"
    fi
fi
