#!/bin/bash
# ============================================================================
# Module: Phylogenetic Tree Visualization
# ============================================================================
# Renders publication-quality tree figures from Newick tree files using ggtree.
#
# Usage:
#   bash visualize_tree.sh \
#       --treedir <phylo_output_dir> \
#       --outdir <figure_output_dir> \
#       --config <merged_config.toml> \
#       [--threads N]
#
# Scans treedir for .treefile (IQ-TREE2), .contree (IQ-TREE2 consensus),
# .raxml.support (RAxML), and .nwk (MEGA_CC) files, then renders each as PNG.
# ============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../utils/logging.sh"

RENDER_SCRIPT="$SCRIPT_DIR/render_tree.R"

# Defaults (overridden by config or flags)
TREE_DIR=""
OUTPUT_DIR=""
CONFIG_FILE=""
MAX_PARALLEL=4
OVERWRITE="true"

# Visualization defaults
LAYOUT="rectangular"
SHOW_BOOTSTRAP="true"
BOOTSTRAP_THRESHOLD=70
BOOTSTRAP_STYLE="text"
LABEL_STYLE="replace"
TIP_LABEL_SIZE=4.5
TIP_LABEL_OFFSET=0.005
BOOTSTRAP_LABEL_SIZE=3.0
NODE_COLOR_HIGH="#B71C1C"
NODE_COLOR_MEDIUM="#E65100"
NODE_SIZE_HIGH=3.0
NODE_SIZE_MEDIUM=2.0
WIDTH=14
HEIGHT=10
DPI=600
BRANCH_WIDTH=0.8
HIGHLIGHT_EGGPLANT="true"
OUTGROUP_PATTERN="Outgroup"
EXCLUDE_TIPS=""
ROOT_OUTGROUP="false"
OPEN_ANGLE=15
BOOTSTRAP_COLOR="#D32F2F"
BRANCH_COLOR="grey20"
TIP_POINT_SIZE=2.5
COLOR_SMELDMP="#4527A0"
COLOR_HAPLOID="#B71C1C"
COLOR_ORTHOLOG="#37474F"
COLOR_OUTGROUP="#78909C"
TREESCALE_FONTSIZE=3.0
TREESCALE_OFFSET=0.3
TREESCALE_COLOR="grey35"

# ======================== Argument Parsing ========================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --treedir)     TREE_DIR="$2"; shift 2 ;;
        --outdir)      OUTPUT_DIR="$2"; shift 2 ;;
        --config)      CONFIG_FILE="$2"; shift 2 ;;
        --threads)     MAX_PARALLEL="$2"; shift 2 ;;
        --overwrite)   OVERWRITE="$2"; shift 2 ;;
        *) log_error "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$TREE_DIR" ]] && { log_error "Missing --treedir"; exit 1; }
[[ -d "$TREE_DIR" ]] || { log_error "Tree directory not found: $TREE_DIR"; exit 1; }

# Default output dir = treedir (figures alongside tree files)
[[ -z "$OUTPUT_DIR" ]] && OUTPUT_DIR="$TREE_DIR"
[[ ! -f "$RENDER_SCRIPT" ]] && { log_error "R script not found: $RENDER_SCRIPT"; exit 1; }

# ======================== Config-driven Overrides ========================

if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
    TOML_PARSER="$SCRIPT_DIR/../utils/parse_toml.py"
    get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

    # cfg <section> <key> <default>
    cfg() {
        local val
        val=$(get_toml "$@" 2>/dev/null) || true
        echo "${val:-}"
    }

    _val=$(cfg visualization layout);              [[ -n "$_val" ]] && LAYOUT="$_val"
    _val=$(cfg visualization show_bootstrap);       [[ -n "$_val" ]] && SHOW_BOOTSTRAP="$_val"
    _val=$(cfg visualization bootstrap_threshold);  [[ -n "$_val" ]] && BOOTSTRAP_THRESHOLD="$_val"
    _val=$(cfg visualization bootstrap_style);      [[ -n "$_val" ]] && BOOTSTRAP_STYLE="$_val"
    _val=$(cfg visualization label_style);          [[ -n "$_val" ]] && LABEL_STYLE="$_val"
    _val=$(cfg visualization tip_label_size);       [[ -n "$_val" ]] && TIP_LABEL_SIZE="$_val"
    _val=$(cfg visualization tip_label_offset);     [[ -n "$_val" ]] && TIP_LABEL_OFFSET="$_val"
    _val=$(cfg visualization bootstrap_label_size); [[ -n "$_val" ]] && BOOTSTRAP_LABEL_SIZE="$_val"
    _val=$(cfg visualization node_color_high);      [[ -n "$_val" ]] && NODE_COLOR_HIGH="$_val"
    _val=$(cfg visualization node_color_medium);    [[ -n "$_val" ]] && NODE_COLOR_MEDIUM="$_val"
    _val=$(cfg visualization node_size_high);       [[ -n "$_val" ]] && NODE_SIZE_HIGH="$_val"
    _val=$(cfg visualization node_size_medium);     [[ -n "$_val" ]] && NODE_SIZE_MEDIUM="$_val"
    _val=$(cfg visualization width);                [[ -n "$_val" ]] && WIDTH="$_val"
    _val=$(cfg visualization height);               [[ -n "$_val" ]] && HEIGHT="$_val"
    _val=$(cfg visualization dpi);                  [[ -n "$_val" ]] && DPI="$_val"
    _val=$(cfg visualization branch_width);         [[ -n "$_val" ]] && BRANCH_WIDTH="$_val"
    _val=$(cfg visualization highlight_eggplant);   [[ -n "$_val" ]] && HIGHLIGHT_EGGPLANT="$_val"
    _val=$(cfg visualization outgroup_pattern);     [[ -n "$_val" ]] && OUTGROUP_PATTERN="$_val"
    _val=$(cfg visualization exclude_tips);         [[ -n "$_val" ]] && EXCLUDE_TIPS="$_val"
    _val=$(cfg visualization root_outgroup);        [[ -n "$_val" ]] && ROOT_OUTGROUP="$_val"
    _val=$(cfg visualization open_angle);           [[ -n "$_val" ]] && OPEN_ANGLE="$_val"
    _val=$(cfg visualization bootstrap_color);      [[ -n "$_val" ]] && BOOTSTRAP_COLOR="$_val"
    _val=$(cfg visualization branch_color);         [[ -n "$_val" ]] && BRANCH_COLOR="$_val"
    _val=$(cfg visualization tip_point_size);       [[ -n "$_val" ]] && TIP_POINT_SIZE="$_val"
    _val=$(cfg visualization color_smeldmp);        [[ -n "$_val" ]] && COLOR_SMELDMP="$_val"
    _val=$(cfg visualization color_haploid);        [[ -n "$_val" ]] && COLOR_HAPLOID="$_val"
    _val=$(cfg visualization color_ortholog);       [[ -n "$_val" ]] && COLOR_ORTHOLOG="$_val"
    _val=$(cfg visualization color_outgroup);       [[ -n "$_val" ]] && COLOR_OUTGROUP="$_val"
    _val=$(cfg visualization treescale_fontsize);   [[ -n "$_val" ]] && TREESCALE_FONTSIZE="$_val"
    _val=$(cfg visualization treescale_offset);     [[ -n "$_val" ]] && TREESCALE_OFFSET="$_val"
    _val=$(cfg visualization treescale_color);      [[ -n "$_val" ]] && TREESCALE_COLOR="$_val"
    unset _val
fi

# ======================== Find Tree Files ========================

TREE_PATTERNS=(
    "*.treefile"        # IQ-TREE2 best tree (with SH-aLRT/UFBoot)
    "*.contree"         # IQ-TREE2 consensus tree
    "*.raxml.support"   # RAxML bootstrap support tree
    "*.nwk"             # MEGA_CC Newick output
)

TREE_FILES=()
for pattern in "${TREE_PATTERNS[@]}"; do
    while IFS= read -r f; do
        TREE_FILES+=("$f")
    done < <(find "$TREE_DIR" -type f -name "$pattern" ! -name "*_rooted*" 2>/dev/null | sort)
done

# Filter MEGA_CC side files: when both <base>.nwk and <base>_consensus.nwk exist
# in the same MEGA_CC folder, drop the consensus variant to avoid double-render.
FILTERED_FILES=()
for f in "${TREE_FILES[@]:-}"; do
    [[ -z "$f" ]] && continue
    if [[ "$f" == *_consensus.nwk ]]; then
        sibling="${f%_consensus.nwk}.nwk"
        if [[ -s "$sibling" ]]; then
            log_info "Skipping consensus variant: $(basename "$f") (using $(basename "$sibling"))"
            continue
        fi
    fi
    FILTERED_FILES+=("$f")
done
TREE_FILES=("${FILTERED_FILES[@]:-}")

if (( ${#TREE_FILES[@]} == 0 )); then
    log_warn "No tree files found in $TREE_DIR"
    exit 0
fi

log_info "Found ${#TREE_FILES[@]} tree file(s) to visualize"

# ======================== Render Trees ========================

wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do
        sleep 0.5
    done
}

rendered=0
skipped=0
pids=()

for treefile in "${TREE_FILES[@]}"; do
    # Derive output PNG path: preserve subdirectory structure
    rel_path="${treefile#$TREE_DIR/}"
    base_name="${rel_path%.*}"
    # Handle double extensions like .raxml.support
    if [[ "$base_name" == *.raxml ]]; then
        base_name="${base_name%.raxml}"
    fi
    out_png="$OUTPUT_DIR/${base_name}.png"

    # Skip if exists and not overwriting
    if [[ -s "$out_png" && "$OVERWRITE" != "true" ]]; then
        log_info "Skipping (exists): $(basename "$out_png")"
        skipped=$(( skipped + 1 ))
        continue
    fi

    mkdir -p "$(dirname "$out_png")"

    # Auto-generate title from filename
    title_base=$(basename "$base_name" | sed 's/_/ /g')

    wait_for_slot "$MAX_PARALLEL"

    # ---- Auto-detect phylo parameters from filename/path/config ----
    PHYLO_SOFTWARE=""
    PHYLO_MODEL=""
    PHYLO_BOOTSTRAP=""
    SEQ_TYPE=""

    # Software — from file extension
    case "$treefile" in
        *.treefile|*.contree) PHYLO_SOFTWARE="IQ-TREE2" ;;
        *.raxml.support)      PHYLO_SOFTWARE="RAxML-NG" ;;
        *.nwk)                PHYLO_SOFTWARE="MEGA_CC"  ;;
    esac

    # Sequence type — match path or filename keywords (lowercased so detection
    # works regardless of casing, including embedded keywords like
    # ".../selected_v3_full_and_HI_DMPs_nucleotide/MAFFT_aligned/IQTREE2/<stem>_NUCLEOTIDE_Sequence_IQTREE2.treefile").
    SEQ_TYPE=""
    _path_lc=$(echo "$treefile" | tr '[:upper:]' '[:lower:]')
    case "$_path_lc" in
        *amino_acid*|*protein*|*polypeptide*|*/aa/*|*_aa_*) SEQ_TYPE="AA" ;;
        *nucleotide*|*_nuc_*|*_nt_*|*_dna_*|*/nt/*) SEQ_TYPE="NT" ;;
    esac
    unset _path_lc

    # Model and bootstrap — from TOML config (if available)
    if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
        case "$PHYLO_SOFTWARE" in
            "IQ-TREE2")
                if [[ "$SEQ_TYPE" == "NT" ]]; then
                    _m=$(cfg phylogenetics iqtree2 nucleotide_model)
                else
                    _m=$(cfg phylogenetics iqtree2 protein_model)
                fi
                [[ -z "$_m" ]] && _m=$(cfg phylogenetics iqtree2 model)
                PHYLO_MODEL="${_m:-}"

                _bs=$(cfg phylogenetics iqtree2 bootstrap)
                _alrt=$(cfg phylogenetics iqtree2 alrt)
                if [[ -n "$_bs" && -n "$_alrt" ]]; then
                    PHYLO_BOOTSTRAP="UFBoot ${_bs} / SH-aLRT ${_alrt}"
                elif [[ -n "$_bs" ]]; then
                    PHYLO_BOOTSTRAP="UFBoot ${_bs}"
                fi
                ;;
            "RAxML-NG")
                if [[ "$SEQ_TYPE" == "NT" ]]; then
                    _m=$(cfg phylogenetics raxml nucleotide model)
                else
                    _m=$(cfg phylogenetics raxml protein model)
                fi
                PHYLO_MODEL="${_m:-}"

                _bs=$(cfg phylogenetics raxml bs_trees)
                [[ -n "$_bs" ]] && PHYLO_BOOTSTRAP="Standard ${_bs}"
                ;;
        esac
    fi

    (
        log_info "Rendering: $rel_path -> $(basename "$out_png")"
        rc=0
        Rscript "$RENDER_SCRIPT" \
            --input              "$treefile" \
            --output             "$out_png" \
            --title              "$title_base" \
            --layout             "$LAYOUT" \
            --show-bootstrap     "$SHOW_BOOTSTRAP" \
            --bootstrap-threshold "$BOOTSTRAP_THRESHOLD" \
            --bootstrap-style    "$BOOTSTRAP_STYLE" \
            --label-style        "$LABEL_STYLE" \
            --tip-label-size     "$TIP_LABEL_SIZE" \
            --tip-label-offset   "$TIP_LABEL_OFFSET" \
            --bootstrap-label-size "$BOOTSTRAP_LABEL_SIZE" \
            --node-color-high    "$NODE_COLOR_HIGH" \
            --node-color-medium  "$NODE_COLOR_MEDIUM" \
            --node-size-high     "$NODE_SIZE_HIGH" \
            --node-size-medium   "$NODE_SIZE_MEDIUM" \
            --width              "$WIDTH" \
            --height             "$HEIGHT" \
            --dpi                "$DPI" \
            --branch-width       "$BRANCH_WIDTH" \
            --highlight-eggplant "$HIGHLIGHT_EGGPLANT" \
            --outgroup-pattern   "$OUTGROUP_PATTERN" \
            --exclude-tips       "$EXCLUDE_TIPS" \
            --root-outgroup      "$ROOT_OUTGROUP" \
            --open-angle         "$OPEN_ANGLE" \
            --bootstrap-color    "$BOOTSTRAP_COLOR" \
            --branch-color       "$BRANCH_COLOR" \
            --tip-point-size     "$TIP_POINT_SIZE" \
            --color-smeldmp      "$COLOR_SMELDMP" \
            --color-haploid      "$COLOR_HAPLOID" \
            --color-ortholog     "$COLOR_ORTHOLOG" \
            --color-outgroup     "$COLOR_OUTGROUP" \
            --treescale-fontsize "$TREESCALE_FONTSIZE" \
            --treescale-offset   "$TREESCALE_OFFSET" \
            --treescale-color    "$TREESCALE_COLOR" \
            --phylo-software     "$PHYLO_SOFTWARE" \
            --phylo-model        "$PHYLO_MODEL" \
            --phylo-bootstrap    "$PHYLO_BOOTSTRAP" \
            --sequence-type      "$SEQ_TYPE" \
            2>&1 || rc=$?

        if (( rc != 0 )); then
            log_error "Failed to render: $rel_path (exit code $rc)"
            exit 1
        fi

        if [[ -s "$out_png" ]]; then
            log_info "Rendered: $(basename "$out_png") ($(du -sh "$out_png" | cut -f1))"
        else
            log_error "Render produced empty output: $out_png"
            exit 1
        fi
    ) &
    pids+=("$!")
    rendered=$(( rendered + 1 ))
done

fail_count=0
for p in "${pids[@]:-}"; do
    [[ -z "$p" ]] && continue
    wait "$p" || fail_count=$(( fail_count + 1 ))
done

log_info "Visualization summary: rendered=$rendered, skipped=$skipped, failed=$fail_count"

if (( fail_count > 0 )); then
    log_error "$fail_count visualization(s) failed"
    exit 1
fi
