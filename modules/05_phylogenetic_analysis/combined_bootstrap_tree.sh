#!/bin/bash
# ============================================================================
# Module: Combined Bootstrap Phylogenetic Tree
# ============================================================================
# Produces a single tree figure (IQ-TREE2 topology) with bootstrap support
# from both IQ-TREE2 (UFBoot) and RAxML-NG overlaid at each internal node
# as "UFBoot/RAX-BS" labels (e.g., "98/95").
#
# Usage (standalone):
#   bash combined_bootstrap_tree.sh \
#       --treedir <phylo_output_dir> \
#       --outdir  <figure_output_dir> \
#       --config  <merged_config.toml> \
#       [--threads N] [--overwrite true|false]
#
# Usage (from orchestrator):
#   bash "$MODULES/05_phylogenetic_analysis/combined_bootstrap_tree.sh" \
#       --treedir  "$PHYLO_DIR" \
#       --outdir   "$PHYLO_DIR" \
#       --config   "$CONFIG_FILE" \
#       --threads  "$CPU" \
#       --overwrite "$OVERWRITE"
#
# Pair-matching logic (same as compare_trees.sh):
#   For every file matching:
#     <treedir>/<...subpath...>/IQTREE2/<stem>_IQTREE2.treefile
#   where <subpath> mirrors the MSA folder hierarchy
#   (e.g. <genome>/<output_subdir>/<set_name>/<METHOD>_aligned/), and falls
#   back to a single <genome>/ component for the legacy flat layout.
#   The script looks for the sibling RAxML output:
#     <treedir>/<...subpath...>/RAXML/<stem>_RAXML.raxml.support
#   Pairs with no .raxml.support counterpart are skipped (bootstrap values
#   required; .raxml.bestTree has no support annotations).
#
# Output: <subpath>/Combined/{Nucleotide|Protein}/<stem>_combined_bootstrap.png
#         (written next to the tree pair, not under --outdir).
# ============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../utils/logging.sh"

RENDER_R="$SCRIPT_DIR/render_combined_bootstrap.R"

# ======================== Defaults ========================

TREE_DIR=""
OUTPUT_DIR=""
CONFIG_FILE=""
MAX_PARALLEL=4
OVERWRITE="true"

# Visualization defaults (overridable via [visualization] TOML section)
LAYOUT="rectangular"
OPEN_ANGLE=15
DPI=600
WIDTH=14
HEIGHT=10
ROOT_OUTGROUP="false"
OUTGROUP_PATTERN="Outgroup"
EXCLUDE_TIPS=""
HIGHLIGHT_EGGPLANT="true"
SHOW_BOOTSTRAP="true"
BOOTSTRAP_THRESHOLD=70
BOOTSTRAP_STYLE="text"
BOOTSTRAP_LABEL_SIZE=3.0
BOOTSTRAP_COLOR="#D32F2F"
NODE_COLOR_HIGH="#B71C1C"
NODE_COLOR_MEDIUM="#E65100"
NODE_SIZE_HIGH=3.0
NODE_SIZE_MEDIUM=2.0
TIP_LABEL_SIZE=4.5
TIP_LABEL_OFFSET=0.005
XLIM_EXPAND=1.0
TIP_POINT_SIZE=2.5
BRANCH_WIDTH=0.8
BRANCH_COLOR="grey20"
LABEL_STYLE="replace"
COMBINED_STYLE="fraction"
COMBINED_SEP="/"
TREE_TITLE=""
TREESCALE_FONTSIZE=3.0
TREESCALE_OFFSET=0.3
TREESCALE_COLOR="grey35"
COLOR_SMELDMP="#4527A0"
COLOR_HAPLOID="#B71C1C"
COLOR_ORTHOLOG="#37474F"
COLOR_OUTGROUP="#78909C"

# ======================== Argument Parsing ========================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --treedir)   TREE_DIR="$2";     shift 2 ;;
        --outdir)    OUTPUT_DIR="$2";   shift 2 ;;
        --config)    CONFIG_FILE="$2";  shift 2 ;;
        --threads)   MAX_PARALLEL="$2"; shift 2 ;;
        --overwrite) OVERWRITE="$2";    shift 2 ;;
        *) log_error "Unknown argument: $1"; exit 1 ;;
    esac
done

[[ -z "$TREE_DIR" ]]   && { log_error "Missing --treedir"; exit 1; }
[[ -d "$TREE_DIR" ]]   || { log_error "Tree directory not found: $TREE_DIR"; exit 1; }
[[ ! -f "$RENDER_R" ]] && { log_error "R script not found: $RENDER_R"; exit 1; }

[[ -z "$OUTPUT_DIR" ]] && OUTPUT_DIR="$TREE_DIR"

# ======================== Config-driven Overrides ========================

if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
    TOML_PARSER="$SCRIPT_DIR/../utils/parse_toml.py"
    get_toml() { python3 "$TOML_PARSER" "$CONFIG_FILE" "$@"; }

    cfg() {
        local val
        val=$(get_toml "$@" 2>/dev/null) || true
        echo "${val:-}"
    }

    _val=$(cfg visualization layout);                          [[ -n "$_val" ]] && LAYOUT="$_val"
    _val=$(cfg visualization open_angle);                      [[ -n "$_val" ]] && OPEN_ANGLE="$_val"
    _val=$(cfg visualization dpi);                             [[ -n "$_val" ]] && DPI="$_val"
    _val=$(cfg visualization width);                           [[ -n "$_val" ]] && WIDTH="$_val"
    _val=$(cfg visualization height);                          [[ -n "$_val" ]] && HEIGHT="$_val"
    _val=$(cfg visualization root_outgroup);                   [[ -n "$_val" ]] && ROOT_OUTGROUP="$_val"
    _val=$(cfg visualization outgroup_pattern);                [[ -n "$_val" ]] && OUTGROUP_PATTERN="$_val"
    _val=$(cfg visualization exclude_tips);                    [[ -n "$_val" ]] && EXCLUDE_TIPS="$_val"
    _val=$(cfg visualization highlight_eggplant);              [[ -n "$_val" ]] && HIGHLIGHT_EGGPLANT="$_val"
    _val=$(cfg visualization show_bootstrap);                  [[ -n "$_val" ]] && SHOW_BOOTSTRAP="$_val"
    _val=$(cfg visualization bootstrap_threshold);             [[ -n "$_val" ]] && BOOTSTRAP_THRESHOLD="$_val"
    _val=$(cfg visualization bootstrap_style);                 [[ -n "$_val" ]] && BOOTSTRAP_STYLE="$_val"
    _val=$(cfg visualization bootstrap_label_size);            [[ -n "$_val" ]] && BOOTSTRAP_LABEL_SIZE="$_val"
    _val=$(cfg visualization bootstrap_color);                 [[ -n "$_val" ]] && BOOTSTRAP_COLOR="$_val"
    _val=$(cfg visualization node_color_high);                 [[ -n "$_val" ]] && NODE_COLOR_HIGH="$_val"
    _val=$(cfg visualization node_color_medium);               [[ -n "$_val" ]] && NODE_COLOR_MEDIUM="$_val"
    _val=$(cfg visualization node_size_high);                  [[ -n "$_val" ]] && NODE_SIZE_HIGH="$_val"
    _val=$(cfg visualization node_size_medium);                [[ -n "$_val" ]] && NODE_SIZE_MEDIUM="$_val"
    _val=$(cfg visualization tip_label_size);                  [[ -n "$_val" ]] && TIP_LABEL_SIZE="$_val"
    _val=$(cfg visualization tip_label_offset);                [[ -n "$_val" ]] && TIP_LABEL_OFFSET="$_val"
    _val=$(cfg visualization tip_point_size);                  [[ -n "$_val" ]] && TIP_POINT_SIZE="$_val"
    _val=$(cfg visualization branch_width);                    [[ -n "$_val" ]] && BRANCH_WIDTH="$_val"
    _val=$(cfg visualization branch_color);                    [[ -n "$_val" ]] && BRANCH_COLOR="$_val"
    _val=$(cfg visualization label_style);                     [[ -n "$_val" ]] && LABEL_STYLE="$_val"
    _val=$(cfg visualization treescale_fontsize);              [[ -n "$_val" ]] && TREESCALE_FONTSIZE="$_val"
    _val=$(cfg visualization treescale_offset);                [[ -n "$_val" ]] && TREESCALE_OFFSET="$_val"
    _val=$(cfg visualization treescale_color);                 [[ -n "$_val" ]] && TREESCALE_COLOR="$_val"
    _val=$(cfg visualization color_smeldmp);                   [[ -n "$_val" ]] && COLOR_SMELDMP="$_val"
    _val=$(cfg visualization color_haploid);                   [[ -n "$_val" ]] && COLOR_HAPLOID="$_val"
    _val=$(cfg visualization color_ortholog);                  [[ -n "$_val" ]] && COLOR_ORTHOLOG="$_val"
    _val=$(cfg visualization color_outgroup);                  [[ -n "$_val" ]] && COLOR_OUTGROUP="$_val"
    _val=$(cfg visualization xlim_expand);                     [[ -n "$_val" ]] && XLIM_EXPAND="$_val"

    # Combined-bootstrap-tree-specific overrides
    _val=$(cfg visualization combined_bootstrap_tree dpi);             [[ -n "$_val" ]] && DPI="$_val"
    _val=$(cfg visualization combined_bootstrap_tree width);           [[ -n "$_val" ]] && WIDTH="$_val"
    _val=$(cfg visualization combined_bootstrap_tree height);          [[ -n "$_val" ]] && HEIGHT="$_val"
    _val=$(cfg visualization combined_bootstrap_tree title);           [[ -n "$_val" ]] && TREE_TITLE="$_val"
    _val=$(cfg visualization combined_bootstrap_tree combined_style);  [[ -n "$_val" ]] && COMBINED_STYLE="$_val"
    _val=$(cfg visualization combined_bootstrap_tree combined_sep);    [[ -n "$_val" ]] && COMBINED_SEP="$_val"
    _val=$(cfg visualization combined_bootstrap_tree xlim_expand);    [[ -n "$_val" ]] && XLIM_EXPAND="$_val"
    unset _val
fi

# ======================== Sequence Type Detection ========================

detect_seq_type() {
    local stem_lower
    stem_lower="$(echo "$1" | tr '[:upper:]' '[:lower:]')"
    if [[ "$stem_lower" =~ (nucleotide|_nuc_|_nt_|_dna_) ]]; then
        echo "Nucleotide"
    elif [[ "$stem_lower" =~ (amino_acid|protein|_aa_) ]]; then
        echo "Protein"
    else
        echo "Unknown"
    fi
}

# ======================== Find IQ-TREE2 / RAxML Pairs ========================

log_info "Scanning for IQ-TREE2/RAxML tree pairs in: $TREE_DIR"

PAIR_IQTREE=()
PAIR_RAXML=()
PAIR_OUTPUT=()
PAIR_SEQTYPE=()

while IFS= read -r iq_tree_file; do
    iq_dir="$(dirname "$iq_tree_file")"        # …/IQTREE2
    genome_dir="$(dirname "$iq_dir")"          # …/<genome>
    iq_fname="$(basename "$iq_tree_file")"     # <stem>_IQTREE2.treefile
    stem="${iq_fname%_IQTREE2.treefile}"       # <stem>

    # Guard: must be inside an IQTREE2 subdirectory
    if [[ "$(basename "$iq_dir")" != "IQTREE2" ]]; then
        log_warn "Skipping unexpected path (not inside IQTREE2/): $iq_tree_file"
        continue
    fi

    raxml_dir="$genome_dir/RAXML"
    raxml_file="$raxml_dir/${stem}_RAXML.raxml.support"

    # Require .raxml.support — bootstrap values are mandatory for combined figure
    if [[ ! -f "$raxml_file" ]]; then
        log_warn "No .raxml.support for stem '$stem' — skipping (looked in $raxml_dir/)"
        continue
    fi

    # genome_dir is the parent of IQTREE2/ — for the new MSA-mirrored layout
    # this is .../<METHOD>_aligned/, not the genome folder. Use the relative
    # path (treedir-stripped) for a clearer log label.
    rel_subpath="${genome_dir#$TREE_DIR/}"
    seq_type="$(detect_seq_type "$stem")"
    # Place combined-bootstrap figures next to the trees they merge (parent of
    # IQTREE2/ and RAXML/). Works for both the MSA-mirrored layout
    # (.../<METHOD>_aligned/Combined/...) and the legacy flat layout
    # (.../<genome>/Combined/...).
    out_png="$genome_dir/Combined/${seq_type}/${stem}_combined_bootstrap.png"

    PAIR_IQTREE+=("$iq_tree_file")
    PAIR_RAXML+=("$raxml_file")
    PAIR_OUTPUT+=("$out_png")
    PAIR_SEQTYPE+=("$seq_type")

    log_info "Matched pair: [$rel_subpath] [$seq_type] $stem"
done < <(find "$TREE_DIR" -type f -name "*_IQTREE2.treefile" 2>/dev/null | sort)

n_pairs=${#PAIR_IQTREE[@]}

if (( n_pairs == 0 )); then
    log_warn "No IQ-TREE2/RAxML pairs found — skipping combined bootstrap tree"
    exit 0
fi

log_info "Found $n_pairs tree pair(s) for combined bootstrap rendering"

# ======================== Render Combined Figures ========================

wait_for_slot() {
    local limit="$1"
    while (( $(jobs -rp | wc -l) >= limit )); do
        sleep 0.5
    done
}

rendered=0
skipped=0
pids=()

for (( i=0; i<n_pairs; i++ )); do
    iq_tree_file="${PAIR_IQTREE[$i]}"
    raxml_file="${PAIR_RAXML[$i]}"
    out_png="${PAIR_OUTPUT[$i]}"
    seq_type="${PAIR_SEQTYPE[$i]}"

    if [[ -s "$out_png" && "$OVERWRITE" != "true" ]]; then
        log_info "Skipping (exists): $(basename "$out_png")"
        skipped=$(( skipped + 1 ))
        continue
    fi

    mkdir -p "$(dirname "$out_png")"

    # Derive substitution model from config for caption
    PHYLO_MODEL=""
    if [[ -n "$CONFIG_FILE" && -f "$CONFIG_FILE" ]]; then
        if [[ "$seq_type" == "Nucleotide" ]]; then
            _m=$(cfg phylogenetics iqtree2 nucleotide_model 2>/dev/null || true)
        else
            _m=$(cfg phylogenetics iqtree2 protein_model 2>/dev/null || true)
        fi
        [[ -n "${_m:-}" ]] && PHYLO_MODEL="$_m"
        unset _m
    fi

    wait_for_slot "$MAX_PARALLEL"

    (
        log_info "Rendering: $(basename "$iq_tree_file") + $(basename "$raxml_file") -> $(basename "$out_png")"
        rc=0
        Rscript "$RENDER_R" \
            --tree-iqtree          "$iq_tree_file" \
            --tree-raxml           "$raxml_file" \
            --output               "$out_png" \
            --layout               "$LAYOUT" \
            --open-angle           "$OPEN_ANGLE" \
            --dpi                  "$DPI" \
            --width                "$WIDTH" \
            --height               "$HEIGHT" \
            --root-outgroup        "$ROOT_OUTGROUP" \
            --outgroup-pattern     "$OUTGROUP_PATTERN" \
            --exclude-tips         "$EXCLUDE_TIPS" \
            --highlight-eggplant   "$HIGHLIGHT_EGGPLANT" \
            --show-bootstrap       "$SHOW_BOOTSTRAP" \
            --bootstrap-threshold  "$BOOTSTRAP_THRESHOLD" \
            --bootstrap-style      "$BOOTSTRAP_STYLE" \
            --bootstrap-label-size "$BOOTSTRAP_LABEL_SIZE" \
            --bootstrap-color      "$BOOTSTRAP_COLOR" \
            --node-color-high      "$NODE_COLOR_HIGH" \
            --node-color-medium    "$NODE_COLOR_MEDIUM" \
            --node-size-high       "$NODE_SIZE_HIGH" \
            --node-size-medium     "$NODE_SIZE_MEDIUM" \
            --tip-label-size       "$TIP_LABEL_SIZE" \
            --tip-label-offset     "$TIP_LABEL_OFFSET" \
            --tip-point-size       "$TIP_POINT_SIZE" \
            --branch-width         "$BRANCH_WIDTH" \
            --branch-color         "$BRANCH_COLOR" \
            --label-style          "$LABEL_STYLE" \
            --combined-style       "$COMBINED_STYLE" \
            --combined-sep         "$COMBINED_SEP" \
            --title                "$TREE_TITLE" \
            --treescale-fontsize   "$TREESCALE_FONTSIZE" \
            --treescale-offset     "$TREESCALE_OFFSET" \
            --treescale-color      "$TREESCALE_COLOR" \
            --color-smeldmp        "$COLOR_SMELDMP" \
            --color-haploid        "$COLOR_HAPLOID" \
            --color-ortholog       "$COLOR_ORTHOLOG" \
            --color-outgroup       "$COLOR_OUTGROUP" \
            --xlim-expand          "$XLIM_EXPAND" \
            --sequence-type        "$seq_type" \
            --phylo-model          "$PHYLO_MODEL" \
            2>&1 || rc=$?

        if (( rc != 0 )); then
            log_error "render_combined_bootstrap.R failed (exit code $rc)"
            exit 1
        fi

        if [[ -s "$out_png" ]]; then
            log_info "Saved: $out_png ($(du -sh "$out_png" | cut -f1))"
        else
            log_error "R script produced no output: $out_png"
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

log_info "Combined bootstrap summary: rendered=$rendered, skipped=$skipped, failed=$fail_count"

if (( fail_count > 0 )); then
    log_error "$fail_count combined bootstrap render(s) failed"
    exit 1
fi
