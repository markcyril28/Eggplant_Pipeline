#!/bin/bash
# ============================================================================
# Module: Phylogenetic Tree Comparison (Tanglegram)
# ============================================================================
# Finds matched tree pairs across IQ-TREE2, RAxML, and MEGA-CC for the same
# alignment stem, then produces a tanglegram comparison figure for each pair.
#
# Pair-matching logic:
#   For every IQ-TREE2 .treefile under <treedir>/<...subpath...>/IQTREE2/,
#   look for the sibling RAxML output in <subpath>/RAXML/ (.raxml.support
#   preferred, .raxml.bestTree as fallback). One IQ vs RAX pair per stem.
#   Additionally, for every MEGA-CC .nwk under <subpath>/MEGA_CC/, pair it
#   with any sibling IQ-TREE2 and/or RAxML tree to emit MEGA vs IQ and
#   MEGA vs RAX comparisons.
#
# Output: <subpath>/Comparison/{Nucleotide|Protein}/<stem>_{IQ_vs_RAX|MEGA_vs_IQ|MEGA_vs_RAX}.png
#         (written next to the tree pair, not under --outdir).
#
# Usage (standalone):
#   bash compare_trees.sh \
#       --treedir  <phylo_output_dir> \
#       --outdir   <comparison_output_dir> \
#       --config   <merged_config.toml> \
#       [--threads N] [--overwrite true|false]
# ============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../utils/logging.sh"

COMPARE_R="$SCRIPT_DIR/compare_trees.R"

# ======================== Defaults ========================

TREE_DIR=""
OUTPUT_DIR=""
CONFIG_FILE=""
MAX_PARALLEL=4
OVERWRITE="true"

# Visualization defaults — aligned with visualize_tree.sh for consistent figures.
# Tanglegram-specific overrides (width, height) come from [visualization.compare_trees] TOML.
# Per-pair labels (IQ-TREE2 / RAxML / MEGA-CC) are assigned during the pair scan.
DPI=600
WIDTH=14
HEIGHT=10
ROOT_OUTGROUP="false"
OUTGROUP_PATTERN="Outgroup"
EXCLUDE_TIPS=""
HIGHLIGHT_EGGPLANT="true"
EGGPLANT_PATTERN="^SMEL5_|^Smel_|^SMELG|^SmelDMP"
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
TIP_POINT_SIZE=2.5
BRANCH_WIDTH=0.8
BRANCH_COLOR="grey20"
LABEL_STYLE="replace"
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
        --treedir)     TREE_DIR="$2";    shift 2 ;;
        --outdir)      OUTPUT_DIR="$2";  shift 2 ;;
        --config)      CONFIG_FILE="$2"; shift 2 ;;
        --threads)     MAX_PARALLEL="$2"; shift 2 ;;
        --overwrite)   OVERWRITE="$2";   shift 2 ;;
        *) log_error "Unknown argument: $1"; exit 1 ;;
    esac
done

[[ -z "$TREE_DIR" ]] && { log_error "Missing --treedir"; exit 1; }
[[ -d "$TREE_DIR" ]] || { log_error "Tree directory not found: $TREE_DIR"; exit 1; }
[[ ! -f "$COMPARE_R" ]] && { log_error "R script not found: $COMPARE_R"; exit 1; }

# Default output dir = treedir
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

    _val=$(cfg visualization dpi);                  [[ -n "$_val" ]] && DPI="$_val"
    _val=$(cfg visualization width);                [[ -n "$_val" ]] && WIDTH="$_val"
    _val=$(cfg visualization height);               [[ -n "$_val" ]] && HEIGHT="$_val"
    _val=$(cfg visualization root_outgroup);        [[ -n "$_val" ]] && ROOT_OUTGROUP="$_val"
    _val=$(cfg visualization outgroup_pattern);     [[ -n "$_val" ]] && OUTGROUP_PATTERN="$_val"
    _val=$(cfg visualization exclude_tips);         [[ -n "$_val" ]] && EXCLUDE_TIPS="$_val"
    _val=$(cfg visualization highlight_eggplant);   [[ -n "$_val" ]] && HIGHLIGHT_EGGPLANT="$_val"
    _val=$(cfg visualization eggplant_pattern);     [[ -n "$_val" ]] && EGGPLANT_PATTERN="$_val"
    _val=$(cfg visualization show_bootstrap);       [[ -n "$_val" ]] && SHOW_BOOTSTRAP="$_val"
    _val=$(cfg visualization bootstrap_threshold);  [[ -n "$_val" ]] && BOOTSTRAP_THRESHOLD="$_val"
    _val=$(cfg visualization bootstrap_style);      [[ -n "$_val" ]] && BOOTSTRAP_STYLE="$_val"
    _val=$(cfg visualization bootstrap_label_size); [[ -n "$_val" ]] && BOOTSTRAP_LABEL_SIZE="$_val"
    _val=$(cfg visualization bootstrap_color);      [[ -n "$_val" ]] && BOOTSTRAP_COLOR="$_val"
    _val=$(cfg visualization node_color_high);      [[ -n "$_val" ]] && NODE_COLOR_HIGH="$_val"
    _val=$(cfg visualization node_color_medium);    [[ -n "$_val" ]] && NODE_COLOR_MEDIUM="$_val"
    _val=$(cfg visualization node_size_high);       [[ -n "$_val" ]] && NODE_SIZE_HIGH="$_val"
    _val=$(cfg visualization node_size_medium);     [[ -n "$_val" ]] && NODE_SIZE_MEDIUM="$_val"
    _val=$(cfg visualization tip_label_size);       [[ -n "$_val" ]] && TIP_LABEL_SIZE="$_val"
    _val=$(cfg visualization tip_label_offset);     [[ -n "$_val" ]] && TIP_LABEL_OFFSET="$_val"
    _val=$(cfg visualization tip_point_size);       [[ -n "$_val" ]] && TIP_POINT_SIZE="$_val"
    _val=$(cfg visualization branch_width);         [[ -n "$_val" ]] && BRANCH_WIDTH="$_val"
    _val=$(cfg visualization branch_color);         [[ -n "$_val" ]] && BRANCH_COLOR="$_val"
    _val=$(cfg visualization label_style);          [[ -n "$_val" ]] && LABEL_STYLE="$_val"
    _val=$(cfg visualization treescale_fontsize);   [[ -n "$_val" ]] && TREESCALE_FONTSIZE="$_val"
    _val=$(cfg visualization treescale_offset);     [[ -n "$_val" ]] && TREESCALE_OFFSET="$_val"
    _val=$(cfg visualization treescale_color);      [[ -n "$_val" ]] && TREESCALE_COLOR="$_val"
    _val=$(cfg visualization color_smeldmp);        [[ -n "$_val" ]] && COLOR_SMELDMP="$_val"
    _val=$(cfg visualization color_haploid);        [[ -n "$_val" ]] && COLOR_HAPLOID="$_val"
    _val=$(cfg visualization color_ortholog);       [[ -n "$_val" ]] && COLOR_ORTHOLOG="$_val"
    _val=$(cfg visualization color_outgroup);       [[ -n "$_val" ]] && COLOR_OUTGROUP="$_val"

    # Compare-trees-specific overrides (tanglegram needs wider canvas)
    _val=$(cfg visualization compare_trees dpi);    [[ -n "$_val" ]] && DPI="$_val"
    _val=$(cfg visualization compare_trees width);  [[ -n "$_val" ]] && WIDTH="$_val"
    _val=$(cfg visualization compare_trees height); [[ -n "$_val" ]] && HEIGHT="$_val"
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

# Resolve a sibling MEGA_CC tree under <genome_dir>/MEGA_CC/, preferring the
# ML best tree over the bootstrap consensus variant. Empty if none.
find_sibling_mega_tree() {
    local genome_dir="$1"
    local mega_dir="$genome_dir/MEGA_CC"
    [[ -d "$mega_dir" ]] || return 0
    local best
    best=$(find "$mega_dir" -maxdepth 1 -type f -name "*.nwk" ! -name "*_consensus.nwk" 2>/dev/null | sort | head -n 1)
    if [[ -z "$best" ]]; then
        best=$(find "$mega_dir" -maxdepth 1 -type f -name "*_consensus.nwk" 2>/dev/null | sort | head -n 1)
    fi
    [[ -n "$best" ]] && echo "$best" || true
}

# ======================== Find Tree Pairs ========================

log_info "Scanning for tree pairs (IQ vs RAX, and MEGA vs IQ/RAX where present) in: $TREE_DIR"

PAIR_TREE1=()
PAIR_TREE2=()
PAIR_LABEL1=()
PAIR_LABEL2=()
PAIR_OUTPUT=()
PAIR_SEQTYPE=()

# --- IQ-TREE2 vs RAxML (the canonical comparison) ---
while IFS= read -r iq_tree; do
    iq_dir="$(dirname "$iq_tree")"          # .../IQTREE2
    genome_dir="$(dirname "$iq_dir")"       # .../<genome>
    iq_fname="$(basename "$iq_tree")"       # <stem>_IQTREE2.treefile
    stem="${iq_fname%_IQTREE2.treefile}"    # <stem>

    if [[ "$(basename "$iq_dir")" != "IQTREE2" ]]; then
        log_warn "Skipping unexpected path (not inside IQTREE2/): $iq_tree"
        continue
    fi

    raxml_dir="$genome_dir/RAXML"
    raxml_tree=""
    if [[ -f "$raxml_dir/${stem}_RAXML.raxml.support" ]]; then
        raxml_tree="$raxml_dir/${stem}_RAXML.raxml.support"
    elif [[ -f "$raxml_dir/${stem}_RAXML.raxml.bestTree" ]]; then
        raxml_tree="$raxml_dir/${stem}_RAXML.raxml.bestTree"
        log_warn "No .raxml.support for stem '$stem' — using .raxml.bestTree (no bootstrap values)"
    fi

    if [[ -z "$raxml_tree" ]]; then
        log_warn "No RAxML counterpart found for: $iq_tree (looked in $raxml_dir/)"
        continue
    fi

    rel_subpath="${genome_dir#$TREE_DIR/}"
    seq_type="$(detect_seq_type "$stem")"
    out_png="$genome_dir/Comparison/${seq_type}/${stem}_IQ_vs_RAX.png"

    PAIR_TREE1+=("$iq_tree")
    PAIR_TREE2+=("$raxml_tree")
    PAIR_LABEL1+=("IQ-TREE2")
    PAIR_LABEL2+=("RAxML")
    PAIR_OUTPUT+=("$out_png")
    PAIR_SEQTYPE+=("$seq_type")

    log_info "Matched pair [IQ vs RAX]: [$rel_subpath] [$seq_type] $stem"
done < <(find "$TREE_DIR" -type f -name "*_IQTREE2.treefile" 2>/dev/null | sort)

# --- MEGA-CC vs IQ-TREE2 and MEGA-CC vs RAxML (when MEGA output exists) ---
while IFS= read -r mega_tree; do
    mega_dir="$(dirname "$mega_tree")"      # .../MEGA_CC
    [[ "$(basename "$mega_dir")" == "MEGA_CC" ]] || continue
    # Skip MEGA consensus variant when sibling ML best tree exists
    mega_fname="$(basename "$mega_tree")"
    if [[ "$mega_fname" == *_consensus.nwk ]]; then
        sibling="${mega_tree%_consensus.nwk}.nwk"
        [[ -s "$sibling" ]] && continue
    fi
    genome_dir="$(dirname "$mega_dir")"
    mega_stem="${mega_fname%.nwk}"
    rel_subpath="${genome_dir#$TREE_DIR/}"
    seq_type="$(detect_seq_type "$mega_stem")"

    # MEGA vs IQ
    iq_dir="$genome_dir/IQTREE2"
    iq_candidate=$(find "$iq_dir" -maxdepth 1 -type f -name "*_IQTREE2.treefile" 2>/dev/null | sort | head -n 1)
    if [[ -n "$iq_candidate" ]]; then
        out_png="$genome_dir/Comparison/${seq_type}/${mega_stem}_MEGA_vs_IQ.png"
        PAIR_TREE1+=("$mega_tree")
        PAIR_TREE2+=("$iq_candidate")
        PAIR_LABEL1+=("MEGA-CC")
        PAIR_LABEL2+=("IQ-TREE2")
        PAIR_OUTPUT+=("$out_png")
        PAIR_SEQTYPE+=("$seq_type")
        log_info "Matched pair [MEGA vs IQ]: [$rel_subpath] [$seq_type] $mega_stem"
    fi

    # MEGA vs RAX (prefer .raxml.support over .raxml.bestTree)
    raxml_dir="$genome_dir/RAXML"
    rax_candidate=$(find "$raxml_dir" -maxdepth 1 -type f -name "*.raxml.support" 2>/dev/null | sort | head -n 1)
    if [[ -z "$rax_candidate" ]]; then
        rax_candidate=$(find "$raxml_dir" -maxdepth 1 -type f -name "*.raxml.bestTree" 2>/dev/null | sort | head -n 1)
    fi
    if [[ -n "$rax_candidate" ]]; then
        out_png="$genome_dir/Comparison/${seq_type}/${mega_stem}_MEGA_vs_RAX.png"
        PAIR_TREE1+=("$mega_tree")
        PAIR_TREE2+=("$rax_candidate")
        PAIR_LABEL1+=("MEGA-CC")
        PAIR_LABEL2+=("RAxML")
        PAIR_OUTPUT+=("$out_png")
        PAIR_SEQTYPE+=("$seq_type")
        log_info "Matched pair [MEGA vs RAX]: [$rel_subpath] [$seq_type] $mega_stem"
    fi
done < <(find "$TREE_DIR" -type f -path "*/MEGA_CC/*.nwk" 2>/dev/null | sort)

n_pairs=${#PAIR_TREE1[@]}

if (( n_pairs == 0 )); then
    log_warn "No tree pairs found — skipping tree comparison"
    exit 0
fi

log_info "Found $n_pairs tree pair(s) to compare"

# ======================== Render Comparisons ========================

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
    tree1="${PAIR_TREE1[$i]}"
    tree2="${PAIR_TREE2[$i]}"
    label1="${PAIR_LABEL1[$i]}"
    label2="${PAIR_LABEL2[$i]}"
    out_png="${PAIR_OUTPUT[$i]}"
    seq_type="${PAIR_SEQTYPE[$i]}"

    if [[ -s "$out_png" && "$OVERWRITE" != "true" ]]; then
        log_info "Skipping (exists): $(basename "$out_png")"
        skipped=$(( skipped + 1 ))
        continue
    fi

    mkdir -p "$(dirname "$out_png")"

    wait_for_slot "$MAX_PARALLEL"

    (
        log_info "Comparing: $label1 ($(basename "$tree1")) vs $label2 ($(basename "$tree2"))"
        rc=0
        Rscript "$COMPARE_R" \
            --tree1                 "$tree1" \
            --tree2                 "$tree2" \
            --label1                "$label1" \
            --label2                "$label2" \
            --output                "$out_png" \
            --dpi                   "$DPI" \
            --width                 "$WIDTH" \
            --height                "$HEIGHT" \
            --root-outgroup         "$ROOT_OUTGROUP" \
            --outgroup-pattern      "$OUTGROUP_PATTERN" \
            --exclude-tips          "$EXCLUDE_TIPS" \
            --highlight-eggplant    "$HIGHLIGHT_EGGPLANT" \
            --eggplant-pattern      "$EGGPLANT_PATTERN" \
            --show-bootstrap        "$SHOW_BOOTSTRAP" \
            --bootstrap-threshold   "$BOOTSTRAP_THRESHOLD" \
            --bootstrap-style       "$BOOTSTRAP_STYLE" \
            --bootstrap-label-size  "$BOOTSTRAP_LABEL_SIZE" \
            --bootstrap-color       "$BOOTSTRAP_COLOR" \
            --node-color-high       "$NODE_COLOR_HIGH" \
            --node-color-medium     "$NODE_COLOR_MEDIUM" \
            --node-size-high        "$NODE_SIZE_HIGH" \
            --node-size-medium      "$NODE_SIZE_MEDIUM" \
            --tip-label-size        "$TIP_LABEL_SIZE" \
            --tip-label-offset      "$TIP_LABEL_OFFSET" \
            --tip-point-size        "$TIP_POINT_SIZE" \
            --branch-width          "$BRANCH_WIDTH" \
            --branch-color          "$BRANCH_COLOR" \
            --label-style           "$LABEL_STYLE" \
            --treescale-fontsize    "$TREESCALE_FONTSIZE" \
            --treescale-offset      "$TREESCALE_OFFSET" \
            --treescale-color       "$TREESCALE_COLOR" \
            --color-smeldmp         "$COLOR_SMELDMP" \
            --color-haploid         "$COLOR_HAPLOID" \
            --color-ortholog        "$COLOR_ORTHOLOG" \
            --color-outgroup        "$COLOR_OUTGROUP" \
            --seq-type              "$seq_type" \
            2>&1 || rc=$?

        if (( rc != 0 )); then
            log_error "compare_trees.R failed for pair: $tree1 vs $tree2 (exit code $rc)"
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

log_info "Comparison summary: rendered=$rendered, skipped=$skipped, failed=$fail_count"

if (( fail_count > 0 )); then
    log_error "$fail_count comparison(s) failed"
    exit 1
fi
