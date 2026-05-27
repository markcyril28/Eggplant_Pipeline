#!/bin/bash
# ============================================================================
# regenerate_motif_locations.sh
# ----------------------------------------------------------------------------
# Re-render the full-view motif-locations PNG (plot_motif_locations.py) from
# an existing meme.xml WITHOUT rerunning MEME / TOMTOM / FIMO.  Useful after
# tweaking plot_motif_locations.py (fonts, colors, layout) or the central
# color TOML when MEME outputs themselves are already on disk.
#
# Two invocation modes:
#
#   1) Explicit:
#        bash regenerate_motif_locations.sh \
#            --meme-xml III_RESULT/DMP/06_Motif_Analysis/.../meme.xml \
#            [--outdir DIR]   [--label NAME]
#
#   2) Auto-discovery (scans all meme.xml under a results root):
#        bash regenerate_motif_locations.sh \
#            --scan-root III_RESULT/DMP/06_Motif_Analysis
#
# When --outdir / --label are omitted, they are inferred from the meme.xml
# location using the standard pipeline layout:
#     <BASE>/04_MEME_Analysis/02_MEME/<alph>/<label>/meme.xml
#         -> outdir = <BASE>/04_MEME_Analysis/05_JPEG/<alph>
#         -> label  = <label>     (folder name containing meme.xml)
#
# Palette / background / domain coloring / font scale / DPI are read from
# 06_motif_analysisCONFIG.toml and config/colors_config/meme_motif_colors.toml
# so the output matches what the pipeline would produce on a full run.
# Per-call overrides are available via --palette, --bg, --font-scale, --dpi.
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
source "$SCRIPT_DIR/../../logging/logging_utils.sh"

PLOT_SCRIPT="$SCRIPT_DIR/plot_motif_locations.py"
TOML_PARSER="$PIPELINE_DIR/modules/utils/parse_toml.py"
MOTIF_COLORS_TOML="$PIPELINE_DIR/config/colors_config/meme_motif_colors.toml"
ROOT_CONFIG="$PIPELINE_DIR/06_motif_analysisCONFIG.toml"

MEME_XML=""
SCAN_ROOT=""
OUTDIR=""
LABEL=""
PALETTE_OVERRIDE=""
BG_OVERRIDE=""
FONT_SCALE_OVERRIDE=""
DPI_OVERRIDE=""
PHYLO_ORDER_OVERRIDE=""
NO_DOMAIN=false

usage() {
    cat <<EOF
Usage: $(basename "$0") (--meme-xml PATH | --scan-root DIR) [options]

Required (one of):
  --meme-xml PATH       Single meme.xml to render
  --scan-root DIR       Directory; auto-discovers every meme.xml underneath
                        and regenerates each one

Optional overrides (defaults read from 06_motif_analysisCONFIG.toml):
  --outdir DIR          Output directory (only with --meme-xml; inferred
                        from layout if omitted)
  --label NAME          Output filename label (only with --meme-xml;
                        inferred from meme.xml folder if omitted)
  --palette NAME        Palette version or comma-separated hex codes
  --bg COLOR            Background: dark | light | hex
  --font-scale NUM      Global font-scale multiplier (default: TOML 1.4)
  --dpi NUM             Output DPI (default: TOML 600)
  --phylo-order FILE    CLUSTAL alignment for row ordering
  --no-domain           Disable structural-domain coloring (palette only)
  -h, --help            Show this help
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --meme-xml)     MEME_XML="$2";              shift 2 ;;
        --scan-root)    SCAN_ROOT="$2";             shift 2 ;;
        --outdir)       OUTDIR="$2";                shift 2 ;;
        --label)        LABEL="$2";                 shift 2 ;;
        --palette)      PALETTE_OVERRIDE="$2";      shift 2 ;;
        --bg)           BG_OVERRIDE="$2";           shift 2 ;;
        --font-scale)   FONT_SCALE_OVERRIDE="$2";   shift 2 ;;
        --dpi)          DPI_OVERRIDE="$2";          shift 2 ;;
        --phylo-order)  PHYLO_ORDER_OVERRIDE="$2";  shift 2 ;;
        --no-domain)    NO_DOMAIN=true;             shift   ;;
        -h|--help)      usage ;;
        *) log_error "Unknown argument: $1"; usage ;;
    esac
done

[[ -z "$MEME_XML" && -z "$SCAN_ROOT" ]] && { log_error "--meme-xml or --scan-root required"; usage; }
[[ -n "$MEME_XML" && -n "$SCAN_ROOT" ]] && { log_error "--meme-xml and --scan-root are mutually exclusive"; exit 1; }
[[ ! -f "$PLOT_SCRIPT" ]] && { log_error "plot_motif_locations.py not found: $PLOT_SCRIPT"; exit 1; }

if ! python3 -c "import matplotlib" &>/dev/null; then
    log_error "matplotlib not available — activate the 'egg' conda env first"
    log_error "  conda activate egg"
    exit 1
fi

# ---------------------------------------------------------------------------
# Resolve defaults from TOML (06_motif_analysisCONFIG.toml -> [meme])
# ---------------------------------------------------------------------------
toml_get() {  # toml_get <file> <section> <key> [default]
    local _file="$1" _section="$2" _key="$3" _default="${4:-}"
    local _val
    _val=$(python3 "$TOML_PARSER" "$_file" "$_section" "$_key" 2>/dev/null) || true
    [[ -n "$_val" ]] && echo "$_val" || echo "$_default"
}

PALETTE_NAME=$(toml_get "$ROOT_CONFIG" meme motif_location_palette "wong")
BG_NAME=$(toml_get      "$ROOT_CONFIG" meme motif_location_bg      "")
FONT_SCALE=$(toml_get   "$ROOT_CONFIG" meme motif_location_font_scale "1.4")
JPEG_DPI=$(toml_get     "$ROOT_CONFIG" meme jpeg_dpi               "600")

[[ -n "$PALETTE_OVERRIDE"    ]] && PALETTE_NAME="$PALETTE_OVERRIDE"
[[ -n "$BG_OVERRIDE"         ]] && BG_NAME="$BG_OVERRIDE"
[[ -n "$FONT_SCALE_OVERRIDE" ]] && FONT_SCALE="$FONT_SCALE_OVERRIDE"
[[ -n "$DPI_OVERRIDE"        ]] && JPEG_DPI="$DPI_OVERRIDE"

# Resolve palette: named version -> hex list, or raw hex passthrough.
RESOLVED_PALETTE=""
case "$PALETTE_NAME" in
    wong|warm|cool|pastel|high_contrast|dark|dmp_interaction)
        if [[ -f "$MOTIF_COLORS_TOML" ]]; then
            _raw=$(python3 "$TOML_PARSER" "$MOTIF_COLORS_TOML" "$PALETTE_NAME" colors 2>/dev/null) || true
            [[ -n "$_raw" ]] && RESOLVED_PALETTE="${_raw//$'\n'/,}"
            unset _raw
        fi
        ;;
    "") ;;
    *)  RESOLVED_PALETTE="$PALETTE_NAME" ;;
esac

# Resolve background colour shorthand.
RESOLVED_BG=""
case "$BG_NAME" in
    dark)  RESOLVED_BG="#111111" ;;
    light) RESOLVED_BG="white"   ;;
    "")    RESOLVED_BG=""        ;;
    *)     RESOLVED_BG="$BG_NAME" ;;
esac

# Resolve domain colors / labels / position-based mapping (dmp_interaction).
DOMAIN_COLORS=""
DOMAIN_LABELS=""
REF_SEQ_IDS=""
DOMAIN_RANGES_JSON=""
if ! $NO_DOMAIN && [[ "$PALETTE_NAME" == "dmp_interaction" ]] && [[ -f "$MOTIF_COLORS_TOML" ]]; then
    DOMAIN_COLORS=$(python3      "$TOML_PARSER" "$MOTIF_COLORS_TOML" dmp_interaction domain_colors_json   2>/dev/null) || true
    DOMAIN_LABELS=$(python3      "$TOML_PARSER" "$MOTIF_COLORS_TOML" dmp_interaction domain_labels_json   2>/dev/null) || true
    DOMAIN_RANGES_JSON=$(python3 "$TOML_PARSER" "$MOTIF_COLORS_TOML" dmp_interaction domain_ranges_json   2>/dev/null) || true
    _rids_raw=$(python3          "$TOML_PARSER" "$MOTIF_COLORS_TOML" dmp_interaction reference_seq_ids    2>/dev/null) || true
    [[ -n "$_rids_raw" ]] && REF_SEQ_IDS="${_rids_raw//$'\n'/,}"
    unset _rids_raw
fi

# ---------------------------------------------------------------------------
# Derive outdir/label from meme.xml location (pipeline-standard layout)
# ---------------------------------------------------------------------------
infer_outdir_label() {  # infer_outdir_label <meme_xml>; sets INFER_OUTDIR INFER_LABEL INFER_ALPH
    local xml="$1"
    local meme_set_dir alph_dir meme_dir alph base_dir
    meme_set_dir="$(dirname "$xml")"          # .../02_MEME/<alph>/<label>
    INFER_LABEL="$(basename "$meme_set_dir")"
    alph_dir="$(dirname "$meme_set_dir")"     # .../02_MEME/<alph>
    INFER_ALPH="$(basename "$alph_dir")"      # amino_acid | nucleotide
    meme_dir="$(dirname "$alph_dir")"         # .../02_MEME
    base_dir="$(dirname "$meme_dir")"         # .../04_MEME_Analysis
    INFER_OUTDIR="$base_dir/05_JPEG/$INFER_ALPH"
}

# ---------------------------------------------------------------------------
# Render a single meme.xml
# ---------------------------------------------------------------------------
render_one() {  # render_one <meme_xml> [outdir] [label]
    local xml="$1" out="${2:-}" lbl="${3:-}"
    [[ ! -f "$xml" ]] && { log_warn "meme.xml not found: $xml — skipping"; return 1; }

    local INFER_OUTDIR INFER_LABEL INFER_ALPH
    infer_outdir_label "$xml"
    [[ -z "$out" ]] && out="$INFER_OUTDIR"
    [[ -z "$lbl" ]] && lbl="$INFER_LABEL"

    mkdir -p "$out"

    # Domain coloring only makes sense for amino-acid datasets.
    local apply_domain=true
    if [[ "$INFER_ALPH" == "nucleotide" ]]; then
        apply_domain=false
    fi

    local -a cmd=(
        python3 "$PLOT_SCRIPT"
        --meme-xml   "$xml"
        --outdir     "$out"
        --label      "$lbl"
        --dpi        "$JPEG_DPI"
        --font-scale "$FONT_SCALE"
    )
    [[ -n "$RESOLVED_PALETTE"     ]] && cmd+=(--palette            "$RESOLVED_PALETTE")
    [[ -n "$RESOLVED_BG"          ]] && cmd+=(--bg-color           "$RESOLVED_BG")
    [[ -n "$PHYLO_ORDER_OVERRIDE" ]] && cmd+=(--phylo-order        "$PHYLO_ORDER_OVERRIDE")
    if $apply_domain; then
        [[ -n "$DOMAIN_COLORS"      ]] && cmd+=(--domain-colors      "$DOMAIN_COLORS")
        [[ -n "$DOMAIN_LABELS"      ]] && cmd+=(--domain-labels      "$DOMAIN_LABELS")
        [[ -n "$REF_SEQ_IDS"        ]] && cmd+=(--ref-seq-ids        "$REF_SEQ_IDS")
        [[ -n "$DOMAIN_RANGES_JSON" ]] && cmd+=(--domain-ranges-json "$DOMAIN_RANGES_JSON")
    fi

    log_info "Render: $lbl  ($INFER_ALPH)"
    log_info "  meme.xml -> $xml"
    log_info "  outdir   -> $out"
    "${cmd[@]}" 2>&1 | tee "$out/${lbl}_fullview.log" || {
        log_warn "Plot failed for $lbl — see $out/${lbl}_fullview.log"
        return 1
    }

    local png="$out/${lbl}_motif_locations_full.png"
    if [[ -f "$png" ]]; then
        local sz
        sz=$(du -k "$png" 2>/dev/null | awk '{print $1}')
        log_info "  Wrote: $png (${sz} KB)"
    else
        log_warn "  Expected PNG not produced: $png"
        return 1
    fi
}

# ---------------------------------------------------------------------------
# Main: dispatch single-xml vs scan-root
# ---------------------------------------------------------------------------
fail_count=0
ok_count=0

if [[ -n "$MEME_XML" ]]; then
    render_one "$MEME_XML" "$OUTDIR" "$LABEL" && ((ok_count++)) || ((fail_count++))
else
    [[ ! -d "$SCAN_ROOT" ]] && { log_error "Scan root not a directory: $SCAN_ROOT"; exit 1; }
    log_info "Scanning $SCAN_ROOT for meme.xml ..."
    mapfile -t xmls < <(find "$SCAN_ROOT" -type f -name meme.xml | sort)
    log_info "Found ${#xmls[@]} meme.xml file(s)"
    [[ ${#xmls[@]} -eq 0 ]] && { log_warn "No meme.xml under $SCAN_ROOT — nothing to do"; exit 0; }
    for xml in "${xmls[@]}"; do
        render_one "$xml" "" "" && ((ok_count++)) || ((fail_count++))
    done
fi

log_info "Done: $ok_count succeeded, $fail_count failed"
(( fail_count == 0 ))
