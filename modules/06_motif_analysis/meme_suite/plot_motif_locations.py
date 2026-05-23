#!/usr/bin/env python3
"""
Full-view motif location diagram - all sequences, all motifs.

Reads a MEME XML file and draws a publication-ready PNG (and optionally SVG)
in the MEME web-suite default layout:

    Header:    Name | p-value | Motif Locations
    Rows:      one per input sequence, with the sequence drawn as a thin
               horizontal line and motif hits as coloured rectangles.
    Footer 1:  Structural-domain legend (when --domain-colors / --ref-seq-ids
               are set, motifs are grouped by structural domain).
    Footer 2:  Motif Consensus table - Motif # | Symbol | Motif Consensus
               (lifted from the <motif name="..."> attribute in meme.xml).

Unlike MEME's default SVG, every input sequence is shown regardless of
dataset size.

Usage:
    python3 plot_motif_locations.py \
        --meme-xml  /path/to/meme.xml \
        --outdir    /path/to/output_dir \
        --label     my_dataset \
        [--dpi      200] \
        [--svg] \
        [--font-scale 1.4] \
        [--phylo-order  /path/to/alignment.aln]
"""

import argparse
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# CLI -------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--meme-xml",  required=True,  help="Path to meme.xml")
    p.add_argument("--outdir",    required=True,  help="Output directory")
    p.add_argument("--label",     default="",     help="Dataset label (used in title)")
    p.add_argument("--dpi",       type=int, default=200, help="Output DPI (default 200)")
    p.add_argument("--svg",       action="store_true",   help="Also write SVG")
    p.add_argument("--font-scale", type=float, default=1.4,
                   help="Global font scaling factor (default 1.4 for readability)")
    p.add_argument("--palette",   default="",
                   help="Motif block colors: comma-separated hex codes, "
                        "or JSON array, e.g. '#E69F00,#56B4E9,...'")
    p.add_argument("--domain-colors", default="",
                   help="JSON object mapping motif number (str) to hex color, "
                        'e.g. \'{"1":"#FFA600","2":"#D9A621",...}\'. '
                        "When set, motifs sharing a color are grouped in the legend "
                        "using --domain-labels.")
    p.add_argument("--domain-labels", default="",
                   help="JSON object mapping hex color to legend label, "
                        'e.g. \'{"#FFA600":"TM Helices","#D9A621":"Amphipathic Helix"}\'. '
                        "Used only together with --domain-colors.")
    p.add_argument("--bar-color", default="",
                   help="Hex color for the sequence backbone line (default #888888).")
    p.add_argument("--bar-edge-color", default="",
                   help="Hex color reserved for backward compat (unused in MEME-default style).")
    p.add_argument("--bg-color", default="",
                   help="Figure background hex color, e.g. '#111111'. "
                        "Dark backgrounds (luminance < 0.4) automatically switch "
                        "all text and decorations to light variants.")
    p.add_argument("--phylo-order", default="",
                   help="Path to a CLUSTAL alignment file whose sequence order "
                        "defines the top-to-bottom order of rows in the figure. "
                        "Sequences not present in the file are appended at the "
                        "bottom in their original MEME order.")
    p.add_argument("--ref-seq-ids", default="",
                   help="Comma-separated reference sequence IDs (sequence names "
                        "as they appear in meme.xml). When set together with "
                        "--domain-ranges-json, motif colors are derived from "
                        "where each motif best matches on this reference, "
                        "overriding any --domain-colors / --domain-labels.")
    p.add_argument("--domain-ranges-json", default="",
                   help="JSON describing protein domain ranges on the reference "
                        "sequence: '{\"ranges\":[{\"start\":N,\"end\":M,"
                        "\"color\":\"#hex\",\"label\":\"...\"},...],"
                        "\"default_color\":\"#hex\",\"default_label\":\"...\"}'. "
                        "Positions are 0-based and inclusive on the reference.")
    return p.parse_args()


# parse XML -------------------------------------------------------------------
def load_meme_xml(xml_path):
    """Parse meme.xml into seq_map, motif_map, sites_by_seq, seq_pvalues.

    seq_map[sid]   = {name, length}
    motif_map[mid] = {num, width, evalue, consensus}
    sites_by_seq[sid] = list of {motif_id, position, pvalue}
    seq_pvalues[sid]  = float (combined p-value from <scanned_sites pvalue=...>)
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    seq_map = {}
    for s in root.findall(".//sequence"):
        seq_map[s.get("id")] = {
            "name":   s.get("name"),
            "length": int(s.get("length")),
        }

    motif_map = {}
    for m in root.findall(".//motif"):
        mid  = m.get("id")
        name = m.get("name", mid)
        try:
            num = int(name.split("-")[-1])
        except ValueError:
            num = int(mid.split("_")[-1])
        # Consensus: prefer <regular_expression>, fall back to MEME name attr.
        re_el = m.find("regular_expression")
        if re_el is not None and (re_el.text or "").strip():
            consensus = re_el.text.strip()
        else:
            consensus = name
        motif_map[mid] = {
            "num":       num,
            "width":     int(m.get("width")),
            "evalue":    m.get("e_value", "?"),
            "consensus": consensus,
        }

    sites_by_seq = {sid: [] for sid in seq_map}
    seq_pvalues  = {sid: None for sid in seq_map}
    for ss in root.findall(".//scanned_sites"):
        sid = ss.get("sequence_id")
        if sid in seq_pvalues:
            try:
                seq_pvalues[sid] = float(ss.get("pvalue"))
            except (TypeError, ValueError):
                seq_pvalues[sid] = None
        for site in ss.findall("scanned_site"):
            sites_by_seq[sid].append({
                "motif_id": site.get("motif_id"),
                "position": int(site.get("position")),
                "pvalue":   float(site.get("pvalue")),
            })

    return seq_map, motif_map, sites_by_seq, seq_pvalues


# phylo ordering --------------------------------------------------------------
def parse_phylo_order(aln_path: str) -> list:
    """Return sequence names in the order they first appear in a CLUSTAL file."""
    seen  = set()
    order = []
    try:
        with open(aln_path, encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.rstrip()
                if not line:
                    continue
                if line.startswith("CLUSTAL") or set(line.strip()).issubset(set("*:. ")):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                name = parts[0]
                if name not in seen:
                    seen.add(name)
                    order.append(name)
    except OSError as exc:
        print(f"WARNING: Cannot read phylo-order file {aln_path}: {exc}",
              file=sys.stderr)
    return order


def apply_phylo_order(seq_map: dict, phylo_names: list) -> dict:
    """Reorder seq_map (OrderedDict semantics) by phylo_names."""
    name_to_sid = {v["name"]: k for k, v in seq_map.items()}

    ordered_sids = []
    used = set()
    for pname in phylo_names:
        if pname in name_to_sid:
            sid = name_to_sid[pname]
            if sid not in used:
                ordered_sids.append(sid)
                used.add(sid)
            continue
        for seq_name, sid in name_to_sid.items():
            if pname in seq_name or seq_name in pname:
                if sid not in used:
                    ordered_sids.append(sid)
                    used.add(sid)
                break

    for sid in seq_map:
        if sid not in used:
            ordered_sids.append(sid)

    return {sid: seq_map[sid] for sid in ordered_sids}


# colour palette (colourblind-safe, up to 10 motifs) --------------------------
DEFAULT_PALETTE = [
    "#E69F00",  # 1  orange
    "#56B4E9",  # 2  sky blue
    "#009E73",  # 3  green
    "#F0E442",  # 4  yellow
    "#0072B2",  # 5  blue
    "#D55E00",  # 6  vermillion
    "#CC79A7",  # 7  pink
    "#999999",  # 8  grey
    "#44AA99",  # 9  teal
    "#882255",  # 10 wine
]


def resolve_palette(palette_arg: str) -> list:
    """Return a list of hex color strings from --palette argument."""
    if not palette_arg:
        return DEFAULT_PALETTE

    palette_arg = palette_arg.strip()

    if palette_arg.startswith("["):
        try:
            colors = json.loads(palette_arg)
            if isinstance(colors, list) and all(isinstance(c, str) for c in colors):
                return colors
        except (json.JSONDecodeError, TypeError):
            pass

    if "," in palette_arg:
        return [c.strip() for c in palette_arg.split(",") if c.strip()]

    if "\n" in palette_arg:
        return [c.strip() for c in palette_arg.split("\n") if c.strip()]

    return [palette_arg]


def derive_domains_from_reference(xml_path: str, ref_seq_ids: list,
                                  ranges_spec: dict) -> tuple:
    """Map each motif to a structural domain by where it best matches on a
    reference sequence in the meme.xml.

    Returns (domain_colors, domain_labels, ref_used) or (None, None, None)
    when no candidate reference is in the file.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    name_to_sid = {s.get("name"): s.get("id") for s in root.findall(".//sequence")}
    ref_sid, ref_name = None, None
    for cand in ref_seq_ids:
        if cand in name_to_sid:
            ref_sid = name_to_sid[cand]
            ref_name = cand
            break
    if ref_sid is None:
        return None, None, None

    best_pos = {}  # motif_id (str) -> (position, pvalue)
    for ss in root.findall(".//scanned_sites"):
        if ss.get("sequence_id") != ref_sid:
            continue
        for site in ss.findall("scanned_site"):
            mid = site.get("motif_id")
            pos = int(site.get("position"))
            try:
                pv = float(site.get("pvalue"))
            except (TypeError, ValueError):
                pv = 1.0
            if mid not in best_pos or pv < best_pos[mid][1]:
                best_pos[mid] = (pos, pv)

    ranges = ranges_spec.get("ranges", [])
    default_color = ranges_spec.get("default_color", "#B3B3B3")
    default_label = ranges_spec.get("default_label", "Loops")

    def domain_at(pos: int):
        for r in ranges:
            if r["start"] <= pos <= r["end"]:
                return r["color"], r["label"]
        return default_color, default_label

    domain_colors = {}
    domain_labels = {}
    for m in root.findall(".//motif"):
        mid = m.get("id")
        name = m.get("name", mid)
        try:
            num = int(name.split("-")[-1])
        except ValueError:
            try:
                num = int(mid.split("_")[-1])
            except ValueError:
                continue
        width = int(m.get("width", 0))
        bp = best_pos.get(mid)
        if bp is None:
            color, label = default_color, default_label
        else:
            center = bp[0] + width // 2
            color, label = domain_at(center)
        domain_colors[num] = color
        domain_labels.setdefault(color, label)

    return domain_colors, domain_labels, ref_name


def _luminance(hex_color: str) -> float:
    """Return the relative luminance (0-1) of a hex colour string."""
    h = hex_color.lstrip("#")
    if len(h) != 6:
        return 1.0
    r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
    def _lin(c):
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * _lin(r) + 0.7152 * _lin(g) + 0.0722 * _lin(b)


def _fmt_pvalue(pv) -> str:
    """Format p-value in MEME style: 1.23e-189, 0.00e+0 for None/0."""
    if pv is None or pv == 0:
        return "0.00e+0"
    try:
        return f"{pv:.2e}"
    except (TypeError, ValueError):
        return "n/a"


def build_figure(seq_map, motif_map, sites_by_seq, seq_pvalues,
                 label="", palette=None,
                 domain_colors=None, domain_labels=None,
                 bar_color="", bar_edge_color="",
                 bg_color="white", font_scale=1.4):
    """Build the MEME web-default motif location figure.

    Layout (top to bottom):
        1. Column header row:  Name  |  p-value  |  Motif Locations
        2. One row per sequence with backbone line + coloured motif blocks
        3. Optional structural-domain legend (when domain_colors is set)
        4. Motif Consensus table:  Motif # | Symbol | Motif Consensus

    Parameters
    ----------
    seq_pvalues : dict
        Maps sequence_id -> combined p-value (float) from MEME's
        <scanned_sites pvalue=...> attribute.
    font_scale : float
        Multiplier applied to every text size.  Default 1.4 (web-default
        sizes are sized for screen, not print).
    """
    dark = _luminance(bg_color) < 0.4 if bg_color and bg_color != "white" else False

    if dark:
        _c = {
            "seq_label":    "#EEEEEE",
            "pvalue_text":  "#DDDDDD",
            "header_text":  "#FFFFFF",
            "header_rule":  "#666666",
            "ruler_text":   "#DDDDDD",
            "ruler_tick":   "#AAAAAA",
            "ruler_axis":   "#AAAAAA",
            "pos_label":    "#DDDDDD",
            "alt_row":      "#1C1C1C",
            "title":        "#FFFFFF",
            "legend_edge":  "#555555",
            "legend_text":  "#EEEEEE",
            "legend_patch_edge": "#CCCCCC",
            "table_grid":   "#555555",
            "table_text":   "#EEEEEE",
            "consensus_text": "#F5F5F5",
        }
        bar_color = bar_color or "#888888"
    else:
        _c = {
            "seq_label":    "#111111",
            "pvalue_text":  "#222222",
            "header_text":  "#000000",
            "header_rule":  "#888888",
            "ruler_text":   "#222222",
            "ruler_tick":   "#333333",
            "ruler_axis":   "#555555",
            "pos_label":    "#222222",
            "alt_row":      "#F8F8F8",
            "title":        "black",
            "legend_edge":  "#888888",
            "legend_text":  "black",
            "legend_patch_edge": "#555555",
            "table_grid":   "#666666",
            "table_text":   "#111111",
            "consensus_text": "#111111",
        }
        bar_color = bar_color or "#666666"

    if palette is None:
        palette = DEFAULT_PALETTE
    motif_ids_sorted = sorted(motif_map.keys(),
                              key=lambda x: int(x.split("_")[-1]))

    # Per-motif colour (domain-override or palette by index)
    if domain_colors:
        color_of = {}
        for mid in motif_ids_sorted:
            num = motif_map[mid]["num"]
            color_of[mid] = domain_colors.get(
                num, domain_colors.get(str(num), "#999999"))
    else:
        color_of = {mid: palette[i % len(palette)]
                    for i, mid in enumerate(motif_ids_sorted)}

    # Scaled font sizes  (web-default sizes multiplied by font_scale) --------
    FS_HEADER    = 13.5 * font_scale
    FS_SEQNAME   = 11.0 * font_scale
    FS_PVALUE    = 10.5 * font_scale
    FS_RULER     =  9.5 * font_scale
    FS_POSLABEL  = 12.0 * font_scale
    FS_TITLE     = 15.0 * font_scale
    FS_LEGEND    = 11.5 * font_scale
    FS_LEGTITLE  = 12.5 * font_scale
    FS_TBL_HDR   = 12.0 * font_scale
    FS_TBL_NUM   = 11.0 * font_scale
    FS_TBL_SEQ   =  9.5 * font_scale   # monospace consensus

    # Layout constants -------------------------------------------------------
    N_SEQ   = len(seq_map)
    N_MOT   = len(motif_ids_sorted)
    MAX_LEN = max(v["length"] for v in seq_map.values())

    # Column widths in normalized figure coords (sum to 1.0)
    NAME_W    = 0.22   # gene name column
    PVAL_W    = 0.11   # p-value column
    GAP_W     = 0.015  # gap between p-value and locations panel
    LOC_LEFT  = NAME_W + PVAL_W + GAP_W
    LOC_W     = 1.0 - LOC_LEFT - 0.02

    BAR_THICK = 0.06   # backbone line thickness (data y-units)
    MOTIF_H   = 0.55   # motif block height (data y-units)
    ROW_H     = 1.00   # one data unit per row

    # Figure dimensions scale with #seq, font, and consensus-table rows
    FIG_W = max(20.0, 18.0 + font_scale * 2.5)
    table_rows_h = N_MOT * 0.32 * font_scale
    legend_extra = 0.6 * font_scale if (domain_colors and domain_labels) else 0.0
    FIG_H = max(11.0, N_SEQ * 0.50 * font_scale + 3.5 + table_rows_h + legend_extra)

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor(bg_color)

    # GridSpec: main plot on top, legend(s) + consensus table at the bottom
    if domain_colors and domain_labels:
        # 3 stacked panels: main / domain legend / consensus table
        legend_h_frac = 0.06 * font_scale / FIG_H
        table_h_frac  = max(0.18, (table_rows_h + 1.2) / FIG_H)
        main_h_frac   = 1.0 - legend_h_frac - table_h_frac - 0.04
        gs = fig.add_gridspec(
            3, 1, height_ratios=[main_h_frac, legend_h_frac, table_h_frac],
            hspace=0.04,
            left=0.015, right=0.985, top=0.965, bottom=0.025)
        ax_main   = fig.add_subplot(gs[0])
        ax_legend = fig.add_subplot(gs[1])
        ax_table  = fig.add_subplot(gs[2])
    else:
        table_h_frac = max(0.20, (table_rows_h + 1.2) / FIG_H)
        main_h_frac  = 1.0 - table_h_frac - 0.04
        gs = fig.add_gridspec(
            2, 1, height_ratios=[main_h_frac, table_h_frac],
            hspace=0.05,
            left=0.015, right=0.985, top=0.965, bottom=0.025)
        ax_main   = fig.add_subplot(gs[0])
        ax_legend = None
        ax_table  = fig.add_subplot(gs[1])

    # ---- MAIN PANEL -------------------------------------------------------
    ax = ax_main
    ax.set_facecolor(bg_color)
    ax.set_xlim(0, 1)
    # +2 rows: 1 for header row, 1 for ruler/position label
    ax.set_ylim(-1.2, N_SEQ + 1.0)
    ax.axis("off")

    def sx(pos):
        return LOC_LEFT + (pos / MAX_LEN) * LOC_W

    def mw(width):
        return (width / MAX_LEN) * LOC_W

    # Header row
    header_y = N_SEQ + 0.45
    ax.text(NAME_W - 0.01, header_y, "Name",
            ha="right", va="center", fontsize=FS_HEADER,
            fontweight="bold", color=_c["header_text"])
    ax.text(NAME_W + PVAL_W - 0.005, header_y, "p-value",
            ha="right", va="center", fontsize=FS_HEADER,
            fontweight="bold", color=_c["header_text"])
    ax.text(LOC_LEFT + LOC_W / 2, header_y, "Motif Locations",
            ha="center", va="center", fontsize=FS_HEADER,
            fontweight="bold", color=_c["header_text"])
    # Underline below header
    rule_y = header_y - 0.32
    ax.plot([0.015, 0.985], [rule_y, rule_y],
            color=_c["header_rule"], lw=1.0, transform=ax.transAxes,
            clip_on=False)
    # Note: transform=ax.transAxes uses axes coords, but our data coords
    # are already 0..1 in x; redo the rule in data coords for x.
    ax.plot([0.005, 0.995], [rule_y, rule_y],
            color=_c["header_rule"], lw=1.0, zorder=2)

    # Position ruler (above first sequence, below header)
    ruler_y = N_SEQ - 0.15
    tick_interval = 100 if MAX_LEN < 1000 else 200
    for tv in range(0, MAX_LEN + 1, tick_interval):
        rx = sx(tv)
        ax.plot([rx, rx], [ruler_y - 0.14, ruler_y],
                color=_c["ruler_tick"], lw=0.8)
        ax.text(rx, ruler_y + 0.10, str(tv), ha="center", va="bottom",
                fontsize=FS_RULER, color=_c["ruler_text"])
    ax.plot([sx(0), sx(MAX_LEN)], [ruler_y - 0.07, ruler_y - 0.07],
            color=_c["ruler_axis"], lw=1.0)

    # Sequence rows
    for row_i, sid in enumerate(seq_map.keys()):
        y   = N_SEQ - 1 - row_i - 0.4
        seq = seq_map[sid]

        # Alternating row stripe across the WHOLE row width for readability
        if row_i % 2 == 0:
            ax.add_patch(mpatches.Rectangle(
                (0.005, y - 0.45), 0.99, 0.90,
                facecolor=_c["alt_row"], edgecolor="none", zorder=-1))

        # Name column (right-aligned)
        ax.text(NAME_W - 0.01, y, seq["name"],
                ha="right", va="center", fontsize=FS_SEQNAME,
                fontfamily="monospace",
                fontweight="bold" if dark else "normal",
                color=_c["seq_label"])

        # p-value column (right-aligned)
        pv_str = _fmt_pvalue(seq_pvalues.get(sid))
        ax.text(NAME_W + PVAL_W - 0.005, y, pv_str,
                ha="right", va="center", fontsize=FS_PVALUE,
                fontfamily="monospace", color=_c["pvalue_text"])

        # Backbone line for the sequence (thin, not a fat bar)
        bx0 = sx(0)
        bx1 = sx(seq["length"])
        ax.add_patch(mpatches.Rectangle(
            (bx0, y - BAR_THICK / 2), bx1 - bx0, BAR_THICK,
            facecolor=bar_color, edgecolor="none", zorder=2))

        # Motif blocks (rectangular, MEME-style)
        for site in sites_by_seq.get(sid, []):
            mid = site["motif_id"]
            bx  = sx(site["position"])
            bw_ = mw(motif_map[mid]["width"])
            col = color_of[mid]
            ax.add_patch(mpatches.Rectangle(
                (bx, y - MOTIF_H / 2), bw_, MOTIF_H,
                facecolor=col, edgecolor="black" if not dark else "white",
                linewidth=0.7, alpha=0.95, zorder=3))

    # Position label below the rows (centered under locations panel)
    ax.text(LOC_LEFT + LOC_W / 2, -0.75,
            "Position (aa)", ha="center", va="center",
            fontsize=FS_POSLABEL, color=_c["pos_label"], style="italic")

    # ---- STRUCTURAL DOMAIN LEGEND (optional) ------------------------------
    if ax_legend is not None:
        ax_legend.set_facecolor(bg_color)
        ax_legend.axis("off")
        seen = {}
        for mid in motif_ids_sorted:
            col = color_of[mid]
            num = motif_map[mid]["num"]
            seen.setdefault(col, []).append(str(num))
        handles = []
        for col, label_text in domain_labels.items():
            if col in seen:
                motif_nums = ", ".join(f"M{n}" for n in sorted(seen[col], key=int))
                handles.append(mpatches.Patch(
                    facecolor=col, edgecolor=_c["legend_patch_edge"],
                    linewidth=0.8,
                    label=f"{label_text}  [{motif_nums}]"))
        ncols = min(len(handles), 4)
        leg = ax_legend.legend(
            handles=handles, title="Structural Domain",
            loc="center", ncol=ncols,
            fontsize=FS_LEGEND, title_fontsize=FS_LEGTITLE,
            frameon=True, framealpha=0.92, edgecolor=_c["legend_edge"],
            handlelength=2.2, handleheight=1.6)
        leg.get_frame().set_facecolor(bg_color)
        for text in leg.get_texts():
            text.set_color(_c["legend_text"])
        if leg.get_title():
            leg.get_title().set_color(_c["legend_text"])
            leg.get_title().set_fontweight("bold")

    # ---- MOTIF CONSENSUS TABLE -------------------------------------------
    ax_t = ax_table
    ax_t.set_facecolor(bg_color)
    ax_t.set_xlim(0, 1)
    ax_t.set_ylim(-0.5, N_MOT + 1.5)
    ax_t.axis("off")

    # Column boundaries within ax_t (normalized 0..1 in x)
    COL_NUM_X  = 0.06   # right edge of "Motif" column (number is centered at COL_NUM_X/2)
    COL_SYM_X  = 0.16   # right edge of "Symbol" column
    CONS_LEFT  = 0.18

    # Background panel
    panel_y0 = -0.45
    panel_h  = N_MOT + 1.30
    ax_t.add_patch(mpatches.Rectangle(
        (0.005, panel_y0), 0.99, panel_h,
        facecolor=bg_color, edgecolor=_c["table_grid"],
        linewidth=1.0, zorder=0))

    # Header
    hdr_y = N_MOT + 0.55
    ax_t.text(COL_NUM_X / 2, hdr_y, "Motif",
              ha="center", va="center", fontsize=FS_TBL_HDR,
              fontweight="bold", color=_c["table_text"])
    ax_t.text((COL_NUM_X + COL_SYM_X) / 2, hdr_y, "Symbol",
              ha="center", va="center", fontsize=FS_TBL_HDR,
              fontweight="bold", color=_c["table_text"])
    ax_t.text(CONS_LEFT, hdr_y, "Motif Consensus",
              ha="left", va="center", fontsize=FS_TBL_HDR,
              fontweight="bold", color=_c["table_text"])

    # Header underline
    ax_t.plot([0.005, 0.995], [N_MOT + 0.10, N_MOT + 0.10],
              color=_c["table_grid"], lw=0.8, zorder=1)
    # Column separators
    for x_sep in (COL_NUM_X, COL_SYM_X):
        ax_t.plot([x_sep, x_sep], [panel_y0, N_MOT + 1.10],
                  color=_c["table_grid"], lw=0.6, zorder=1)

    # Symbol swatch geometry (centered between COL_NUM_X and COL_SYM_X)
    sw_cx = (COL_NUM_X + COL_SYM_X) / 2
    sw_w  = (COL_SYM_X - COL_NUM_X) * 0.55
    sw_h  = 0.50

    for i, mid in enumerate(motif_ids_sorted):
        m = motif_map[mid]
        y = N_MOT - 1 - i

        # Row stripe
        if i % 2 == 0:
            ax_t.add_patch(mpatches.Rectangle(
                (0.005, y - 0.5), 0.99, 1.0,
                facecolor=_c["alt_row"], edgecolor="none", zorder=-1))

        # Motif number
        ax_t.text(COL_NUM_X / 2, y, f"{m['num']}.",
                  ha="center", va="center", fontsize=FS_TBL_NUM,
                  color=_c["table_text"], fontweight="bold")

        # Symbol swatch
        ax_t.add_patch(mpatches.Rectangle(
            (sw_cx - sw_w / 2, y - sw_h / 2), sw_w, sw_h,
            facecolor=color_of[mid], edgecolor=_c["legend_patch_edge"],
            linewidth=0.7, zorder=2))

        # Consensus sequence (monospace)
        ax_t.text(CONS_LEFT, y, m["consensus"],
                  ha="left", va="center", fontsize=FS_TBL_SEQ,
                  fontfamily="monospace", color=_c["consensus_text"])

    # ---- TITLE -----------------------------------------------------------
    title = (f"MEME Motif Locations - {label}  ({N_SEQ} sequences, {N_MOT} motifs)"
             if label else
             f"MEME Motif Locations  ({N_SEQ} sequences, {N_MOT} motifs)")
    fig.suptitle(title, y=0.992, fontsize=FS_TITLE, fontweight="bold",
                 ha="center", color=_c["title"])
    return fig


def main():
    args = parse_args()

    xml_path = Path(args.meme_xml)
    if not xml_path.is_file():
        print(f"ERROR: meme.xml not found: {xml_path}", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    seq_map, motif_map, sites_by_seq, seq_pvalues = load_meme_xml(xml_path)

    if not seq_map:
        print("ERROR: no sequences found in meme.xml", file=sys.stderr)
        sys.exit(1)

    if args.phylo_order:
        phylo_names = parse_phylo_order(args.phylo_order)
        if phylo_names:
            seq_map = apply_phylo_order(seq_map, phylo_names)
            print(f"Applied phylo order from {args.phylo_order} "
                  f"({len(phylo_names)} names parsed)")
        else:
            print(f"WARNING: no names parsed from {args.phylo_order}; "
                  "using original MEME order", file=sys.stderr)

    palette = resolve_palette(args.palette)

    domain_colors = None
    domain_labels = None
    if args.domain_colors:
        try:
            raw = json.loads(args.domain_colors)
            domain_colors = {int(k): v for k, v in raw.items()}
        except (json.JSONDecodeError, ValueError) as exc:
            print(f"ERROR: --domain-colors is not valid JSON: {exc}", file=sys.stderr)
            sys.exit(1)
    if args.domain_labels:
        try:
            domain_labels = json.loads(args.domain_labels)
        except (json.JSONDecodeError, ValueError) as exc:
            print(f"ERROR: --domain-labels is not valid JSON: {exc}", file=sys.stderr)
            sys.exit(1)

    if args.ref_seq_ids and args.domain_ranges_json:
        try:
            ranges_spec = json.loads(args.domain_ranges_json)
        except (json.JSONDecodeError, ValueError) as exc:
            print(f"ERROR: --domain-ranges-json is not valid JSON: {exc}",
                  file=sys.stderr)
            sys.exit(1)
        ref_ids = [s.strip() for s in args.ref_seq_ids.split(",") if s.strip()]
        auto_dc, auto_dl, ref_name = derive_domains_from_reference(
            str(xml_path), ref_ids, ranges_spec)
        if auto_dc:
            domain_colors = auto_dc
            domain_labels = auto_dl
            print(f"Position-based domain colors: mapped {len(auto_dc)} motifs "
                  f"via reference {ref_name!r}")
        else:
            print(f"WARNING: no reference sequence in {ref_ids} found in "
                  f"{xml_path}; falling back to static domain mapping",
                  file=sys.stderr)

    bar_color      = args.bar_color      if args.bar_color      else ""
    bar_edge_color = args.bar_edge_color if args.bar_edge_color else ""
    bg_color       = args.bg_color       if args.bg_color       else "white"

    fig = build_figure(seq_map, motif_map, sites_by_seq, seq_pvalues,
                       label=args.label, palette=palette,
                       domain_colors=domain_colors, domain_labels=domain_labels,
                       bar_color=bar_color, bar_edge_color=bar_edge_color,
                       bg_color=bg_color, font_scale=args.font_scale)

    stem = f"{args.label}_motif_locations_full" if args.label \
           else "motif_locations_full"

    png_out = out_dir / f"{stem}.png"
    fig.savefig(png_out, dpi=args.dpi, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print(f"Saved: {png_out}")

    if args.svg:
        svg_out = out_dir / f"{stem}.svg"
        fig.savefig(svg_out, format="svg", bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"Saved: {svg_out}")

    plt.close(fig)


if __name__ == "__main__":
    main()
