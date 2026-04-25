#!/usr/bin/env python3
"""
Full-view motif location diagram — all sequences, all motifs.

Reads a MEME XML file and draws a publication-ready PNG (and optionally SVG)
where every input sequence is shown as a horizontal bar with coloured blocks
for each motif hit.  Unlike the default meme-motif-locations.svg, this script
never truncates the sequence list regardless of dataset size.

Usage:
    python3 plot_motif_locations.py \
        --meme-xml  /path/to/meme.xml \
        --outdir    /path/to/output_dir \
        --label     my_dataset \
        [--dpi      200] \
        [--svg] \
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


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--meme-xml",  required=True,  help="Path to meme.xml")
    p.add_argument("--outdir",    required=True,  help="Output directory")
    p.add_argument("--label",     default="",     help="Dataset label (used in title)")
    p.add_argument("--dpi",       type=int, default=200, help="Output DPI (default 200)")
    p.add_argument("--svg",       action="store_true",   help="Also write SVG")
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
                   help="Hex color for the sequence backbone bar (default #DDEEFF).")
    p.add_argument("--bar-edge-color", default="",
                   help="Hex color for the sequence backbone bar edge (default #99AACC).")
    p.add_argument("--bg-color", default="",
                   help="Figure background hex color, e.g. '#111111'. "
                        "Dark backgrounds (luminance < 0.4) automatically switch "
                        "all text and decorations to light variants.")
    p.add_argument("--phylo-order", default="",
                   help="Path to a CLUSTAL alignment file whose sequence order "
                        "defines the top-to-bottom order of rows in the figure. "
                        "Sequences not present in the file are appended at the "
                        "bottom in their original MEME order.")
    return p.parse_args()


# ── parse XML ─────────────────────────────────────────────────────────────────
def load_meme_xml(xml_path):
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
        motif_map[mid] = {
            "num":    num,
            "width":  int(m.get("width")),
            "evalue": m.get("e_value", "?"),
        }

    sites_by_seq = {sid: [] for sid in seq_map}
    for ss in root.findall(".//scanned_sites"):
        sid = ss.get("sequence_id")
        for site in ss.findall("scanned_site"):
            sites_by_seq[sid].append({
                "motif_id": site.get("motif_id"),
                "position": int(site.get("position")),
                "pvalue":   float(site.get("pvalue")),
            })

    return seq_map, motif_map, sites_by_seq


# ── phylo ordering ────────────────────────────────────────────────────────────
def parse_phylo_order(aln_path: str) -> list:
    """Return sequence names in the order they first appear in a CLUSTAL file.

    Handles CLUSTAL W / CLUSTAL OMEGA format.  Only non-empty lines that do
    not start with 'CLUSTAL' or '*:.' (conservation lines) and that contain
    a whitespace-separated name are treated as sequence rows.
    """
    seen  = set()
    order = []
    try:
        with open(aln_path, encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.rstrip()
                if not line:
                    continue
                # Skip CLUSTAL header and conservation rows
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
    """Reorder seq_map (OrderedDict semantics) by phylo_names.

    Sequences whose names match (exact or substring) are placed first in
    phylo order; unmatched sequences are appended at the end.
    """
    # Build name → sid lookup
    name_to_sid = {v["name"]: k for k, v in seq_map.items()}

    ordered_sids = []
    used = set()
    for pname in phylo_names:
        # Exact match
        if pname in name_to_sid:
            sid = name_to_sid[pname]
            if sid not in used:
                ordered_sids.append(sid)
                used.add(sid)
            continue
        # Substring match (handle truncation in some CLUSTAL files)
        for seq_name, sid in name_to_sid.items():
            if pname in seq_name or seq_name in pname:
                if sid not in used:
                    ordered_sids.append(sid)
                    used.add(sid)
                break

    # Append any sequences not matched by the phylo order
    for sid in seq_map:
        if sid not in used:
            ordered_sids.append(sid)

    return {sid: seq_map[sid] for sid in ordered_sids}


# ── colour palette (colourblind-safe, up to 10 motifs) ───────────────────────
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
    """Return a list of hex color strings from --palette argument.

    Accepts:
      - Empty string           → DEFAULT_PALETTE
      - Comma-separated hex    → split on commas
      - JSON array string      → json.loads
    """
    if not palette_arg:
        return DEFAULT_PALETTE

    palette_arg = palette_arg.strip()

    # JSON array
    if palette_arg.startswith("["):
        try:
            colors = json.loads(palette_arg)
            if isinstance(colors, list) and all(isinstance(c, str) for c in colors):
                return colors
        except (json.JSONDecodeError, TypeError):
            pass

    # Comma-separated hex codes
    if "," in palette_arg:
        return [c.strip() for c in palette_arg.split(",") if c.strip()]

    # Newline-separated hex codes (parse_toml.py list output piped via shell)
    if "\n" in palette_arg:
        return [c.strip() for c in palette_arg.split("\n") if c.strip()]

    # Single color (unlikely but safe)
    return [palette_arg]


def _luminance(hex_color: str) -> float:
    """Return the relative luminance (0–1) of a hex colour string."""
    h = hex_color.lstrip("#")
    if len(h) != 6:
        return 1.0  # assume light on bad input
    r, g, b = int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255
    def _lin(c):
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * _lin(r) + 0.7152 * _lin(g) + 0.0722 * _lin(b)


def build_figure(seq_map, motif_map, sites_by_seq, label="", palette=None,
                 domain_colors=None, domain_labels=None,
                 bar_color="", bar_edge_color="",
                 bg_color="white"):
    """Build the full-view motif location figure.

    Parameters
    ----------
    domain_colors : dict or None
        Maps motif number (int) → hex color string.  When provided, motif
        blocks are coloured by structural domain rather than palette index.
    domain_labels : dict or None
        Maps hex color string → human-readable domain label.  Used to build
        a grouped domain legend when domain_colors is set.
    bar_color : str
        Fill colour for the sequence backbone bar.  Auto-selected by bg_color
        when empty.
    bar_edge_color : str
        Edge colour for the sequence backbone bar.  Auto-selected by bg_color
        when empty.
    bg_color : str
        Figure background colour.  Dark values (luminance < 0.4) trigger
        automatic light-text mode for all decorations.
    """
    dark = _luminance(bg_color) < 0.4 if bg_color and bg_color != "white" else False

    # Colour token dict — keys used throughout the function
    if dark:
        _c = {
            "seq_label":    "#EEEEEE",
            "len_label":    "#CCCCCC",
            "ruler_text":   "#DDDDDD",
            "ruler_tick":   "#AAAAAA",
            "ruler_axis":   "#AAAAAA",
            "pos_label":    "#DDDDDD",
            "alt_row":      "#1C1C1C",
            "title":        "#FFFFFF",
            "legend_edge":  "#555555",
            "legend_text":  "#EEEEEE",
            "legend_patch_edge": "#CCCCCC",
        }
        bar_color      = bar_color      or "#1A2A3A"
        bar_edge_color = bar_edge_color or "#2C4A6E"
    else:
        _c = {
            "seq_label":    "#111111",
            "len_label":    "#555555",
            "ruler_text":   "#222222",
            "ruler_tick":   "#333333",
            "ruler_axis":   "#555555",
            "pos_label":    "#222222",
            "alt_row":      "#F8F8F8",
            "title":        "black",
            "legend_edge":  "#CCCCCC",
            "legend_text":  "black",
            "legend_patch_edge": "#555555",
        }
        bar_color      = bar_color      or "#DDEEFF"
        bar_edge_color = bar_edge_color or "#99AACC"
    if palette is None:
        palette = DEFAULT_PALETTE
    motif_ids_sorted = sorted(motif_map.keys(), key=lambda x: int(x.split("_")[-1]))

    # Resolve per-motif color ─────────────────────────────────────────────────
    if domain_colors:
        color_of = {}
        for mid in motif_ids_sorted:
            num = motif_map[mid]["num"]
            color_of[mid] = domain_colors.get(num, domain_colors.get(str(num), "#999999"))
    else:
        color_of = {mid: palette[i % len(palette)] for i, mid in enumerate(motif_ids_sorted)}

    N_SEQ   = len(seq_map)
    BAR_H   = 0.30
    MOTIF_H = 0.40
    MAX_LEN = max(v["length"] for v in seq_map.values())
    LABEL_W = 0.32          # wider label area for long sequence names
    SEQ_W   = 1 - LABEL_W - 0.04

    ROW_H   = 0.65          # vertical spacing per row (in data units)
    FIG_W   = 22            # wider figure for readability
    FIG_H   = max(14, N_SEQ * ROW_H + 4.5)

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor(bg_color)
    ax.set_facecolor(bg_color)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, N_SEQ - 0.5)
    ax.axis("off")

    def sx(pos):
        return LABEL_W + (pos / MAX_LEN) * SEQ_W

    def mw(width):
        return (width / MAX_LEN) * SEQ_W

    # ruler
    ruler_y = N_SEQ - 0.15
    tick_interval = 100 if MAX_LEN < 1000 else 200
    for tv in range(0, MAX_LEN + 1, tick_interval):
        rx = sx(tv)
        ax.plot([rx, rx], [ruler_y - 0.14, ruler_y], color=_c["ruler_tick"], lw=0.7)
        ax.text(rx, ruler_y + 0.06, str(tv), ha="center", va="bottom",
                fontsize=8, color=_c["ruler_text"])
    ax.plot([sx(0), sx(MAX_LEN)], [ruler_y - 0.07, ruler_y - 0.07],
            color=_c["ruler_axis"], lw=1.0)
    ax.text(LABEL_W + SEQ_W / 2, ruler_y + 0.35,
            "Position (aa)", ha="center", va="bottom", fontsize=10,
            color=_c["pos_label"])

    # sequence rows
    for row_i, sid in enumerate(seq_map.keys()):
        y   = N_SEQ - 1 - row_i
        seq = seq_map[sid]

        # Sequence name — readable, monospace, slightly bold on dark bg
        ax.text(LABEL_W - 0.007, y, seq["name"],
                ha="right", va="center", fontsize=9,
                fontfamily="monospace",
                fontweight="bold" if dark else "normal",
                color=_c["seq_label"])

        bx0 = sx(0)
        bx1 = sx(seq["length"])
        ax.add_patch(mpatches.FancyBboxPatch(
            (bx0, y - BAR_H / 2), bx1 - bx0, BAR_H,
            boxstyle="round,pad=0.002",
            facecolor=bar_color, edgecolor=bar_edge_color, linewidth=0.6))

        ax.text(bx1 + 0.006, y, f"{seq['length']} aa",
                ha="left", va="center", fontsize=7.5, color=_c["len_label"])

        for site in sites_by_seq.get(sid, []):
            mid = site["motif_id"]
            bx  = sx(site["position"])
            bw_ = mw(motif_map[mid]["width"])
            col = color_of[mid]
            ax.add_patch(mpatches.FancyBboxPatch(
                (bx, y - MOTIF_H / 2), bw_, MOTIF_H,
                boxstyle="round,pad=0.001",
                facecolor=col, edgecolor="white" if dark else "#333333",
                linewidth=0.5, alpha=0.93))
            if bw_ > 0.010:
                lum = _luminance(col)
                txt_col = "#FFFFFF" if lum < 0.45 else "#111111"
                ax.text(bx + bw_ / 2, y, str(motif_map[mid]["num"]),
                        ha="center", va="center",
                        fontsize=6.5, color=txt_col, fontweight="bold")

        if row_i % 2 == 0:
            ax.add_patch(mpatches.FancyBboxPatch(
                (0, y - 0.5), 1, 1,
                boxstyle="square,pad=0",
                facecolor=_c["alt_row"], edgecolor="none", zorder=-1))

    # legend ──────────────────────────────────────────────────────────────────
    if domain_colors and domain_labels:
        seen = {}
        for mid in motif_ids_sorted:
            col = color_of[mid]
            num = motif_map[mid]["num"]
            if col not in seen:
                seen[col] = []
            seen[col].append(str(num))
        handles = []
        for col, label_text in domain_labels.items():
            if col in seen:
                motif_nums = ", ".join(f"M{n}" for n in sorted(seen[col], key=int))
                handles.append(mpatches.Patch(
                    facecolor=col, edgecolor=_c["legend_patch_edge"], linewidth=0.7,
                    label=f"{label_text}  [{motif_nums}]"))
        legend_title = "Structural Domain"
        ncols = min(len(handles), 4)
    else:
        handles = [
            mpatches.Patch(
                facecolor=color_of[mid], edgecolor=_c["legend_patch_edge"], linewidth=0.6,
                label=f"Motif {motif_map[mid]['num']}  "
                      f"(w={motif_map[mid]['width']}, E={motif_map[mid]['evalue']})")
            for mid in motif_ids_sorted
        ]
        legend_title = "Motifs"
        ncols = 5
    leg = ax.legend(handles=handles, title=legend_title,
                    loc="lower center", bbox_to_anchor=(0.5, -0.06),
                    ncol=ncols, fontsize=9.5, title_fontsize=10.5,
                    frameon=True, framealpha=0.92, edgecolor=_c["legend_edge"],
                    handlelength=2.0, handleheight=1.5)
    leg.get_frame().set_facecolor(bg_color)
    for text in leg.get_texts():
        text.set_color(_c["legend_text"])
    if leg.get_title():
        leg.get_title().set_color(_c["legend_text"])
        leg.get_title().set_fontweight("bold")

    title = (f"MEME Motif Locations — {label}  ({N_SEQ} sequences)" if label
             else f"MEME Motif Locations  ({N_SEQ} sequences)")
    fig.suptitle(title, y=0.997, fontsize=13, fontweight="bold", ha="center",
                 color=_c["title"])
    plt.tight_layout(rect=[0, 0.05, 1, 0.995])
    return fig


def main():
    args = parse_args()

    xml_path = Path(args.meme_xml)
    if not xml_path.is_file():
        print(f"ERROR: meme.xml not found: {xml_path}", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    seq_map, motif_map, sites_by_seq = load_meme_xml(xml_path)

    if not seq_map:
        print("ERROR: no sequences found in meme.xml", file=sys.stderr)
        sys.exit(1)

    # Apply phylo ordering if provided
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

    bar_color      = args.bar_color      if args.bar_color      else ""
    bar_edge_color = args.bar_edge_color if args.bar_edge_color else ""
    bg_color       = args.bg_color       if args.bg_color       else "white"

    fig = build_figure(seq_map, motif_map, sites_by_seq, label=args.label,
                       palette=palette,
                       domain_colors=domain_colors, domain_labels=domain_labels,
                       bar_color=bar_color, bar_edge_color=bar_edge_color,
                       bg_color=bg_color)

    stem = f"{args.label}_motif_locations_full" if args.label \
           else "motif_locations_full"

    png_out = out_dir / f"{stem}.png"
    fig.savefig(png_out, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved: {png_out}")

    if args.svg:
        svg_out = out_dir / f"{stem}.svg"
        fig.savefig(svg_out, format="svg", bbox_inches="tight")
        print(f"Saved: {svg_out}")

    plt.close(fig)


if __name__ == "__main__":
    main()
