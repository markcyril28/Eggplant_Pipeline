#!/usr/bin/env python3
"""
visualize_blast_results.py
==========================
Generates two publication-quality figures from BLASTn curated results:

  Figure 1  - Cross-species % Identity Heatmap
  Figure 1b - Cross-species % Identity + E-value Heatmap
  Figure 3  - Per-gene Top-Hits Lollipop (ranked by Bit Score)

Usage (auto-discover CSVs from curated_results directory):
    python3 visualize_blast_results.py

Usage (explicit CSV path):
    python3 visualize_blast_results.py \\
        --plant-csv  path/to/*_plant_only.csv \\
        --out-dir    path/to/figures/

Options:
    --results-dir   Directory containing *_plant_only.csv (default: auto-detected)
    --plant-csv     Explicit path to *_plant_only.csv (overrides --results-dir)
    --out-dir       Output directory for figures (default: <results-dir>/figures/)
    --top-n         Top-N hits to show in the Fig 3 lollipop (default: 10)
"""

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class VizConfig:
    """Visualization settings read from [blast_visualize] TOML section."""
    colormap: str = "RdYlGn"
    figure_dpi: int = 150
    save_dpi: int = 300
    heatmap_vmin: float = 65.0
    heatmap_vmax: float = 100.0
    heatmap_w_scale: float = 1.20
    heatmap_h_scale: float = 0.92
    lollipop_ncols: int = 2
    lollipop_x_pad: float = 1.60
    lollipop_dot_size: int = 100
    lollipop_dot_size_hi: int = 150
    hi_stem_color: str = "#c7920a"
    stem_color: str = "#d1d5db"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpe
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import FancyBboxPatch

# ── Global style ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

# ─────────────────────────────────────────────────────────────────────────────
# Gene metadata: update GENE_LABELS if the thesis naming convention differs
# ─────────────────────────────────────────────────────────────────────────────

# Display order: DMP clade ascending (DMP2 → DMP3 → DMP4 → DMP6 → DMP7a → DMP7b → DMP8/9)
GENE_ORDER = [
    "SMEL5_01g026030.1",
    "SMEL5_01g008730.1",
    "SMEL5_04g005390.1",
    "SMEL5_02g013320.1",
    "SMEL5_10g003660.1",
    "SMEL5_12g005350.1",
    "SMEL5_10g017610.1",
]

# Two-line label used in panel titles (update to match thesis convention if needed)
GENE_LABELS = {
    "SMEL5_01g026030.1": "SmelDMP01.030\n(DMP2-like)",
    "SMEL5_01g008730.1": "SmelDMP01.730\n(DMP3-like)",
    "SMEL5_04g005390.1": "SmelDMP04.390\n(DMP4-like)",
    "SMEL5_02g013320.1": "SmelDMP02.320\n(DMP6-like)",
    "SMEL5_10g003660.1": "SmelDMP10.660\n(DMP7-like)",
    "SMEL5_12g005350.1": "SmelDMP12.350\n(DMP7b-like)",
    "SMEL5_10g017610.1": "SmelDMP10.610\n(DMP8/9-like)",
}

# Single-line short name for axis tick labels
GENE_SHORT = {
    "SMEL5_01g026030.1": "SmelDMP01.030",
    "SMEL5_01g008730.1": "SmelDMP01.730",
    "SMEL5_04g005390.1": "SmelDMP04.390",
    "SMEL5_02g013320.1": "SmelDMP02.320",
    "SMEL5_10g003660.1": "SmelDMP10.660",
    "SMEL5_12g005350.1": "SmelDMP12.350",
    "SMEL5_10g017610.1": "SmelDMP10.610",
}

GENE_CLADE = {
    "SMEL5_01g026030.1": "DMP2",
    "SMEL5_01g008730.1": "DMP3",
    "SMEL5_04g005390.1": "DMP4",
    "SMEL5_02g013320.1": "DMP6",
    "SMEL5_10g003660.1": "DMP7a",
    "SMEL5_12g005350.1": "DMP7b",
    "SMEL5_10g017610.1": "DMP8/9",
}

# ─────────────────────────────────────────────────────────────────────────────
# Species metadata  (taxonomic order: Solanaceae → Asterids → Rosids → other)
# ─────────────────────────────────────────────────────────────────────────────

SPECIES_ORDER = [
    "Sl", "Ca", "Nt", "Ib", "Cs", "Pt", "Gr", "At", "Br", "Gm", "Mt",
    "Ma", "Si", "Bd", "Pp",  # added with the older 2024_Paper_DMPs panel
]

SPECIES_LABELS = {
    "Sl": "S. lycopersicum",
    "Ca": "C. annuum",
    "Nt": "N. tabacum",
    "Ib": "I. batatas",
    "Cs": "C. sativus",
    "Pt": "P. trichocarpa",
    "Gr": "G. raimondii",
    "At": "A. thaliana",
    "Br": "B. rapa",
    "Gm": "G. max",
    "Mt": "M. truncatula",
    "Ma": "M. acuminata",
    "Si": "S. italica",
    "Bd": "B. distachyon",
    "Pp": "P. patens",
}

# Family grouping for the heatmap annotation bar
# Each tuple: (family_name, [species_codes_in_order])
FAMILY_GROUPS = [
    ("Solanaceae",      ["Sl", "Ca", "Nt"]),
    ("Convolvulaceae",  ["Ib"]),
    ("Cucurbitaceae",   ["Cs"]),
    ("Salicaceae",      ["Pt"]),
    ("Malvaceae",       ["Gr"]),
    ("Brassicaceae",    ["At", "Br"]),
    ("Fabaceae",        ["Gm", "Mt"]),
    ("Musaceae",        ["Ma"]),
    ("Poaceae",         ["Si", "Bd"]),
    ("Funariaceae",     ["Pp"]),
]

FAMILY_COLORS = [
    "#d62839",  # Solanaceae - warm red
    "#e0a458",  # Convolvulaceae - amber
    "#60b5d1",  # Cucurbitaceae - sky blue
    "#1a8faa",  # Salicaceae - teal
    "#0a3d62",  # Malvaceae - deep navy
    "#36b37e",  # Brassicaceae - green
    "#6baa3d",  # Fabaceae - olive-green
    "#8952d4",  # Musaceae - purple
    "#c45a8a",  # Poaceae - rose
    "#5c6b73",  # Funariaceae - slate (moss outgroup)
]

# Disambiguation notes for the haploid-inducer label set below:
#   "CsDMP9"/XP_006482605 is Citrus sinensis (NC_068561), NOT Cucumis sativus.
#     Yin et al. 2024 validated cucumber CsDMP (CsaV3_1G028660); not in query set.
#     The "Cs" prefix collision caused an earlier misclassification.
#   FASTA-labeled SlDMP3 (Solyc05g007920, chr5) IS the verified inducer per
#     Deng et al. 2025 (despite NCBI annotating XP_004239396 as "DMP8-like" -
#     same locus, so _short_label() yields "SlDMP8-like" for the protein hit).
# Citations and full label-to-locus mapping are in the shared root config
# 02_blast_alignmentCONFIG.toml under [blast_visualize].haploid_inducer_labels.
# Default fallback used when --hi-labels is not provided. The authoritative
# list lives in [blast_visualize].haploid_inducer_labels of the shared root
# config 02_blast_alignmentCONFIG.toml; the orchestrator reads it, joins with
# commas, and passes via --hi-labels. Keep this default in sync so standalone
# invocations still produce correct DMP figures.
_DEFAULT_HAPLOID_INDUCER_LABELS = frozenset({
    # Mirrors [blast_visualize].haploid_inducer_labels in the shared root
    # 02_blast_alignmentCONFIG.toml. See that file for citations.
    "AtDMP8", "AtDMP9", "AtDMP8+AtDMP9",            # A. thaliana
    "ZmDMP",                                         # Z. mays
    "SlDMP3", "SlDMP8-like",                         # S. lycopersicum (same locus)
    "StDMP",                                         # S. tuberosum
    "NtDMP", "NtDMP2-like", "NtDMP3-like",           # N. tabacum
    "CDX74441", "CDX81135", "CDY30259", "CDY56548",  # B. napus
    "BoDMP_LOC106333617", "BoDMP_LOC106333853",      # B. oleracea
    "BjuA03g54090S", "BjuA04g10430S", "BjuB08g57390S",  # B. juncea
    "ClDMP3",                                        # C. lanatus
    "CsDMP",                                         # C. sativus
    "MtDMP8", "MtDMP9",                              # M. truncatula
    "GmDMP1", "GmDMP2",                              # G. max
    "GhDMPa", "GhDMPd",                              # G. hirsutum
    "OsDMP1", "OsDMP3",                              # O. sativa
})

# Mutable module global referenced by plot_lollipop(); main() may replace this
# with the parsed --hi-labels set at runtime.
HAPLOID_INDUCER_LABELS = set(_DEFAULT_HAPLOID_INDUCER_LABELS)

# ─── Dynamic metadata (resolved at runtime for non-DMP gene groups) ──────────
_GENE_GROUP = "DMP"  # updated by _init_gene_metadata()


def _auto_short_label(gene_id):
    """Generate a readable short label from a Subject ID.

    SMEL4.1_01g000730.1.01 → SMEL4.1_01g000730
    SMEL5_10g003660.1      → SMEL5_10g003660
    """
    m = re.match(r"(SMEL\d+(?:\.\d+)?_\d+g\d+)", gene_id)
    return m.group(1) if m else gene_id.rsplit(".", 1)[0]


def _detect_gene_group(csv_path):
    """Extract gene group name from a CSV path under 3_RESULT/<GROUP>/..."""
    parts = Path(csv_path).resolve().parts
    for i, part in enumerate(parts):
        if part == "3_RESULT" and i + 1 < len(parts):
            return parts[i + 1]
    return None


def _init_gene_metadata(df, gene_group=None):
    """Configure gene metadata from CSV data.  Updates module-level globals.

    Returns the (possibly filtered) DataFrame.

    Routing rule (in priority order):
      1. If gene_group is explicitly "DMP" → use hardcoded DMP metadata.
      2. If gene_group is unknown (None) AND subject IDs overlap with the
         hardcoded DMP list → infer DMP and use DMP metadata.
      3. Any other known gene group (e.g. "HAP2") → auto-generate metadata
         from the CSV data, regardless of subject ID content.  This prevents
         cross-contamination when a HAP2 BLAST result incidentally contains
         eggplant gene IDs that also appear in the DMP gene set.
    """
    global GENE_ORDER, GENE_SHORT, GENE_CLADE, GENE_LABELS, _GENE_GROUP

    dmp_ids = set(GENE_ORDER)  # initial hardcoded DMP list
    subject_ids = set(df["Subject ID"].unique())

    is_dmp = (gene_group == "DMP") or (
        gene_group is None and bool(subject_ids & dmp_ids)
    )

    if is_dmp:
        # Keep hardcoded DMP metadata and filter to known DMP subject IDs
        _GENE_GROUP = "DMP"
        return df[df["Subject ID"].isin(dmp_ids)].copy()

    # ── Auto-generate metadata for non-DMP gene groups ────────────────────
    sorted_ids = sorted(subject_ids)
    GENE_ORDER = sorted_ids
    GENE_SHORT = {g: _auto_short_label(g) for g in sorted_ids}
    GENE_CLADE = {}
    GENE_LABELS = {g: GENE_SHORT[g] for g in sorted_ids}
    _GENE_GROUP = gene_group or "Unknown"
    return df.copy()


# ─────────────────────────────────────────────────────────────────────────────
# Parsing helpers
# ─────────────────────────────────────────────────────────────────────────────

# Hardcoded mapping for raw NCBI/genome-DB query IDs that appear in the curated
# CSV without a species-prefixed name. Resolved from the FASTA headers under
# 1_RefSeqs/d_DMP_Query_Fasta/<species>/*_merged*.fa (protein description in
# square brackets) plus the species inferred from the parent directory.
# Format: accession_stem -> friendly short label.
ACCESSION_LABELS = {
    # Nicotiana tabacum (NtDMPs_merged.fa) - RefSeq XM_016*
    "XM_016586331": "NtDMP2-like",
    "XM_016642301": "NtDMP3-like",
    "XM_016578562": "NtDMP4-like",
    "XM_016591727": "NtDMP7-like",
    "XM_016604032": "NtDMP9-like",
    # Cucumis sativus (CsDMPs_merged.fa)
    "XM_004146681": "CsDMP9",
    # Medicago truncatula (MtDMPs_merged.fa)
    "XM_003621193": "MtDMP9",   # haploid inducer (N. Wang et al., 2022)
    "XM_003614037": "MtDMP",
    # Brassica oleracea (BoDMPs_merged.fa) - uncharacterized LOC*
    "XM_013772041": "BoDMP_LOC106333617",
    "XM_013772244": "BoDMP_LOC106333853",
}


# Protein accession (XP_*) → short label, derived from header scans of the 8
# gap-filler species under II_INPUTS/DMP_query_fasta_file/. Used when BLAST
# emits the secondary `lcl|NC_*_cds_XP_*` token form as Query ID, or when a
# Query ID is just a bare XP accession.
#
# To regenerate after wiring new species, run a header scan that emits
# (XP_accession, short_label) pairs from each species' *_merged_fasta.fa.
# The short_label is the first whitespace token of the FASTA header passed
# through _short_label() (regex-split on _XP).
PROTEIN_ACC_LABELS = {
    # Pp - Physcomitrella patens (1)
    "XP_024392000": "PpDMP5-like",
    # Si - Setaria italica (14)
    "XP_004952503": "SiDMP6", "XP_004958514": "SiDMP1", "XP_004961228": "SiDMP4",
    "XP_004965304": "SiDMP2", "XP_004965424": "SiDMP4", "XP_004966913": "SiDMP10",
    "XP_004968801": "SiDMP3", "XP_004969121": "SiDMP3", "XP_004969493": "SiDMP6",
    "XP_004970780": "SiDMP7", "XP_004971846": "SiDMP2", "XP_004972417": "SiDMP3",
    "XP_004984268": "SiDMP2", "XP_022681673": "SiDMP9-like",
    # Ma - Musa acuminata (42)
    "XP_009380593": "MaDMP6-like", "XP_009380594": "MaDMP4",
    "XP_009381725": "MaDMP8-like", "XP_009385032": "MaDMP2",
    "XP_009387184": "MaDMP5-like", "XP_009394826": "MaDMP4",
    "XP_009395997": "MaDMP4",      "XP_009402533": "MaDMP3-like",
    "XP_009402791": "MaDMP7-like", "XP_009403823": "MaDMP2-like",
    "XP_009408744": "MaDMP3-like", "XP_009410685": "MaDMP3",
    "XP_009413807": "MaDMP4",      "XP_064937510": "MaDMP4-like",
    "XP_064937589": "MaDMP7-like", "XP_064943435": "MaDMP8-like",
    "XP_064943850": "MaDMP4-like", "XP_064943851": "MaDMP6-like",
    "XP_064944811": "MaDMP2-like", "XP_064958017": "MaDMP4-like",
    "XP_064968105": "MaDMP2-like", "XP_064972103": "MaDMP5-like",
    "XP_064979449": "MaDMP7-like", "XP_065001710": "MaDMP4-like",
    "XP_065003221": "MaDMP5-like", "XP_065007018": "MaDMP3-like",
    "XP_065007094": "MaDMP7-like", "XP_065012787": "MaDMP6-like",
    "XP_065015983": "MaDMP3-like", "XP_065016215": "MaDMP5-like",
    "XP_065017841": "MaDMP7-like", "XP_065019000": "MaDMP4-like",
    "XP_065027107": "MaDMP4-like", "XP_065027109": "MaDMP6-like",
    "XP_065027149": "MaDMP8-like", "XP_065036817": "MaDMP4-like",
    "XP_065036951": "MaDMP4-like", "XP_065037018": "MaDMP5-like",
    "XP_065041527": "MaDMP7-like", "XP_065041528": "MaDMP3-like",
    "XP_065042776": "MaDMP2-like", "XP_065049239": "MaDMP3-like",
    # Pt - Populus trichocarpa (11)
    "XP_002305111": "PtDMP4", "XP_002312247": "PtDMP8", "XP_002312376": "PtDMP2",
    "XP_002315049": "PtDMP8", "XP_002315542": "PtDMP2", "XP_006372758": "PtDMP7",
    "XP_006376329": "PtDMP2", "XP_006377414": "PtDMP4", "XP_006385029": "PtDMP3",
    "XP_006389536": "PtDMP3", "XP_024439153": "PtDMP10",
    # Bd - Brachypodium distachyon (13)
    "XP_003557584": "BdDMP3", "XP_003557899": "BdDMP2", "XP_003565771": "BdDMP2",
    "XP_003566981": "BdDMP4", "XP_003567905": "BdDMP6", "XP_003572244": "BdDMP2",
    "XP_003572527": "BdDMP6", "XP_010228847": "BdDMP2", "XP_010236366": "BdDMP2",
    "XP_024314432": "BdDMP10-like-X1", "XP_024314433": "BdDMP10-like-X2",
    "XP_024316011": "BdDMP7", "XP_024316012": "BdDMP7",
    # Gr - Gossypium raimondii (9)
    "XP_012453761": "GrDMP2", "XP_012457653": "GrDMP4", "XP_012457725": "GrDMP2",
    "XP_012458528": "GrDMP7", "XP_012467832": "GrDMP3", "XP_012478015": "GrDMP10",
    "XP_012481184": "GrDMP2", "XP_012487405": "GrDMP2", "XP_012489581": "GrDMP9",
    # Ca - Capsicum annuum (7)
    "XP_016539005": "CaDMP3",      "XP_016540695": "CaDMP7-like",
    "XP_016557530": "CaDMP2",      "XP_016571689": "CaDMP9-like",
    "XP_016582229": "CaDMP10",     "XP_016582285": "CaDMP3-like",
    "XP_047261276": "CaDMP6-like",
    # Ib - Ipomoea batatas uses BSXM accessions only (no XP_*); Query IDs
    # of the form 'IbDMP3_BSXM01000036.1' resolve via the regex split below.
}


def _extract_species(query_id):
    """Return two-letter species prefix, e.g. 'Sl' from 'SlDMP3_XP_...'.

    Handles raw-ID styles seen in the curated CSV. The XP_* and lcl|NC_*_cds_XP_*
    branches resolve via PROTEIN_ACC_LABELS to recover the species code (e.g.
    'XP_024392000' → 'PpDMP5-like' → 'Pp') for the older 8-species panel under
    II_INPUTS/DMP_query_fasta_file/.
    """
    qid = query_id.strip().removeprefix("lcl|")
    if qid.startswith(("CDX", "CDY")) or qid.startswith("mRNA.Bju"):
        return "Br"
    if qid.startswith("XM_"):
        m = re.match(r"(XM_\d+)", qid)
        if m and m.group(1) in ACCESSION_LABELS:
            return _extract_species(ACCESSION_LABELS[m.group(1)])
    # NCBI secondary-token form: NC_<contig>.<v>_cds_XP_<acc>_<n>
    m = re.match(r"NC_\d+\.\d+_cds_(XP_\d+)", qid)
    if m and m.group(1) in PROTEIN_ACC_LABELS:
        return _extract_species(PROTEIN_ACC_LABELS[m.group(1)])
    # Bare protein accession: 'XP_024392000.1' or 'XP_024392000_*'
    if qid.startswith("XP_"):
        m = re.match(r"(XP_\d+)", qid)
        if m and m.group(1) in PROTEIN_ACC_LABELS:
            return _extract_species(PROTEIN_ACC_LABELS[m.group(1)])
    m = re.match(r"^([A-Z][a-z]+)", qid)
    return m.group(1) if m else None


def _short_label(query_id):
    """Return a publication-ready label for a BLAST query ID.

    Resolves five raw-ID styles, falling back to the prefix-strip regex:
      'lcl|XM_003621193.2_cds_XP_003621241.1_1'    → 'MtDMP9'        (ACCESSION_LABELS, XM key)
      'lcl|NC_037265.1_cds_XP_024392000.1_26917'   → 'PpDMP5-like'   (PROTEIN_ACC_LABELS, secondary-token form)
      'XP_024392000.1'                              → 'PpDMP5-like'   (PROTEIN_ACC_LABELS, bare-XP form)
      'CDX74441.'                                   → 'CDX74441'      (strip trailing '.')
      'mRNA.BjuA04g10430S.'                         → 'BjuA04g10430S' (strip 'mRNA.' prefix and dot)
      'PpDMP5-like_XP_024392000.1_26917'            → 'PpDMP5-like'   (regex split on '_XP')
      'AtDMP8+AtDMP9'                               → 'AtDMP8+AtDMP9' (pass-through)
    """
    qid = query_id.strip().removeprefix("lcl|")
    qid = qid.rstrip(".")  # CDX74441. → CDX74441

    if qid.startswith("XM_"):
        m = re.match(r"(XM_\d+)", qid)
        if m:
            stem = m.group(1)
            if stem in ACCESSION_LABELS:
                return ACCESSION_LABELS[stem]
            return stem  # fallback: bare accession without _cds_XP_* tail

    # NCBI secondary-token form (e.g. headers retain lcl|NC_*_cds_XP_* as the
    # second whitespace token; appears as Query ID if the first token is dropped).
    m = re.match(r"NC_\d+\.\d+_cds_(XP_\d+)", qid)
    if m and m.group(1) in PROTEIN_ACC_LABELS:
        return PROTEIN_ACC_LABELS[m.group(1)]

    # Bare protein accession.
    if qid.startswith("XP_"):
        m = re.match(r"(XP_\d+)", qid)
        if m and m.group(1) in PROTEIN_ACC_LABELS:
            return PROTEIN_ACC_LABELS[m.group(1)]

    if qid.startswith("mRNA.Bju"):
        return qid[len("mRNA."):]  # mRNA.BjuA04g10430S → BjuA04g10430S

    label = re.split(r"[_\s](?:XP|XM|BSXM|NM|AT[0-9]|Glyma\.|OZ[0-9]|NC_[0-9]|LOC[0-9])", qid)[0]
    return label.rstrip("_.")


def _load_plant_csv(path):
    """Load the plant-only CSV and attach helper columns."""
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    df["Species"]    = df["Query ID"].apply(_extract_species)
    df["ShortLabel"] = df["Query ID"].apply(_short_label)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 - Cross-species % Identity Heatmap
# ─────────────────────────────────────────────────────────────────────────────

# ── Shared helpers for heatmaps ───────────────────────────────────────────────

SHORT_FAMILY = {
    "Solanaceae":     "Solanaceae",
    "Convolvulaceae": "Conv.",
    "Cucurbitaceae":  "Cucurb.",
    "Salicaceae":     "Salic.",
    "Malvaceae":      "Malv.",
    "Brassicaceae":   "Brassicaceae",
    "Fabaceae":       "Fabaceae",
    "Musaceae":       "Musa.",
}

def _txt_col(hex_color):
    """Return 'white' or '#1a1a2e' based on WCAG relative luminance."""
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i+2], 16) / 255 for i in (0, 2, 4))
    def _lin(c): return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    lum = 0.2126 * _lin(r) + 0.7152 * _lin(g) + 0.0722 * _lin(b)
    return "white" if lum < 0.40 else "#1a1a2e"

def _cell_text_color(val):
    """Return appropriate text color for a heatmap cell value."""
    if val >= 91:
        return "white"
    if val >= 85:
        return "#1a1a2e"
    return "#1a1a2e"

def _draw_family_bar(ax, ordered_species, pad=0.08, rounding=0.15):
    """Draw rounded family annotation badges on *ax*."""
    ax.set_xlim(0, len(ordered_species))
    ax.set_ylim(0, 1)
    ax.axis("off")

    for (family, members), fcolor in zip(FAMILY_GROUPS, FAMILY_COLORS):
        idxs = [i for i, s in enumerate(ordered_species) if s in members]
        if not idxs:
            continue
        x0, x1 = min(idxs), max(idxs) + 1
        badge = FancyBboxPatch(
            (x0 + pad, 0.10), (x1 - x0) - 2 * pad, 0.80,
            boxstyle=f"round,pad={rounding}",
            facecolor=fcolor, edgecolor="white", linewidth=1.8,
        )
        ax.add_patch(badge)
        fname = SHORT_FAMILY.get(family, family)
        txt = ax.text(
            (x0 + x1) / 2, 0.50, fname,
            ha="center", va="center",
            fontsize=7.8 if len(fname) > 6 else 9.0,
            fontweight="bold",
            color=_txt_col(fcolor),
        )
        txt.set_path_effects([
            mpe.withStroke(linewidth=0.3, foreground=_txt_col(fcolor)),
        ])

def _draw_colorbar(fig, label, cfg: "VizConfig" = None, rect=None):
    """Draw a polished horizontal colorbar at the bottom of *fig*."""
    if cfg is None:
        cfg = VizConfig()
    if rect is None:
        rect = [0.21, 0.055, 0.65, 0.020]
    cbar_ax = fig.add_axes(rect)
    norm = mcolors.Normalize(vmin=cfg.heatmap_vmin, vmax=cfg.heatmap_vmax)
    cb = mcolorbar.ColorbarBase(
        cbar_ax, cmap=plt.get_cmap(cfg.colormap),
        norm=norm, orientation="horizontal",
    )
    cb.set_label(label, fontsize=10, labelpad=5)
    cb.ax.tick_params(labelsize=8.5, length=3)
    cb.outline.set_linewidth(0.6)
    return cb


def plot_heatmap(df, out_dir, cfg: VizConfig = None):
    if cfg is None:
        cfg = VizConfig()
    from matplotlib.gridspec import GridSpec
    print("  Generating Fig 1: Cross-species % Identity Heatmap …")

    # ── Data prep ─────────────────────────────────────────────────────────────
    pivot = (
        df.groupby(["Subject ID", "Species"])["Percent Identity"]
        .max()
        .unstack()
    )
    ordered_species = [s for s in SPECIES_ORDER if s in pivot.columns]
    ordered_species += sorted(s for s in pivot.columns if s not in set(SPECIES_ORDER))
    ordered_genes   = [g for g in GENE_ORDER if g in pivot.index]
    pivot = pivot.reindex(index=ordered_genes, columns=ordered_species)

    row_labels = [GENE_SHORT.get(g, g)     for g in pivot.index]
    col_labels = [SPECIES_LABELS.get(s, s) for s in pivot.columns]

    n_sp   = len(ordered_species)
    n_gene = len(ordered_genes)
    fig_w  = max(14, n_sp * cfg.heatmap_w_scale + 3.5)
    fig_h  = max(6.0, n_gene * cfg.heatmap_h_scale + 3.5)

    # ── Figure: annotation bar (top) + heatmap (bottom) ──────────────────────
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    gs  = GridSpec(
        2, 1, figure=fig,
        height_ratios=[0.10, 1],
        hspace=0.04,
        left=0.18, right=0.88,
        top=0.89, bottom=0.17,
    )
    ax_ann  = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1])

    # ── Heatmap ───────────────────────────────────────────────────────────────
    cmap = plt.get_cmap(cfg.colormap)
    sns.heatmap(
        pivot, ax=ax_heat,
        cmap=cmap, vmin=cfg.heatmap_vmin, vmax=cfg.heatmap_vmax,
        mask=pivot.isna(),
        linewidths=1.2, linecolor="white",
        annot=False, cbar=False,
    )

    # Subtle no-data cells with soft hatching
    for ri, gene in enumerate(pivot.index):
        for ci, sp in enumerate(pivot.columns):
            if pd.isna(pivot.loc[gene, sp]):
                ax_heat.add_patch(mpatches.Rectangle(
                    (ci + 0.02, ri + 0.02), 0.96, 0.96,
                    facecolor="#eef0f2", edgecolor="#d5d8dc",
                    linewidth=0.4, zorder=2,
                ))

    # Cell text with subtle shadow for readability
    for ri, gene in enumerate(pivot.index):
        for ci, sp in enumerate(pivot.columns):
            val = pivot.loc[gene, sp]
            if pd.isna(val):
                ax_heat.text(
                    ci + 0.5, ri + 0.5, "n/d",
                    ha="center", va="center",
                    fontsize=7.5, color="#9ca3af", style="italic", zorder=3,
                )
            else:
                tc = _cell_text_color(val)
                txt = ax_heat.text(
                    ci + 0.5, ri + 0.5, f"{val:.1f}",
                    ha="center", va="center", fontsize=9.5, zorder=3,
                    color=tc,
                    fontweight="bold" if val >= 88 else "medium",
                )
                if tc == "white":
                    txt.set_path_effects([
                        mpe.withStroke(linewidth=1.5, foreground="#00000033"),
                    ])

    ax_heat.set_xticklabels(
        col_labels, rotation=45, ha="right", fontsize=10, style="italic"
    )
    ax_heat.set_yticklabels(row_labels, rotation=0, fontsize=10, fontweight="medium")
    ax_heat.set_xlabel("")
    ax_heat.set_ylabel("")
    ax_heat.tick_params(left=False, bottom=False, pad=4)

    # ── Family annotation bar ─────────────────────────────────────────────────
    _draw_family_bar(ax_ann, ordered_species)

    # ── Colorbar ──────────────────────────────────────────────────────────────
    _draw_colorbar(fig, "Max % Identity (%)", cfg=cfg)

    # ── Title ─────────────────────────────────────────────────────────────────
    fig.suptitle(
        f"BLASTn Percent Identity: Smel{_GENE_GROUP} Genes vs. Plant Orthologs",
        fontsize=14, fontweight="bold", y=0.97,
        color="#1a1a2e",
    )

    _save(fig, out_dir, "fig1_blastn_identity_heatmap", save_dpi=cfg.save_dpi)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1b - Heatmap with % Identity colour + E-value annotation
# ─────────────────────────────────────────────────────────────────────────────

def _fmt_ev(ev):
    """Format E-value as a compact exponent string.

    0          → '0'
    1.15e-176  → 'e-176'
    5.43e-45   → 'e-45'
    """
    if ev == 0:
        return "0"
    import math
    return f"e{int(math.floor(math.log10(ev)))}"


def plot_heatmap_evalue(df, out_dir, cfg: VizConfig = None):
    """Heatmap with % Identity colour scale and E-value exponent in each cell."""
    if cfg is None:
        cfg = VizConfig()
    from matplotlib.gridspec import GridSpec
    print("  Generating Fig 1b: Heatmap with % Identity + E-value …")

    # Best hit per (Subject Gene × Species): highest % identity; record its E-value
    idx_max  = df.groupby(["Subject ID", "Species"])["Percent Identity"].idxmax()
    best     = df.loc[idx_max].set_index(["Subject ID", "Species"])
    pivot_id = best["Percent Identity"].unstack()
    pivot_ev = best["E-value"].unstack()

    ordered_species = [s for s in SPECIES_ORDER if s in pivot_id.columns]
    ordered_species += sorted(s for s in pivot_id.columns if s not in set(SPECIES_ORDER))
    ordered_genes   = [g for g in GENE_ORDER if g in pivot_id.index]
    pivot_id = pivot_id.reindex(index=ordered_genes, columns=ordered_species)
    pivot_ev = pivot_ev.reindex(index=ordered_genes, columns=ordered_species)

    row_labels = [GENE_SHORT.get(g, g)     for g in pivot_id.index]
    col_labels = [SPECIES_LABELS.get(s, s) for s in pivot_id.columns]

    n_sp   = len(ordered_species)
    n_gene = len(ordered_genes)
    fig_w  = max(14, n_sp * cfg.heatmap_w_scale + 3.5)
    fig_h  = max(7.0, n_gene * (cfg.heatmap_h_scale * 1.20) + 3.5)

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    gs  = GridSpec(
        2, 1, figure=fig,
        height_ratios=[0.10, 1],
        hspace=0.04,
        left=0.18, right=0.88,
        top=0.89, bottom=0.17,
    )
    ax_ann  = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1])

    # ── Heatmap (colour = % identity) ────────────────────────────────────────
    cmap = plt.get_cmap(cfg.colormap)
    sns.heatmap(
        pivot_id, ax=ax_heat,
        cmap=cmap, vmin=cfg.heatmap_vmin, vmax=cfg.heatmap_vmax,
        mask=pivot_id.isna(),
        linewidths=1.2, linecolor="white",
        annot=False, cbar=False,
    )

    # No-data cells
    for ri, gene in enumerate(pivot_id.index):
        for ci, sp in enumerate(pivot_id.columns):
            if pd.isna(pivot_id.loc[gene, sp]):
                ax_heat.add_patch(mpatches.Rectangle(
                    (ci + 0.02, ri + 0.02), 0.96, 0.96,
                    facecolor="#eef0f2", edgecolor="#d5d8dc",
                    linewidth=0.4, zorder=2,
                ))

    # Two-line cell annotation: % identity (upper) + E-value exponent (lower)
    for ri, gene in enumerate(pivot_id.index):
        for ci, sp in enumerate(pivot_id.columns):
            val = pivot_id.loc[gene, sp]
            if pd.isna(val):
                ax_heat.text(
                    ci + 0.5, ri + 0.5, "n/d",
                    ha="center", va="center",
                    fontsize=7.5, color="#9ca3af", style="italic", zorder=3,
                )
            else:
                ev      = pivot_ev.loc[gene, sp]
                tc      = _cell_text_color(val)
                ev_tc   = "#ffffffbb" if val >= 91 else "#6b7280"
                # % identity: upper half of cell
                id_txt = ax_heat.text(
                    ci + 0.5, ri + 0.36, f"{val:.1f}",
                    ha="center", va="center", fontsize=9.0, zorder=3,
                    color=tc,
                    fontweight="bold" if val >= 88 else "medium",
                )
                if tc == "white":
                    id_txt.set_path_effects([
                        mpe.withStroke(linewidth=1.5, foreground="#00000033"),
                    ])
                # E-value exponent: lower half of cell
                ax_heat.text(
                    ci + 0.5, ri + 0.67, _fmt_ev(ev),
                    ha="center", va="center", fontsize=7.0, zorder=3,
                    color=ev_tc, style="italic",
                )

    ax_heat.set_xticklabels(
        col_labels, rotation=45, ha="right", fontsize=10, style="italic"
    )
    ax_heat.set_yticklabels(row_labels, rotation=0, fontsize=10, fontweight="medium")
    ax_heat.set_xlabel("")
    ax_heat.set_ylabel("")
    ax_heat.tick_params(left=False, bottom=False, pad=4)

    # ── Family annotation bar ─────────────────────────────────────────────────
    _draw_family_bar(ax_ann, ordered_species)

    # ── Colorbar + legend note ───────────────────────────────────────────────
    cb = _draw_colorbar(
        fig,
        "Max % Identity (%)     \u2502     cell bottom: E-value exponent  "
        "(e.g. e\u221250 = 1\u00d710\u207b\u2075\u2070)",
        cfg=cfg,
    )
    cb.set_label(
        "Max % Identity (%)     \u2502     cell bottom: E-value exponent  "
        "(e.g. e\u221250 = 1\u00d710\u207b\u2075\u2070)",
        fontsize=9.0, labelpad=5,
    )

    fig.suptitle(
        f"BLASTn Percent Identity + E-value: Smel{_GENE_GROUP} Genes vs. Plant Orthologs",
        fontsize=14, fontweight="bold", y=0.97,
        color="#1a1a2e",
    )

    _save(fig, out_dir, "fig1b_blastn_identity_evalue_heatmap", save_dpi=cfg.save_dpi)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 - Per-gene Top-Hits Lollipop (Bit Score)
# ─────────────────────────────────────────────────────────────────────────────

def _lollipop_row_count(gene_id, df, top_n):
    """Return the number of rows that the lollipop panel for *gene_id* will
    display: top-N + force-included haploid-inducer hits (DMP only).
    Mirrors the row-selection logic in plot_lollipop / _plot_lollipop_splits
    so grid heights can be sized proportionally before rendering."""
    gene_df = df[df["Subject ID"] == gene_id]
    if gene_df.empty:
        return 1
    top = (gene_df.sort_values("Bit Score", ascending=False)
                  .drop_duplicates(subset=["Query ID"]).head(top_n))
    if _GENE_GROUP != "DMP":
        return len(top)
    hi_rows = (gene_df[gene_df["ShortLabel"].isin(HAPLOID_INDUCER_LABELS)]
                      .drop_duplicates(subset=["Query ID"]))
    extra = hi_rows[~hi_rows["Query ID"].isin(set(top["Query ID"]))]
    return len(top) + len(extra)


def plot_lollipop(df, out_dir, top_n=10, cfg: VizConfig = None):
    """Top-N hits per gene ranked by Bit Score.

    Known haploid-inducer genes (AtDMP8, AtDMP9) are always shown even when
    they fall outside the top-N, marked with a gold stem and ★ prefix in the
    y-tick label.  Source: phylogenetic input FASTA
    (3_RESULT/DMP/04_MSA/.../input_fasta.fa).
    """
    print(f"  Generating Fig 3: Top-{top_n} Hits Lollipop (Bit Score) …")

    if cfg is None:
        cfg = VizConfig()
    HI_STEM_COLOR = cfg.hi_stem_color
    HI_EDGE_COLOR = "#9e7608"   # darker gold dot edge (fixed relative to stem)
    STEM_COLOR    = cfg.stem_color
    DOT_EDGE      = "#374151"   # dark gray dot edge

    n_genes = len(GENE_ORDER)
    ncols   = cfg.lollipop_ncols
    nrows   = math.ceil(n_genes / ncols)

    norm    = mcolors.Normalize(vmin=cfg.heatmap_vmin, vmax=cfg.heatmap_vmax)
    cmap    = plt.get_cmap(cfg.colormap)

    # Per-gene row count (same logic as the panel render loop), used to size
    # each grid row so every DMP gene gets the same vertical pixel space.
    _row_counts = [_lollipop_row_count(g, df, top_n) for g in GENE_ORDER]
    # Pad to fill the grid (legend slot uses the row average so it doesn't
    # look squashed).
    while len(_row_counts) < nrows * ncols:
        _row_counts.append(max(1, sum(_row_counts) // max(1, len(_row_counts))))
    # Per-grid-row max determines that row's height ratio.
    _height_ratios = [
        max(_row_counts[r * ncols : (r + 1) * ncols]) for r in range(nrows)
    ]
    _row_pixel_h = 0.32  # inches per gene row, tuned to match figures_v1
    fig_height   = sum(_height_ratios) * _row_pixel_h + 2.5
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(8.0 * ncols, fig_height),
        facecolor="white",
        gridspec_kw={"height_ratios": _height_ratios},
    )
    axes_flat = axes.flatten()

    for idx, gene in enumerate(GENE_ORDER):
        ax      = axes_flat[idx]
        gene_df = df[df["Subject ID"] == gene].copy()

        if gene_df.empty:
            clade = GENE_CLADE.get(gene, "")
            short = GENE_SHORT.get(gene, gene)
            title = f"{short}  ({clade})" if clade else short
            ax.set_title(
                title, fontsize=10.5, fontweight="bold",
                loc="left", color="#1a1a2e", pad=8,
            )
            ax.text(
                0.5, 0.5, "No plant ortholog hits",
                transform=ax.transAxes,
                ha="center", va="center",
                fontsize=10, color="#9ca3af", style="italic",
            )
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            continue

        # Top-N hits by Bit Score (deduplicated per Query ID)
        top = (
            gene_df
            .sort_values("Bit Score", ascending=False)
            .drop_duplicates(subset=["Query ID"])
            .head(top_n)
        )

        # Force-include known haploid-inducer hits even if outside top-N (DMP only)
        if _GENE_GROUP == "DMP":
            hi_rows = (
                gene_df[gene_df["ShortLabel"].isin(HAPLOID_INDUCER_LABELS)]
                .sort_values("Bit Score", ascending=False)
                .drop_duplicates(subset=["Query ID"])
            )
            already_in = set(top["Query ID"])
            hi_extra   = hi_rows[~hi_rows["Query ID"].isin(already_in)]
            sub = pd.concat([top, hi_extra], ignore_index=True)
        else:
            sub = top

        # Sort: lowest bit score at bottom (so best hit is at top of panel)
        sub = sub.sort_values("Bit Score", ascending=True).reset_index(drop=True)

        y_pos  = np.arange(len(sub))
        is_hi  = (sub["ShortLabel"].isin(HAPLOID_INDUCER_LABELS).values
                  if _GENE_GROUP == "DMP" else np.zeros(len(sub), dtype=bool))
        colors = [cmap(norm(v)) for v in sub["Percent Identity"]]

        # Subtle horizontal gridlines behind everything
        ax.set_axisbelow(True)
        ax.yaxis.grid(False)
        ax.xaxis.grid(True, color="#f3f4f6", linewidth=0.6, zorder=0)

        # Stems: gold dashed for HI genes, soft gray for others
        for yi, bs, hi in zip(y_pos, sub["Bit Score"], is_hi):
            ax.hlines(yi, 0, bs,
                      color=HI_STEM_COLOR if hi else STEM_COLOR,
                      linewidth=2.2 if hi else 1.4,
                      linestyle="--" if hi else "-",
                      zorder=1)

        # Dots: polished with subtle shadow
        for yi, bs, col, hi in zip(y_pos, sub["Bit Score"], colors, is_hi):
            # Shadow dot (slightly offset)
            ax.scatter(bs + 2, yi - 0.06,
                       c=["#00000015"],
                       s=(cfg.lollipop_dot_size_hi - 10) if hi else cfg.lollipop_dot_size,
                       zorder=2, edgecolors="none")
            # Main dot
            ax.scatter(bs, yi,
                       c=[col],
                       s=cfg.lollipop_dot_size_hi if hi else cfg.lollipop_dot_size,
                       zorder=3,
                       edgecolors=HI_EDGE_COLOR if hi else DOT_EDGE,
                       linewidths=1.8 if hi else 0.8)

        # Inline annotation: % identity | coverage (with white stroke for readability)
        x_max = sub["Bit Score"].max()
        for yi, (_, row) in zip(y_pos, sub.iterrows()):
            ann_txt = ax.text(
                row["Bit Score"] + x_max * 0.03, yi,
                f"{row['Percent Identity']:.0f}%  \u2502  {row['Query Coverage']:.0f}% cov",
                va="center", fontsize=7.5, color="#374151",
                fontweight="medium",
            )
            ann_txt.set_path_effects([
                mpe.withStroke(linewidth=2.5, foreground="white"),
            ])

        # Y-tick labels: prefix ★ and color gold for HI genes
        tick_labels = []
        tick_colors = []
        for lbl, hi in zip(sub["ShortLabel"], is_hi):
            tick_labels.append(f"\u2605 {lbl}" if hi else lbl)
            tick_colors.append(HI_STEM_COLOR if hi else "#374151")

        ax.set_yticks(y_pos)
        ax.set_yticklabels(tick_labels, fontsize=8.5)
        for tick_lbl, col in zip(ax.get_yticklabels(), tick_colors):
            tick_lbl.set_color(col)
            if col == HI_STEM_COLOR:
                tick_lbl.set_fontweight("bold")

        ax.set_xlabel("Bit Score", fontsize=9, color="#4b5563")
        ax.set_xlim(0, x_max * cfg.lollipop_x_pad)
        clade = GENE_CLADE.get(gene, "")
        short = GENE_SHORT.get(gene, gene)
        title = f"{short}  ({clade})" if clade else short
        ax.set_title(
            title, fontsize=10.5, fontweight="bold",
            loc="left", color="#1a1a2e", pad=8,
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#e5e7eb")
        ax.spines["bottom"].set_color("#e5e7eb")
        ax.tick_params(axis="x", labelsize=8.5, colors="#6b7280")
        ax.tick_params(axis="y", length=0)

    # Use empty panel slot(s) for a unified HI gene legend (DMP only)
    for idx in range(n_genes, len(axes_flat)):
        ax_empty = axes_flat[idx]
        ax_empty.axis("off")

        if _GENE_GROUP != "DMP":
            continue

        ax_empty.set_visible(True)

        # ── Unified haploid-inducer legend ──
        # Labels mirror the actual ShortLabels that appear in the panels (resolved
        # by _short_label + ACCESSION_LABELS). When a query FASTA was merged or
        # uses a non-standard header, the displayed label may differ from the
        # locus's literature name (e.g. SlDMP3 here = Solyc05g007920, the
        # validated haploid inducer per Deng et al., 2025).
        hi_genes_table = [
            ("AtDMP8+AtDMP9", "A. thaliana",     "Zhong et al., 2020"),
            ("NtDMP",         "N. tabacum",      "X. Zhang et al., 2022"),
            ("SlDMP3",        "S. lycopersicum", "Deng et al., 2025"),
            ("MtDMP9",        "M. truncatula",   "N. Wang et al., 2022"),
            ("GmDMP2",        "G. max",          "Zhong et al., 2024"),
        ]

        # Title
        ax_empty.text(
            0.05, 0.92, "\u2605  Known Haploid-Inducer DMP Genes",
            transform=ax_empty.transAxes, fontsize=10.5,
            fontweight="bold", color=HI_STEM_COLOR, va="top",
        )

        # Visual key: stem + dot sample
        ax_empty.plot([0.06, 0.14], [0.82, 0.82], color=HI_STEM_COLOR,
                      linewidth=2.2, linestyle="--", transform=ax_empty.transAxes,
                      clip_on=False, zorder=2)
        ax_empty.scatter([0.14], [0.82], c=[cmap(norm(80))], s=100, zorder=3,
                         edgecolors=HI_EDGE_COLOR, linewidths=1.8,
                         transform=ax_empty.transAxes, clip_on=False)
        ax_empty.text(
            0.17, 0.82, "Gold dashed stem  +  \u2605 prefix",
            transform=ax_empty.transAxes, fontsize=8.5, va="center", color="#374151",
        )

        # Table of HI genes
        header = f"{'Gene':<10}{'Species':<20}{'Reference'}"
        lines  = [header, "\u2500" * 48]
        for gene_name, species, ref in hi_genes_table:
            lines.append(f"{gene_name:<10}{species:<20}{ref}")

        ax_empty.text(
            0.05, 0.72, "\n".join(lines),
            transform=ax_empty.transAxes, fontsize=7.5,
            fontfamily="monospace", va="top", color="#374151",
            bbox=dict(boxstyle="round,pad=0.5", fc="#fffbeb",
                      ec=HI_STEM_COLOR, lw=1.0, alpha=0.95),
        )

    # Shared colorbar at bottom
    cbar_ax = fig.add_axes([0.30, 0.015, 0.40, 0.012])
    cb = mcolorbar.ColorbarBase(
        cbar_ax, cmap=cmap, norm=norm, orientation="horizontal"
    )
    cb.set_label("% Identity", fontsize=9.5)
    cb.ax.tick_params(labelsize=8.5, length=3)
    cb.outline.set_linewidth(0.6)

    _lollipop_title = f"Top-{top_n} BLASTn Hits per Smel{_GENE_GROUP} Gene  (ranked by Bit Score)"
    if _GENE_GROUP == "DMP":
        _lollipop_title += (
            "\n\u2605 Gold = known haploid-inducer DMP "
            "(AtDMP8+AtDMP9, NtDMP, SlDMP3, MtDMP9, GmDMP2)"
        )
    fig.suptitle(
        _lollipop_title,
        fontsize=12.5, fontweight="bold",
        color="#1a1a2e",
    )
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout(rect=[0, 0.05, 1, 0.95])

    _save(fig, out_dir, "fig3_blastn_lollipop", save_dpi=cfg.save_dpi)

    # Also emit split variants (split_1 = first ceil(N/2) panels with the
    # suptitle; split_2 = remaining panels + HI legend + colorbar).
    _plot_lollipop_splits(df, out_dir, top_n, cfg)


def _plot_lollipop_splits(df, out_dir, top_n, cfg):
    """Render fig3_blastn_lollipop_split_{1,2}.{png,svg}.

    Splits Fig 3 into two equal-height figures:
      Split 1: first ceil(N/2) panels with the suptitle.
      Split 2: remaining panels + HI legend + % identity colorbar.

    Styling mirrors plot_lollipop() exactly. The rendering loop is duplicated
    to avoid touching the working combined-figure code path.
    """
    HI_STEM_COLOR = cfg.hi_stem_color
    HI_EDGE_COLOR = "#9e7608"
    STEM_COLOR    = cfg.stem_color
    DOT_EDGE      = "#374151"

    norm = mcolors.Normalize(vmin=cfg.heatmap_vmin, vmax=cfg.heatmap_vmax)
    cmap = plt.get_cmap(cfg.colormap)
    row_pixel_h = 0.32

    genes = list(GENE_ORDER)
    half  = math.ceil(len(genes) / 2)
    parts = [genes[:half], genes[half:]]

    for split_idx, gene_subset in enumerate(parts, start=1):
        ncols = 2
        nrows = math.ceil(max(len(gene_subset), 4) / ncols)
        is_last = split_idx == len(parts)

        # Per-grid-row max row count -> height ratio so every DMP gene gets
        # the same vertical pixel space across the split figure.
        row_counts = [_lollipop_row_count(g, df, top_n) for g in gene_subset]
        while len(row_counts) < nrows * ncols:
            row_counts.append(max(1, sum(row_counts) // max(1, len(row_counts))))
        height_ratios = [
            max(row_counts[r * ncols : (r + 1) * ncols]) for r in range(nrows)
        ]
        fig_height = sum(height_ratios) * row_pixel_h + 2.5

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(8.0 * ncols, fig_height),
            facecolor="white",
            gridspec_kw={"height_ratios": height_ratios},
        )
        axes_flat = axes.flatten()

        for idx, gene in enumerate(gene_subset):
            ax = axes_flat[idx]
            gene_df = df[df["Subject ID"] == gene].copy()

            if gene_df.empty:
                clade = GENE_CLADE.get(gene, "")
                short = GENE_SHORT.get(gene, gene)
                title = f"{short}  ({clade})" if clade else short
                ax.set_title(title, fontsize=10.5, fontweight="bold",
                             loc="left", color="#1a1a2e", pad=8)
                ax.text(0.5, 0.5, "No plant ortholog hits",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=10, color="#9ca3af", style="italic")
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                continue

            top = (gene_df.sort_values("Bit Score", ascending=False)
                          .drop_duplicates(subset=["Query ID"]).head(top_n))

            if _GENE_GROUP == "DMP":
                hi_rows = (gene_df[gene_df["ShortLabel"].isin(HAPLOID_INDUCER_LABELS)]
                                  .sort_values("Bit Score", ascending=False)
                                  .drop_duplicates(subset=["Query ID"]))
                already_in = set(top["Query ID"])
                hi_extra   = hi_rows[~hi_rows["Query ID"].isin(already_in)]
                sub = pd.concat([top, hi_extra], ignore_index=True)
            else:
                sub = top

            sub = sub.sort_values("Bit Score", ascending=True).reset_index(drop=True)

            y_pos  = np.arange(len(sub))
            is_hi  = (sub["ShortLabel"].isin(HAPLOID_INDUCER_LABELS).values
                      if _GENE_GROUP == "DMP" else np.zeros(len(sub), dtype=bool))
            colors = [cmap(norm(v)) for v in sub["Percent Identity"]]

            ax.set_axisbelow(True)
            ax.yaxis.grid(False)
            ax.xaxis.grid(True, color="#f3f4f6", linewidth=0.6, zorder=0)

            for yi, bs, hi in zip(y_pos, sub["Bit Score"], is_hi):
                ax.hlines(yi, 0, bs,
                          color=HI_STEM_COLOR if hi else STEM_COLOR,
                          linewidth=2.2 if hi else 1.4,
                          linestyle="--" if hi else "-", zorder=1)

            for yi, bs, col, hi in zip(y_pos, sub["Bit Score"], colors, is_hi):
                ax.scatter(bs + 2, yi - 0.06, c=["#00000015"],
                           s=(cfg.lollipop_dot_size_hi - 10) if hi else cfg.lollipop_dot_size,
                           zorder=2, edgecolors="none")
                ax.scatter(bs, yi, c=[col],
                           s=cfg.lollipop_dot_size_hi if hi else cfg.lollipop_dot_size,
                           zorder=3,
                           edgecolors=HI_EDGE_COLOR if hi else DOT_EDGE,
                           linewidths=1.8 if hi else 0.8)

            x_max = sub["Bit Score"].max()
            for yi, (_, row) in zip(y_pos, sub.iterrows()):
                ann_txt = ax.text(
                    row["Bit Score"] + x_max * 0.03, yi,
                    f"{row['Percent Identity']:.0f}%  │  {row['Query Coverage']:.0f}% cov",
                    va="center", fontsize=7.5, color="#374151",
                    fontweight="medium")
                ann_txt.set_path_effects([
                    mpe.withStroke(linewidth=2.5, foreground="white")])

            tick_labels = [f"★ {lbl}" if hi else lbl
                           for lbl, hi in zip(sub["ShortLabel"], is_hi)]
            tick_colors = [HI_STEM_COLOR if hi else "#374151" for hi in is_hi]
            ax.set_yticks(y_pos)
            ax.set_yticklabels(tick_labels, fontsize=8.5)
            for tick_lbl, col in zip(ax.get_yticklabels(), tick_colors):
                tick_lbl.set_color(col)
                if col == HI_STEM_COLOR:
                    tick_lbl.set_fontweight("bold")

            ax.set_xlabel("Bit Score", fontsize=9, color="#4b5563")
            ax.set_xlim(0, x_max * cfg.lollipop_x_pad)
            clade = GENE_CLADE.get(gene, "")
            short = GENE_SHORT.get(gene, gene)
            title = f"{short}  ({clade})" if clade else short
            ax.set_title(title, fontsize=10.5, fontweight="bold",
                         loc="left", color="#1a1a2e", pad=8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_color("#e5e7eb")
            ax.spines["bottom"].set_color("#e5e7eb")
            ax.tick_params(axis="x", labelsize=8.5, colors="#6b7280")
            ax.tick_params(axis="y", length=0)

        for idx in range(len(gene_subset), len(axes_flat)):
            ax_e = axes_flat[idx]
            ax_e.axis("off")
            if not (is_last and _GENE_GROUP == "DMP"):
                continue
            ax_e.set_visible(True)
            hi_genes_table = [
                ("AtDMP8+AtDMP9", "A. thaliana",     "Zhong et al., 2020"),
                ("NtDMP",         "N. tabacum",      "X. Zhang et al., 2022"),
                ("SlDMP3",        "S. lycopersicum", "Deng et al., 2025"),
                ("MtDMP9",        "M. truncatula",   "N. Wang et al., 2022"),
                ("GmDMP2",        "G. max",          "Zhong et al., 2024"),
            ]
            ax_e.text(0.05, 0.92, "★  Known Haploid-Inducer DMP Genes",
                      transform=ax_e.transAxes, fontsize=10.5,
                      fontweight="bold", color=HI_STEM_COLOR, va="top")
            ax_e.plot([0.06, 0.14], [0.82, 0.82], color=HI_STEM_COLOR,
                      linewidth=2.2, linestyle="--", transform=ax_e.transAxes,
                      clip_on=False, zorder=2)
            ax_e.scatter([0.14], [0.82], c=[cmap(norm(80))], s=100, zorder=3,
                         edgecolors=HI_EDGE_COLOR, linewidths=1.8,
                         transform=ax_e.transAxes, clip_on=False)
            ax_e.text(0.17, 0.82, "Gold dashed stem  +  ★ prefix",
                      transform=ax_e.transAxes, fontsize=8.5, va="center",
                      color="#374151")
            header = f"{'Gene':<10}{'Species':<20}{'Reference'}"
            lines  = [header, "─" * 48]
            for gene_name, species, ref in hi_genes_table:
                lines.append(f"{gene_name:<10}{species:<20}{ref}")
            ax_e.text(0.05, 0.72, "\n".join(lines),
                      transform=ax_e.transAxes, fontsize=7.5,
                      fontfamily="monospace", va="top", color="#374151",
                      bbox=dict(boxstyle="round,pad=0.5", fc="#fffbeb",
                                ec=HI_STEM_COLOR, lw=1.0, alpha=0.95))

        if is_last:
            cbar_ax = fig.add_axes([0.30, 0.015, 0.40, 0.012])
            cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm,
                                        orientation="horizontal")
            cb.set_label("% Identity", fontsize=9.5)
            cb.ax.tick_params(labelsize=8.5, length=3)
            cb.outline.set_linewidth(0.6)

        if split_idx == 1:
            title = f"Top-{top_n} BLASTn Hits per Smel{_GENE_GROUP} Gene  (ranked by Bit Score)"
            if _GENE_GROUP == "DMP":
                title += ("\n★ Gold = known haploid-inducer DMP "
                          "(AtDMP8+AtDMP9, NtDMP, SlDMP3, MtDMP9, GmDMP2)")
            fig.suptitle(title, fontsize=12.5, fontweight="bold", color="#1a1a2e")
            rect = [0, 0.02, 1, 0.95]
        else:
            rect = [0, 0.05, 1, 1.0]

        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig.tight_layout(rect=rect)

        _save(fig, out_dir, f"fig3_blastn_lollipop_split_{split_idx}",
              save_dpi=cfg.save_dpi)


# ─────────────────────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────────────────────

def _save(fig, out_dir, stem, save_dpi=300):
    for ext in ("svg", "png"):
        p = out_dir / f"{stem}.{ext}"
        fig.savefig(p, dpi=save_dpi, bbox_inches="tight")
        print(f"    Saved: {p.name}")
    plt.close(fig)


def _discover_plant_csv(results_dir):
    hits = sorted(Path(results_dir).glob("*_plant_only.csv"))
    if not hits:
        sys.exit(f"ERROR: No *_plant_only.csv found in {results_dir}")
    if len(hits) > 1:
        # Pick the CSV with the most data rows (not just the newest by name)
        best = max(hits, key=lambda p: sum(1 for _ in open(p, encoding="utf-8")))
        print(f"  WARNING: Multiple plant CSVs found; using largest: {best.name}")
        return best
    return hits[0]


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    _script_dir     = Path(__file__).resolve().parent
    _pipeline_root  = _script_dir.parent.parent
    _default_results = (
        _pipeline_root
        / "III_RESULT" / "DMP" / "02_BLAST_Alignment"
        / "GPE001970_SMEL5" / "curated_results"
    )

    parser = argparse.ArgumentParser(
        description="Generate BLASTn visualisation figures (heatmap, lollipop)."
    )
    parser.add_argument(
        "--results-dir", default=str(_default_results),
        help="Directory containing *_plant_only.csv (default: auto-detected from script location)",
    )
    parser.add_argument(
        "--plant-csv", default=None,
        help="Explicit path to *_plant_only.csv (overrides --results-dir)",
    )
    parser.add_argument(
        "--gene-group", default=None,
        help=(
            "Gene group name (e.g. DMP, HAP2).  When provided, overrides path-based "
            "auto-detection and ensures the correct gene metadata is applied."
        ),
    )
    parser.add_argument(
        "--out-dir", default=None,
        help="Output directory for figures (default: <results-dir>/figures/)",
    )
    parser.add_argument(
        "--top-n", type=int, default=10,
        help="Top-N hits per gene in Fig 3 lollipop (default: 10)",
    )
    parser.add_argument(
        "--figures", default="heatmap,heatmap_evalue,lollipop",
        help=(
            "Comma-separated list of figures to generate "
            "(default: heatmap,heatmap_evalue,lollipop). "
            "Valid values: heatmap, heatmap_evalue, lollipop."
        ),
    )
    # ── Visualization configuration ────────────────────────────────────────────
    parser.add_argument(
        "--colormap", default="RdYlGn",
        help="Matplotlib colormap for heatmaps and lollipop colorbar (default: RdYlGn)",
    )
    parser.add_argument(
        "--figure-dpi", type=int, default=150,
        help="Screen/preview DPI for figures (default: 150)",
    )
    parser.add_argument(
        "--save-dpi", type=int, default=300,
        help="Export DPI for PNG output (default: 300)",
    )
    parser.add_argument(
        "--heatmap-vmin", type=float, default=65.0,
        help="Lower bound of percent identity color scale (default: 65)",
    )
    parser.add_argument(
        "--heatmap-vmax", type=float, default=100.0,
        help="Upper bound of percent identity color scale (default: 100)",
    )
    parser.add_argument(
        "--heatmap-w-scale", type=float, default=1.20,
        help="Figure width per species column in inches (default: 1.20)",
    )
    parser.add_argument(
        "--heatmap-h-scale", type=float, default=0.92,
        help="Figure height per gene row in inches (default: 0.92)",
    )
    parser.add_argument(
        "--lollipop-ncols", type=int, default=2,
        help="Number of columns in the gene panel grid (default: 2)",
    )
    parser.add_argument(
        "--lollipop-x-pad", type=float, default=1.60,
        help="X-axis limit = max_bit_score × this factor (default: 1.60)",
    )
    parser.add_argument(
        "--lollipop-dot-size", type=int, default=100,
        help="Dot marker area (pt²) for regular hits (default: 100)",
    )
    parser.add_argument(
        "--lollipop-dot-size-hi", type=int, default=150,
        help="Dot marker area for haploid-inducer hits (default: 150)",
    )
    parser.add_argument(
        "--hi-stem-color", default="#c7920a",
        help="Stem color for haploid-inducer hits (default: #c7920a gold)",
    )
    parser.add_argument(
        "--stem-color", default="#d1d5db",
        help="Stem color for regular hits (default: #d1d5db gray)",
    )
    parser.add_argument(
        "--hi-labels", default=None,
        help=(
            "Comma-separated haploid-inducer ShortLabels (override the built-in default). "
            "Authoritative source: [blast_visualize].haploid_inducer_labels in the DMP TOML."
        ),
    )
    args = parser.parse_args()

    if args.hi_labels is not None:
        new_labels = {s.strip() for s in args.hi_labels.split(",") if s.strip()}
        if new_labels:
            global HAPLOID_INDUCER_LABELS
            HAPLOID_INDUCER_LABELS = new_labels

    figures = {f.strip() for f in args.figures.split(",") if f.strip()}

    cfg = VizConfig(
        colormap=args.colormap,
        figure_dpi=args.figure_dpi,
        save_dpi=args.save_dpi,
        heatmap_vmin=args.heatmap_vmin,
        heatmap_vmax=args.heatmap_vmax,
        heatmap_w_scale=args.heatmap_w_scale,
        heatmap_h_scale=args.heatmap_h_scale,
        lollipop_ncols=args.lollipop_ncols,
        lollipop_x_pad=args.lollipop_x_pad,
        lollipop_dot_size=args.lollipop_dot_size,
        lollipop_dot_size_hi=args.lollipop_dot_size_hi,
        hi_stem_color=args.hi_stem_color,
        stem_color=args.stem_color,
    )
    plt.rcParams["figure.dpi"] = cfg.figure_dpi

    plant_csv = (
        Path(args.plant_csv) if args.plant_csv
        else _discover_plant_csv(args.results_dir)
    )
    if not plant_csv.exists():
        sys.exit(f"ERROR: CSV not found: {plant_csv}")

    out_dir = Path(args.out_dir) if args.out_dir else plant_csv.parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input CSV : {plant_csv}")
    print(f"Output dir: {out_dir}")
    print(f"Figures   : {', '.join(sorted(figures))}\n")

    df = _load_plant_csv(plant_csv)
    # Prefer explicitly passed --gene-group; fall back to path-based detection.
    gene_group = args.gene_group or _detect_gene_group(plant_csv)
    df = _init_gene_metadata(df, gene_group)
    n_genes = df["Subject ID"].nunique()
    print(f"Gene group: {_GENE_GROUP}")
    print(f"Loaded {len(df)} hits across {n_genes} subject gene(s)\n")

    if "heatmap" in figures:
        plot_heatmap(df, out_dir, cfg)
    if "heatmap_evalue" in figures:
        plot_heatmap_evalue(df, out_dir, cfg)
    if "lollipop" in figures:
        plot_lollipop(df, out_dir, top_n=args.top_n, cfg=cfg)

    print(f"\nAll figures saved to: {out_dir}/")


if __name__ == "__main__":
    main()
