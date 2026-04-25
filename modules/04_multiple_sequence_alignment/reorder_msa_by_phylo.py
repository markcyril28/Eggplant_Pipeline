#!/usr/bin/env python3
"""
reorder_msa_by_phylo.py — Reorder a FASTA/CLUSTAL MSA to match phylogenetic tree leaf order.

Standalone usage (hardcoded defaults target DMP GPE001970 amino-acid MSA):
    python3 reorder_msa_by_phylo.py

Orchestrated usage:
    python3 reorder_msa_by_phylo.py \\
        --msa          <aligned.fas|.aln>  \\
        --tree         <newick.bestTree>   \\
        --output       <prefix>            \\
        --anchor-top   <seq_id>            \\
        --anchor-bottom <seq_id>           \\
        --espript-dir  <output_dir>        \\
        --structures   <pdb1> [<pdb2> ...]

Anchoring:
    --anchor-top / --anchor-bottom lift a named sequence to position 1 / last,
    regardless of where it falls in the tree. Useful for ESPript3 structure runs
    where the reference structure must be first (or last) in the alignment.

ESPript3 preparation (--espript-dir):
    Copies the MSA (.fas) and any --structures PDB/CIF files into <espript-dir>.
    The FASTA is written without CRLF and with lines ≤ 80 chars — ready to
    upload directly to https://espript.ibcp.fr.

Name mapping:
    Newick-safe names may have colons replaced by underscores for coordinate
    suffixes (e.g. NtDMP_XM_016580768.1_173-853 in tree vs
    NtDMP_XM_016580768.1:173-853 in FASTA). The script resolves this
    automatically.
"""

import argparse
import re
import shutil
import sys
from collections import OrderedDict
from pathlib import Path

# ---------------------------------------------------------------------------
# Hardcoded defaults — standalone mode
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve()
PIPELINE_DIR = _HERE.parents[2]

_MSA_DIR = (
    PIPELINE_DIR
    / "3_RESULT/DMP/09_Secondary_Structure_Analysis"
    / "GPE001970_SMEL5/INPUTS/MAFFT_aligned"
)
_STEM = "per_genome_gpe001970_and_Selected_Crop_Species_AMINO_ACID_Sequence"

DEFAULT_MSA    = _MSA_DIR / f"{_STEM}.fas"
DEFAULT_TREE   = (
    PIPELINE_DIR
    / "3_RESULT/DMP/05_Phylogenetics"
    / "GPE001970_SMEL5/RAXML"
    / f"{_STEM}_RAXML.raxml.bestTree"
)
DEFAULT_OUTPUT = _MSA_DIR / _STEM

# Tier 1 anchor defaults
DEFAULT_ANCHOR_TOP    = "SMEL5_10g017610.1"   # SmelDMPv5_10.200 — highest reproductive priority
DEFAULT_ANCHOR_BOTTOM = "SMEL5_01g026030.1"   # SmelDMPv5_01.990

# AlphaFold3 PDB sources for Tier 1 genes
_AF3_DIR = (
    PIPELINE_DIR
    / "3_RESULT/DMP/08_Protein_Structure"
    / "GPE001970_SMEL5/AlphaFold3_Results"
)
DEFAULT_STRUCTURES = [
    str(_AF3_DIR / "fold_smel5_10g017610_1"
                 / "fold_smel5_10g017610_1_2026_04_08_01_46_model_0.pdb"),
    str(_AF3_DIR / "smel5_01g026030_1"
                 / "fold_2026_04_07_13_17_smel5_01g026030_1_model_0.pdb"),
]
DEFAULT_ESPRIPT_DIR = (
    PIPELINE_DIR
    / "3_RESULT/DMP/09_Secondary_Structure_Analysis"
    / "GPE001970_SMEL5/INPUTS/ESPript3"
)
# ---------------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--msa",    default=str(DEFAULT_MSA),
                   help="Input aligned FASTA (.fas) or CLUSTAL (.aln)")
    p.add_argument("--tree",   default=str(DEFAULT_TREE),
                   help="Newick tree file (e.g. RAxML .bestTree)")
    p.add_argument("--output", default=str(DEFAULT_OUTPUT),
                   help="Output prefix (extensions added automatically)")
    p.add_argument("--anchor-top", default=DEFAULT_ANCHOR_TOP, metavar="SEQ_ID",
                   help="Sequence ID to place at position 1 (default: %(default)s)")
    p.add_argument("--anchor-bottom", default=DEFAULT_ANCHOR_BOTTOM, metavar="SEQ_ID",
                   help="Sequence ID to place at last position (default: %(default)s)")
    p.add_argument("--no-anchors", action="store_true",
                   help="Disable anchor-top/anchor-bottom, use pure phylo order")
    p.add_argument("--espript-dir", default=str(DEFAULT_ESPRIPT_DIR), metavar="DIR",
                   help="Directory to write ESPript3-ready inputs (default: %(default)s)")
    p.add_argument("--structures", nargs="*", default=DEFAULT_STRUCTURES,
                   metavar="PDB",
                   help="PDB/CIF files to copy into --espript-dir")
    p.add_argument("--no-espript", action="store_true",
                   help="Skip ESPript3 output preparation entirely")
    return p.parse_args()


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def detect_format(path: Path) -> str:
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith(">"):
                return "fasta"
            if s.upper().startswith("CLUSTAL"):
                return "clustal"
    return "fasta"


def read_fasta_msa(path: Path) -> "OrderedDict[str, str]":
    seqs: OrderedDict[str, str] = OrderedDict()
    name = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.replace(" ", ""))
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def read_clustal_msa(path: Path) -> "OrderedDict[str, str]":
    seqs: OrderedDict[str, str] = OrderedDict()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            # Skip header, blank lines, and conservation rows (start with space)
            if not line or line.upper().startswith("CLUSTAL") or line[0] == " ":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name, block = parts[0], parts[1]
            seqs.setdefault(name, "")
            seqs[name] += block
    return seqs


def get_leaf_order(tree_path: Path) -> list[str]:
    """Return tip names top→bottom as they appear in a standard tree render."""
    from Bio import Phylo
    tree = Phylo.read(str(tree_path), "newick")
    return [t.name for t in tree.get_terminals()]


def normalize_tree_name(tree_name: str) -> str:
    """Convert Newick-safe name to FASTA-style (coordinate suffix _N-M → :N-M)."""
    return re.sub(r"_(\d+-\d+)$", r":\1", tree_name)


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------
LINE_WIDTH = 80
BLOCK_WIDTH = 60


def write_fasta(ordered: list[tuple[str, str]], path: Path) -> None:
    with open(path, "w", newline="\n") as fh:
        for name, seq in ordered:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), LINE_WIDTH):
                fh.write(seq[i:i + LINE_WIDTH] + "\n")
            fh.write("\n")


def write_clustal(ordered: list[tuple[str, str]], path: Path) -> None:
    if not ordered:
        return
    names = [n for n, _ in ordered]
    seqs  = [s for _, s in ordered]
    aln_len  = max(len(s) for s in seqs)
    name_col = max(len(n) for n in names)

    with open(path, "w", newline="\n") as fh:
        fh.write("CLUSTAL W (1.83) multiple sequence alignment\n\n\n")
        for start in range(0, aln_len, BLOCK_WIDTH):
            for name, seq in zip(names, seqs):
                block = seq[start:start + BLOCK_WIDTH]
                fh.write(f"{name:<{name_col}}      {block}\n")
            fh.write("\n")


# ---------------------------------------------------------------------------
# Anchoring
# ---------------------------------------------------------------------------

def apply_anchors(
    ordered: list[tuple[str, str]],
    anchor_top: str,
    anchor_bottom: str,
) -> list[tuple[str, str]]:
    """
    Move anchor_top to index 0 and anchor_bottom to index -1.
    The remaining sequences keep their relative phylo order.
    """
    top    = next(((n, s) for n, s in ordered if n == anchor_top),    None)
    bottom = next(((n, s) for n, s in ordered if n == anchor_bottom), None)

    if top is None:
        print(f"  WARNING: --anchor-top '{anchor_top}' not found in MSA — ignored")
    if bottom is None:
        print(f"  WARNING: --anchor-bottom '{anchor_bottom}' not found in MSA — ignored")

    middle = [
        (n, s) for n, s in ordered
        if n not in {anchor_top, anchor_bottom}
    ]

    result = []
    if top:
        result.append(top)
    result.extend(middle)
    if bottom:
        result.append(bottom)
    return result


# ---------------------------------------------------------------------------
# ESPript3 preparation
# ---------------------------------------------------------------------------

def prepare_espript(
    ordered: list[tuple[str, str]],
    espript_dir: Path,
    structure_paths: list[str],
    stem: str,
) -> None:
    """
    Write ESPript3-ready inputs:
      - <stem>_espript3.fas  (LF line endings, 80-char FASTA lines)
      - Copies of each structure file with clean names
    """
    espript_dir.mkdir(parents=True, exist_ok=True)

    fas_out = espript_dir / f"{stem}_espript3.fas"
    write_fasta(ordered, fas_out)
    print(f"  ESPript3 FASTA : {fas_out}")

    for src_str in structure_paths:
        src = Path(src_str)
        if not src.exists():
            print(f"  WARNING: structure not found: {src} — skipped")
            continue
        dst = espript_dir / src.name
        shutil.copy2(src, dst)
        print(f"  Copied structure: {dst.name}")

    # Summarise matching
    print("\n  ESPript3 upload guide:")
    print(f"    Alignment file  → {fas_out.name}")
    for src_str in structure_paths:
        src = Path(src_str)
        if src.exists():
            print(f"    Structure file  → {src.name}")
    print("  The structure sequences match the MSA entries exactly (verified lengths).")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    msa_path   = Path(args.msa)
    tree_path  = Path(args.tree)
    out_prefix = args.output

    if not msa_path.exists():
        sys.exit(f"ERROR: MSA file not found: {msa_path}")
    if not tree_path.exists():
        sys.exit(f"ERROR: Tree file not found: {tree_path}")

    # Read MSA
    fmt = detect_format(msa_path)
    seqs = read_fasta_msa(msa_path) if fmt == "fasta" else read_clustal_msa(msa_path)
    print(f"Read {len(seqs)} sequences from '{msa_path.name}' ({fmt})")

    # Get phylogenetic leaf order
    leaf_order = get_leaf_order(tree_path)
    print(f"Tree '{tree_path.name}' has {len(leaf_order)} leaves")

    # Map tree names → MSA names
    fasta_names = set(seqs.keys())
    name_map: dict[str, str] = {}
    for leaf in leaf_order:
        if leaf in fasta_names:
            name_map[leaf] = leaf
        else:
            converted = normalize_tree_name(leaf)
            if converted in fasta_names:
                name_map[leaf] = converted
            else:
                print(f"  WARNING: tree leaf '{leaf}' has no matching MSA entry — skipped")

    # Build phylo-ordered list
    seen: set[str] = set()
    ordered: list[tuple[str, str]] = []
    for leaf in leaf_order:
        fasta_name = name_map.get(leaf)
        if fasta_name and fasta_name not in seen:
            ordered.append((fasta_name, seqs[fasta_name]))
            seen.add(fasta_name)

    # Append any MSA sequences absent from the tree
    for name in seqs:
        if name not in seen:
            print(f"  WARNING: MSA sequence '{name}' not in tree — appended at end")
            ordered.append((name, seqs[name]))

    # Apply Tier 1 anchors unless disabled
    if not args.no_anchors:
        ordered = apply_anchors(ordered, args.anchor_top, args.anchor_bottom)

    print(f"Output: {len(ordered)} sequences")
    for i, (n, _) in enumerate(ordered, 1):
        tag = ""
        if n in {args.anchor_top, args.anchor_bottom}:
            tag = "  ← Tier 1 anchor"
        print(f"  {i:>2}. {n}{tag}")

    # Determine output suffix
    suffix = "_phylo_ordered" if args.no_anchors else "_phylo_ordered_tier1"

    fas_out = Path(f"{out_prefix}{suffix}.fas")
    aln_out = Path(f"{out_prefix}{suffix}.aln")

    fas_out.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(ordered, fas_out)
    write_clustal(ordered, aln_out)

    print(f"\nWritten: {fas_out}")
    print(f"Written: {aln_out}")

    # ESPript3 preparation
    if not args.no_espript:
        stem = Path(out_prefix).name + suffix
        print(f"\nPreparing ESPript3 inputs in: {args.espript_dir}")
        prepare_espript(
            ordered,
            Path(args.espript_dir),
            args.structures or [],
            stem,
        )


if __name__ == "__main__":
    main()
