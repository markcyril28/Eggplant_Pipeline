#!/usr/bin/env python3
"""
build_af3_job_json.py -- emit a SINGLE AlphaFold Server-compatible JSON job
file covering every (HAP2 variant, DMP variant, stoichiometry) cell of a
Stage-14 experiment.

The AlphaFold Server (https://alphafoldserver.com/) accepts a single JSON
file as job input. The top-level value is an ARRAY of job objects; each job
has a `name`, a `sequences` list (chain-count-encoded), and the literal
fields `"dialect": "alphafoldserver"`, `"version": 1`, `"modelSeeds": []`.

Custom templates are optional. When --hap2-template-mmcif is provided, every
HAP2 chain gets that mmCIF embedded with per-variant queryIndices /
templateIndices arrays computed from the variant's deletion string. The
template typically covers a fixed residue range of the original (full-length)
HAP2 sequence (e.g. SWISS-MODEL of AtHAP2 22-530 from PDB 5OW3); only
template positions whose original residue survives the deletion(s) and falls
inside [coverage_start, coverage_end] are emitted.

Job name format: "{pair_label}_{stoich_label}" (e.g. "WT__WT_monomeric",
"delC__WT_postfusion_like") so each downloaded result zip carries both axes
in its filename for the orchestrator's SCATTER step.
"""

import argparse
import json
import sys
from pathlib import Path


def read_first_fasta(path):
    """Return (header_without_gt, sequence) for the first FASTA record."""
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:]
            else:
                parts.append(line)
    if header is None:
        raise ValueError(f"No FASTA record found in {path}")
    seq = "".join(parts).upper()
    if not seq:
        raise ValueError(f"FASTA record in {path} has no sequence body")
    return header, seq


def parse_deletions(s):
    """Parse '25-530,596-705' (1-indexed inclusive) into [(25, 530), (596, 705)].
    Empty / 'WT' / None returns [].
    """
    if not s or s.strip() in ("", "WT"):
        return []
    ranges = []
    for chunk in s.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        a, b = chunk.split("-")
        ranges.append((int(a), int(b)))
    return ranges


def position_after_deletions(orig_pos, deletions):
    """Map 1-indexed original position to 1-indexed post-deletion position.
    Returns None if orig_pos falls inside a deletion range.
    """
    for s, e in deletions:
        if s <= orig_pos <= e:
            return None
    deleted_before = sum(e - s + 1 for s, e in deletions if e < orig_pos)
    return orig_pos - deleted_before


def build_template_alignment(deletions, coverage_start, coverage_end):
    """Build parallel (queryIndices, templateIndices) arrays (both 0-indexed)
    for a HAP2 variant. Template residue at coverage_start has templateIndex 0.
    Query residues are renumbered to reflect the post-deletion variant sequence.
    Positions that fall inside a deletion are skipped (no entry emitted).
    """
    q, t = [], []
    for orig_pos in range(coverage_start, coverage_end + 1):
        new_pos = position_after_deletions(orig_pos, deletions)
        if new_pos is None:
            continue
        q.append(new_pos - 1)
        t.append(orig_pos - coverage_start)
    return q, t


def build_protein_chain(sequence, count, template=None):
    """Return a proteinChain dict; attach the template only if one was built
    and has at least one anchor residue (an empty alignment makes the template
    a no-op and just bloats the JSON).
    """
    chain = {"sequence": sequence, "count": int(count)}
    if template is not None and template["queryIndices"]:
        chain["templates"] = [template]
    return chain


def build_job(name, hap2_chain, dmp_chain):
    return {
        "name": name,
        "modelSeeds": [],
        "sequences": [
            {"proteinChain": hap2_chain},
            {"proteinChain": dmp_chain},
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pair-labels", required=True,
                   help="space-separated pair labels (e.g. 'WT__WT delC__WT WT__delN')")
    p.add_argument("--hap2-fastas", required=True,
                   help="space-separated HAP2 variant FASTA paths, parallel to --pair-labels")
    p.add_argument("--dmp-fastas", required=True,
                   help="space-separated DMP variant FASTA paths, parallel to --pair-labels")
    p.add_argument("--stoich-labels", required=True,
                   help="space-separated stoich names, e.g. 'monomeric dimeric postfusion_like'")
    p.add_argument("--stoich-counts", required=True,
                   help="space-separated HAP2 copy counts parallel to --stoich-labels, e.g. '1 2 3'")
    p.add_argument("--dmp-copies", type=int, default=1,
                   help="DMP copies per job (default 1)")
    # Optional HAP2 template. When --hap2-template-mmcif is supplied, every
    # HAP2 chain in the output JSON gets the same mmCIF embedded plus a
    # per-variant queryIndices/templateIndices alignment derived from the
    # variant's deletion string.
    p.add_argument("--hap2-template-mmcif", default="",
                   help="path to mmCIF template for HAP2 (e.g. SWISS-MODEL output); "
                        "leave empty to skip template embedding")
    p.add_argument("--hap2-coverage-start", type=int, default=0,
                   help="1-indexed first residue of the original (full-length) HAP2 "
                        "sequence covered by the template")
    p.add_argument("--hap2-coverage-end", type=int, default=0,
                   help="1-indexed last residue covered by the template")
    p.add_argument("--hap2-deletions", default="",
                   help="pipe-separated deletion strings parallel to --pair-labels "
                        "(e.g. '|25-530|' for WT/delEcto/WT). Each entry uses commas "
                        "for multi-range, e.g. '25-530,596-705'. Empty entry = WT.")
    p.add_argument("--output", required=True,
                   help="output JSON path; parent dirs created if missing")
    args = p.parse_args()

    pair_labels = args.pair_labels.split()
    hap2_paths = args.hap2_fastas.split()
    dmp_paths = args.dmp_fastas.split()
    stoich_labels = args.stoich_labels.split()
    stoich_counts = args.stoich_counts.split()

    if not pair_labels:
        sys.exit("ERROR: --pair-labels is empty; nothing to emit")
    if len(pair_labels) != len(hap2_paths) or len(pair_labels) != len(dmp_paths):
        sys.exit(
            f"ERROR: parallel arrays must be the same length: "
            f"pair_labels={len(pair_labels)}, hap2_fastas={len(hap2_paths)}, "
            f"dmp_fastas={len(dmp_paths)}"
        )
    if len(stoich_labels) != len(stoich_counts):
        sys.exit(
            f"ERROR: --stoich-labels ({len(stoich_labels)}) and --stoich-counts "
            f"({len(stoich_counts)}) must be parallel arrays"
        )
    if not stoich_labels:
        sys.exit("ERROR: --stoich-labels is empty; nothing to emit")

    # Template loading. The mmCIF text is read once and reused across every
    # HAP2 chain to keep the helper fast; the JSON serializer still
    # duplicates the string per job in the output file (the AF3 Server JSON
    # has no shared-template construct).
    template_mmcif = ""
    if args.hap2_template_mmcif:
        tpath = Path(args.hap2_template_mmcif)
        if not tpath.is_file():
            sys.exit(f"ERROR: --hap2-template-mmcif not found: {tpath}")
        if not (args.hap2_coverage_start and args.hap2_coverage_end):
            sys.exit(
                "ERROR: when --hap2-template-mmcif is set, both "
                "--hap2-coverage-start and --hap2-coverage-end are required "
                "(1-indexed residue range of the original HAP2 covered by the template)"
            )
        if args.hap2_coverage_start > args.hap2_coverage_end:
            sys.exit(
                f"ERROR: --hap2-coverage-start ({args.hap2_coverage_start}) > "
                f"--hap2-coverage-end ({args.hap2_coverage_end})"
            )
        template_mmcif = tpath.read_text(encoding="utf-8")

    # Pipe-separated to keep bash quoting sane (deletion strings already use
    # commas internally for multi-range entries).
    hap2_deletions = []
    if args.hap2_deletions:
        hap2_deletions = args.hap2_deletions.split("|")
        if len(hap2_deletions) != len(pair_labels):
            sys.exit(
                f"ERROR: --hap2-deletions ({len(hap2_deletions)}) must be parallel "
                f"to --pair-labels ({len(pair_labels)}); pipe-separate one entry per "
                "pair (empty entry for WT)"
            )
    else:
        hap2_deletions = ["" for _ in pair_labels]

    # Cache FASTA reads + per-deletion-string template alignments so repeated
    # WT references (orthogonal mode pins one side to WT) don't redo the work.
    seq_cache = {}
    def seq_of(path):
        if path not in seq_cache:
            _, seq_cache[path] = read_first_fasta(path)
        return seq_cache[path]

    align_cache = {}
    def template_for(deletion_str):
        if not template_mmcif:
            return None
        key = deletion_str or ""
        if key not in align_cache:
            dels = parse_deletions(deletion_str)
            q, t = build_template_alignment(
                dels, args.hap2_coverage_start, args.hap2_coverage_end
            )
            align_cache[key] = {
                "mmcif": template_mmcif,
                "queryIndices": q,
                "templateIndices": t,
            }
        return align_cache[key]

    jobs = []
    for pair, hap2_fa, dmp_fa, hap2_del in zip(
        pair_labels, hap2_paths, dmp_paths, hap2_deletions
    ):
        hap2_seq = seq_of(hap2_fa)
        dmp_seq = seq_of(dmp_fa)
        tmpl = template_for(hap2_del)
        for slab, n_hap2 in zip(stoich_labels, stoich_counts):
            hap2_chain = build_protein_chain(hap2_seq, n_hap2, tmpl)
            dmp_chain = build_protein_chain(dmp_seq, args.dmp_copies, None)
            jobs.append(build_job(f"{pair}_{slab}", hap2_chain, dmp_chain))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as fh:
        json.dump(jobs, fh, indent=2)
        fh.write("\n")

    if template_mmcif:
        n_anchored = sum(1 for k, v in align_cache.items() if v["queryIndices"])
        n_total = len(align_cache)
        print(
            f"  -> wrote {len(jobs)} AF3 job(s) "
            f"({len(pair_labels)} pair(s) x {len(stoich_labels)} stoich(s)); "
            f"HAP2 template attached to {n_anchored}/{n_total} unique variant(s) "
            f"-> {out}",
            file=sys.stderr,
        )
    else:
        print(
            f"  -> wrote {len(jobs)} AF3 job(s) "
            f"({len(pair_labels)} pair(s) x {len(stoich_labels)} stoich(s)) to {out}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
