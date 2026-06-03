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
HAP2 chain *that matches the template* gets that mmCIF embedded with a
per-variant queryIndices / templateIndices alignment.

Two correctness rules are enforced before a template is embedded (the
AlphaFold Server rejects the upload otherwise, which surfaces as the per-job
"error uploading files" message):

  1. SINGLE CHAIN. AlphaFold3 refuses any template mmCIF containing more than
     one chain ("it is not possible to determine which of the chains should be
     used as the template"). SWISS-MODEL homo-oligomer models (chains A/B/C)
     are reduced in-memory to the first chain only: rows in every
     chain-keyed category (_atom_site, _struct_asym, _pdbx_poly_seq_scheme,
     _ma_qa_metric_local, _ma_target_entity_instance) are filtered to the kept
     chain, and _entity.pdbx_number_of_molecules / _entity_poly.pdbx_strand_id
     are rewritten to reflect a monomer.

  2. RESIDUE-IDENTITY GATE. queryIndices are 0-based into the query
     (proteinChain) sequence; templateIndices are 0-based into the template's
     _entity_poly_seq (the full-length polymer, including residues with no
     coordinates). Anchors are only emitted for original positions inside
     [coverage_start, coverage_end] that (a) survive the variant's deletions,
     (b) are resolved (have coordinates) in the template, and (c) exist in the
     template polymer. The fraction of those anchors whose query residue
     matches the template residue is computed per variant; if it falls below
     --template-min-identity the template is dropped for that variant (the JSON
     uploads template-free and AlphaFold3 runs its own PDB template search).
     This makes a mis-paired template (e.g. an eggplant SmelHAP2 model attached
     to an Arabidopsis AtHAP2 query) silently no-op instead of breaking the
     upload, while a correctly-paired same-protein model attaches cleanly.

Job name format: "{pair_label}_{stoich_label}" (e.g. "WT__WT_monomeric",
"delC__WT_postfusion_like") so each downloaded result zip carries both axes
in its filename for the orchestrator's SCATTER step.
"""

import argparse
import json
import sys
from pathlib import Path

# Three-letter -> one-letter amino-acid map. Non-standard / unknown monomers
# map to "X" so they never count as an identity match in the gate.
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}


def three_to_one(mon):
    return _THREE_TO_ONE.get((mon or "").upper(), "X")


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


# ---------------------------------------------------------------------------
# Minimal, dependency-free mmCIF handling (stdlib only -- the generator runs
# in the pure-stdlib prepare_variants env). We never build a full structure;
# we only need to (a) reduce a multi-chain model to a single chain and (b)
# read the resolved-residue set and the entity_poly_seq sequence for the
# identity gate.
# ---------------------------------------------------------------------------

# category -> the column tag that carries the chain id; rows whose chain id is
# not the kept chain are dropped during single-chain reduction.
_CHAIN_KEYED = {
    "_atom_site": "_atom_site.label_asym_id",
    "_atom_site_anisotrop": "_atom_site_anisotrop.label_asym_id",
    "_struct_asym": "_struct_asym.id",
    "_pdbx_poly_seq_scheme": "_pdbx_poly_seq_scheme.asym_id",
    "_pdbx_nonpoly_scheme": "_pdbx_nonpoly_scheme.asym_id",
    "_ma_qa_metric_local": "_ma_qa_metric_local.label_asym_id",
    "_ma_target_entity_instance": "_ma_target_entity_instance.asym_id",
}


def _category_of(tag):
    return tag.split(".", 1)[0]


def split_cif_line(line):
    """Tokenize one mmCIF data line, honouring '...' and \"...\" quoting.
    Adequate for single-line data rows (the only rows we parse/rewrite); it is
    NOT used on multi-line ';'-delimited fields.
    """
    tokens = []
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c in " \t":
            i += 1
            continue
        if c == "#":          # comment to end of line (only at token start)
            break
        if c in ("'", '"'):
            q = c
            i += 1
            buf = []
            while i < n:
                if line[i] == q and (i + 1 >= n or line[i + 1] in " \t"):
                    break
                buf.append(line[i])
                i += 1
            tokens.append("".join(buf))
            i += 1            # skip closing quote
        else:
            start = i
            while i < n and line[i] not in " \t":
                i += 1
            tokens.append(line[start:i])
    return tokens


def _format_value(v):
    """Re-quote a token for emission if it contains whitespace or starts with a
    quote/semicolon. Bare tokens (the common case) pass through unchanged.
    """
    if v == "":
        return "."
    if any(ch in v for ch in " \t") or v[0] in "'\";":
        return '"' + v + '"' if '"' not in v else "'" + v + "'"
    return v


def iter_loop_rows(text, wanted):
    """Yield (category, headers, tokens) for single-line data rows of the
    `wanted` categories. ';'-delimited multi-line fields are tracked and
    skipped so they are never mistaken for data rows.
    """
    in_multiline = False
    reading_headers = False
    cur_cat = None
    cur_headers = []
    for ln in text.splitlines():
        if in_multiline:
            if ln.startswith(";"):
                in_multiline = False
            continue
        if ln.startswith(";"):
            in_multiline = True
            continue
        s = ln.strip()
        if s == "loop_":
            reading_headers = True
            cur_headers = []
            cur_cat = None
            continue
        if reading_headers:
            if s.startswith("_"):
                cur_headers.append(s)
                continue
            reading_headers = False
            cur_cat = _category_of(cur_headers[0]) if cur_headers else None
            # fall through: this line may already be the first data row
        if s == "" or s.startswith("#"):
            cur_cat = None
            continue
        if s.startswith("_") or s.startswith("data_"):
            cur_cat = None
            continue
        if cur_cat in wanted:
            yield cur_cat, cur_headers, split_cif_line(ln)


def read_mmcif_metadata(text):
    """Return chain list, kept chain, resolved-residue set, and the
    num -> one-letter map for the kept chain's entity. None if no chain found.
    """
    struct_asym = []                  # [(chain_id, entity_id)]
    atom_chains = []                  # ordered unique label_asym_id
    seen_atom_chains = set()
    resolved_by_chain = {}            # chain -> set(int label_seq_id)
    polyseq = []                      # [(entity_id, num, mon_id)]

    wanted = {"_struct_asym", "_atom_site", "_entity_poly_seq"}
    for cat, headers, toks in iter_loop_rows(text, wanted):
        if len(toks) != len(headers):
            continue
        row = dict(zip(headers, toks))
        if cat == "_struct_asym":
            cid = row.get("_struct_asym.id")
            if cid is not None:
                struct_asym.append((cid, row.get("_struct_asym.entity_id")))
        elif cat == "_atom_site":
            ch = row.get("_atom_site.label_asym_id")
            if ch is None:
                continue
            if ch not in seen_atom_chains:
                seen_atom_chains.add(ch)
                atom_chains.append(ch)
            seq = row.get("_atom_site.label_seq_id")
            if seq not in (None, ".", "?"):
                try:
                    resolved_by_chain.setdefault(ch, set()).add(int(seq))
                except ValueError:
                    pass
        elif cat == "_entity_poly_seq":
            polyseq.append((
                row.get("_entity_poly_seq.entity_id"),
                row.get("_entity_poly_seq.num"),
                row.get("_entity_poly_seq.mon_id"),
            ))

    chains = [c for c, _ in struct_asym] or atom_chains
    if not chains:
        return None
    kept = chains[0]
    kept_entity = next((eid for cid, eid in struct_asym if cid == kept), None)

    resolved = resolved_by_chain.get(kept, set())
    num_to_one = {}
    for eid, num, mon in polyseq:
        # If we could not resolve the kept chain's entity, fall back to taking
        # every entity_poly_seq row (single-entity templates are the norm).
        if kept_entity is not None and eid != kept_entity:
            continue
        try:
            num_to_one[int(num)] = three_to_one(mon)
        except (TypeError, ValueError):
            continue

    return {
        "chains": chains,
        "kept_chain": kept,
        "resolved": resolved,
        "num_to_one": num_to_one,
    }


def reduce_to_single_chain(text, kept_chain, all_chains):
    """Return `text` with every chain other than `kept_chain` removed so the
    AlphaFold3 single-chain template requirement is met. Non-chain-keyed
    categories pass through verbatim; molecule-count metadata is rewritten to
    reflect a monomer.
    """
    chains_joined = ",".join(all_chains)
    out = []
    in_multiline = False
    reading_headers = False
    cur_cat = None
    cur_headers = []
    chain_idx = None      # column index of the chain id in a chain-keyed loop
    mol_idx = None        # column index of pdbx_number_of_molecules in _entity

    for ln in text.splitlines():
        if in_multiline:
            out.append(ln)
            if ln.startswith(";"):
                in_multiline = False
            continue
        if ln.startswith(";"):
            in_multiline = True
            out.append(ln)
            continue
        s = ln.strip()
        if s == "loop_":
            reading_headers = True
            cur_headers = []
            cur_cat = None
            chain_idx = None
            mol_idx = None
            out.append(ln)
            continue
        if reading_headers:
            if s.startswith("_"):
                cur_headers.append(s)
                out.append(ln)
                continue
            reading_headers = False
            cur_cat = _category_of(cur_headers[0]) if cur_headers else None
            if cur_cat in _CHAIN_KEYED:
                try:
                    chain_idx = cur_headers.index(_CHAIN_KEYED[cur_cat])
                except ValueError:
                    chain_idx = None
            elif cur_cat == "_entity":
                try:
                    mol_idx = cur_headers.index("_entity.pdbx_number_of_molecules")
                except ValueError:
                    mol_idx = None
            # fall through to handle this line as data / terminator
        if s == "" or s.startswith("#"):
            cur_cat = None
            chain_idx = None
            mol_idx = None
            out.append(ln)
            continue
        if s.startswith("_") or s.startswith("data_"):
            cur_cat = None
            chain_idx = None
            mol_idx = None
            out.append(ln)
            continue

        # --- data row of cur_cat ---
        if cur_cat in _CHAIN_KEYED and chain_idx is not None:
            toks = split_cif_line(ln)
            if len(toks) == len(cur_headers) and toks[chain_idx] != kept_chain:
                continue                          # drop non-kept chain
            out.append(ln)
            continue
        if cur_cat == "_entity" and mol_idx is not None:
            toks = split_cif_line(ln)
            if len(toks) == len(cur_headers):
                toks[mol_idx] = "1"
                out.append(" ".join(_format_value(t) for t in toks))
                continue
            out.append(ln)
            continue
        if cur_cat == "_entity_poly" and chains_joined and chains_joined in ln:
            toks = split_cif_line(ln)
            if chains_joined in toks:
                toks = [kept_chain if t == chains_joined else t for t in toks]
                out.append(" ".join(_format_value(t) for t in toks))
                continue
            out.append(ln)
            continue
        out.append(ln)

    return "\n".join(out) + ("\n" if text.endswith("\n") else "")


def build_template_alignment(deletions, coverage_start, coverage_end,
                             query_seq, resolved, num_to_one):
    """Build parallel (queryIndices, templateIndices) arrays plus identity
    stats for one HAP2 variant.

    queryIndices: 0-based into `query_seq` (the post-deletion variant).
    templateIndices: 0-based into the template _entity_poly_seq (orig_pos - 1).

    Only original positions in [coverage_start, coverage_end] that survive the
    variant's deletions, are resolved in the template, and exist in the
    template polymer become anchors. Returns (q, t, n_match, n_considered).
    """
    q, t = [], []
    n_match = 0
    n_considered = 0
    qlen = len(query_seq)
    for orig_pos in range(coverage_start, coverage_end + 1):
        new_pos = position_after_deletions(orig_pos, deletions)
        if new_pos is None:                 # deleted in this variant
            continue
        if orig_pos not in resolved:        # no coordinates in template
            continue
        tmpl_res = num_to_one.get(orig_pos)
        if tmpl_res is None:                # outside template polymer
            continue
        qi = new_pos - 1
        if qi < 0 or qi >= qlen:
            continue
        n_considered += 1
        if query_seq[qi] == tmpl_res:
            n_match += 1
        q.append(qi)
        t.append(orig_pos - 1)
    return q, t, n_match, n_considered


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
    # Optional HAP2 template. When --hap2-template-mmcif is supplied, the mmCIF
    # is reduced to a single chain and embedded into every HAP2 chain whose
    # sequence matches the template (see --template-min-identity).
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
    p.add_argument("--template-min-identity", type=float, default=0.7,
                   help="minimum fraction of anchored positions whose query residue "
                        "must match the template residue for the template to be "
                        "embedded for a given variant (default 0.7). Below this the "
                        "template is dropped and AlphaFold3 does its own search.")
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

    # Template loading + single-chain reduction. Done once; the reduced mmCIF
    # text is reused across every HAP2 chain that passes the identity gate.
    single_chain_mmcif = ""
    tmpl_resolved = set()
    tmpl_num_to_one = {}
    tmpl_kept_chain = ""
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
        raw_mmcif = tpath.read_text(encoding="utf-8")
        meta = read_mmcif_metadata(raw_mmcif)
        if meta is None:
            print(f"  [warn] no chain found in template {tpath}; skipping template embedding.",
                  file=sys.stderr)
        else:
            tmpl_kept_chain = meta["kept_chain"]
            tmpl_resolved = meta["resolved"]
            tmpl_num_to_one = meta["num_to_one"]
            if len(meta["chains"]) > 1:
                single_chain_mmcif = reduce_to_single_chain(
                    raw_mmcif, tmpl_kept_chain, meta["chains"]
                )
                print(
                    f"  [template] reduced multi-chain model "
                    f"({','.join(meta['chains'])}) to single chain "
                    f"'{tmpl_kept_chain}' for AF3 upload.",
                    file=sys.stderr,
                )
            else:
                single_chain_mmcif = raw_mmcif

    # Pipe-separated to keep bash quoting sane (deletion strings already use
    # commas internally for multi-range entries).
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

    # Cache FASTA reads + per-(sequence, deletion) template alignments so the
    # repeated WT references in orthogonal mode don't redo the work.
    seq_cache = {}

    def seq_of(path):
        if path not in seq_cache:
            _, seq_cache[path] = read_first_fasta(path)
        return seq_cache[path]

    align_cache = {}
    template_stats = []   # (variant_key, n_match, n_considered, attached)

    def template_for(hap2_path, deletion_str):
        """Return a template dict for this variant, or None if no template is
        configured or the variant fails the residue-identity gate.
        """
        if not single_chain_mmcif:
            return None
        key = (hap2_path, deletion_str or "")
        if key in align_cache:
            return align_cache[key]
        dels = parse_deletions(deletion_str)
        query_seq = seq_of(hap2_path)
        q, t, n_match, n_considered = build_template_alignment(
            dels, args.hap2_coverage_start, args.hap2_coverage_end,
            query_seq, tmpl_resolved, tmpl_num_to_one,
        )
        frac = (n_match / n_considered) if n_considered else 0.0
        attached = n_considered > 0 and frac >= args.template_min_identity
        template_stats.append((key, n_match, n_considered, attached))
        result = None
        if attached:
            result = {
                "mmcif": single_chain_mmcif,
                "queryIndices": q,
                "templateIndices": t,
            }
        align_cache[key] = result
        return result

    jobs = []
    for pair, hap2_fa, dmp_fa, hap2_del in zip(
        pair_labels, hap2_paths, dmp_paths, hap2_deletions
    ):
        hap2_seq = seq_of(hap2_fa)
        dmp_seq = seq_of(dmp_fa)
        tmpl = template_for(hap2_fa, hap2_del)
        for slab, n_hap2 in zip(stoich_labels, stoich_counts):
            hap2_chain = build_protein_chain(hap2_seq, n_hap2, tmpl)
            dmp_chain = build_protein_chain(dmp_seq, args.dmp_copies, None)
            jobs.append(build_job(f"{pair}_{slab}", hap2_chain, dmp_chain))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as fh:
        json.dump(jobs, fh, indent=2)
        fh.write("\n")

    if single_chain_mmcif:
        n_attached = sum(1 for _, _, _, a in template_stats if a)
        n_dropped = len(template_stats) - n_attached
        print(
            f"  -> wrote {len(jobs)} AF3 job(s) "
            f"({len(pair_labels)} pair(s) x {len(stoich_labels)} stoich(s)); "
            f"HAP2 template (chain '{tmpl_kept_chain}') attached to "
            f"{n_attached}/{len(template_stats)} unique variant(s) -> {out}",
            file=sys.stderr,
        )
        for (hap2_path, deletion_str), n_match, n_considered, attached in template_stats:
            if attached:
                continue
            frac = (n_match / n_considered) if n_considered else 0.0
            print(
                f"     [template dropped] {Path(hap2_path).name} "
                f"(del='{deletion_str or 'WT'}'): identity "
                f"{n_match}/{n_considered} = {frac:.0%} "
                f"< {args.template_min_identity:.0%}; AlphaFold3 will run its own "
                f"template search.",
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
