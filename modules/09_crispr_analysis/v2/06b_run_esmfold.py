#!/usr/bin/env python3
"""
Module: 06b_run_esmfold.py
Stage [6b] — Local ESMFold structure prediction for guides flagged in [6].

Reads every <stem>.protein.tsv produced by stage [6] in --protein-tsv-dir.
For each row whose `structure_flag == "recommend_structure"` and whose
indel{i}_seq points to an existing CDS FASTA, translates the CDS to protein
and queues it for folding under the id

    <stem>_<guideId>_indel{i}

All queued jobs are written to a single JSONL and folded in one
`conda run -n <esmfold-env> python _esmfold_worker.py ...` call so the
~8 GB ESMFold model is loaded only once. PDBs land in --outdir.

Each <stem>.protein.tsv is then rewritten with three new columns per indel:
    indel{i}_pdb_path       relative path to PDB (or "")
    indel{i}_pdb_status     ok | skipped | error | not_flagged
    indel{i}_pdb_msg        worker message (e.g. peak VRAM / skip reason)

This module runs under the main `crispr_v2` env. The folding subprocess
runs under a dedicated `esmfold` env (PyTorch 2.8 + cu128 + openfold).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _log import _log  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Local ESMFold for stage [6] flagged guides")
    p.add_argument("--protein-tsv-dir", required=True,
                   help="Stage [6] output dir containing *.protein.tsv")
    p.add_argument("--transcripts-dir", required=True,
                   help="Stage [5] output dir; root for indel{i}_seq relative paths")
    p.add_argument("--outdir",          required=True,
                   help="Where ESMFold PDBs are written (typically <stage6>/esmfold_structures)")
    p.add_argument("--esmfold-env",     default="esmfold",
                   help="Conda env that has fair-esm + openfold installed")
    p.add_argument("--worker",          default="",
                   help="Path to _esmfold_worker.py (defaults to sibling of this script)")
    p.add_argument("--chunk-size",      type=int, default=64)
    p.add_argument("--max-protein-length", type=int, default=400,
                   help="Skip proteins longer than this (RTX 5050 / 8 GB VRAM cap)")
    p.add_argument("--overwrite",       action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()


def translate_cds(fasta_path: Path) -> str:
    """Translate first record of a CDS FASTA to protein (stops at first stop)."""
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
        rec = next(SeqIO.parse(str(fasta_path), "fasta"), None)
        if rec is None:
            return ""
        nt = str(rec.seq).upper()
        nt = nt[: len(nt) - (len(nt) % 3)]
        return str(Seq(nt).translate(to_stop=True))
    except Exception as exc:
        _log(f"[06b_esmfold] translate failed for {fasta_path}: {exc}", level="WARN")
        return ""


def collect_jobs(protein_tsv_dir: Path,
                 transcripts_dir: Path,
                 outdir: Path,
                 max_len: int,
                 overwrite: bool):
    """Return (jobs_list, per_tsv_plan) where:

    jobs_list = [{"id": ..., "sequence": ...}, ...]      -> JSONL for the worker
    per_tsv_plan = {tsv_path: [(row_index, indel_index, job_id, status, msg)]}
                                                        -> for TSV rewrite
    """
    jobs = []
    plan: dict[Path, list[tuple]] = {}
    seen_ids: set[str] = set()

    for tsv in sorted(protein_tsv_dir.glob("*.protein.tsv")):
        stem = tsv.name.split(".")[0]
        with open(tsv, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
            fields = list(reader.fieldnames or [])
        n_indels = sum(1 for f in fields if f.startswith("indel") and f.endswith("_seq"))
        per_rows: list[tuple] = []
        for r_idx, row in enumerate(rows):
            flag = (row.get("structure_flag") or "").strip()
            guide_id = (row.get("guideId") or row.get("guide_id") or f"row{r_idx}").strip()
            guide_id = guide_id.replace("/", "_").replace(" ", "_")
            for i in range(1, n_indels + 1):
                fa_rel = (row.get(f"indel{i}_seq") or "NA").strip()
                if flag != "recommend_structure":
                    per_rows.append((r_idx, i, "", "not_flagged", ""))
                    continue
                if not fa_rel or fa_rel == "NA":
                    per_rows.append((r_idx, i, "", "skipped", "no_cds_fasta"))
                    continue
                fa_path = transcripts_dir / fa_rel
                if not fa_path.exists():
                    per_rows.append((r_idx, i, "", "error", f"missing:{fa_rel}"))
                    continue
                protein = translate_cds(fa_path)
                if not protein:
                    per_rows.append((r_idx, i, "", "error", "translate_failed"))
                    continue
                if len(protein) > max_len:
                    per_rows.append((r_idx, i, "", "skipped",
                                     f"too_long:{len(protein)}>{max_len}"))
                    continue
                jid = f"{stem}_{guide_id}_indel{i}"
                # Disambiguate if two guideIds collapse to the same id
                if jid in seen_ids:
                    jid = f"{jid}_r{r_idx}"
                seen_ids.add(jid)
                pdb_path = outdir / f"{jid}.pdb"
                if not overwrite and pdb_path.exists():
                    per_rows.append((r_idx, i, jid, "skipped", "exists"))
                    continue
                jobs.append({"id": jid, "sequence": protein})
                per_rows.append((r_idx, i, jid, "queued", f"len={len(protein)}"))
        plan[tsv] = per_rows
    return jobs, plan


def run_worker(esmfold_env: str,
               worker_script: Path,
               jobs_jsonl: Path,
               outdir: Path,
               chunk_size: int,
               max_len: int,
               overwrite: bool) -> dict[str, dict]:
    """Invoke _esmfold_worker.py via `conda run -n <env>` and collect results."""
    cmd = [
        "conda", "run", "-n", esmfold_env, "--no-capture-output",
        "python3", str(worker_script),
        "--jobs",       str(jobs_jsonl),
        "--outdir",     str(outdir),
        "--chunk-size", str(chunk_size),
        "--max-length", str(max_len),
        "--overwrite" if overwrite else "--no-overwrite",
    ]
    _log(f"[06b_esmfold] {' '.join(shlex.quote(c) for c in cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.stderr:
        for line in proc.stderr.splitlines():
            print(line, file=sys.stderr, flush=True)
    results: dict[str, dict] = {}
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line or not line.startswith("{"):
            continue
        try:
            rec = json.loads(line)
            results[rec["id"]] = rec
        except (json.JSONDecodeError, KeyError):
            continue
    if proc.returncode not in (0, 3):  # 3 = some jobs errored, others may be ok
        _log(f"[06b_esmfold] worker exited with code {proc.returncode}", level="ERROR")
    return results


def rewrite_tsv(tsv: Path,
                per_rows: list[tuple],
                results: dict[str, dict],
                outdir: Path) -> tuple[int, int, int]:
    """Add indel{i}_pdb_path/_status/_msg columns. Return (ok, skipped, error)."""
    with open(tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    n_indels = sum(1 for f in fields if f.startswith("indel") and f.endswith("_seq"))
    new_cols: list[str] = []
    for i in range(1, n_indels + 1):
        for suffix in ("pdb_path", "pdb_status", "pdb_msg"):
            col = f"indel{i}_{suffix}"
            if col not in fields:
                new_cols.append(col)
    new_fields = fields + new_cols

    # Initialize new cells
    for row in rows:
        for col in new_cols:
            row.setdefault(col, "")

    n_ok = n_skip = n_err = 0
    for r_idx, indel_i, jid, status, msg in per_rows:
        if not (0 <= r_idx < len(rows)):
            continue
        row = rows[r_idx]
        if jid and jid in results:
            rec = results[jid]
            status = rec.get("status", status)
            msg    = rec.get("msg", msg)
            pdb    = rec.get("pdb", "")
        elif status == "queued":
            # Queued but worker returned no record → treat as error
            status, msg, pdb = "error", "no_worker_response", ""
        else:
            pdb = ""
        if pdb:
            try:
                pdb_rel = str(Path(pdb).relative_to(outdir.parent))
            except ValueError:
                pdb_rel = pdb
        else:
            pdb_rel = ""
        row[f"indel{indel_i}_pdb_path"]   = pdb_rel
        row[f"indel{indel_i}_pdb_status"] = status
        row[f"indel{indel_i}_pdb_msg"]    = msg
        if status == "ok":
            n_ok += 1
        elif status == "skipped" or status == "not_flagged":
            n_skip += 1
        elif status == "error":
            n_err += 1

    with open(tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return n_ok, n_skip, n_err


def main() -> int:
    args = parse_args()

    protein_dir = Path(args.protein_tsv_dir)
    transcripts_dir = Path(args.transcripts_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not protein_dir.is_dir():
        _log(f"[06b_esmfold] protein-tsv-dir does not exist: {protein_dir}",
             level="ERROR")
        return 1

    worker = Path(args.worker) if args.worker else \
             Path(__file__).parent / "_esmfold_worker.py"
    if not worker.exists():
        _log(f"[06b_esmfold] worker script not found: {worker}", level="ERROR")
        return 1

    jobs, plan = collect_jobs(protein_dir, transcripts_dir, outdir,
                              args.max_protein_length, args.overwrite)
    n_tsvs = len(plan)
    n_jobs = len(jobs)
    _log(f"[06b_esmfold] {n_jobs} sequences queued from {n_tsvs} protein.tsv files")

    results: dict[str, dict] = {}
    if n_jobs > 0:
        with tempfile.NamedTemporaryFile("w", suffix=".jsonl",
                                         delete=False, dir=str(outdir)) as tf:
            for job in jobs:
                tf.write(json.dumps(job) + "\n")
            jobs_path = Path(tf.name)
        try:
            results = run_worker(args.esmfold_env, worker, jobs_path,
                                 outdir, args.chunk_size,
                                 args.max_protein_length, args.overwrite)
        finally:
            try:
                os.remove(jobs_path)
            except OSError:
                pass
    else:
        _log("[06b_esmfold] No flagged guides — skipping ESMFold invocation.")

    tot_ok = tot_skip = tot_err = 0
    for tsv, per_rows in plan.items():
        ok, skip, err = rewrite_tsv(tsv, per_rows, results, outdir)
        tot_ok   += ok
        tot_skip += skip
        tot_err  += err
        _log(f"[06b_esmfold] {tsv.name}: ok={ok} skipped={skip} error={err}")

    _log(f"[06b_esmfold] TOTAL ok={tot_ok} skipped={tot_skip} error={tot_err} "
         f"(PDBs in {outdir})")
    return 0 if tot_err == 0 else 4


if __name__ == "__main__":
    sys.exit(main())
