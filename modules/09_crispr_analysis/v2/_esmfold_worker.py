#!/usr/bin/env python3
"""
Worker invoked by 06b_run_esmfold.py via `conda run -n esmfold python ...`.
Loads ESMFold v1 once, then folds every (id, sequence) record in --jobs (JSONL)
to <outdir>/<id>.pdb. Stays GPU-resident across all jobs to amortize the
multi-GB model load.

Stdout is reserved for one JSON line per job:
    {"id": "...", "pdb": "<path>", "status": "ok|skipped|error", "msg": "..."}
Logs and progress go to stderr.

Required env (set by the conda activation hook in the esmfold env):
    CUDA_HOME, LD_LIBRARY_PATH (with /usr/lib/wsl/lib for WSL libcuda).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path


def log(msg: str, level: str = "INFO") -> None:
    print(f"[esmfold_worker] [{level}] {msg}", file=sys.stderr, flush=True)


def parse_args():
    p = argparse.ArgumentParser(description="ESMFold worker (single-GPU)")
    p.add_argument("--jobs",       required=True, help="JSONL of {id, sequence} records")
    p.add_argument("--outdir",     required=True, help="Directory to write <id>.pdb files")
    p.add_argument("--chunk-size", type=int, default=64,
                   help="model.set_chunk_size(N); lower = less VRAM, slower")
    p.add_argument("--max-length", type=int, default=400,
                   help="Skip sequences longer than this many residues (VRAM cap)")
    p.add_argument("--overwrite",  action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()


def emit(record: dict) -> None:
    print(json.dumps(record), flush=True)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    jobs_path = Path(args.jobs)
    if not jobs_path.exists():
        log(f"jobs file not found: {jobs_path}", level="ERROR")
        return 1

    with open(jobs_path) as fh:
        jobs = [json.loads(line) for line in fh if line.strip()]

    if not jobs:
        log("No folding jobs supplied — nothing to do.")
        return 0

    log(f"Importing torch / esm ...")
    import torch
    import warnings
    warnings.filterwarnings("ignore")
    import esm

    if not torch.cuda.is_available():
        log("CUDA not available — refusing to fold on CPU (would take hours).",
            level="ERROR")
        for job in jobs:
            emit({"id": job["id"], "pdb": "",
                  "status": "error", "msg": "no_cuda"})
        return 2

    log(f"GPU: {torch.cuda.get_device_name(0)} "
        f"({torch.cuda.get_device_properties(0).total_memory/1e9:.1f} GB)")

    log(f"Loading ESMFold v1 (~8 GB weights cached) ...")
    t0 = time.time()
    model = esm.pretrained.esmfold_v1().eval().cuda()
    model.set_chunk_size(args.chunk_size)
    log(f"Model ready in {time.time()-t0:.1f}s; "
        f"VRAM after load: {torch.cuda.memory_allocated()/1e9:.2f} GB")

    n_ok = n_skip = n_err = 0
    for i, job in enumerate(jobs, 1):
        jid = job["id"]
        seq = job["sequence"].strip().upper().replace("*", "")
        out_pdb = outdir / f"{jid}.pdb"

        if not args.overwrite and out_pdb.exists():
            n_skip += 1
            emit({"id": jid, "pdb": str(out_pdb),
                  "status": "skipped", "msg": "exists"})
            continue

        if not seq:
            n_err += 1
            emit({"id": jid, "pdb": "",
                  "status": "error", "msg": "empty_sequence"})
            continue

        if len(seq) > args.max_length:
            n_skip += 1
            emit({"id": jid, "pdb": "",
                  "status": "skipped",
                  "msg": f"too_long ({len(seq)} > {args.max_length})"})
            continue

        log(f"[{i}/{len(jobs)}] Folding {jid} ({len(seq)} aa) ...")
        try:
            t0 = time.time()
            with torch.no_grad():
                pdb = model.infer_pdb(seq)
            out_pdb.write_text(pdb)
            n_ok += 1
            peak = torch.cuda.max_memory_allocated() / 1e9
            log(f"    -> {out_pdb.name}  ({time.time()-t0:.1f}s, peak {peak:.2f} GB)")
            emit({"id": jid, "pdb": str(out_pdb),
                  "status": "ok",
                  "msg": f"len={len(seq)} time_s={time.time()-t0:.1f} peak_gb={peak:.2f}"})
        except torch.cuda.OutOfMemoryError as exc:
            n_err += 1
            torch.cuda.empty_cache()
            log(f"    OOM at {len(seq)} aa: {exc}", level="ERROR")
            emit({"id": jid, "pdb": "",
                  "status": "error", "msg": f"oom_at_{len(seq)}aa"})
        except Exception as exc:
            n_err += 1
            log(f"    failed: {type(exc).__name__}: {exc}", level="ERROR")
            emit({"id": jid, "pdb": "",
                  "status": "error", "msg": f"{type(exc).__name__}: {exc}"})

    log(f"Done. ok={n_ok} skipped={n_skip} error={n_err}")
    return 0 if n_err == 0 else 3


if __name__ == "__main__":
    sys.exit(main())
