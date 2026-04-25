#!/usr/bin/env python3
# PEP 563: keep type hints as strings so this module loads under the legacy
# crispr_v2_indelphi env (Python 3.7, needed for scikit-learn 0.20 compatibility
# with the bundled inDelphi pickles). Without this, `tuple[...]` / `X | None`
# crash at import time with `TypeError: 'type' object is not subscriptable`.
from __future__ import annotations

"""
Module: 04_predict_indels.py
Stage [4] — Indel outcome prediction (inDelphi / Lindel).

For each guide in the curated TSV (stage [3]), predicts the distribution of
repair outcomes and annotates:
  - frameshift_fraction  : fraction of alleles predicted to cause a frameshift
  - top_indels           : JSON list of the top N (sequence, frequency) pairs
  - frameshift_flag      : "likely_KO" if frameshift_fraction >= threshold

Usage:
    python3 04_predict_indels.py \\
        --input              <curated.tsv>      \\
        --outdir             <output_dir>       \\
        --predictors         inDelphi Lindel    \\
        --indelphi-cell-type HEK293             \\
        --frameshift-threshold 0.5             \\
        --top-outcomes       5
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _log import _log  # noqa: E402

# Add the bundled inDelphi directory to sys.path so the module is importable
# regardless of whether it was pip-installed. The bundled copy lives alongside
# this script under tools/inDelphi/ and includes model pickles for both
# sklearn 0.18.1 and 0.20.0.
_INDELPHI_TOOLS = Path(__file__).parent / "tools" / "inDelphi"
if _INDELPHI_TOOLS.exists() and str(_INDELPHI_TOOLS) not in sys.path:
    sys.path.insert(0, str(_INDELPHI_TOOLS))

# Bundled Lindel clone (contains Model_weights.pkl + model_prereq.pkl).
# Prepend to sys.path so `import Lindel` resolves to this copy instead of a
# site-packages install that may ship the code without the pickles.
_LINDEL_TOOLS = Path(__file__).parent / "tools" / "Lindel"
if _LINDEL_TOOLS.exists() and str(_LINDEL_TOOLS) not in sys.path:
    sys.path.insert(0, str(_LINDEL_TOOLS))


# ─── inDelphi wrapper ─────────────────────────────────────────────────────────
# Module-level init so we don't re-initialise (and re-log errors) per guide.
# _INDELPHI_STATE is one of:
#   "ok"           — model loaded, run_inDelphi is functional
#   "unavailable"  — import or init failed once; skip silently thereafter
#   None           — not yet attempted
_INDELPHI_STATE = None

def _init_indelphi_once(cell_type: str) -> bool:
    """Attempt to import + initialise inDelphi exactly once per process.

    Returns True if the model is ready, False otherwise. Logs a single
    diagnostic line on the first failure; subsequent calls are silent.
    """
    global _INDELPHI_STATE
    if _INDELPHI_STATE in ("ok", "predict_failed"):
        return True
    if _INDELPHI_STATE == "unavailable":
        return False

    try:
        import inDelphi
    except ImportError:
        _log("[04_indels] inDelphi not installed — frameshift prediction will "
              "fall back to Lindel only.", level="WARN")
        _INDELPHI_STATE = "unavailable"
        return False

    # Upstream inDelphi uses `celltype` (no underscore). Older forks used
    # `cell_type`. Accept either so the config key remains stable.
    # Capture init_model's unprefixed stdout chatter ("Initializing model
    # aax/aag, HEK293...") so it doesn't pollute the orchestrator log; emit
    # a single timestamped line instead.
    import contextlib, io
    init_stdout = io.StringIO()
    try:
        with contextlib.redirect_stdout(init_stdout):
            try:
                inDelphi.init_model(celltype=cell_type)
            except TypeError:
                inDelphi.init_model(cell_type=cell_type)
    except AssertionError as exc:
        # Upstream init_model asserts on unsupported sklearn versions and on
        # missing bundled model pickles. Both produce AssertionError.
        _log(f"[04_indels] inDelphi unavailable: {exc}. The bundled models "
              "require scikit-learn 0.18.1 or 0.20.0 — falling back to Lindel.", level="ERROR")
        _INDELPHI_STATE = "unavailable"
        return False
    except Exception as exc:
        _log(f"[04_indels] inDelphi init failed ({type(exc).__name__}: {exc}) "
              "— falling back to Lindel.", level="ERROR")
        _INDELPHI_STATE = "unavailable"
        return False

    _INDELPHI_STATE = "ok"
    _chatter = init_stdout.getvalue().strip().replace("\n", " ")
    if _chatter:
        _log(f"[04_indels] inDelphi init: {_chatter}")
    return True


def run_inDelphi(seq30: str, cell_type: str, top_n: int) -> tuple[float, list]:
    """
    seq30: 30-nt sequence centred on the cut site
           (positions 17-20 = guide, cut between 17 and 18 from 5' end).
    Returns (frameshift_fraction, top_outcomes).
    """
    if not _init_indelphi_once(cell_type):
        return float("nan"), []

    # inDelphi rejects any character outside ACGT. CRISPOR's 23-nt targetSeq
    # padded to 30 nt with 'N' (see get_context_seq) always fails validation,
    # which historically produced a "too many values to unpack" ValueError
    # because inDelphi.predict() returns a string on validation errors and a
    # tuple on success. Reject N-padded / short contexts up front so the
    # script falls back to Lindel cleanly without a noisy per-guide error.
    if len(seq30) < 30 or any(ch not in "ACGT" for ch in seq30):
        return float("nan"), []

    import inDelphi
    try:
        # SpCas9 cuts between nt 17 and 18 of the 20-nt guide (0-based).
        # For a 30-nt context centred on the cut, cutsite = 17.
        result = inDelphi.predict(seq30, 17)
    except Exception as exc:
        global _INDELPHI_STATE
        if _INDELPHI_STATE != "predict_failed":
            _log(f"[04_indels] inDelphi.predict failed: {exc} "
                  "(further per-guide errors suppressed).", level="ERROR")
            _INDELPHI_STATE = "predict_failed"
        return float("nan"), []

    # inDelphi.predict returns a string on validation errors and a
    # (pred_df, stats) tuple on success. Guard the unpack.
    if not isinstance(result, tuple) or len(result) != 2:
        return float("nan"), []
    pred_df, stats = result

    # Upstream stats dict already contains the frameshift fraction (0-100).
    # Prefer it; fall back to computing from pred_df if absent.
    fs_pct = stats.get("Frameshift frequency") if isinstance(stats, dict) else None
    if fs_pct is None:
        # Derive: Category=='del' with Length % 3 != 0, plus 1-bp insertions
        mask = (((pred_df["Category"] == "del") & (pred_df["Length"] % 3 != 0)) |
                 (pred_df["Category"] == "ins"))
        fs_pct = float(pred_df.loc[mask, "Predicted frequency"].sum())
    fs_frac = float(fs_pct) / 100.0

    # Build a human-readable Indel label from available columns
    def _label(row):
        if row["Category"] == "ins":
            return f"+1{row.get('Inserted Bases', '')}"
        return f"-{int(row['Length'])}"
    top_rows = pred_df.nlargest(top_n, "Predicted frequency")
    top_list = [{"indel": _label(r), "freq": round(float(r["Predicted frequency"]), 3)}
                for _, r in top_rows.iterrows()]
    return fs_frac, top_list


# ─── Lindel wrapper ───────────────────────────────────────────────────────────
# The real Lindel API is Lindel.Predictor.gen_prediction(seq, weights, prereq),
# where `weights` and `prereq` come from the bundled pickles. The upstream repo
# never exposes a top-level Lindel.run_lindel() or Lindel.predict(), so probing
# for those attributes always fails. Load the pickles once per process and
# call gen_prediction directly.
_LINDEL_STATE = None     # "ok" | "unavailable"
_LINDEL_WEIGHTS = None
_LINDEL_PREREQ = None


def _init_lindel_once() -> bool:
    """Load Lindel weights + prereq pickles exactly once per process."""
    global _LINDEL_STATE, _LINDEL_WEIGHTS, _LINDEL_PREREQ
    if _LINDEL_STATE == "ok":
        return True
    if _LINDEL_STATE == "unavailable":
        return False

    try:
        import pickle
        import Lindel  # noqa: F401 — needed to resolve Lindel.__path__
        from Lindel.Predictor import gen_prediction  # noqa: F401 — API smoke test
        # Locate the model pickles. Prefer the import's own __path__ (site-
        # packages install may carry the pickles); fall back to the bundled
        # clone under tools/Lindel/ which is guaranteed to ship them.
        candidate_roots = [Path(p) for p in Lindel.__path__]
        bundled_pkg = _LINDEL_TOOLS / "Lindel"
        if bundled_pkg.exists() and bundled_pkg not in candidate_roots:
            candidate_roots.append(bundled_pkg)

        weights_path = next(
            (r / "Model_weights.pkl" for r in candidate_roots
             if (r / "Model_weights.pkl").exists()), None)
        prereq_path = next(
            (r / "model_prereq.pkl" for r in candidate_roots
             if (r / "model_prereq.pkl").exists()), None)
        if weights_path is None or prereq_path is None:
            raise FileNotFoundError(
                str(candidate_roots[0] / "Model_weights.pkl"))
        with open(weights_path, "rb") as fh:
            _LINDEL_WEIGHTS = pickle.load(fh)
        with open(prereq_path, "rb") as fh:
            _LINDEL_PREREQ = pickle.load(fh)
    except ImportError:
        _log("[04_indels] Lindel not installed — skipping.", level="WARN")
        _LINDEL_STATE = "unavailable"
        return False
    except FileNotFoundError as exc:
        _log(f"[04_indels] Lindel pickles missing ({exc.filename}) — skipping.",
             level="WARN")
        _LINDEL_STATE = "unavailable"
        return False
    except Exception as exc:
        _log(f"[04_indels] Lindel init failed ({type(exc).__name__}: {exc}) "
              "— skipping.", level="WARN")
        _LINDEL_STATE = "unavailable"
        return False

    _LINDEL_STATE = "ok"
    return True


def run_Lindel(seq60: str, top_n: int) -> tuple[float, list]:
    """
    seq60: 60-nt sequence (30 nt upstream + 30 nt downstream of cut site).
           Lindel requires a valid PAM (NGG) at positions 33-36; sequences
           lacking this are rejected silently.
    Returns (frameshift_fraction, top_outcomes).
    """
    if not _init_lindel_once():
        return float("nan"), []

    # Lindel's gen_prediction rejects non-ACGT characters and requires a
    # valid PAM at positions 33-36. Reject up front so N-padded contexts
    # don't generate misleading scores.
    if len(seq60) < 60 or any(ch not in "ACGT" for ch in seq60[:60]):
        return float("nan"), []
    if seq60[33:36] not in ("AGG", "TGG", "CGG", "GGG"):
        return float("nan"), []

    try:
        from Lindel.Predictor import gen_prediction
        y_hat, fs = gen_prediction(seq60[:60], _LINDEL_WEIGHTS, _LINDEL_PREREQ)
    except Exception as exc:
        global _LINDEL_STATE
        if _LINDEL_STATE != "predict_failed":
            _log(f"[04_indels] Lindel.gen_prediction failed: {exc} "
                  "(further per-guide errors suppressed).", level="ERROR")
            _LINDEL_STATE = "predict_failed"
        return float("nan"), []

    # Map indices back to indel labels via the reverse index in prereq[1],
    # mirroring the layout used by the bundled Lindel_prediction.py helper.
    rev_index = _LINDEL_PREREQ[1]
    pred_freq = {rev_index[i]: float(y_hat[i]) for i in range(len(y_hat))
                 if y_hat[i] != 0}
    if not pred_freq:
        return float(fs), []

    sorted_items = sorted(pred_freq.items(), key=lambda kv: kv[1], reverse=True)[:top_n]
    top_list = [{"indel": k, "freq": round(v, 3)} for k, v in sorted_items]
    return float(fs), top_list


# ─── Sequence helpers ─────────────────────────────────────────────────────────

def get_context_seq(row: dict, length: int = 30) -> str:
    """Try to get cut-site context from CRISPOR columns; return N-padded if absent.

    NOTE: CRISPR-P v2.0 exports only the 23-nt protospacer+PAM, so this
    fallback always N-pads, which both inDelphi and Lindel reject. When a
    genome FASTA + GTF are available, prefer build_cut_context() which
    returns a true genomic flank centered on the cut site.
    """
    for col in ("targetSeq", "guideSeq", "sequence", "Sequence"):
        if col in row and row[col]:
            seq = row[col].strip().upper()
            if len(seq) < length:
                seq = seq.ljust(length, "N")
            return seq[:length]
    return "N" * length


# ─── Genome-aware cut-site context (CRISPR-P v2.0) ────────────────────────────
# CRISPR-P exports only the 23-nt protospacer+PAM, which is too short for
# inDelphi (30 nt) or Lindel (60 nt). Both predictors silently return NaN on
# N-padded input, so every guide ends up with fs_frac=NA and the composite
# KO score collapses to a flat baseline. To recover real predictions, load
# the gene's genomic CDS region once per stage-[4] invocation and locate
# each guide by sequence match (forward + reverse-complement).

def _reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def _load_gene_cds(genome_fasta: str, gtf: str, gene_group: str,
                   gene_id: str) -> tuple[str, str]:
    """Return (cds_sense_seq, strand) for the target gene, or ("", "")."""
    if not (genome_fasta and gtf and gene_id):
        return "", ""
    try:
        # Reuse stage-[5]'s GTF + CDS extraction so both stages stay in sync.
        sys.path.insert(0, str(Path(__file__).parent))
        rebuild = __import__("05_rebuild_transcripts")
    except Exception as exc:
        _log(f"[04_indels] Could not import 05_rebuild_transcripts for context "
              f"extraction: {exc}", level="WARN")
        return "", ""

    coords = rebuild.load_cds_coords(gtf, gene_group, gene_id)
    if not coords:
        return "", ""
    # Pick the first transcript; all transcripts of one gene share strand.
    first_tx = next(iter(coords))
    exons = coords[first_tx]
    strand = exons[0][3] if exons else "+"
    cds = rebuild.extract_cds_sequence(genome_fasta, exons, strand)
    return cds, strand


def build_cut_context(rows: list[dict], length: int,
                      cds_sense: str) -> list[str]:
    """Return one length-nt sense-strand context per row, centred on the cut.

    Locates each guide's 20-nt protospacer in the sense-strand CDS (forward
    or reverse-complemented), computes the SpCas9 cut at +17 from the
    protospacer start, and slices ±length/2 around it. Returns an empty
    string when the guide cannot be unambiguously located or the flanking
    window would run off the CDS — both predictors treat empty strings as
    "skip" and emit NaN, which is preferable to fabricating context.

    The 30-nt window for inDelphi places the cut at index 17 (17 nt upstream
    + 13 nt downstream); the 60-nt window for Lindel places the cut at
    index 30 (30 nt upstream + 30 nt downstream, with PAM at 33-35).
    """
    if not cds_sense:
        return ["" for _ in rows]
    half_up = 17 if length == 30 else 30
    half_dn = length - half_up
    contexts: list[str] = []
    for row in rows:
        proto = (row.get("protospacer") or "").strip().upper()
        if not proto:
            # Fallback to the 23-nt Sequence column (strip PAM).
            for col in ("Sequence", "sequence", "targetSeq", "guideSeq"):
                v = (row.get(col) or "").strip().upper()
                if len(v) >= 20:
                    proto = v[:20]
                    break
        if len(proto) < 20:
            contexts.append("")
            continue
        guide = proto[:20]
        idx = cds_sense.find(guide)
        cut = -1
        if idx != -1:
            cut = idx + 17
        else:
            rc = _reverse_complement(guide)
            idx = cds_sense.find(rc)
            if idx != -1:
                # Guide is on the antisense strand: cut on sense strand sits
                # 3 nt from the 5' end of the rc-match.
                cut = idx + 3
        if cut < half_up or cut + half_dn > len(cds_sense):
            contexts.append("")
            continue
        contexts.append(cds_sense[cut - half_up : cut + half_dn])
    return contexts


# ─── Legacy-env inDelphi dispatcher ──────────────────────────────────────────

def _run_indelphi_in_legacy_env(input_path: str, cell_type: str,
                                 top_n: int, conda_env: str,
                                 ctx30: list[str] | None = None) -> list[dict] | None:
    """Spawn this script under conda_env as a worker to run inDelphi only.

    Returns a list of per-row dicts with keys inDelphi_fs_frac and
    inDelphi_top_indels, in the same row order as input_path. Returns None
    if the subprocess fails.

    When ctx30 is provided (one 30-nt sense-strand context per row), writes
    a temp TSV with an extra __ctx30__ column and passes that to the worker
    instead of input_path. The worker prefers __ctx30__ over the fallback
    get_context_seq path so the legacy env doesn't need to import the GTF
    parser (which uses PEP-585 generics incompatible with Python 3.7).
    """
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
        tmp_json = tmp.name

    # When genome-derived contexts are available, materialize a temp TSV with
    # an extra __ctx30__ column so the worker can use real cut-site flanks
    # without re-parsing GTF/genome under Python 3.7.
    worker_input = input_path
    tmp_input: str | None = None
    if ctx30 is not None:
        with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tsv", delete=False, newline="") as tin:
            with open(input_path, newline="") as fin:
                rows = list(csv.DictReader(fin, delimiter="\t"))
            fields = list(rows[0].keys()) if rows else []
            new_fields = fields + ["__ctx30__"]
            writer = csv.DictWriter(
                tin, fieldnames=new_fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for i, row in enumerate(rows):
                row["__ctx30__"] = ctx30[i] if i < len(ctx30) else ""
                writer.writerow(row)
            tmp_input = tin.name
        worker_input = tmp_input

    cmd = [
        "conda", "run", "-n", conda_env, "--no-capture-output",
        "python3", __file__,
        "--input",               worker_input,
        "--outdir",              "/tmp",    # unused in worker mode
        "--predictors",          "inDelphi",
        "--indelphi-cell-type",  cell_type,
        "--top-outcomes",        str(top_n),
        "--indelphi-json-out",   tmp_json,
    ]
    # Pass the bundled inDelphi path through PYTHONPATH so the worker process
    # (which runs under a different conda env) can also import it without a
    # pip install step.
    proc_env = os.environ.copy()
    proc_env["PYTHONPATH"] = (
        str(_INDELPHI_TOOLS) + os.pathsep + proc_env.get("PYTHONPATH", "")
    )
    result = subprocess.run(cmd, capture_output=False, env=proc_env)
    if result.returncode != 0:
        _log(f"[04_indels] Legacy inDelphi subprocess failed (rc={result.returncode}). "
              "inDelphi results will be NA.", level="ERROR")
        if tmp_input:
            Path(tmp_input).unlink(missing_ok=True)
        return None

    try:
        with open(tmp_json) as fh:
            return json.load(fh)
    except Exception as exc:
        _log(f"[04_indels] Could not read inDelphi worker output: {exc}", level="ERROR")
        return None
    finally:
        Path(tmp_json).unlink(missing_ok=True)
        if tmp_input:
            Path(tmp_input).unlink(missing_ok=True)


# ─── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Indel outcome prediction")
    p.add_argument("--input",                  required=True)
    p.add_argument("--outdir",                 required=True)
    p.add_argument("--predictors",             nargs="+", default=["inDelphi", "Lindel"])
    p.add_argument("--indelphi-cell-type",     default="HEK293")
    p.add_argument("--frameshift-threshold",   type=float, default=0.5)
    p.add_argument("--top-outcomes",           type=int,   default=5)
    # Python 3.7 (legacy inDelphi env) lacks argparse.BooleanOptionalAction, so
    # expose --overwrite / --no-overwrite as two separate actions on the same
    # dest. Behaviour matches BooleanOptionalAction for both envs.
    p.add_argument("--workers",      type=int, default=1,
                   help="Reserved for future parallel indel prediction (accepted but unused).")
    p.add_argument("--overwrite",    dest="overwrite", action="store_true",  default=True)
    p.add_argument("--no-overwrite", dest="overwrite", action="store_false")
    # Legacy-env support: when inDelphi is unavailable in the current env,
    # re-invoke this script under the specified conda env (which ships
    # scikit-learn 0.18/0.20) and merge the results back.
    p.add_argument("--indelphi-conda-env",     default="",
                   help="Conda env name with legacy scikit-learn for inDelphi.")
    # Internal worker flag: write only inDelphi JSON results to this path, then exit.
    p.add_argument("--indelphi-json-out",      default="",
                   help=argparse.SUPPRESS)
    # Optional genome context for CRISPR-P inputs (which only carry the 23-nt
    # protospacer+PAM). When provided, real cut-site flanks are extracted and
    # passed to inDelphi/Lindel; otherwise both predictors fall back to the
    # N-padded path and emit NaN.
    p.add_argument("--genome-fasta",           default="",
                   help="Genome FASTA for cut-site context extraction.")
    p.add_argument("--gtf",                    default="",
                   help="GTF annotation matching --genome-fasta.")
    p.add_argument("--gene-group",             default="DMP",
                   help="Family tag used as a fallback regex when --gene-id misses.")
    p.add_argument("--gene-id",                default="",
                   help="Specific target gene id (e.g. SMEL5_01g008730).")
    return p.parse_args()


def main():
    args = parse_args()

    inpath = Path(args.input)

    # ── Worker mode: run inDelphi only, write JSON results, exit ─────────────
    if args.indelphi_json_out:
        # Exit rc=1 immediately if inDelphi cannot initialize — the caller
        # (_run_indelphi_in_legacy_env) interprets a non-zero return code as a
        # clean failure and falls back to Lindel-only in the main env. This
        # prevents the caller from logging "inDelphi results loaded" when the
        # actual predictions would all be NA (silent sklearn mismatch).
        if not _init_indelphi_once(args.indelphi_cell_type):
            sys.exit(1)
        with open(inpath, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        results = []
        for row in rows:
            # Prefer the precomputed __ctx30__ column (genomic flank from the
            # main env) over the N-padded fallback, which both predictors
            # silently reject for CRISPR-P inputs.
            ctx = (row.get("__ctx30__") or "").strip().upper()
            seq30 = ctx if len(ctx) == 30 else get_context_seq(row, 30)
            fs, top = run_inDelphi(seq30, args.indelphi_cell_type, args.top_outcomes)
            results.append({
                "inDelphi_fs_frac":    f"{fs:.4f}" if fs == fs else "NA",
                "inDelphi_top_indels": json.dumps(top),
            })
        with open(args.indelphi_json_out, "w") as fh:
            json.dump(results, fh)
        return

    # ── Normal mode ───────────────────────────────────────────────────────────
    outdir  = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    gene_stem = inpath.stem.split(".")[0]
    outpath = outdir / f"{gene_stem}.indels.tsv"

    if not args.overwrite and outpath.exists():
        _log(f"[04_indels] Skipping (overwrite=false): {outpath}", level="WARN")
        return

    with open(inpath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows   = list(reader)
        fields = list(reader.fieldnames or [])

    new_cols = []
    for predictor in args.predictors:
        new_cols += [f"{predictor}_fs_frac", f"{predictor}_top_indels"]
    new_cols += ["best_fs_frac", "frameshift_flag"]
    new_fields = fields + new_cols

    # Build per-row genomic cut-site contexts when genome+GTF are provided.
    # This is the only source of real flanking sequence for CRISPR-P inputs;
    # without it both predictors silently reject the N-padded 23-nt fallback
    # and every guide ends up with fs_frac=NA. ctx30/ctx60 are sense-strand
    # 30-nt and 60-nt windows centred on the SpCas9 cut site (3 nt from PAM).
    ctx30_list: list[str] = []
    ctx60_list: list[str] = []
    if args.genome_fasta and args.gtf and args.gene_id:
        cds_sense, _strand = _load_gene_cds(
            args.genome_fasta, args.gtf, args.gene_group, args.gene_id)
        if cds_sense:
            ctx30_list = build_cut_context(rows, 30, cds_sense)
            ctx60_list = build_cut_context(rows, 60, cds_sense)
            n30 = sum(1 for s in ctx30_list if s)
            n60 = sum(1 for s in ctx60_list if s)
            _log(f"[04_indels] Cut-site context built for {args.gene_id}: "
                  f"30-nt={n30}/{len(rows)}, 60-nt={n60}/{len(rows)}.", level="INFO")
        else:
            _log(f"[04_indels] CDS extraction failed for {args.gene_id}; "
                  f"falling back to N-padded context (predictors will return NA).",
                  level="WARN")

    # Pre-fetch inDelphi results.
    # When --indelphi-conda-env is provided the orchestrator has already
    # verified that env has the pinned sklearn 0.20 that inDelphi requires,
    # so go there directly.  Trying the main env first only produces a noisy
    # ERROR for every parallel gene process when the main env has an
    # incompatible sklearn (e.g. 1.7.x).
    indelphi_prefetch: list[dict] | None = None
    if "inDelphi" in args.predictors:
        if args.indelphi_conda_env:
            indelphi_prefetch = _run_indelphi_in_legacy_env(
                str(inpath), args.indelphi_cell_type,
                args.top_outcomes, args.indelphi_conda_env,
                ctx30=ctx30_list if ctx30_list else None,
            )
            if indelphi_prefetch is not None:
                _log(f"[04_indels] inDelphi results loaded from legacy env "
                      f"({len(indelphi_prefetch)} guides).", level="INFO")
            else:
                _log(f"[04_indels] inDelphi unavailable in legacy env "
                      f"'{args.indelphi_conda_env}'; Lindel-only scores used.", level="WARN")
        else:
            available = _init_indelphi_once(args.indelphi_cell_type)
            if not available:
                _log("[04_indels] inDelphi unavailable in current env and no "
                      "legacy env configured; Lindel-only scores used.", level="WARN")

    for i, row in enumerate(rows):
        best_fs = float("nan")

        if "inDelphi" in args.predictors:
            if indelphi_prefetch is not None:
                # Results came from the legacy subprocess.
                pf = indelphi_prefetch[i] if i < len(indelphi_prefetch) else {}
                row["inDelphi_fs_frac"]    = pf.get("inDelphi_fs_frac", "NA")
                row["inDelphi_top_indels"] = pf.get("inDelphi_top_indels", "[]")
                fs_val = row["inDelphi_fs_frac"]
                if fs_val != "NA":
                    try:
                        fs = float(fs_val)
                        if fs == fs:
                            best_fs = fs if (best_fs != best_fs or fs > best_fs) else best_fs
                    except ValueError:
                        pass
            else:
                ctx30 = ctx30_list[i] if i < len(ctx30_list) else ""
                seq30 = ctx30 if len(ctx30) == 30 else get_context_seq(row, 30)
                fs, top = run_inDelphi(seq30, args.indelphi_cell_type, args.top_outcomes)
                row["inDelphi_fs_frac"]    = f"{fs:.4f}" if fs == fs else "NA"
                row["inDelphi_top_indels"] = json.dumps(top)
                if fs == fs:
                    best_fs = fs if (best_fs != best_fs or fs > best_fs) else best_fs

        if "Lindel" in args.predictors:
            # Prefer the Lindel-Score column CRISPOR computes with proper
            # genomic flanking context. Fall back to standalone Lindel only
            # when the column is absent or flagged NotEnoughFlankSeq (which
            # happens when --skipAlign was used or the locus is near contig
            # edges). Standalone Lindel needs 60-nt context with PAM at
            # positions 33-36; the 23-nt targetSeq from CRISPOR produces
            # unreliable results when padded with N.
            crispor_lindel = row.get("Lindel-Score", "")
            ctx60 = ctx60_list[i] if i < len(ctx60_list) else ""
            if crispor_lindel and crispor_lindel not in (
                    "", "NA", "N/A", "NotEnoughFlankSeq"):
                try:
                    # CRISPOR reports Lindel-Score as a percentage (0-100).
                    fs = float(crispor_lindel) / 100.0
                    top = []
                except ValueError:
                    seq60 = ctx60 if len(ctx60) == 60 else get_context_seq(row, 60)
                    fs, top = run_Lindel(seq60, args.top_outcomes)
            else:
                seq60 = ctx60 if len(ctx60) == 60 else get_context_seq(row, 60)
                fs, top = run_Lindel(seq60, args.top_outcomes)
            row["Lindel_fs_frac"]    = f"{fs:.4f}" if fs == fs else "NA"
            row["Lindel_top_indels"] = json.dumps(top)
            if fs == fs:
                best_fs = fs if (best_fs != best_fs or fs > best_fs) else best_fs

        row["best_fs_frac"]    = f"{best_fs:.4f}" if best_fs == best_fs else "NA"
        row["frameshift_flag"] = "likely_KO" if (best_fs == best_fs and
                                                  best_fs >= args.frameshift_threshold) else ""

    with open(outpath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    ko_n = sum(1 for r in rows if r.get("frameshift_flag") == "likely_KO")
    _log(f"[04_indels] {len(rows)} guides; {ko_n} likely_KO (fs>={args.frameshift_threshold}) -> {outpath}", level="INFO")


if __name__ == "__main__":
    main()
