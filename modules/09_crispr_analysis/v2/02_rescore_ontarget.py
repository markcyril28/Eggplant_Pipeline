#!/usr/bin/env python3
"""
Module: 02_rescore_ontarget.py
Stage [2] — On-target rescoring with CRISPRon and/or DeepSpCas9.

Reads the filtered CRISPOR TSV from stage [1], adds rescored columns
(crisprOn_score, DeepSpCas9_score), and writes an updated TSV.

DeepSpCas9 (Kim et al. 2019) requires a 30-nt context:
    4 nt upstream + 20-nt guide + 3-nt PAM + 3 nt downstream.
CRISPOR emits only a 23-nt `targetSeq` (guide + PAM). To recover the
flanking 4 + 3 = 7 nt, this module looks up the target nucleotide FASTA
that was fed to CRISPOR in stage [1]. When `--target-fasta` is supplied,
DeepSpCas9 receives the real genomic/target context; when it is absent or
the lookup fails, the module falls back to N-padding AND logs a loud
warning because N-padding systematically depresses DeepSpCas9 scores.

Usage:
    python3 02_rescore_ontarget.py \\
        --input         <crispor_filtered.tsv> \\
        --outdir        <output_dir>           \\
        --target-fasta  <gene_or_locus.fa>     \\
        --predictors    crisprOn DeepSpCas9    \\
        --flag-threshold 0.3
"""

import argparse
import csv
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# Ensure the bundled crisprOn shim is importable regardless of PYTHONPATH
_TOOLS_DIR = Path(__file__).parent / "tools"
if str(_TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(_TOOLS_DIR))

sys.path.insert(0, str(Path(__file__).parent))
from _log import _log  # noqa: E402


# ─── Helper: extract 20-nt protospacer from CRISPOR row ──────────────────────

def get_protospacer(row: dict) -> str:
    """Return the 20-nt guide sequence from a CRISPOR output row."""
    for col in ("targetSeq", "guideSeq", "protospacer", "sequence", "Sequence"):
        if col in row and row[col]:
            seq = row[col].strip().upper()
            # CRISPOR often returns 23-nt (20-nt guide + 3-nt PAM); take first 20
            return seq[:20]
    return ""


def _get_targetseq(row: dict) -> str:
    """Return the 23-nt (guide + PAM) sequence from a CRISPOR row."""
    for col in ("targetSeq", "guideSeq", "sequence", "Sequence"):
        if col in row and row[col]:
            return row[col].strip().upper().replace("U", "T")
    return ""


# ─── Target FASTA (genomic/CDS context for DeepSpCas9) ───────────────────────

def _load_target_fasta(path: str) -> str:
    """Concatenate every record in a FASTA file into one uppercase string.

    For plant CRISPR, each stage-1 input FASTA is a single gene/CDS record
    (e.g. SMEL5_01g008730.fa), so concatenation is effectively a no-op. When
    the file holds multiple records the concatenated string still serves as
    a valid search substrate — the target-seq match must be unique to hit.
    """
    try:
        seq_chunks: list[str] = []
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    # Separator so records do not fuse across headers, which
                    # could create a spurious match spanning two records.
                    seq_chunks.append("N")
                    continue
                seq_chunks.append(line.strip().upper().replace("U", "T"))
        return "".join(seq_chunks)
    except FileNotFoundError:
        return ""
    except Exception as exc:
        _log(f"[02_rescore] Could not read target FASTA {path}: {exc}", level="WARN")
        return ""


def _rc(seq: str) -> str:
    """Minimal reverse-complement (A/T/G/C/N only)."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def build_deepspcas9_context(row: dict, target_fasta_seq: str) -> tuple[str, str]:
    """Build the 30-nt DeepSpCas9 context for a single CRISPOR row.

    Returns (context_30nt, source) where source ∈ {'fasta', 'n_padded'}.
    The fasta source is used when the 23-nt targetSeq can be uniquely located
    in the target FASTA on either strand; otherwise we fall back to N-padding
    the 23-nt column (systematically depresses the score — logged as WARN
    by the caller so the thesis Methods can report the fallback rate).
    """
    target_23 = _get_targetseq(row)
    if len(target_23) < 23:
        return (target_23.ljust(30, "N")[:30], "n_padded")
    if not target_fasta_seq:
        return (target_23.ljust(30, "N")[:30], "n_padded")

    # Forward strand — expect exactly one hit for a CRISPOR-emitted guide.
    pos = target_fasta_seq.find(target_23)
    if pos != -1 and pos >= 4 and pos + 23 + 3 <= len(target_fasta_seq):
        return (target_fasta_seq[pos - 4 : pos + 23 + 3], "fasta")

    # Reverse strand — guide was called on the minus strand of the input.
    rc_target = _rc(target_23)
    pos = target_fasta_seq.find(rc_target)
    if pos != -1 and pos >= 3 and pos + 23 + 4 <= len(target_fasta_seq):
        # Extract on the plus strand then reverse-complement so the returned
        # context is in the canonical (protospacer + PAM) orientation.
        plus_window = target_fasta_seq[pos - 3 : pos + 23 + 4]
        return (_rc(plus_window), "fasta")

    # Match found but too close to the FASTA edge for the full 30-nt window,
    # or no unique match. Fall back to N-padding.
    return (target_23.ljust(30, "N")[:30], "n_padded")


# ─── CRISPRon wrapper ─────────────────────────────────────────────────────────

def score_crisprOn(seqs: list[str], workers: int = 1) -> list[float]:
    try:
        import crisprOn

        def _predict_one(seq: str) -> float:
            try:
                sc = crisprOn.predict(seq)
                return float(sc) if sc is not None else float("nan")
            except Exception:
                return float("nan")

        if workers > 1 and len(seqs) > 1:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                return list(pool.map(_predict_one, seqs))
        return [_predict_one(seq) for seq in seqs]
    except ImportError:
        _log("[02_rescore] crisprOn not installed — skipping.", level="WARN")
        return [float("nan")] * len(seqs)


# ─── DeepCRISPR wrapper (plant-usable; mammalian-trained, data-efficient) ──
# Reference: Chuai et al. 2018, doi:10.1186/s13059-018-1459-4
# Upstream repo: https://github.com/bm2-lab/DeepCRISPR
#
# DeepCRISPR expects 23-nt guide+PAM input and returns an activity score in
# [0, 1]. The upstream package pins TensorFlow 1.x and sonnet 1.x, which do
# not coexist with the modern `crispr_v3` env. This wrapper therefore
# prefers a subprocess into an isolated DeepCRISPR conda env
# (CONDA_ENV_DEEPCRISPR, if defined) and only falls back to direct import
# as a last resort. If neither path is available, scores are NaN.
#
# To wire DeepCRISPR fully:
#   1. Clone into modules/09_crispr_analysis/v2/tools/DeepCRISPR.
#   2. Create a dedicated env `crispr_v3_deepcrispr`:
#        conda create -n crispr_v3_deepcrispr python=3.6 tensorflow=1.8 \
#            dm-sonnet=1.13 numpy pandas
#      then register the clone on its sys.path via a .pth file.
#   3. Export CONDA_ENV_DEEPCRISPR=crispr_v3_deepcrispr before invoking the
#      v3 orchestrator, or set [crispr_v3.plant_scorer].deepcrispr_env in
#      the TOML (the orchestrator threads it through).

_TOOLS_ROOT = Path(__file__).parent / "tools"


def _deepcrispr_dir() -> Path | None:
    """Return the path to a cloned DeepCRISPR tree, or None if absent.

    Accepts both the canonical "DeepCRISPR" name and "DeepCRISPR-master"
    (GitHub ZIP-extraction default) so either install path works.
    """
    for name in ("DeepCRISPR", "DeepCRISPR-master"):
        d = _TOOLS_ROOT / name
        if d.is_dir() and (d / "run_examples.py").exists():
            return d
    return None


def score_DeepCRISPR(seqs23: list[str], workers: int = 1,
                      conda_env: str = "") -> list[float]:
    """Plant-usable on-target rescorer (Chuai 2018).

    Inputs:  23-nt protospacer+PAM sequences.
    Returns: one activity score in [0, 1] per input (NaN on any failure).

    Invocation preference:
      1. Subprocess into `conda_env` (typically `crispr_v3_deepcrispr`)
         so TF 1.x / sonnet 1.x dependencies stay isolated.
      2. Direct Python import (only works if the main env somehow has
         tensorflow 1.x + dm-sonnet 1.x).
    """
    tool_dir = _deepcrispr_dir()
    if tool_dir is None:
        _log("[02_rescore] DeepCRISPR clone not found under "
             f"{_TOOLS_ROOT}/DeepCRISPR — returning NaN. "
             "Install by running setup_conda_crispr_v3.sh --with-plant-scorers "
             "or clone the repo manually.", level="WARN")
        return [float("nan")] * len(seqs23)

    # Locate the v3-bundled inference helper. It is installed by
    # setup_conda_crispr_v3.sh --with-plant-scorers alongside the clone,
    # and exposes a simple text-in / text-out CLI suitable for subprocess.
    infer_script = tool_dir / "deepcrispr_infer.py"
    if not infer_script.exists():
        _log(f"[02_rescore] DeepCRISPR clone present at {tool_dir} but the v3 "
             "inference helper deepcrispr_infer.py is missing. Run "
             "setup_conda_crispr_v3.sh --with-plant-scorers to install it, "
             "or copy it manually from the setup script's 'Step 4c' block.",
             level="WARN")
        return [float("nan")] * len(seqs23)

    # ── Path 1: subprocess under an isolated conda env (recommended) ────
    if conda_env:
        try:
            import subprocess
            import tempfile
            with tempfile.NamedTemporaryFile(suffix=".txt", mode="w",
                                              delete=False) as tmp_in:
                for s in seqs23:
                    tmp_in.write(s.strip().upper() + "\n")
                tmp_in_path = tmp_in.name
            tmp_out_path = tmp_in_path + ".scores"

            cmd = ["conda", "run", "-n", conda_env, "--no-capture-output",
                    "python", str(infer_script),
                    "--input", tmp_in_path,
                    "--output", tmp_out_path,
                    "--model-dir", str(tool_dir / "trained_models" /
                                          "ontar_cnn_reg_seq"),
                    "--seq-only"]
            try:
                r = subprocess.run(cmd, capture_output=True, text=True,
                                    timeout=1200, cwd=str(tool_dir))
            except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
                _log(f"[02_rescore] DeepCRISPR subprocess failed to launch: "
                     f"{exc}", level="WARN")
                r = None

            scores: list[float] = [float("nan")] * len(seqs23)
            if r is not None and r.returncode == 0 and Path(tmp_out_path).exists():
                with open(tmp_out_path) as fh:
                    parsed = [float(l.strip()) for l in fh if l.strip()]
                if len(parsed) == len(seqs23):
                    scores = parsed
                else:
                    _log(f"[02_rescore] DeepCRISPR output row count "
                         f"({len(parsed)}) != input ({len(seqs23)}); "
                         "returning NaN.", level="WARN")
            elif r is not None and r.returncode != 0:
                _log(f"[02_rescore] deepcrispr_infer.py exited rc={r.returncode}. "
                     "TF-1.x env likely not set up (need tensorflow=1.3, "
                     "dm-sonnet=1.9, python=3.6). stderr head: "
                     f"{r.stderr[:200]}", level="WARN")

            try:
                Path(tmp_in_path).unlink(missing_ok=True)
                Path(tmp_out_path).unlink(missing_ok=True)
            except Exception:
                pass
            return scores
        except Exception as exc:
            _log(f"[02_rescore] DeepCRISPR subprocess dispatch failed: {exc} — "
                 "falling back to direct import.", level="WARN")

    # ── Path 2: direct Python import (only works if the main env somehow
    #            has TF 1.x + sonnet 1.x — unlikely under crispr_v3). Skip
    #            this path entirely when conda_env is empty: under the
    #            modern crispr_v3 env the import always fails, which just
    #            produces a noisy WARN every time the stage runs. One
    #            INFO at dispatch time is enough; the setup instructions
    #            are in the TOML and stage-[1] README.
    if not conda_env:
        _log("[02_rescore] DeepCRISPR skipped: no [crispr_v3.plant_scorer]."
             "deepcrispr_env configured. Create an isolated TF-1.x env and "
             "set the TOML key to enable plant-trained rescoring.",
             level="INFO")
        return [float("nan")] * len(seqs23)

    try:
        if str(tool_dir) not in sys.path:
            sys.path.insert(0, str(tool_dir))
        import deepcrispr  # noqa: F401  — upstream namespace varies

        def _predict_one(seq: str) -> float:
            try:
                sc = deepcrispr.predict(seq.strip().upper()[:23])
                return float(sc) if sc is not None else float("nan")
            except Exception:
                return float("nan")

        if workers > 1 and len(seqs23) > 1:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                return list(pool.map(_predict_one, seqs23))
        return [_predict_one(s) for s in seqs23]
    except ImportError:
        _log("[02_rescore] DeepCRISPR clone present but not importable "
             f"from env '{conda_env}'. Verify the env has tensorflow=1.3 "
             "and dm-sonnet=1.9 installed.", level="WARN")
        return [float("nan")] * len(seqs23)


# ─── CRISPR-Local wrapper (plant-trained, rice-calibrated) ─────────────────
# Reference: Sun et al. 2019, doi:10.1093/bib/bbz110
# Upstream: NOT on GitHub — distributed as a tarball from
#           http://crispr.hzau.edu.cn/CRISPR-Local/
#
# CRISPR-Local is a Perl + R toolkit with a Python scoring helper. Its
# `sgRNA_Efficiency.py` (or `SgRNA_Efficiency.py` in older releases)
# reads a list of 20/23-nt protospacers and writes scores to stdout. The
# wrapper locates that script under tools/CRISPR-Local/ and invokes it
# via subprocess so its Perl/R runtime is unaffected by conda.

def _crispr_local_dir() -> Path | None:
    """Return the path to a CRISPR-Local tree, or None if absent.

    Accepts several canonical folder names because Windows NTFS can
    transiently block renaming from the ZIP-extraction default
    ("CRISPR-Local.new" during an `unzip -> mv` sequence) to the clean
    "CRISPR-Local" name. Check in preference order.
    """
    for name in ("CRISPR-Local", "CRISPR-Local.new", "CRISPR-Local-master",
                  "crispr_local"):
        d = _TOOLS_ROOT / name
        # rs2_score_calculator.py is CRISPR-Local's canonical scoring entry
        # point (Doench Rule Set 2 / Azimuth); its presence indicates a
        # complete clone. Use it as the detection fingerprint.
        if d.is_dir() and (d / "Rule_Set_2_scoring_v1" /
                            "analysis" / "rs2_score_calculator.py").exists():
            return d
    return None


def _find_crispr_local_scorer(tool_dir: Path) -> Path | None:
    """Locate CRISPR-Local's on-target scorer.

    The canonical entry point is Rule_Set_2_scoring_v1/analysis/
    rs2_score_calculator.py (Doench Rule Set 2, Python 2). Returns None
    if the expected layout is missing; caller handles the NaN fallback.

    IMPORTANT (training-species provenance): Rule Set 2 is MAMMALIAN-
    trained (Doench 2016). CRISPR-Local's plant-specific value comes
    from its paralog-aware PL-search.pl + a pre-built species sgRNA
    database — NOT from a plant-trained efficacy model. The v3 pipeline
    exposes this scorer for completeness; treat its scores as a local-
    install drop-in for Azimuth, not as a plant-calibrated metric.
    """
    rs2 = tool_dir / "Rule_Set_2_scoring_v1" / "analysis" / "rs2_score_calculator.py"
    if rs2.exists():
        return rs2
    # Fallback — try any legacy filename under the tree.
    for name in ("sgRNA_Efficiency.py", "SgRNA_Efficiency.py",
                  "on_target_score.py", "score_sgrna.py",
                  "Efficiency_Score.py"):
        for p in tool_dir.rglob(name):
            return p
    return None


def score_CRISPR_Local(seqs23: list[str], workers: int = 1,
                        conda_env: str = "") -> list[float]:
    """CRISPR-Local local on-target rescorer.

    Training data provenance: Doench Rule Set 2 (mammalian). CRISPR-Local
    ships Rule Set 2 as its scoring engine; treat this as a local-install
    drop-in for Azimuth, not as a plant-calibrated metric. See the
    docstring of `_find_crispr_local_scorer` for the full note.

    Inputs:  23-nt protospacer+PAM sequences. A 30-mer of the form
             NNNN<23-mer>NNN is built from each (4 nt upstream + 23 nt
             guide+PAM + 3 nt downstream) because rs2_score_calculator.py
             expects a 30-mer context.
    Returns: one Rule Set 2 score per input, or NaN on any failure.

    rs2_score_calculator.py is Python 2. This wrapper invokes it via
    `conda run -n <conda_env>` when `conda_env` is set; otherwise it
    attempts the current env (will fail on Python 3 with syntax errors
    and a clear WARN).
    """
    # Skip dispatch entirely when no legacy env is configured. The bundled
    # rs2_score_calculator.py pulls in model_comparison.py, which in turn
    # imports pylab/matplotlib and other heavy deps that are absent from
    # the modern crispr_v3 env. Without the dedicated env the subprocess
    # always fails with an ImportError at rc=1 — one INFO at dispatch time
    # is enough, the setup instructions live in the TOML.
    if not conda_env:
        _log("[02_rescore] CRISPR-Local skipped: no [crispr_v3.plant_scorer]."
             "crispr_local_env configured. Create the dedicated env and "
             "set the TOML key to enable Rule Set 2 rescoring.",
             level="INFO")
        return [float("nan")] * len(seqs23)

    tool_dir = _crispr_local_dir()
    if tool_dir is None:
        _log("[02_rescore] CRISPR-Local clone not found. Expected one of: "
             "CRISPR-Local/, CRISPR-Local-master/, crispr_local/ under "
             f"{_TOOLS_ROOT}. Clone from "
             "https://github.com/sunjiamin0824/CRISPR-Local.git (or "
             "download the tarball from http://crispr.hzau.edu.cn/CRISPR-Local/). "
             "Returning NaN.", level="WARN")
        return [float("nan")] * len(seqs23)

    scorer = _find_crispr_local_scorer(tool_dir)
    if scorer is None:
        _log(f"[02_rescore] CRISPR-Local clone at {tool_dir} has no "
             "Rule_Set_2_scoring_v1/analysis/rs2_score_calculator.py. "
             "Upstream layout may have changed — returning NaN.",
             level="WARN")
        return [float("nan")] * len(seqs23)

    try:
        import subprocess
        import tempfile

        # rs2_score_calculator.py reads a FASTA-like file
        # (alternating >name\n30mer\n lines) and writes >name#score\n20mer
        # for each entry. Build that input here.
        with tempfile.NamedTemporaryFile(suffix=".fa", mode="w",
                                          delete=False) as tmp_in:
            for i, s in enumerate(seqs23):
                s23 = s.strip().upper().ljust(23, "N")[:23]
                mer30 = "NNNN" + s23 + "NNN"  # 4 + 23 + 3 = 30
                tmp_in.write(f">g{i}\n{mer30}\n")
            tmp_in_path = tmp_in.name
        tmp_out_path = tmp_in_path + ".scores"

        cmd = ["python", str(scorer), "--input", tmp_in_path,
                "--output", tmp_out_path]
        if conda_env:
            cmd = ["conda", "run", "-n", conda_env,
                    "--no-capture-output"] + cmd
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=600, cwd=str(tool_dir))
        except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
            _log(f"[02_rescore] CRISPR-Local subprocess failed to launch: "
                 f"{exc}", level="WARN")
            return [float("nan")] * len(seqs23)

        scores: list[float] = [float("nan")] * len(seqs23)
        if r.returncode != 0:
            _log(f"[02_rescore] rs2_score_calculator.py exited rc={r.returncode}. "
                 "Python-2 env likely not set — pass its name via "
                 "[crispr_v3.plant_scorer].crispr_local_env. "
                 f"stderr: {r.stderr[:200]}", level="WARN")
        elif Path(tmp_out_path).exists():
            # Parse ">gN#SCORE\n20mer" — one entry per guide.
            with open(tmp_out_path) as fh:
                current_idx = None
                for line in fh:
                    line = line.strip()
                    if line.startswith(">"):
                        m = line[1:].split("#", 1)
                        if len(m) == 2 and m[0].startswith("g"):
                            try:
                                current_idx = int(m[0][1:])
                                sc = float(m[1])
                                if 0 <= current_idx < len(scores):
                                    scores[current_idx] = sc
                            except ValueError:
                                pass

        # Cleanup
        try:
            Path(tmp_in_path).unlink(missing_ok=True)
            Path(tmp_out_path).unlink(missing_ok=True)
        except Exception:
            pass
        return scores
    except Exception as exc:
        _log(f"[02_rescore] CRISPR-Local invocation failed: {exc}",
             level="WARN")
        return [float("nan")] * len(seqs23)


# ─── DeepSpCas9 wrapper ───────────────────────────────────────────────────────

def score_DeepSpCas9(seqs30: list[str], workers: int = 1) -> list[float]:
    """
    DeepSpCas9 (Kim et al. 2019) requires a 30-nt context:
        4 nt upstream + 20-nt guide + 3-nt PAM + 3 nt downstream.

    Input sequences MUST already be 30 nt; build them with
    `build_deepspcas9_context` before calling this function. Passing
    N-padded CRISPOR `targetSeq` (23 nt) silently tanks the score — the
    upstream model treats N's as unknowns and the 4 nt upstream flank is
    particularly informative for the prediction.

    Falls back to NaN if the package is unavailable.
    """
    try:
        import DeepSpCas9 as dsc

        def _predict_one(seq: str) -> float:
            try:
                # Defensive: trust but verify the caller.
                if len(seq) < 30:
                    seq = seq.ljust(30, "N")
                elif len(seq) > 30:
                    seq = seq[:30]
                sc = dsc.predict(seq)
                return float(sc) if sc is not None else float("nan")
            except Exception:
                return float("nan")

        if workers > 1 and len(seqs30) > 1:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                return list(pool.map(_predict_one, seqs30))
        return [_predict_one(seq) for seq in seqs30]
    except ImportError:
        _log("[02_rescore] DeepSpCas9 not installed — skipping.", level="WARN")
        return [float("nan")] * len(seqs30)


# ─── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="On-target rescoring for gRNA candidates")
    p.add_argument("--input",           required=True,  help="Filtered CRISPOR TSV from stage [1]")
    p.add_argument("--outdir",          required=True,  help="Output directory")
    p.add_argument("--target-fasta",    default="",
                   help="Target nucleotide FASTA fed to CRISPOR in stage [1]. Required for a "
                        "correctly scaled DeepSpCas9 score (provides the 4 nt upstream + 3 nt "
                        "downstream flank CRISPOR does not emit). Missing ⇒ N-padded fallback "
                        "with a loud warning.")
    p.add_argument("--predictors",      nargs="+",      default=["crisprOn", "DeepSpCas9"],
                   help="Predictors to run. Mammalian: crisprOn, DeepSpCas9. "
                        "Plant-trained: DeepCRISPR, CRISPR-Local.")
    p.add_argument("--flag-threshold",  type=float,     default=0.3,
                   help="Guides whose best rescore < threshold get flagged")
    p.add_argument("--workers",         type=int,       default=1,
                   help="Worker threads for per-guide scoring")
    p.add_argument("--deepcrispr-env",  default="",
                   help="Name of an isolated conda env holding TensorFlow 1.x "
                        "+ dm-sonnet 1.x for DeepCRISPR. When set, the module "
                        "subprocess-invokes DeepCRISPR under that env so its "
                        "deps don't collide with the modern v3 env.")
    p.add_argument("--crispr-local-env", default="",
                   help="Name of a Python-2 conda env (Python 2.7 + "
                        "scikit-learn 0.16/0.17) for CRISPR-Local's "
                        "rs2_score_calculator.py. Leave empty to attempt "
                        "invocation in the current env (will fail with a "
                        "SyntaxError on Python 3).")
    p.add_argument("--overwrite",       action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()


def main():
    args = parse_args()

    inpath  = Path(args.input)
    outdir  = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Gene stem is the prefix before the first dot (e.g. SMEL5_01g008730.filtered
    # → SMEL5_01g008730), so the filename stays flat across stages.
    gene_stem = inpath.stem.split(".")[0]
    outpath = outdir / f"{gene_stem}.rescored.tsv"

    if not args.overwrite and outpath.exists():
        _log(f"[02_rescore] Skipping (overwrite=false): {outpath}")
        return

    # Read input
    with open(inpath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows   = list(reader)
        fields = list(reader.fieldnames) if reader.fieldnames else []

    if not rows:
        # Write an empty-but-headered output so downstream stages see a fresh
        # 0-row file instead of reading a stale *_rescored.tsv from a prior run.
        empty_fields = fields + ["crisprOn_score", "DeepSpCas9_score", "rescore_flag"]
        with open(outpath, "w", newline="") as fh:
            csv.DictWriter(fh, fieldnames=empty_fields, delimiter="\t").writeheader()
        _log(f"[02_rescore] Input is empty: {inpath} — wrote empty output: {outpath}")
        return

    seqs = [get_protospacer(r) for r in rows]

    # Load the target FASTA once up front so per-row DeepSpCas9 context
    # lookups are cheap. This is the genomic/CDS sequence CRISPOR received
    # as input; its flanks are what DeepSpCas9 actually needs.
    target_fasta_seq = ""
    if args.target_fasta:
        target_fasta_seq = _load_target_fasta(args.target_fasta)
        if not target_fasta_seq:
            _log(f"[02_rescore] --target-fasta supplied but could not be read: "
                 f"{args.target_fasta} — DeepSpCas9 will use N-padded fallback.",
                 level="WARN")
    elif "DeepSpCas9" in args.predictors:
        _log("[02_rescore] No --target-fasta provided; DeepSpCas9 scores will be "
             "computed from N-padded 23-nt targetSeq and will be systematically "
             "depressed. Pass the stage-[1] target FASTA to get calibrated scores.",
             level="WARN")

    # Per-row DeepSpCas9 context and source bookkeeping (so the thesis can
    # report the fraction of scores that used the real genomic flank vs the
    # N-padded fallback).
    ds_contexts: list[str] = []
    ds_sources:  list[str] = []
    if "DeepSpCas9" in args.predictors:
        for row in rows:
            ctx, src = build_deepspcas9_context(row, target_fasta_seq)
            ds_contexts.append(ctx)
            ds_sources.append(src)
        n_fasta    = sum(1 for s in ds_sources if s == "fasta")
        n_fallback = sum(1 for s in ds_sources if s == "n_padded")
        level = "INFO" if n_fallback == 0 else "WARN"
        _log(f"[02_rescore] DeepSpCas9 context: {n_fasta} from target FASTA, "
             f"{n_fallback} N-padded fallback (N-padded scores systematically low).",
             level=level)

    # Build 23-nt protospacer+PAM list once for the plant-trained scorers.
    # Prefer the CRISPOR-emitted targetSeq (23 nt = 20 guide + 3 PAM); fall
    # back to the protospacer padded with NGG if targetSeq is absent.
    seqs23 = []
    for row in rows:
        t23 = _get_targetseq(row)
        if len(t23) >= 23:
            seqs23.append(t23[:23])
        else:
            seqs23.append((get_protospacer(row) + "NGG")[:23])

    # Run predictors
    score_cols: dict[str, list[float]] = {}
    for predictor in args.predictors:
        if predictor == "crisprOn":
            score_cols["crisprOn_score"] = score_crisprOn(seqs, workers=args.workers)
        elif predictor == "DeepSpCas9":
            # DeepSpCas9 gets the full 30-nt context (fasta-sourced when possible).
            score_cols["DeepSpCas9_score"] = score_DeepSpCas9(ds_contexts, workers=args.workers)
        elif predictor == "DeepCRISPR":
            # Plant-usable rescorer (Chuai 2018). Dispatches via subprocess
            # to an isolated TF-1.x env when --deepcrispr-env is provided.
            score_cols["DeepCRISPR_score"] = score_DeepCRISPR(
                seqs23, workers=args.workers, conda_env=args.deepcrispr_env)
        elif predictor == "CRISPR-Local" or predictor == "CRISPR_Local":
            # Local rescorer (Sun 2019). Scoring engine is Doench Rule Set 2
            # — see score_CRISPR_Local docstring re. training provenance.
            score_cols["CRISPR_Local_score"] = score_CRISPR_Local(
                seqs23, workers=args.workers,
                conda_env=args.crispr_local_env)
        else:
            _log(f"[02_rescore] Unknown predictor '{predictor}' — skipping.", level="WARN")

    # Add flag: guides where ALL available scores are < threshold. An extra
    # DeepSpCas9_context_source column records whether each row's 30-nt
    # context came from the target FASTA or the N-padded fallback, so the
    # calibration status is auditable per guide.
    extra_cols = list(score_cols.keys()) + ["rescore_flag"]
    if "DeepSpCas9" in args.predictors:
        extra_cols.append("DeepSpCas9_context_source")
    # Alias the CRISPR-P v2.0 raw "Score" column to "crisprP_score" so the
    # downstream rank script (08_rank_guides.py) and comparison scatter
    # (09_comparison_scatter.py) recognise it as a 0-100 on-target scorer.
    # Without the alias the original CRISPR-P score is silently dropped from
    # the composite KO score, leaving every guide with c_ontarget=0 and a
    # flat composite of ~0.15, and the scatter falls back to c_ontarget × 100.
    alias_crisprp = "Score" in fields and "crisprP_score" not in fields
    if alias_crisprp:
        extra_cols.append("crisprP_score")
    new_fields = fields + extra_cols
    for i, row in enumerate(rows):
        best = None
        for col, scores in score_cols.items():
            sc = scores[i]
            if sc == sc:  # not NaN
                row[col] = f"{sc:.4f}"
                best = sc if (best is None or sc > best) else best
            else:
                row[col] = "NA"
        row["rescore_flag"] = "low_efficiency" if (best is not None and best < args.flag_threshold) else ""
        if "DeepSpCas9" in args.predictors:
            row["DeepSpCas9_context_source"] = ds_sources[i] if i < len(ds_sources) else "n_padded"
        if alias_crisprp:
            row["crisprP_score"] = row.get("Score", "")

    with open(outpath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=new_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    flagged = sum(1 for r in rows if r.get("rescore_flag") == "low_efficiency")
    _log(f"[02_rescore] {len(rows)} guides rescored; {flagged} flagged low_efficiency -> {outpath}")


if __name__ == "__main__":
    main()
