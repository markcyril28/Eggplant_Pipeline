#!/usr/bin/env python3
"""
DeepCRISPR sequence-only inference helper for v3.

This wrapper takes a plain text file of 23-nt protospacer+PAM sequences (one
per line), runs them through DeepCRISPR's sequence-only on-target regression
model (`trained_models/ontar_cnn_reg_seq`), and writes the predicted
activity score (one per line, in input order) to the output file.

v3's 02_rescore_ontarget.py subprocess-invokes this script under the
isolated `crispr_v3_deepcrispr` conda env (TF 1.3 / sonnet 1.9 / Python 3.6).
Do not rename or relocate this file — the dispatch expects it at
`modules/09_crispr_analysis/v2/tools/DeepCRISPR/deepcrispr_infer.py`.

Training provenance: DeepCRISPR was trained on HUMAN CRISPR screens
(HCT116, HEK293, HeLa, HL60). It is commonly used as a plant on-target
baseline in the absence of a plant-native deep model, but its scores
should be treated as relative ordering, not absolute plant editing
activity.

Usage:
    python deepcrispr_infer.py \\
        --input   protospacers.txt \\
        --output  scores.txt \\
        --model-dir trained_models/ontar_cnn_reg_seq \\
        --seq-only
"""

from __future__ import print_function

import argparse
import os
import sys


# sgRNA one-hot encoding matches DeepCRISPR's training convention:
#   A -> [1, 0, 0, 0]
#   C -> [0, 1, 0, 0]
#   G -> [0, 0, 1, 0]
#   T -> [0, 0, 0, 1]
#   N / other -> [0, 0, 0, 0]
_ONEHOT = {
    "A": [1, 0, 0, 0],
    "C": [0, 1, 0, 0],
    "G": [0, 0, 1, 0],
    "T": [0, 0, 0, 1],
}


def one_hot_23mers(sequences):
    """Return a numpy array of shape [N, 4, 1, 23] from 23-nt strings."""
    import numpy as np
    n = len(sequences)
    arr = np.zeros((n, 4, 1, 23), dtype=np.float32)
    for i, seq in enumerate(sequences):
        s = (seq or "").strip().upper().ljust(23, "N")[:23]
        for j, ch in enumerate(s):
            vec = _ONEHOT.get(ch, [0, 0, 0, 0])
            for k in range(4):
                arr[i, k, 0, j] = vec[k]
    return arr


def parse_args():
    p = argparse.ArgumentParser(
        description="DeepCRISPR sequence-only on-target inference (v3 helper)")
    p.add_argument("--input", required=True,
                    help="Text file — one 23-nt protospacer+PAM per line.")
    p.add_argument("--output", required=True,
                    help="Text file — one float score per line, input order.")
    p.add_argument("--model-dir", required=True,
                    help="Directory containing the DeepCRISPR on-target "
                         "sequence-only regression model (e.g. "
                         "trained_models/ontar_cnn_reg_seq).")
    p.add_argument("--seq-only", action="store_true",
                    help="Use the sequence-only model (recommended for plant "
                         "work when no epigenetic tracks are available).")
    return p.parse_args()


def main():
    args = parse_args()

    # Read input sequences
    with open(args.input) as fh:
        seqs = [line.strip() for line in fh if line.strip()]
    if not seqs:
        open(args.output, "w").close()
        return

    # Validate TF / sonnet are importable (these only install in the
    # isolated crispr_v3_deepcrispr env).
    try:
        import numpy as np        # noqa: F401
        import tensorflow as tf
        import deepcrispr as dc
    except ImportError as exc:
        sys.stderr.write(
            "ERROR: DeepCRISPR inference requires tensorflow 1.3 + sonnet 1.9 "
            "+ python 3.6. Import failed: {0}\n".format(exc))
        sys.stderr.write(
            "       Create the env with:  conda create -n crispr_v3_deepcrispr "
            "python=3.6 && conda activate crispr_v3_deepcrispr && "
            "pip install tensorflow==1.3.0 dm-sonnet==1.9\n")
        sys.exit(1)

    x = one_hot_23mers(seqs)

    # TF 1.x graph-mode inference. This matches run_examples.py's
    # "On-target Seq-only Regression Task" block line-for-line.
    sess = tf.InteractiveSession() if hasattr(tf, "InteractiveSession") \
        else tf.compat.v1.InteractiveSession()
    model_dir = args.model_dir
    if not os.path.isdir(model_dir):
        sys.stderr.write(
            "ERROR: model dir not found: {0}\n".format(model_dir))
        sys.exit(2)

    try:
        dcmodel = dc.DCModelOntar(sess, model_dir, is_reg=True,
                                    seq_feature_only=True)
    except Exception as exc:
        sys.stderr.write(
            "ERROR: DCModelOntar init failed: {0}\n".format(exc))
        sys.exit(3)

    # DeepCRISPR expects shape [N, 4, 1, 23] for sequence-only regression.
    preds = dcmodel.ontar_predict(x)

    # `ontar_predict` returns a numpy array of shape [N] or [N, 1].
    flat = [float(v) for v in (preds.flatten() if hasattr(preds, "flatten")
                                 else preds)]
    if len(flat) != len(seqs):
        sys.stderr.write(
            "WARN: output length ({0}) != input length ({1}); padding with "
            "NaN.\n".format(len(flat), len(seqs)))
        flat = flat + [float("nan")] * (len(seqs) - len(flat))
        flat = flat[:len(seqs)]

    with open(args.output, "w") as fh:
        for v in flat:
            fh.write("{0:.6f}\n".format(v))


if __name__ == "__main__":
    main()
