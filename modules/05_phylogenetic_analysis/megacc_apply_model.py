#!/usr/bin/env python3
"""Apply best-fit substitution model from a MEGA-CC model-selection CSV to an
inference .mao template.

Usage:
    megacc_apply_model.py --csv <model_select.csv> \\
                          --template <inference.mao> \\
                          --output <derived.mao> \\
                          [--threads N]

Reads the best-fit model from the top data row of the CSV (MEGA ranks by BIC
ascending, so row 1 = best), parses +G / +I / +F suffixes, and substitutes
"Model/Method" and "Rates among Sites" in the inference .mao. The output is a
new .mao ready for `megacc -a`.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

NUC_MODEL_MAP = {
    "GTR":  "General Time Reversible model",
    "TN":   "Tamura-Nei model",
    "TN93": "Tamura-Nei model",
    "HKY":  "Hasegawa-Kishino-Yano model",
    "T92":  "Tamura 3-parameter model",
    "K2":   "Kimura 2-parameter model",
    "K80":  "Kimura 2-parameter model",
    "JC":   "Jukes-Cantor model",
    "JC69": "Jukes-Cantor model",
}

PROT_MODEL_MAP = {
    "LG":      "LG model",
    "JTT":     "Jones-Taylor-Thornton (JTT) model",
    "WAG":     "WAG model",
    "Dayhoff": "Dayhoff model",
    "Poisson": "Poisson model",
    "mtREV":   "mtREV24 model",
    "mtREV24": "mtREV24 model",
    "cpREV":   "cpREV model",
    "rtREV":   "rtREV model",
}

PROT_MODEL_F_MAP = {
    "LG":      "LG with Freqs. (+F) model",
    "JTT":     "JTT with Freqs. (+F) model",
    "WAG":     "WAG with Freqs. (+F) model",
    "Dayhoff": "Dayhoff with Freqs. (+F) model",
}

RATES_MAP = {
    ("G", "I"): "Gamma Distributed with Invariant Sites (G+I)",
    ("G",):     "Gamma Distributed (G)",
    ("I",):     "Has Invariant Sites (I)",
    ():         "Uniform rates",
}


def parse_model_code(code: str):
    parts = re.split(r"\+", code.strip())
    base = parts[0]
    # Strip trailing digits so ModelTest-NG style "G4"/"G8" normalize to "G".
    suffixes = {re.sub(r"\d+$", "", p.strip().upper()) for p in parts[1:] if p.strip()}
    suffixes.discard("")
    return base, suffixes


def best_model_from_csv(csv_path: Path) -> str:
    """Return the model code from the top data row of the MEGA model-selection CSV."""
    text = csv_path.read_text(encoding="utf-8", errors="replace")
    sample = text[:2048]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
    except csv.Error:
        dialect = csv.excel
    reader = csv.reader(text.splitlines(), dialect)
    rows = [r for r in reader if any(cell.strip() for cell in r)]
    if not rows:
        raise ValueError(f"Empty CSV: {csv_path}")
    header = [h.strip().lower() for h in rows[0]]
    model_idx = next((i for i, h in enumerate(header) if h in ("model", "model name", "name")), 0)
    for row in rows[1:]:
        if len(row) > model_idx and row[model_idx].strip():
            return row[model_idx].strip()
    raise ValueError(f"No data rows with model name in: {csv_path}")


def map_model_to_mao_keys(code: str, alphabet: str):
    base, suffixes = parse_model_code(code)
    has_g, has_i, has_f = "G" in suffixes, "I" in suffixes, "F" in suffixes

    if alphabet == "nucleotide":
        model_method = NUC_MODEL_MAP.get(base)
    elif alphabet == "protein":
        model_method = (PROT_MODEL_F_MAP.get(base) if has_f else None) or PROT_MODEL_MAP.get(base)
    else:
        raise ValueError(f"Unknown alphabet: {alphabet}")
    if not model_method:
        raise ValueError(f"Unknown {alphabet} model '{base}' (from '{code}')")

    rates_key = tuple(s for s, present in (("G", has_g), ("I", has_i)) if present)
    return model_method, RATES_MAP[rates_key]


def detect_alphabet(template: Path) -> str:
    text = template.read_text(encoding="utf-8", errors="replace")
    if re.search(r"datatype\s*=\s*snNucleotide", text):
        return "nucleotide"
    if re.search(r"datatype\s*=\s*snProtein", text):
        return "protein"
    raise ValueError(f"Cannot detect alphabet from template: {template}")


def _replace_value(line: str, new_value: str) -> str:
    m = re.match(r"^([^=]*=)\s*(.*)$", line)
    if not m:
        return line
    return f"{m.group(1)} {new_value}"


def substitute_mao(template: Path, output: Path, model_method: str, rates: str, threads):
    lines = template.read_text(encoding="utf-8").splitlines()
    out = []
    seen = {"model": False, "rates": False, "threads": False}

    for line in lines:
        s = line.lstrip()
        if s.startswith("Model/Method"):
            out.append(_replace_value(line, model_method)); seen["model"] = True
        elif s.startswith("Rates among Sites"):
            out.append(_replace_value(line, rates)); seen["rates"] = True
        elif threads is not None and s.startswith("Number of Threads"):
            out.append(_replace_value(line, str(threads))); seen["threads"] = True
        else:
            out.append(line)

    needs_inject = (not seen["model"] or not seen["rates"]
                    or (threads is not None and not seen["threads"]))
    if needs_inject:
        injected = False
        merged = []
        for line in out:
            merged.append(line)
            if not injected and line.strip().startswith("[ AnalysisSettings ]"):
                if not seen["model"]:
                    merged.append(f"Model/Method                         = {model_method}")
                if not seen["rates"]:
                    merged.append(f"Rates among Sites                    = {rates}")
                if threads is not None and not seen["threads"]:
                    merged.append(f"Number of Threads                    = {threads}")
                injected = True
        out = merged

    output.write_text("\n".join(out) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--csv", type=Path, help="MEGA-CC model-selection CSV (best model = top data row)")
    src.add_argument("--model-code", type=str,
                     help="Direct model code (e.g. 'JTT+I+G4+F'); bypasses CSV parsing")
    ap.add_argument("--template", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=None)
    args = ap.parse_args()

    if args.csv is not None and not args.csv.is_file():
        sys.exit(f"CSV not found: {args.csv}")
    if not args.template.is_file():
        sys.exit(f"Template .mao not found: {args.template}")

    alphabet = detect_alphabet(args.template)
    code = args.model_code if args.model_code else best_model_from_csv(args.csv)
    model_method, rates = map_model_to_mao_keys(code, alphabet)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    substitute_mao(args.template, args.output, model_method, rates, args.threads)
    print(f"Selected model: {code}", file=sys.stderr)
    print(f"  Model/Method      = {model_method}", file=sys.stderr)
    print(f"  Rates among Sites = {rates}", file=sys.stderr)
    print(f"  Output .mao       = {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
