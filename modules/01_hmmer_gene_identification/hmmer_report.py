#!/usr/bin/env python3
"""
HMMER Identification Report Generator
======================================
Scans all *_hits_filtered.domtbl, *_nhmmer_hits.tbl, and *.clstr files under an
01_Identification directory tree, then produces:
  d_REPORT/
    ├── summary.csv                        consolidated hit table (with search_source)
    └── report.md                          Markdown report with tables

Usage (standalone):
    python3 hmmer_report.py /path/to/01_Identification --gene-group DMP

Usage (called by orchestrator):
    python3 hmmer_report.py "$IDENT_DIR" --gene-group "$GENE_GROUP" --evalue "$E_VALUE"
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from datetime import datetime


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _derive_genome(ident_dir, path):
    """Derive genome label from a result file path under ident_dir.

    Handles both layouts:
      {ident_dir}/{genome}/e-value_.../...          -> genome
      {ident_dir}/shared/{genome}/e-value_.../...   -> genome (not "shared")
    """
    rel = path.relative_to(ident_dir)
    if len(rel.parts) < 2:
        return "unknown"
    if rel.parts[0] == "shared" and len(rel.parts) > 2:
        return rel.parts[1]
    return rel.parts[0]


def _index_result_files(ident_dir):
    """Single-pass directory walk to collect result files by type.

    Avoids repeated rglob() calls which are slow on WSL/network filesystems.
    Returns dict with keys: domtbl, nhmmer_tbl, prot_hit_ids, nucl_hit_ids,
    clstr, cdhit_fa — each a sorted list of Path objects.
    """
    import os as _os
    index = {
        "domtbl": [], "nhmmer_tbl": [], "prot_hit_ids": [],
        "nucl_hit_ids": [], "clstr": [], "cdhit_fa": [],
    }
    for root, _dirs, files in _os.walk(ident_dir):
        root_path = Path(root)
        for fname in files:
            fp = root_path / fname
            if fname.endswith("_hits_filtered.domtbl"):
                index["domtbl"].append(fp)
            elif fname.endswith("_nhmmer_hits.tbl"):
                index["nhmmer_tbl"].append(fp)
            elif fname.endswith("_prot_hit_ids.txt"):
                index["prot_hit_ids"].append(fp)
            elif fname.endswith("_nucl_hit_ids.txt"):
                index["nucl_hit_ids"].append(fp)
            elif fname.endswith(".clstr"):
                index["clstr"].append(fp)
            elif fname.endswith("_cdhit.fa"):
                index["cdhit_fa"].append(fp)
    return {k: sorted(v) for k, v in index.items()}


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_domtbl_line(line):
    """Parse a single non-comment line from HMMER domtblout format."""
    parts = line.strip().split()
    if len(parts) < 22:
        return None
    desc = " ".join(parts[22:]) if len(parts) > 22 else ""
    gene_match = re.search(r"gene=([^\s]+)", desc)
    name_match = re.search(r'Name:"([^"]*)"', desc)
    return {
        "gene_id":          parts[0],
        "target_acc":       parts[1],
        "protein_length":   int(parts[2]),
        "domain_name":      parts[3],
        "pfam_id":          parts[4],
        "model_length":     int(parts[5]),
        "full_seq_evalue":  float(parts[6]),
        "full_seq_score":   float(parts[7]),
        "full_seq_bias":    float(parts[8]),
        "domain_num":       int(parts[9]),
        "total_domains":    int(parts[10]),
        "c_evalue":         float(parts[11]),
        "i_evalue":         float(parts[12]),
        "domain_score":     float(parts[13]),
        "domain_bias":      float(parts[14]),
        "hmm_from":         int(parts[15]),
        "hmm_to":           int(parts[16]),
        "ali_from":         int(parts[17]),
        "ali_to":           int(parts[18]),
        "env_from":         int(parts[19]),
        "env_to":           int(parts[20]),
        "accuracy":         float(parts[21]),
        "gene_name":        gene_match.group(1) if gene_match else "",
        "annotation":       name_match.group(1) if name_match else "",
    }


def scan_domtbl_files(ident_dir, file_list=None):
    """Walk the identification dir tree and collect all filtered domtbl hits.

    Returns list of dicts, each augmented with 'genome' key.
    If *file_list* is given, skip rglob and iterate over that instead.
    """
    records = []
    for domtbl in (file_list if file_list is not None
                   else sorted(ident_dir.rglob("*_hits_filtered.domtbl"))):
        genome = _derive_genome(ident_dir, domtbl)
        with open(domtbl) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                rec = parse_domtbl_line(line)
                if rec:
                    rec["genome"] = genome
                    records.append(rec)
    return records


def parse_nhmmer_tbl_line(line):
    """Parse a single non-comment line from nhmmer --tblout format.

    nhmmer tblout columns:
      0: target name, 1: accession, 2: query name, 3: query accession,
      4: hmmfrom, 5: hmmto, 6: alifrom, 7: alito, 8: envfrom, 9: envto,
      10: sq len, 11: strand, 12: E-value, 13: score, 14: bias,
      15+: description
    """
    parts = line.strip().split()
    if len(parts) < 15:
        return None
    try:
        return {
            "target_name":  parts[0],
            "query_name":   parts[2],
            "hmm_from":     int(parts[4]),
            "hmm_to":       int(parts[5]),
            "ali_from":     int(parts[6]),
            "ali_to":       int(parts[7]),
            "env_from":     int(parts[8]),
            "env_to":       int(parts[9]),
            "seq_length":   int(parts[10]),
            "strand":       parts[11],
            "evalue":       float(parts[12]),
            "score":        float(parts[13]),
            "bias":         float(parts[14]),
            "description":  " ".join(parts[15:]) if len(parts) > 15 else "",
        }
    except (ValueError, IndexError):
        return None


def scan_nhmmer_tbl_files(ident_dir, file_list=None):
    """Walk the identification dir tree for nhmmer tblout files.

    Returns list of dicts, each augmented with 'genome' and 'profile' keys.
    If *file_list* is given, skip rglob and iterate over that instead.
    """
    records = []
    for tbl in (file_list if file_list is not None
                else sorted(ident_dir.rglob("*_nhmmer_hits.tbl"))):
        genome = _derive_genome(ident_dir, tbl)
        profile = tbl.parent.name  # e.g. PF05078_DMP
        with open(tbl) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                rec = parse_nhmmer_tbl_line(line)
                if rec:
                    rec["genome"] = genome
                    rec["profile"] = profile
                    records.append(rec)
    return records


def scan_hit_source_files(ident_dir, prot_files=None, nucl_files=None):
    """Scan *_prot_hit_ids.txt and *_nucl_hit_ids.txt to compute hit-source breakdown.

    Returns dict: (genome, profile) -> {"prot_only": set, "nucl_only": set, "both": set}
    If *prot_files* / *nucl_files* are given, skip rglob for those file types.
    """
    # Collect protein hit IDs
    prot_ids = {}  # (genome, profile) -> set
    for f in (prot_files if prot_files is not None
              else sorted(ident_dir.rglob("*_prot_hit_ids.txt"))):
        genome = _derive_genome(ident_dir, f)
        profile = f.parent.name
        ids = {l.strip() for l in open(f) if l.strip()}
        prot_ids[(genome, profile)] = ids

    # Collect nucleotide hit IDs
    nucl_ids = {}  # (genome, profile) -> set
    for f in (nucl_files if nucl_files is not None
              else sorted(ident_dir.rglob("*_nucl_hit_ids.txt"))):
        genome = _derive_genome(ident_dir, f)
        profile = f.parent.name
        ids = {l.strip() for l in open(f) if l.strip()}
        nucl_ids[(genome, profile)] = ids

    all_keys = set(prot_ids) | set(nucl_ids)
    breakdown = {}
    for key in all_keys:
        p = prot_ids.get(key, set())
        n = nucl_ids.get(key, set())
        breakdown[key] = {
            "prot_only": p - n,
            "nucl_only": n - p,
            "both":      p & n,
        }
    return breakdown


def scan_clstr_files(ident_dir, clstr_files=None, cdhit_files=None):
    """Parse CD-HIT .clstr files to get before/after counts per profile × genome.

    Returns separate dicts for protein and transcript clustering.
    If *clstr_files* / *cdhit_files* are given, skip rglob for those file types.
    """
    prot_stats = {}   # (genome, profile) -> {"before": int, "after": int}
    trans_stats = {}  # (genome, profile) -> {"before": int, "after": int}

    for clstr in (clstr_files if clstr_files is not None
                  else sorted(ident_dir.rglob("*.clstr"))):
        genome = _derive_genome(ident_dir, clstr)
        profile = clstr.parent.name
        total_seqs = 0
        n_clusters = 0
        with open(clstr) as fh:
            for line in fh:
                if line.startswith(">Cluster"):
                    n_clusters += 1
                elif line.strip():
                    total_seqs += 1
        is_transcript = "_transcripts_cdhit" in clstr.name
        target = trans_stats if is_transcript else prot_stats
        target[(genome, profile)] = {"before": total_seqs, "after": n_clusters}

    # Definitive "after" counts from output FASTAs
    for cdhit_fa in (cdhit_files if cdhit_files is not None
                     else sorted(ident_dir.rglob("*_cdhit.fa"))):
        if cdhit_fa.suffix == ".clstr":
            continue
        genome = _derive_genome(ident_dir, cdhit_fa)
        profile = cdhit_fa.parent.name
        n_seqs = sum(1 for l in open(cdhit_fa) if l.startswith(">"))
        is_transcript = "_transcripts_cdhit" in cdhit_fa.name
        target = trans_stats if is_transcript else prot_stats
        key = (genome, profile)
        if key in target:
            target[key]["after"] = n_seqs
        else:
            target[key] = {"before": n_seqs, "after": n_seqs}

    return prot_stats, trans_stats


# ---------------------------------------------------------------------------
# Per-genome stats CSV  (one row per genome — designed for 30+ genome runs)
# ---------------------------------------------------------------------------

PGS_FIELDS = [
    "genome", "total_domain_hits", "unique_genes", "unique_domains",
    "best_evalue", "best_score",
    "prot_hits_before_cdhit", "prot_hits_after_cdhit",
    "nhmmer_hits", "nhmmer_unique_targets",
    "source_prot_only", "source_nucl_only", "source_both",
]


def write_per_genome_stats_csv(records, nhmmer_records, prot_stats, trans_stats,
                                source_breakdown, out_path):
    """Write per_genome_stats.csv — one aggregate row per genome."""
    all_genomes = sorted(
        {r["genome"] for r in records} | {r["genome"] for r in nhmmer_records}
    )
    rows = []
    for g in all_genomes:
        g_recs       = [r for r in records       if r["genome"] == g]
        g_nhmmer     = [r for r in nhmmer_records if r["genome"] == g]
        unique_genes  = len({r["gene_id"]     for r in g_recs})
        unique_doms   = len({r["domain_name"] for r in g_recs})
        total_hits    = len(g_recs)
        best_ev    = min((r["full_seq_evalue"] for r in g_recs), default=None)
        best_score = max((r["full_seq_score"]  for r in g_recs), default=None)
        nhmmer_hits = len(g_nhmmer)
        nhmmer_tgts = len({r["target_name"] for r in g_nhmmer})
        prot_before = sum(s["before"] for (gg, _), s in prot_stats.items()      if gg == g)
        prot_after  = sum(s["after"]  for (gg, _), s in prot_stats.items()      if gg == g)
        prot_only   = sum(len(v["prot_only"]) for (gg, _), v in source_breakdown.items() if gg == g)
        nucl_only   = sum(len(v["nucl_only"]) for (gg, _), v in source_breakdown.items() if gg == g)
        both        = sum(len(v["both"])      for (gg, _), v in source_breakdown.items() if gg == g)
        rows.append({
            "genome":                 g,
            "total_domain_hits":      total_hits,
            "unique_genes":           unique_genes,
            "unique_domains":         unique_doms,
            "best_evalue":            f"{best_ev:.2e}"    if best_ev    is not None else "",
            "best_score":             f"{best_score:.1f}" if best_score is not None else "",
            "prot_hits_before_cdhit": prot_before,
            "prot_hits_after_cdhit":  prot_after,
            "nhmmer_hits":            nhmmer_hits,
            "nhmmer_unique_targets":  nhmmer_tgts,
            "source_prot_only":       prot_only,
            "source_nucl_only":       nucl_only,
            "source_both":            both,
        })
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=PGS_FIELDS)
        w.writeheader()
        w.writerows(rows)
    return rows


# ---------------------------------------------------------------------------
# Per-genome gene list CSVs  (keeps report.md readable for 30+ genomes)
# ---------------------------------------------------------------------------

def write_gene_list_csvs(records, gene_list_dir):
    """Write one CSV per genome under gene_list_dir.

    Returns list of (genome_label, csv_path, unique_gene_count).
    """
    gene_list_dir.mkdir(parents=True, exist_ok=True)
    fields = [
        "gene_id", "domain_name", "full_seq_evalue", "full_seq_score",
        "protein_length", "ali_from", "ali_to", "accuracy",
        "gene_name", "annotation", "search_source",
    ]
    written = []
    for g in sorted({r["genome"] for r in records}):
        g_recs = [r for r in records if r["genome"] == g]
        # Deduplicate by gene_id — keep highest-scoring domain hit per gene
        seen: dict = {}
        for r in sorted(g_recs, key=lambda x: x["full_seq_score"], reverse=True):
            if r["gene_id"] not in seen:
                seen[r["gene_id"]] = r
        safe_name = re.sub(r"[^\w.-]", "_", g)
        out_path = gene_list_dir / f"{safe_name}_genes.csv"
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(seen.values())
        _write_tsv_copy(out_path)
        written.append((g, out_path, len(seen)))
    return written


# ---------------------------------------------------------------------------
# Summary CSV
# ---------------------------------------------------------------------------

CSV_FIELDS = [
    "genome", "gene_id", "search_source", "gene_name", "annotation",
    "protein_length", "domain_name", "pfam_id", "model_length",
    "full_seq_evalue", "full_seq_score", "full_seq_bias",
    "domain_num", "total_domains",
    "c_evalue", "i_evalue", "domain_score", "domain_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    "accuracy",
]


def write_summary_csv(records, out_path):
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(records, key=lambda r: (r["genome"], r["gene_id"], r["full_seq_evalue"])))



# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def _md_table(headers, rows):
    """Build a simple Markdown table string."""
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines)


def _md_to_html_simple(md_text):
    """Minimal Markdown→HTML converter (no external deps).

    Handles: ATX headings, GFM pipe tables, horizontal rules, bold, italic,
    inline code, <br> pass-through, and paragraph wrapping.
    """
    import html as html_mod

    def _inline(s):
        s = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', s)
        s = re.sub(r'\*(.+?)\*',     r'<em>\1</em>',         s)
        s = re.sub(r'`(.+?)`',       r'<code>\1</code>',      s)
        return s

    lines = md_text.splitlines()
    out = []
    i = 0
    while i < len(lines):
        line = lines[i]

        # Pass-through for HTML block tags (e.g. our landscape <div>) and
        # for HTML comments — do not wrap in <p>.
        stripped = line.strip()
        if stripped.startswith('<!--') and stripped.endswith('-->'):
            out.append(stripped)
            i += 1; continue
        if re.match(r'^\s*</?(div|section|article)\b', line):
            out.append(stripped)
            i += 1; continue

        # Heading
        m = re.match(r'^(#{1,6})\s+(.*)', line)
        if m:
            lvl = len(m.group(1))
            out.append(f"<h{lvl}>{_inline(m.group(2))}</h{lvl}>")
            i += 1; continue

        # Horizontal rule
        if re.match(r'^---+$', line.strip()):
            out.append("<hr>")
            i += 1; continue

        # GFM pipe table: collect consecutive | lines
        if line.strip().startswith('|'):
            table_lines = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                table_lines.append(lines[i])
                i += 1
            # Identify separator row (all dashes/colons)
            sep_idx = None
            for ti, tl in enumerate(table_lines):
                if re.match(r'^[\|\-\: ]+$', tl):
                    sep_idx = ti; break
            out.append('<table>')
            for ti, tl in enumerate(table_lines):
                if sep_idx is not None and ti == sep_idx:
                    continue
                cells = [c.strip() for c in tl.strip().strip('|').split('|')]
                tag = 'th' if (sep_idx is not None and ti < sep_idx) else 'td'
                row = ''.join(f'<{tag}>{_inline(c)}</{tag}>' for c in cells)
                out.append(f'<tr>{row}</tr>')
            out.append('</table>')
            continue

        # Blank line
        if not line.strip():
            out.append('')
            i += 1; continue

        # Paragraph: group consecutive non-special lines into one <p> so that
        # embedded <br> tags (as in the metadata block) produce visible breaks.
        para_lines = []
        while i < len(lines):
            l = lines[i]
            if (not l.strip()
                    or re.match(r'^(#{1,6})\s', l)
                    or re.match(r'^---+$', l.strip())
                    or l.strip().startswith('|')):
                break
            para_lines.append(_inline(l))
            i += 1
        if para_lines:
            out.append('<p>' + '\n'.join(para_lines) + '</p>')

    return '\n'.join(out)


def _write_report_pdf(md_path):
    """Render a Markdown report to PDF alongside the .md file.

    Tries (in order):
      1. pandoc --pdf-engine=xelatex/pdflatex/lualatex
      2. weasyprint + markdown  (pure-Python, needs GTK on Windows)
      3. Chrome / Edge headless --print-to-pdf  (same HTML+CSS as weasyprint)
      4. reportlab               (pure-Python, no system libs required)

    Returns the PDF Path on success, None if no tool is available.
    """
    import subprocess
    import shutil
    import tempfile

    md_path = Path(md_path)
    pdf_path = md_path.with_suffix(".pdf")

    # Shared HTML+CSS — used by both weasyprint and Chrome headless so the
    # output style is identical regardless of which tool succeeds.
    # The @page rules give the ortholog section its own landscape page.
    CSS = (
        "@page{size:A4 portrait;margin:1.5cm}"
        "@page landscape{size:A4 landscape;margin:1.2cm}"
        "body{font-family:DejaVu Sans,Arial,sans-serif;max-width:960px;"
        "margin:auto;font-size:10pt}"
        "p{margin:0.4em 0;line-height:1.5}"
        "table{border-collapse:collapse;width:100%;margin:1em 0}"
        "th,td{border:1px solid #bbb;padding:3px 7px;font-size:9pt}"
        "th{background:#f0f0f0}"
        "code{background:#f5f5f5;padding:1px 4px;font-family:monospace}"
        "pre{background:#f5f5f5;padding:8px;overflow-x:auto}"
        "h1{font-size:16pt}h2{font-size:13pt}h3{font-size:11pt}"
        f".{_LANDSCAPE_CSS_CLASS}{{page:landscape;"
        "page-break-before:always;page-break-after:always;"
        "break-before:page;break-after:page;max-width:none}}"
        f".{_LANDSCAPE_CSS_CLASS} table{{font-size:8pt}}"
        f".{_LANDSCAPE_CSS_CLASS} th,.{_LANDSCAPE_CSS_CLASS} td{{font-size:7pt;padding:2px 4px}}"
    )

    def _build_html(md_text):
        md_text = _md_with_landscape_html(md_text)
        try:
            import markdown as md_lib          # type: ignore
            body = md_lib.markdown(md_text, extensions=["tables", "fenced_code"])
        except ImportError:
            body = _md_to_html_simple(md_text)
        return (
            f"<html><head><meta charset='utf-8'>"
            f"<style>{CSS}</style></head><body>{body}</body></html>"
        )

    # --- 1. pandoc ---
    pandoc = shutil.which("pandoc")
    if pandoc:
        # Preprocess Markdown: replace landscape markers with raw LaTeX so
        # pandoc emits \begin{landscape}...\end{landscape} via the pdflscape
        # package. We switch --from=markdown (not gfm) to enable raw_tex.
        md_for_latex = _md_with_landscape_latex(md_path.read_text(encoding="utf-8"))
        tmp_md = md_path.with_name(md_path.stem + ".__landscape_tmp__.md")
        tmp_md.write_text(md_for_latex, encoding="utf-8")
        # pdflscape alone fails to emit /Rotate 90 on xelatex+xdvipdfmx
        # (the PDF backend used by xelatex on many distros). Adding a
        # shipout hook that issues \special{pdf:put @thispage <</Rotate 90>>}
        # for every page inside a landscape environment forces the
        # rotation on every affected page, regardless of driver quirks.
        header_tex = (
            r"\usepackage{pdflscape}"  "\n"
            r"\usepackage{etoolbox}"   "\n"
            r"\usepackage{atbegshi}"   "\n"
            r"\makeatletter"           "\n"
            r"\newif\ifLS@inlandscape" "\n"
            r"\AtBeginShipout{\ifLS@inlandscape"
            r"\special{pdf:put @thispage <</Rotate 90>>}\fi}" "\n"
            r"\AtBeginEnvironment{landscape}{\LS@inlandscapetrue}" "\n"
            r"\AfterEndEnvironment{landscape}{\LS@inlandscapefalse}" "\n"
            r"\makeatother"            "\n"
        )
        tmp_hdr = md_path.with_name(md_path.stem + ".__landscape_hdr__.tex")
        tmp_hdr.write_text(header_tex, encoding="utf-8")
        try:
            for engine in ("xelatex", "pdflatex", "lualatex"):
                try:
                    r = subprocess.run(
                        [pandoc, str(tmp_md), "-o", str(pdf_path),
                         "--from=markdown+raw_tex+pipe_tables",
                         f"--pdf-engine={engine}",
                         "-V", "geometry:margin=1in", "-V", "fontsize=11pt",
                         f"--include-in-header={tmp_hdr}",
                         "--standalone"],
                        capture_output=True, text=True, timeout=120,
                    )
                    if r.returncode == 0 and pdf_path.exists():
                        return pdf_path
                except (subprocess.TimeoutExpired, FileNotFoundError):
                    continue
        finally:
            for p in (tmp_md, tmp_hdr):
                try:
                    p.unlink()
                except OSError:
                    pass

    # --- 2. weasyprint ---
    try:
        import weasyprint                  # type: ignore
        html = _build_html(md_path.read_text(encoding="utf-8"))
        weasyprint.HTML(string=html, base_url=str(md_path.parent)).write_pdf(str(pdf_path))
        if pdf_path.exists():
            return pdf_path
    except (ImportError, Exception):
        pass

    # --- 3. Chrome / Edge headless (same HTML+CSS → identical visual style) ---
    try:
        html = _build_html(md_path.read_text(encoding="utf-8"))
        html_tmp = Path(tempfile.mktemp(suffix=".html"))
        html_tmp.write_text(html, encoding="utf-8")
        browser = None
        for candidate in [
            shutil.which("google-chrome"),
            shutil.which("google-chrome-stable"),
            shutil.which("chromium"),
            shutil.which("chromium-browser"),
            r"C:\Program Files\Google\Chrome\Application\chrome.exe",
            r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        ]:
            if candidate and Path(candidate).exists():
                browser = candidate
                break
        if browser:
            r = subprocess.run(
                [browser,
                 "--headless", "--disable-gpu", "--no-sandbox",
                 f"--print-to-pdf={pdf_path}",
                 "--print-to-pdf-no-header",
                 html_tmp.as_uri()],
                capture_output=True, timeout=60,
            )
            html_tmp.unlink(missing_ok=True)
            if pdf_path.exists():
                return pdf_path
        html_tmp.unlink(missing_ok=True)
    except (ImportError, Exception):
        pass

    # --- 4. reportlab (pure-Python, no system libs required) ---
    try:
        import re as _re
        import markdown as md_lib          # type: ignore  # noqa: F401
        from reportlab.lib.pagesizes import A4, landscape
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import cm
        from reportlab.lib import colors
        from reportlab.platypus import (
            BaseDocTemplate, PageTemplate, Frame, NextPageTemplate, PageBreak,
            Paragraph, Spacer, Table, TableStyle, HRFlowable,
        )
        text = md_path.read_text(encoding="utf-8")
        s = getSampleStyleSheet()
        h1s  = ParagraphStyle("rH1", parent=s["Heading1"], fontSize=14, spaceAfter=8)
        h2s  = ParagraphStyle("rH2", parent=s["Heading2"], fontSize=11,
                              spaceAfter=6, spaceBefore=12)
        norm = ParagraphStyle("rNorm", parent=s["Normal"], fontSize=8, leading=11)

        # Two page templates: 'portrait' for the bulk of the report,
        # 'landscape' for the ortholog table section.
        portrait_size = A4
        landscape_size = landscape(A4)
        margin = 1.5 * cm
        portrait_frame = Frame(margin, margin,
                               portrait_size[0] - 2*margin,
                               portrait_size[1] - 2*margin,
                               id='p_frame')
        landscape_frame = Frame(margin, margin,
                                landscape_size[0] - 2*margin,
                                landscape_size[1] - 2*margin,
                                id='l_frame')
        doc = BaseDocTemplate(
            str(pdf_path),
            pagesize=portrait_size,
            leftMargin=margin, rightMargin=margin,
            topMargin=margin, bottomMargin=margin,
        )
        doc.addPageTemplates([
            PageTemplate(id='portrait', frames=[portrait_frame], pagesize=portrait_size),
            PageTemplate(id='landscape', frames=[landscape_frame], pagesize=landscape_size),
        ])

        def _table_for(rows, page_width):
            cw = (page_width - 3*cm) / max(len(rows[0]), 1)
            tbl = Table(rows, colWidths=[cw]*len(rows[0]), repeatRows=1)
            tbl.setStyle(TableStyle([
                ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#f0f0f0")),
                ("FONTSIZE",(0,0),(-1,-1),7),
                ("GRID",(0,0),(-1,-1),0.3,colors.HexColor("#bbbbbb")),
                ("VALIGN",(0,0),(-1,-1),"TOP"),
                ("LEFTPADDING",(0,0),(-1,-1),3),
                ("RIGHTPADDING",(0,0),(-1,-1),3),
                ("TOPPADDING",(0,0),(-1,-1),2),
                ("BOTTOMPADDING",(0,0),(-1,-1),2),
            ]))
            return tbl

        story, lines, i = [], text.split("\n"), 0
        # Track whether we are currently rendering inside the landscape block
        in_landscape = False
        current_page_w = portrait_size[0]
        while i < len(lines):
            ln = lines[i]
            stripped = ln.strip()

            # Landscape region markers: switch page template + force a page break.
            if stripped == _LANDSCAPE_START_MARKER:
                story.append(NextPageTemplate('landscape'))
                story.append(PageBreak())
                in_landscape = True
                current_page_w = landscape_size[0]
                i += 1
                continue
            if stripped == _LANDSCAPE_END_MARKER:
                story.append(NextPageTemplate('portrait'))
                story.append(PageBreak())
                in_landscape = False
                current_page_w = portrait_size[0]
                i += 1
                continue

            if ln.startswith("# ") and not ln.startswith("## "):
                story.append(Paragraph(ln[2:].strip(), h1s))
            elif ln.startswith("## "):
                story.append(Spacer(1, 6))
                story.append(Paragraph(ln[3:].strip(), h2s))
            elif stripped == "---":
                story.append(HRFlowable(width="100%", thickness=0.5, color=colors.grey))
            elif ln.startswith("| "):
                rows = []
                while i < len(lines) and lines[i].startswith("|"):
                    if "---" not in lines[i]:
                        cells = [c.strip().replace("**","")
                                 for c in lines[i].strip("|").split("|")]
                        rows.append(cells)
                    i += 1
                if rows:
                    story.append(_table_for(rows, current_page_w))
                    story.append(Spacer(1, 4))
                continue
            elif stripped.startswith(("- ","* ")):
                story.append(Paragraph("• " + stripped[2:].replace("**",""), norm))
            elif stripped:
                story.append(Paragraph(
                    _re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", ln), norm))
            i += 1
        doc.build(story)
        if pdf_path.exists():
            return pdf_path
    except (ImportError, Exception):
        pass

    return None


def _write_tsv_copy(csv_path):
    """Write a TSV copy alongside a CSV file (same stem, .tsv extension).

    Uses csv.reader/writer so embedded commas/quotes are handled correctly
    rather than a naive comma→tab substitution.
    """
    csv_path = Path(csv_path)
    tsv_path = csv_path.with_suffix(".tsv")
    with open(csv_path, newline="", encoding="utf-8") as src, \
         open(tsv_path, "w", newline="", encoding="utf-8") as dst:
        reader = csv.reader(src)
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        for row in reader:
            writer.writerow(row)
    return tsv_path


# ---------------------------------------------------------------------------
# GRF / GIF family classification
# ---------------------------------------------------------------------------

FAMILY_FIELDS = ["genome", "GRF", "QLQ", "WRC", "GIF", "SSXT"]


# ---------------------------------------------------------------------------
# Ortholog/Gene Family Mapping Table (BLAST-based reciprocal best hits)
# ---------------------------------------------------------------------------

# BLAST reference / columns: v4.1 (so all 8 v4.1 SmelDMPs are column positions).
ORTHO_REF_GENOME = "Solanum_melongena_v4.1"
# Row rendering order: V3 as second row (first data row), v4.1 as third row,
# remaining genomes alphabetical afterwards. Decoupled from BLAST ref so we
# can pin row positions without changing the column reference genome.
ORTHO_PINNED_ROWS = ("Solanum_melongena_V3", "Solanum_melongena_v4.1")
# Kept for back-compat with memory/docs that reference it.
ORTHO_SECOND_GENOME = "Solanum_melongena_V3"
# Safe field separator for combined BLAST FASTA IDs (must not occur in genome
# labels or gene IDs). "|" is reserved by some BLAST versions; use "___".
_ORTHO_SEP = "___"

# Invisible Markdown markers used to denote the region that should render as
# a dedicated landscape-oriented page in the PDF. In report.md these appear
# as HTML comments (harmless in MD viewers); PDF backends translate them
# into LaTeX pdflscape blocks (pandoc), @page CSS (weasyprint/chrome), or a
# PageBreak + landscape template switch (reportlab).
_LANDSCAPE_START_MARKER = "<!-- PAGE:LANDSCAPE:START -->"
_LANDSCAPE_END_MARKER = "<!-- PAGE:LANDSCAPE:END -->"
_LANDSCAPE_CSS_CLASS = "landscape-page"


def _md_with_landscape_html(md_text):
    """Replace landscape markers with an HTML div wrapper for HTML backends."""
    return (
        md_text
        .replace(_LANDSCAPE_START_MARKER,
                 f'<div class="{_LANDSCAPE_CSS_CLASS}">')
        .replace(_LANDSCAPE_END_MARKER, '</div>')
    )


def _md_with_landscape_latex(md_text):
    """Replace landscape markers with raw LaTeX blocks for pandoc."""
    return (
        md_text
        .replace(_LANDSCAPE_START_MARKER,
                 '\n\\clearpage\n\\begin{landscape}\n')
        .replace(_LANDSCAPE_END_MARKER,
                 '\n\\end{landscape}\n\\clearpage\n')
    )


def _collect_gene_group_fastas(ident_dir, gene_group):
    """Find per-genome transcript FASTA files for *gene_group*.

    Uses the pre-CD-HIT transcript file ``{gene_group}_transcripts.fa`` so
    near-identical paralogs (e.g. the two V3 SmelDMPs that collapse into a
    single CD-HIT cluster) remain distinct reference genes. Duplicates
    between row genomes are resolved afterwards by gene-name similarity
    so each row gene still appears in at most one column.

    Returns dict: genome_label -> Path.
    """
    import os as _os

    target = f"{gene_group}_transcripts.fa"
    fastas = {}
    for root, _dirs, files in _os.walk(ident_dir):
        if target not in files:
            continue
        root_path = Path(root)
        fpath = root_path / target
        try:
            if fpath.stat().st_size == 0:
                continue
        except OSError:
            continue
        genome = _derive_genome(ident_dir, root_path)
        fastas[genome] = fpath
    return fastas


def _write_combined_blast_input(genome_fastas, combined_path):
    """Concatenate per-genome FASTAs with ``{genome}{SEP}{gene_id}`` IDs.

    Returns dict: combined_id -> (genome, original_gene_id).
    """
    id_map = {}
    with open(combined_path, "w", encoding="utf-8") as out:
        for genome, fasta in genome_fastas.items():
            with open(fasta, encoding="utf-8") as f:
                for line in f:
                    if line.startswith(">"):
                        orig = line[1:].strip().split()[0]
                        combined_id = f"{genome}{_ORTHO_SEP}{orig}"
                        id_map[combined_id] = (genome, orig)
                        out.write(f">{combined_id}\n")
                    else:
                        out.write(line)
    return id_map


def _run_cross_genome_blastn(combined_fa, work_dir, threads):
    """Build a BLAST DB from *combined_fa* and run blastn all-vs-all.

    Returns Path to the outfmt-6 results file, or None if BLAST is
    unavailable or any step fails.
    """
    import subprocess
    import shutil as _shutil

    if not _shutil.which("makeblastdb") or not _shutil.which("blastn"):
        return None

    db_prefix = work_dir / "ortho_db"
    try:
        subprocess.run(
            ["makeblastdb", "-in", str(combined_fa),
             "-dbtype", "nucl", "-out", str(db_prefix)],
            check=True, capture_output=True, text=True, timeout=600,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
        return None

    results = work_dir / "ortho_blastn.tsv"
    try:
        subprocess.run(
            ["blastn", "-query", str(combined_fa), "-db", str(db_prefix),
             "-outfmt", "6 qseqid sseqid pident evalue bitscore length",
             "-out", str(results),
             "-evalue", "1e-5",
             "-num_threads", str(max(1, int(threads))),
             "-max_target_seqs", "1000"],
            check=True, capture_output=True, text=True, timeout=1800,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
        return None
    return results


def _parse_all_hits(blast_tsv, id_map):
    """Parse BLAST outfmt-6 into all cross-genome hits (not just the top one).

    Keeping every hit (instead of only the single best per subject genome)
    lets the bipartite matcher find secondary matches for near-identical
    V3 paralogs: when two V3 genes both have the same top hit in v4.1,
    one gets the top and the other can fall through to its next-best hit.

    Returns dict: (q_genome, q_gene) -> {s_genome: [(s_gene, pident, evalue, bitscore), ...]}
    sorted by bit-score descending.
    """
    from collections import defaultdict as _dd
    all_hits = _dd(lambda: _dd(list))
    with open(blast_tsv, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            q_comb, s_comb, pident, evalue, bitscore = parts[:5]
            if q_comb == s_comb:
                continue
            q_info = id_map.get(q_comb)
            s_info = id_map.get(s_comb)
            if q_info is None or s_info is None:
                continue
            q_gen, q_id = q_info
            s_gen, s_id = s_info
            if q_gen == s_gen:  # skip within-genome hits
                continue
            try:
                bits = float(bitscore)
                pid = float(pident)
                ev = float(evalue)
            except ValueError:
                continue
            all_hits[(q_gen, q_id)][s_gen].append((s_id, pid, ev, bits))

    for qk in all_hits:
        for sg in all_hits[qk]:
            all_hits[qk][sg].sort(key=lambda h: h[3], reverse=True)
    return all_hits


def build_ortholog_matching(all_hits, ref_genome):
    """Build 1:1 ortholog assignment between ref_genome and each other genome.

    Uses greedy bipartite matching on BLAST bit-score across all hits (not
    just the top), with reciprocal-best-hit edges prioritised over one-way
    best matches. This guarantees:
      * each row genome's gene appears in at most one V3 column, and
      * near-identical V3 paralogs (e.g. SMEL_010g339680 / SMEL_010g339690)
        each receive a distinct ortholog when the target genome has enough
        genes to support it, instead of one being dropped by strict RBH.

    Returns dict: ref_gene -> {genome -> ortholog_gene}.
    """
    orthologs = {}
    ref_genes = sorted({gene for (gen, gene) in all_hits.keys() if gen == ref_genome})
    other_genomes = sorted({g for (g, _) in all_hits.keys() if g != ref_genome})

    for rg in ref_genes:
        orthologs[rg] = {ref_genome: rg}

    for og_genome in other_genomes:
        # Collect every candidate edge (rg, og) from both directions, keeping
        # the highest bit-score seen. Both directions are needed because a
        # true ortholog pair may have one side's top hit "stolen" by a paralog
        # while the reverse side still points correctly.
        edges = {}  # (rg, og) -> max bit-score

        # Forward: all rg -> og_genome hits
        for rg in ref_genes:
            for (og, _pid, _ev, bits) in all_hits.get((ref_genome, rg), {}).get(og_genome, []):
                key = (rg, og)
                edges[key] = max(edges.get(key, 0.0), bits)

        # Reverse: all og -> ref_genome hits
        for (gen, og) in [k for k in all_hits.keys() if k[0] == og_genome]:
            for (rg, _pid, _ev, bits) in all_hits[(gen, og)].get(ref_genome, []):
                if rg not in orthologs:
                    continue
                key = (rg, og)
                edges[key] = max(edges.get(key, 0.0), bits)

        # Flag reciprocal-best edges: rg's top hit in og_genome == og AND
        # og's top hit in ref_genome == rg.
        def _is_rbh(rg, og):
            fwd = all_hits.get((ref_genome, rg), {}).get(og_genome, [])
            rev = all_hits.get((og_genome, og), {}).get(ref_genome, [])
            return bool(fwd) and bool(rev) and fwd[0][0] == og and rev[0][0] == rg

        edge_list = [(rg, og, bits, _is_rbh(rg, og)) for (rg, og), bits in edges.items()]
        # RBH edges first, then by bit-score descending; ref_gene + other_gene
        # break remaining ties deterministically.
        edge_list.sort(key=lambda e: (1 if e[3] else 0, e[2], e[0], e[1]), reverse=True)

        claimed_ref, claimed_other = set(), set()
        for rg, og, _bits, _rbh in edge_list:
            if rg in claimed_ref or og in claimed_other:
                continue
            orthologs[rg][og_genome] = og
            claimed_ref.add(rg)
            claimed_other.add(og)

    return orthologs


def build_ortholog_mapping(ident_dir, gene_group, threads=4):
    """Build cross-genome ortholog mapping using BLASTn + reciprocal best hits.

    Returns dict: ref_gene -> {genome -> ortholog_gene}. Empty dict when
    BLAST is unavailable, when the reference genome is missing, or when no
    gene-group FASTAs are found. Prints status lines on stdout.
    """
    import tempfile
    import shutil as _shutil

    fasta_name = f"{gene_group}_transcripts.fa"
    genome_fastas = _collect_gene_group_fastas(ident_dir, gene_group)
    if not genome_fastas:
        print(f"  ! Ortholog table skipped: no {fasta_name} files found under {ident_dir}")
        return {}
    if ORTHO_REF_GENOME not in genome_fastas:
        print(f"  ! Ortholog table skipped: reference genome '{ORTHO_REF_GENOME}' "
              f"has no {fasta_name}")
        return {}
    if ORTHO_SECOND_GENOME not in genome_fastas:
        print(f"  ! Ortholog table warning: expected second row genome "
              f"'{ORTHO_SECOND_GENOME}' has no {fasta_name}")

    work_dir = Path(tempfile.mkdtemp(prefix=f"{gene_group}_ortho_"))
    try:
        combined = work_dir / f"{gene_group}_combined.fa"
        id_map = _write_combined_blast_input(genome_fastas, combined)
        if not id_map:
            print("  ! Ortholog table skipped: combined FASTA is empty")
            return {}

        blast_tsv = _run_cross_genome_blastn(combined, work_dir, threads)
        if blast_tsv is None:
            print("  ! Ortholog table skipped: blastn/makeblastdb not available or BLAST run failed")
            return {}

        all_hits = _parse_all_hits(blast_tsv, id_map)
        orthologs = build_ortholog_matching(all_hits, ORTHO_REF_GENOME)
        n_other = len(genome_fastas) - 1
        print(f"  ✓ Ortholog table: {len(orthologs)} {ORTHO_REF_GENOME} reference gene(s), "
              f"1:1 bipartite matching across {n_other} other genome(s)")
        return orthologs
    finally:
        _shutil.rmtree(work_dir, ignore_errors=True)


def _order_ortholog_genomes(all_genomes):
    """Return genomes in report order per ORTHO_PINNED_ROWS, then alphabetical.

    Currently: V3 is the first data row, v4.1 is the second data row, every
    other genome follows alphabetically. Note this is independent of which
    genome was used as the BLAST reference (ORTHO_REF_GENOME).
    """
    remaining = set(all_genomes)
    ordered = []
    for pinned in ORTHO_PINNED_ROWS:
        if pinned in remaining:
            ordered.append(pinned)
            remaining.discard(pinned)
    ordered.extend(sorted(remaining))
    return ordered


def generate_ortholog_table(orthologs):
    """Generate ortholog table rows with V3 first, v4.1 second, then alphabetical.

    Columns are numbered 1..N for compactness — the first data row (V3)
    carries the full reference gene IDs, so readers can still map each
    column back to its reference gene by looking at that first row.

    Returns (headers, rows) suitable for _md_table.
    """
    ref_genes = sorted(orthologs.keys())

    all_genomes = set()
    for ortho in orthologs.values():
        all_genomes.update(ortho.keys())
    ordered_genomes = _order_ortholog_genomes(all_genomes)

    rows = []
    for genome in ordered_genomes:
        row = [genome]
        for ref_gene in ref_genes:
            row.append(orthologs.get(ref_gene, {}).get(genome, "-"))
        rows.append(row)

    headers = ["Genome"] + [str(i + 1) for i in range(len(ref_genes))]
    return headers, rows


def write_ortholog_table_csv(headers, rows, out_path):
    """Write the ortholog table to CSV (companion TSV produced by caller)."""
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(headers)
        w.writerows(rows)
    return out_path


def compute_family_classification(records):
    """Classify genes per genome into GRF (QLQ+WRC) and GIF (SSXT) families.

    Column definitions:
      - GRF:  unique genes containing BOTH QLQ and WRC domains (canonical GRF)
      - QLQ:  unique genes containing a QLQ domain hit
      - WRC:  unique genes containing a WRC domain hit
      - GIF:  unique genes containing an SSXT domain hit (canonical GIF)
      - SSXT: unique genes containing an SSXT domain hit (= GIF)

    Returns list of dicts keyed by FAMILY_FIELDS, sorted by (GRF+GIF) desc.
    """
    gene_doms = defaultdict(set)  # (genome, gene_id) -> {domain_name, ...}
    for r in records:
        gene_doms[(r["genome"], r["gene_id"])].add(r["domain_name"])

    per_genome = defaultdict(lambda: {"QLQ": set(), "WRC": set(), "SSXT": set()})
    for (genome, gid), doms in gene_doms.items():
        if "QLQ" in doms:
            per_genome[genome]["QLQ"].add(gid)
        if "WRC" in doms:
            per_genome[genome]["WRC"].add(gid)
        if "SSXT" in doms:
            per_genome[genome]["SSXT"].add(gid)

    rows = []
    for genome in sorted(per_genome):
        sets = per_genome[genome]
        grf = sets["QLQ"] & sets["WRC"]
        ssxt = sets["SSXT"]
        rows.append({
            "genome": genome,
            "GRF":  len(grf),
            "QLQ":  len(sets["QLQ"]),
            "WRC":  len(sets["WRC"]),
            "GIF":  len(ssxt),
            "SSXT": len(ssxt),
        })
    rows.sort(key=lambda r: (-(r["GRF"] + r["GIF"]), r["genome"]))
    return rows


def write_family_classification(records, out_path):
    """Write GRF/GIF family classification CSV. Returns (rows, csv_path) or (None, None)."""
    rows = compute_family_classification(records)
    if not rows:
        return None, None
    # Only emit when at least one GRF/GIF-defining domain was observed
    domains_seen = {r["domain_name"] for r in records}
    if not (domains_seen & {"QLQ", "WRC", "SSXT"}):
        return None, None
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FAMILY_FIELDS)
        w.writeheader()
        w.writerows(rows)
    return rows, out_path


def generate_report(records, prot_stats, trans_stats, nhmmer_records,
                    source_breakdown, gene_group, evalue, report_dir,
                    gene_list_written=None, family_rows=None, orthologs=None):
    """Write report.md summarising everything (protein + nucleotide)."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    genomes = sorted({r["genome"] for r in records})
    domains = sorted({r["domain_name"] for r in records})

    # --- Per-genome per-domain summary ---
    genome_domain = defaultdict(lambda: defaultdict(list))
    for r in records:
        genome_domain[r["genome"]][r["domain_name"]].append(r)

    overview_rows = []
    for g in genomes:
        for d in domains:
            hits = genome_domain[g].get(d, [])
            unique_genes = len({h["gene_id"] for h in hits})
            best_ev = min((h["full_seq_evalue"] for h in hits), default="—")
            best_score = max((h["full_seq_score"] for h in hits), default="—")
            overview_rows.append([g, d, len(hits), unique_genes,
                                  f"{best_ev:.2e}" if isinstance(best_ev, float) else best_ev,
                                  f"{best_score:.1f}" if isinstance(best_score, float) else best_score])

    # --- Top hits table (top 20 by score) ---
    # Scale top hits with genome count so large runs stay informative
    top_k = min(50, max(20, len(genomes)))
    top = sorted(records, key=lambda r: r["full_seq_score"], reverse=True)[:top_k]
    top_rows = []
    for r in top:
        top_rows.append([
            r["genome"], r["gene_id"], r["domain_name"],
            f"{r['full_seq_evalue']:.2e}", f"{r['full_seq_score']:.1f}",
            r["protein_length"],
            f"{r['ali_from']}–{r['ali_to']}",
            f"{r['accuracy']:.2f}",
            r["annotation"][:60] + ("…" if len(r["annotation"]) > 60 else ""),
        ])

    # --- CD-HIT summary (protein) ---
    prot_cdhit_rows = []
    for (g, p) in sorted(prot_stats):
        s = prot_stats[(g, p)]
        reduced = s["before"] - s["after"]
        pct = f"{reduced / s['before'] * 100:.1f}" if s["before"] > 0 else "0"
        prot_cdhit_rows.append([g, p, s["before"], s["after"], reduced, f"{pct}%"])

    # --- CD-HIT-EST summary (transcripts) ---
    trans_cdhit_rows = []
    for (g, p) in sorted(trans_stats):
        s = trans_stats[(g, p)]
        reduced = s["before"] - s["after"]
        pct = f"{reduced / s['before'] * 100:.1f}" if s["before"] > 0 else "0"
        trans_cdhit_rows.append([g, p, s["before"], s["after"], reduced, f"{pct}%"])

    # --- nhmmer summary ---
    nhmmer_genomes = sorted({r["genome"] for r in nhmmer_records}) if nhmmer_records else []
    nhmmer_profiles = sorted({r["profile"] for r in nhmmer_records}) if nhmmer_records else []
    nhmmer_overview_rows = []
    if nhmmer_records:
        ng_prof = defaultdict(lambda: defaultdict(list))
        for r in nhmmer_records:
            ng_prof[r["genome"]][r["profile"]].append(r)
        for g in nhmmer_genomes:
            for p in nhmmer_profiles:
                hits = ng_prof[g].get(p, [])
                unique_targets = len({h["target_name"] for h in hits})
                best_ev = min((h["evalue"] for h in hits), default=0)
                best_sc = max((h["score"] for h in hits), default=0)
                nhmmer_overview_rows.append([
                    g, p, len(hits), unique_targets,
                    f"{best_ev:.2e}" if best_ev else "—",
                    f"{best_sc:.1f}" if best_sc else "—",
                ])

    # --- Multi-domain genes ---
    gene_doms = defaultdict(set)
    for r in records:
        gene_doms[(r["genome"], r["gene_id"])].add(r["domain_name"])
    multi = {k: v for k, v in gene_doms.items() if len(v) > 1}
    multi_rows = [[g, gid, ", ".join(sorted(ds))] for (g, gid), ds in sorted(multi.items())]

    # --- Cross-genome hit count matrix (genome × domain pivot) ---
    all_domains = sorted({r["domain_name"] for r in records})
    pivot_rows = []
    for g in genomes:
        g_recs = [r for r in records if r["genome"] == g]
        _total_g = len({r["gene_id"] for r in g_recs})
        # After CD-HIT: max per-profile representative count for this genome.
        # Using max correctly handles redundant profiles (e.g. DMP + PF05078_DMP
        # that target the same gene set and produce the same after count).
        g_after = max(
            (s["after"] for (gg, _), s in prot_stats.items() if gg == g),
            default=_total_g,
        )
        row = [g, _total_g, g_after]
        for d in all_domains:
            row.append(len({r["gene_id"] for r in g_recs if r["domain_name"] == d}))
        pivot_rows.append((_total_g, row))
    pivot_rows.sort(key=lambda r: r[0], reverse=True)
    pivot_rows = [r[1] for r in pivot_rows]
    pivot_headers = ["Genome", "HMMER Search", "After CD-HIT"] + all_domains

    # --- Compact per-genome gene list summary (full lists written to gene_lists/ CSVs) ---
    compact_genome_rows = []
    for g in genomes:
        g_recs = [r for r in records if r["genome"] == g]
        unique = len({r["gene_id"] for r in g_recs})
        top_r  = max(g_recs, key=lambda r: r["full_seq_score"]) if g_recs else None
        sanitized = re.sub(r'[^\w.-]', '_', g)
        csv_ref = f"`gene_lists/{sanitized}_genes.csv`"
        compact_genome_rows.append([
            g, unique,
            f"{top_r['full_seq_evalue']:.2e}" if top_r else "—",
            f"{top_r['full_seq_score']:.1f}"  if top_r else "—",
            top_r["gene_id"]    if top_r else "—",
            top_r["domain_name"] if top_r else "—",
            csv_ref,
        ])
    compact_genome_rows.sort(key=lambda r: -r[1])  # sort by unique gene count desc

    total_hits = len(records)
    total_unique = len({r["gene_id"] for r in records})
    total_nhmmer_hits = len(nhmmer_records) if nhmmer_records else 0
    total_nhmmer_targets = len({r["target_name"] for r in nhmmer_records}) if nhmmer_records else 0

    # --- Source annotation ---
    has_source = bool(source_breakdown)
    search_mode = "dual (protein + nucleotide)" if total_nhmmer_hits > 0 else "protein-only"

    # Family classification comes first when available.
    # fam_offset shifts all subsequent section numbers by 1.
    fam_offset = 1 if family_rows else 0

    # Sort family rows: Solanum_melongena_v4.1 always first, then by (GRF+GIF) desc.
    def _fam_sort_key(r):
        is_smel = 0 if "Solanum_melongena" in r["genome"] else 1
        return (is_smel, -(r["GRF"] + r["GIF"]), r["genome"])

    report = f"""\
# HMMER Identification Report — {gene_group}

| Field | Value |
| --- | --- |
| Generated | {now} |
| E-value threshold | {evalue} |
| Search mode | {search_mode} |
| Genomes screened | {len(genomes)} |
| Total domain hits (protein) | {total_hits} |
| Total nhmmer hits (nucleotide) | {total_nhmmer_hits} |
| Unique genes identified (protein) | {total_unique} |
| Unique transcripts matched (nucleotide) | {total_nhmmer_targets} |

---

"""

    if family_rows:
        fam_sorted = sorted(family_rows, key=_fam_sort_key)
        family_md_rows = [
            [r["genome"], r["GRF"], r["QLQ"], r["WRC"], r["GIF"], r["SSXT"]]
            for r in fam_sorted
        ]
        report += f"""## 1. GRF / GIF Family Classification

Per-genome counts of unique genes classified into the **GRF** (QLQ + WRC
domains) and **GIF** (SSXT domain) families.

- **GRF**: unique genes containing **both** QLQ **and** WRC domains (canonical GRF)
- **QLQ**: unique genes containing a QLQ domain hit
- **WRC**: unique genes containing a WRC domain hit
- **GIF**: unique genes containing an SSXT domain hit (canonical GIF)
- **SSXT**: unique genes containing an SSXT domain hit

Genomes: *Solanum melongena* v4.1 and Unito Genomics Accessions.
Full table also written to `family_classification.csv` and `family_classification.tsv`.

{_md_table(["Genome", "GRF", "QLQ", "WRC", "GIF", "SSXT"], family_md_rows)}

"""

    # Ortholog mapping table rendered below; prepare rows here so a failure
    # in table assembly does not abort the whole report.
    ortholog_headers = None
    ortholog_rows = None
    if orthologs:
        ortholog_headers, ortholog_rows = generate_ortholog_table(orthologs)

    # Running counter ensures sequential section numbers regardless of which
    # optional sections (nhmmer, source breakdown) are present.
    sec = 1 + fam_offset

    report += f"""## {sec}. Cross-Genome Domain Count Matrix

**HMMER Search**: unique gene IDs from hmmsearch output (before CD-HIT).
**After CD-HIT**: representative sequences retained after redundancy reduction (max across profiles).
Domain columns: unique gene IDs per domain profile.

{_md_table(pivot_headers, pivot_rows)}

"""
    sec += 1

    if ortholog_headers and ortholog_rows:
        report += f"""{_LANDSCAPE_START_MARKER}

## {sec}. Gene Orthologs Across Genomes

Cross-genome ortholog mapping with columns rooted on *Solanum melongena*
v4.1 (BLAST reference). Rows: *Solanum melongena* V3 appears as the
second row (first data row), *Solanum melongena* v4.1 appears as the
third row, and remaining genomes follow in alphabetical order. Columns
are numbered 1..N for compactness — the v4.1 row (third row) carries the
reference gene IDs that each numbered column stands for. Each cell
reports the best 1:1 BLASTn-based ortholog assignment of the row genome
to the column v4.1 reference gene (E-value ≤ 1e-5), prioritising
reciprocal best hits and then filling remaining columns via bipartite
matching on bit-score so near-identical v4.1 paralogs each get a distinct
ortholog when the target genome has enough genes. A dash (-) indicates
no candidate match was found in that genome.

{_md_table(ortholog_headers, ortholog_rows)}

{_LANDSCAPE_END_MARKER}

"""
        sec += 1

    report += f"""## {sec}. Overview — Protein Search (hmmsearch)

{_md_table(
    ["Genome", "Domain", "Total Hits", "Unique Genes", "Best E-value", "Best Score"],
    overview_rows,
)}

"""
    sec += 1

    if nhmmer_overview_rows:
        report += f"""## {sec}. Overview — Nucleotide Search (nhmmer)

{_md_table(
    ["Genome", "Profile", "Total Hits", "Unique Targets", "Best E-value", "Best Score"],
    nhmmer_overview_rows,
)}

"""
        sec += 1

    section_num = sec

    report += f"""## {section_num}. Multi-domain Genes

{_md_table(["Genome", "Gene ID", "Domains"], multi_rows) if multi_rows else "_No genes with multiple domain types detected._"}

---
"""
    report_path = report_dir / "report.md"
    report_path.write_text(report, encoding="utf-8")
    return report_path


# ---------------------------------------------------------------------------
# Cross-gene-group report
# ---------------------------------------------------------------------------

def generate_cross_group_report(result_base, gene_groups, target_labels, evalue,
                                output_dir):
    """Scan multiple gene group result trees and produce a consolidated report.

    Outputs:
      {output_dir}/cross_group_summary.csv     one row per genome
      {output_dir}/cross_group_genes.csv       one row per discovered gene
      {output_dir}/cross_group_report.md       Markdown summary
    """
    result_base = Path(result_base)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Scan each gene group's identification directory
    all_records = {}  # group -> list of record dicts
    for group in gene_groups:
        ident_dir = result_base / group / "01_Identification"
        if ident_dir.is_dir():
            recs = scan_domtbl_files(ident_dir)
            all_records[group] = recs
            print(f"  {group}: {len(recs)} domain hits in "
                  f"{len({r['genome'] for r in recs})} genome(s)")
        else:
            all_records[group] = []
            print(f"  {group}: identification directory not found")

    # Collect all genomes observed + target labels (for zero-hit genomes)
    all_genomes = set()
    for recs in all_records.values():
        all_genomes.update(r["genome"] for r in recs)
    if target_labels:
        all_genomes.update(target_labels)
    all_genomes = sorted(all_genomes)

    # ---- CSV 1: Cross-group genome summary ----
    summary_fields = ["genome", "gene_groups_found"]
    for g in gene_groups:
        summary_fields.extend([f"{g}_unique_genes", f"{g}_best_evalue", f"{g}_best_score"])
    summary_fields.append("total_unique_genes")

    summary_rows = []
    for genome in all_genomes:
        row = {"genome": genome, "gene_groups_found": 0, "total_unique_genes": 0}
        total = 0
        found = 0
        for group in gene_groups:
            g_recs = [r for r in all_records[group] if r["genome"] == genome]
            unique = len({r["gene_id"] for r in g_recs})
            if unique > 0:
                found += 1
            best_ev = min((r["full_seq_evalue"] for r in g_recs), default=None)
            best_sc = max((r["full_seq_score"] for r in g_recs), default=None)
            row[f"{group}_unique_genes"] = unique
            row[f"{group}_best_evalue"] = f"{best_ev:.2e}" if best_ev is not None else ""
            row[f"{group}_best_score"] = f"{best_sc:.1f}" if best_sc is not None else ""
            total += unique
        row["gene_groups_found"] = found
        row["total_unique_genes"] = total
        summary_rows.append(row)
    summary_rows.sort(key=lambda r: (-r["total_unique_genes"], r["genome"]))

    summary_csv = output_dir / "cross_group_summary.csv"
    with open(summary_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=summary_fields)
        w.writeheader()
        w.writerows(summary_rows)
    _write_tsv_copy(summary_csv)
    print(f"  ✓ cross_group_summary.csv (+ .tsv, {len(summary_rows)} genomes)")

    # ---- CSV 2: Detailed gene list (one row per discovered gene) ----
    gene_fields = [
        "genome", "gene_group", "gene_id", "domain_name", "pfam_id",
        "full_seq_evalue", "full_seq_score", "domain_score",
        "protein_length", "accuracy", "gene_name", "annotation",
    ]
    gene_rows = []
    for group in gene_groups:
        # Deduplicate: keep highest-scoring domain hit per (genome, gene_id)
        seen = {}
        for r in sorted(all_records[group], key=lambda x: x["full_seq_score"],
                        reverse=True):
            key = (r["genome"], r["gene_id"])
            if key not in seen:
                seen[key] = r
        for (genome, gene_id), r in sorted(seen.items()):
            gene_rows.append({
                "genome": genome,
                "gene_group": group,
                "gene_id": gene_id,
                "domain_name": r["domain_name"],
                "pfam_id": r.get("pfam_id", ""),
                "full_seq_evalue": f"{r['full_seq_evalue']:.2e}",
                "full_seq_score": f"{r['full_seq_score']:.1f}",
                "domain_score": f"{r['domain_score']:.1f}",
                "protein_length": r["protein_length"],
                "accuracy": f"{r['accuracy']:.2f}",
                "gene_name": r.get("gene_name", ""),
                "annotation": r.get("annotation", ""),
            })

    genes_csv = output_dir / "cross_group_genes.csv"
    with open(genes_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=gene_fields)
        w.writeheader()
        w.writerows(gene_rows)
    _write_tsv_copy(genes_csv)
    print(f"  ✓ cross_group_genes.csv (+ .tsv, {len(gene_rows)} genes)")

    # ---- Markdown report ----
    md_headers = ["Genome", "Groups Found"]
    for g in gene_groups:
        md_headers.extend([f"{g} Genes", f"{g} Best E-value", f"{g} Best Score"])
    md_headers.append("**Total Genes**")

    md_rows = []
    for sr in summary_rows:
        row = [sr["genome"], sr["gene_groups_found"]]
        for g in gene_groups:
            row.append(sr[f"{g}_unique_genes"])
            row.append(sr[f"{g}_best_evalue"] or "—")
            row.append(sr[f"{g}_best_score"] or "—")
        row.append(sr["total_unique_genes"])
        md_rows.append(row)

    # Top genes table (top 50 by score across all groups)
    all_flat = []
    for group in gene_groups:
        for r in all_records[group]:
            all_flat.append({**r, "_group": group})
    all_flat.sort(key=lambda x: x["full_seq_score"], reverse=True)
    top_k = min(50, len(all_flat))
    top_rows = []
    seen_genes = set()
    for r in all_flat:
        if len(top_rows) >= top_k:
            break
        key = (r["genome"], r["gene_id"], r["_group"])
        if key in seen_genes:
            continue
        seen_genes.add(key)
        top_rows.append([
            r["genome"], r["_group"], r["gene_id"], r["domain_name"],
            f"{r['full_seq_evalue']:.2e}", f"{r['full_seq_score']:.1f}",
            r["protein_length"],
            f"{r['accuracy']:.2f}",
            r.get("annotation", "")[:50],
        ])

    genomes_with_hits = sum(1 for sr in summary_rows if sr["gene_groups_found"] > 0)

    report_text = f"""\
# Cross-Gene-Group HMMER Identification Report

**Generated:** {now}
**E-value threshold:** {evalue}
**Gene groups analyzed:** {', '.join(gene_groups)}
**Target genomes:** {len(all_genomes)}
**Genomes with hits:** {genomes_with_hits}

---

## 1. Genome × Gene Group Summary

All target genomes with per-group unique gene counts and best scores.
Sorted by total unique genes (highest first). Genomes with zero hits
across all groups are included at the bottom for completeness.

{_md_table(md_headers, md_rows)}

## 2. Top {top_k} Hits by Bit-Score (All Groups)

Deduplicated by (genome, gene, group). Sorted by full-sequence bit-score.

{_md_table(
    ["Genome", "Gene Group", "Gene ID", "Domain", "E-value", "Score",
     "Length", "Accuracy", "Annotation"],
    top_rows,
) if top_rows else "_No hits found across any gene group._"}

---
"""
    report_path = output_dir / "cross_group_report.md"
    report_path.write_text(report_text, encoding="utf-8")
    print(f"  ✓ cross_group_report.md")
    return report_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate HMMER identification report with tables (summary.csv + report.md).",
    )
    parser.add_argument("ident_dir", nargs="?",
                        help="Path to 01_Identification directory (per-group mode)")
    parser.add_argument("--gene-group", default="unknown", help="Gene group label")
    parser.add_argument("--evalue", default="auto", help="E-value used (for report text)")
    # Cross-group mode
    parser.add_argument("--cross-group", action="store_true",
                        help="Generate cross-gene-group summary report")
    parser.add_argument("--result-base",
                        help="Path to 3_RESULT base directory (cross-group mode)")
    parser.add_argument("--gene-groups", nargs="+",
                        help="Gene group names to include (cross-group mode)")
    parser.add_argument("--output",
                        help="Output directory for cross-group report")
    parser.add_argument("--target-labels", nargs="*", default=[],
                        help="All target genome labels (includes zero-hit genomes)")
    parser.add_argument("--threads", type=int, default=0,
                        help="Threads for BLAST-based ortholog table (default: auto)")
    args = parser.parse_args()

    # --- Cross-group mode ---
    if args.cross_group:
        if not args.result_base or not args.gene_groups:
            parser.error("--cross-group requires --result-base and --gene-groups")
        result_base = Path(args.result_base)
        output_dir = (Path(args.output) if args.output
                      else result_base / "cross_group_hmmer_report")
        evalue = args.evalue if args.evalue != "auto" else "N/A"
        print(f"Cross-group report for: {', '.join(args.gene_groups)}")
        generate_cross_group_report(
            result_base, args.gene_groups, args.target_labels,
            evalue, output_dir,
        )
        print(f"Cross-group report complete: {output_dir}")
        return

    # --- Per-group mode (existing behavior) ---
    if not args.ident_dir:
        parser.error("ident_dir is required in per-group mode")

    ident_dir = Path(args.ident_dir)
    if not ident_dir.is_dir():
        print(f"ERROR: {ident_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning {ident_dir} for HMMER results...")
    file_index = _index_result_files(ident_dir)

    records = scan_domtbl_files(ident_dir, file_index["domtbl"])
    if not records:
        print("WARNING: No filtered domtbl records found. Skipping report generation.")
        sys.exit(0)

    nhmmer_records = scan_nhmmer_tbl_files(ident_dir, file_index["nhmmer_tbl"])
    source_breakdown = scan_hit_source_files(
        ident_dir, file_index["prot_hit_ids"], file_index["nucl_hit_ids"])
    prot_stats, trans_stats = scan_clstr_files(
        ident_dir, file_index["clstr"], file_index["cdhit_fa"])

    # Annotate each protein record with search_source by cross-referencing hit ID files
    for rec in records:
        rec["search_source"] = "protein"  # default
    if source_breakdown:
        # Build a fast lookup: gene_id -> source
        source_lookup = {}  # (genome, gene_id) -> "protein" | "nucleotide" | "both"
        for (g, _prof), bd in source_breakdown.items():
            for gid in bd["prot_only"]:
                source_lookup[(g, gid)] = "protein"
            for gid in bd["nucl_only"]:
                source_lookup[(g, gid)] = "nucleotide"
            for gid in bd["both"]:
                source_lookup[(g, gid)] = "both"
        for rec in records:
            rec["search_source"] = source_lookup.get(
                (rec["genome"], rec["gene_id"]), "protein")

    # Auto-detect e-value from paths if not given
    if args.evalue == "auto":
        ev_dirs = [d.name for d in ident_dir.rglob("e-value_*") if d.is_dir()]
        args.evalue = ev_dirs[0].replace("e-value_", "") if ev_dirs else "N/A"

    report_dir = ident_dir / "d_REPORT"
    report_dir.mkdir(parents=True, exist_ok=True)

    print(f"  Found {len(records)} protein domain hits across {len({r['genome'] for r in records})} genome(s)")
    if nhmmer_records:
        print(f"  Found {len(nhmmer_records)} nhmmer nucleotide hits across {len({r['genome'] for r in nhmmer_records})} genome(s)")
    if source_breakdown:
        total_prot = sum(len(v["prot_only"]) for v in source_breakdown.values())
        total_nucl = sum(len(v["nucl_only"]) for v in source_breakdown.values())
        total_both = sum(len(v["both"]) for v in source_breakdown.values())
        print(f"  Hit sources: {total_prot} protein-only, {total_nucl} nucleotide-only, {total_both} both")
    print(f"  Writing to {report_dir}/")

    # summary.csv — all domain hits consolidated
    summary_csv = report_dir / "summary.csv"
    write_summary_csv(records, summary_csv)
    _write_tsv_copy(summary_csv)
    print("  ✓ summary.csv (+ summary.tsv)")

    # per_genome_stats.csv — one aggregate row per genome (useful for 30+ genome runs)
    pgs_csv = report_dir / "per_genome_stats.csv"
    write_per_genome_stats_csv(
        records, nhmmer_records, prot_stats, trans_stats, source_breakdown,
        pgs_csv,
    )
    _write_tsv_copy(pgs_csv)
    print("  ✓ per_genome_stats.csv (+ per_genome_stats.tsv)")

    # family_classification.csv — GRF (QLQ+WRC) / GIF (SSXT) per genome
    family_csv = report_dir / "family_classification.csv"
    family_rows, family_path = write_family_classification(records, family_csv)
    if family_path:
        _write_tsv_copy(family_csv)
        print("  ✓ family_classification.csv (+ family_classification.tsv)")

    # gene_lists/ — one CSV per genome (keeps report.md readable with many genomes)
    gene_list_dir = report_dir / "gene_lists"
    gene_list_written = write_gene_list_csvs(records, gene_list_dir)
    print(f"  ✓ gene_lists/ ({len(gene_list_written)} genome CSV(s) + TSV(s))")

    # ortholog mapping — cross-genome BLASTn + reciprocal best hits
    import os as _os
    if args.threads and args.threads > 0:
        _ortho_threads = args.threads
    else:
        _ortho_threads = min(8, max(1, (_os.cpu_count() or 4)))
    orthologs = build_ortholog_mapping(ident_dir, args.gene_group, threads=_ortho_threads)

    # orthologs_table.csv + .tsv — matches the Markdown section 2 layout
    if orthologs:
        ortho_headers, ortho_rows = generate_ortholog_table(orthologs)
        ortho_csv = report_dir / "orthologs_table.csv"
        write_ortholog_table_csv(ortho_headers, ortho_rows, ortho_csv)
        _write_tsv_copy(ortho_csv)
        print("  ✓ orthologs_table.csv (+ orthologs_table.tsv)")

    # report.md — Markdown report with cross-genome matrix and compact gene list summary
    rpath = generate_report(records, prot_stats, trans_stats, nhmmer_records,
                            source_breakdown, args.gene_group, args.evalue, report_dir,
                            gene_list_written=gene_list_written,
                            family_rows=family_rows, orthologs=orthologs)
    print(f"  ✓ {rpath.name}")

    # report.pdf — optional; skipped gracefully when no PDF tool is installed
    pdf_path = _write_report_pdf(rpath)
    if pdf_path:
        print(f"  ✓ {pdf_path.name}")
    else:
        print("  ! report.pdf skipped — install pandoc or weasyprint+markdown to enable PDF output")

    print(f"Report complete: {report_dir}")


if __name__ == "__main__":
    main()
