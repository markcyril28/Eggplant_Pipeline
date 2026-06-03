#!/usr/bin/env python3
"""Compose the four-priority SmelDMP AlphaFold 3 figure.

Builds a vertical stack of four per-protein side-by-side composites
(pLDDT confidence | interaction-domain) for the priority paralogs only:

  (ii)  SmelDMPv5_02       Tier 2
  (iv)  SmelDMPv5_01.030   Tier 1
  (v)   SmelDMPv5_01.730   Tier 2
  (vii) SmelDMPv5_10.610   Tier 1, primary HI candidate

Order mirrors Table 10 (descending mean pLDDT). Panel numerals match the
seven-protein numbering used by ``compose_confidence_interaction.py`` so
cross-figure references (Appendix K, Appendix L) stay consistent.

Output:
  AlphaFold3_Results/Combined_Confidence_Interaction/
      four_priority_smeldmp_confidence_vs_interaction_composite.png
"""

from __future__ import annotations
import sys
from pathlib import Path

from PIL import Image, ImageDraw

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from compose_confidence_interaction import (  # noqa: E402
    AF3_DIR,
    BG_COLOR,
    FONT_SIZE,
    GRID_GAP,
    TEXT_COLOR,
    get_font,
    make_pair,
)

PRIORITY_PANELS = [
    {  # (ii) highest pLDDT
        "confidence_dir":  "2026_04_07_13_17_smel5_02g013320_1",
        "confidence_stem": "2026_04_07_13_17_smel5_02g013320_1",
        "interaction_dir":  "2026_04_07_13_17_smel5_02g013320_1",
        "interaction_stem": "2026_04_07_13_17_smel5_02g013320_1",
        "name": "SmelDMPv5_02",
        "gene": "SMEL5_02g013320.1",
        "panel": "ii",
        "tier": "Tier 2",
    },
    {  # (iv)
        "confidence_dir":  "2026_04_07_13_17_smel5_01g026030_1",
        "confidence_stem": "2026_04_07_13_17_smel5_01g026030_1",
        "interaction_dir":  "2026_04_07_13_17_smel5_01g026030_1",
        "interaction_stem": "2026_04_07_13_17_smel5_01g026030_1",
        "name": "SmelDMPv5_01.030",
        "gene": "SMEL5_01g026030.1",
        "panel": "iv",
        "tier": "Tier 1",
    },
    {  # (v)
        "confidence_dir":  "2026_04_07_13_16_smel5_01g008730_1",
        "confidence_stem": "2026_04_07_13_16_smel5_01g008730_1",
        "interaction_dir":  "smel5_01g008730_1",
        "interaction_stem": "smel5_01g008730_1",
        "name": "SmelDMPv5_01.730",
        "gene": "SMEL5_01g008730.1",
        "panel": "v",
        "tier": "Tier 2",
    },
    {  # (vii) primary HI candidate
        "confidence_dir":  "fold_smel5_10g017610_1",
        "confidence_stem": "fold_smel5_10g017610_1",
        "interaction_dir":  "smel5_10g017610_1",
        "interaction_stem": "smel5_10g017610_1",
        "name": "SmelDMPv5_10.610",
        "gene": "SMEL5_10g017610.1",
        "panel": "vii",
        "tier": "Tier 1, primary",
    },
]

OUTPUT_NAME = "four_priority_smeldmp_confidence_vs_interaction_composite.png"


def main() -> None:
    pipeline_dir = SCRIPT_DIR.parent.parent
    out_dir = pipeline_dir / AF3_DIR / "Combined_Confidence_Interaction"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Pipeline root: {pipeline_dir}")
    print(f"Output dir:    {out_dir}")
    print()

    composites = []
    for p in PRIORITY_PANELS:
        print(f"Processing {p['name']} ({p['gene']})...")
        comp = make_pair(pipeline_dir, p, out_dir, "black", "interaction")
        if comp is None:
            print(f"  ERROR: composite generation failed for {p['name']}")
            sys.exit(1)
        composites.append(comp)

    cw, ch = composites[0].size
    title_h = 110
    total_w = cw
    total_h = title_h + len(composites) * ch + (len(composites) - 1) * GRID_GAP

    canvas = Image.new("RGB", (total_w, total_h), BG_COLOR)
    draw = ImageDraw.Draw(canvas)
    font = get_font(FONT_SIZE + 4)
    title = (
        "Figure 12.  AlphaFold 3 SmelDMP Priority Structures: "
        "pLDDT Confidence vs Interaction Domain Coloring"
    )
    bbox = draw.textbbox((0, 0), title, font=font)
    draw.text(((total_w - (bbox[2] - bbox[0])) / 2, 12), title,
              fill=TEXT_COLOR, font=font)

    y = title_h
    for comp in composites:
        canvas.paste(comp, (0, y))
        y += ch + GRID_GAP

    out_path = out_dir / OUTPUT_NAME
    canvas.save(out_path)
    print()
    print(f"Saved: {out_path.relative_to(pipeline_dir)}")


if __name__ == "__main__":
    main()
