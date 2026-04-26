#!/usr/bin/env python3
"""
Shared DMP query-ID -> publication-friendly label mappings.

This is the single source of truth for the mapping tables and the
``short_label()`` / ``extract_species()`` helpers used to convert raw
NCBI/Phytozome FASTA header IDs (as they appear in the BLAST query files)
into compact, human-readable names like ``CaDMP3-like`` or ``MtDMP9``.

Used by:
  - ``modules/02_blast_ortholog_alignment/visualize_blast_results.py``
    (lollipop / heatmap y-axis labels)
  - ``modules/04_multiple_sequence_alignment/build_v4_blast_groups.py``
    (per-paralog FASTA headers)
  - ``modules/utils/prettify_fasta_headers.py``
    (CLI to rewrite headers in any FASTA)

Add new mappings here, not in the BLAST viz module.
"""
from __future__ import annotations

import re

# ---------------------------------------------------------------------------
# Anchor SmelDMP single-line labels (axis tick / FASTA header form)
# ---------------------------------------------------------------------------
GENE_SHORT: dict[str, str] = {
    "SMEL5_01g026030.1": "SmelDMP01.030",
    "SMEL5_01g008730.1": "SmelDMP01.730",
    "SMEL5_04g005390.1": "SmelDMP04.390",
    "SMEL5_02g013320.1": "SmelDMP02.320",
    "SMEL5_10g003660.1": "SmelDMP10.660",
    "SMEL5_12g005350.1": "SmelDMP12.350",
    "SMEL5_10g017610.1": "SmelDMP10.610",
}

# ---------------------------------------------------------------------------
# Mapping tables for cross-species BLAST query IDs
# ---------------------------------------------------------------------------
# Hardcoded mapping for raw NCBI/genome-DB query IDs that appear in the curated
# CSV without a species-prefixed name. Resolved from the FASTA headers under
# I_RefSeqs/d_DMP_Query_Fasta/<species>/*_merged*.fa (protein description in
# square brackets) plus the species inferred from the parent directory.
# Format: accession_stem -> friendly short label.
ACCESSION_LABELS: dict[str, str] = {
    # Nicotiana tabacum (NtDMPs_merged.fa) - RefSeq XM_016*
    "XM_016586331": "NtDMP2-like",
    "XM_016642301": "NtDMP3-like",
    "XM_016578562": "NtDMP4-like",
    "XM_016591727": "NtDMP7-like",
    "XM_016604032": "NtDMP9-like",
    # Cucumis sativus (CsDMPs_merged.fa)
    "XM_004146681": "CsDMP9",
    # Medicago truncatula (MtDMPs_merged.fa)
    "XM_003621193": "MtDMP9",   # haploid inducer (N. Wang et al., 2022)
    "XM_003614037": "MtDMP",
    # Brassica oleracea (BoDMPs_merged.fa) - uncharacterized LOC*
    "XM_013772041": "BoDMP_LOC106333617",
    "XM_013772244": "BoDMP_LOC106333853",
    # Solanum lycopersicum RefSeq NM_/NP_ accessions
    "NM_001150103": "SlDMP2",
    "NM_001150831": "SlDMP3-like",
    "NM_001151404": "SlDMP8",
    "NM_001151895": "SlDMP4-like",
    "NM_001371746": "SlDMP1-like",
}

# Protein accession (XP_*) -> short label, derived from header scans of the 8
# gap-filler species under II_INPUTS/DMP_query_fasta_file/. Used when BLAST
# emits the secondary `lcl|NC_*_cds_XP_*` token form as Query ID, or when a
# Query ID is just a bare XP accession.
PROTEIN_ACC_LABELS: dict[str, str] = {
    # Pp - Physcomitrella patens (1)
    "XP_024392000": "PpDMP5-like",
    # Si - Setaria italica (14)
    "XP_004952503": "SiDMP6", "XP_004958514": "SiDMP1", "XP_004961228": "SiDMP4",
    "XP_004965304": "SiDMP2", "XP_004965424": "SiDMP4", "XP_004966913": "SiDMP10",
    "XP_004968801": "SiDMP3", "XP_004969121": "SiDMP3", "XP_004969493": "SiDMP6",
    "XP_004970780": "SiDMP7", "XP_004971846": "SiDMP2", "XP_004972417": "SiDMP3",
    "XP_004984268": "SiDMP2", "XP_022681673": "SiDMP9-like",
    # Ma - Musa acuminata (42)
    "XP_009380593": "MaDMP6-like", "XP_009380594": "MaDMP4",
    "XP_009381725": "MaDMP8-like", "XP_009385032": "MaDMP2",
    "XP_009387184": "MaDMP5-like", "XP_009394826": "MaDMP4",
    "XP_009395997": "MaDMP4",      "XP_009402533": "MaDMP3-like",
    "XP_009402791": "MaDMP7-like", "XP_009403823": "MaDMP2-like",
    "XP_009408744": "MaDMP3-like", "XP_009410685": "MaDMP3",
    "XP_009413807": "MaDMP4",      "XP_064937510": "MaDMP4-like",
    "XP_064937589": "MaDMP7-like", "XP_064943435": "MaDMP8-like",
    "XP_064943850": "MaDMP4-like", "XP_064943851": "MaDMP6-like",
    "XP_064944811": "MaDMP2-like", "XP_064958017": "MaDMP4-like",
    "XP_064968105": "MaDMP2-like", "XP_064972103": "MaDMP5-like",
    "XP_064979449": "MaDMP7-like", "XP_065001710": "MaDMP4-like",
    "XP_065003221": "MaDMP5-like", "XP_065007018": "MaDMP3-like",
    "XP_065007094": "MaDMP7-like", "XP_065012787": "MaDMP6-like",
    "XP_065015983": "MaDMP3-like", "XP_065016215": "MaDMP5-like",
    "XP_065017841": "MaDMP7-like", "XP_065019000": "MaDMP4-like",
    "XP_065027107": "MaDMP4-like", "XP_065027109": "MaDMP6-like",
    "XP_065027149": "MaDMP8-like", "XP_065036817": "MaDMP4-like",
    "XP_065036951": "MaDMP4-like", "XP_065037018": "MaDMP5-like",
    "XP_065041527": "MaDMP7-like", "XP_065041528": "MaDMP3-like",
    "XP_065042776": "MaDMP2-like", "XP_065049239": "MaDMP3-like",
    # Pt - Populus trichocarpa (11)
    "XP_002305111": "PtDMP4", "XP_002312247": "PtDMP8", "XP_002312376": "PtDMP2",
    "XP_002315049": "PtDMP8", "XP_002315542": "PtDMP2", "XP_006372758": "PtDMP7",
    "XP_006376329": "PtDMP2", "XP_006377414": "PtDMP4", "XP_006385029": "PtDMP3",
    "XP_006389536": "PtDMP3", "XP_024439153": "PtDMP10",
    # Bd - Brachypodium distachyon (13)
    "XP_003557584": "BdDMP3", "XP_003557899": "BdDMP2", "XP_003565771": "BdDMP2",
    "XP_003566981": "BdDMP4", "XP_003567905": "BdDMP6", "XP_003572244": "BdDMP2",
    "XP_003572527": "BdDMP6", "XP_010228847": "BdDMP2", "XP_010236366": "BdDMP2",
    "XP_024314432": "BdDMP10-like-X1", "XP_024314433": "BdDMP10-like-X2",
    "XP_024316011": "BdDMP7", "XP_024316012": "BdDMP7",
    # Gr - Gossypium raimondii (9)
    "XP_012453761": "GrDMP2", "XP_012457653": "GrDMP4", "XP_012457725": "GrDMP2",
    "XP_012458528": "GrDMP7", "XP_012467832": "GrDMP3", "XP_012478015": "GrDMP10",
    "XP_012481184": "GrDMP2", "XP_012487405": "GrDMP2", "XP_012489581": "GrDMP9",
    # Ca - Capsicum annuum (7)
    "XP_016539005": "CaDMP3",      "XP_016540695": "CaDMP7-like",
    "XP_016557530": "CaDMP2",      "XP_016571689": "CaDMP9-like",
    "XP_016582229": "CaDMP10",     "XP_016582285": "CaDMP3-like",
    "XP_047261276": "CaDMP6-like",
    # NCBI NP_* (Solanum lycopersicum, in NM_/NP_ pairs)
    "NP_001143575": "SlDMP2",
    "NP_001144303": "SlDMP3-like",
    "NP_001144876": "SlDMP8",
    "NP_001145367": "SlDMP4-like",
    "NP_001358675": "SlDMP1-like",
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def short_label(query_id: str) -> str:
    """Return a publication-ready label for a BLAST query ID or FASTA header.

    Resolves the common raw-ID styles, falling back to a prefix-strip regex:

      'lcl|XM_003621193.2_cds_XP_003621241.1_1'    -> 'MtDMP9'         (ACCESSION_LABELS, XM key)
      'lcl|NM_001151404.2_cds_NP_001144876.2_1'    -> 'SlDMP8'         (ACCESSION_LABELS, NM key)
      'lcl|NC_037265.1_cds_XP_024392000.1_26917'   -> 'PpDMP5-like'    (PROTEIN_ACC_LABELS, secondary-token)
      'XP_024392000.1'                              -> 'PpDMP5-like'    (PROTEIN_ACC_LABELS, bare-XP)
      'CDX74441.'                                   -> 'CDX74441'       (strip trailing '.')
      'mRNA.BjuA04g10430S.'                         -> 'BjuA04g10430S'  (strip 'mRNA.' prefix)
      'PpDMP5-like_XP_024392000.1_26917'            -> 'PpDMP5-like'    (regex split on '_XP')
      'AtDMP8+AtDMP9'                               -> 'AtDMP8+AtDMP9'  (pass-through)
      'SMEL5_10g017610.1'                           -> 'SmelDMP10.610'  (GENE_SHORT)
    """
    qid = query_id.strip().removeprefix("lcl|")
    qid = qid.rstrip(".")

    if qid in GENE_SHORT:
        return GENE_SHORT[qid]

    if qid.startswith("XM_"):
        m = re.match(r"(XM_\d+)", qid)
        if m:
            stem = m.group(1)
            if stem in ACCESSION_LABELS:
                return ACCESSION_LABELS[stem]
            return stem

    if qid.startswith("NM_"):
        m = re.match(r"(NM_\d+)", qid)
        if m and m.group(1) in ACCESSION_LABELS:
            return ACCESSION_LABELS[m.group(1)]

    # NCBI secondary-token form (e.g. headers retain lcl|NC_*_cds_XP_* as the
    # second whitespace token; appears as Query ID if the first token is dropped).
    m = re.match(r"NC_\d+\.\d+_cds_(XP_\d+)", qid)
    if m and m.group(1) in PROTEIN_ACC_LABELS:
        return PROTEIN_ACC_LABELS[m.group(1)]

    if qid.startswith("XP_"):
        m = re.match(r"(XP_\d+)", qid)
        if m and m.group(1) in PROTEIN_ACC_LABELS:
            return PROTEIN_ACC_LABELS[m.group(1)]

    if qid.startswith("NP_"):
        m = re.match(r"(NP_\d+)", qid)
        if m and m.group(1) in PROTEIN_ACC_LABELS:
            return PROTEIN_ACC_LABELS[m.group(1)]

    if qid.startswith("mRNA.Bju"):
        return qid[len("mRNA."):]

    label = re.split(r"[_\s](?:XP|XM|BSXM|NM|AT[0-9]|Glyma\.|OZ[0-9]|NC_[0-9]|LOC[0-9])", qid)[0]
    return label.rstrip("_.")


def extract_species(query_id: str) -> str | None:
    """Return the two-letter species prefix, e.g. 'Sl' from 'SlDMP3_XP_...'.
    Returns None when the species cannot be inferred.
    """
    qid = query_id.strip().removeprefix("lcl|")
    if qid.startswith(("CDX", "CDY")) or qid.startswith("mRNA.Bju"):
        return "Br"
    if qid.startswith("XM_"):
        m = re.match(r"(XM_\d+)", qid)
        if m and m.group(1) in ACCESSION_LABELS:
            return extract_species(ACCESSION_LABELS[m.group(1)])
    if qid.startswith("NM_"):
        m = re.match(r"(NM_\d+)", qid)
        if m and m.group(1) in ACCESSION_LABELS:
            return extract_species(ACCESSION_LABELS[m.group(1)])
    m = re.match(r"NC_\d+\.\d+_cds_(XP_\d+)", qid)
    if m and m.group(1) in PROTEIN_ACC_LABELS:
        return extract_species(PROTEIN_ACC_LABELS[m.group(1)])
    if qid.startswith("XP_"):
        m = re.match(r"(XP_\d+)", qid)
        if m and m.group(1) in PROTEIN_ACC_LABELS:
            return extract_species(PROTEIN_ACC_LABELS[m.group(1)])
    if qid.startswith("NP_"):
        m = re.match(r"(NP_\d+)", qid)
        if m and m.group(1) in PROTEIN_ACC_LABELS:
            return extract_species(PROTEIN_ACC_LABELS[m.group(1)])
    m = re.match(r"^([A-Z][a-z]+)", qid)
    return m.group(1) if m else None
