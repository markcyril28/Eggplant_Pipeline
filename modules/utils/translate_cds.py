#!/usr/bin/env python3
"""Translate a multi-record CDS FASTA to protein in frame 1 (stdlib only).

Used by `00_RefSeq_setup_and_download.sh` (operation `MERGE_DMP_HI_GENES`) to
derive `*_merged_protein.fa` sister files alongside the CDS `*_merged*.fa`,
so that the BLASTp stage can consume the verified haploid-inducer DMP panel
without an external translator (Biopython, EMBOSS, gffread).

Behaviour:
- Frame 1 only (records arrive as canonical NCBI/Ensembl mRNA CDS).
- Trailing stop codon trimmed; internal stops counted and reported.
- Records < 30 nt or yielding < 30 aa are skipped.
- Unknown codons (N-containing or ambiguous) translate to 'X'.

Usage:
    python3 translate_cds.py <in_cds.fa> <out_protein.fa>
"""
import sys

CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def translate(nuc):
    nuc = nuc.upper().replace("U", "T")
    aa = []
    for i in range(0, len(nuc) - 2, 3):
        c = nuc[i:i+3]
        aa.append(CODON.get(c, "X"))
    s = "".join(aa)
    if s.endswith("*"):
        s = s[:-1]
    return s

def parse_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
        if name:
            yield name, "".join(seq)

def main(in_fa, out_fa):
    n_in = n_out = warn = 0
    with open(out_fa, "w") as out:
        for hdr, nuc in parse_fasta(in_fa):
            n_in += 1
            if len(nuc) < 30:
                continue
            aa = translate(nuc)
            if "*" in aa:
                warn += 1
            if len(aa) >= 30:
                out.write(f">{hdr}\n")
                for i in range(0, len(aa), 60):
                    out.write(aa[i:i+60] + "\n")
                n_out += 1
    print(f"  translate_cds: {in_fa} -> {out_fa}  ({n_in} CDS -> {n_out} protein, {warn} internal-stop)", flush=True)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: translate_cds.py <in_cds.fa> <out_protein.fa>")
    main(sys.argv[1], sys.argv[2])
