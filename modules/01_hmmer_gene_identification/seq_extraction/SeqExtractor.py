#!/usr/bin/env python3
"""
Extract FASTA records for a list of accessions from multiple FASTA files.

- Reads an accession list from a text file (one per line).
- Searches for each accession across all FASTA files (up to 4 levels deep)
  within a source directory.
- Writes matching records to an output FASTA file.
- If an accession is not found, appends only the header (>accession) to the output.

Usage:
    python SeqExtractor.py --source-folder /path/to/genomes [--phylo Debernardi Bull]
"""

import argparse
import sys
import os
from pathlib import Path


def read_fasta(filename):
    """Yield (header, sequence) tuples from a FASTA file."""
    header, seq_lines = None, []
    with open(filename, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header, seq_lines = line[1:], []
            else:
                seq_lines.append(line)
        if header is not None:
            yield header, "".join(seq_lines)


def write_fasta_record(out_handle, header, sequence):
    """Write a single FASTA record to the output handle."""
    out_handle.write(f"\n>{header}\n")
    if sequence:
        for start in range(0, len(sequence), 60):
            out_handle.write(sequence[start:start + 60] + "\n")
    out_handle.write("\n")
    out_handle.flush()


def find_fasta_files(source_folder, genome_type="cds"):
    """Return a list of FASTA files up to 4 levels deep with allowed extensions."""
    exts = (".fa", ".fasta", ".fna", ".faa")
    result = []
    for root, dirs, files in os.walk(source_folder):
        depth = Path(root).relative_to(source_folder).parts
        if len(depth) > 4:
            dirs[:] = []
            continue
        for fname in files:
            fname_lower = fname.lower()
            if genome_type == "protein":
                if ("protein" in fname_lower or "pep" in fname_lower or "peptide" in fname_lower) and fname_lower.endswith(exts):
                    result.append(os.path.join(root, fname))
            else:
                if "cds" in fname_lower and fname_lower.endswith(exts):
                    result.append(os.path.join(root, fname))
    return sorted(result)


# Per-phylo input/output path definitions
PHYLO_CONFIG = {
    "Debernardi": {
        "GRF": {
            "input": "LIST/Debernardi_List_for_extract_GRF.txt",
            "nucleotide": "OUTPUT_FASTA/Debernardi_Phylo_GRF_NucSeq.fasta",
            "protein": "OUTPUT_FASTA/Debernardi_Phylo_GRF_ProtSeq.fasta",
        },
        "GIF": {
            "input": "LIST/Debernardi_List_for_extract_GIF.txt",
            "nucleotide": "OUTPUT_FASTA/Debernardi_Phylo_GIF_NucSeq.fasta",
            "protein": "OUTPUT_FASTA/Debernardi_Phylo_GIF_ProtSeq.fasta",
        },
    },
    "Bull": {
        "GRF": {
            "input": "LIST/Bull_List_for_extract_GRF.txt",
            "nucleotide": "OUTPUT_FASTA/Bull_Phylo_GRF_NucSeq.fasta",
            "protein": "OUTPUT_FASTA/Bull_Phylo_GRF_ProtSeq.fasta",
        },
        "GIF": {
            "input": "LIST/Bull_List_for_extract_GIF.txt",
            "nucleotide": "OUTPUT_FASTA/Bull_Phylo_GIF_NucSeq.fasta",
            "protein": "OUTPUT_FASTA/Bull_Phylo_GIF_ProtSeq.fasta",
        },
    },
}


def main():
    parser = argparse.ArgumentParser(
        description="Extract FASTA records for accession lists from genome FASTA files."
    )
    parser.add_argument(
        "--source-folder", required=True,
        help="Root directory containing genome FASTA files."
    )
    parser.add_argument(
        "--phylo", nargs="+", default=["Debernardi", "Bull"],
        choices=["Debernardi", "Bull"],
        help="Which phylogeny sets to process (default: both)."
    )
    args = parser.parse_args()

    source_folder = args.source_folder

    for phylo in args.phylo:
        print(f"\n=== Processing for {phylo} Phylogeny ===")
        for gene_type in ["GIF", "GRF"]:
            cfg = PHYLO_CONFIG[phylo][gene_type]
            input_query_list = cfg["input"]

            with open(input_query_list, "r", encoding="utf-8") as f:
                accessions = [line.strip() for line in f if line.strip()]

            print(f"\nProcessing for {gene_type}:")
            for genome in ["cds", "protein"]:
                fasta_files = find_fasta_files(source_folder, genome_type=genome)
                if not fasta_files:
                    print("Error: No FASTA files found in the source folder.", file=sys.stderr)
                    sys.exit(1)

                output_fasta = cfg["nucleotide"] if genome == "cds" else cfg["protein"]
                with open(output_fasta, "w", encoding="utf-8") as out_handle:
                    for query_accession in accessions:
                        print(f"\nProcessing accession: {query_accession}")
                        found = False
                        for fasta in fasta_files:
                            for header, sequence in read_fasta(fasta):
                                if query_accession.split()[0] == header.split()[0]:
                                    write_fasta_record(out_handle, query_accession, sequence)
                                    print(f"\n    Found in {fasta}")
                                    found = True
                                    break
                            if found:
                                break
                        if not found:
                            print(f"Warning: Accession {query_accession} not found.", file=sys.stderr)
                            out_handle.write(f">{query_accession}\n")
                print(f"Extraction complete. Results saved in {output_fasta}.")


if __name__ == "__main__":
    main()
