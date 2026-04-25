#!/usr/bin/env python3
"""Fix malformed FASTA files where > headers appear mid-line
or header text is concatenated with sequence data.

Usage:
    python fix_fasta_format.py input.fa [output.fa]
    python fix_fasta_format.py input.fa --inplace

If output is omitted, prints to stdout. Use --inplace to overwrite.
"""
import re
import sys


def fix_fasta(input_path, output_path=None):
    with open(input_path, "r") as f:
        raw = f.read()

    # Step 1: Ensure every '>' starts on a new line
    fixed = re.sub(r"(?<!\n)(?<!^)>", "\n>", raw)

    lines = fixed.split("\n")
    entries = []
    current_header = None
    current_seq_lines = []

    for line in lines:
        line = line.rstrip()
        if not line:
            continue

        if line.startswith(">"):
            # Save previous entry
            if current_header is not None:
                entries.append((current_header, current_seq_lines))
            current_header = line
            current_seq_lines = []
        elif current_header is not None and not current_seq_lines:
            # First line after header — could be sequence or header continuation
            if re.match(r"^[ACGTNacgtn]+$", line):
                current_seq_lines.append(line)
            else:
                # Contains non-DNA chars → header continuation or mixed
                m = re.search(r"([ACGTacgt]{15,})$", line)
                if m:
                    header_part = line[: m.start()]
                    seq_part = line[m.start() :]
                    current_header += header_part
                    current_seq_lines.append(seq_part)
                else:
                    # Pure header continuation (no long DNA suffix)
                    current_header += line
        else:
            # Normal sequence line
            current_seq_lines.append(line)

    # Last entry
    if current_header is not None:
        entries.append((current_header, current_seq_lines))

    # Build output
    out_lines = []
    for header, seqs in entries:
        out_lines.append(header)
        out_lines.extend(seqs)

    output = "\n".join(out_lines) + "\n"

    if output_path:
        with open(output_path, "w") as f:
            f.write(output)
    else:
        sys.stdout.write(output)

    return len(entries)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    inp = sys.argv[1]
    if len(sys.argv) > 2:
        if sys.argv[2] == "--inplace":
            n = fix_fasta(inp, inp)
        else:
            n = fix_fasta(inp, sys.argv[2])
    else:
        n = fix_fasta(inp)

    print(f"Fixed {n} entries", file=sys.stderr)
