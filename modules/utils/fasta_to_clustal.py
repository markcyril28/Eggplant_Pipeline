#!/usr/bin/env python3
"""Convert aligned FASTA (.fas) to Clustal (.aln) format."""
import sys


def fasta_to_clustal(fasta_path, aln_path, line_width=60):
    sequences = {}
    order = []
    current = None
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                # Handle duplicate IDs (e.g. merged cross-species files)
                if name in sequences:
                    n = 2
                    while f"{name}_{n}" in sequences:
                        n += 1
                    name = f"{name}_{n}"
                current = name
                order.append(current)
                sequences[current] = []
            elif current:
                sequences[current].append(line)

    for name in order:
        sequences[name] = "".join(sequences[name])

    lengths = set(len(sequences[n]) for n in order)
    if len(lengths) != 1:
        print(f"Error: sequences have different lengths: {lengths}", file=sys.stderr)
        sys.exit(1)

    seq_len = lengths.pop()
    pad = max(max(len(n) for n in order) + 4, 16)

    with open(aln_path, "w") as out:
        out.write("CLUSTAL W (1.83) multiple sequence alignment\n\n\n")
        for start in range(0, seq_len, line_width):
            end = min(start + line_width, seq_len)
            for name in order:
                out.write(f"{name:<{pad}}{sequences[name][start:end]}\n")
            cons = []
            for i in range(start, end):
                col = {sequences[n][i] for n in order} - {"-"}
                cons.append("*" if len(col) == 1 and col else " ")
            out.write(f"{'':<{pad}}{''.join(cons)}\n\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.fas> <output.aln>", file=sys.stderr)
        sys.exit(1)
    fasta_to_clustal(sys.argv[1], sys.argv[2])
