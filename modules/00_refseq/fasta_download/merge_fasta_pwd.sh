#!/bin/bash
# ============================================================================
# Generic merge: concatenate all *.fasta files in PWD into one .fa
# Usage: bash merge_fasta_pwd.sh <output_filename.fa>
# Example: bash merge_fasta_pwd.sh SlDMPs_merged_fasta.fa
# ============================================================================
set -euo pipefail

OUTPUT_FILE="${1:-merged_fasta.fa}"

# Skip the merge output itself when globbing
shopt -s nullglob
> "$OUTPUT_FILE"
count=0
for fasta_file in "$PWD"/*.fasta; do
    [[ -f "$fasta_file" ]] || continue
    [[ "$(basename "$fasta_file")" == "$OUTPUT_FILE" ]] && continue
    cat "$fasta_file" >> "$OUTPUT_FILE"
    echo >> "$OUTPUT_FILE"
    echo "  added: $(basename "$fasta_file")"
    count=$((count + 1))
done

if (( count == 0 )); then
    echo "  [WARN] No *.fasta files found in $PWD" >&2
    rm -f "$OUTPUT_FILE"
    exit 0
fi

echo "Merged $count file(s) into $OUTPUT_FILE"
