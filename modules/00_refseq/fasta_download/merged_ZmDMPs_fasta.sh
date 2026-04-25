#!/bin/bash
set -euo pipefail

# Output file to store the combined fasta sequences
touch ZmDMPs_v2_merged_fasta.fa

output_file="ZmDMPs_v2_merged_fasta.fa"

# Initialize (or clear) the output file
> "$output_file"

# Loop through all fasta files in the directory
for fasta_file in "$PWD"/*.fasta; do
  if [ -f "$fasta_file" ]; then
    cat "$fasta_file" >> "$output_file"
    echo >> "$output_file"  
    echo "Added: $fasta_file"
  fi
done

echo "All fasta files have been merged into $output_file."
