#!/bin/bash
set -euo pipefail

# Assign command-line arguments to variables
input_fasta=$1
header='>'
output_fasta=$3

# Extract the sequence using awk
awk -v header=">$header" '
BEGIN {found = 0}
# Match the header
$0 ~ header {
    found = 1
    print $0 > "'"$output_fasta"'"
    next
}
# Print lines as part of the sequence
found == 1 && $0 !~ /^>/ {
    print $0 > "'"$output_fasta"'"
}
# Stop when another header is found
found == 1 && $0 ~ /^>/ && $0 != header {
    exit
}' "$input_fasta"

# Check if the output file was created successfully
if [[ -s $output_fasta ]]; then
    echo "Sequence extracted to $output_fasta."
else
    echo "Header '$header' not found in $input_fasta or no sequence extracted."
    rm -f "$output_fasta" # Remove empty file
    exit 2
fi
