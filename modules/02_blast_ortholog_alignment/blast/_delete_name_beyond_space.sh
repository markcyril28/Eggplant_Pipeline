#!/bin/bash
set -euo pipefail

for file in *.fasta *.fa *.fas; do
    if [[ -f "$file" ]]; then
        sed -i 's/^\(>[^ ]*\) .*/\1/' "$file"
    fi
done