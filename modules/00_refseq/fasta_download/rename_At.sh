#!/bin/bash
set -euo pipefail


# ===================== IMPORTANT VARIABLES =====================
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

# Strings to Remove
string_to_remove="A_thaliana_"
directory="$PIPELINE_DIR/2_INPUTS/DMP/query_fasta/Arabidopsis_thaliana"
# ===============================================================

# Loop through all files in the directory
for file in "$directory"/*; do
  # Extract the file name from the full path
  filename=$(basename "$file")
  echo "$filename"

  # Remove the string from the filename
  new_filename="${filename//$string_to_remove/}"
  
  # Rename the file only if the new name is different
  if [ "$filename" != "$new_filename" ]; then
    mv "$file" "$directory/$new_filename"
    echo "Renamed: $filename -> $new_filename"
  fi
done

echo "Renaming complete."
