#!/bin/bash
set -euo pipefail
# Add Major Group column to tab files based on function mapping

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAPPING_FILE="Major_Groups_to_Function_mapping.csv"

# Load mapping into associative array
declare -A function_to_group
while IFS=',' read -r function major_group; do
    [[ "$function" == "Function" ]] && continue
    # Trim whitespace and carriage returns
    function=$(echo "$function" | sed 's/\r$//' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    major_group=$(echo "$major_group" | sed 's/\r$//' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    function_to_group["$function"]="$major_group"
done < "$MAPPING_FILE"

# Process each .tab file
for tab_file in "$SCRIPT_DIR"/*.tab; do
    [[ ! -f "$tab_file" ]] && continue
    
    temp_file="${tab_file}.tmp"
    
    while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 rest; do
        # Trim whitespace from col8
        col8=$(echo "$col8" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        
        # Map function (col8) to Major Group only if col8 is not empty
        # This will overwrite col9 if it exists
        major_group=""
        if [[ -n "$col8" ]]; then
            major_group="${function_to_group[$col8]}"
        fi
        
        # Output with new/overwritten column 9
        echo -e "${col1}\t${col2}\t${col3}\t${col4}\t${col5}\t${col6}\t${col7}\t${col8}\t${major_group}"
    done < "$tab_file" > "$temp_file"
    
    mv "$temp_file" "$tab_file"
    echo "Processed: $(basename "$tab_file")"
done

echo "Done."
