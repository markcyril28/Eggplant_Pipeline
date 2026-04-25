#!/bin/bash
set -euo pipefail

# Script to unzip all zip files in the current directory
# - Output folder name is the same as the zip filename (without .zip extension)
# - Skips unzipping if output folder already exists
# - Keeps the original zip file

echo "Starting to process zip files..."

# Find all zip files in the current directory
for zipfile in *.zip; do
    # Check if any zip files exist
    if [ ! -f "$zipfile" ]; then
        echo "No zip files found in the current directory."
        break
    fi
    
    # Get the filename without the .zip extension
    folder_name="${zipfile%.zip}"
    
    # Check if the output folder already exists
    if [ -d "$folder_name" ]; then
        echo "Skipping '$zipfile' - folder '$folder_name' already exists."
    else
        echo "Unzipping '$zipfile' to '$folder_name'..."
        unzip -q "$zipfile" -d "$folder_name"
        
        if [ $? -eq 0 ]; then
            echo "Successfully unzipped '$zipfile'"
        else
            echo "Error unzipping '$zipfile'"
        fi
    fi
done

echo "Done processing all zip files."
