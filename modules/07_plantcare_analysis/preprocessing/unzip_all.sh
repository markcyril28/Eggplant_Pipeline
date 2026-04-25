#!/bin/bash
set -euo pipefail

# Script to extract archive files in the current directory
# Supports: .zip, .tar.gz, .tgz
# - Output folder name is the same as the archive filename (without extension)
# - Skips extraction if output folder already exists
# - Keeps the original archive file

echo "Starting to process archive files..."

# Look for these archive patterns
archives=( *.zip *.tar.gz *.tgz )

found=false
for archive in "${archives[@]}"; do
    # If pattern didn't match any file, it will remain literal; skip those
    if [ ! -f "$archive" ]; then
        continue
    fi

    found=true

    # Determine output folder name by stripping the archive extension
    case "$archive" in
        *.zip)
            folder_name="${archive%.zip}"
            ;;
        *.tar.gz)
            folder_name="${archive%.tar.gz}"
            ;;
        *.tgz)
            folder_name="${archive%.tgz}"
            ;;
        *)
            folder_name="${archive%.*}"
            ;;
    esac

    # Skip extraction if the target folder already exists
    if [ -d "$folder_name" ]; then
        echo "Skipping '$archive' - folder '$folder_name' already exists."
        continue
    fi

    echo "Extracting '$archive' to '$folder_name'..."
    # Create the output folder and extract into it
    mkdir -p "$folder_name"

    if [[ "$archive" == *.zip ]]; then
        unzip -q "$archive" -d "$folder_name"
    else
        # For tar.gz and tgz
        tar -xzf "$archive" -C "$folder_name"
    fi

    if [ $? -eq 0 ]; then
        echo "Successfully extracted '$archive'"
    else
        echo "Error extracting '$archive'"
    fi
done

if [ "$found" = false ]; then
    echo "No archive files found in the current directory."
fi

echo "Done processing all archive files."
