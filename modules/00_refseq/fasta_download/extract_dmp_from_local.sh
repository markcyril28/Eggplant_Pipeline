#!/bin/bash
# ============================================================================
# Read II_INPUTS/DMP_HI_registry.tsv and, for each row with source_mode=local,
# extract the matching FASTA records from the configured local genome file
# into the species output directory.
#
# Args:
#   $1  Registry TSV path (default: $PROJECT_ROOT/II_INPUTS/DMP_HI_registry.tsv)
#   $2  Output base dir   (default: $PROJECT_ROOT/I_RefSeqs/d_DMP_Query_Fasta)
#
# OVERWRITE=true forces re-extraction even if the per-species *_local.fasta
# already exists.
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
REGISTRY="${1:-$PIPELINE_DIR/II_INPUTS/DMP_HI_registry.tsv}"
OUT_BASE="${2:-$PIPELINE_DIR/I_RefSeqs/d_DMP_Query_Fasta}"
EXTRACTOR="$SCRIPT_DIR/extract_dmp_from_local.py"
OVERWRITE="${OVERWRITE:-false}"

if [[ ! -f "$REGISTRY" ]]; then
    echo "[ERROR] registry not found: $REGISTRY" >&2
    exit 1
fi

local_count=0
ncbi_count=0
extracted_count=0

while IFS=$'\t' read -r species_full species_short gene_name locus_id hir_pct citation \
                       output_dir source_mode local_cds local_prot patterns; do
    # Skip comments and header
    [[ "$species_full" =~ ^# ]] && continue
    [[ "$species_full" == "species_full" ]] && continue
    [[ -z "${species_full:-}" ]] && continue

    if [[ "$source_mode" != "local" ]]; then
        ncbi_count=$((ncbi_count + 1))
        continue
    fi
    local_count=$((local_count + 1))

    species_out="$OUT_BASE/$output_dir"
    mkdir -p "$species_out"
    out_file="$species_out/${species_short}DMPs_local_extract.fasta"

    if [[ "$OVERWRITE" != "true" && -s "$out_file" ]]; then
        echo "  [SKIP] $species_short ($species_full) -> $(basename "$out_file") exists"
        extracted_count=$((extracted_count + 1))
        continue
    fi

    # Truncate so re-runs do not double-append
    : > "$out_file"

    # Try CDS first
    if [[ -n "$local_cds" && "$local_cds" != "-" ]]; then
        cds_path="$PIPELINE_DIR/$local_cds"
        if [[ -f "$cds_path" ]]; then
            python3 "$EXTRACTOR" --fasta "$cds_path" --patterns "$patterns" \
                --out "$out_file" --name "${gene_name//|/+}" || true
        else
            echo "  [WARN] local CDS not found: $cds_path" >&2
        fi
    fi
    # Then protein (appended to same out_file)
    if [[ -n "$local_prot" && "$local_prot" != "-" ]]; then
        prot_path="$PIPELINE_DIR/$local_prot"
        if [[ -f "$prot_path" ]]; then
            python3 "$EXTRACTOR" --fasta "$prot_path" --patterns "$patterns" \
                --out "$out_file" --name "${gene_name//|/+}" || true
        else
            echo "  [WARN] local protein not found: $prot_path" >&2
        fi
    fi

    if [[ -s "$out_file" ]]; then
        extracted_count=$((extracted_count + 1))
        echo "  [OK]   $species_short ($species_full) -> $(basename "$out_file")"
    else
        echo "  [MISS] $species_short ($species_full): no records matched patterns; will fall back to NCBI"
        rm -f "$out_file"
    fi
done < "$REGISTRY"

echo ""
echo "Summary: $extracted_count/$local_count local extractions succeeded; $ncbi_count rows deferred to NCBI"
