#!/bin/bash
# ============================================================================
# Module: GTF/GFF Extraction from Reference Annotations
# ============================================================================
# Filters a reference GFF3 annotation to extract entries for HMMER-identified
# genes. Produces:
#   1. A filtered GFF3 (standard format, all features for matched genes)
#   2. A custom GTF compatible with extract_sequences.sh
#
# Custom GTF format (tab-separated):
#   GeneID  TranscriptID  Chromosome  Start  End  Strand  .  Attributes
#
# Usage:
#   bash extract_gtf.sh \
#       --annotation <reference.gff> \
#       --ids <hit_ids.txt> \
#       --outdir <output_dir> \
#       [--label <genome_label>] \
#       [--overwrite]
# ============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../logging/logging_utils.sh"

# ===================== IMPORTANT VARIABLES =====================
ANNOTATION=""
IDS_FILE=""
OUTDIR=""
LABEL="identified"
OVERWRITE=false
# ===============================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --annotation) ANNOTATION="$2"; shift 2 ;;
        --ids)        IDS_FILE="$2"; shift 2 ;;
        --outdir)     OUTDIR="$2"; shift 2 ;;
        --label)      LABEL="$2"; shift 2 ;;
        --overwrite)  OVERWRITE=true; shift ;;
        *) log_error "Unknown option: $1"; exit 1 ;;
    esac
done

[[ -z "$ANNOTATION" ]] && { log_error "--annotation required"; exit 1; }
[[ -z "$IDS_FILE" ]]   && { log_error "--ids required"; exit 1; }
[[ -z "$OUTDIR" ]]     && { log_error "--outdir required"; exit 1; }
[[ ! -f "$ANNOTATION" ]] && { log_error "Annotation file not found: $ANNOTATION"; exit 1; }
[[ ! -f "$IDS_FILE" ]]   && { log_error "IDs file not found: $IDS_FILE"; exit 1; }

mkdir -p "$OUTDIR"

FILTERED_GFF="$OUTDIR/${LABEL}_filtered.gff3"
CUSTOM_GTF="$OUTDIR/${LABEL}.gtf"

if [[ -f "$CUSTOM_GTF" && "$OVERWRITE" != true ]]; then
    log_info "Output exists, skipping: $CUSTOM_GTF"
    exit 0
fi

# ---------------------------------------------------------------------------
# Step 1: Derive gene-level IDs from HMMER hit IDs
#   SMEL4.1_01g000730.1.01  -->  SMEL4.1_01g000730.1   (strip .01 isoform)
#   SMEL5_01g008730.1       -->  SMEL5_01g008730        (strip .1 isoform)
# ---------------------------------------------------------------------------
TMP_GENE_IDS=$(mktemp "${TMPDIR:-/tmp}/gene_ids_XXXXXX.txt")
TMP_HIT_MAP=$(mktemp "${TMPDIR:-/tmp}/hit_map_XXXXXX.txt")
trap "rm -f '$TMP_GENE_IDS' '$TMP_HIT_MAP'" EXIT

while IFS= read -r hit_id; do
    [[ -z "$hit_id" || "$hit_id" =~ ^# ]] && continue
    gene_id="${hit_id%.*}"
    echo "$gene_id" >> "$TMP_GENE_IDS"
    printf '%s\t%s\n' "$gene_id" "$hit_id" >> "$TMP_HIT_MAP"
done < "$IDS_FILE"
sort -u -o "$TMP_GENE_IDS" "$TMP_GENE_IDS"

NUM_GENES=$(wc -l < "$TMP_GENE_IDS")
log_info "Found $NUM_GENES unique gene IDs from HMMER results"

# ---------------------------------------------------------------------------
# Step 2: Filter GFF3 — keep features where ID or Parent matches a gene ID
#   gene   → ID=<gene_id>
#   mRNA   → Parent=<gene_id>
#   CDS/exon → Parent=<mRNA_id>  (parent stripped = gene_id)
# ---------------------------------------------------------------------------
log_info "Filtering annotation: $(basename "$ANNOTATION")"

awk -F'\t' -v GENE_IDS="$TMP_GENE_IDS" '
BEGIN {
    while ((getline id < GENE_IDS) > 0)
        gene_ids[id] = 1
}
/^#/ { print; next }
{
    # Extract ID and Parent from GFF3 attributes (column 9)
    id = ""; parent = ""
    n = split($9, attrs, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^ +| +$/, "", attrs[i])
        if (substr(attrs[i], 1, 3) == "ID=")
            id = substr(attrs[i], 4)
        else if (substr(attrs[i], 1, 7) == "Parent=")
            parent = substr(attrs[i], 8)
    }

    # Direct gene-level match (gene features)
    if (id in gene_ids) { print; next }

    # Direct parent match (mRNA features whose Parent is a gene_id)
    if (parent in gene_ids) { print; next }

    # Child features (CDS/exon) whose Parent is a transcript_id
    # Strip last .suffix to get gene_id
    parent_gene = parent
    sub(/\.[^.]+$/, "", parent_gene)
    if (parent_gene in gene_ids) { print; next }
}' "$ANNOTATION" > "$FILTERED_GFF"

FILTERED_LINES=$(grep -cv '^#' "$FILTERED_GFF" 2>/dev/null || echo "0")
log_info "Filtered GFF3: $FILTERED_LINES feature lines -> $FILTERED_GFF"

# ---------------------------------------------------------------------------
# Step 3: Create custom GTF for extract_sequences.sh
#   Format: GeneID <TAB> TranscriptID <TAB> Chr <TAB> Start <TAB> End <TAB> Strand <TAB> . <TAB> Attributes
#   One row per gene (uses gene-level coordinates for full span)
# ---------------------------------------------------------------------------
log_info "Generating custom GTF for sequence extraction"

awk -F'\t' -v HIT_MAP="$TMP_HIT_MAP" '
BEGIN {
    while ((getline line < HIT_MAP) > 0) {
        split(line, a, "\t")
        hit_ids[a[1]] = a[2]
    }
}
!/^#/ && $3 == "gene" {
    # Extract ID from attributes
    id = ""
    n = split($9, attrs, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^ +| +$/, "", attrs[i])
        if (substr(attrs[i], 1, 3) == "ID=") {
            id = substr(attrs[i], 4)
            break
        }
    }
    if (!(id in hit_ids)) next

    transcript_id = hit_ids[id]
    chr = $1; start = $4; end = $5; strand = $7

    printf "%s\t%s\t%s\t%s\t%s\t%s\t.\t\"transcript_id \"\"%s\"\"; gene_id \"\"%s\"\";\"\n", \
        id, transcript_id, chr, start, end, strand, transcript_id, id
}' "$FILTERED_GFF" > "$CUSTOM_GTF"

GTF_LINES=$(wc -l < "$CUSTOM_GTF")
log_info "Custom GTF: $GTF_LINES genes -> $CUSTOM_GTF"
