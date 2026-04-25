#!/bin/bash
# ============================================================================
# Module: Genomic Sequence Extraction
# ============================================================================
# Extracts gene sequences with upstream/downstream flanking regions from
# reference genome using GTF coordinates.
#
# FASTA OUTPUT NAMING CONVENTION:
#   {GeneID}_{Chromosome}_{CDS_Start}_{CDS_End}_{Strand}_genomic_{UPSTREAMbp}up_{DOWNSTREAMbp}down.fa
#   Example: AtGRF1_chr1_15230_17891_+_genomic_1000up_1000down.fa
#
# FASTA HEADER FORMAT:
#   >{GeneID} chr={Chromosome} cds={CDS_Start}-{CDS_End} strand={Strand} region={Chr}:{ExtractStart}-{ExtractEnd} upstream={UPSTREAMbp}bp downstream={DOWNSTREAMbp}bp
#   Example: >AtGRF1 chr=chr1 cds=15230-17891 strand=+ region=chr1:14230-18891 upstream=1000bp downstream=1000bp
#   NOTE: Sequence ID (first token) is the gene name only — PlantCARE requires alphanumeric-only IDs.
#
# Conventions:
#   - Strand: '+' (forward) or '-' (reverse), sequence is rev-comped if '-'
#   - All coordinates are 1-based, inclusive
#   - CDS coordinates represent the original annotated gene boundaries
#   - region coordinates include upstream/downstream extensions
#
# Usage: bash extract_sequences.sh --genome <ref.fa> --gtf <annotations.gtf> --outdir <dir> [--upstream N] [--downstream N] [--threads N]
# ============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../utils/logging.sh"

GENOME_FASTA=""
GTF_FILE=""
OUTPUT_DIR="."
UPSTREAM=1000
DOWNSTREAM=1000
MAX_PARALLEL=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome)     GENOME_FASTA="$2"; shift 2 ;;
        --gtf)        GTF_FILE="$2"; shift 2 ;;
        --outdir)     OUTPUT_DIR="$2"; shift 2 ;;
        --upstream)   UPSTREAM="$2"; shift 2 ;;
        --downstream) DOWNSTREAM="$2"; shift 2 ;;
        --threads)    MAX_PARALLEL="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

[[ -z "$GENOME_FASTA" ]] && { log_error "Missing --genome"; exit 1; }
[[ -z "$GTF_FILE" ]]     && { log_error "Missing --gtf"; exit 1; }

# Organize into a subdirectory by parameter type (e.g. 2000up_0down/)
OUTPUT_DIR="$OUTPUT_DIR/${UPSTREAM}up_${DOWNSTREAM}down"
mkdir -p "$OUTPUT_DIR"

# Ensure genome is indexed
if [[ ! -f "${GENOME_FASTA}.fai" ]]; then
    log_info "Indexing genome..."
    samtools faidx "$GENOME_FASTA"
fi

log_step "Extracting sequences (upstream=$UPSTREAM, downstream=$DOWNSTREAM, parallel=$MAX_PARALLEL)"

# Wait until a parallel slot opens
wait_for_slot() {
    local limit="${1:-$MAX_PARALLEL}"
    while (( $(jobs -rp | wc -l) >= limit )); do
        sleep 0.3
    done
}

# GTF columns: GeneName, Col2, ChromosomeNumber, LowerNumber, UpperNumber, Strand, Rest
PIDS=()
while IFS=$'\t' read -r GeneName Col2 ChromosomeNumber LowerNumber UpperNumber Strand Rest; do
    [[ -z "$GeneName" || "$GeneName" =~ ^# ]] && continue

    wait_for_slot "$MAX_PARALLEL"
    (
        Start=$((LowerNumber - UPSTREAM))
        End=$((UpperNumber + DOWNSTREAM))
        [[ "$Start" -lt 1 ]] && Start=1
        if [[ "$Start" -gt "$End" ]]; then
            tmp=$Start; Start=$End; End=$tmp
        fi

        Region="${ChromosomeNumber}:${Start}-${End}"
        # Filename: {gene}_{chr}_{cds_start}_{cds_end}_{strand}_genomic_{upstream}up_{downstream}down.fa
        OutFile="$OUTPUT_DIR/${GeneName}_${ChromosomeNumber}_${LowerNumber}_${UpperNumber}_${Strand}_genomic_${UPSTREAM}up_${DOWNSTREAM}down.fa"

        # Extract sequence to temp file
        samtools faidx "$GENOME_FASTA" "$Region" > "${OutFile}.tmp"

        # Reverse complement for negative strand
        if [[ "$Strand" == "-" ]]; then
            seqtk seq -r "${OutFile}.tmp" > "${OutFile}.revcomp"
            mv "${OutFile}.revcomp" "${OutFile}.tmp"
        fi

        # Reformat header with gene name and metadata.
        # Sequence ID (first token) = gene name only — PlantCARE requires alphanumeric IDs.
        # Coordinates and strand are kept as description fields after the first space.
        # Format: >GeneID chr=CHR cds=START-END strand=STRAND region=CHR:start-end upstream=Nbp downstream=Mbp
        awk -v gene="$GeneName" -v chr="$ChromosomeNumber" -v cds_start="$LowerNumber" \
            -v cds_end="$UpperNumber" -v strand="$Strand" -v start="$Start" \
            -v end="$End" -v up="$UPSTREAM" -v down="$DOWNSTREAM" \
            'NR==1 {
                print ">" gene " chr=" chr " cds=" cds_start "-" cds_end " strand=" strand " region=" chr ":" start "-" end " upstream=" up "bp downstream=" down "bp"
                next
            }
            {
                gsub(/[[:space:]]/, "", $0)
                while (length($0) > 80) {
                    print substr($0, 1, 80)
                    $0 = substr($0, 81)
                }
                if (length($0) > 0) {
                    print $0
                }
            }' "${OutFile}.tmp" > "$OutFile"
        rm "${OutFile}.tmp"

        log_info "Extracted: $GeneName -> $OutFile"
    ) &
    PIDS+=($!)
done < "$GTF_FILE"

# Wait for all and collect failures
FAILED=0
for pid in "${PIDS[@]}"; do
    wait "$pid" || ((FAILED++))
done
(( FAILED > 0 )) && { log_error "$FAILED gene extraction(s) failed"; exit 1; }

log_info "Sequence extraction complete -> $OUTPUT_DIR"
