#!/bin/bash
set -euo pipefail
# NOTE: This is a legacy script. The active pipeline uses 2_BLAST_Alignment.sh + modules/01_identification/blastn.sh

# ===================== IMPORTANT VARIABLES =====================
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
REFSEQS_DIR="$PIPELINE_DIR/1_RefSeqs/a_Smel_RefSeqs"

Eggplant_V4_1_genome="$REFSEQS_DIR/Solanum_melongena_v4.1/Eggplant_V4.1.fa"
Eggplant_V4_1_transcripts="$REFSEQS_DIR/Solanum_melongena_v4.1/Eggplant_V4.1_transcripts.function.fa"
# ===============================================================
query_sequence="Gene_Control_Query.fasta"  # Query sequence file
e="${2:-1e-10}"      # E-value threshold (default: 1e-10)
ws="${3:-11}"        # Word size (default: 11)

# Check if input files exist
if [ ! -f "$Eggplant_V4_1_transcripts" ]; then
    echo "Error: Reference transcript file not found: $Eggplant_V4_1_transcripts"
    exit 1
fi
if [ ! -f "$query_sequence" ]; then
    echo "Error: Query sequence file not found: $query_sequence"
    exit 1
fi

output_folder="002_Gene_Control_BLASTn_Search_Results"
rm -rf "$output_folder"
mkdir -p "$output_folder"

echo "=============================== Creating BLAST DBs ==============================="
# Create BLAST DBs (skip if already exist)
mkdir -p "DB"
gbase=$(basename "$Eggplant_V4_1_genome")
gbase="${gbase%%.*}"
dbpath="DB/${gbase}_DB"
if [ -f "${dbpath}.nsq" ]; then
    echo "BLAST DB for $gbase exists. Skipping."
    echo ""
else
    echo "Creating BLAST DB for $gbase..."
    makeblastdb -in "$Eggplant_V4_1_transcripts" -dbtype nucl -out "$dbpath"
fi

echo "=============================== Running BLASTn ==============================="


query_base=$(basename "$query_sequence")
query_base="${query_base%%.*}"

csv="$output_folder/blastn_${gbase}_VS_${query_base}_${e}_${ws}.csv"
txt="$output_folder/blastn_${gbase}_VS_${query_base}_${e}_${ws}.txt"

# Write CSV header
echo "Subject ID,Query ID,E-value,Percent Identity,Alignment Length,Query Length,Target Length,Query Coverage,Mismatches,Gaps,Bit Score,Max Score,Total Score" > "$csv"

# Run BLASTn and output results
blastn -query "$query_sequence" -db "$dbpath" -evalue "$e" -word_size "$ws" \
    -outfmt "10 sseqid qseqid evalue pident length qlen slen qcovs mismatch gaps bitscore score score" >> "$csv"

blastn -query "$query_sequence" -db "$dbpath" -evalue "$e" -word_size "$ws" > "$txt"
echo "BLASTn search completed. Results saved in $output_folder."
