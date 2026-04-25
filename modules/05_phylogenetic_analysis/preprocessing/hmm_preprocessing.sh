#!/bin/bash
set -euo pipefail

RAW_HMM_DIR="a_RAW_HMM"
PRE_PROCESSED_HMM_DIR="b_PRE_PROCESSING_for_T_COFFEE"
ALIGNMENT_DIR="c_ALIGNMENT"
CPU=4

# Create directories if they don't exist
mkdir -p "$RAW_HMM_DIR" "$PRE_PROCESSED_HMM_DIR" "$ALIGNMENT_DIR"

# Check if required HMM files exist
if [[ ! -f "$RAW_HMM_DIR/PF05030_SSXT_Profile.hmm" ]]; then
    echo "Error: $RAW_HMM_DIR/PF05030_SSXT_Profile.hmm not found"
    exit 1
fi

if [[ ! -f "$RAW_HMM_DIR/PF08879_WRC_Profile.hmm" ]]; then
    echo "Error: $RAW_HMM_DIR/PF08879_WRC_Profile.hmm not found"
    exit 1
fi

if [[ ! -f "$RAW_HMM_DIR/PF08880_QLQ_Profile.hmm" ]]; then
    echo "Error: $RAW_HMM_DIR/PF08880_QLQ_Profile.hmm not found"
    exit 1
fi

# FOR SSXT DOMAIN - Generate consensus alignment from HMM
echo "Generating SSXT domain alignment..."
hmmemit -c "$RAW_HMM_DIR/PF05030_SSXT_Profile.hmm" > "$ALIGNMENT_DIR/PF05030_SSXT_Profile_alignment.aln"

# FOR WRC and QLQ DOMAINS - Generate consensus sequences from HMMs
echo "Generating WRC domain consensus..."
hmmemit -c "$RAW_HMM_DIR/PF08879_WRC_Profile.hmm" > "$PRE_PROCESSED_HMM_DIR/GRF_PF08879_WRC_Profile_consensus.fasta"

echo "Generating QLQ domain consensus..."
hmmemit -c "$RAW_HMM_DIR/PF08880_QLQ_Profile.hmm" > "$PRE_PROCESSED_HMM_DIR/GRF_PF08880_QLQ_Profile_consensus.fasta"

# Combine the consensus sequences by concatenating the sequences (not just files)
echo "Creating combined consensus sequence..."
python3 << 'EOF'
import os

# Read QLQ consensus
with open("b_PRE_PROCESSING_for_T_COFFEE/GRF_PF08880_QLQ_Profile_consensus.fasta", "r") as f:
    qlq_lines = f.readlines()
qlq_seq = "".join(line.strip() for line in qlq_lines[1:])  # Skip header

# Read WRC consensus  
with open("b_PRE_PROCESSING_for_T_COFFEE/GRF_PF08879_WRC_Profile_consensus.fasta", "r") as f:
    wrc_lines = f.readlines()
wrc_seq = "".join(line.strip() for line in wrc_lines[1:])  # Skip header

# Create combined sequence with a linker (you may want to adjust this)
linker = "GGGGG"  # 5 glycines as linker, adjust as needed
combined_seq = qlq_seq + linker + wrc_seq

# Write combined sequence
with open("b_PRE_PROCESSING_for_T_COFFEE/GRF_PF08880_QLQ_PF08879_WRC_combined.fasta", "w") as f:
    f.write(">GRF_QLQ_WRC_combined\n")
    f.write(combined_seq + "\n")

print(f"QLQ length: {len(qlq_seq)}")
print(f"WRC length: {len(wrc_seq)}")
print(f"Combined length: {len(combined_seq)} (including {len(linker)} linker residues)")
EOF

# Convert to Stockholm format for hmmbuild (single sequence alignment)
echo "Converting to Stockholm format..."
python3 << 'EOF'
# Read the combined sequence
with open("b_PRE_PROCESSING_for_T_COFFEE/GRF_PF08880_QLQ_PF08879_WRC_combined.fasta", "r") as f:
    lines = f.readlines()
    sequence = "".join(line.strip() for line in lines[1:])  # Skip header, join all sequence lines

# Write proper Stockholm format
with open("b_PRE_PROCESSING_for_T_COFFEE/GRF_PF08880_QLQ_PF08879_WRC_alignment.sto", "w") as f:
    f.write("# STOCKHOLM 1.0\n")
    f.write(f"GRF_QLQ_WRC_combined    {sequence}\n")
    f.write("//\n")

print(f"Stockholm file created with sequence length: {len(sequence)}")
EOF

# Build a combined HMM profile from the combined alignment
echo "Building combined HMM profile..."
hmmbuild "$PRE_PROCESSED_HMM_DIR/GRF_PF08880_QLQ_PF08879_WRC_alignment.hmm" "$PRE_PROCESSED_HMM_DIR/GRF_PF08880_QLQ_PF08879_WRC_alignment.sto"

# Convert the combined HMM profile to a multiple sequence alignment in FASTA format
echo "Generating final combined alignment..."
hmmemit -c "$PRE_PROCESSED_HMM_DIR/GRF_PF08880_QLQ_PF08879_WRC_alignment.hmm" > "$ALIGNMENT_DIR/GRF_PF08880_QLQ_PF08879_WRC_alignment.aln"

echo ""
echo "=== PREPROCESSING COMPLETE ==="
echo "Generated alignment files:"
echo "SSXT domain: $(realpath "$ALIGNMENT_DIR/PF05030_SSXT_Profile_alignment.aln")"
echo "Combined QLQ+WRC: $(realpath "$ALIGNMENT_DIR/GRF_PF08880_QLQ_PF08879_WRC_alignment.aln")"
echo ""

# REFERENCE:
#t_coffee -seq GRF_candidates.fasta -profile GRF_QLQ_WRC.aln -outfile GRF_QLQ_WRC_aligned.fasta -output=fasta_aln
