#!/bin/bash
# =============================================================================
# Fix: Re-run MAFFT + IQ-TREE2 + RAxML-NG + combined bootstrap visualization
#      for the v4_BLAST_Groups_bitscore200 NUCLEOTIDE merged set after
#      replacing the wrong IbDMP3 genomic sequence with the correct CDS
#      (reverse complement of BSXM01000036.1:complement(3330812-3331456),
#       protein GMD11507.1, "protein DMP6-like").
#
# Run on the HPC (IBSONE) after conda activate egg:
#   bash rerun_IbDMP3_tree.sh
#
# Estimated time: 30-90 min (MAFFT ~1 min, IQ-TREE2 ~20 min, RAxML ~60 min
#                            at 10000 bootstraps each on 12 threads)
# =============================================================================

set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES="$PIPELINE_DIR/modules"
CPU=12

# ---- Paths ------------------------------------------------------------------
GENE_GROUP="DMP"
GENOME="GPE001970_SMEL5"
SET="v4_BLAST_Groups_bitscore200"

MSA_DIR="$PIPELINE_DIR/III_RESULT/$GENE_GROUP/04_MSA/$GENOME/Selected_Result/$SET"
MERGED_NUC_DIR="$MSA_DIR/selected_v4_blast_groups_merged_nucleotide/MAFFT_aligned"
ALIGNED_FAS="$MERGED_NUC_DIR/per_genome_gpe001970_${SET}_NUCLEOTIDE_Sequence.fas"
INPUT_FA="$MERGED_NUC_DIR/input_fasta.fa"

PHYLO_DIR="$PIPELINE_DIR/III_RESULT/$GENE_GROUP/05_Phylogenetics/$GENOME/$SET/MAFFT_aligned"
IQ_PREFIX="$PHYLO_DIR/IQTREE2/per_genome_gpe001970_${SET}_NUCLEOTIDE_Sequence_IQTREE2"
RAX_PREFIX="$PHYLO_DIR/RAXML/per_genome_gpe001970_${SET}_NUCLEOTIDE_Sequence_RAXML"

CONFIG_FILE="$PIPELINE_DIR/05_phyloCONFIG.toml"

# ---- Sanity check -----------------------------------------------------------
[[ -f "$INPUT_FA" ]] || { echo "ERROR: input FASTA not found: $INPUT_FA"; exit 1; }
python3 - <<PYCHECK
import re
with open("$INPUT_FA") as f:
    content = f.read()
entries = re.split(r'(?=^>)', content, flags=re.MULTILINE)
for e in entries:
    if 'IbDMP3' in e:
        lines = e.strip().splitlines()
        seq = "".join(lines[1:])
        assert seq[:3] == "ATG", f"IbDMP3 does not start with ATG: {seq[:10]}"
        assert seq[-3:] == "TAA", f"IbDMP3 does not end with TAA: {seq[-10:]}"
        print(f"[OK] IbDMP3 CDS confirmed: starts ATG, ends TAA, {len(seq)} bp")
        break
else:
    print("[WARN] IbDMP3 not found in merged input - proceeding anyway")
PYCHECK

# =============================================================================
# STEP 1: MAFFT - re-align the merged nucleotide set
# =============================================================================
echo "=== STEP 1: MAFFT alignment ==="
mkdir -p "$MERGED_NUC_DIR"

cleaned_fa="$(mktemp --suffix=.fa)"
awk '/^>/{print; next} {gsub(/[^A-Za-z*\-]/, ""); print}' "$INPUT_FA" > "$cleaned_fa"

mafft --thread "$CPU" --localpair --maxiterate 1000 "$cleaned_fa" > "$ALIGNED_FAS"
rm -f "$cleaned_fa"

# Generate .aln (Clustal format) for compatibility
ALN_FILE="${ALIGNED_FAS%.fas}.aln"
python3 "$MODULES/utils/fasta_to_clustal.py" "$ALIGNED_FAS" "$ALN_FILE" \
    && echo "[OK] .aln generated" || echo "[WARN] .aln conversion failed (non-fatal)"

echo "[OK] MAFFT done: $ALIGNED_FAS"

# =============================================================================
# STEP 2: IQ-TREE2 - rebuild nucleotide tree (GTR+F+R4, 10000 UFBoot + alrt)
# =============================================================================
echo "=== STEP 2: IQ-TREE2 ==="
mkdir -p "$(dirname "$IQ_PREFIX")"

# Deduplicate sequences (pipeline standard)
dedup_fa="$(mktemp --suffix=_dedup.fas)"
python3 - <<PYDEDUP
import re
with open("$ALIGNED_FAS") as f:
    content = f.read()
entries = [e for e in re.split(r'(?=^>)', content, flags=re.MULTILINE) if e.strip()]
seen_seqs, out = set(), []
for e in entries:
    lines = e.strip().splitlines()
    seq = "".join(lines[1:])
    if seq not in seen_seqs:
        seen_seqs.add(seq)
        out.append(e.rstrip())
with open("$dedup_fa", 'w') as f:
    f.write("\n".join(out) + "\n")
print(f"[OK] Dedup: {len(entries)} -> {len(out)} sequences")
PYDEDUP

iqtree -s "$dedup_fa" \
    -m GTR+F+R4 --seqtype DNA \
    -T "$CPU" \
    -bb 10000 -alrt 10000 \
    --allnni --safe --bnni -pers 0.05 \
    --redo \
    -pre "$IQ_PREFIX"

rm -f "$dedup_fa"
echo "[OK] IQ-TREE2 done: ${IQ_PREFIX}.treefile"

# =============================================================================
# STEP 3: RAxML-NG - rebuild nucleotide tree (GTR+FC+R4, 10000 bootstraps)
# =============================================================================
echo "=== STEP 3: RAxML-NG ==="
mkdir -p "$(dirname "$RAX_PREFIX")"

# Reuse same dedup FASTA (create fresh copy)
dedup_fa2="$(mktemp --suffix=_dedup_rax.fas)"
python3 - <<PYDEDUP2
import re
with open("$ALIGNED_FAS") as f:
    content = f.read()
entries = [e for e in re.split(r'(?=^>)', content, flags=re.MULTILINE) if e.strip()]
seen_seqs, out = set(), []
for e in entries:
    lines = e.strip().splitlines()
    seq = "".join(lines[1:])
    if seq not in seen_seqs:
        seen_seqs.add(seq)
        out.append(e.rstrip())
with open("$dedup_fa2", 'w') as f:
    f.write("\n".join(out) + "\n")
PYDEDUP2

raxml-ng \
    --msa "$dedup_fa2" \
    --model GTR+FC+R4 \
    --threads "$CPU" \
    --force perf_threads \
    --seed 12345 \
    --prefix "$RAX_PREFIX" \
    --all \
    --bs-trees 10000 \
    --redo

rm -f "$dedup_fa2"
echo "[OK] RAxML-NG done: ${RAX_PREFIX}.raxml.support"

# =============================================================================
# STEP 4: Combined bootstrap visualization
# =============================================================================
echo "=== STEP 4: Combined bootstrap tree ==="

bash "$MODULES/05_phylogenetic_analysis/combined_bootstrap_tree.sh" \
    --treedir "$PHYLO_DIR" \
    --outdir  "$PHYLO_DIR" \
    --config  "$CONFIG_FILE" \
    --threads "$CPU" \
    --overwrite true

echo "=== All done ==="
echo "Combined bootstrap PNG:"
find "$PHYLO_DIR" -name "*combined_bootstrap*.png" | sort
