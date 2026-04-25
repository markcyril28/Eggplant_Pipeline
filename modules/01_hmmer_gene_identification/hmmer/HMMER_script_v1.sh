#!/bin/bash
set -euo pipefail

# ================================================================
# Profile-HMM Based Gene Identification Workflow using HMMER
# ================================================================
# Dependencies: cd-hit, mafft, trimal, hmmer, wget, gunzip
# Ensure these are installed in your conda environment before running
# ================================================================

# -----------------------------
# (Optional) Install required tools
# -----------------------------
# conda install -y -c bioconda cd-hit mafft trimal hmmer

# NOTE: This is a legacy script. The active pipeline uses 1_HMMER_identify.sh + modules/01_identification/hmmer.sh

# ===================== IMPORTANT VARIABLES =====================
readonly PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
readonly REFSEQS_DIR="$PIPELINE_DIR/1_RefSeqs/a_Smel_RefSeqs"

readonly EGGPLANT_V4_1_GENOME="$REFSEQS_DIR/Solanum_melongena_v4.1/Eggplant_V4.1.fa"
readonly EGGPLANT_V4_1_TRANSCRIPTS="$REFSEQS_DIR/Solanum_melongena_v4.1/Eggplant_V4.1_transcripts.function.fa"
readonly EGGPLANT_V4_1_PROTEINS="$REFSEQS_DIR/Solanum_melongena_v4.1/Eggplant_V4.1_protein.function.fa"
# ===============================================================
# tail -n 10 "$EGGPLANT_V4_1_PROTEINS"

CPU=4   # Number of threads to use

echo "----------------------------------------------------------"
echo "Download"
echo "----------------------------------------------------------"

#wget https://www.ebi.ac.uk/interpro/wwwapi/protein/reviewed/entry/pfam/PF08879/ || \
#    echo "[WARNING] Download may fail; verify URL manually."

# -----------------------------
# Working Directories
# -----------------------------
#rm -rf "3_PREPROCESSING" "4_HMMER_RESULTS"
mkdir -p "2_CONFIGS" "3_PREPROCESSING" "4_HMMER_RESULTS"

SEED_HMM="3_PREPROCESSING/PF08879_profile_from_seed.hmm"        # Custom profile HMM built from seeds
PROFILE_HMM="1_INPUTS/a_Raw/PF08879_profile.hmm"                      # Pre-existing profile HMM from Pfam
ALIGNMENT="1_INPUTS/a_Raw/PF08879.alignment.seed"                     # Curated seed protein sequences (Stockholm format from Pfam)
# NOTE: cd-hit expects FASTA input. If ALIGNMENT is Stockholm format, convert to FASTA first.
# Point SEED_FASTA to the FASTA version of the seed sequences (may need a separate export step).
SEED_FASTA="$ALIGNMENT"  # TODO: verify this file is in FASTA format for cd-hit compatibility

CD_HIT_REDUCED_FASTA="3_PREPROCESSING/seeds_nr_cdhit.faa"
ALIGNMENT_FILE="3_PREPROCESSING/seed_alignment.faa"
TRIMMED_FILE="3_PREPROCESSING/seed_alignment_trimmed.faa"
HMMER_OUTDIR="4_HMMER_RESULTS"

hmmpress "$PROFILE_HMM"


echo "----------------------------------------------------------"
echo "Pfam HMM Database Setup"
echo "----------------------------------------------------------"
echo "[INFO] Downloading Pfam-A HMM database..."
if [ ! -f 2_CONFIGS/Pfam-A.hmm ]; then
    wget -O 2_CONFIGS/Pfam-A.hmm.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip -f 2_CONFIGS/Pfam-A.hmm.gz
else
    echo "[INFO] Pfam-A.hmm already exists, skipping download."
fi

echo "[INFO] Previewing Pfam-A HMM entries for PF08879, PF08880, PF05030..."
# grep -A2 "PF08879" 2_CONFIGS/Pfam-A.hmm | head
# grep -A2 "PF08880" 2_CONFIGS/Pfam-A.hmm | head
# grep -A2 "PF05030" 2_CONFIGS/Pfam-A.hmm | head

hmmpress 2_CONFIGS/Pfam-A.hmm


echo "[INFO] Indexing Pfam-A HMM database..."
if [ -f 2_CONFIGS/Pfam-A.hmm.h3m ]; then
    echo "[INFO] Pfam-A.hmm already indexed, skipping hmmpress."
else
    hmmpress 2_CONFIGS/Pfam-A.hmm
fi
# IGNORE


echo "[INFO] Extracting specific HMM models (PF08879, PF08880, PF05030) from Pfam-A database using hmmfetch..."
mkdir -p "2_CONFIGS/hmm_fetch/"
for HMM_ID in PF08879 PF08880 PF05030; do
    HMM_FILE="2_CONFIGS/hmm_fetch/${HMM_ID}.hmm"
    if [ -f "$HMM_FILE" ]; then
        echo "[INFO] $HMM_FILE already exists, skipping fetch."
    else
        echo "[INFO] Fetching HMM model $HMM_ID..."
        hmmfetch 2_CONFIGS/Pfam-A.hmm "$HMM_ID" > "$HMM_FILE"
    fi
done


echo "[INFO] Extracting specific HMM models (PF08879, PF08880, PF05030) from Pfam-A database using wget..."
mkdir -p "2_CONFIGS/downloaded/"
for HMM_ID in PF08879 PF08880 PF05030; do
    HMM_FILE="2_CONFIGS/downloaded/${HMM_ID}.hmm"
    if [ -f "$HMM_FILE" ]; then
        echo "[INFO] $HMM_FILE already exists, skipping download."
    else
        echo "[INFO] Downloading $HMM_ID HMM from InterPro..."
        wget "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/${HMM_ID}/hmmer" -O "$HMM_FILE"
    fi
done

# Assign Pfam HMM model variables
QLQ_HMM="2_CONFIGS/downloaded/PF08879.hmm"   # QLQ domain (GRFs)
WRC_HMM="2_CONFIGS/downloaded/PF08880.hmm"   # WRC domain (GRFs)
SSXT_HMM="2_CONFIGS/downloaded/PF05030.hmm"  # SSXT domain (GIFs)

echo "----------------------------------------------------------"
echo  "Step 1: CD-HIT Redundancy Reduction"
echo "----------------------------------------------------------"
cd-hit -i "$SEED_FASTA" \
       -o "$CD_HIT_REDUCED_FASTA" \
       -c 0.95 -n 5

echo "----------------------------------------------------------"
echo  "Step 2: Multiple Sequence Alignment and Trimming"
echo "----------------------------------------------------------"
mafft --maxiterate 1000 --localpair \
      "$CD_HIT_REDUCED_FASTA" > "$ALIGNMENT_FILE"

trimal -in "$ALIGNMENT_FILE" \
       -out "$TRIMMED_FILE" -automated1


echo "----------------------------------------------------------"
echo  "Step 3: Build Custom Profile-HMM and Search Proteome"
echo "----------------------------------------------------------"
#hmmbuild "$PROFILE_HMM" "$TRIMMED_FILE"

rm -rf "${HMMER_OUTDIR:?}"/*
hmmbuild "$SEED_HMM" "$ALIGNMENT"
hmmsearch --cpu "$CPU" \
          --tblout "$HMMER_OUTDIR/seed_hits.tbl" \
          --domtblout "$HMMER_OUTDIR/seed_dom.tbl" \
          "$PROFILE_HMM" \
          "$EGGPLANT_V4_1_PROTEINS" > "$HMMER_OUTDIR/hmmsearch_seed.log"

echo "----------------------------------------------------------"
echo  "Step 4: Post-processing"
echo "----------------------------------------------------------"
awk '$1 !~ /^#/ && $7 <= 1e-5 {print $0}' \
    "$HMMER_OUTDIR/seed_dom.tbl" > "$HMMER_OUTDIR/seed_hits_filtered.domtbl"

: << 'DEBUG'
# -----------------------------
# Optional: Run domain-specific searches (Pfam HMMs)
# -----------------------------
for DOMAIN_HMM in "$QLQ_HMM" "$WRC_HMM" "$SSXT_HMM"; do
    DOMAIN_NAME=$(basename "$DOMAIN_HMM" .hmm)
    echo "[INFO] Running hmmsearch with Pfam model: $DOMAIN_NAME"
    hmmsearch --cpu "$CPU" \
              --tblout "$HMMER_OUTDIR/${DOMAIN_NAME}_hits.tbl" \
              --domtblout "$HMMER_OUTDIR/${DOMAIN_NAME}_dom.tbl" \
              "$DOMAIN_HMM" \
              "$EGGPLANT_V4_1_PROTEINS" > "$HMMER_OUTDIR/hmmsearch_${DOMAIN_NAME}.log"
done

DEBUG
# -----------------------------
# Done
# -----------------------------
echo "[INFO] Workflow complete!"
echo "Results (custom HMM):"
echo "- Log:              $HMMER_OUTDIR/hmmsearch_seed.log"
echo "- Raw hits:         $HMMER_OUTDIR/seed_hits.tbl"
echo "- Domain results:   $HMMER_OUTDIR/seed_dom.tbl"
echo "- Filtered hits:    $HMMER_OUTDIR/seed_hits_filtered.domtbl"
echo
echo "Results (Pfam domains: QLQ, WRC, SSXT) are in $HMMER_OUTDIR/"
