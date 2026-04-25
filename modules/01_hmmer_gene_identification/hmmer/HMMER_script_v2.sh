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

CPU=12   # Number of threads to use

echo "----------------------------------------------------------"
echo "Download"
echo "----------------------------------------------------------"

#wget https://www.ebi.ac.uk/interpro/wwwapi/protein/reviewed/entry/pfam/PF08879/ || \
#    echo "[WARNING] Download may fail; verify URL manually."

# -----------------------------
# Working Directories
# -----------------------------
#rm -rf "3_PREPROCESSING" "4_HMMER_RESULTS" "2_CONFIGS"
#mkdir -p "2_CONFIGS" "3_PREPROCESSING" "4_HMMER_RESULTS"

HMM_PROFILES=(
    #"1_INPUTS/a_HMM/PF08879_WRC_Profile.hmm"
    #"1_INPUTS/a_HMM/PF08880_QLQ_Profile.hmm"
    #"1_INPUTS/a_HMM/PF05030_SSXT_Profile.hmm" # or SNH domain
    "1_INPUTS/a_HMM/PF05078_DMP.hmm"
)

ALIGNMENTS=(
    #"1_INPUTS/b_Alignments/PF08879_WRC_clustal2.alignment.seed"
    #"1_INPUTS/b_Alignments/PF08880_QLQ_clustal2.alignment.seed"
    #"1_INPUTS/b_Alignments/PF05030_SSXT_clustal2.alignment.seed"
    # PF05078_DMP - pre-built HMM, no alignment needed
)

#e_value=1e-15 #1e-5 #1e-15 ##1e-10

e_value_LIST=(
    1e-5 
    #1e-10 
    #1e-15
    #1e-25
)

MAIN_OUTDIR="2_MAIN_RESULTS_per_evalue"
mkdir -p "$MAIN_OUTDIR"
rm -rf "${MAIN_OUTDIR:?}/"*

: << 'OFF'
HMMER_OUTDIR="2_MAIN_RESULTS_per_evalue/a_HMMER_RESULTS"
rm -rf "$HMMER_OUTDIR"
mkdir -p "$HMMER_OUTDIR"

RAW_EXTRACTED_FOLDER="2_MAIN_RESULTS_per_evalue/b_RAW_EXTRACTED"
rm -rf "$RAW_EXTRACTED_FOLDER"
mkdir -p "$RAW_EXTRACTED_FOLDER"

CD_HIT_OUTDIR="2_MAIN_RESULTS_per_evalue/c_CD_HIT_Reduced"
rm -rf "$CD_HIT_OUTDIR"
mkdir -p "$CD_HIT_OUTDIR"
OFF

for e_value in "${e_value_LIST[@]}"; do
    mkdir -p "$MAIN_OUTDIR/e-value_$e_value"
    HMMER_OUTDIR="$MAIN_OUTDIR/e-value_$e_value/a_HMMER_RESULTS"
    RAW_EXTRACTED_FOLDER="$MAIN_OUTDIR/e-value_$e_value/b_RAW_EXTRACTED"
    CD_HIT_OUTDIR="$MAIN_OUTDIR/e-value_$e_value/c_CD_HIT_Reduced"
    rm -rf "${HMMER_OUTDIR:?}" "${RAW_EXTRACTED_FOLDER:?}" "${CD_HIT_OUTDIR:?}"
    mkdir -p "$HMMER_OUTDIR"
    mkdir -p "$RAW_EXTRACTED_FOLDER"
    mkdir -p "$CD_HIT_OUTDIR"

    echo "--------------------------------------------------------------------------------------------------------------------"
    echo "[INFO] Starting HMMER workflow with E-value threshold: $e_value"
    echo "--------------------------------------------------------------------------------------------------------------------"
    for i in "${!HMM_PROFILES[@]}"; do
        base=$(basename "${HMM_PROFILES[$i]}" _Profile.hmm)
        # Handle DMP profile which doesn't have _Profile suffix
        if [[ "$base" == "${HMM_PROFILES[$i]##*/}" ]]; then
            base=$(basename "${HMM_PROFILES[$i]}" .hmm)
        fi

        # Only run hmmbuild if alignment file exists and is not empty
        if [[ -n "${ALIGNMENTS[$i]}" && -f "${ALIGNMENTS[$i]}" ]]; then
            echo "[INFO] Building HMM profile for $base..."
            hmmbuild "${HMM_PROFILES[$i]}" "${ALIGNMENTS[$i]}"
        else
            echo "[INFO] Using pre-built HMM profile for $base (no alignment provided)..."
        fi

        echo "[INFO] Pressing HMM profile for $base..."
        hmmpress -f "${HMM_PROFILES[$i]}"
        mkdir -p "$HMMER_OUTDIR/$base"

        echo "[INFO] Running hmmsearch for $base..."
        hmmsearch --cpu "$CPU" \
            --tblout "$HMMER_OUTDIR/$base/${base}_hits.tbl" \
            --domtblout "$HMMER_OUTDIR/$base/${base}_dom.tbl" \
            "${HMM_PROFILES[$i]}" \
            "$EGGPLANT_V4_1_PROTEINS" > "$HMMER_OUTDIR/hmmsearch_${base}.log"

        echo "[INFO] Filtering hmmsearch domain hits for $base (E-value <= $e_value)..."
        awk -v ev="$e_value" '$1 !~ /^#/ && $7 <= ev {print $0}' \
            "$HMMER_OUTDIR/$base/${base}_dom.tbl" > "$HMMER_OUTDIR/$base/${base}_hits_filtered.domtbl"
        
        echo "[INFO] Done with $base."

        echo "[INFO] Extracting matched protein and transcript sequences..."
        hmm_result="$HMMER_OUTDIR/$base/${base}_hits_filtered.domtbl"
        if [[ -s "$hmm_result" ]]; then
            while IFS= read -r first_col; do
                echo "[INFO] Extracting protein sequences for ID:       $first_col"
                awk -v id="$first_col" '
                $0 ~ "^>" && $0 ~ id {print_header=1; print $0; next}
                print_header && $0 !~ "^>" {gsub(/\.$/, "", $0); print $0}
                $0 ~ "^>" && $0 !~ id {
                    if (print_header) print "";  # Add blank line after sequence
                    print_header=0
                }
                ' "$EGGPLANT_V4_1_PROTEINS" >> "$RAW_EXTRACTED_FOLDER/${base}_matched_proteins.fa"

                echo "[INFO] Extracting transcript sequence for ID:     $first_col"
                awk -v id="$first_col" '
                $0 ~ "^>" && $0 ~ id {print_header=1; print $0; next}
                print_header && $0 !~ "^>" {gsub(/\.$/, "", $0); print $0}
                $0 ~ "^>" && $0 !~ id {
                    if (print_header) print "";  # Add blank line after sequence
                    print_header=0
                }
                ' "$EGGPLANT_V4_1_TRANSCRIPTS" >> "$RAW_EXTRACTED_FOLDER/${base}_matched_transcripts.fa"
            done < <(awk '{print $1}' "$hmm_result" | sort | uniq)
        
            echo "[INFO] Reducing redundancy with CD-HIT for $base..."
            mkdir -p "$CD_HIT_OUTDIR/$base"
            cd-hit -i "$RAW_EXTRACTED_FOLDER/${base}_matched_proteins.fa" \
                -o "$CD_HIT_OUTDIR/$base/${base}_matched_proteins_cdhit_reduced.fa" \
                -c 0.9 -n 5 -T "$CPU" -M 16000
            cd-hit -i "$RAW_EXTRACTED_FOLDER/${base}_matched_transcripts.fa" \
                -o "$CD_HIT_OUTDIR/$base/${base}_matched_transcripts_cdhit_reduced.fa" \
                -c 0.9 -n 5 -T "$CPU" -M 16000
        fi
    done
    # -----------------------------
    # NEW: QLQ + WRC Intersection
    # -----------------------------
    echo "[INFO] Identifying proteins with BOTH QLQ and WRC domains..."

    WRC_DOMTBL="$HMMER_OUTDIR/PF08879_WRC/PF08879_WRC_hits_filtered.domtbl"
    QLQ_DOMTBL="$HMMER_OUTDIR/PF08880_QLQ/PF08880_QLQ_hits_filtered.domtbl"

    if [[ ! -f "$WRC_DOMTBL" || ! -f "$QLQ_DOMTBL" ]]; then
        echo "[WARN] QLQ+WRC intersection skipped: domtbl missing (WRC: $WRC_DOMTBL, QLQ: $QLQ_DOMTBL)"
        echo "[WARN]   Ensure PF08879_WRC and PF08880_QLQ are in HMM_PROFILES and HMMER ran successfully."
    else
    WRC_IDS="$HMMER_OUTDIR/PF08879_WRC/PF08879_WRC_hits.ids"
    QLQ_IDS="$HMMER_OUTDIR/PF08880_QLQ/PF08880_QLQ_hits.ids"
    INTERSECT_IDS="$HMMER_OUTDIR/QLQ_WRC_common.ids"

    awk '{print $1}' "$WRC_DOMTBL" | sort -u > "$WRC_IDS"
    awk '{print $1}' "$QLQ_DOMTBL" | sort -u > "$QLQ_IDS"
    comm -12 "$QLQ_IDS" "$WRC_IDS" > "$INTERSECT_IDS"

    RAW_EXTRACTED_COMBINED_PROTEIN="$RAW_EXTRACTED_FOLDER/PF08879_PF08880_QLQ_WRC_matched_proteins.fa"
    RAW_EXTRACTED_COMBINED_TRANSCRIPT="$RAW_EXTRACTED_FOLDER/PF08879_PF08880_QLQ_WRC_matched_transcripts.fa"
    > "$RAW_EXTRACTED_COMBINED_PROTEIN"
    > "$RAW_EXTRACTED_COMBINED_TRANSCRIPT"

    while read -r id; do
        echo "[INFO] Extracting protein sequence for ID:        $id"
        awk -v id="$id" '
        $0 ~ "^>" && $0 ~ id {print_header=1; print $0; next}
        print_header && $0 !~ "^>" {gsub(/\.$/, "", $0); print $0}
        $0 ~ "^>" && $0 !~ id {
            if (print_header) print "";
            print_header=0
        }
        ' "$EGGPLANT_V4_1_PROTEINS" >> "$RAW_EXTRACTED_COMBINED_PROTEIN"

        echo "[INFO] Extracting transcript sequence for ID:     $id"
        awk -v id="$id" '
        $0 ~ "^>" && $0 ~ id {print_header=1; print $0; next}
        print_header && $0 !~ "^>" {gsub(/\.$/, "", $0); print $0}
        $0 ~ "^>" && $0 !~ id {
            if (print_header) print "";  # Add blank line after sequence
            print_header=0
        }
        ' "$EGGPLANT_V4_1_TRANSCRIPTS" >> "$RAW_EXTRACTED_COMBINED_TRANSCRIPT"
    done < "$INTERSECT_IDS"

    echo "[INFO] Reducing redundancy with CD-HIT for QLQ+WRC intersection..."
    cd-hit -i "$RAW_EXTRACTED_COMBINED_PROTEIN" \
        -o "$CD_HIT_OUTDIR/QLQ_WRC_matched_proteins_cdhit_reduced.fa" \
        -c 0.9 -n 5 -T "$CPU" -M 16000

    cd-hit -i "$RAW_EXTRACTED_COMBINED_TRANSCRIPT" \
        -o "$CD_HIT_OUTDIR/QLQ_WRC_matched_transcripts_cdhit_reduced.fa" \
        -c 0.9 -n 5 -T "$CPU" -M 16000
    fi  # end QLQ+WRC intersection guard
done


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
