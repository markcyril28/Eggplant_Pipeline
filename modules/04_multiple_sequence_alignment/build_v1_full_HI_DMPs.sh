#!/bin/bash
# ============================================================================
# Build pre-curated MSA Selected/v1_Full HI DMP input FASTAs
# ============================================================================
# Concatenates the eggplant per-genome HMMER hits (CD-HIT reduced) with the
# 14 verified haploid-inducer DMP merged FASTAs from
# I_RefSeqs/d_DMP_Query_Fasta/<species>/*_merged*.fa  and writes:
#
#   <out_dir>/per_genome_<egg>_and_HI_DMPs_NUCLEOTIDE_Sequence.fasta
#   <out_dir>/per_genome_<egg>_and_HI_DMPs_AMINO_ACID_Sequence.fasta
#
# Defaults to the Selected/v1_Full curated location for the DMP gene group.
# Idempotent: re-running with overwrite=true regenerates; otherwise skips.
# ============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DMP_QUERY_DIR="${DMP_QUERY_DIR:-$PIPELINE_DIR/I_RefSeqs/d_DMP_Query_Fasta}"
EGG_GENOME="${EGG_GENOME:-GPE001970_SMEL5}"
EGG_TAG="${EGG_TAG:-gpe001970}"
EVAL_TAG="${EVAL_TAG:-e-value_1e-5}"
PFAM_TAG="${PFAM_TAG:-PF05078_DMP}"
OUT_DIR="${OUT_DIR:-$PIPELINE_DIR/III_RESULT/DMP/04_MSA/merged_input/Selected/v1_Full}"
OVERWRITE="${OVERWRITE:-true}"

EGG_NUC="$PIPELINE_DIR/III_RESULT/DMP/01_Identification/$EGG_GENOME/$EVAL_TAG/c_CD_HIT_Reduced/$PFAM_TAG/${PFAM_TAG}_transcripts_cdhit.fa"
EGG_AA="$PIPELINE_DIR/III_RESULT/DMP/01_Identification/$EGG_GENOME/$EVAL_TAG/c_CD_HIT_Reduced/$PFAM_TAG/${PFAM_TAG}_proteins_cdhit.fa"

OUT_NUC="$OUT_DIR/per_genome_${EGG_TAG}_and_HI_DMPs_NUCLEOTIDE_Sequence.fasta"
OUT_AA="$OUT_DIR/per_genome_${EGG_TAG}_and_HI_DMPs_AMINO_ACID_Sequence.fasta"

# Index-aligned with [dmp_hi_genes] in 00_RefSeq_setup_and_downloadCONFIG.toml
HI_DIRS=(Arabidopsis_thaliana Zea_mays Solanum_lycopersicum Solanum_tuberosum Nicotiana_tabacum \
         Brassica_napus Brassica_oleracea Brassica_juncea Citrullus_lanatus Cucumis_sativus \
         Medicago_truncatula Glycine_max Gossypium_hirsutum Oryza_sativa)
HI_NUCS=(AtDMPs_merged_fasta.fa ZmDMPs_v2_merged_fasta.fa SlDMPs_merged.fa StDMPs_merged.fa NtDMPs_merged.fa \
         BnaDMPs_merged.fa BoDMPs_merged.fa BjuDMPs_merged.fa ClDMPs_merged.fa CsDMPs_merged.fa \
         MtDMPs_merged.fa GmDMPs_merged.fa GhDMPs_merged.fa OsDMPs_merged.fa)
HI_AAS=(AtDMPs_merged_protein.fa ZmDMPs_v2_merged_protein.fa SlDMPs_merged_protein.fa StDMPs_merged_protein.fa \
        NtDMPs_merged_protein.fa BnaDMPs_merged_protein.fa BoDMPs_merged_protein.fa BjuDMPs_merged_protein.fa \
        ClDMPs_merged_protein.fa CsDMPs_merged_protein.fa MtDMPs_merged_protein.fa GmDMPs_merged_protein.fa \
        GhDMPs_merged_protein.fa OsDMPs_merged_protein.fa)

mkdir -p "$OUT_DIR"

_concat_with_egg() {
    local egg_src="$1" out="$2"; shift 2
    if [[ -s "$out" && "$OVERWRITE" != "true" ]]; then
        echo "[SKIP-EXISTS] $(basename "$out") -- set OVERWRITE=true to regenerate"
        return 0
    fi
    : > "$out"
    if [[ -s "$egg_src" ]]; then
        cat "$egg_src" >> "$out"
        printf '\n' >> "$out"
    else
        echo "[WARN] eggplant source missing: $egg_src" >&2
    fi
    local added=0
    for path in "$@"; do
        if [[ -s "$path" ]]; then
            cat "$path" >> "$out"
            printf '\n' >> "$out"
            added=$((added + 1))
        else
            echo "[WARN] HI source missing: $path" >&2
        fi
    done
    # Strip blank-line clutter and dedup by header (first word) using stdlib python
    python3 - "$out" "$out.tmp" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
seen = set(); kept = 0; dup = 0
out_lines = []
hdr, seq = None, []
def flush():
    global kept, dup
    if hdr is None: return
    key = hdr[1:].split()[0]
    if key in seen:
        dup += 1
        return
    seen.add(key)
    out_lines.append(hdr)
    out_lines.extend(seq)
    kept += 1
with open(src) as f:
    for line in f:
        line = line.rstrip()
        if not line: continue
        if line.startswith(">"):
            flush()
            hdr, seq = line, []
        else:
            seq.append(line)
    flush()
with open(dst, "w") as f:
    f.write("\n".join(out_lines) + "\n")
print(f"  -> {dst}  kept={kept}  dup={dup}", flush=True)
PY
    mv "$out.tmp" "$out"
    echo "  built: $(basename "$out")  ($(grep -c '^>' "$out") seqs, $(wc -c < "$out") bytes; +$added HI species)"
}

NUC_INPUTS=()
for i in "${!HI_DIRS[@]}"; do
    NUC_INPUTS+=("$DMP_QUERY_DIR/${HI_DIRS[i]}/${HI_NUCS[i]}")
done
AA_INPUTS=()
for i in "${!HI_DIRS[@]}"; do
    AA_INPUTS+=("$DMP_QUERY_DIR/${HI_DIRS[i]}/${HI_AAS[i]}")
done

echo "=== build_v1_full_HI_DMPs ==="
echo "  PIPELINE_DIR : $PIPELINE_DIR"
echo "  DMP_QUERY_DIR: $DMP_QUERY_DIR"
echo "  EGG_GENOME   : $EGG_GENOME ($EVAL_TAG / $PFAM_TAG)"
echo "  OUT_DIR      : $OUT_DIR"
echo "  OVERWRITE    : $OVERWRITE"
echo ""
_concat_with_egg "$EGG_NUC" "$OUT_NUC" "${NUC_INPUTS[@]}"
_concat_with_egg "$EGG_AA"  "$OUT_AA"  "${AA_INPUTS[@]}"
