#!/bin/bash
# Download BjuDMP1-4 CDS for Brassica juncea (Chu et al., 2025).
#
# Source paper: Chu et al. (2025) Horticulture Research, 12(7):uhaf094.
#   "In vivo maternal haploid induction in Brassica juncea". HIR 0.64-1.51%.
#
# The Chu paper does not publish per-paralog locus IDs (the supplementary
# kdocs.cn link cited in the paper hosts no downloadable sequence file at
# verification time). Per-paralog identities resolved here by k-mer overlap of
# AtDMP8 (UniProt O80493) + AtDMP9 (NP_198781.1) against the Ensembl Plants
# B. juncea proteome (assembly ASM1870372v1, T84-66 var. tumida; release 62).
# The four highest-scoring proteins are clear orthologs and span both AABB
# subgenomes (2 BjuA + 2 BjuB), matching the Chu nomenclature of four BjuDMPs:
#
#   shared k-mers (k=7) vs AtDMP8/9 union | gene id          | chr | candidate
#   ------------------------------------- | ---------------- | --- | ---------
#   191 / 238   (78.3% normalised)        | BjuA04g10430S    | A04 | BjuDMP1
#   153 / 238   (63.2%)                   | BjuA03g54090S    | A03 | BjuDMP2
#   147 / 238   (60.7%)                   | BjuB08g57390S    | B08 | BjuDMP3
#   120 / 238   (43.5%)                   | BjuB01g27600S    | B01 | BjuDMP4
#   --- sharp drop to 2 shared k-mers in next hits ---
#
# Source: Ensembl Plants REST API (https://rest.ensembl.org/sequence/id/<id>)
# returns the canonical CDS in FASTA format.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REST="https://rest.ensembl.org"
GENE_IDS=(
    "BjuA04g10430S"
    "BjuA03g54090S"
    "BjuB08g57390S"
    "BjuB01g27600S"
)
GENE_NAMES=(
    "BjuDMP1_AA_subgenome_A04"
    "BjuDMP2_AA_subgenome_A03"
    "BjuDMP3_BB_subgenome_B08"
    "BjuDMP4_BB_subgenome_B01"
)

# Prefer wget (pipeline default on Linux/WSL); fall back to curl.
_fetch() {
    local url="$1" out="$2"
    if command -v wget >/dev/null 2>&1; then
        wget -qO "$out" --header="Accept: text/x-fasta" --timeout=30 --tries=3 "$url"
    elif command -v curl >/dev/null 2>&1; then
        curl -sf --max-time 30 --retry 3 -H "Accept: text/x-fasta" -o "$out" "$url"
    else
        echo "    [ERROR] need wget or curl in PATH" >&2
        return 1
    fi
}

for ((i=0; i<${#GENE_IDS[@]}; i++)); do
    gid="${GENE_IDS[i]}"
    name="${GENE_NAMES[i]}"
    out_file="${name}_${gid}.fasta"
    if [[ -s "$out_file" && "${OVERWRITE:-false}" != "true" ]]; then
        echo "    [SKIP] $out_file exists"
        continue
    fi
    url="${REST}/sequence/id/${gid}?type=cds"
    echo ">>> $name  ($gid, Ensembl Plants)"
    if _fetch "$url" "$out_file" && [[ -s "$out_file" ]] && head -1 "$out_file" | grep -q "^>"; then
        echo "    -> $out_file ($(wc -c < "$out_file") bytes)"
    else
        echo "    [WARN] Ensembl REST fetch failed for $gid" >&2
        rm -f "$out_file"
    fi
    sleep 0.3
done
