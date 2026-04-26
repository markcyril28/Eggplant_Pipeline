#!/bin/bash
# Download BnaDMP 4 paralogs for Brassica napus (Y. Li et al., 2022).
#
# Source paper: Y. Li et al. (2022) Plant Biotechnology J., 20:2052-2063.
#   "An in planta haploid induction system in Brassica napus". HIR up to 2.53%.
#
# The four cited Bna* loci are Darmor-bzh v4.1 IDs (Phytozome/BnPIR style),
# not registered in NCBI nuccore. Searching NCBI by these symbols falls
# through to whole-chromosome WGS records (>40 MB each), so we resolve them
# via Ensembl Plants which hosts the same Darmor-bzh assembly
# (AST_PRJEB5043_v1) with the Bna* gene_symbol attribute. Mapping verified
# 2026-04-26 via k-mer overlap of AtDMP8 + AtDMP9 against the Ensembl B. napus
# proteome (top 4 hits == the four cited paralogs, with sharp cliff to rank 5):
#
#   gene_symbol       | Ensembl gene id          | shared k-mers vs AtDMP8/9
#   ----------------  | ------------------------ | ------------------------
#   BnaA04g09480D     | GSBRNA2T00114309001      | 184/238 (77.5%)
#   BnaA03g55920D     | GSBRNA2T00019729001      | 153/238 (63.2%)
#   BnaC03g03890D     | GSBRNA2T00134587001      | 145/238 (59.9%)
#   BnaC04g31700D     | GSBRNA2T00044898001      |  66/238 (49.6%, partial)
#   --- next hit drops to 2 shared k-mers ---
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REST="https://rest.ensembl.org"
GENE_IDS=(
    "GSBRNA2T00019729001"
    "GSBRNA2T00134587001"
    "GSBRNA2T00114309001"
    "GSBRNA2T00044898001"
)
GENE_NAMES=(
    "BnaDMP_A03_BnaA03g55920D"
    "BnaDMP_C03_BnaC03g03890D"
    "BnaDMP_A04_BnaA04g09480D"
    "BnaDMP_C04_BnaC04g31700D"
)

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
    out_file="${name}.fasta"
    if [[ -s "$out_file" && "${OVERWRITE:-false}" != "true" ]]; then
        echo "    [SKIP] $out_file exists"
        continue
    fi
    url="${REST}/sequence/id/${gid}?type=cds"
    echo ">>> $name  ($gid, Ensembl Plants Darmor-bzh)"
    if _fetch "$url" "$out_file" && [[ -s "$out_file" ]] && head -1 "$out_file" | grep -q "^>"; then
        echo "    -> $out_file ($(wc -c < "$out_file") bytes)"
    else
        echo "    [WARN] Ensembl REST fetch failed for $gid" >&2
        rm -f "$out_file"
    fi
    sleep 0.3
done
