#!/bin/bash
# Download IpDMP sequences for Ipomoea batatas (sweet potato).
# Strategy 1: NCBI Gene DB + elink (gene symbol search).
# Strategy 2: Phytozome Ibatatas_v1_1 CDS with DMP/DUF679 pattern (requires JGI credentials).
#
# Set JGI_USER and JGI_PASSWORD in the environment to enable Phytozome fallback.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib_ncbi_fetch.sh"
source "$SCRIPT_DIR/_lib_phytozome.sh"

ORGANISM="Ipomoea batatas"
PHYTOZOME_ORGANISM="Ibatatas"
# Pattern matches any DMP / DUF679 gene in the Phytozome CDS FASTA
PHYTOZOME_PATTERN="DMP|DUF679|dmp|duf679"

# NCBI Gene DB search: try known/expected gene symbols for sweet potato DMP
SYMBOL_LIST=("IbDMP" "IbDMP1" "IbDMP2" "IbDMP3")
GENE_NAME_LIST=("IbDMP" "IbDMP1" "IbDMP2" "IbDMP3")

ncbi_success=false
for ((i=0; i<${#SYMBOL_LIST[@]}; i++)); do
    if ncbi_fetch_via_gene_db "${SYMBOL_LIST[i]}" "${GENE_NAME_LIST[i]}" "$ORGANISM"; then
        ncbi_success=true
    fi
done

# Phytozome fallback: download primary-transcript CDS and extract by DMP/DUF679 pattern
if [[ "$ncbi_success" == "false" ]]; then
    echo "    [INFO] NCBI gene DB returned no hits; trying Phytozome ($PHYTOZOME_ORGANISM)"
    phytozome_fetch_gene_cds \
        "$PHYTOZOME_ORGANISM" \
        "IbDMPs" \
        "$PHYTOZOME_PATTERN" \
        "IbDMPs_phytozome.fasta" || true
fi
