#!/bin/bash
# ============================================================================
# NCBI E-utilities fetch library
# ============================================================================
# Sourced by per-species DMP download modules. Provides:
#   ncbi_fetch_by_locus  <locus> <gene_name> <organism>
#   ncbi_fetch_by_symbol <gene_symbol> <gene_name> <organism>
#
# Both resolve a query to nuccore UIDs via esearch, then pull FASTA via efetch.
# Output: <gene_name>_<locus_or_symbol>.fasta in PWD.
# Requires: wget, python3 (no extra bioinformatics tooling).
# Set NCBI_API_KEY in env to raise the eutils rate limit.
# ============================================================================

[[ "${LIB_NCBI_FETCH_SOURCED:-}" == "true" ]] && return 0
LIB_NCBI_FETCH_SOURCED="true"

EUTILS_BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EUTILS_DB="${EUTILS_DB:-nuccore}"
EUTILS_RETMAX="${EUTILS_RETMAX:-10}"
EUTILS_DELAY="${EUTILS_DELAY:-0.4}"   # seconds between requests; >=0.34 keeps under 3 req/s

_urlenc() {
    python3 -c 'import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1]))' "$1"
}

_extract_ids() {
    python3 -c 'import json,sys; d=json.load(sys.stdin); print(",".join(d["esearchresult"].get("idlist",[])))'
}

_eutils_url() {
    local endpoint="$1" qs="$2"
    if [[ -n "${NCBI_API_KEY:-}" ]]; then
        printf '%s/%s?%s&api_key=%s' "$EUTILS_BASE" "$endpoint" "$qs" "$NCBI_API_KEY"
    else
        printf '%s/%s?%s' "$EUTILS_BASE" "$endpoint" "$qs"
    fi
}

# Internal: query → write FASTA
_fetch_query() {
    local query="$1" out_file="$2"
    local q_enc ids esearch_url efetch_url
    q_enc=$(_urlenc "$query")
    esearch_url=$(_eutils_url esearch.fcgi "db=${EUTILS_DB}&term=${q_enc}&retmax=${EUTILS_RETMAX}&retmode=json")
    ids=$(wget -qO- "$esearch_url" | _extract_ids 2>/dev/null) || ids=""

    if [[ -z "$ids" ]]; then
        echo "    [WARN] No NCBI ${EUTILS_DB} hits for: $query" >&2
        return 1
    fi

    efetch_url=$(_eutils_url efetch.fcgi "db=${EUTILS_DB}&id=${ids}&rettype=fasta&retmode=text")
    if wget -qO "$out_file" "$efetch_url" && [[ -s "$out_file" ]]; then
        echo "    -> $out_file"
    else
        echo "    [ERROR] efetch failed for: $query" >&2
        rm -f "$out_file"
        return 1
    fi
    sleep "$EUTILS_DELAY"
}

# Public: fetch by exact locus identifier (e.g. Solyc05g007920)
ncbi_fetch_by_locus() {
    local locus="$1" gene_name="$2" organism="$3"
    local out_file="${gene_name}_${locus}.fasta"
    if [[ -s "$out_file" && "${OVERWRITE:-false}" != "true" ]]; then
        echo "    [SKIP] $out_file exists"
        return 0
    fi
    echo ">>> $gene_name  ($locus)  in  $organism"
    _fetch_query "${locus} AND \"${organism}\"[ORGN]" "$out_file" || return 0
}

# Public: fetch by gene symbol when no per-paralog locus is published
ncbi_fetch_by_symbol() {
    local symbol="$1" gene_name="$2" organism="$3"
    local out_file="${gene_name}_${symbol}.fasta"
    if [[ -s "$out_file" && "${OVERWRITE:-false}" != "true" ]]; then
        echo "    [SKIP] $out_file exists"
        return 0
    fi
    echo ">>> $gene_name  ($symbol)  in  $organism"
    _fetch_query "${symbol}[Gene Name] AND \"${organism}\"[ORGN]" "$out_file" || return 0
}
