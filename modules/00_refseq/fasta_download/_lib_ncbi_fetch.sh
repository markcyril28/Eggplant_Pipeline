#!/bin/bash
# ============================================================================
# NCBI E-utilities fetch library
# ============================================================================
# Sourced by per-species DMP download modules. Provides:
#   ncbi_fetch_by_locus    <locus> <gene_name> <organism>
#   ncbi_fetch_by_symbol   <gene_symbol> <gene_name> <organism>
#   ncbi_fetch_via_gene_db <id> <gene_name> <organism>
#
# ncbi_fetch_by_locus / ncbi_fetch_by_symbol resolve directly in nuccore.
# ncbi_fetch_via_gene_db routes through db=gene + elink for cross-database IDs
# (Phytozome, CottonGen, RAP-DB, CucurBit) not indexed directly in nuccore.
# Output: <gene_name>_<locus_or_symbol>.fasta in PWD.
# Requires: wget, python3 (no extra bioinformatics tooling).
# Set NCBI_API_KEY in env to raise the eutils rate limit.
# ============================================================================

[[ "${LIB_NCBI_FETCH_SOURCED:-}" == "true" ]] && return 0
LIB_NCBI_FETCH_SOURCED="true"

EUTILS_BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EUTILS_DB="${EUTILS_DB:-nuccore}"
EUTILS_RETMAX="${EUTILS_RETMAX:-10}"
EUTILS_TIMEOUT="${EUTILS_TIMEOUT:-30}"   # wget --timeout in seconds (per-attempt)
EUTILS_TRIES="${EUTILS_TRIES:-3}"        # wget --tries
# Default rettype for nuccore efetch:
#   fasta         = full mRNA (5'UTR + CDS + 3'UTR) -- legacy default; produces
#                   internal stops when frame-1 translated
#   fasta_cds_na  = pure CDS nucleotide (no UTRs)   -- preferred for DMP queries
#   fasta_cds_aa  = pure CDS protein                -- canonical RefSeq protein
EUTILS_RETTYPE="${EUTILS_RETTYPE:-fasta_cds_na}"
# Rate-limit-aware default delay between requests within ONE worker:
#   anonymous NCBI limit = 3 req/s -> safe single-worker delay = 0.4s
#   API-key NCBI limit   = 10 req/s -> safe single-worker delay = 0.12s
# Caller may override EUTILS_DELAY explicitly; otherwise we pick from API key.
if [[ -z "${EUTILS_DELAY:-}" ]]; then
    if [[ -n "${NCBI_API_KEY:-}" ]]; then
        EUTILS_DELAY="0.12"
    else
        EUTILS_DELAY="0.4"
    fi
fi

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
    ids=$(wget -qO- --timeout="$EUTILS_TIMEOUT" --tries="$EUTILS_TRIES" "$esearch_url" \
        | _extract_ids 2>/dev/null) || ids=""

    if [[ -z "$ids" ]]; then
        echo "    [WARN] No NCBI ${EUTILS_DB} hits for: $query" >&2
        return 1
    fi

    efetch_url=$(_eutils_url efetch.fcgi "db=${EUTILS_DB}&id=${ids}&rettype=${EUTILS_RETTYPE}&retmode=text")
    if wget -qO "$out_file" --timeout="$EUTILS_TIMEOUT" --tries="$EUTILS_TRIES" "$efetch_url" \
       && [[ -s "$out_file" ]]; then
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

# Internal: given comma-separated NCBI Gene IDs, return linked nuccore IDs (up to 10)
_elink_gene_to_nuccore() {
    local gene_ids="$1"
    local nuccore_ids="" elink_url

    # Prefer RefSeq RNA links (fewer off-target records)
    elink_url=$(_eutils_url elink.fcgi \
        "dbfrom=gene&db=nuccore&id=${gene_ids}&linkname=gene_nuccore_refseqrna&retmode=json")
    nuccore_ids=$(wget -qO- --timeout="$EUTILS_TIMEOUT" --tries="$EUTILS_TRIES" \
        "$elink_url" | python3 -c '
import json,sys
d=json.load(sys.stdin)
ids=[]
for ls in d.get("linksets",[]):
    for ld in ls.get("linksetdbs",[]):
        for l in ld.get("links",[]):
            ids.append(str(l["id"]) if isinstance(l,dict) else str(l))
print(",".join(ids[:10]))
' 2>/dev/null) || nuccore_ids=""

    if [[ -z "$nuccore_ids" ]]; then
        # Fallback: all mRNA links (broader but may include genomic)
        elink_url=$(_eutils_url elink.fcgi \
            "dbfrom=gene&db=nuccore&id=${gene_ids}&linkname=gene_nuccore&retmode=json")
        nuccore_ids=$(wget -qO- --timeout="$EUTILS_TIMEOUT" --tries="$EUTILS_TRIES" \
            "$elink_url" | python3 -c '
import json,sys
d=json.load(sys.stdin)
ids=[]
for ls in d.get("linksets",[]):
    for ld in ls.get("linksetdbs",[]):
        for l in ld.get("links",[]):
            ids.append(str(l["id"]) if isinstance(l,dict) else str(l))
print(",".join(ids[:10]))
' 2>/dev/null) || nuccore_ids=""
    fi

    echo "$nuccore_ids"
    sleep "$EUTILS_DELAY"
}

# Public: search NCBI Gene DB then elink to nuccore.
# Use for cross-database locus IDs (Phytozome, CottonGen, RAP-DB, CucurBit, etc.)
# where the identifier is indexed in the gene DB but not directly in nuccore.
ncbi_fetch_via_gene_db() {
    local id="$1" gene_name="$2" organism="$3"
    local out_file="${gene_name}_${id}.fasta"
    if [[ -s "$out_file" && "${OVERWRITE:-false}" != "true" ]]; then
        echo "    [SKIP] $out_file exists"
        return 0
    fi
    echo ">>> $gene_name  ($id)  in  $organism  [via gene DB + elink]"

    local gene_ids nuccore_ids q_enc esearch_url efetch_url
    q_enc=$(_urlenc "\"${id}\"[All Fields] AND \"${organism}\"[ORGN]")
    esearch_url=$(_eutils_url esearch.fcgi "db=gene&term=${q_enc}&retmax=5&retmode=json")
    gene_ids=$(wget -qO- --timeout="$EUTILS_TIMEOUT" --tries="$EUTILS_TRIES" "$esearch_url" \
        | _extract_ids 2>/dev/null) || gene_ids=""
    sleep "$EUTILS_DELAY"

    if [[ -z "$gene_ids" ]]; then
        echo "    [WARN] No NCBI gene hits for: $id ($organism)" >&2
        return 1
    fi

    nuccore_ids=$(_elink_gene_to_nuccore "$gene_ids")

    if [[ -z "$nuccore_ids" ]]; then
        echo "    [WARN] No nuccore records linked from gene entry: $id ($organism)" >&2
        return 1
    fi

    efetch_url=$(_eutils_url efetch.fcgi "db=nuccore&id=${nuccore_ids}&rettype=${EUTILS_RETTYPE}&retmode=text")
    if wget -qO "$out_file" --timeout="$EUTILS_TIMEOUT" --tries="$EUTILS_TRIES" "$efetch_url" \
       && [[ -s "$out_file" ]]; then
        echo "    -> $out_file"
    else
        echo "    [ERROR] efetch failed for: $id" >&2
        rm -f "$out_file"
        return 1
    fi
    sleep "$EUTILS_DELAY"
}
