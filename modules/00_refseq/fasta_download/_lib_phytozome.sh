#!/bin/bash
# ============================================================================
# Phytozome / JGI Genome Portal sequence fetch library
# ============================================================================
# Sourced by per-species download modules when NCBI does not carry the locus
# identifier (Phytozome-specific IDs such as Glyma.18G, CsaV3, Cla97, Ibatatas,
# Soltu.DM, LOC_Os).
#
# Portal: https://phytozome-next.jgi.doe.gov/
# Requires: curl, python3, gzip
# Env vars (must be set before sourcing):
#   JGI_USER         JGI / Phytozome account login email
#   JGI_PASSWORD     JGI / Phytozome account password
#   JGI_COOKIE_JAR   (optional) path for session cookie; default /tmp/jgi_<PID>
#
# Usage example in a download script:
#   source "$SCRIPT_DIR/_lib_phytozome.sh"
#   phytozome_fetch_gene_cds "Gmax" "GmDMP1" \
#       "Glyma\.18G097400" "GmDMP1_Glyma.18G097400.fasta"
#
# Phytozome JGI portal organism codes (verified empirically — these are the
# *short* names accepted by the get-directory ext-api, not the assembly slugs
# shown in the report URL):
#
#   Species              Portal code   Used for
#   ────────────────────────────────────────────────────────────────────
#   Glycine max          Gmax          Glyma.18G IDs
#   Cucumis sativus      Csativus      CsaV3_ IDs
#   Citrullus lanatus    Clanatus      Cla97 IDs
#   Solanum tuberosum    Stuberosum    Soltu.DM IDs
#   Ipomoea batatas      Ibatatas      Ibatatas IDs
#   Oryza sativa         Osativa       LOC_Os IDs
#   Gossypium hirsutum   Ghirsutum     Gh_ IDs (CottonGen mirror)
#
# The portal returns XML (not JSON), so this lib parses XML via python.
# ============================================================================

[[ "${LIB_PHYTOZOME_SOURCED:-}" == "true" ]] && return 0
LIB_PHYTOZOME_SOURCED="true"

_PHYTO_SIGNON="https://signon.jgi.doe.gov/signon/create"
_PHYTO_FILES_API="https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism="
_PHYTO_DL_BASE="https://genome.jgi.doe.gov"

_SCRIPT_DIR_PHYTO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_PHYTO_EXTRACTOR="${_PHYTO_EXTRACTOR:-$_SCRIPT_DIR_PHYTO/extract_dmp_from_local.py}"

# Returns cookie-jar path on success, empty string on failure.
_phytozome_login() {
    if [[ -z "${JGI_USER:-}" || -z "${JGI_PASSWORD:-}" ]]; then
        echo "    [WARN] JGI_USER / JGI_PASSWORD not set; Phytozome download skipped" >&2
        echo ""
        return 1
    fi
    local jar="${JGI_COOKIE_JAR:-/tmp/jgi_cookies_${$}_${BASHPID}.txt}"
    curl -s --max-time 30 --connect-timeout 10 \
        -c "$jar" \
        --data-urlencode "login=${JGI_USER}" \
        --data-urlencode "password=${JGI_PASSWORD}" \
        "$_PHYTO_SIGNON" -o /dev/null
    if [[ -s "$jar" ]] && grep -q "jgi_session" "$jar"; then
        echo "$jar"
    else
        echo "    [WARN] JGI login failed (check JGI_USER / JGI_PASSWORD)" >&2
        rm -f "$jar"
        echo ""
        return 1
    fi
}

# Return the relative download URL for the primary-transcript CDS file of an organism.
# Picks the highest-version *primaryTranscriptOnly* CDS available; falls back to any CDS.
_phytozome_cds_url() {
    local organism_code="$1" cookie_jar="$2"
    local raw_xml
    raw_xml=$(curl -s --max-time 60 --connect-timeout 10 \
        -b "$cookie_jar" "${_PHYTO_FILES_API}${organism_code}")
    if [[ -z "$raw_xml" ]]; then
        echo "    [WARN] Phytozome files API returned empty response for: $organism_code" >&2
        echo ""
        return 0
    fi
    if [[ "$raw_xml" == "Portal does not exist" ]]; then
        echo "    [WARN] Phytozome portal does not exist for: $organism_code (try the short form, e.g., 'Gmax' not 'Gmax_Wm82.a2.v1')" >&2
        echo ""
        return 0
    fi
    echo "$raw_xml" | python3 -c '
import sys, re
import xml.etree.ElementTree as ET

try:
    root = ET.fromstring(sys.stdin.read())
except ET.ParseError as e:
    print(f"    [WARN] Phytozome API returned non-XML (terms not accepted? wrong org code?): {e}", file=sys.stderr)
    sys.exit(0)

best_primary = ""
best_primary_ver = ""
best_any = ""

for f in root.iter("file"):
    fn = f.attrib.get("filename", "")
    url = f.attrib.get("url", "")
    if not (url and fn.endswith((".fa.gz", ".fasta.gz"))):
        continue
    if "cds" not in fn.lower():
        continue
    # Extract version-like token to prefer newer assemblies
    m = re.search(r"_v?(\d+(?:\.\d+)*)", fn)
    ver = m.group(1) if m else "0"
    if "primaryTranscriptOnly" in fn:
        if ver > best_primary_ver:
            best_primary_ver = ver
            best_primary = url
    elif not best_any:
        best_any = url

print(best_primary or best_any)
' 2>/dev/null
}

# Public: download CDS FASTA from Phytozome, extract sequences matching pattern,
# write to out_file.
#
# phytozome_fetch_gene_cds <organism_code> <gene_name> <pattern> <out_file>
#   organism_code  Phytozome JGI portal code (Gmax, Csativus, Stuberosum, ...)
#   gene_name      label embedded in the extracted FASTA header
#   pattern        pipe-separated regex for FASTA-header matching
#   out_file       destination FASTA (truncated and rewritten on each run)
phytozome_fetch_gene_cds() {
    local organism_code="$1" gene_name="$2" pattern="$3" out_file="$4"
    if [[ -s "$out_file" && "${OVERWRITE:-false}" != "true" ]]; then
        echo "    [SKIP] $out_file exists (Phytozome)"
        return 0
    fi
    echo ">>> $gene_name  from Phytozome  ($organism_code)"

    local cookie_jar
    cookie_jar=$(_phytozome_login) || return 1
    [[ -z "$cookie_jar" ]] && return 1

    local cds_url
    cds_url=$(_phytozome_cds_url "$organism_code" "$cookie_jar")
    if [[ -z "$cds_url" ]]; then
        echo "    [WARN] No CDS file found on Phytozome for: $organism_code" >&2
        [[ -z "${JGI_COOKIE_JAR:-}" ]] && rm -f "$cookie_jar"
        return 1
    fi

    # XML attribute values come URL-encoded with &amp;
    cds_url="${cds_url//&amp;/&}"

    local tmp_gz="/tmp/phyto_cds_${$}_${BASHPID}.fa.gz"
    local tmp_fa="/tmp/phyto_cds_${$}_${BASHPID}.fa"
    curl -sL --max-time 600 -b "$cookie_jar" "${_PHYTO_DL_BASE}${cds_url}" -o "$tmp_gz"
    [[ -z "${JGI_COOKIE_JAR:-}" ]] && rm -f "$cookie_jar"

    if [[ ! -s "$tmp_gz" ]]; then
        echo "    [ERROR] Phytozome CDS download empty for: $organism_code" >&2
        rm -f "$tmp_gz"
        return 1
    fi

    if ! gzip -t "$tmp_gz" 2>/dev/null; then
        echo "    [ERROR] Phytozome CDS download is not gzip (auth/terms issue?) for: $organism_code" >&2
        rm -f "$tmp_gz"
        return 1
    fi

    gzip -dc "$tmp_gz" > "$tmp_fa"
    rm -f "$tmp_gz"

    # Truncate before extract so reruns don't append duplicates.
    : > "$out_file"
    python3 "$_PHYTO_EXTRACTOR" \
        --fasta "$tmp_fa" \
        --patterns "$pattern" \
        --out "$out_file" \
        --name "$gene_name" || true
    rm -f "$tmp_fa"

    if [[ -s "$out_file" ]]; then
        echo "    -> $out_file  (Phytozome)"
    else
        echo "    [WARN] No records matched pattern '$pattern' for $organism_code" >&2
        rm -f "$out_file"
        return 1
    fi
}
