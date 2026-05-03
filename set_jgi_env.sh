#!/bin/bash
# ============================================================================
# set_jgi_env.sh  --  Configure JGI/Phytozome credentials and run Stage 0
# ============================================================================
# Exports JGI_USER and JGI_PASSWORD so the Phytozome download fallback in
# Stage 0 (0_RefSeq_setup_and_download.sh) can authenticate with JGI.
#
# Credential lookup order (first match wins):
#   1. Already exported in the calling shell  ($JGI_USER / $JGI_PASSWORD)
#   2. .jgi_credentials file in the pipeline root  (git-ignored; see template)
#   3. Interactive prompt (password is hidden)
#
# Usage — two modes:
#
#   a) SOURCE into your current shell, then run the pipeline separately:
#        source set_jgi_env.sh
#        bash 0_RefSeq_setup_and_download.sh
#
#   b) Run directly (auto-sources creds, then launches Stage 0):
#        bash set_jgi_env.sh
#        bash set_jgi_env.sh --only DOWNLOAD_DMP_HI_GENES   # pass-through args
#        bash set_jgi_env.sh --dry-run
#
# Register your JGI account at: https://signon.jgi.doe.gov/signon/create
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CREDS_FILE="$SCRIPT_DIR/.jgi_credentials"

# ── Step 1: load from file if present ────────────────────────────────────────
if [[ -f "$CREDS_FILE" ]]; then
    # File must contain:  JGI_USER="..."  and  JGI_PASSWORD="..."
    # (same syntax as bash variable assignments, one per line)
    # shellcheck disable=SC1090
    source "$CREDS_FILE"
fi

# ── Step 2: prompt for anything still missing ─────────────────────────────────
if [[ -z "${JGI_USER:-}" ]]; then
    read -r -p "JGI email: " JGI_USER
fi

if [[ -z "${JGI_PASSWORD:-}" ]]; then
    read -r -s -p "JGI password: " JGI_PASSWORD
    echo   # newline after hidden input
fi

# ── Step 3: export ────────────────────────────────────────────────────────────
export JGI_USER
export JGI_PASSWORD

echo "[set_jgi_env] JGI_USER  = $JGI_USER"
echo "[set_jgi_env] JGI_PASSWORD = (set, ${#JGI_PASSWORD} chars)"

# ── Step 4: if executed (not sourced), launch Stage 0 ─────────────────────────
# When sourced, BASH_SOURCE[0] == the sourcing script; skip auto-launch.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "[set_jgi_env] Launching: bash 0_RefSeq_setup_and_download.sh $*"
    bash "$SCRIPT_DIR/0_RefSeq_setup_and_download.sh" "$@"
fi
