#!/usr/bin/env bash
set -euo pipefail
#######################################################################
# TOML Configuration Parser
#
# Reads TOML config files using Python3's tomllib (3.11+) or tomli.
# Source this file in scripts to load configuration from TOML files.
#
# Usage:
#   source "$(dirname "$0")/modules/config_parser.sh"
#   load_config "/path/to/config.toml"
#   value=$(toml_get "section.key" "default_value")
#   array=($(toml_get_array "section.list"))
#######################################################################

# Resolve a single TOML value via Python
# Usage: toml_get <dotted.key> [default]
toml_get() {
    local key="$1"
    local default="${2:-}"
    local config_file="${_TOML_CONFIG_FILE:-}"

    if [[ -z "$config_file" || ! -f "$config_file" ]]; then
        echo "$default"
        return
    fi

    local fallback_file="${_TOML_FALLBACK_FILE:-}"

    local _raw
    _raw=$(python3 - "$key" "$config_file" "$fallback_file" "$default" <<'PYEOF'
import sys, os
try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        print(sys.argv[4])
        sys.exit(0)

def _resolve(config_path, keys):
    with open(config_path, 'rb') as f:
        data = tomllib.load(f)
    val = data
    for k in keys:
        if isinstance(val, dict) and k in val:
            val = val[k]
        else:
            return None
    return val

key, cfg, fb, default = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
keys = key.split('.')
val = _resolve(cfg, keys)
if val is None and fb:
    val = _resolve(fb, keys)
if val is None:
    print(default)
    sys.exit(0)

if isinstance(val, bool):
    print('true' if val else 'false')
elif isinstance(val, list):
    for item in val:
        print(item)
else:
    print(val)
PYEOF
) || _raw="$default"
    # Expand leading ~ to $HOME for portability
    echo "${_raw/#\~/$HOME}"
}

# Resolve a TOML array into bash array (one element per line)
# Usage: mapfile -t my_array < <(toml_get_array "section.list")
toml_get_array() {
    local key="$1"
    local config_file="${_TOML_CONFIG_FILE:-}"

    if [[ -z "$config_file" || ! -f "$config_file" ]]; then
        return
    fi

    python3 - "$key" "$config_file" <<'PYEOF'
import sys
try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        sys.exit(0)

key, cfg = sys.argv[1], sys.argv[2]
with open(cfg, 'rb') as f:
    data = tomllib.load(f)

keys = key.split('.')
val = data
for k in keys:
    if isinstance(val, dict) and k in val:
        val = val[k]
    else:
        sys.exit(0)

if isinstance(val, list):
    for item in val:
        s = str(item).strip()
        if not s.startswith('#'):
            print(s)
PYEOF
}

# Load a TOML config file (sets the file for subsequent toml_get calls)
# Usage: load_config "/path/to/config.toml"
load_config() {
    local config_file="$1"
    if [[ ! -f "$config_file" ]]; then
        echo "[config_parser] WARNING: Config file not found: $config_file" >&2
        return 1
    fi
    export _TOML_CONFIG_FILE="$config_file"
}

# Load a secondary TOML config that overrides the primary one
# Values from the override file take precedence.
# Usage: load_config_override "/path/to/override.toml"
load_config_with_fallback() {
    local primary="$1"
    local override="$2"

    if [[ -f "$override" ]]; then
        export _TOML_CONFIG_FILE="$override"
        export _TOML_FALLBACK_FILE="$primary"
    elif [[ -f "$primary" ]]; then
        export _TOML_CONFIG_FILE="$primary"
    else
        echo "[config_parser] WARNING: No config files found" >&2
        return 1
    fi
}

# Resolve the config directory relative to a script
# Usage: CONFIG_DIR=$(resolve_config_dir)
resolve_config_dir() {
    local script_dir="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)}"
    echo "${script_dir}/config/PPI"
}

# Batch-resolve multiple TOML keys in a single Python call (avoids N subprocess forks)
# Usage: eval "$(toml_get_bulk VAR1=section.key1 VAR2=section.key2:default ...)"
# Outputs: VAR1='value1'\nVAR2='value2'\n... — suitable for eval.
# Each item is VAR=dotted.key[:default]
toml_get_bulk() {
    local config_file="${_TOML_CONFIG_FILE:-}"
    local fallback_file="${_TOML_FALLBACK_FILE:-}"

    if [[ -z "$config_file" || ! -f "$config_file" ]]; then
        # Emit defaults only
        local spec
        for spec in "$@"; do
            local varname="${spec%%=*}"
            local rest="${spec#*=}"
            local default="${rest#*:}"
            [[ "$rest" == *:* ]] || default=""
            echo "${varname}='${default}'"
        done
        return
    fi

    # Pass specs as positional args after config/fallback
    python3 - "$config_file" "$fallback_file" "$@" <<'PYEOF'
import sys
try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        sys.exit(1)

def _resolve(path, keys):
    with open(path, 'rb') as f:
        data = tomllib.load(f)
    val = data
    for k in keys:
        if isinstance(val, dict) and k in val:
            val = val[k]
        else:
            return None
    return val

cfg, fb = sys.argv[1], sys.argv[2]

for spec in sys.argv[3:]:
    varname, rest = spec.split('=', 1)
    if ':' in rest:
        key, default = rest.split(':', 1)
    else:
        key, default = rest, ''
    keys = key.split('.')
    val = _resolve(cfg, keys)
    if val is None and fb:
        val = _resolve(fb, keys)
    if val is None:
        out = default
    elif isinstance(val, bool):
        out = 'true' if val else 'false'
    elif isinstance(val, list):
        out = ' '.join(str(x) for x in val)
    else:
        out = str(val)
    # Shell-safe single-quote escaping
    out = out.replace("'", "'\\''")
    print(f"{varname}='{out}'")
PYEOF
}
