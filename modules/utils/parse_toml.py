#!/usr/bin/env python3
"""
TOML Configuration Parser for Eggplant Pipeline
Parses a TOML config file and outputs shell-compatible variable exports.

Usage:
    python3 parse_toml.py <config.toml>                  # Export all variables
    python3 parse_toml.py <config.toml> <section>        # Export one section
    python3 parse_toml.py <config.toml> <section> <key>  # Get single value

Output is eval-safe shell assignments:
    SECTION__KEY="value"
    SECTION__LIST=("val1" "val2")
"""

import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

# Force LF-only output even on Windows — otherwise Python's text-mode stdout
# appends \r to every line, and bash `mapfile -t` leaves those \r in array
# elements (breaking `[[ -d path ]]` and every downstream use).
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(newline="\n")
    except AttributeError:
        pass


def flatten_toml(data: dict, prefix: str = "") -> list[str]:
    """Convert nested TOML dict into flat shell variable assignments."""
    lines = []
    for key, value in data.items():
        var_name = f"{prefix}__{key}".upper() if prefix else key.upper()
        # Sanitize variable name: replace hyphens/dots with underscores
        var_name = var_name.replace("-", "_").replace(".", "_")

        if isinstance(value, dict):
            lines.extend(flatten_toml(value, var_name))
        elif isinstance(value, list):
            # Shell array: VARNAME=("v1" "v2" ...)
            items = " ".join(f'"{v}"' for v in value)
            lines.append(f'{var_name}=({items})')
        elif isinstance(value, bool):
            lines.append(f'{var_name}={"true" if value else "false"}')
        elif isinstance(value, (int, float)):
            lines.append(f'{var_name}="{value}"')
        else:
            lines.append(f'{var_name}="{value}"')
    return lines


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <config.toml> [key ...]", file=sys.stderr)
        sys.exit(1)

    config_path = Path(sys.argv[1])
    if not config_path.exists():
        print(f"Error: Config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)

    with open(config_path, "rb") as f:
        config = tomllib.load(f)

    keys = sys.argv[2:]

    if not keys:
        for line in flatten_toml(config):
            print(line)
        return

    # Navigate arbitrarily deep: get_toml section subsection ... key
    data = config
    for i, key in enumerate(keys):
        if not isinstance(data, dict) or key not in data:
            path = " → ".join(keys[: i + 1])
            print(f"Error: Key path '{path}' not found in config", file=sys.stderr)
            sys.exit(1)
        data = data[key]

    if isinstance(data, dict):
        for line in flatten_toml(data, keys[-1]):
            print(line)
    elif isinstance(data, list):
        for item in data:
            print(item)
    elif isinstance(data, bool):
        print("true" if data else "false")
    else:
        print(data)


if __name__ == "__main__":
    main()
