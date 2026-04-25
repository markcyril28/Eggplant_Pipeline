#!/usr/bin/env python3
"""
Deep-merge multiple TOML files and output valid TOML to stdout.
Later files override earlier ones at the leaf level (scalars and lists
replace; dicts are recursively merged).

Missing files are silently skipped so orchestrators can list shared +
group files without checking existence first.

Usage:
    python3 merge_toml.py shared/00.toml shared/01.toml group/00.toml group/01.toml > merged.toml
"""

import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


def deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge *override* into *base*.

    - Dicts: recurse into matching keys.
    - Everything else (scalars, lists): override replaces base.
    """
    result = dict(base)
    for key, val in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = deep_merge(result[key], val)
        else:
            result[key] = val
    return result


# ---------------------------------------------------------------------------
# Minimal TOML serializer (supports the types used in this pipeline)
# ---------------------------------------------------------------------------

def _fmt_value(val) -> str:
    """Format a single TOML value."""
    if isinstance(val, bool):
        return "true" if val else "false"
    if isinstance(val, int):
        return str(val)
    if isinstance(val, float):
        return str(val)
    if isinstance(val, str):
        # Escape backslashes and double-quotes inside the string
        escaped = val.replace("\\", "\\\\").replace('"', '\\"')
        return f'"{escaped}"'
    if isinstance(val, list):
        if len(val) == 0:
            return "[]"
        items_str = ", ".join(_fmt_value(v) for v in val)
        # Multi-line for long or many-element arrays
        if len(items_str) > 80 or len(val) > 3:
            nl = "\n".join(f"    {_fmt_value(v)}," for v in val)
            return f"[\n{nl}\n]"
        return f"[{items_str}]"
    # Fallback
    return f'"{val}"'


def _serialize(data: dict, section_path: str = "") -> list[str]:
    """Recursively serialize a dict to TOML lines."""
    lines: list[str] = []
    # Pass 1: emit scalar / list keys at this level
    for key, val in data.items():
        if not isinstance(val, dict):
            lines.append(f"{key} = {_fmt_value(val)}")
    # Pass 2: emit sub-tables
    for key, val in data.items():
        if isinstance(val, dict):
            path = f"{section_path}.{key}" if section_path else key
            lines.append("")
            lines.append(f"[{path}]")
            lines.extend(_serialize(val, path))
    return lines


def main() -> None:
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} file1.toml [file2.toml ...]", file=sys.stderr)
        sys.exit(1)

    merged: dict = {}
    for path_str in sys.argv[1:]:
        path = Path(path_str)
        if not path.exists():
            # Silently skip — allows listing shared + group files even when
            # a gene group hasn't defined a stage-specific override file yet.
            continue
        with open(path, "rb") as f:
            data = tomllib.load(f)
        merged = deep_merge(merged, data)

    if not merged:
        print("Error: no valid TOML files found", file=sys.stderr)
        sys.exit(1)

    for line in _serialize(merged):
        print(line)


if __name__ == "__main__":
    main()
