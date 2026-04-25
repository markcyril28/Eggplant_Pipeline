"""Stamped log helper shared by v2 CRISPR pipeline modules.

Matches the format used by modules/logging/logging_utils.sh so stdout from
stage scripts interleaves cleanly with orchestrator logs:

    [2026-04-19 09:59:17] [INFO] [02_rescore] 4 guides rescored

Use `_log(msg)` for INFO, `_log(msg, level="WARN")` / `"ERROR"` for others.
Routes to stderr so it merges into the orchestrator's error/full log streams
without interfering with any data the script might print to stdout.
"""
from __future__ import annotations

import sys
from datetime import datetime


def _log(msg: str, level: str = "INFO", file=sys.stderr) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] [{level}] {msg}", file=file, flush=True)
