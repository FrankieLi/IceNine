"""Identity check of two BFS pilot runs (item 14, rule R2).

  uv run python scripts/autoresearch/item14/same_run.py RUN_A RUN_B

Equal means: ``recon.mic`` byte for byte; ``records.json`` and ``run.json`` equal once every clock
field is dropped (keys starting with ``wall``, ``setup_s``, ``seconds``, ``full_search_seconds``,
``timing_label``). Exit 0 if equal, 1 with the first difference named.
"""

import json
import sys
from pathlib import Path
from typing import Any, List, Optional

CLOCK_KEYS = {"setup_s", "seconds", "full_search_seconds", "timing_label"}


def strip_clock(obj: Any) -> Any:
    """``obj`` without the clock fields, at any depth."""
    if isinstance(obj, dict):
        return {
            k: strip_clock(v)
            for k, v in obj.items()
            if not (str(k).startswith("wall") or k in CLOCK_KEYS)
        }
    if isinstance(obj, list):
        return [strip_clock(v) for v in obj]
    return obj


def first_difference(a: Any, b: Any, path: str = "") -> Optional[str]:
    """Path and values of the first difference of two JSON trees (None if equal)."""
    if type(a) is not type(b):
        return f"{path or '.'}: {a!r} vs {b!r}"
    if isinstance(a, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                return f"{path}/{k}: present in only one run"
            d = first_difference(a[k], b[k], f"{path}/{k}")
            if d:
                return d
        return None
    if isinstance(a, list):
        if len(a) != len(b):
            return f"{path}: length {len(a)} vs {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            d = first_difference(x, y, f"{path}[{i}]")
            if d:
                return d
        return None
    return None if a == b else f"{path}: {a!r} vs {b!r}"


def diff_runs(a: Path, b: Path) -> Optional[str]:
    """First difference between two run folders, or None if they are identical."""
    ma, mb = (a / "recon.mic").read_bytes(), (b / "recon.mic").read_bytes()
    if ma != mb:
        return f"recon.mic differs ({len(ma)} vs {len(mb)} bytes)"
    for name in ("records.json", "run.json"):
        ja = strip_clock(json.loads((a / name).read_text()))
        jb = strip_clock(json.loads((b / name).read_text()))
        d = first_difference(ja, jb)
        if d:
            return f"{name}: {d}"
    return None


def main(argv: Optional[List[str]] = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print("usage: same_run.py RUN_A RUN_B", file=sys.stderr)
        return 2
    d = diff_runs(Path(args[0]), Path(args[1]))
    if d:
        print(f"RUNS DIFFER (runs differ): {d}")
        return 1
    print(f"runs identical: {args[0]} == {args[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
