"""Print (and optionally save) the machine state before a timing run.

Usage: timing_preflight.py [--require] [--json OUT.json] [--max-load 1.5] [--allow-battery]
With --require, exits 1 if the machine is not quiet (see preflight.require_quiet).
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Optional, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "icenine_py" / "scripts" / "common"))
import preflight  # noqa: E402


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--require", action="store_true")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--max-load", type=float, default=1.5)
    ap.add_argument("--allow-battery", action="store_true")
    args = ap.parse_args(argv)

    info = preflight.preflight()
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=2))
    if args.require:
        try:
            preflight.require_quiet(args.max_load, args.allow_battery, info=info)
        except preflight.MachineBusyError as e:
            print(f"error: {e}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
