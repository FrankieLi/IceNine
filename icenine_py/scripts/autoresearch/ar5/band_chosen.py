"""AR5 step 15: print the form the gate chose (``ar5/band/chosen.json``, committed and clean).

  uv run python scripts/autoresearch/ar5/band_chosen.py        # prints PB1000 or PB500, exit 0
                                                              # exit 2 (and a reason) otherwise

No import of the study's other modules (the run-all script calls it before anything starts).
"""

import json
import subprocess
import sys
from pathlib import Path

ICE = Path(__file__).resolve().parents[3]
CHOSEN = ICE / "benchmarks" / "autoresearch" / "ar5" / "band" / "chosen.json"


def main() -> int:
    if not CHOSEN.exists():
        print("chosen.json does not exist (the gate has not chosen a form)", file=sys.stderr)
        return 2
    log = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", str(CHOSEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    st = subprocess.run(
        ["git", "status", "--porcelain", "--", str(CHOSEN)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not log or st:
        print("chosen.json is not committed or differs from its commit", file=sys.stderr)
        return 2
    print(json.loads(CHOSEN.read_text())["chosen"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
