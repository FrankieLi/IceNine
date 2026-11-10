"""Print the headline of a paint run folder: uv run python scripts/autoresearch/ar3/peek.py NAME"""

import json
import sys
from pathlib import Path

d = (
    Path(__file__).resolve().parents[3]
    / "benchmarks"
    / "autoresearch"
    / "ar3"
    / "runs"
    / sys.argv[1]
)
r = json.loads((d / "run.json").read_text())
for k in (
    "n_entries",
    "evals",
    "evals_total",
    "counts",
    "wall_total_s",
    "wall_by_stage_s",
    "index_pairing_s",
    "setup_s",
    "side",
):
    print(k, r.get(k))
p = json.loads((d / "paint.json").read_text())
for e in p["entries"][:8]:
    print({k: v for k, v in e.items() if k != "R"})
print(p["log"][:3])
