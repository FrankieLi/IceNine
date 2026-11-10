"""Print the recall table and the misses of a recall_<kind>_<tag>.json (a quick look)."""

import json
import sys
from pathlib import Path

import numpy as np

d = json.loads(
    (Path(__file__).resolve().parents[3] / "benchmarks/autoresearch/ar1" / sys.argv[1]).read_text()
)
rows = d["rows"]
print("n", len(rows))
for t in (0.5, 2, 5):
    print(
        "within",
        t,
        {K: sum(1 for r in rows if r["dist"] and min(r["dist"][:K]) < t) for K in (1, 5, 20)},
    )
print("median wall s", np.median([r["wall_s"] for r in rows]))
nmiss = int(sys.argv[2]) if len(sys.argv) > 2 else 40
for r in rows:
    if not (r["dist"] and min(r["dist"][:20]) < 2):
        if nmiss <= 0:
            break
        nmiss -= 1
        print(
            r["voxel"],
            r["n_att"],
            r["n_true_explained"],
            r["votes_at_truth"],
            [round(x, 1) for x in r["dist"][:6]],
            r["explained"][:6],
            len(r["dist"]),
        )
