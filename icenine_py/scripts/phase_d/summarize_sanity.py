"""Print the sanity JSONs (results/sanity_*.json) as a compact table."""

import json
import sys
from pathlib import Path

import numpy as np

res = Path(__file__).parent / "results"
for f in sorted(res.glob("sanity_*.json")):
    d = json.loads(f.read_text())
    print(f.name)
    if "containment" in d:
        r = d["containment"]["rows"]
        print("  serial-path containment:", [(x["lit"], x["frac"]) for x in r])
    if "cost" in d:
        rows = d["cost"]["rows"]
        keys = ["q_true", "q_0.25deg_mean", "q_0.5deg_mean", "q_1.0deg_mean", "q_1.0deg_max"]
        for k in keys:
            v = np.array([x[k] for x in rows])
            print(f"  {k:16s} min {v.min():.3f}  mean {v.mean():.3f}  max {v.max():.3f}")
        gap = min(x["q_true"] - x["q_1.0deg_max"] for x in rows)
        print(f"  min (q_true - q_1deg_max) over {len(rows)} voxels: {gap:.3f}")
