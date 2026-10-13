#!/usr/bin/env python3
"""Review-fix extras: paired level-0 rank counts (30 voxels) and E0 medians over the voxels
right under all arms.

  uv run python scripts/phase_d/ray_overlap_review_extras.py

Writes benchmarks/phase_d_seed_diag/ray_overlap/review_extras.json.
Rank tie rule: equal ranks are "tied".
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "ray_overlap"
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
import common as C  # noqa: E402

CACHE = HERE / "cache"
DIRS = {"old": "off_ms0", "cpp": "on_ms0", "R2": "rayR2_ms0", "A": "rayA_ms0"}


def main() -> None:
    res: Dict[str, Any] = {}
    for var in ("clean", "all"):
        label = "realistic" if var == "all" else var
        rk: Dict[int, Dict[str, int]] = {}
        data: Dict[int, Dict[str, Any]] = {}
        for p in (OUT / "coarse30").glob(f"{label}_v*.json"):
            r = json.loads(p.read_text())
            data.setdefault(r["voxel"], {}).update(r["arms"])
        vox = sorted(
            v for v, a in data.items() if all(k in a for k in ("old", "cpp", "f0_R2", "f0_A"))
        )
        for arm in ("f0_R2", "f0_A", "old"):
            c = dict(better=0, worse=0, tied=0)
            for v in vox:
                a = data[v][arm]["levels"][0]["rank_closest_disc"]
                b = data[v]["cpp"]["levels"][0]["rank_closest_disc"]
                c["better" if a < b else "worse" if a > b else "tied"] += 1
            res[f"{label}|rank_L0|{arm}_vs_cpp"] = dict(n=len(vox), **c)
        names = sorted(p.name for p in (CACHE / DIRS["old"]).glob(f"*_{var}.npz"))
        ok = [n for n in names if all((CACHE / d / n).exists() for d in DIRS.values())]
        errs = {}
        for a, d in DIRS.items():
            z = [np.load(CACHE / d / n) for n in ok]
            errs[a] = C.err_deg(np.stack([x["R"] for x in z]), np.stack([x["R_true"] for x in z]))
        common = np.all([e <= C.WRONG_DEG for e in errs.values()], axis=0)
        res[f"{label}|common_right"] = dict(
            n=int(common.sum()),
            **{f"median_{a}": float(np.median(e[common])) for a, e in errs.items()},
        )
        for a, e in errs.items():
            res[f"{label}|own_right_n|{a}"] = int((e <= C.WRONG_DEG).sum())
    (OUT / "review_extras.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
