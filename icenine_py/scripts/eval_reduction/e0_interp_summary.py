"""Task A step 3: E0 clean (200 voxels, per-voxel images, seed 0, `mc` arm, min_sin_eta 0), the
opt-in screen switches against the stored arms, paired by voxel.

Arms (scripts/phase_d/e0_cpprule.py, run with --msin 0 --tag a --variant clean):
  edge    C++ overlap rule + CppEdgeRules                      vs  on  (C++ rule, centroid)
  interp  C++ rule + CppEdgeRules + CppOmegaInterp + BatchedScreen   vs  on, and vs  edge
  lenint  lenient rule + CppEdgeRules + CppOmegaInterp + BatchedScreen   vs  off (lenient, centroid)
Stored arms: scripts/phase_d/cache/{on,off}_ms0_centroid (Task 0 re-measure). Wrong = more than
1 degree from the truth (symmetry-reduced). Paired counts are by voxel: fixed = wrong before and
right after, broken = right before and wrong after.

  uv run python scripts/eval_reduction/e0_interp_summary.py --out benchmarks/eval_reduction/A/e0
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
CACHE = ICE / "scripts" / "phase_d" / "cache"
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import common as C  # noqa: E402
import doc_tables as DT  # noqa: E402
import stats as S  # noqa: E402

PAIRS = [
    ("edge", "on", "C++ rule + edge rules vs C++ rule"),
    ("interp", "on", "C++ rule + edge rules + omega interpolation + batched vs C++ rule"),
    ("interp", "edge", "omega interpolation + batched vs edge rules alone (C++ rule)"),
    ("lenint", "off", "lenient + edge rules + interpolation + batched vs lenient (default)"),
]
DIRS = {
    "on": "on_ms0_centroid",
    "off": "off_ms0_centroid",
    "edge": "edge_ms0_a",
    "interp": "interp_ms0_a",
    "lenint": "lenint_ms0_a",
}


def load(arm: str) -> Dict[str, Dict[str, Any]]:
    d = CACHE / DIRS[arm]
    return {p.name: dict(np.load(p)) for p in sorted(d.glob("v*_clean.npz"))} if d.exists() else {}


def arm_stats(z: Dict[str, Dict[str, Any]], names: List[str]) -> Dict[str, Any]:
    Rt = np.stack([z[n]["R_true"] for n in names])
    e = C.err_deg(np.stack([z[n]["R"] for n in names]), Rt)
    ev = np.array([int(z[n]["evals_global"]) + int(z[n]["evals_local"]) for n in names])
    return dict(
        err=e,
        wrong=e > C.WRONG_DEG,
        evals=ev,
        n_l0=np.array([int(z[n]["n_l0"]) for n in names]),
        runtime=np.array([float(z[n]["runtime"]) for n in names]),
        basin_miss=np.array([float(z[n]["min_dist_l0"]) > 3.0 for n in names]),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    data = {k: load(k) for k in DIRS}
    summ: Dict[str, Any] = {}
    rows = []
    for new, base, label in PAIRS:
        if not data[new] or not data[base]:
            continue
        names = sorted(set(data[new]) & set(data[base]))
        A, B = arm_stats(data[new], names), arm_stats(data[base], names)
        n = len(names)
        broken, fixed = S.paired_discordant(A["wrong"], B["wrong"])
        lo_a, hi_a = S.wilson(int(A["wrong"].sum()), n)
        lo_b, hi_b = S.wilson(int(B["wrong"].sum()), n)
        both_ok = ~A["wrong"] & ~B["wrong"]
        s = dict(
            label=label,
            n=n,
            wrong_base=int(B["wrong"].sum()),
            wrong_new=int(A["wrong"].sum()),
            broken=broken,
            fixed=fixed,
            both_wrong=int((A["wrong"] & B["wrong"]).sum()),
            evals_base_mean=float(B["evals"].mean()),
            evals_new_mean=float(A["evals"].mean()),
            l0_base_median=float(np.median(B["n_l0"])),
            l0_new_median=float(np.median(A["n_l0"])),
            basin_miss_base=int(B["basin_miss"].sum()),
            basin_miss_new=int(A["basin_miss"].sum()),
            med_err_both_right_base=float(np.median(B["err"][both_ok])),
            med_err_both_right_new=float(np.median(A["err"][both_ok])),
            runtime_base_mean_s=float(B["runtime"].mean()),
            runtime_new_mean_s=float(A["runtime"].mean()),
            wall_note="multi-worker, contended: not a timing claim",
        )
        summ[f"{new}_vs_{base}"] = s
        rows.append(
            dict(
                new=new,
                base=base,
                n=n,
                wrong_base=f"{s['wrong_base']}/{n} [{lo_b:.3f}, {hi_b:.3f}]",
                wrong_new=f"{s['wrong_new']}/{n} [{lo_a:.3f}, {hi_a:.3f}]",
                fixed=fixed,
                broken=broken,
                both_wrong=s["both_wrong"],
                evals_base=s["evals_base_mean"],
                evals_new=s["evals_new_mean"],
                l0_base=s["l0_base_median"],
                l0_new=s["l0_new_median"],
                basin_miss_base=s["basin_miss_base"],
                basin_miss_new=s["basin_miss_new"],
                runtime_base=s["runtime_base_mean_s"],
                runtime_new=s["runtime_new_mean_s"],
            )
        )
    (out / "summary.json").write_text(json.dumps(summ, indent=1, default=float))
    tbl = DT.markdown_table(
        rows,
        [
            "new",
            "base",
            "n",
            "wrong_base",
            "wrong_new",
            "fixed",
            "broken",
            "both_wrong",
            "evals_base",
            "evals_new",
            "l0_base",
            "l0_new",
            "basin_miss_base",
            "basin_miss_new",
            "runtime_base",
            "runtime_new",
        ],
        formats={
            "evals_base": ",.0f",
            "evals_new": ",.0f",
            "l0_base": ".0f",
            "l0_new": ".0f",
            "runtime_base": ".1f",
            "runtime_new": ".1f",
        },
    )
    DT.write_tables(out / "tables.md", {"e0_interp_paired": tbl})
    print(tbl)


if __name__ == "__main__":
    main()
