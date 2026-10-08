"""Rank analysis of the level-0 candidates saved by seed_diag.py --save-cands: where do the
candidates near the truth (symmetry-reduced misorientation < T deg) sit in the discrete-stage
ordering (what a top-K cap before the quick MC would use) and in the post-quick-MC ordering (what
the reconstructor keeps the top quarter of)?

Usage (from icenine_py/): uv run python scripts/phase_d/seed_diag_ranks.py --tag cap2k
"""

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
from csl import reduced_misorientation_deg  # noqa: E402

OUT = ICE / "benchmarks" / "phase_d_seed_diag"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="cap2k")
    ap.add_argument("--near-deg", type=float, nargs="+", default=[2.0, 3.0])
    a = ap.parse_args()
    rows = []
    for f in sorted(glob.glob(str(OUT / "runs" / f"{a.tag}_full_v*_cands.npz"))):
        z = np.load(f)
        v = int(Path(f).name.split("_v")[1].split("_")[0])
        Rt = z["R_true"]
        R_disc = z["L0_discrete_R"]  # (N,3,3) discrete order (sorted by discrete cost)
        score = z["L0_discrete_score"]  # 1 - confidence, ascending order
        perm = z["L0_quick_mc_perm"]  # sorted-after-MC -> discrete index
        cost_mc = z["L0_quick_mc_cost"]
        n = len(R_disc)
        err = reduced_misorientation_deg(R_disc, Rt)
        for T in a.near_deg:
            near = err < T
            if not near.any():
                rows.append(dict(voxel=v, near_deg=T, n=n, n_near=0))
                continue
            rank_disc = np.argsort(np.argsort(score, kind="stable"), kind="stable")  # 0 = best
            best_disc = int(rank_disc[near].min())
            # post-MC rank of each discrete index
            rank_mc = np.empty(n, int)
            rank_mc[perm] = np.arange(n)
            best_mc = int(rank_mc[near].min())
            keep = max(1, int(n * 0.25))
            rows.append(
                dict(
                    voxel=v,
                    near_deg=T,
                    n=n,
                    n_near=int(near.sum()),
                    best_near_rank_discrete=best_disc,
                    best_near_rank_after_quickmc=best_mc,
                    n_near_in_kept_quarter=int((rank_mc[near] < keep).sum()),
                    min_err_deg=float(err.min()),
                )
            )
    for r in rows:
        print(r)
    (OUT / f"rank_analysis_{a.tag}.json").write_text(json.dumps(rows, indent=1) + "\n")


if __name__ == "__main__":
    main()
