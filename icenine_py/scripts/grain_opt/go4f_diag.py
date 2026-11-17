"""GO4f: why do rand voxels stay wrong after the best G1b pipeline (L/366/3,G/15/0.15,N/30/0.3)?

For each level, the voxels still more than 1 deg from the ground truth: gated by the seed path
(strict cost after the ladder above 0.15), improved by it (strict cost lowered), and whether the
strict cost at the ground truth is lower than at the end (search failure) or not (objective
failure: the wrong end looks at least as good to the strict cost).

Usage (from icenine_py/): uv run python scripts/grain_opt/go4f_diag.py --run go4e_rand_G15_N30 --sample rand
"""

import argparse
import gzip
import json

import common as CM  # noqa: I001
import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--thr", type=float, default=0.15)
    a = ap.parse_args()
    rows = json.load(gzip.open(CM.OUT / "runs" / a.run / "voxels.json.gz"))
    S = CM.Sample(CM.OUT / "go4f" / f"run_{a.sample}", start="centre", intensity=False, sample=a.sample)
    out = {}
    for lv in sorted({r["level"] for r in rows}):
        sel = [r for r in rows if r["level"] == lv]
        last = max(int(k[5:]) for k in sel[0] if k.startswith("err_s") and k[5:].isdigit())
        wrong = [r for r in sel if r[f"err_s{last}"] > 1.0]
        near = [r for r in sel if 0.1 < r[f"err_s{last}"] <= 1.0]
        res = dict(n=len(sel), wrong=len(wrong), near=len(near))
        if wrong:
            v = [r["voxel"] for r in wrong]
            st_gt = S.strict(S.T.R[v], v)
            st_end = np.array([r[f"strict_s{last}"] for r in wrong])
            gated = np.array([r["strict_s0"] > a.thr for r in wrong])
            improved = np.array([r["strict_s1"] < r["strict_s0"] for r in wrong])
            res.update(
                gated=int(gated.sum()),
                improved_by_seed=int(improved.sum()),
                gt_lower=int((st_gt < st_end).sum()),
                strict_gt_med=float(np.median(st_gt)),
                strict_end_med=float(np.median(st_end)),
                strict_after_ladder_med=float(np.median([r["strict_s0"] for r in wrong])),
                ungated_gt_lower=int(((st_gt < st_end) & ~gated).sum()),
                gated_gt_lower=int(((st_gt < st_end) & gated).sum()),
            )
        out[str(lv)] = res
        print(lv, json.dumps(res), flush=True)
    (CM.OUT / "go4f").mkdir(parents=True, exist_ok=True)
    (CM.OUT / "go4f" / f"{a.run}.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
