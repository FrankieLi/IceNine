"""GO7f-d: where does the rocking damage to the strict cost come from?

For a fixed set of G20 voxels (every voxel of grain 197's pick plus a seeded random 200 of the
other G20 voxels), evaluates the strict cost at the ground-truth orientation on three image sets:
clean ``grad``, ``grad_rock10``, and ``grad_rock10`` with the strict cost's rocking window. Per arm:
mean strict cost and the share of predicted peaks on the detector with no overlapping lit pixel
(``1 - peak_overlap / peak_on_detector``, summed over the voxels).

Usage (from icenine_py/): uv run python scripts/grain_opt/go7f_diag.py [--n 200] [--rocking 1.0]
"""

import argparse
import gzip
import json
import time

import common as CM  # noqa: I001
import numpy as np
import torch


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--rocking", type=float, default=1.0)
    ap.add_argument("--stacks", nargs="*", default=None, help="GO7g: stack tags (no window) in place of the GO7f arms")
    ap.add_argument("--out", default="diag_d.json")
    a = ap.parse_args()
    rows = json.load(gzip.open(CM.OUT / "runs" / "go7_clean_As100" / "voxels.json.gz"))
    g197 = [r["voxel"] for r in rows if r["pick"] == 197]
    rest = [r["voxel"] for r in rows if r["pick"] != 197]
    rng = np.random.default_rng(0)
    pick = sorted(rng.choice(rest, size=min(a.n, len(rest)), replace=False).tolist())
    arms = [("clean", None, 0.0), ("rock10", "grad_rock10_clean.npy", 0.0),
            ("rock10_window", "grad_rock10_clean.npy", a.rocking)]
    if a.stacks:
        arms = [(s, None if s == "grad" else f"{s}_clean.npy", 0.0) for s in a.stacks]
    out = {"n_g197": len(g197), "n_rest": len(pick), "rocking_deg": a.rocking, "arms": {}}
    for name, stack, rock in arms:
        t0 = time.time()
        S = CM.Sample(CM.OUT / "runs" / f"go7f_diag_{name}", sample="grad", stack=stack, rocking_deg=rock)
        res = {}
        for sel_name, sel in (("g197", g197), ("rest", pick)):
            cost, po, pd = [], 0, 0
            for v in sel:
                info = S.cost_fn.evaluate(
                    np.asarray(S.T.R[v], np.float64), torch.as_tensor(S.V[v], dtype=torch.float32), S.phase
                )
                cost.append(float(info.cost))
                po += int(info.peak_overlap)
                pd += int(info.peak_on_detector)
            res[sel_name] = dict(mean_strict=float(np.mean(cost)), median_strict=float(np.median(cost)),
                                 unlit_peak_share=float(1 - po / max(pd, 1)), peaks=pd)
        res["t_s"] = time.time() - t0
        out["arms"][name] = res
        print(json.dumps({name: res}), flush=True)
    (CM.OUT / "go7f").mkdir(exist_ok=True)
    json.dump(out, open(CM.OUT / "go7f" / a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
