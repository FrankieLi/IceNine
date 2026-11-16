"""Probe the objectives on a few grains: value at the truth, at the painted start and at the truth
turned by a fixed angle about random axes (whole grain), with timings.

Usage (from icenine_py/): uv run python scripts/grain_opt/probe_objective.py --grains 64 381
"""

import argparse
import json
import time

import common as CM  # noqa: I001
import numpy as np
import torch

from icenine import grain_fit as GF


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grains", type=int, nargs="+", default=[64, 381])
    ap.add_argument("--n-sub", type=int, default=2)
    ap.add_argument("--blur", type=float, default=1.0)
    ap.add_argument("--diameter", type=float, default=12.0)
    ap.add_argument("--reach", type=float, default=0.0)
    ap.add_argument("--lam", type=float, default=10.0)
    ap.add_argument("--keep", type=int, default=3)
    a = ap.parse_args()
    S = CM.Sample(CM.OUT / "probe")
    T = S.T
    gm = np.asarray(T.gmap)
    bary = torch.as_tensor(GF.subtriangle_bary(a.n_sub))
    rng = np.random.default_rng(0)
    reach = a.reach or None
    for g in a.grains:
        members = [int(v) for v in np.nonzero(gm == g)[0]]
        shell = GF.shell_of(members, S.nb)
        reg = GF.build_region(S.model, S.images, S.R0, S.V, members, shell, bary, keep_px=a.keep)
        y, b, _ = GF.pack_observed(reg)
        obs_g = GF.observed_grids(reg, intensity=True)

        def objective(Rm):
            with torch.no_grad():
                pr = S.model.predict(torch.as_tensor(Rm), torch.as_tensor(S.V[members]))
                x, am = GF.pack_predicted(reg, pr, bary, S.model.n_det)
                t = time.time()
                sd = GF.sinkhorn_divergence(
                    x, am, y, b, blur_px=a.blur, lam=a.lam, diameter_px=a.diameter, reach_px=reach
                )
                ts = time.time() - t
                t = time.time()
                ls = GF.ls_objective(GF.soft_raster(reg, x, am, 0.5, 0.1), obs_g)
                tl = time.time() - t
            return float(sd.mean()), float(ls.mean()), ts, tl

        res = {"truth": objective(T.R[members]), "paint": objective(S.R0[members])}
        for ang in (0.02, 0.05, 0.1, 0.2, 0.5):
            vals = []
            for _ in range(3):
                ax = rng.standard_normal(3)
                ax /= np.linalg.norm(ax)
                dR = GF.rotvec_to_matrix(torch.as_tensor(ax * np.radians(ang))).numpy()
                vals.append(objective(dR @ T.R[members])[:2])
            res[f"truth+{ang}"] = np.mean(vals, 0).tolist()
        print(
            json.dumps(
                dict(grain=g, level=T.level_of[g], n_vox=len(members), S=len(reg.spots),
                     P=int(a.n_sub ** 2 * len(members)), O=int(y.shape[1]),
                     obs_med=float(np.median([len(s.obs) for s in reg.spots])), res=res)
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
