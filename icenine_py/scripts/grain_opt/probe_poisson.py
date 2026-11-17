"""GO9 (I18): bias of a Poisson likelihood on pixel counts as an orientation objective.

GO1's protocol (``probe_bias.py``): level-0 ``grad`` grains, CMA over the grain rotation only,
started at the ground truth; an unbiased objective stays there. Objective: per spot window,
mu = s_j m + b over the union of observed lit pixels and pixels with predicted coverage mass m
(grain plus the fixed neighbour shell), y = clean intensity / gain, s_j = sum(y) / sum(m),
Poisson deviance summed over spots.

Usage (from icenine_py/): uv run python scripts/grain_opt/probe_poisson.py --grains 381 302
"""

import argparse
import json
import time

import common as CM  # noqa: I001
import numpy as np
import torch

from icenine import grain_fit as GF

W = 1 << 12


def enc(j: np.ndarray, p: np.ndarray) -> np.ndarray:
    return ((j.astype(np.int64) * W + p[..., 2]) * W + p[..., 1]) * W + p[..., 0]


class Poisson:
    def __init__(self, prob: object, gain: float, b: float) -> None:
        self.prob, self.b = prob, b
        ks, ys = [], []
        for j, sp in enumerate(prob.region.spots):
            p = np.floor(np.asarray(sp.obs)).astype(np.int64)
            ks.append(enc(np.full(len(p), j), p))
            ys.append(np.asarray(sp.obs_val, float) / gain)
        self.ok = np.concatenate(ks)
        self.oy = np.concatenate(ys)
        self.S = len(prob.region.spots)

    def __call__(self, R: torch.Tensor) -> float:
        prob = self.prob
        with torch.no_grad():
            rows = prob.model.predict(R, prob.V)
            x, a = GF.pack_predicted(prob.region, rows, prob.bary, prob.model.n_det, with_shell=True)
        x, a = x.numpy(), a.numpy()
        p = np.clip(np.floor(x).astype(np.int64), 0, W - 1)
        k = enc(np.broadcast_to(np.arange(self.S)[:, None], a.shape), p)[a > 0]
        m = a[a > 0]
        keys = np.concatenate([k, self.ok])
        uk, inv = np.unique(keys, return_inverse=True)
        mm = np.bincount(inv[: len(k)], weights=m, minlength=len(uk))
        yy = np.bincount(inv[len(k):], weights=self.oy, minlength=len(uk))
        spot = uk // (W * W * W)
        sm = np.bincount(spot, weights=mm, minlength=self.S)
        sy = np.bincount(spot, weights=yy, minlength=self.S)
        s = np.where(sm > 0, sy / np.maximum(sm, 1e-12), 0.0)
        mu = s[spot] * mm + self.b
        dev = mu - yy + np.where(yy > 0, yy * np.log(np.maximum(yy, 1e-12) / mu), 0.0)
        return float(2.0 * dev.sum())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grains", type=int, nargs="+", default=[381, 302])
    ap.add_argument("--evals", type=int, default=150)
    a = ap.parse_args()
    import cma

    S = CM.Sample(CM.OUT / "go9" / "run")
    T = S.T
    gm = np.asarray(T.gmap)
    allv = np.concatenate([S.images.store["val"]]) if getattr(S.images, "store", None) else None
    gain = float(np.percentile(allv, 5)) / 50.0 if allv is not None else 1244.0
    out = []
    for g in a.grains:
        members = [int(v) for v in np.nonzero(gm == g)[0]]
        shell = GF.shell_of(members, S.nb)
        Rt = T.R[members]
        prob = GF.GrainProblem(S.model, S.images, S.R0, S.V, S.xy, members, shell, S.sym,
                               spec=GF.ObjectiveSpec(kind="ls"), strict_fn=S.strict)
        Rt_al = GF.align_to(Rt, prob.field.R_ref, S.sym)
        w_t = GF.matrix_to_rotvec(Rt_al[:1] @ prob.field.R_ref.T)[0]
        th0 = np.zeros(prob.field.n_params)
        th0[:3] = w_t
        for name, b in (("pois", 0.5), ("pois_b5", 5.0)):
            obj = Poisson(prob, gain, b)

            def f(w: np.ndarray) -> float:
                th = th0.copy()
                th[:3] = w
                with torch.no_grad():
                    return obj(prob.field.orientations(torch.as_tensor(th)))

            t0 = time.time()
            es = cma.CMAEvolutionStrategy(list(w_t), np.radians(0.03),
                                          {"seed": 1, "maxfevals": a.evals, "verbose": -9})
            es.optimize(f)
            w = np.asarray(es.result.xbest)
            th = th0.copy()
            th[:3] = w
            with torch.no_grad():
                R_opt = prob.field.orientations(torch.as_tensor(th)).numpy()
            row = dict(grain=g, variant=name, gain=gain, n_spots=prob.n_spots,
                       offset_deg=float(np.median(T.err(R_opt, members))), f_truth=f(w_t),
                       f_opt=f(w), strict_truth=float(S.strict(Rt, members).mean()),
                       strict_opt=float(S.strict(R_opt, members).mean()), t_s=time.time() - t0)
            print(json.dumps(row), flush=True)
            out.append(row)
    (CM.OUT / "go9" / "poisson.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
