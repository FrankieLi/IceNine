"""GO5b (idea I38): per-grain strain fit on sample ``strain`` with the coverage cost.

  uv run python scripts/grain_opt/fit_strain.py --arm E@GT|FE@GT|FE@P [--out go5b]
      [--evals 600] [--only <grain id>]

Benchmark S12: 4 grains per strain level (1e-4, 5e-4, 1e-3 tensor RMS), sizes at quantiles 0.2,
0.4, 0.6, 0.8 of the level; true grain ids as membership. Cost: 1 - share of predicted coverage
mass on lit pixels, mean over the grain's spots (GO5 "coverage"). Spot windows are built without
strain. Parameters: grain field (FieldA order 1, 9 numbers) and one crystal-frame strain per
grain (6). CMA with per-coordinate steps 0.05 deg (field) and 5e-4 (strain), strain from 0.

GO5e (I40): ``--cost sd|int`` replaces the coverage cost by a sub-pixel cost (``SubPix``: signed
distance to the lit set, or clean intensity) in any arm; ``--seed`` offsets every CMA seed.

GO5d arms (I41), strain only (``--evals_e``) at a fixed field: ``E@F`` (field fitted to the
ground truth); ``O@P`` (field leg ``--evals_f`` on the coverage cost at zero strain); ``Ost@P``
(field leg on the strict cost, no strain model).

GO5c arms (I39): ``J@P`` (field leg, then one joint leg of ``--evals``); ``S@P`` (field leg,
strain-only leg ``--evals_e`` at the fixed field, joint leg ``--evals`` at steps 0.02 deg and
2e-4); ``SR@P`` (``S@P`` plus two joint legs of 300, steps halved each leg); ``S@GT`` (``S@P``
from the field fitted to the ground truth). GO5b arms: ``E@GT`` (orientations at the ground truth, strain only, 300 evaluations); ``FE@GT``
(field + strain from the field fitted to the ground truth); ``FE@P`` (painted start: field-only
CMA, 200 evaluations, then field + strain; the start field ignores painted outliers,
``robust_theta0``). One JSON line per grain; all rows in
``benchmarks/grain_opt/<out>/<arm>.json``.
"""

import argparse
import json
import math
import time
from typing import Any, Dict, List

import common as CM  # noqa: I001
import numpy as np
import torch

from icenine import grain_fit as GF

QUANTILES = (0.2, 0.4, 0.6, 0.8)
LEVELS = (1e-4, 5e-4, 1e-3)
STEP_ROT = math.radians(0.05)
STEP_EPS = 5e-4


def tensor_rms(e: np.ndarray) -> float:
    e = np.asarray(e, float)
    return float(np.sqrt(((e[:3] ** 2).sum() + 2 * (e[3:] ** 2).sum()) / 3.0))


def pick_grains(T: Any) -> List[int]:
    gm = np.asarray(T.gmap)
    out = []
    for lv in LEVELS:
        gids = [int(g) for g, v in T.strain_level_of.items() if abs(float(v) - lv) < 1e-12]
        gids = sorted(gids, key=lambda g: (int((gm == g).sum()), g))
        out += [gids[int(q * (len(gids) - 1))] for q in QUANTILES]
    return out


class Coverage:
    """Vectorised GO5 coverage cost of one grain problem (grain rows only, no shell)."""

    W = 1 << 12  # column / row / frame ranges are below 4096

    def __init__(self, prob: Any) -> None:
        self.prob = prob
        keys = [
            self._enc(j, np.floor(np.asarray(sp.obs)).astype(np.int64))
            for j, sp in enumerate(prob.region.spots)
        ]
        self.lit = np.unique(np.concatenate(keys)) if keys else np.zeros(0, np.int64)
        self.S = len(prob.region.spots)

    def _enc(self, j: Any, p: np.ndarray) -> np.ndarray:
        W = self.W
        return ((np.asarray(j, np.int64) * W + p[..., 2]) * W + p[..., 1]) * W + p[..., 0]

    def __call__(self, R: torch.Tensor, strain: torch.Tensor) -> float:
        prob = self.prob
        prob.n_evals += 1
        with torch.no_grad():
            rows = prob.model.predict(R, prob.V, strain)
            x, a = GF.pack_predicted(
                prob.region, rows, prob.bary, prob.model.n_det, with_shell=False
            )
        x, a = x.numpy(), a.numpy()
        p = np.clip(np.floor(x).astype(np.int64), 0, self.W - 1)
        on = np.isin(self._enc(np.arange(self.S)[:, None], p), self.lit)
        tot = a.sum(1)
        ok = tot > 0
        if not ok.any():
            return 1.0
        return 1.0 - float(((a * on).sum(1)[ok] / tot[ok]).mean())


class SubPix:
    """GO5e sub-pixel costs on the same predicted coverage points as ``Coverage``. Per spot a map
    over its window (frames, rows, columns), read by trilinear interpolation at each point;
    mass-weighted mean per spot, mean over spots. ``sd``: signed distance to the lit set (pixel
    / frame units, negative inside, boundary at 0, clipped to [-3, 5]). ``int``: 1 - clean
    intensity / the spot's top intensity (intensity 0 off the lit set)."""

    def __init__(self, prob: Any, kind: str) -> None:
        from scipy import ndimage

        self.prob, self.kind = prob, kind
        self.maps, self.lo = [], []
        for sp in prob.region.spots:
            shape = (sp.f1 - sp.f0, sp.r1 - sp.r0, sp.c1 - sp.c0)
            o = np.floor(np.asarray(sp.obs)).astype(np.int64)
            fi = np.clip(o[:, 2] - sp.f0, 0, shape[0] - 1)
            ri = np.clip(o[:, 1] - sp.r0, 0, shape[1] - 1)
            ci = np.clip(o[:, 0] - sp.c0, 0, shape[2] - 1)
            if kind == "sd":
                lit = np.zeros(shape, bool)
                lit[fi, ri, ci] = True
                if lit.any():
                    d_out = ndimage.distance_transform_edt(~lit)
                    d_in = ndimage.distance_transform_edt(lit)
                    m = np.where(lit, -(d_in - 0.5), d_out - 0.5)
                else:
                    m = np.full(shape, 5.0)
                m = np.clip(m, -3.0, 5.0)
            else:
                m = np.zeros(shape)
                v = np.asarray(sp.obs_val, float)
                m[fi, ri, ci] = v / max(float(v.max()), 1e-12) if len(v) else 0.0
            self.maps.append(m)
            self.lo.append(np.array([sp.f0, sp.r0, sp.c0], float) + 0.5)

    def __call__(self, R: torch.Tensor, strain: torch.Tensor) -> float:
        from scipy import ndimage

        prob = self.prob
        prob.n_evals += 1
        with torch.no_grad():
            rows = prob.model.predict(R, prob.V, strain)
            x, a = GF.pack_predicted(
                prob.region, rows, prob.bary, prob.model.n_det, with_shell=False
            )
        x, a = x.numpy(), a.numpy()
        vals = []
        for j, m in enumerate(self.maps):
            ok = a[j] > 0
            if not ok.any():
                continue
            c = x[j][ok][:, ::-1] - self.lo[j]  # (frame, row, col) in map index units
            if self.kind == "sd":
                v = ndimage.map_coordinates(m, c.T, order=1, mode="nearest")
            else:
                v = ndimage.map_coordinates(m, c.T, order=1, mode="constant", cval=0.0)
            w = a[j][ok]
            vals.append(float((w * v).sum() / w.sum()))
        if not vals:
            return 5.0 if self.kind == "sd" else 1.0
        f = float(np.mean(vals))
        return f if self.kind == "sd" else 1.0 - f


def robust_theta0(prob: Any, sym: np.ndarray, tol_deg: float = 1.5) -> Any:
    """Field fit of the start that ignores painted outliers: the medoid voxel (most members within
    ``tol_deg``) is the anchor; members farther than ``tol_deg`` from it are set to it before the
    least-squares field fit (true-grain membership keeps voxels the paint got wrong)."""
    from icenine.ray_indexing import misorientation_deg

    R = prob.R_start
    n = len(R)
    D = np.stack(
        [np.asarray(misorientation_deg(R, np.repeat(R[i : i + 1], n, 0), sym)) for i in range(n)]
    )
    m = int(np.argmax((D <= tol_deg).sum(1)))
    out = D[m] > tol_deg
    Rc = R.copy()
    Rc[out] = R[m]
    return prob.field.fit_theta(Rc, sym), int(out.sum())


SEED = 0  # GO5e: offset added to every CMA seed (--seed)


def cma_min(f: Any, x0: np.ndarray, stds: np.ndarray, evals: int, seed: int = 0) -> Dict[str, Any]:
    import cma

    f0 = f(x0)
    es = cma.CMAEvolutionStrategy(
        list(x0),
        1.0,
        {
            "seed": seed + 1 + 1000 * SEED,
            "maxfevals": evals,
            "verbose": -9,
            "CMA_stds": list(stds),
            "tolfun": 1e-7,
            "tolx": 1e-9,
        },
    )
    best_x, best_f, n = np.asarray(x0, float), f0, 1
    while not es.stop():
        X = es.ask()
        F = [f(np.asarray(x)) for x in X]
        n += len(X)
        es.tell(X, F)
        i = int(np.argmin(F))
        if F[i] < best_f:
            best_f, best_x = F[i], np.asarray(X[i], float)
    return dict(x=best_x, f0=f0, f=best_f, evals=n, stop=str(dict(es.stop())))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--arm",
        required=True,
        choices=["E@GT", "FE@GT", "FE@P", "J@P", "S@P", "SR@P", "S@GT", "E@F", "O@P", "Ost@P"],
    )
    ap.add_argument("--out", default="go5b")
    ap.add_argument("--evals", type=int, default=600)
    ap.add_argument("--evals_e", type=int, default=300)
    ap.add_argument("--evals_f", type=int, default=200)
    ap.add_argument("--only", type=int, default=None)
    ap.add_argument("--cost", default="cov", choices=["cov", "sd", "int"])
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    global SEED
    SEED = a.seed
    out = CM.OUT / a.out
    tag = a.arm.replace("@", "_")
    if a.cost != "cov":
        tag += f"_{a.cost}"
    if a.seed:
        tag += f"_s{a.seed}"
    start = "paint" if a.arm.endswith("@P") else "truth"
    S = CM.Sample(out / f"run_{tag}", start=start, intensity=True, sample="strain")
    T = S.T
    gm = np.asarray(T.gmap)
    label = CM.C.preflight_label(out / f"run_{tag}")
    grains = pick_grains(T) if a.only is None else [a.only]
    rows_out: List[Dict[str, Any]] = []
    for g in grains:
        t0 = time.time()
        members = [int(v) for v in np.nonzero(gm == g)[0]]
        shell = GF.shell_of(members, S.nb)
        n = len(members)
        prob = GF.GrainProblem(
            S.model,
            S.images,
            S.R0,
            S.V,
            S.xy,
            members,
            shell,
            S.sym,
            n_sub=2,
            order=1,
            spec=GF.ObjectiveSpec(kind="sinkhorn"),
            strain=np.zeros((S.n, 6)),
        )
        cov = Coverage(prob) if a.cost == "cov" else SubPix(prob, a.cost)
        zero = torch.zeros((n, 6), dtype=GF.DT)
        # ground truth in the fitted symmetry representative (that of prob.R_start), the same
        # re-expression GrainProblem applies to a strain when it swaps representatives
        ops = GF.align_ops(T.R[members], prob.R_start[0], S.sym)
        R_gt = GF.align_to(T.R[members], prob.R_start[0], S.sym)
        eps_gt = GF.rotate_strain(T.strain[members], ops)[0]
        Rg = torch.as_tensor(R_gt, dtype=GF.DT)

        def eps_of(e6: np.ndarray) -> torch.Tensor:
            return torch.as_tensor(np.tile(e6, (n, 1)), dtype=GF.DT)

        f_gt = cov(Rg, eps_of(eps_gt))
        f_gt_unstr = cov(Rg, zero)
        field = prob.field
        info: Dict[str, Any] = {}
        if a.arm == "E@GT":
            r = cma_min(lambda e: cov(Rg, eps_of(e)), np.zeros(6), np.full(6, STEP_EPS), a.evals_e)
            eps_fit, R_fit = r["x"], R_gt
        else:
            th0 = prob.theta0
            if a.arm.endswith("@P"):
                th0, info["n_outliers"] = robust_theta0(prob, S.sym)
                if a.arm == "Ost@P":  # GO5d: field leg on the strict cost (no strain model)

                    def f_field(th: np.ndarray) -> float:
                        with torch.no_grad():
                            R = field.orientations(torch.as_tensor(th, dtype=GF.DT)).numpy()
                        return float(np.mean(S.strict(R, members)))

                else:

                    def f_field(th: np.ndarray) -> float:
                        return cov(field.orientations(torch.as_tensor(th, dtype=GF.DT)), zero)

                rf = cma_min(
                    f_field,
                    th0,
                    np.full(len(th0), STEP_ROT),
                    a.evals_f,
                )
                th0 = rf["x"]
                info["field"] = {k: rf[k] for k in ("f0", "f", "evals", "stop")}
            k = len(th0)

            def fj(x: np.ndarray) -> float:
                return cov(field.orientations(torch.as_tensor(x[:k], dtype=GF.DT)), eps_of(x[k:]))

            def stds(rot: float, eps: float) -> np.ndarray:
                return np.concatenate([np.full(k, rot), np.full(6, eps)])

            if a.arm in ("E@F", "O@P", "Ost@P"):  # GO5d: strain only at the fixed field
                th_t = torch.as_tensor(th0, dtype=GF.DT)
                rs = cma_min(
                    lambda e: cov(field.orientations(th_t), eps_of(e)),
                    np.zeros(6),
                    np.full(6, STEP_EPS),
                    a.evals_e,
                )
                r = dict(rs, x=np.concatenate([th0, rs["x"]]))
            elif a.arm in ("FE@GT", "FE@P", "J@P"):  # one joint leg
                r = cma_min(
                    fj, np.concatenate([th0, np.zeros(6)]), stds(STEP_ROT, STEP_EPS), a.evals
                )
            else:  # GO5c staged: strain only at the fixed field, then joint legs
                th_t = torch.as_tensor(th0, dtype=GF.DT)
                rs = cma_min(
                    lambda e: cov(field.orientations(th_t), eps_of(e)),
                    np.zeros(6),
                    np.full(6, STEP_EPS),
                    a.evals_e,
                )
                info["strain"] = {k_: rs[k_] for k_ in ("f0", "f", "evals", "stop")}
                rot, ep = math.radians(0.02), 2e-4
                legs = [a.evals] + ([300, 300] if a.arm == "SR@P" else [])
                x = np.concatenate([th0, rs["x"]])
                info["legs"] = []
                for i, ev in enumerate(legs):
                    r = cma_min(fj, x, stds(rot, ep), ev, seed=i)
                    info["legs"].append({k_: r[k_] for k_ in ("f0", "f", "evals", "stop")})
                    x = r["x"]
                    rot, ep = rot / 2, ep / 2
            eps_fit = r["x"][k:]
            with torch.no_grad():
                R_fit = field.orientations(torch.as_tensor(r["x"][:k], dtype=GF.DT)).numpy()
        info["joint"] = {k_: r[k_] for k_ in ("f0", "f", "evals", "stop")}
        eps_fit = np.asarray(eps_fit, float)
        err_vox = T.err(R_fit, members)
        err_start = T.err(prob.R_start, members)
        row = dict(
            arm=a.arm,
            cost=a.cost,
            seed=a.seed,
            grain=int(g),
            level=float(T.strain_level_of[g]),
            n_vox=n,
            n_spots=prob.n_spots,
            eps_gt=eps_gt.tolist(),
            eps_fit=eps_fit.tolist(),
            err_rms=tensor_rms(eps_fit - eps_gt),
            err_comp=(eps_fit - eps_gt).tolist(),
            f_fit=float(r["f"]),
            f_gt=f_gt,
            f_gt_unstrained=f_gt_unstr,
            mis_start_med=float(np.median(err_start)),
            mis_fit_med=float(np.median(err_vox)),
            mis_fit_p95=float(np.percentile(err_vox, 95)),
            evals=prob.n_evals,
            info=info,
            t_s=time.time() - t0,
        )
        rows_out.append(row)
        print(
            json.dumps(
                {
                    k_: row[k_]
                    for k_ in (
                        "grain",
                        "level",
                        "n_vox",
                        "err_rms",
                        "f_fit",
                        "f_gt",
                        "f_gt_unstrained",
                        "mis_start_med",
                        "mis_fit_med",
                        "evals",
                        "t_s",
                    )
                },
                default=float,
            ),
            flush=True,
        )
    name = tag + ("" if a.only is None else f"_g{a.only}")
    (out / f"{name}.json").write_text(
        json.dumps(dict(arm=a.arm, timing=label, rows=rows_out), indent=1, default=float) + "\n"
    )


if __name__ == "__main__":
    main()
