"""GO5 (ideas I16, I30): how much each strain component matters to the cost, at the ground truth.

  uv run python scripts/grain_opt/strain_sensitivity.py [--sample strain] [--out go5]

Five grains of the ``strain`` sample (one at 1e-4, two at 5e-4, two at 1e-3 tensor RMS; middle
sizes). Each grain is evaluated at its ground-truth orientations and ground-truth strain, with
one strain component moved by delta at a time (all voxels of the grain together, six components
[e11 e22 e33 e23 e13 e12]). Two costs, both through ``SpotModel`` with strain:

* ``coverage``: 1 - (fraction of predicted coverage mass on lit pixels), mean over the grain's
  spots (the strict-like comparison: ``VoxelCostFunction`` has no strain). A point counts as lit
  when its pixel (floor of column, row, omega bin) is lit in the clean image of the grain's
  spot window.
* ``sinkhorn``: the intensity-weighted Sinkhorn divergence of ``GrainProblem`` (reach 3 px,
  blur 1 px), mean over spots, with ``strain`` passed to the spot model and the region.

Per grain, component and cost: cost at delta = 0 (``cost0``), the symmetrised curve
(f(+d) + f(-d)) / 2 on a grid of d, the second difference at d = 2e-4 (curvature), and
``step10``: the smallest d at which the symmetrised curve reaches 1.1 x cost0 (linear
interpolation between grid points; null if not reached up to the largest d). Also the mean spot
shift in pixels per 1e-3 strain of each component (central difference of the predicted
vertices), and the strain that shifts a spot by one pixel, to compare with pixel / distance
(4.4e-4 and 2.7e-4 for the two detectors). Writes ``benchmarks/grain_opt/<out>/sensitivity.json``.
"""

import argparse
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import common as CM  # noqa: I001
import numpy as np
import torch

from icenine import grain_fit as GF

COMPONENTS = ("e11", "e22", "e33", "e23", "e13", "e12")
GRID = np.array([2e-5, 5e-5, 1e-4, 2e-4, 3e-4, 5e-4, 8e-4, 1.2e-3, 2e-3, 3e-3, 5e-3])
H_CURV = 2e-4
# (strain level, size quantile within the level) of the five grains
PICKS: List[Tuple[float, float]] = [(1e-4, 0.5), (5e-4, 0.3), (5e-4, 0.7), (1e-3, 0.4), (1e-3, 0.6)]


def pick_grains(T: Any) -> List[int]:
    gm = np.asarray(T.gmap)
    levels: Dict[float, List[int]] = {}
    for g, lv in T.strain_level_of.items():
        levels.setdefault(float(lv), []).append(int(g))
    out = []
    for lv, q in PICKS:
        gids = sorted(levels[lv], key=lambda g: (int((gm == g).sum()), g))
        out.append(gids[int(q * (len(gids) - 1))])
    return out


def coverage_cost(prob: Any, rows: Any) -> float:
    """1 - mean over spots of (predicted coverage mass on lit pixels / predicted mass)."""
    x, a = GF.pack_predicted(prob.region, rows, prob.bary, prob.model.n_det, with_shell=False)
    fr = []
    for j, sp in enumerate(prob.region.spots):
        lit = {tuple(p) for p in np.floor(np.asarray(sp.obs)).astype(np.int64).tolist()}
        m = a[j].numpy()
        pts = np.floor(x[j].numpy()).astype(np.int64)
        on = np.array([tuple(p) in lit for p in pts.tolist()], dtype=float)
        tot = m.sum()
        if tot > 0:
            fr.append(float((m * on).sum() / tot))
    return 1.0 - float(np.mean(fr))


def make_costs(S: Any, prob: Any, R: torch.Tensor, V: torch.Tensor, base: np.ndarray) -> Any:
    n_all = base.shape[0]
    members = prob.members

    def at(delta: float, comp: int) -> Dict[str, float]:
        full = np.zeros((n_all, 6))
        full[:] = base
        full[members, comp] += delta
        st = torch.as_tensor(full[members], dtype=GF.DT)
        prob.strain = st
        with torch.no_grad():
            rows = S.model.predict(R, V, st)
            cov = coverage_cost(prob, rows)
            sk = float(prob.objective(R))
        return dict(coverage=cov, sinkhorn=sk)

    return at


def shift_px_per_1e3(S: Any, R: torch.Tensor, V: torch.Tensor, base_m: torch.Tensor, comp: int
                     ) -> Dict[str, float]:
    """Mean |shift| (pixels) of the predicted vertices per 1e-3 of one strain component, per
    detector (central difference at 1e-5, rows present in both)."""
    h = 1e-5
    sel = torch.zeros(6, dtype=GF.DT)
    sel[comp] = h
    with torch.no_grad():
        rp = S.model.predict(R, V, base_m + sel)
        rm = S.model.predict(R, V, base_m - sel)
    assert torch.equal(rp.refl, rm.refl) and torch.equal(rp.vox, rm.vox)
    out: Dict[str, float] = {}
    for det in range(S.model.n_det):
        d = torch.sqrt((rp.col[:, det] - rm.col[:, det]) ** 2 + (rp.row[:, det] - rm.row[:, det]) ** 2)
        ok = (rp.hit[:, det].all(1)) & (rm.hit[:, det].all(1))
        out[f"det{det}"] = float(d[ok].mean() / (2 * h) * 1e-3)
    return out


def step_for_rise(curve: np.ndarray, cost0: float, frac: float = 0.10) -> Any:
    tgt = cost0 * (1.0 + frac)
    for i, v in enumerate(curve):
        if v >= tgt:
            if i == 0:
                return float(GRID[0])
            lo, hi = curve[i - 1], curve[i]
            w = (tgt - lo) / (hi - lo) if hi > lo else 1.0
            return float(GRID[i - 1] + w * (GRID[i] - GRID[i - 1]))
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default="strain")
    ap.add_argument("--out", default="go5")
    a = ap.parse_args()
    out = CM.OUT / a.out
    S = CM.Sample(out / "run", start="truth", intensity=True, sample=a.sample)
    T = S.T
    assert T.strain is not None
    gm = np.asarray(T.gmap)
    spec = GF.ObjectiveSpec(kind="sinkhorn", intensity=True, reach_px=3.0, blur_px=1.0)
    res: List[Dict[str, Any]] = []
    for g in pick_grains(T):
        t0 = time.time()
        members = [int(v) for v in np.nonzero(gm == g)[0]]
        shell = GF.shell_of(members, S.nb)
        prob = GF.GrainProblem(
            S.model, S.images, S.T.R, S.V, S.xy, members, shell, S.sym,
            n_sub=2, order=1, spec=spec, strain=T.strain,
        )
        R = torch.as_tensor(prob.R_start, dtype=GF.DT)
        # the strain is in the crystal frame of the ground-truth orientations: they must not
        # have been moved by a symmetry operator
        assert np.abs(prob.R_start - T.R[members]).max() < 1e-6
        V = prob.V
        base = T.strain
        at = make_costs(S, prob, R, V, base)
        zero = np.zeros_like(base)
        prob.strain = torch.as_tensor(zero[members], dtype=GF.DT)
        with torch.no_grad():
            unstr = dict(
                coverage=coverage_cost(prob, S.model.predict(R, V, None)),
                sinkhorn=float(prob.objective(R)),
            )
        c0 = at(0.0, 0)
        eps_g = base[members[0]]
        rec: Dict[str, Any] = dict(
            grain=int(g), strain_level=float(np.sqrt((eps_g[:3] ** 2).sum() + 2 * (eps_g[3:] ** 2).sum()) / np.sqrt(3)),
            eps=eps_g.tolist(), n_vox=len(members), n_spots=prob.n_spots,
            cost0=c0, cost_unstrained_model=unstr, components={},
        )
        st_gt = torch.as_tensor(base[members], dtype=GF.DT)
        for k, name in enumerate(COMPONENTS):
            curves: Dict[str, List[float]] = {"coverage": [], "sinkhorn": []}
            for d in GRID:
                p, m = at(float(d), k), at(-float(d), k)
                for c in curves:
                    curves[c].append(0.5 * (p[c] + m[c]))
            pm = at(H_CURV, k), at(-H_CURV, k)
            entry: Dict[str, Any] = {}
            for c in curves:
                cv = np.array(curves[c])
                entry[c] = dict(
                    curve_sym=cv.tolist(),
                    curvature=float((pm[0][c] + pm[1][c] - 2 * c0[c]) / H_CURV ** 2),
                    step10=step_for_rise(cv, c0[c]),
                )
            entry["shift_px_per_1e-3"] = shift_px_per_1e3(S, R, V, st_gt, k)
            rec["components"][name] = entry
        rec["t_s"] = time.time() - t0
        res.append(rec)
        print(json.dumps(dict(
            grain=g, level=rec["strain_level"], n_vox=len(members), n_spots=prob.n_spots,
            cost0=c0, unstrained=unstr,
            step10={n: (e["coverage"]["step10"], e["sinkhorn"]["step10"])
                    for n, e in rec["components"].items()},
        ), default=float), flush=True)
    doc = dict(sample=a.sample, grid=GRID.tolist(), h_curv=H_CURV, components=list(COMPONENTS),
               spec=dict(kind="sinkhorn", intensity=True, reach_px=3.0, blur_px=1.0),
               reference_pixel_over_distance=[4.4e-4, 2.7e-4], grains=res,
               timing=CM.C.preflight_label(out / "run"))
    (out / "sensitivity.json").write_text(json.dumps(doc, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
