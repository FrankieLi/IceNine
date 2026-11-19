"""GO5b (idea I38): per-grain strain fit on sample ``strain`` with the coverage cost.

  uv run python scripts/grain_opt/fit_strain.py --arm E@GT|FE@GT|FE@P [--out go5b]
      [--evals 600] [--only <grain id>]

Benchmark S12: 4 grains per strain level (1e-4, 5e-4, 1e-3 tensor RMS), sizes at quantiles 0.2,
0.4, 0.6, 0.8 of the level; true grain ids as membership. Cost: 1 - share of predicted coverage
mass on lit pixels, mean over the grain's spots (GO5 "coverage"). Spot windows are built without
strain. Parameters: grain field (FieldA order 1, 9 numbers) and one crystal-frame strain per
grain (6). CMA with per-coordinate steps 0.05 deg (field) and 5e-4 (strain), strain from 0.

GO7d (I53): ``--cen_split``: observed lit pixels nearer the fixed neighbour-shell prediction
than the grain's start prediction are left out of the observed centroid.

GO7c (I52): ``--psf_sigma S``: the coverage cost scores K = 4 Gaussian copies (sigma S px) of each
predicted coverage point (a PSF on the prediction).

GO5n (I49): ``--cen_w bin|binclean``: observed centroids with binary weights (and isolated lit
pixels removed) in every centroid use; ``--stack`` (GO5m) reads another image stack.

GO5l (I47): ``GP4@P`` / ``GP6@P``: ``GFF@P`` plus a fallback when every first-round stage is
skipped, and a joint CMA polish (400 / 600, 0.02 deg and 1e-4) when the strain RMS exceeds 6e-4.

GO5k (I46): ``GF2@P`` (``GA2@P`` with per-stage frozen Tukey weights, a stage skipped below
half inliers) and ``GFF@P`` (plus the GO5h CMA hybrid for grains with fewer than 40 spots).

GO5j (I45): ``GA1@P`` / ``GA2@P``: ``GN@P`` plus one or two rounds of (gradient-only CMA 100 on
coverage at the Gauss-Newton strain, then Gauss-Newton again).

GO5i (I44): arms ``GN@P`` (coverage field leg, then Levenberg-Marquardt on per-spot centroid
residuals over the mean rotation and strain, ``gn_fit``) and ``GNH@P`` (then one joint CMA leg of
300 on the main cost at 0.01 deg and 1e-4).

GO5g (I43): ``--cost hyb --lam L`` (``cen`` + L x coverage); arm ``Sx@P`` (field leg on coverage,
strain-only leg on ``cen``, joint leg on ``hyb`` with L = 2).

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


def pick_grains(T: Any, per_level: int = 4) -> List[int]:
    gm = np.asarray(T.gmap)
    out = []
    for lv in LEVELS:
        gids = [int(g) for g, v in T.strain_level_of.items() if abs(float(v) - lv) < 1e-12]
        gids = sorted(gids, key=lambda g: (int((gm == g).sum()), g))
        qs = QUANTILES if per_level == 4 else [(i + 1) / (per_level + 1) for i in range(per_level)]
        out += [gids[int(q * (len(gids) - 1))] for q in qs]
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
        if PSF_SIGMA > 0:  # GO7c: K Gaussian copies of each point (column, row), fixed seed
            rng = np.random.default_rng(0)
            off = rng.normal(0.0, PSF_SIGMA, size=(PSF_K,) + x.shape[:2] + (2,))
            xk = np.repeat(x[None], PSF_K, 0)
            xk[..., :2] += off
            x = np.concatenate(list(xk), axis=1)
            a = np.concatenate([a / PSF_K] * PSF_K, axis=1)
        p = np.clip(np.floor(x).astype(np.int64), 0, self.W - 1)
        on = np.isin(self._enc(np.arange(self.S)[:, None], p), self.lit)
        tot = a.sum(1)
        ok = tot > 0
        if not ok.any():
            return 1.0
        return 1.0 - float(((a * on).sum(1)[ok] / tot[ok]).mean())


CEN_SPLIT = False  # GO7d: split merged spots by the predicted positions (--cen_split)
PSF_SIGMA, PSF_K = 0.0, 4  # GO7c: PSF-aware coverage (--psf_sigma)
CEN_W = "int"  # GO5n/o: observed-centroid weights, "int" | "bin" | "binclean" | "intclean"


def isolated_mask(ob: np.ndarray) -> np.ndarray:
    """GO5n: True for lit pixels with at least 2 lit 8-neighbours among the spot's pixels in the
    same frame (pixel centres (col, row, frame) + 0.5)."""
    p = np.floor(ob).astype(np.int64)
    lit = {(int(c), int(r), int(f)) for c, r, f in p}
    nb = [(dc, dr) for dc in (-1, 0, 1) for dr in (-1, 0, 1) if dc or dr]
    cnt = np.array([sum((c + dc, r + dr, f) in lit for dc, dr in nb) for c, r, f in p])
    return cnt >= 2


class SubPix:
    """GO5e sub-pixel costs on the same predicted coverage points as ``Coverage``. Per spot a map
    over its window (frames, rows, columns), read by trilinear interpolation at each point;
    mass-weighted mean per spot, mean over spots. ``sd``: signed distance to the lit set (pixel
    / frame units, negative inside, boundary at 0, clipped to [-3, 5]). ``int``: 1 - clean
    intensity / the spot's top intensity (intensity 0 off the lit set)."""

    def __init__(self, prob: Any, kind: str) -> None:
        from scipy import ndimage

        self.prob, self.kind = prob, kind
        self.maps, self.lo, self.obs, self.w = [], [], [], []
        self.n_split = 0
        gpts = None
        if kind == "cen" and CEN_SPLIT:  # GO7d: grain prediction at the start, no strain
            with torch.no_grad():
                rows = prob.model.predict(torch.as_tensor(prob.R_start, dtype=GF.DT), prob.V)
                xg, ag = GF.pack_predicted(
                    prob.region, rows, prob.bary, prob.model.n_det, with_shell=False
                )
            xg, ag = xg.numpy(), ag.numpy()
            gpts = [xg[j][ag[j] > 0] for j in range(len(prob.region.spots))]
        for sp in prob.region.spots:
            shape = (sp.f1 - sp.f0, sp.r1 - sp.r0, sp.c1 - sp.c0)
            o = np.floor(np.asarray(sp.obs)).astype(np.int64)
            fi = np.clip(o[:, 2] - sp.f0, 0, shape[0] - 1)
            ri = np.clip(o[:, 1] - sp.r0, 0, shape[1] - 1)
            ci = np.clip(o[:, 0] - sp.c0, 0, shape[2] - 1)
            w = np.asarray(sp.obs_val, float)
            ob = np.asarray(sp.obs, float)
            if kind == "cen" and CEN_W != "int":  # GO5n/GO5o: noise-robust observed centroid
                if CEN_W in ("bin", "binclean"):
                    w = np.ones(len(ob))
                if CEN_W in ("binclean", "intclean") and len(ob):
                    keep = isolated_mask(ob)
                    if keep.any():
                        ob, w = ob[keep], w[keep]
            if gpts is not None and len(ob):  # GO7d: drop pixels nearer a neighbour's prediction
                from scipy.spatial import cKDTree

                jj = len(self.obs)
                sh = np.asarray(prob.region.shell_pts[jj], float)
                if len(sh) and len(gpts[jj]):
                    dg, _ = cKDTree(gpts[jj]).query(ob)
                    ds, _ = cKDTree(sh).query(ob)
                    keep = dg <= ds
                    if keep.any() and not keep.all():
                        self.n_split += int((~keep).sum())
                        ob, w = ob[keep], w[keep]
            self.obs.append(ob)
            self.w.append(w / max(float(w.sum()), 1e-12))
            if kind in ("sd", "ch"):
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
            w = a[j][ok]
            if self.kind == "cen":  # GO5f: predicted mass centroid vs observed intensity centroid
                cp = (w[:, None] * x[j][ok]).sum(0) / w.sum()
                co = (self.w[j][:, None] * self.obs[j]).sum(0)
                vals.append(min(float(((cp - co) ** 2).sum()), 4.0))
                continue
            c = x[j][ok][:, ::-1] - self.lo[j]  # (frame, row, col) in map index units
            if self.kind == "ch":  # GO5f: predicted mass off the lit set + lit pixels uncovered
                from scipy.spatial import cKDTree

                v = np.clip(ndimage.map_coordinates(m, c.T, order=1, mode="nearest"), 0.0, 5.0)
                t1 = float((w * v).sum() / w.sum())
                if len(self.obs[j]):
                    dist, _ = cKDTree(x[j][ok]).query(self.obs[j])
                    t1 += float((self.w[j] * np.minimum(dist, 5.0)).sum())
                vals.append(t1)
                continue
            if self.kind == "sd":
                v = ndimage.map_coordinates(m, c.T, order=1, mode="nearest")
            else:
                v = ndimage.map_coordinates(m, c.T, order=1, mode="constant", cval=0.0)
            vals.append(float((w * v).sum() / w.sum()))
        if not vals:
            return 5.0 if self.kind == "sd" else 1.0
        f = float(np.mean(vals))
        return 1.0 - f if self.kind == "int" else f


def centroid_residuals(cen: Any, R: torch.Tensor, strain: torch.Tensor) -> Any:
    """GO5i: per-spot (predicted mass centroid - observed intensity centroid), (S, 3) in column,
    row and frame units, and the mask of spots that have predicted mass and observed pixels."""
    prob = cen.prob
    prob.n_evals += 1
    with torch.no_grad():
        rows = prob.model.predict(R, prob.V, strain)
        x, a = GF.pack_predicted(prob.region, rows, prob.bary, prob.model.n_det, with_shell=False)
    x, a = x.numpy(), a.numpy()
    S = len(cen.obs)
    res, has = np.zeros((S, 3)), np.zeros(S, bool)
    for j in range(S):
        ok = a[j] > 0
        if not ok.any() or not len(cen.obs[j]):
            continue
        w = a[j][ok]
        res[j] = (w[:, None] * x[j][ok]).sum(0) / w.sum() - (cen.w[j][:, None] * cen.obs[j]).sum(0)
        has[j] = True
    return res, has


def robust_w(d: np.ndarray, c: float) -> np.ndarray:
    """Per-spot weights: Huber at 2 px when c == 0, else Tukey biweight with cut-off c px."""
    if c <= 0:
        return np.where(d <= 2.0, 1.0, 2.0 / np.maximum(d, 1e-12))
    return np.where(d < c, (1 - (d / c) ** 2) ** 2, 0.0)


def gn_fit(
    cen: Any,
    field: Any,
    th: np.ndarray,
    eps0: np.ndarray,
    n: int,
    iters: int = 15,
    c: float = 0.0,
    freeze: bool = False,
    min_share: float = 0.0,
) -> Any:
    """GO5i: Levenberg-Marquardt on the centroid residuals over the field's mean rotation
    (theta[:3]) and one strain (6); the field gradient theta[3:] stays fixed. Finite-difference
    Jacobian, Huber weights per spot (1 up to 2 px, 2/|d| beyond), Marquardt scaling."""
    h = np.full(9, 1e-5)

    def resid(q: np.ndarray) -> Any:
        t = np.array(th, float)
        t[:3] = q[:3]
        R = field.orientations(torch.as_tensor(t, dtype=GF.DT))
        e = torch.as_tensor(np.tile(q[3:], (n, 1)), dtype=GF.DT)
        return centroid_residuals(cen, R, e)

    q = np.concatenate([np.asarray(th[:3], float), np.asarray(eps0, float)])
    res, active = resid(q)
    w0 = robust_w(np.linalg.norm(res, axis=1), c)  # GO5k: frozen weights of the stage start
    n_in = int(((w0 > 0) & active).sum())
    if freeze and n_in < min_share * active.sum():  # GO5k: too few inliers, skip the stage
        return q, dict(iters=0, skipped=True, n_in=n_in, n_active=int(active.sum()))
    mu, it, c_start = 1e-3, 0, None
    for it in range(1, iters + 1):
        d = np.linalg.norm(res, axis=1)
        sw = np.sqrt(w0 if freeze else robust_w(d, c))[:, None]

        def rv(r: np.ndarray) -> np.ndarray:
            return (sw * r)[active].ravel()

        r0 = rv(res)
        c0 = float(r0 @ r0)
        c_start = c0 if c_start is None else c_start
        J = np.empty((len(r0), 9))
        for i in range(9):
            qq = q.copy()
            qq[i] += h[i]
            J[:, i] = (rv(resid(qq)[0]) - r0) / h[i]
        A, g = J.T @ J, J.T @ r0
        D = np.diag(np.maximum(np.diag(A), 1e-30))
        accepted = False
        for _ in range(6):
            delta = -np.linalg.solve(A + mu * D, g)
            res_n, _h = resid(q + delta)
            rn = rv(res_n)
            if float(rn @ rn) < c0:
                q, res, accepted = q + delta, res_n, True
                mu = max(mu / 3.0, 1e-9)
                break
            mu *= 10.0
        if not accepted or np.all(np.abs(delta) < 1e-7):
            break
    d = np.linalg.norm(res, axis=1)
    w = w0 if freeze else robust_w(d, c)
    c_end = float((w * d**2)[active].sum())
    return q, dict(
        iters=it, cost0=c_start, cost=c_end, n_active=int(active.sum()), n_in=n_in, mu=mu
    )


class Hybrid:
    """GO5g: ``cen`` + lam x coverage (centroids for strain, coverage for the orientation field)."""

    def __init__(self, prob: Any, lam: float) -> None:
        self.prob, self.lam = prob, lam
        self.cov, self.cen = Coverage(prob), SubPix(prob, "cen")

    def __call__(self, R: torch.Tensor, strain: torch.Tensor) -> float:
        f = self.cen(R, strain) + self.lam * self.cov(R, strain)
        self.prob.n_evals -= 1  # one hybrid evaluation, two predictions
        return f


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


GN_ARMS = ("GN@P", "GNH@P", "GA1@P", "GA2@P", "GF2@P", "GFF@P", "GP4@P", "GP6@P")
GP_ARMS = {"GP4@P": 400, "GP6@P": 600}  # GO5l: polish budget
SEED = 0  # GO5e: offset added to every CMA seed (--seed)


def grad_errors(fit_vox: np.ndarray, gt_vox: np.ndarray) -> Dict[str, float]:
    """GO5q: error of the intragranular part (deviation from the grain mean): RMS over voxels of
    the tensor RMS of fitted minus ground-truth deviation, absolute and relative to the
    ground-truth deviation RMS; and the per-voxel strain error RMS."""
    fd = fit_vox - fit_vox.mean(0)
    gd = gt_vox - gt_vox.mean(0)

    def rms(x: np.ndarray) -> float:
        return float(np.sqrt(np.mean([tensor_rms(r) ** 2 for r in x])))

    g = rms(gd)
    return dict(
        grad_err=rms(fd - gd),
        grad_gt=g,
        grad_rel=rms(fd - gd) / g if g > 0 else 0.0,
        vox_err=rms(fit_vox - gt_vox),
    )


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
        choices=[
            "E@GT",
            "FE@GT",
            "FE@P",
            "J@P",
            "S@P",
            "SR@P",
            "S@GT",
            "E@F",
            "O@P",
            "Ost@P",
            "Sx@P",
            "GN@P",
            "GNH@P",
            "GA1@P",
            "GA2@P",
            "GF2@P",
            "GFF@P",
            "GP4@P",
            "GP6@P",
        ],
    )
    ap.add_argument("--out", default="go5b")
    ap.add_argument("--evals", type=int, default=600)
    ap.add_argument("--evals_e", type=int, default=300)
    ap.add_argument("--evals_f", type=int, default=200)
    ap.add_argument("--only", type=int, default=None)
    ap.add_argument("--cost", default="cov", choices=["cov", "sd", "int", "cen", "ch", "hyb"])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--lam", type=float, default=2.0)
    ap.add_argument("--per_level", type=int, default=4)
    ap.add_argument("--psf_sigma", type=float, default=0.0)
    ap.add_argument("--cen_split", action="store_true")
    ap.add_argument("--ratio_hi", type=float, default=None)
    ap.add_argument("--grad", action="store_true", help="GO5q: 6 + 12 strain numbers per grain")
    ap.add_argument("--evals_g", type=int, default=600)
    ap.add_argument("--sample", default="strain")
    ap.add_argument("--cen_w", default="int", choices=["int", "bin", "binclean", "intclean"])
    ap.add_argument("--stack", default=None, help="GO5m: stack file in place of strain_clean.npy")
    a = ap.parse_args()
    global SEED, CEN_W, PSF_SIGMA, CEN_SPLIT
    SEED = a.seed
    CEN_W = a.cen_w
    PSF_SIGMA = a.psf_sigma
    CEN_SPLIT = a.cen_split
    out = CM.OUT / a.out
    tag = a.arm.replace("@", "_")
    if a.cost != "cov":
        tag += f"_{a.cost}" + (f"{a.lam:g}" if a.cost == "hyb" else "")
    if a.seed:
        tag += f"_s{a.seed}"
    if a.per_level != 4:
        tag += f"_S{3 * a.per_level}"
    if a.cen_w != "int":
        tag += f"_{a.cen_w}"
    if a.psf_sigma > 0:
        tag += f"_psf{a.psf_sigma:g}"
    if a.cen_split:
        tag += "_split"
    if a.ratio_hi is not None:
        tag += f"_r{a.ratio_hi:g}"
    if a.grad:
        tag += "_grad"
    start = "paint" if a.arm.endswith("@P") else "truth"
    S = CM.Sample(out / f"run_{tag}", start=start, intensity=True, sample=a.sample, stack=a.stack)
    T = S.T
    gm = np.asarray(T.gmap)
    label = CM.C.preflight_label(out / f"run_{tag}")
    grains = pick_grains(T, a.per_level) if a.only is None else [a.only]
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
            region_kw=None if a.ratio_hi is None else dict(ratio=(0.8, a.ratio_hi)),
        )
        if a.cost == "hyb":
            cov: Any = Hybrid(prob, a.lam)
        else:
            cov = Coverage(prob) if a.cost == "cov" else SubPix(prob, a.cost)
        cov_f = cov_e = cov_j = cov  # field leg, strain-only leg, joint legs
        if a.arm == "Sx@P":  # GO5g by role: field on coverage, strain on centroids, joint hybrid
            cov_f, cov_e, cov_j = Coverage(prob), SubPix(prob, "cen"), Hybrid(prob, 2.0)
        if a.arm in GN_ARMS:  # GO5i: field leg on coverage, then Gauss-Newton
            cov_f = Coverage(prob)
        zero = torch.zeros((n, 6), dtype=GF.DT)
        # ground truth in the fitted symmetry representative (that of prob.R_start), the same
        # re-expression GrainProblem applies to a strain when it swaps representatives
        ops = GF.align_ops(T.R[members], prob.R_start[0], S.sym)
        R_gt = GF.align_to(T.R[members], prob.R_start[0], S.sym)
        eps_gt_vox = GF.rotate_strain(T.strain[members], ops)
        eps_gt = eps_gt_vox.mean(0) if a.grad else eps_gt_vox[0]
        Rg = torch.as_tensor(R_gt, dtype=GF.DT)
        xy_m = np.asarray(S.xy[members], float)  # GO5q: normalised in-plane voxel offsets
        d_m = xy_m - xy_m.mean(0)
        d_m = d_m / max(float(np.sqrt((d_m**2).sum(1).mean())), 1e-12)
        NP = 18 if a.grad else 6

        def eps_np(e: np.ndarray) -> np.ndarray:
            e = np.asarray(e, float)
            if e.ndim == 2:
                return e
            if len(e) == 18:  # mean + in-plane linear gradient D (6 x 2)
                return e[:6] + d_m @ e[6:].reshape(6, 2).T
            return np.tile(e, (n, 1))

        def eps_of(e6: np.ndarray) -> torch.Tensor:
            return torch.as_tensor(eps_np(e6), dtype=GF.DT)

        f_gt = cov(Rg, eps_of(eps_gt_vox))
        f_gt_unstr = cov(Rg, zero)
        field = prob.field
        info: Dict[str, Any] = {}
        if a.arm == "E@GT":
            r = cma_min(
                lambda e: cov(Rg, eps_of(e)), np.zeros(NP), np.full(NP, STEP_EPS), a.evals_e
            )
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
                        return cov_f(field.orientations(torch.as_tensor(th, dtype=GF.DT)), zero)

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
                return cov_j(field.orientations(torch.as_tensor(x[:k], dtype=GF.DT)), eps_of(x[k:]))

            def stds(rot: float, eps: float, ne: int = 6) -> np.ndarray:
                return np.concatenate([np.full(k, rot), np.full(ne, eps)])

            if a.arm in GN_ARMS:  # GO5i: Gauss-Newton on centroid residuals (GO5j: + rounds)
                cen_o, th_cur = SubPix(prob, "cen"), np.array(th0, float)
                q, info["gn"], info["grad"] = np.concatenate([th0[:3], np.zeros(6)]), [], []
                frz = a.arm in ("GF2@P", "GFF@P") or a.arm in GP_ARMS  # GO5k: frozen weights
                th0_t = torch.as_tensor(th0, dtype=GF.DT)
                n_act0 = int(centroid_residuals(cen_o, field.orientations(th0_t), zero)[1].sum())
                info["n_active0"] = n_act0
                fallback = (a.arm == "GFF@P" or a.arm in GP_ARMS) and n_act0 < 40
                info["fallback"] = fallback
                n_rounds = (
                    0
                    if fallback
                    else {
                        "GA1@P": 1,
                        "GA2@P": 2,
                        "GF2@P": 2,
                        "GFF@P": 2,
                        "GP4@P": 2,
                        "GP6@P": 2,
                    }.get(a.arm, 0)
                )
                for rd in range(n_rounds + (0 if fallback else 1)):
                    if rd > 0:  # GO5j: gradient refit on coverage at the Gauss-Newton strain
                        eps_t, base = eps_of(q[3:]), th_cur.copy()

                        def fgrad(gv: np.ndarray, base: Any = base, eps_t: Any = eps_t) -> float:
                            t = base.copy()
                            t[3:] = gv
                            return cov_f(field.orientations(torch.as_tensor(t, dtype=GF.DT)), eps_t)

                        rg = cma_min(fgrad, base[3:], np.full(k - 3, math.radians(0.01)), 100)
                        th_cur[3:] = rg["x"]
                        info["grad"].append({k_: rg[k_] for k_ in ("f0", "f", "evals", "stop")})
                    for c_px in (3.0, 1.5, 1.0):  # Tukey cut-off annealed (GO5i amendment)
                        th_c = th_cur.copy()
                        th_c[:3] = q[:3]
                        q, gi = gn_fit(
                            cen_o, field, th_c, q[3:], n, c=c_px, freeze=frz, min_share=0.5
                        )
                        info["gn"].append(dict(gi, c=c_px, round=rd))
                    th_cur[:3] = q[:3]
                    if (
                        rd == 0
                        and a.arm in GP_ARMS
                        and all(g_.get("skipped") for g_ in info["gn"] if g_["round"] == 0)
                    ):  # GO5l: Gauss-Newton stalled, fall back
                        fallback = info["fallback"] = info["stall"] = True
                        break
                xg = np.concatenate([th_cur, q[3:]])
                if fallback:  # GO5k: GO5h CMA hybrid protocol from the coverage field leg
                    rs = cma_min(
                        lambda e: cov(field.orientations(th0_t), eps_of(e)),
                        np.zeros(6),
                        np.full(6, STEP_EPS),
                        a.evals_e,
                    )
                    r = cma_min(
                        fj,
                        np.concatenate([th0, rs["x"]]),
                        stds(math.radians(0.02), 2e-4),
                        a.evals,
                    )
                elif a.arm in GP_ARMS and tensor_rms(q[3:]) > 6e-4:  # GO5l: polish
                    info["polish"] = True
                    r = cma_min(fj, xg, stds(math.radians(0.02), 1e-4), GP_ARMS[a.arm])
                elif a.arm == "GNH@P":
                    r = cma_min(fj, xg, stds(math.radians(0.01), 1e-4), 300)
                else:
                    fg = fj(xg)
                    r = dict(x=xg, f0=fg, f=fg, evals=0, stop="gn")
                if a.grad:  # GO5q: strain leg over mean + gradient at the fitted field
                    th_f = torch.as_tensor(r["x"][:k], dtype=GF.DT)
                    e0 = np.concatenate([r["x"][k : k + 6], np.zeros(12)])
                    rg_ = cma_min(
                        lambda e: cov(field.orientations(th_f), eps_of(e)),
                        e0,
                        np.full(18, 2e-4),
                        a.evals_g,
                    )
                    info["grad_leg"] = {k_: rg_[k_] for k_ in ("f0", "f", "evals", "stop")}
                    r = dict(rg_, x=np.concatenate([r["x"][:k], rg_["x"]]))
            elif a.arm in ("E@F", "O@P", "Ost@P"):  # GO5d: strain only at the fixed field
                th_t = torch.as_tensor(th0, dtype=GF.DT)
                rs = cma_min(
                    lambda e: cov_e(field.orientations(th_t), eps_of(e)),
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
                    lambda e: cov_e(field.orientations(th_t), eps_of(e)),
                    np.zeros(NP),
                    np.full(NP, STEP_EPS),
                    a.evals_e,
                )
                info["strain"] = {k_: rs[k_] for k_ in ("f0", "f", "evals", "stop")}
                rot, ep = math.radians(0.02), 2e-4
                legs = [a.evals] + ([300, 300] if a.arm == "SR@P" else [])
                x = np.concatenate([th0, rs["x"]])
                info["legs"] = []
                for i, ev in enumerate(legs):
                    r = cma_min(fj, x, stds(rot, ep, NP), ev, seed=i)
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
            lam=a.lam if a.cost == "hyb" else None,
            cen_w=a.cen_w,
            psf_sigma=a.psf_sigma,
            cen_split=a.cen_split,
            ratio_hi=a.ratio_hi,
            grad=a.grad,
            sample=a.sample,
            stack=a.stack,
            seed=a.seed,
            grain=int(g),
            level=float(T.strain_level_of[g]),
            n_vox=n,
            n_spots=prob.n_spots,
            eps_gt=eps_gt.tolist(),
            eps_fit=eps_fit.tolist(),
            err_rms=tensor_rms(eps_fit[:6] - eps_gt),
            err_comp=(eps_fit[:6] - eps_gt).tolist(),
            **grad_errors(eps_np(eps_fit), eps_gt_vox),
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
                        "grad_rel",
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
