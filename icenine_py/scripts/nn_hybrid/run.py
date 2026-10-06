#!/usr/bin/env python3
"""
Hybrid network -> FindOptimal on the single-voxel perturbation-sweep cases (Task 1 of
hybrid NN / Q_max-8 proxy / run-time profiling; see MIGRATION_HISTORY.md).

The cases are exactly those of scripts/findoptimal_sweep.py experiment B (50 voxels, 10 radii,
20 directions, clean and realistic ["all"] data, per-voxel images; never full-sample renders).
A finisher is always run on the case's experiment-B image with experiment B's MC seed
(realism_seed + 1000003 vpos + 1009 ri + 7919 j + 17 vi), so its result is paired with
findoptimal_b_raw.

Pipelines (--pipeline):
  H0   FindOptimal alone from the perturbed nominal (re-run for timing / exactness check)
  H1   net x1 -> refine_from_candidates([estimate]) (default box)   [strictly paired with B: the
       image and the net's pass 1 see identical noise]
  H3   net x3 -> refine_from_candidates (default box)
  H3c  net x3 -> FindOptimal with box b = clip(3 sigma_max, default, 2 deg), diameter = 3 b;
       where b is not above the default box the H3 result is copied (H3 cache if present)
  H3m  net x3 -> MC (optimizer_sweep.run_mc_adam, 3500 steps x 2 restarts) with box 3 sigma_max
  HG   Huber GN x3 (c = 1) -> net x3 -> FindOptimal      (radii 1.5, 2, 3, 5 by default)
sigma_max = sqrt(lambda_max(L L^T)) in degrees of the LAST net pass's predicted covariance factor
(pass 1's is stored as sig1 for the post-hoc fall-back rule).
The net-only results N1 / N3 are stored in every output file (R_x1, R_x3) and are checked against
perturbation_sweep_raw by tests/test_nn_hybrid.py.

Usage (from icenine_py/):
  uv run python scripts/nn_hybrid/run.py pilot
  uv run python scripts/nn_hybrid/run.py run --pipeline H3 --workers 10
  uv run python scripts/nn_hybrid/summary.py            (metrics, criteria; see summary.py)
"""

import argparse
import contextlib
import io
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import torch
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(ICENINE_PY / "benchmarks"))

import findoptimal_sweep as fs  # noqa: E402
import optimizer_sweep as osw  # noqa: E402
import perturbation_sweep as ps  # noqa: E402
from stage_timer import StageTimer  # noqa: E402

VARIANTS = ps.VARIANTS  # ["clean", "all"]; "all" is the realistic variant
SWEEP_DIR = ICENINE_PY / "benchmarks" / "toy_orientation_sweep"
OUT_DIR = ICENINE_PY / "benchmarks" / "nn_hybrid"
CACHE_DIR = HERE / "cache"
PIPELINES = ["H0", "H1", "H3", "H3c", "H3m", "HG"]
HG_RADII_IDX = [6, 7, 8, 9]  # 1.5, 2, 3, 5 deg
SIGMA_BOX_FACTOR = 3.0
BOX_MAX_DEG = 2.0
MC_STEPS, MC_RESTARTS, MC_STEP_FRAC = 3500, 2, 0.5  # the optimizer sweep's MC protocol

# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------


def default_box_deg(params: Any) -> float:
    """The search box refine_from_candidates uses by default, max(d/3, 0.2 deg) with
    d = local_grid_radius / 1.5^(max_local_resolution + 1) (0.329 deg for ReconstructQ8)."""
    d = params.local_grid_radius / 1.5 ** (params.max_local_resolution + 1)
    return float(np.degrees(max(d / 3.0, np.radians(0.2))))


def sigma_max_deg(L: np.ndarray) -> np.ndarray:
    """sqrt(lambda_max(L L^T)) (deg) of Cholesky factors (..., 3, 3); NaN where L is NaN."""
    S = np.asarray(L, dtype=np.float64)
    S = S @ np.swapaxes(S, -1, -2)
    lam = np.linalg.eigvalsh(np.nan_to_num(S, nan=0.0))[..., -1]
    out = np.sqrt(np.clip(lam, 0.0, None))
    return np.where(np.isnan(np.asarray(L, dtype=np.float64)).any(axis=(-1, -2)), np.nan, out)


def covariance_box_deg(sigma_max: float, default_deg: float, hi_deg: float = BOX_MAX_DEG) -> float:
    """H3c box b = clip(3 sigma_max, default, hi) in degrees."""
    return float(np.clip(SIGMA_BOX_FACTOR * sigma_max, default_deg, hi_deg))


def box_to_diameter_rad(box_deg: float) -> float:
    """Value of refine_from_candidates(diameter=...) giving search radius box_deg, valid for
    box_deg >= 0.2 deg (final radius = max(diameter / 3, 0.2 deg))."""
    assert box_deg >= 0.2, "refine_from_candidates cannot produce boxes below 0.2 deg"
    return float(3.0 * np.radians(box_deg))


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------

_T = StageTimer()
_NETS: Dict[str, torch.nn.Module] = {}


def init_worker(wargs: Dict[str, Any]) -> None:
    from icenine.cost_functions import VoxelCostFunction
    from icenine.orientation_search import MCOptimizer

    fs.init_worker(wargs)
    assert fs._W is not None
    ctx = fs._W.ctx
    # the checks the plan asks for: what ps.prepare_nominal / render_batch need is on the context
    for attr in ("args", "example_dir", "setup"):
        assert hasattr(ctx, attr), f"worker context lacks {attr}"
    paths = dict((n, p) for n, p in ctx.args.models)
    for name in wargs["use_models"]:
        _NETS[name] = ps.load_model(paths[name])[0]
    _T.count_evaluations(VoxelCostFunction)
    _T.patch(MCOptimizer, "optimize", "find_optimal")
    _T.patch(MCOptimizer, "variance_minimizing_optimize", "variance_min")
    _T.label(fs._W.local_fn, "local")


def net_stage(
    ctx: SimpleNamespace,
    net: torch.nn.Module,
    vidx: int,
    prep1: List[Optional[Any]],
    b1: Dict[str, torch.Tensor],
    ok1: np.ndarray,
    delta0: np.ndarray,
    R_nom0: np.ndarray,
    R_true: np.ndarray,
    draws: List[Tuple[np.ndarray, np.ndarray]],
    sources: List[Any],
    variant: str,
    seed: int,
    timer: StageTimer = _T,
) -> Tuple[np.ndarray, np.ndarray]:
    """perturbation_sweep.sweep_voxel's pass loop for one (voxel, radius, variant, model), copied
    (gating, re-centring at the estimate, re-rendering the same true orientation with the same
    draws and seed, carry-forward when a case stops). Returns the estimate of every pass,
    R (D, P, 3, 3) and the predicted Cholesky factor L (D, P, 3, 3) (deg), NaN for cases that
    never had an estimate; entries of later passes of a stopped case repeat its last estimate."""
    a = ctx.args
    K, P, D = a.frame_half_width, a.passes, len(prep1)
    alive = ok1.copy()
    R_nom, delta_t = R_nom0.copy(), delta0.copy()
    R_out = np.full((D, P, 3, 3), np.nan)
    L_out = np.full((D, P, 3, 3), np.nan)
    est: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    for p_i in range(P):
        if p_i == 0:
            batch = b1
        else:
            preps: List[Optional[Any]] = [None] * D
            for j in np.nonzero(alive)[0]:
                with timer.stage("prepare"):
                    pj, _why = ps.prepare_nominal(ctx, vidx, R_nom[j])
                if pj is None:
                    alive[j] = False
                preps[j] = pj
            with timer.stage("render"):
                batch = ps.render_batch(preps, delta_t, draws, sources, variant, seed, a)
            alive &= ~(alive & (batch["n_present"].numpy() < a.min_present))
        rows = np.nonzero(alive)[0]
        if len(rows):
            with timer.stage("forward"):
                dh, L = ps.run_net(net, batch, rows, K)
            Rn, dt = ps.recentre(R_nom[rows], dh, R_true)
            for k, j in enumerate(rows):
                est[int(j)] = (Rn[k], L[k])
            R_nom[rows], delta_t[rows] = Rn, dt
        for j, (R, L) in est.items():
            R_out[j, p_i], L_out[j, p_i] = R, L
    return R_out, L_out


class VariantBatch(SimpleNamespace):
    """Pass-1 state of one (voxel, radius, variant): see variant_batches."""


def variant_batches(
    ctx: SimpleNamespace,
    vidx: int,
    vpos: int,
    ri: int,
    ref_nroi: np.ndarray,
    ref_fail: np.ndarray,
    variants: Sequence[str],
) -> Iterator[VariantBatch]:
    """The pass-1 batches of experiment B's cases, built as findoptimal_sweep.case_images builds
    them (same assertions that the case is aligned with the network sweep), but yielding the batch
    itself so the net stage can use it."""
    a = ctx.args
    vctx = osw.voxel_context(ctx, vidx)
    D, r = a.n_dirs, a.radii[ri]
    sigma_comp = a.neighbor_sigma_deg / np.sqrt(3.0)
    delta0, draws = ps.case_draws(
        a.sweep_seed, vidx, ri, r, D, len(vctx.sources), sigma_comp, a.neighbor_p
    )
    R_nom0 = ps.perturbed_nominal(vctx.R_true, delta0)
    prep1: List[Optional[Any]] = []
    reason1 = np.zeros(D, dtype=np.int8)
    t0 = time.perf_counter()
    for j in range(D):
        with _T.stage("prepare"):
            p, why = ps.prepare_nominal(ctx, vidx, R_nom0[j])
        prep1.append(p)
        reason1[j] = why
    prep1_s = time.perf_counter() - t0  # shared by the variants
    assert (np.array([p.n if p is not None else 0 for p in prep1]) == ref_nroi).all()
    for vi, variant in enumerate(VARIANTS):
        if variant not in variants:
            continue
        seed = a.realism_seed + 1000003 * vpos + 1009 * ri
        layers: Dict[str, torch.Tensor] = {}
        t0 = time.perf_counter()
        with _T.stage("render"):
            b1 = ps.render_batch(
                prep1, delta0, draws, vctx.sources, variant, seed, a, layers=layers
            )
        render1_s = time.perf_counter() - t0
        have = np.array([p is not None for p in prep1])
        ok1 = have & (b1["n_present"].numpy() >= a.min_present)
        fail = np.where(reason1 > 0, reason1, np.where(ok1, 0, 4)).astype(np.int8)
        assert (fail == ref_fail[:, vi]).all(), "pass-1 failures differ from the sweep"
        yield VariantBatch(
            prep1_s=prep1_s / len(variants), render1_s=render1_s,
            vi=vi, variant=variant, seed=seed, vctx=vctx, delta0=delta0, draws=draws,
            R_nom0=R_nom0, prep1=prep1, b1=b1, layers=layers, have=have, ok1=ok1, fail=fail,
        )  # fmt: skip


def case_keys(ctx: SimpleNamespace, vb: VariantBatch, j: int) -> np.ndarray:
    """Pixel keys of case j's detector images (as findoptimal_sweep.case_images)."""
    a = ctx.args
    edit = None
    if vb.variant != "clean":
        p = vb.prep1[j]
        n = p.n
        edit = osw.realism_edit(
            vb.layers["clean"][j, :n].numpy(),
            vb.layers["dis"][j, :n].numpy(),
            vb.b1["windows"][j, :n].numpy(),
            p.spec,
            p.obs.det_idx.numpy(),
            a.frame_half_width,
            ctx.geo,
        )
    return osw.case_image_keys(ctx, vb.vctx, vb.variant, vb.draws[j], edit)


def b_seed(a: Any, vpos: int, ri: int, j: int, vi: int) -> int:
    return int(a.realism_seed + 1000003 * vpos + 1009 * ri + 7919 * j + 17 * vi)


def refine_fo(
    start: np.ndarray, vctx: SimpleNamespace, seed: int, diameter: Optional[float] = None
) -> Tuple[np.ndarray, float, bool]:
    """refine_from_candidates([start]) exactly as findoptimal_sweep.b_task (optionally with a
    search diameter). Returns (R, cost, converged)."""
    from icenine.orientation_search import MCOptimizer, SearchCandidate

    assert fs._W is not None
    rec, lf = fs._W.rec, fs._W.local_fn
    mc = MCOptimizer(
        cost_fn=lf, voxel_vertices=vctx.vertices, phase_index=vctx.voxel.phase,
        rng=np.random.default_rng(seed),
    )  # fmt: skip
    with contextlib.redirect_stdout(io.StringIO()):
        res = rec.refine_from_candidates(
            [SearchCandidate(orientation=start.astype(np.float32), cost=1.0)],
            vctx.vertices,
            vctx.voxel.phase,
            diameter=diameter,
            local_cost_fn=lf,
            mc_optimizer=mc,
        )
    return np.asarray(res.orientation, dtype=np.float64), float(res.cost), bool(
        rec.last_find_optimal.get("converged", False)
    )


def gn_estimates(
    ctx: SimpleNamespace, vidx: int, vb: VariantBatch
) -> Tuple[np.ndarray, np.ndarray]:
    """Huber (c = 1) Gauss-Newton, up to 3 re-centring passes (optimizer_sweep.run_gn). Returns
    the pass-3 estimate R (D, 3, 3) (the perturbed nominal where GN produced none) and a (D,) mask
    of cases where it did."""
    with _T.stage("gn"):
        g = osw.run_gn(
            ctx, vidx, vb.prep1, vb.b1, vb.ok1, vb.delta0, vb.R_nom0, vb.vctx.R_true, vb.draws,
            vb.vctx.sources, vb.variant, vb.seed, 1.0,
        )  # fmt: skip
    e = np.stack([g["err_x"][:, -1], g["err_y"][:, -1], g["err_z"][:, -1]], axis=-1)
    good = np.isfinite(e).all(axis=-1)
    R = vb.R_nom0.copy()
    if good.any():
        R[good] = (
            Rotation.from_rotvec(np.radians(e[good].astype(np.float64))).as_matrix()
            @ vb.vctx.R_true
        )
    return R, good


def task(item: Tuple[Any, ...]) -> Tuple[int, int, float]:
    """One (voxel, radius) of one pipeline: both variants, all directions. Caches the result."""
    vidx, vpos, ri, path, ref_nroi, ref_fail, spec = item
    assert fs._W is not None
    ctx, lf = fs._W.ctx, fs._W.local_fn
    a = ctx.args
    t_start = time.time()
    pipe, model, variants = spec["pipeline"], spec["model"], spec["variants"]
    net = _NETS[model]
    default_deg = default_box_deg(fs._W.rec.params)
    D, V = a.n_dirs, len(VARIANTS)
    f = lambda *s: np.full((D, V) + s, np.nan)  # noqa: E731
    out = dict(
        R_final=f(3, 3), cost_final=f(), cost_true=f(), t_finish=f(), t_fo=f(), t_vm=f(),
        t_eval=f(), evals=np.zeros((D, V), dtype=np.int32), converged=np.zeros((D, V), bool),
        ran=np.zeros((D, V), bool), fallback=np.zeros((D, V), bool), copied=np.zeros((D, V), bool),
        start=f(3, 3), R_x1=f(3, 3), R_x3=f(3, 3), sig1=f(), sig3=f(), box_deg=f(),
        R_gn=f(3, 3), t_prep=f(), t_render=f(), t_forward=f(), t_gn=f(),
        fail_pass1=np.zeros((D, V), dtype=np.int8),
    )  # fmt: skip
    h3_cache = None
    if pipe == "H3c":
        p3 = Path(str(path).replace("H3c_", "H3_"))
        h3_cache = np.load(p3) if p3.exists() else None
    for vb in variant_batches(ctx, vidx, vpos, ri, ref_nroi, ref_fail, variants):
        vi, vctx = vb.vi, vb.vctx
        out["fail_pass1"][:, vi] = vb.fail
        n_have = max(int(vb.have.sum()), 1)
        snap = _T.snapshot()
        R_gn = None
        if pipe == "HG":
            R_gn, gn_ok = gn_estimates(ctx, vidx, vb)
            out["R_gn"][:, vi] = R_gn
            # net from the GN estimate: re-render the true orientation at that nominal
            dstart = ps.relative_offset_deg(vctx.R_true, R_gn)
            preps2: List[Optional[Any]] = [None] * D
            for j in np.nonzero(vb.ok1)[0]:
                with _T.stage("prepare"):
                    preps2[j], _ = ps.prepare_nominal(ctx, vidx, R_gn[j])
            with _T.stage("render"):
                b2 = ps.render_batch(preps2, dstart, vb.draws, vctx.sources, vb.variant, vb.seed, a)
            ok2 = np.array([p is not None for p in preps2]) & (
                b2["n_present"].numpy() >= a.min_present
            )
            Rn, Ln = net_stage(
                ctx, net, vidx, preps2, b2, ok2, dstart, R_gn, vctx.R_true, vb.draws,
                vctx.sources, vb.variant, vb.seed,
            )  # fmt: skip
            # a case with no net estimate keeps the GN estimate
            for p_i in range(Rn.shape[1]):
                miss = np.isnan(Rn[:, p_i, 0, 0]) & vb.ok1
                Rn[miss, p_i] = R_gn[miss]
            Rx1, Rx3 = Rn[:, 0], Rn[:, -1]
            L1, L3 = Ln[:, 0], Ln[:, -1]
        else:
            Rn, Ln = net_stage(
                ctx, net, vidx, vb.prep1, vb.b1, vb.ok1, vb.delta0, vb.R_nom0, vctx.R_true,
                vb.draws, vctx.sources, vb.variant, vb.seed,
            )  # fmt: skip
            Rx1, Rx3, L1, L3 = Rn[:, 0], Rn[:, -1], Ln[:, 0], Ln[:, -1]
        d = _T.delta_since(snap, "inclusive")
        out["R_x1"][:, vi], out["R_x3"][:, vi] = Rx1, Rx3
        out["sig1"][:, vi], out["sig3"][:, vi] = sigma_max_deg(L1), sigma_max_deg(L3)
        for key, st in (("t_prep", "prepare"), ("t_render", "render"), ("t_forward", "forward"),
                        ("t_gn", "gn")):  # fmt: skip
            extra = {"prepare": vb.prep1_s, "render": vb.render1_s}.get(st, 0.0)
            out[key][vb.have, vi] = (d.get(st, 0.0) + extra) / n_have  # batch time per case
        for j in range(D):
            if not vb.have[j]:
                continue
            R_nom = vb.R_nom0[j]
            usable = bool(vb.ok1[j]) and np.isfinite(Rx3[j]).all()
            fallback = pipe != "H0" and not usable
            start = {
                "H0": R_nom, "H1": Rx1[j], "H3": Rx3[j], "H3c": Rx3[j], "H3m": Rx3[j], "HG": Rx3[j],
            }[pipe]  # fmt: skip
            if fallback:
                start = R_nom
            out["fallback"][j, vi] = fallback
            out["start"][j, vi] = start
            seed = b_seed(a, vpos, ri, j, vi)
            attach = fs.attach_images
            attach(case_keys(ctx, vb, j))
            diameter: Optional[float] = None
            box = np.nan
            copied = False
            if pipe == "H3c" and not fallback:
                box = covariance_box_deg(out["sig3"][j, vi], default_deg)
                if box <= default_deg * (1 + 1e-9):
                    if h3_cache is not None:
                        copied = True
                    else:
                        box = default_deg
                else:
                    diameter = box_to_diameter_rad(box)
            if pipe == "H3m" and not fallback:
                box = SIGMA_BOX_FACTOR * float(out["sig3"][j, vi])
            out["box_deg"][j, vi] = box
            out["cost_true"][j, vi] = lf.evaluate(
                vctx.R_true.astype(np.float32), vctx.vertices, vctx.voxel.phase
            ).cost
            if copied:
                for k in ("R_final", "cost_final", "t_finish", "t_fo", "t_vm", "t_eval", "evals"):
                    out[k][j, vi] = h3_cache[k][j, vi]
                out["converged"][j, vi] = h3_cache["converged"][j, vi]
                out["ran"][j, vi], out["copied"][j, vi] = True, True
                continue
            snap = _T.snapshot()
            n0 = lf.eval_count
            t0 = time.perf_counter()
            if pipe == "H3m" and not fallback:
                # run_one_mc uses the search box 1.5 r: r = box / 1.5 gives a box of 3 sigma_max
                o = osw.run_mc_adam(
                    ctx, vctx, case_keys(ctx, vb, j), start, float(box) / 1.5, ["mc"], seed
                )["mc"]
                R_f = np.asarray(o["R_final"], dtype=np.float64)
                cost_f = float(
                    lf.evaluate(R_f.astype(np.float32), vctx.vertices, vctx.voxel.phase).cost
                )
                conv = False
                out["evals"][j, vi] = int(o["n_evals"])
            else:
                R_f, cost_f, conv = refine_fo(start, vctx, seed, diameter)
                out["evals"][j, vi] = lf.eval_count - n0
            out["t_finish"][j, vi] = time.perf_counter() - t0
            dd = _T.delta_since(snap, "inclusive")  # find_optimal / variance_min / evaluate overlap
            out["t_fo"][j, vi] = dd.get("find_optimal", 0.0)
            out["t_vm"][j, vi] = dd.get("variance_min", 0.0)
            out["t_eval"][j, vi] = dd.get("evaluate", 0.0)
            out["R_final"][j, vi], out["cost_final"][j, vi] = R_f, cost_f
            out["converged"][j, vi], out["ran"][j, vi] = conv, True
    np.savez(path, **out)
    return vidx, ri, time.time() - t_start


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def worker_args(models: Sequence[str]) -> Tuple[Dict[str, Any], Dict[str, np.ndarray]]:
    cfg, sweep = osw.load_sweep_config(SWEEP_DIR / "perturbation_sweep_raw.npz")
    w = dict(cfg)
    w.update(methods=[], opt_dirs=0, mc_steps=MC_STEPS, mc_restarts=MC_RESTARTS,
             mc_step_frac=MC_STEP_FRAC, huber_c=1.0, use_models=list(models))  # fmt: skip
    return w, sweep


def make_items(
    args: argparse.Namespace, sweep: Dict[str, np.ndarray], cache: Path
) -> List[Tuple[Any, ...]]:
    voxels = [int(v) for v in sweep["voxel_indices"]][: args.n_voxels or None]
    radii_idx = args.radii_idx or (HG_RADII_IDX if args.pipeline == "HG" else list(range(10)))
    variants = args.variants
    spec = dict(pipeline=args.pipeline, model=args.model, variants=variants)
    tag = f"{args.pipeline}_{args.model}" + ("" if len(variants) == 2 else "_" + variants[0])
    return [
        (v, vpos, ri, str(cache / f"{tag}_v{v}_r{ri}_d{args.n_dirs}.npz"),
         sweep["n_roi"][vpos, ri], sweep["fail_pass1"][vpos, ri], spec)
        for vpos, v in enumerate(voxels)
        for ri in radii_idx
    ]  # fmt: skip


def run_pool(todo: List[Any], workers: int, wargs: Dict[str, Any], label: str) -> None:
    import multiprocessing as mp

    t0 = time.time()
    with mp.get_context("spawn").Pool(workers, initializer=init_worker, initargs=(wargs,)) as pool:
        for k, (vidx, ri, secs) in enumerate(pool.imap_unordered(task, todo), 1):
            print(
                f"  [{label} {k}/{len(todo)}] voxel {vidx} r#{ri} done in {secs:.0f}s "
                f"({time.time() - t0:.0f}s total)",
                flush=True,
            )


def assemble(items: List[Any], args: argparse.Namespace, out_dir: Path) -> Optional[Path]:
    if any(not Path(it[3]).exists() for it in items):
        print("not all tasks cached; nothing assembled")
        return None
    voxels = sorted({it[0] for it in items})
    ris = sorted({it[2] for it in items})
    by = {(it[0], it[2]): np.load(it[3]) for it in items}
    raw = {k: np.stack([np.stack([by[(v, r)][k] for r in ris]) for v in voxels])
           for k in by[(voxels[0], ris[0])].files}  # fmt: skip
    raw["voxel_indices"], raw["radii_idx"] = np.array(voxels), np.array(ris)
    raw["radii"] = np.array([ps.RADII_DEG[i] for i in ris])
    raw["variants"] = np.array(VARIANTS)
    spec = items[0][6]
    name = out_dir / (
        f"{spec['pipeline']}_{spec['model']}"
        + ("" if len(spec["variants"]) == 2 else "_" + spec["variants"][0])
        + "_raw.npz"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(name, **raw)
    print(f"saved {name}", flush=True)
    return name


def do_run(args: argparse.Namespace, cache: Path, out_dir: Path) -> float:
    cache.mkdir(parents=True, exist_ok=True)
    wargs, sweep = worker_args([args.model])
    items = make_items(args, sweep, cache)
    todo = [it for it in items if not Path(it[3]).exists()]
    print(
        f"{args.pipeline}/{args.model}: {len(items)} tasks ({len(items) - len(todo)} cached), "
        f"{args.workers} workers",
        flush=True,
    )
    t0 = time.time()
    run_pool(todo, args.workers, wargs, args.pipeline)
    wall = time.time() - t0
    print(f"{args.pipeline} finished; wall {wall:.0f}s", flush=True)
    import json

    out_dir.mkdir(parents=True, exist_ok=True)
    wt = out_dir / "wall_times.json"
    d = json.loads(wt.read_text()) if wt.exists() else {}
    d[f"{args.pipeline}_{args.model}" + ("" if len(args.variants) == 2 else "_" + args.variants[0])] = dict(
        wall_s=wall, tasks=len(todo), workers=args.workers)
    wt.write_text(json.dumps(d, indent=1))
    assemble(items, args, out_dir)
    return wall


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "pilot"):
        r = sub.add_parser(name)
        r.add_argument("--pipeline", choices=PIPELINES, default="H3")
        r.add_argument("--model", default="realistic_s0")
        r.add_argument("--variants", nargs="+", default=list(VARIANTS), choices=VARIANTS)
        r.add_argument("--workers", type=int, default=10)
        r.add_argument("--n-voxels", type=int, default=0, help="first N sweep voxels (0 = all)")
        r.add_argument("--n-dirs", type=int, default=20)
        r.add_argument("--radii-idx", type=int, nargs="*", default=[])
        r.add_argument("--cache-dir", default=str(CACHE_DIR))
        r.add_argument("--out-dir", default=str(OUT_DIR))
    return ap


def main() -> None:
    args = build_parser().parse_args()
    if args.cmd == "pilot":
        pilot_pipes = [args.pipeline] if "--pipeline" in sys.argv else PIPELINES
        args.n_voxels = 2
        for p in pilot_pipes:
            args.pipeline = p
            do_run(args, Path(args.cache_dir) / "pilot", Path(args.out_dir) / "pilot")
    else:
        do_run(args, Path(args.cache_dir), Path(args.out_dir))


if __name__ == "__main__":
    main()
