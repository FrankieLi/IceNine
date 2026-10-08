#!/usr/bin/env python3
"""
Phase B3 (finisher/MC study): head-to-head of local finishers on the T5 cases and on the
perturbation-sweep cases, all on the same VoxelCostFunction with exactly counted evaluations.

Cases (per variant clean / realistic; images exactly as findoptimal_sweep.case_images builds them):
  T5   the 300 cases of Phase A / B1 (diagnose.build_items): 200 H3 starts (net x3) and 100 H0
       starts (perturbed nominal), the realistic start used for both variants
  SW   50 sweep voxels x radii {0.05, 0.1, 0.25, 0.5, 1} deg x the first 4 directions valid in
       both variants; start = the perturbed nominal (the H0 start)

Methods (finisher_bench.optimizers; every one starts by evaluating the start):
  mc_deployed  MCOptimizer, box 0.3292 / step0 0.1317 deg, 200 steps, 2 restarts (natural stop)
  mc_april     MCOptimizer as in the April sweep: 3500 steps, 2 restarts, box 1.5 r, step 0.5 box,
               r = the case's sweep radius (told to the method, unlike the finisher)
  mc_local     MC with local restarts (fixed 31 stuck steps, restart near the best at half the step)
  es_box, es_002   (1+1)-ES with a success-rate step rule, step0 = 0.1317 / 0.02 deg
  nm           Nelder-Mead on the rotation vector, simplex edge 0.1317 deg, shrinking restarts
  cma_005, cma_02  local CMA-ES, sigma0 = 0.05 / 0.2 deg
  vm_small     VarianceMinimizing with a quarter box (0.0823 deg)
  gn           centroid Huber (c = 1) Gauss-Newton, 3 re-centring passes on the net's windows
  adam         RiemannianAdamOptimizer (April settings: 100 steps, lr 1e-4, scale 2, 2 restarts)
mc_april, mc_local, es_*, nm, cma_*, vm_small are run once to 10000 evaluations; the result at a
smaller budget is the best orientation after that many evaluations (a budgeted run is a prefix
of the long run: tests/test_finisher_bench.py). gn and adam are run once in their own units.

Usage (from icenine_py/):
  uv run python scripts/finisher_bench/bench.py pilot
  uv run python scripts/finisher_bench/bench.py run --workers 10 [--n-voxels 50]
"""

import argparse
import math
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
for sub in (HERE, HERE.parent / "finisher_diagnosis", HERE.parent / "nn_hybrid", HERE.parent):
    sys.path.insert(0, str(sub))
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(ICENINE_PY / "benchmarks"))

import diagnose as D  # noqa: E402
import findoptimal_sweep as fs  # noqa: E402
import optimizer_sweep as osw  # noqa: E402
import optimizers as O  # noqa: E402
import perturbation_sweep as ps  # noqa: E402
import run as nnrun  # noqa: E402

CACHE_DIR = HERE / "cache"
VARIANTS = ["clean", "all"]  # "all" = realistic
VI_REAL = 1
CKPTS = [100, 250, 500, 1000, 2600, 5000, 10000]
BUDGETS = [250, 1000, 2600, 10000]
MAX_BUDGET = 10000
SW_RADII_IDX = [0, 1, 2, 3, 5]  # 0.05, 0.1, 0.25, 0.5, 1 deg
SW_DIRS = 4
ANYTIME = [
    "mc_deployed", "mc_april", "mc_local", "es_box", "es_002", "nm", "cma_005", "cma_02",
    "vm_small",
]  # fmt: skip
CFG_SEED = {n: 7919 * (i + 1) for i, n in enumerate(ANYTIME)}
ES_SMALL_DEG = 0.02


def make_methods(box: float, step0: float, params: Any, r_told_rad: float) -> Dict[str, Any]:
    """name -> function(cc, R0, seed). box / step0 in rad (the finisher's default box and step)."""
    stuck = max(1, int(2.0 * (box / step0) ** 3))  # the deployed initial min_ergodic (31)
    return {
        "mc_deployed": lambda cc, R0, s: O.mc_plain(
            cc, R0, s, box, step0, params.max_mc_steps, params.successive_restarts,
            params.max_convergence_cost,
        ),
        "mc_april": lambda cc, R0, s: O.mc_plain(
            cc, R0, s, 1.5 * r_told_rad, 0.5 * 1.5 * r_told_rad, 3500, 2, 0.0
        ),
        "mc_local": lambda cc, R0, s: O.mc_local_restarts(
            cc, R0, np.random.default_rng(s), step0, stuck
        ),
        "es_box": lambda cc, R0, s: O.one_plus_one_es(cc, R0, np.random.default_rng(s), step0),
        "es_002": lambda cc, R0, s: O.one_plus_one_es(
            cc, R0, np.random.default_rng(s), math.radians(ES_SMALL_DEG)
        ),
        "nm": lambda cc, R0, s: O.nelder_mead_rot(cc, R0, math.degrees(step0)),
        "cma_005": lambda cc, R0, s: O.cma_local(cc, R0, s, 0.05),
        "cma_02": lambda cc, R0, s: O.cma_local(cc, R0, s, 0.2),
        "vm_small": lambda cc, R0, s: O.variance_min_small_box(cc, R0, s, box / 4.0),
    }  # fmt: skip


# ---------------------------------------------------------------------------
# GN and Adam
# ---------------------------------------------------------------------------


def gn_from_starts(
    ctx: Any, vidx: int, vb: Any, starts: Dict[int, np.ndarray]
) -> Tuple[Dict[int, np.ndarray], Dict[int, Dict[str, float]], float]:
    """Centroid Huber GN (3 re-centring passes, c = 1) from each start, on the net's window
    pipeline: windows of the true orientation rendered at the start (same draws / realism seed as
    the case), exactly as nn_hybrid HG's gn_estimates does from the perturbed nominal. Returns
    ({j: estimate}, {j: info}, batch seconds). The estimate is the start where GN gave none."""
    a = ctx.args
    Dn = a.n_dirs
    t0 = time.perf_counter()
    R_s = vb.R_nom0.copy()
    for j, s in starts.items():
        R_s[j] = s
    delta_s = ps.relative_offset_deg(vb.vctx.R_true, R_s)
    preps: List[Any] = [None] * Dn
    for j in starts:
        preps[j], _why = ps.prepare_nominal(ctx, vidx, R_s[j])
    batch = ps.render_batch(preps, delta_s, vb.draws, vb.vctx.sources, vb.variant, vb.seed, a)
    ok = np.array([p is not None for p in preps]) & (batch["n_present"].numpy() >= a.min_present)
    g = osw.run_gn(
        ctx, vidx, preps, batch, ok, delta_s, R_s, vb.vctx.R_true, vb.draws, vb.vctx.sources,
        vb.variant, vb.seed, 1.0,
    )  # fmt: skip
    from scipy.spatial.transform import Rotation

    e = np.stack([g["err_x"][:, -1], g["err_y"][:, -1], g["err_z"][:, -1]], axis=-1)
    good = np.isfinite(e).all(axis=-1)
    R_est = R_s.copy()
    if good.any():
        R_est[good] = (
            Rotation.from_rotvec(np.radians(e[good].astype(np.float64))).as_matrix()
            @ vb.vctx.R_true
        )
    secs = time.perf_counter() - t0
    out, info = {}, {}
    for j in starts:
        out[j] = R_est[j]
        info[j] = dict(
            ok=float(good[j]), n_used=float(g["aux"][j, -1]), status=float(g["status"][j, -1]),
            solve_s=float(g["runtime"][j, -1]),
        )  # fmt: skip
    return out, info, secs


def run_adam(
    ctx: Any, vctx: Any, keys: np.ndarray, R0: np.ndarray, params: Any, box: float, seed: int
) -> Dict[str, Any]:
    """The hybrid RiemannianAdamOptimizer with the April settings on the case's coarse (8x) image
    stack. Counts hard evaluations (CountingCost) and differentiable evaluations (forward+backward
    steps) separately."""
    import bench_hp_sweep as hp
    from icenine.orientation_search import RiemannianAdamOptimizer

    lf = fs._W.local_fn
    diff = ctx.diff_fn
    groups = osw.group_pixels(keys, ctx.geo)
    diff.image_stack = osw.CoarseStack(groups, ctx.geo, osw.FACTOR, hp.FIXED_OW, osw.SCALE_INDEX)
    n_diff = [0]
    orig = diff.evaluate

    def counted(*a: Any, **k: Any) -> Any:
        n_diff[0] += 1
        return orig(*a, **k)

    diff.evaluate = counted  # type: ignore[assignment]
    cc = O.CountingCost(lf, vctx.vertices, vctx.voxel.phase, 10**6, ())
    opt = RiemannianAdamOptimizer(
        hard_cost_fn=cc, diff_cost_fn=diff, voxel_vertices=vctx.vertices,
        phase_index=vctx.voxel.phase, rng=np.random.default_rng(seed),
    )  # fmt: skip
    t0 = time.perf_counter()
    try:
        opt.optimize(
            np.asarray(R0, dtype=np.float64), box, n_steps=params.adam_n_steps, lr=params.adam_lr,
            scale=params.adam_scale, max_restarts=params.successive_restarts,
            max_convergence_cost=params.max_convergence_cost,
        )  # fmt: skip
    finally:
        secs = time.perf_counter() - t0
        diff.evaluate = orig  # type: ignore[assignment]
        diff.image_stack = None
    assert cc.best_R is not None
    return dict(R=cc.best_R, cost=cc.best_cost, hard=cc.n, diff=n_diff[0], secs=secs)


# ---------------------------------------------------------------------------
# One case, one task
# ---------------------------------------------------------------------------


def run_case(
    vb: Any,
    j: int,
    start: np.ndarray,
    seed: int,
    r_told_deg: float,
    keys: np.ndarray,
    gn_R: np.ndarray,
    raw_R: Any,
) -> Dict[str, np.ndarray]:
    ctx = fs._W.ctx
    rec, lf = fs._W.rec, fs._W.local_fn
    vctx = vb.vctx
    box, step0 = D.final_box(rec)
    R_true = np.asarray(vctx.R_true, dtype=np.float64)
    R0 = np.asarray(start, dtype=np.float64)
    methods = make_methods(box, step0, rec.params, math.radians(r_told_deg))
    nck, nm = len(CKPTS), len(ANYTIME)
    res_R = np.zeros((nm, nck, 3, 3))
    res_cost = np.zeros((nm, nck))
    res_used = np.zeros((nm, nck), dtype=np.int64)
    total = np.zeros(nm, dtype=np.int64)
    secs = np.zeros(nm)
    for mi, name in enumerate(ANYTIME):
        cc = O.CountingCost(lf, vctx.vertices, vctx.voxel.phase, MAX_BUDGET, CKPTS)
        n0 = lf.eval_count
        t0 = time.perf_counter()
        O.run_budgeted(lambda c, f=methods[name]: f(c, R0, seed + CFG_SEED[name]), cc)
        secs[mi] = time.perf_counter() - t0
        assert lf.eval_count - n0 == cc.n <= MAX_BUDGET, (name, lf.eval_count - n0, cc.n)
        total[mi] = cc.n
        for ki, c in enumerate(CKPTS):
            res_cost[mi, ki], res_R[mi, ki], res_used[mi, ki] = cc.at(c)
    out: Dict[str, np.ndarray] = dict(
        R_start=R0, R_true=R_true, r_told_deg=np.array(r_told_deg),
        cost_true=np.array(lf.evaluate(R_true, vctx.vertices, vctx.voxel.phase).cost),
        cost_start=np.array(lf.evaluate(R0, vctx.vertices, vctx.voxel.phase).cost),
        res_R=res_R, res_cost=res_cost, res_used=res_used, total_evals=total, secs=secs,
    )  # fmt: skip
    # the default finisher (FindOptimal MC + VarianceMinimizing), unmodified, for reference
    n0, t0 = lf.eval_count, time.perf_counter()
    R_f, c_f, _conv = nnrun.refine_fo(R0.astype(np.float32), vctx, seed)
    out.update(
        fin_R=R_f, fin_cost=np.array(c_f), fin_evals=np.array(lf.eval_count - n0),
        fin_secs=np.array(time.perf_counter() - t0),
    )  # fmt: skip
    # the Task 1 / T5 result, bit for bit (realistic variant of the T5 cases; NaN elsewhere)
    ref_ok = raw_R is not None and vb.vi == VI_REAL
    out["fin_vs_raw_maxabs"] = np.array(np.abs(R_f - raw_R).max() if ref_ok else np.nan)
    out["gn_cost"] = np.array(lf.evaluate(gn_R, vctx.vertices, vctx.voxel.phase).cost)
    ad = run_adam(ctx, vctx, keys, R0, rec.params, box, seed + 104729)
    out.update(
        adam_R=ad["R"], adam_cost=np.array(ad["cost"]), adam_hard=np.array(ad["hard"]),
        adam_diff=np.array(ad["diff"]), adam_secs=np.array(ad["secs"]),
    )  # fmt: skip
    return out


def task(item: Dict[str, Any]) -> Tuple[int, int, float]:
    assert fs._W is not None
    ctx = fs._W.ctx
    a = ctx.args
    vidx, vpos, ri, dirs = item["vidx"], item["vpos"], item["ri"], item["dirs"]
    t_start = time.time()
    r_told = float(ps.RADII_DEG[ri])
    res: Dict[Tuple[int, int], Dict[str, np.ndarray]] = {}
    gn: Dict[Tuple[int, int], Tuple[np.ndarray, Dict[str, float]]] = {}
    gn_secs = np.zeros((len(VARIANTS), len(dirs)))
    for vb in nnrun.variant_batches(
        ctx, vidx, vpos, ri, item["ref_nroi"], item["ref_fail"], VARIANTS
    ):
        starts = {
            j: (item["raw_start"][j] if item["raw_start"] is not None else vb.R_nom0[j])
            for j in dirs
        }
        est, info, secs = gn_from_starts(ctx, vidx, vb, starts)
        for jj, j in enumerate(dirs):
            gn[(vb.vi, jj)] = (est[j], info[j])
            gn_secs[vb.vi, jj] = secs / len(dirs)  # batch time per case
        for jj, j in enumerate(dirs):
            keys = nnrun.case_keys(ctx, vb, j)
            fs.attach_images(keys)
            seed = nnrun.b_seed(a, vpos, ri, j, VI_REAL)  # the same stream for both variants
            raw_R = item["raw_R"][j] if item["raw_R"] is not None else None
            res[(vb.vi, jj)] = run_case(vb, j, starts[j], seed, r_told, keys, est[j], raw_R)
    out: Dict[str, Any] = dict(
        dirs=np.array(dirs), vidx=np.array(vidx), ri=np.array(ri), kind=np.array(item["kind"]),
        methods=np.array(ANYTIME), ckpts=np.array(CKPTS), gn_secs=gn_secs,
    )  # fmt: skip
    for k in res[(0, 0)]:
        out[k] = np.stack(
            [np.stack([res[(vi, jj)][k] for jj in range(len(dirs))]) for vi in range(len(VARIANTS))]
        )  # (variant, case, ...)
    out["gn_R"] = np.stack(
        [np.stack([gn[(vi, jj)][0] for jj in range(len(dirs))]) for vi in range(len(VARIANTS))]
    )
    for key in ("ok", "n_used", "status", "solve_s"):
        out["gn_" + key] = np.array(
            [[gn[(vi, jj)][1][key] for jj in range(len(dirs))] for vi in range(len(VARIANTS))]
        )
    np.savez(item["path"], **out)
    return vidx, ri, time.time() - t_start


# ---------------------------------------------------------------------------
# Items and driver
# ---------------------------------------------------------------------------


def build_task_items(
    cache: Path, pilot: bool, n_voxels: int
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    t5, wargs = D.build_items(cache / "t5_tmp", only_first=2 if pilot else 0)
    items: List[Dict[str, Any]] = []
    for vidx, vpos, ri, dirs, pipe, _p, nroi, fail, Rf, St in t5:
        items.append(
            dict(kind=pipe, vidx=vidx, vpos=vpos, ri=ri, dirs=list(dirs), ref_nroi=nroi,
                 ref_fail=fail, raw_start=St, raw_R=Rf,
                 path=str(cache / f"T5{pipe}_v{vidx}_r{ri}.npz"))
        )  # fmt: skip
    _w, sweep = nnrun.worker_args(["realistic_s0"])
    voxels = [int(v) for v in sweep["voxel_indices"]]
    if pilot:
        voxels = voxels[:2]
    elif n_voxels:
        voxels = voxels[:n_voxels]
    for vpos, v in enumerate(voxels):
        for ri in SW_RADII_IDX[:2] if pilot else SW_RADII_IDX:
            valid = np.nonzero((sweep["fail_pass1"][vpos, ri] == 0).all(axis=-1))[0]
            dirs = [int(j) for j in valid[: 2 if pilot else SW_DIRS]]
            if dirs:
                items.append(
                    dict(kind="SW", vidx=v, vpos=vpos, ri=ri, dirs=dirs,
                         ref_nroi=sweep["n_roi"][vpos, ri], ref_fail=sweep["fail_pass1"][vpos, ri],
                         raw_start=None, raw_R=None, path=str(cache / f"SW_v{v}_r{ri}.npz"))
                )  # fmt: skip
    return items, wargs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pilot", "run"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--n-voxels", type=int, default=0, help="first N sweep voxels (0 = all 50)")
    ap.add_argument("--cache", default="")
    args = ap.parse_args()
    import multiprocessing as mp

    cache = CACHE_DIR / (args.cache or args.cmd)
    cache.mkdir(parents=True, exist_ok=True)
    (cache / "t5_tmp").mkdir(exist_ok=True)
    items, wargs = build_task_items(cache, args.cmd == "pilot", args.n_voxels)
    # longest tasks first
    items.sort(key=lambda it: -len(it["dirs"]))
    todo = [it for it in items if not Path(it["path"]).exists()]
    n_cases = sum(len(it["dirs"]) for it in items) * len(VARIANTS)
    print(
        f"{len(items)} tasks ({len(items) - len(todo)} cached), {n_cases} cases, "
        f"{args.workers} workers",
        flush=True,
    )
    t0 = time.time()
    with mp.get_context("spawn").Pool(
        min(args.workers, max(len(todo), 1)), initializer=fs.init_worker, initargs=(wargs,)
    ) as pool:
        for k, (vidx, ri, secs) in enumerate(pool.imap_unordered(task, todo), 1):
            print(
                f"  [{k}/{len(todo)}] voxel {vidx} r#{ri} {secs:.0f}s ({time.time()-t0:.0f}s)",
                flush=True,
            )
    print(f"finished; wall {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
