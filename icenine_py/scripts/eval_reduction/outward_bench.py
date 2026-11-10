"""Task F: searching outward from a start orientation, per method, start offset and direction.

Cases. A voxel (the E0 set: per-voxel clean images, or the six interior voxels of the clean
full-sample images), the true orientation R_true, and starts at ``OFFSETS`` degrees from it in
``--dirs`` random directions (the same directions at every offset). Every method gets the start and
the images and nothing else, except F1 which is told the offset (its box is the offset). The strict
cost is the library's local cost function (lenient rule, as the E0 and Phase D runs use).

Methods (``icenine.outward_search`` for f4, f5):
  f1     FindOptimal (MC + VarianceMinimizing) with the box widened to max(default box, offset)
  f2_05, f2_1, f2_2   CMA-ES, sigma0 0.5, 1, 2 deg, stops on its own tolerance or 4000 evaluations
  f3     centroid Gauss-Newton with Huber loss, three re-centring passes, windows of the images
         at the start (E0 only; the HG path without the network)
  f4     Gauss-Newton on the observed-ray residual (Task E arm E4), trust region
  f5     expanding shells (screen, then E4, then strict), stops when a strict cost passes
  f54    f5 followed by f4 from the point f5 accepted (or from its best point if none passed)

Counted per run: strict evaluations, screen evaluations (batched pixel-radius screen), smooth
evaluations (E4 values and E4 value-and-Jacobian evaluations), wall time. Diagnostics from the
truth: the error of each method's current best against the cumulative evaluations (a curve on the
checkpoints ``CKPT`` and the evaluations until the error first falls below 0.05 and 1 deg).

  uv run python scripts/eval_reduction/outward_bench.py pilot
  nohup uv run python scripts/eval_reduction/outward_bench.py run --kind e0 --workers 5 &
  uv run python scripts/eval_reduction/outward_bench.py timing   # single worker, preflight-gated
  uv run python scripts/eval_reduction/outward_bench.py timing-f5 --kind e0 --dirs 5  # F5 batches

Outputs: benchmarks/eval_reduction/F/cache/<set>/v<voxel>_o<offset index>.npz (gitignored).
"""

import argparse
import contextlib
import io
import json
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
for sub in (
    HERE,
    ICE / "scripts",
    ICE / "scripts" / "common",
    ICE / "scripts" / "finisher_bench",
    ICE / "scripts" / "findoptimal_robustness",
    ICE / "benchmarks",
):
    sys.path.insert(0, str(sub))
import cost_edges as CE  # noqa: E402
import cost_edges_run as RUN  # noqa: E402

OUT = ICE / "benchmarks" / "eval_reduction" / "F"
OFFSETS = [0.1, 0.5, 1.0, 2.0, 3.0, 5.0]
CKPT = [3, 10, 30, 100, 300, 1000, 3000, 10000, 30000]
THR = (0.05, 0.1, 1.0)
METHODS = ["f1", "f2_05", "f2_1", "f2_2", "f3", "f4", "f5", "f54"]
EXTRA = ["f4w"]  # run separately (--methods f4w --tag <tag>_f4w), merged by outward_summary
F4W_FRAMES = 5  # frame window of f4w, chosen on the pilot voxels (RUN_LOG 3)
SIGMA0 = {"f2_05": 0.5, "f2_1": 1.0, "f2_2": 2.0}
CMA_BUDGET = 4000
CMA_TOLX_DEG = 0.005  # CMA-ES stops when its step-size x scale falls below this
ACCEPT = 0.8  # strict cost that passes the acceptance test of f4 and f5 (see RUN_LOG, section 3)
SPACING_DEG = 0.3  # f5 shell spacing (see RUN_LOG, section 3)
R_MAX_DEG = 6.0
TAU_SCHEDULE = (0.25,)  # a pilot on four other voxels beat the annealed schedules (RUN_LOG 3)
NEAR_UM = 120.0  # near-axis: r_perp below this (the crossover of the HP sweep findings)

_E: Dict[str, Any] = {}  # worker state


# --------------------------------------------------------------------------------------------
# voxel sets
# --------------------------------------------------------------------------------------------
def e0_sets(n_near: int, n_far: int, n_pilot: int = 4) -> Dict[str, List[int]]:
    """Near-axis (r_perp < 120 um) and far-from-axis (>= 250 um) voxels of the E0 200 that have
    clean images, spread evenly; pilot voxels are taken from the rest and not used afterwards."""
    import common as C

    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]
    rp = np.asarray(info["voxel_r_perp_um"])[: C.N_VOXELS]
    ok = np.array([(C.CACHE_DIR / "images" / f"v{v}_clean.npz").exists() for v in vox])
    near = [v for v, r, o in zip(vox, rp, ok) if o and r < NEAR_UM]
    far = [v for v, r, o in zip(vox, rp, ok) if o and r >= 250.0]
    take = lambda L, n: (
        [L[int(i)] for i in np.linspace(0, len(L) - 1, n)] if n < len(L) else L
    )  # noqa: E731
    pilot_near = near[::-1][: n_pilot // 2]
    pilot_far = far[::-1][: n_pilot - len(pilot_near)]
    near = [v for v in near if v not in pilot_near]
    far = [v for v in far if v not in pilot_far]
    # the pilot voxels were removed from the pool; the sets are spread over what is left
    return dict(
        near=take(near, n_near), far=take(far, n_far), pilot=pilot_near + pilot_far,
        n_near_all=near, n_far_all=far,
    )  # fmt: skip


# --------------------------------------------------------------------------------------------
# environment of one voxel
# --------------------------------------------------------------------------------------------
class Tracer:
    """Counts the calls of the strict cost function ``lf.evaluate`` and records the error to the
    truth of its best-cost orientation each time the best improves (a trace of (n, error deg))."""

    def __init__(self, lf: Any, R_true: np.ndarray, verts: Any = None, phase: int = 0) -> None:
        self.lf, self.R_true, self.verts, self.phase = lf, R_true, verts, phase
        self.n = 0
        self.best_c = math.inf
        self.best_R: Optional[np.ndarray] = None
        self.trace: List[Tuple[int, float]] = []
        self._orig = lf.evaluate

    def __enter__(self) -> "Tracer":
        from icenine.outward_search import angle_deg

        orig = self._orig

        def evaluate(orientation: np.ndarray, voxel_vertices: Any = None, phase_index: int = 0):
            info = orig(orientation, voxel_vertices, phase_index)
            self.n += 1
            c = float(info.cost)
            if c < self.best_c:
                self.best_c = c
                self.best_R = np.asarray(orientation, dtype=np.float64).copy()
                self.trace.append((self.n, angle_deg(self.best_R, self.R_true)))
            return info

        self.lf.evaluate = evaluate  # type: ignore[method-assign]
        return self

    def __exit__(self, *a: Any) -> None:
        self.lf.__dict__.pop("evaluate", None)

    def strict(self, R: np.ndarray) -> float:
        return float(self.lf.evaluate(np.asarray(R, dtype=np.float64), self.verts, self.phase).cost)


def init_worker(kind: str) -> None:
    import torch

    torch.set_num_threads(1)
    _E["kind"] = kind
    (RUN.init_full if kind == "full" else RUN.init_e0)()
    if kind == "full":
        _E["rec"] = RUN._ENV["ctx"].rec
        _E["lf"] = _E["rec"]._make_local_cost_fn()


def voxel_env(voxel: int) -> Dict[str, Any]:
    """Everything one voxel needs; E0 attaches the voxel's images (cached while it stays the
    current voxel). The SmoothCost, the screen cost and the observed rays are built lazily."""
    if _E.get("voxel") == voxel:
        return _E["env"]
    base, Rt, verts, phase = RUN.case_data(voxel)
    if _E["kind"] == "e0":
        import common as C
        from optimizer_sweep import voxel_context

        W = C.get_worker()
        rec, lf = W.rec, W.local_fn
        vctx = voxel_context(W.ctx, voxel)
        geo = dict(ctx=W.ctx, vctx=vctx)
    else:
        rec, lf, geo = _E["rec"], _E["lf"], {}
    env = dict(base=base, R_true=Rt, verts=verts, phase=phase, rec=rec, lf=lf, voxel=voxel, **geo)
    if _E["kind"] == "full" and _E.get("shared") is not None:
        env.update(_E["shared"])
    _E.update(voxel=voxel, env=env)
    return env


def lazy(env: Dict[str, Any], what: str) -> Any:
    if what in env:
        return env[what]
    from icenine.smooth_cost import SmoothCost

    b = env["base"]
    if what == "sc":
        env[what] = SmoothCost(
            b.simulator, b.detector_list, b.range_map, b.exp_data, b.sample, b.structure_list,
            eta_limit=b.eta_limit,
        )  # fmt: skip
        env["sc"].observed_rays()
    elif what == "screen_cf":
        env[what] = CE.make_cf(b, "screen", "old")
    if _E["kind"] == "full":
        _E.setdefault("shared", {})[what] = env[what]
    return env[what]


# --------------------------------------------------------------------------------------------
# methods
# --------------------------------------------------------------------------------------------
def curve_of(trace: List[Tuple[int, float]], err0: float) -> np.ndarray:
    """Error of the method's current best at each checkpoint (the start's error before the first
    trace point); ``trace`` is [(cumulative evaluations, error)] with increasing n."""
    out = np.empty(len(CKPT))
    for i, c in enumerate(CKPT):
        e = err0
        for n, er in trace:
            if n <= c:
                e = er
            else:
                break
        out[i] = e
    return out


CAPTURE_TOL = 1e-3  # deg. R_true is float32: a start exactly 1 deg out measures 0.9999 deg


def evals_to(trace: List[Tuple[int, float]], thr: float, err0: float) -> float:
    if err0 < thr - CAPTURE_TOL:
        return 0.0
    for n, er in trace:
        if er < thr - CAPTURE_TOL:
            return float(n)
    return math.nan


def pack(R: np.ndarray, R_true: np.ndarray, n_strict: int, n_screen: int, n_smooth: int,
         n_jac: int, secs: float, trace: List[Tuple[int, float]], err0: float, strict_cost: float,
         flag: float) -> Dict[str, Any]:  # fmt: skip
    from icenine.outward_search import angle_deg

    return dict(
        err=angle_deg(R, R_true), cost=strict_cost, n_strict=n_strict, n_screen=n_screen,
        n_smooth=n_smooth, n_jac=n_jac, secs=secs, curve=curve_of(trace, err0),
        cap005=evals_to(trace, THR[0], err0), cap01=evals_to(trace, THR[1], err0),
        cap1=evals_to(trace, THR[2], err0), flag=flag,
    )  # fmt: skip


def run_f1(env: Dict[str, Any], R0: np.ndarray, offset: float, seed: int) -> Dict[str, Any]:
    """FindOptimal as the reconstructor's final stage, box = max(default box, offset)."""
    from icenine.orientation_search import MCOptimizer, SearchCandidate
    from icenine.outward_search import angle_deg

    rec, lf, verts, phase = env["rec"], env["lf"], env["verts"], env["phase"]
    p = rec.params
    box0 = max(
        (p.local_grid_radius / 1.5 ** (p.max_local_resolution + 1)) / 3.0, math.radians(0.2)
    ) / (2**p.min_local_resolution)
    box = max(box0, math.radians(offset))
    diameter = 3.0 * box * (2**p.min_local_resolution)
    err0 = angle_deg(R0, env["R_true"])
    with Tracer(lf, env["R_true"]) as tr:
        t0 = time.perf_counter()
        mc = MCOptimizer(
            cost_fn=lf, voxel_vertices=verts, phase_index=phase, rng=np.random.default_rng(seed)
        )
        with contextlib.redirect_stdout(io.StringIO()):
            res = rec.refine_from_candidates(
                [SearchCandidate(orientation=R0.astype(np.float32), cost=1.0)],
                verts, phase, diameter=diameter, local_cost_fn=lf, mc_optimizer=mc,
            )  # fmt: skip
        secs = time.perf_counter() - t0
    R = np.asarray(res.orientation, dtype=np.float64)
    return pack(R, env["R_true"], tr.n, 0, 0, 0, secs, tr.trace, err0, float(res.cost), 1.0)


def cma_run(f: Any, R0: np.ndarray, seed: int, sigma0_deg: float) -> None:
    """finisher_bench.optimizers.cma_local with the x-tolerance as a parameter (``CMA_TOLX_DEG``):
    local CMA-ES on a rotation vector in degrees about R0, the package's flat-fitness tests off
    (the cost is quantised); stops on the x-tolerance or when the budget ends it."""
    import cma
    import optimizers as FO

    f(R0)
    opts = dict(seed=int(seed) + 1, verbose=-9, tolfun=0.0, tolfunhist=0.0,
                tolflatfitness=10**9, tolstagnation=10**9, tolx=CMA_TOLX_DEG)  # fmt: skip
    es = cma.CMAEvolutionStrategy(np.zeros(3), float(sigma0_deg), opts)
    while not es.stop():
        X = es.ask()
        es.tell(X, [f(FO.rot_about(R0, x)) for x in X])


def run_f2(env: Dict[str, Any], R0: np.ndarray, sigma0: float, seed: int) -> Dict[str, Any]:
    import optimizers as FO
    from icenine.outward_search import angle_deg

    lf, verts, phase = env["lf"], env["verts"], env["phase"]
    err0 = angle_deg(R0, env["R_true"])
    with Tracer(lf, env["R_true"]) as tr:
        cc = FO.CountingCost(lf, verts, phase, CMA_BUDGET, ())
        t0 = time.perf_counter()
        FO.run_budgeted(lambda c: cma_run(c, R0, seed, sigma0), cc)
        secs = time.perf_counter() - t0
    return pack(cc.best_R, env["R_true"], tr.n, 0, 0, 0, secs, tr.trace, err0, cc.best_cost, 1.0)


def run_f3_batch(env: Dict[str, Any], starts: np.ndarray, offset: float) -> List[Dict[str, Any]]:
    """Centroid Huber GN, three passes, on the windows of the images at each start (one batch)."""
    import optimizer_sweep as osw
    import perturbation_sweep as ps
    from scipy.spatial.transform import Rotation
    from icenine.outward_search import angle_deg

    ctx, vctx = env["ctx"], env["vctx"]
    a = ctx.args
    D = len(starts)
    R_true = env["R_true"]
    t0 = time.perf_counter()
    delta_s = ps.relative_offset_deg(R_true, starts)
    preps = [ps.prepare_nominal(ctx, env["voxel"], starts[j])[0] for j in range(D)]
    batch = ps.render_batch(preps, delta_s, None, vctx.sources, "clean", 0, a)
    ok = np.array([p is not None for p in preps]) & (batch["n_present"].numpy() >= a.min_present)
    g = osw.run_gn(
        ctx, env["voxel"], preps, batch, ok, delta_s, starts, R_true, None, vctx.sources,
        "clean", 0, 1.0,
    )  # fmt: skip
    e = np.stack([g["err_x"][:, -1], g["err_y"][:, -1], g["err_z"][:, -1]], axis=-1)
    good = np.isfinite(e).all(axis=-1)
    R_est = starts.copy()
    if good.any():
        R_est[good] = (
            Rotation.from_rotvec(np.radians(e[good].astype(np.float64))).as_matrix() @ R_true
        )
    secs = (time.perf_counter() - t0) / D
    out = []
    lf = env["lf"]
    for j in range(D):
        err0 = angle_deg(starts[j], R_true)
        c = float(lf.evaluate(R_est[j], env["verts"], env["phase"]).cost)
        npass = int(np.isfinite(g["err_angle"][j]).sum())
        tr = [(max(npass, 1), angle_deg(R_est[j], R_true))] if good[j] else []
        out.append(pack(R_est[j], R_true, 0, 0, 0, 0, secs, tr, err0, c, float(good[j])))
        out[-1]["npass"] = float(npass)
    return out


def make_gn(
    env: Dict[str, Any], tr: Tracer, taus: Tuple[float, ...] = TAU_SCHEDULE, window: int = 1
) -> Any:
    from icenine.outward_search import RayGaussNewton

    return RayGaussNewton(
        lazy(env, "sc"), env["verts"], tr.strict, env["phase"], tau_schedule_deg=taus,
        accept_cost=ACCEPT, frame_window=window,
    )  # fmt: skip


def run_f4(env: Dict[str, Any], R0: np.ndarray, window: int = 1) -> Dict[str, Any]:
    from icenine.outward_search import angle_deg

    err0 = angle_deg(R0, env["R_true"])
    trace: List[Tuple[int, float]] = []
    with Tracer(env["lf"], env["R_true"], env["verts"], env["phase"]) as tr:
        gn = make_gn(env, tr, window=window)
        t0 = time.perf_counter()
        res = gn.solve(
            R0, callback=lambda R, ns, nt: trace.append((ns + nt, angle_deg(R, env["R_true"])))
        )
        secs = time.perf_counter() - t0
    trace.append((res.n_evals, angle_deg(res.R, env["R_true"])))
    return pack(res.R, env["R_true"], res.n_strict, 0, res.n_smooth, res.n_jac, secs, trace,
                err0, res.strict_cost, float(res.accepted))  # fmt: skip


def run_f5(
    env: Dict[str, Any], R0: np.ndarray, spacing: float = SPACING_DEG
) -> Tuple[Dict[str, Any], Any]:
    from icenine.outward_search import ShellSearch, angle_deg, screen_fn_from, smooth_fn_from

    err0 = angle_deg(R0, env["R_true"])
    with Tracer(env["lf"], env["R_true"], env["verts"], env["phase"]) as tr:
        scr = screen_fn_from(lazy(env, "screen_cf"), env["verts"], env["phase"])
        smo = smooth_fn_from(lazy(env, "sc"), env["verts"], env["phase"], "E4,tau=0.25")
        ss = ShellSearch(scr, tr.strict, smo, spacing_deg=spacing, r_max_deg=R_MAX_DEG,
                         accept_cost=ACCEPT)  # fmt: skip
        t0 = time.perf_counter()
        res = ss.search(R0)
        secs = time.perf_counter() - t0
    trace = [(n, angle_deg(Rb, env["R_true"])) for _, _, n, Rb in res.trace]
    out = pack(res.R, env["R_true"], res.n_strict, res.n_screen, res.n_smooth, 0, secs, trace,
               err0, res.strict_cost, float(res.found))  # fmt: skip
    out["radius"] = float(res.radius_deg)
    return out, res


def run_f54(
    env: Dict[str, Any], R0: np.ndarray, f5: Dict[str, Any], f5res: Any, window: int = 1
) -> Dict[str, Any]:
    """f5, then f4 from f5's result; counts and time are the sums, the trace is f5's followed by
    f4's shifted by f5's evaluations. ``window`` is f4's frame window (Fcompare uses 5)."""
    from icenine.outward_search import angle_deg

    err0 = angle_deg(R0, env["R_true"])
    base_n = f5res.n_evals
    trace = [(n, angle_deg(Rb, env["R_true"])) for _, _, n, Rb in f5res.trace]
    with Tracer(env["lf"], env["R_true"], env["verts"], env["phase"]) as tr:
        gn = make_gn(env, tr, window=window)
        t0 = time.perf_counter()
        res = gn.solve(
            f5res.R,
            callback=lambda R, ns, nt: trace.append(
                (base_n + ns + nt, angle_deg(R, env["R_true"]))
            ),
        )
        secs = time.perf_counter() - t0
    trace.append((base_n + res.n_evals, angle_deg(res.R, env["R_true"])))
    return pack(res.R, env["R_true"], f5res.n_strict + res.n_strict, f5res.n_screen,
                f5res.n_smooth + res.n_smooth, res.n_jac, f5["secs"] + secs, trace, err0,
                res.strict_cost, float(res.accepted))  # fmt: skip


# --------------------------------------------------------------------------------------------
# one (voxel, offset)
# --------------------------------------------------------------------------------------------
def starts_for(env: Dict[str, Any], offset: float, nd: int) -> np.ndarray:
    dirs = CE.unit_directions(nd, 7000 + env["voxel"])
    return np.stack([CE.rotate_about(env["R_true"], np.array([offset]), d)[0] for d in dirs])


def run_group(env: Dict[str, Any], oi: int, nd: int, methods: List[str]) -> Dict[str, np.ndarray]:
    from icenine.outward_search import angle_deg

    offset = OFFSETS[oi]
    starts = starts_for(env, offset, nd)
    lf, R_true = env["lf"], env["R_true"]
    rows: Dict[str, List[Dict[str, Any]]] = {m: [] for m in methods}
    if "f3" in methods:
        if _E["kind"] == "e0":
            rows["f3"] = run_f3_batch(env, starts, offset)
        else:
            methods = [m for m in methods if m != "f3"]
            rows.pop("f3")
    for d in range(nd):
        R0 = starts[d]
        seed = int(np.random.default_rng([env["voxel"], oi, d]).integers(2**31))
        if "f1" in methods:
            rows["f1"].append(run_f1(env, R0, offset, seed))
        for m, sg in SIGMA0.items():
            if m in methods:
                rows[m].append(run_f2(env, R0, sg, seed))
        if "f4" in methods:
            rows["f4"].append(run_f4(env, R0))
        if "f4w" in methods:
            rows["f4w"].append(run_f4(env, R0, F4W_FRAMES))
        if "f5" in methods or "f54" in methods:
            r5, res5 = run_f5(env, R0)
            if "f5" in methods:
                rows["f5"].append(r5)
            if "f54" in methods:
                rows["f54"].append(run_f54(env, R0, r5, res5))
    out: Dict[str, np.ndarray] = dict(
        err_start=np.array([angle_deg(s, R_true) for s in starts]),
        cost_start=np.array(
            [float(lf.evaluate(s, env["verts"], env["phase"]).cost) for s in starts]
        ),
        cost_true=np.array(float(lf.evaluate(R_true, env["verts"], env["phase"]).cost)),
        offset=np.array(offset),
    )
    for m, rs in rows.items():
        if not rs:
            continue
        for k in rs[0]:
            out[f"{m}|{k}"] = np.array([r[k] for r in rs])
    return out


def task(item: Tuple[str, int, str, int, List[str]]) -> str:
    kind, voxel, tag, nd, methods = item
    t0 = time.time()
    done = 0
    for oi in range(len(OFFSETS)):
        path = OUT / "cache" / tag / f"v{voxel}_o{oi}.npz"
        if path.exists():
            continue
        if _E.get("only_offsets") and oi not in _E["only_offsets"]:
            continue
        env = voxel_env(voxel)
        res = run_group(env, oi, nd, methods)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, **res)
        done += 1
    return f"{kind} v{voxel}: {done} offsets written, {time.time() - t0:.0f}s"


def _init(kind: str, only: Optional[List[int]]) -> None:
    init_worker(kind)
    _E["only_offsets"] = only


def run_pool(kind: str, voxels: List[int], tag: str, nd: int, workers: int,
             methods: List[str], only: Optional[List[int]] = None) -> None:  # fmt: skip
    (OUT / "cache" / "logs").mkdir(parents=True, exist_ok=True)
    items = [(kind, v, tag, nd, methods) for v in voxels]
    print(f"{tag}: {len(items)} voxels, {workers} workers, {nd} dirs", flush=True)
    t0 = time.time()
    if workers <= 1:
        _init(kind, only)
        for i, it in enumerate(items):
            print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {task(it)}", flush=True)
    else:
        with mp.get_context("spawn").Pool(workers, initializer=_init, initargs=(kind, only)) as p:
            for i, msg in enumerate(p.imap_unordered(task, items)):
                print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {msg}", flush=True)
    (OUT / "cache" / "logs" / f"{tag}.done").write_text(f"{time.time() - t0:.0f}s\n")


def timing(kind: str, voxels: List[int], dirs: int) -> None:
    """Single-worker timing (preflight-gated): per-evaluation costs of each kind on the first
    voxel, then a small serial run of every method. Writes F/timing/<kind>.json (+ preflight)."""
    import preflight
    from icenine.outward_search import RayGaussNewton, screen_fn_from

    pf = preflight.preflight()
    (OUT / "timing").mkdir(parents=True, exist_ok=True)
    (OUT / "timing" / f"preflight_{kind}.json").write_text(json.dumps(pf, indent=1) + "\n")
    preflight.require_quiet(info=pf)
    _init(kind, [2, 3, 5])
    env = voxel_env(voxels[0])
    R_true, verts, phase, lf = env["R_true"], env["verts"], env["phase"], env["lf"]
    starts = starts_for(env, 1.0, 64)
    ms: Dict[str, float] = {}
    t0 = time.perf_counter()
    for R in starts[: 500 if len(starts) >= 500 else len(starts)]:
        lf.evaluate(R, verts, phase)
    ms["strict"] = 1e3 * (time.perf_counter() - t0) / min(500, len(starts))
    scr = screen_fn_from(lazy(env, "screen_cf"), verts, phase)
    sc = lazy(env, "sc")
    scr(starts[:8])
    t0 = time.perf_counter()
    for _ in range(10):
        scr(starts)
    ms["screen"] = 1e3 * (time.perf_counter() - t0) / (10 * len(starts))
    with Tracer(lf, R_true, verts, phase) as tr:
        gn = make_gn(env, tr)
        gn.e4_value(starts[0], 0.004)
        t0 = time.perf_counter()
        for R in starts[:30]:
            gn.e4_value(R, 0.004)
        ms["e4_value"] = 1e3 * (time.perf_counter() - t0) / 30
        t0 = time.perf_counter()
        for R in starts[:30]:
            gn.linearize(R, 0.004)
        ms["e4_value_and_jacobian"] = 1e3 * (time.perf_counter() - t0) / 30
    t0 = time.time()
    run_pool_serial = [(kind, v, "timing_" + kind, dirs, METHODS + EXTRA) for v in voxels]
    for it in run_pool_serial:
        task(it)
    ms["serial_run_s"] = time.time() - t0
    (OUT / "timing" / f"{kind}.json").write_text(
        json.dumps(
            dict(label="single-worker", ms_per_evaluation=ms, voxels=voxels, dirs=dirs), indent=1
        )
        + "\n"
    )
    print(json.dumps(ms, indent=1))


def timing_f5(kind: str, voxels: List[int], dirs: int) -> None:
    """Single-worker timing (preflight-gated) of F5's evaluation kinds at the batch sizes F5 uses
    them: ``ShellSearch`` is run from the timing run's starts (offsets 1, 2, 5 deg, ``dirs``
    directions, every voxel of ``voxels``) with timers around its callables. The screen is called
    in chunks of up to 256 rotations (a whole shell), E4 on the best 32 of a shell through
    ``SmoothCost.evaluate_many``, the strict cost one rotation at a time. Also times F4's single
    E4 value and value-and-Jacobian on every voxel (the per-evaluation weights of ``timing`` come
    from the first voxel only). Writes F/timing/<kind>_f5.json (+ preflight_<kind>_f5.json)."""
    import preflight
    from icenine.outward_search import ShellSearch, screen_fn_from, smooth_fn_from

    pf = preflight.preflight()
    (OUT / "timing").mkdir(parents=True, exist_ok=True)
    (OUT / "timing" / f"preflight_{kind}_f5.json").write_text(json.dumps(pf, indent=1) + "\n")
    preflight.require_quiet(info=pf)
    _init(kind, [2, 3, 5])
    per_voxel: Dict[str, Any] = {}
    for v in voxels:
        env = voxel_env(v)
        verts, phase, lf = env["verts"], env["phase"], env["lf"]
        acc = {k: [0.0, 0] for k in ("screen", "smooth", "strict")}  # ms, orientations
        calls = {k: 0 for k in acc}

        def timed(fn: Any, k: str, many: bool) -> Any:
            def g(x: np.ndarray) -> Any:
                t0 = time.perf_counter()
                out = fn(x)
                acc[k][0] += 1e3 * (time.perf_counter() - t0)
                acc[k][1] += len(x) if many else 1
                calls[k] += 1
                return out

            return g

        scr = screen_fn_from(lazy(env, "screen_cf"), verts, phase)
        smo = smooth_fn_from(lazy(env, "sc"), verts, phase, "E4,tau=0.25")

        def strict(R: np.ndarray) -> float:
            return float(lf.evaluate(np.asarray(R, dtype=np.float64), verts, phase).cost)

        warm = starts_for(env, 1.0, 8)
        scr(warm), smo(warm), strict(warm[0])
        ss = ShellSearch(timed(scr, "screen", True), timed(strict, "strict", False),
                         timed(smo, "smooth", True), spacing_deg=SPACING_DEG,
                         r_max_deg=R_MAX_DEG, accept_cost=ACCEPT)  # fmt: skip
        search_ms = 0.0
        for oi in (2, 3, 5):
            for R0 in starts_for(env, OFFSETS[oi], dirs):
                t0 = time.perf_counter()
                ss.search(R0)
                search_ms += 1e3 * (time.perf_counter() - t0)
        rec = {
            f"{k}_ms_per_orientation": (acc[k][0] / acc[k][1]) if acc[k][1] else math.nan
            for k in acc
        }
        rec.update({f"{k}_orientations": acc[k][1] for k in acc})
        rec.update({f"{k}_calls": calls[k] for k in calls})
        rec["search_ms_total"] = search_ms
        rec["callables_ms_total"] = sum(a[0] for a in acc.values())
        with Tracer(lf, env["R_true"], verts, phase) as tr:
            gn = make_gn(env, tr)
            st = starts_for(env, 1.0, 30)
            gn.e4_value(st[0], 0.004)
            t0 = time.perf_counter()
            for R in st:
                gn.e4_value(R, 0.004)
            rec["f4_e4_value_ms"] = 1e3 * (time.perf_counter() - t0) / len(st)
            t0 = time.perf_counter()
            for R in st:
                gn.linearize(R, 0.004)
            rec["f4_e4_value_and_jacobian_ms"] = 1e3 * (time.perf_counter() - t0) / len(st)
        per_voxel[str(v)] = rec
        print(v, json.dumps(rec, indent=1), flush=True)
    (OUT / "timing" / f"{kind}_f5.json").write_text(
        json.dumps(
            dict(
                label="single-worker",
                voxels=voxels,
                dirs=dirs,
                offsets_deg=[1.0, 2.0, 5.0],
                per_voxel=per_voxel,
            ),  # fmt: skip
            indent=1,
        )
        + "\n"
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pilot", "run", "sets", "timing", "timing-f5"])
    ap.add_argument("--kind", choices=["e0", "full"], default="e0")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--dirs", type=int, default=20)
    ap.add_argument("--n-near", type=int, default=8)
    ap.add_argument("--n-far", type=int, default=8)
    ap.add_argument("--tag", default="")
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--voxels", default="")
    ap.add_argument("--offsets", default="", help="indices into OFFSETS, comma separated")
    a = ap.parse_args()
    methods = a.methods.split(",")
    only = [int(x) for x in a.offsets.split(",")] if a.offsets else None
    if a.cmd == "sets":
        s = e0_sets(a.n_near, a.n_far)
        print(json.dumps({k: v for k, v in s.items() if not k.endswith("_all")}))
        print("pool sizes near/far:", len(s["n_near_all"]), len(s["n_far_all"]))
        return
    if a.cmd == "timing":
        s = e0_sets(a.n_near, a.n_far)
        timing(a.kind, [s["near"][0], s["far"][0]] if a.kind == "e0" else [15901], a.dirs)
        return
    if a.cmd == "timing-f5":
        s = e0_sets(a.n_near, a.n_far)
        timing_f5(a.kind, [s["near"][0], s["far"][0]] if a.kind == "e0" else [15901], a.dirs)
        return
    if a.kind == "full":
        from parity_common import SIX_VOXELS

        vox = list(SIX_VOXELS)
    elif a.voxels:
        vox = [int(x) for x in a.voxels.split(",")]
    else:
        s = e0_sets(a.n_near, a.n_far)
        vox = s["pilot"] if a.cmd == "pilot" else s["near"] + s["far"]
    tag = a.tag or (f"{a.kind}_pilot" if a.cmd == "pilot" else a.kind)
    run_pool(a.kind, vox, tag, a.dirs, a.workers, methods, only)


if __name__ == "__main__":
    main()
