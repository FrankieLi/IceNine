#!/usr/bin/env python3
"""
T5 finisher diagnosis: why does refine_from_candidates stop above the truth's cost?

Cases: realistic ("all") single-voxel perturbation-sweep cases of the Task 1 pipelines H3 (net x3
start) and H0 (perturbed nominal start), radii <= 3 deg, stratified by radius. The finisher is run
again with the same seed through the UNMODIFIED refine_from_candidates, with LoggedMC (a
MCOptimizer subclass whose optimize / variance_minimizing_optimize are copies that also record why
and when they stopped; same random stream) as its mc_optimizer. The re-run result is checked
against the Task 1 raw R_final.

Per case (see MIGRATION_HISTORY.md "T5 results"):
  (a) cost along the straight geodesic result -> truth
  (b) continuations from the result: FindOptimal MC (long / small step / small box) and
      VarianceMinimizing (long / small box) called as components, plus a default finisher re-run
      from the result and from the truth
  (c) cost granularity: cost at truth and at result +- tiny rotations
  (d) stopping rule and iteration of the finisher run

Usage (from icenine_py/):
  uv run python scripts/finisher_diagnosis/diagnose.py pilot
  uv run python scripts/finisher_diagnosis/diagnose.py run --workers 10
"""

import argparse
import contextlib
import io
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "nn_hybrid"))
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(ICENINE_PY / "benchmarks"))

import findoptimal_sweep as fs  # noqa: E402
import run as nnrun  # noqa: E402  (scripts/nn_hybrid/run.py)
from stats import reorder  # noqa: E402

from icenine.orientation_search import (  # noqa: E402
    MCOptimizer,
    SearchCandidate,
    _quat_multiply,
    matrix_to_quaternion,
    quaternion_to_matrix,
)

CACHE_DIR = HERE / "cache"
OUT_DIR = ICENINE_PY / "benchmarks" / "finisher_diagnosis"
RAW_DIR = ICENINE_PY / "benchmarks" / "nn_hybrid"
VI = 1  # the realistic ("all") variant
N_RADII = 9  # radii index 0..8 = 0.05 .. 3 deg
N_CASES = {"H3": 200, "H0": 100}
VOX_PER_RADIUS = {"H3": 4, "H0": 2}
SELECT_SEED = 20261006
N_GEO = 41  # uniform points along the geodesic (plus near-start points)
GEO_NEAR = np.geomspace(1e-3, 0.1, 7)
GRAN_ANGLES_DEG = np.array([0.001, 0.003, 0.01, 0.02, 0.05, 0.1, 0.2])
GRAN_DIRS = 12
LONG_FACTOR = 25  # long continuations: 25x the config max_mc_steps (200)
SMALL_STEP_DIV = 10.0
SMALL_BOX_DIV = 4.0

STOP_NAMES = ["step_budget", "restarts_exhausted", "cost_converged", "not_run"]
CONT_NAMES = [
    "mc_long", "mc_smallstep", "mc_smallbox", "vm_long", "vm_smallbox", "rerun_default",
    "from_truth",
]  # fmt: skip
LOG_KEYS = [
    "stop", "steps_run", "n_accept", "last_accept", "n_restarts", "final_step_deg",
    "min_ergodic", "since_improve", "cost_start", "cost_end",
]  # fmt: skip
VM_KEYS = [
    "n_sub", "n_improve", "last_improve", "final_radius_deg", "final_variance",
    "max_steps_end", "steps_taken", "cost_start", "cost_end", "n_restart", "capped",
]  # fmt: skip


class LoggedMC(MCOptimizer):
    """MCOptimizer whose optimize / variance_minimizing_optimize record why and when they stop.
    optimize is the library's (it logs MCOptimizer.last_run); variance_minimizing_optimize is a
    line-for-line copy of orientation_search.MCOptimizer's (same RNG calls, same arithmetic), so
    the results are identical; tests/test_finisher_diagnosis.py verifies this."""

    def __init__(self, *a: Any, **k: Any) -> None:
        super().__init__(*a, **k)
        self.mc_logs: List[Dict[str, float]] = []
        self.vm_logs: List[Dict[str, float]] = []
        # Safety cap on total VM steps (None = off, as in the reconstructor). The VM budget
        # grows by n_subregion_steps whenever the subregion variance exceeds the threshold, so
        # it has no upper bound; the cap is only used for the continuation runs.
        self.vm_step_cap: Optional[int] = None

    def optimize(  # type: ignore[override]
        self,
        initial_orientation: np.ndarray,
        angular_box_side: float,
        angular_step: float,
        max_mc_steps: int,
        max_restarts: int,
        max_convergence_cost: float = 0.0,
        trajectory: Optional[List[Dict]] = None,
    ) -> SearchCandidate:
        """MCOptimizer.optimize itself (the C++-faithful block loop); the stop reason and counts
        it leaves in last_run are appended to mc_logs (stop 0 budget, 1 restarts, 2 converged)."""
        res = super().optimize(
            initial_orientation,
            angular_box_side,
            angular_step,
            max_mc_steps,
            max_restarts,
            max_convergence_cost,
            trajectory,
        )
        self.mc_logs.append(dict(self.last_run))
        return res

    def variance_minimizing_optimize(  # type: ignore[override]
        self,
        initial_orientation: np.ndarray,
        search_box_side: float,
        max_mc_steps: int,
        successive_restarts: int,
        max_convergence_cost: float,
        convergence_variance: float,
        cost_fn_angular_resolution: float = math.radians(0.5),
    ) -> SearchCandidate:
        global_best_q = matrix_to_quaternion(initial_orientation)
        initial_q = global_best_q.copy()
        current_q = global_best_q.copy()
        global_best_info = self.cost_fn.evaluate(
            initial_orientation, self.voxel_vertices, self.phase_index
        )
        global_min_cost = global_best_info.cost
        cost_start = global_min_cost
        subregion_radius = math.tan(search_box_side) / math.sqrt(48.0)
        total_steps_taken = 0
        max_steps = max_mc_steps
        n_global_restarts = 0
        n_sub = n_improve = 0
        last_improve = -1
        variance = float("nan")
        capped = 0
        while total_steps_taken < max_steps:
            if self.vm_step_cap is not None and total_steps_taken >= self.vm_step_cap:
                capped = 1
                break
            if cost_fn_angular_resolution > 0:
                n_subregion_steps = int(
                    math.ceil(subregion_radius / cost_fn_angular_resolution) ** 2.7
                )
            else:
                n_subregion_steps = 10
            n_subregion_steps = max(n_subregion_steps, 10)
            current_mat = quaternion_to_matrix(current_q)
            new_orient, new_cost, variance, new_info = self._zero_temp_with_variance(
                current_mat, subregion_radius, n_subregion_steps
            )
            n_sub += 1
            if variance > convergence_variance:
                max_steps += n_subregion_steps
            total_steps_taken += n_subregion_steps
            if new_cost >= global_min_cost:
                subregion_radius = min(2.0 * subregion_radius, search_box_side)
                # restart as C++ and MCOptimizer: +-SubregionRadius about the initial orientation
                rx = self._rng.uniform(-subregion_radius, subregion_radius)
                ry = self._rng.uniform(-subregion_radius, subregion_radius)
                rz = self._rng.uniform(-subregion_radius, subregion_radius)
                restart_q = self._grid_gen.get_near_identity_point(rx, ry, rz)
                current_q = _quat_multiply(restart_q, initial_q)
                n_global_restarts += 1
            else:
                global_min_cost = new_cost
                global_best_q = matrix_to_quaternion(new_orient)
                global_best_info = new_info
                current_q = global_best_q.copy()
                subregion_radius *= 0.5
                n_improve += 1
                last_improve = n_sub - 1
            if global_min_cost < max_convergence_cost and abs(variance) < convergence_variance:
                break
        self.vm_logs.append(
            dict(
                n_sub=n_sub, n_improve=n_improve, last_improve=last_improve,
                final_radius_deg=math.degrees(subregion_radius), final_variance=variance,
                max_steps_end=max_steps, steps_taken=total_steps_taken, cost_start=cost_start,
                cost_end=global_min_cost, n_restart=n_global_restarts, capped=capped,
            )
        )  # fmt: skip
        return SearchCandidate(
            orientation=quaternion_to_matrix(global_best_q),
            cost=global_min_cost,
            overlap_info=global_best_info,
        )


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------


def geodesic_points(R_from: np.ndarray, R_to: np.ndarray, ts: np.ndarray) -> np.ndarray:
    """Rotations R(t) = exp(t log(R_to R_from^T)) R_from, shape (len(ts), 3, 3): the straight
    (geodesic) path, R(0) = R_from and R(1) = R_to."""
    rv = Rotation.from_matrix(R_to @ R_from.T).as_rotvec()
    return Rotation.from_rotvec(np.asarray(ts)[:, None] * rv[None]).as_matrix() @ R_from


def angle_deg(R1: np.ndarray, R2: np.ndarray) -> float:
    """Angle (deg) of the rotation R1 R2^T (symmetry-agnostic)."""
    return float(np.degrees(np.linalg.norm(Rotation.from_matrix(R1 @ R2.T).as_rotvec())))


def tiny_rotations(
    R: np.ndarray, angles_deg: np.ndarray, n_dirs: int, rng: np.random.Generator
) -> np.ndarray:
    """R rotated by each angle about n_dirs random axes: (len(angles), n_dirs, 3, 3)."""
    ax = rng.normal(size=(n_dirs, 3))
    ax /= np.linalg.norm(ax, axis=1, keepdims=True)
    rv = np.radians(angles_deg)[:, None, None] * ax[None]
    rots = (
        Rotation.from_rotvec(rv.reshape(-1, 3)).as_matrix().reshape(len(angles_deg), n_dirs, 3, 3)
    )
    return rots @ R


# ---------------------------------------------------------------------------
# Case selection
# ---------------------------------------------------------------------------


def select_cases(
    pipe: str, raw: Dict[str, np.ndarray], sweep_vox: np.ndarray
) -> List[Tuple[int, int, List[int]]]:
    """Stratified by radius (index 0..8): N_CASES[pipe] cases split as evenly as possible over the
    radii, each radius over VOX_PER_RADIUS[pipe] random voxels (valid cases only). Returns
    [(vpos, ri, [directions])], vpos indexing the sweep's voxel order."""
    ran = reorder(raw["ran"], raw["voxel_indices"], sweep_vox)[..., VI]
    fb = reorder(raw["fallback"], raw["voxel_indices"], sweep_vox)[..., VI]
    fail = reorder(raw["fail_pass1"], raw["voxel_indices"], sweep_vox)[..., VI]
    valid = ran & ~fb & (fail == 0)  # (voxel, radius, dir)
    rng = np.random.default_rng(SELECT_SEED + (pipe == "H0"))
    n_tot, nv = N_CASES[pipe], VOX_PER_RADIUS[pipe]
    out: List[Tuple[int, int, List[int]]] = []
    for ri in range(N_RADII):
        n_r = n_tot // N_RADII + (1 if ri < n_tot % N_RADII else 0)
        per = [n_r // nv + (1 if k < n_r % nv else 0) for k in range(nv)]
        ok = [v for v in range(valid.shape[0]) if valid[v, ri].sum() >= max(per)]
        for k, v in enumerate(rng.choice(ok, size=nv, replace=False)):
            dirs = rng.choice(np.nonzero(valid[v, ri])[0], size=per[k], replace=False)
            out.append((int(v), ri, sorted(int(d) for d in dirs)))
    assert sum(len(d) for _, _, d in out) == n_tot
    return out


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------


def init_worker(wargs: Dict[str, Any]) -> None:
    fs.init_worker(wargs)


def _cost(R: np.ndarray, vctx: Any) -> float:
    lf = fs._W.local_fn
    return float(lf.evaluate(R.astype(np.float32), vctx.vertices, vctx.voxel.phase).cost)


def _pad(log: Optional[Dict[str, float]], keys: List[str]) -> np.ndarray:
    if log is None:
        return np.full(len(keys), np.nan)
    return np.array([log[k] for k in keys], dtype=np.float64)


def run_finisher(
    start: np.ndarray, vctx: Any, seed: int
) -> Tuple[np.ndarray, float, LoggedMC, Dict[str, Any]]:
    """The unmodified refine_from_candidates([start]) with a LoggedMC (stdout silenced)."""
    rec, lf = fs._W.rec, fs._W.local_fn
    mc = LoggedMC(
        cost_fn=lf, voxel_vertices=vctx.vertices, phase_index=vctx.voxel.phase,
        rng=np.random.default_rng(seed),
    )  # fmt: skip
    with contextlib.redirect_stdout(io.StringIO()):
        res = rec.refine_from_candidates(
            [SearchCandidate(orientation=start.astype(np.float32), cost=1.0)],
            vctx.vertices,
            vctx.voxel.phase,
            local_cost_fn=lf,
            mc_optimizer=mc,
        )
    return np.asarray(res.orientation, dtype=np.float64), float(res.cost), mc, rec.last_find_optimal


def final_box(rec: Any) -> Tuple[float, float]:
    """(box width, MC step) in rad exactly as refine_from_candidates computes them by default."""
    p = rec.params
    diameter = p.local_grid_radius / 1.5 ** (p.max_local_resolution + 1)
    radius = max(diameter / 3.0, math.radians(0.2))
    box = radius / (2**p.min_local_resolution)
    return box, box * p.mc_radius_scale_factor


def continuation(
    kind: str, R0: np.ndarray, vctx: Any, seed: int, R_true: np.ndarray
) -> Dict[str, Any]:
    """One continuation from R0; the optimizers are called as components with changed parameters."""
    rec, lf = fs._W.rec, fs._W.local_fn
    p = rec.params
    box, step = final_box(rec)
    if kind in ("rerun_default", "from_truth"):
        R, c, mc, _ = run_finisher(R0, vctx, seed)
        mc_log, vm_log = mc.mc_logs[0], mc.vm_logs[0]
        c = _cost(R, vctx)
    else:
        mc = LoggedMC(
            cost_fn=lf, voxel_vertices=vctx.vertices, phase_index=vctx.voxel.phase,
            rng=np.random.default_rng(seed),
        )  # fmt: skip
        n_long = p.max_mc_steps * LONG_FACTOR
        mc.vm_step_cap = 10 * n_long  # VM only: 10x the continuation budget
        mc_log = vm_log = None
        if kind.startswith("mc_"):
            b, s, n, rs = {
                "mc_long": (box, step, n_long, 10 * p.successive_restarts),
                "mc_smallstep": (box, step / SMALL_STEP_DIV, 10000, 5),
                "mc_smallbox": (box / SMALL_BOX_DIV, step / SMALL_BOX_DIV, 10000, 5),
            }[kind]
            r = mc.optimize(R0.astype(np.float32), b, s, n, rs, 0.0)
            mc_log = mc.mc_logs[0]
        else:
            b, n = {"vm_long": (box, n_long), "vm_smallbox": (box / SMALL_BOX_DIV, n_long)}[kind]
            r = mc.variance_minimizing_optimize(
                R0.astype(np.float32), b, n, p.successive_restarts, 0.0, 0.02**2
            )
            vm_log = mc.vm_logs[0]
        R = np.asarray(r.orientation, dtype=np.float64)
        c = _cost(R, vctx)
    return dict(R=R, cost=c, ang_true=angle_deg(R, R_true), mc_log=mc_log, vm_log=vm_log)


def diagnose_case(
    vctx: Any, start: np.ndarray, R_res_raw: np.ndarray, seed: int, case_rng_seed: int
) -> Dict[str, np.ndarray]:
    lf = fs._W.local_fn
    R_true = np.asarray(vctx.R_true, dtype=np.float64)
    out: Dict[str, np.ndarray] = {}
    n0 = lf.eval_count
    # --- (d) re-run of the finisher with the Task 1 seed ---
    R_res, cost_rep, mc, fo_info = run_finisher(start, vctx, seed)
    out["rerun_maxabs"] = np.array(np.abs(R_res - R_res_raw).max())
    out["rerun_evals"] = np.array(lf.eval_count - n0)
    out["cost_true"] = np.array(_cost(R_true, vctx))
    out["cost_res"] = np.array(_cost(R_res, vctx))
    out["cost_start"] = np.array(_cost(start, vctx))
    out["ang_start_true"] = np.array(angle_deg(start, R_true))
    out["ang_res_true"] = np.array(angle_deg(R_res, R_true))
    out["fo_log"] = _pad(mc.mc_logs[0], LOG_KEYS)
    out["vm_log"] = _pad(mc.vm_logs[0], VM_KEYS)
    out["fo_converged_flag"] = np.array(bool(fo_info.get("converged", False)))
    for tag, R in (("truth", R_true), ("res", R_res)):
        info = lf.evaluate(R.astype(np.float32), vctx.vertices, vctx.voxel.phase)
        out[f"counts_{tag}"] = np.array(
            [info.pixel_overlap, info.pixel_on_detector, info.peak_overlap,
             info.peak_on_detector, info.n_quality_points]
        )  # fmt: skip
    print(f"    case seed {seed}: finisher re-run done, evals {out['rerun_evals']}", flush=True)
    # --- (a) geodesic result -> truth ---
    ts = np.unique(np.concatenate([np.linspace(0, 1, N_GEO), GEO_NEAR]))
    pts = geodesic_points(R_res, R_true, ts)
    out["geo_t"] = ts
    out["geo_cost"] = np.array([_cost(R, vctx) for R in pts])
    # --- (c) granularity at the truth and at the result ---
    rng = np.random.default_rng(case_rng_seed)
    for tag, R in (("truth", R_true), ("res", R_res)):
        rots = tiny_rotations(R, GRAN_ANGLES_DEG, GRAN_DIRS, rng)
        out[f"gran_{tag}"] = np.array(
            [[_cost(rots[a, k], vctx) for k in range(GRAN_DIRS)] for a in range(len(rots))]
        )
    out["eval_repeat_diff"] = np.array(abs(_cost(R_true, vctx) - float(out["cost_true"])))
    # --- (b) continuations ---
    for ci, kind in enumerate(CONT_NAMES):
        print(f"    case seed {seed}: continuation {kind} ({time.time():.0f})", flush=True)
        R0 = R_true if kind == "from_truth" else R_res
        r = continuation(kind, R0, vctx, seed + 100003 * (ci + 1), R_true)
        out[f"c_{kind}_cost"] = np.array(r["cost"])
        out[f"c_{kind}_ang"] = np.array(r["ang_true"])
        out[f"c_{kind}_moved"] = np.array(angle_deg(r["R"], R0))
        out[f"c_{kind}_mc"] = _pad(r["mc_log"], LOG_KEYS)
        out[f"c_{kind}_vm"] = _pad(r["vm_log"], VM_KEYS)
    return out


def task(item: Tuple[Any, ...]) -> Tuple[int, int, float]:
    vidx, vpos, ri, dirs, pipe, path, ref_nroi, ref_fail, raw_R, raw_start = item
    assert fs._W is not None
    ctx = fs._W.ctx
    a = ctx.args
    t_start = time.time()
    res: Dict[int, Dict[str, np.ndarray]] = {}
    for vb in nnrun.variant_batches(ctx, vidx, vpos, ri, ref_nroi, ref_fail, ["all"]):
        for j in dirs:
            keys = nnrun.case_keys(ctx, vb, j)
            fs.attach_images(keys)
            seed = nnrun.b_seed(a, vpos, ri, j, vb.vi)
            res[j] = diagnose_case(vb.vctx, raw_start[j], raw_R[j], seed, seed + 7)
    out: Dict[str, np.ndarray] = {"dirs": np.array(dirs)}
    for k in res[dirs[0]]:
        out[k] = np.stack([res[j][k] for j in dirs])
    out.update(vidx=np.array(vidx), ri=np.array(ri), pipe=np.array(pipe))
    np.savez(path, **out)
    return vidx, ri, time.time() - t_start


def build_items(cache: Path, only_first: int = 0) -> Tuple[List[Tuple[Any, ...]], Dict[str, Any]]:
    wargs, sweep = nnrun.worker_args(["realistic_s0"])
    wargs["use_models"] = []
    sweep_vox = np.asarray(sweep["voxel_indices"])
    items: List[Tuple[Any, ...]] = []
    for pipe in ("H3", "H0"):
        raw = dict(np.load(RAW_DIR / f"{pipe}_realistic_s0_raw.npz"))
        Rf = reorder(raw["R_final"], raw["voxel_indices"], sweep_vox)[:, :, :, VI]
        St = reorder(raw["start"], raw["voxel_indices"], sweep_vox)[:, :, :, VI]
        cases = select_cases(pipe, raw, sweep_vox)
        if only_first:
            cases = cases[:: max(1, len(cases) // only_first)][:only_first]
        for vpos, ri, dirs in cases:
            if only_first:
                dirs = dirs[:2]
            items.append(
                (int(sweep_vox[vpos]), vpos, ri, dirs, pipe,
                 str(cache / f"{pipe}_v{sweep_vox[vpos]}_r{ri}.npz"),
                 sweep["n_roi"][vpos, ri], sweep["fail_pass1"][vpos, ri],
                 Rf[vpos, ri], St[vpos, ri])
            )  # fmt: skip
    return items, wargs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pilot", "run"])
    ap.add_argument("--workers", type=int, default=10)
    args = ap.parse_args()
    import multiprocessing as mp

    cache = CACHE_DIR / ("pilot" if args.cmd == "pilot" else "full")
    cache.mkdir(parents=True, exist_ok=True)
    items, wargs = build_items(cache, only_first=3 if args.cmd == "pilot" else 0)
    todo = [it for it in items if not Path(it[5]).exists()]
    print(
        f"{len(items)} tasks ({len(items) - len(todo)} cached), {args.workers} workers", flush=True
    )
    t0 = time.time()
    with mp.get_context("spawn").Pool(
        min(args.workers, max(len(todo), 1)), initializer=init_worker, initargs=(wargs,)
    ) as pool:
        for k, (vidx, ri, secs) in enumerate(pool.imap_unordered(task, todo), 1):
            print(
                f"  [{k}/{len(todo)}] voxel {vidx} r#{ri} {secs:.0f}s ({time.time()-t0:.0f}s)",
                flush=True,
            )
    print(f"finished; wall {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
