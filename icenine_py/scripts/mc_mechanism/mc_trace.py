#!/usr/bin/env python3
"""
Phase B1 (finisher/MC study): per-step MC trace and improvement-probability curves on the T5 cases.

Cases are exactly the T5 / Phase A cases (diagnose.build_items: 200 H3 + 100 H0), both variants
(clean, realistic) of the same voxel/radius/direction, the same starts, and the realistic variant's
seed for both (so the clean run is a paired re-run of the same random stream).

Per case and variant:
  (1) the unmodified refine_from_candidates with TracedMC (a LoggedMC whose optimize also
      records, per step, the step size, min_ergodic, steps since the last global improvement, the
      event, the trial rotation angle and the distance of the best orientation to the truth). The
      realistic variant's result is asserted bit-identical to the T5 result.
  (2) the improvement-probability curve: at the finisher's result, the truth, and the truth
      rotated by 0.005/0.01/0.02/0.05 deg (2 directions each), for each MC step size s, N
      proposals drawn exactly as MCOptimizer draws them; the cost of each is evaluated.

Usage (from icenine_py/):
  uv run python scripts/mc_mechanism/mc_trace.py pilot
  uv run python scripts/mc_mechanism/mc_trace.py run --workers 10
"""

import argparse
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "finisher_diagnosis"))
sys.path.insert(0, str(HERE.parent / "nn_hybrid"))
sys.path.insert(0, str(HERE.parent / "common"))
sys.path.insert(0, str(HERE.parent))

import diagnose as D  # noqa: E402
import findoptimal_sweep as fs  # noqa: E402
import mc_helpers as H  # noqa: E402
import run as nnrun  # noqa: E402

from icenine.orientation_search import (  # noqa: E402
    SearchCandidate,
    _quat_multiply,
    matrix_to_quaternion,
    quaternion_to_matrix,
)

CACHE_DIR = HERE / "cache"
# The T5 result (stored with the pre-port MC and VarianceMinimizing) is bit-identical only for the
# original B1 cache; a tagged rerun (C++-faithful MC and VM) records the difference instead.
TAGGED = os.environ.get("MC_TRACE_TAG", "") != ""
VARIANTS = ["clean", "all"]  # "all" = realistic
VI_REAL = 1
STEPS_DEG = np.array(
    [0.0005, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3]
)
OFFSETS_DEG = np.array([0.005, 0.01, 0.02, 0.05])
DIRS_PER_OFFSET = 2
N_PROP = 200
PT_NAMES = ["result", "truth"] + [
    f"d{o:g}_{k}" for o in OFFSETS_DEG for k in range(DIRS_PER_OFFSET)
]


class TracedMC(D.LoggedMC):
    """LoggedMC whose optimize also records a per-step trace (arrays of length max_mc_steps,
    event = EV_NOT_RUN after the run ended). The body is a copy of LoggedMC.optimize (same RNG
    calls, same arithmetic); mc_trace.py asserts the realistic result equals the T5 result."""

    def __init__(self, *a: Any, **k: Any) -> None:
        super().__init__(*a, **k)
        self.R_true: Optional[np.ndarray] = None
        self.traces: List[Dict[str, np.ndarray]] = []

    def _ang(self, q: np.ndarray) -> float:
        return (
            D.angle_deg(quaternion_to_matrix(q), self.R_true)
            if self.R_true is not None
            else math.nan
        )

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
        """MCOptimizer.optimize (library code) with per-step arrays filled by the _mc_block and
        _on_block_end hooks. Step t of a block is trial t; the block's last step carries the
        block outcome (EV_GLOBAL: global best improved and step halved, EV_RESTART: failed block,
        EV_EXHAUSTED: failed block that ended the run); other steps are EV_LOCAL when the trial
        lowered the block's own best, else EV_NONE. n_since is the successive-failure count after
        the step's block (constant inside a block), best_* the global best at the block's start
        (updated on the block's last step)."""
        n = max_mc_steps
        self._max_restarts = max_restarts
        self._tr = dict(
            step_deg=np.full(n, np.nan),
            min_erg=np.full(n, np.nan),
            n_since=np.full(n, np.nan),
            event=np.full(n, H.EV_NOT_RUN, dtype=np.int8),
            trial_ang=np.full(n, np.nan),
            trial_cost=np.full(n, np.nan),
            cur_cost=np.full(n, np.nan),
            best_ang=np.full(n, np.nan),
            best_cost=np.full(n, np.nan),
        )
        self._n_min = H.min_ergodic(angular_box_side, angular_step, max_mc_steps)
        self._best_ang = self._ang(matrix_to_quaternion(initial_orientation))
        self._best_cost = math.nan
        self._blk = (0, 0)
        res = super().optimize(
            initial_orientation,
            angular_box_side,
            angular_step,
            max_mc_steps,
            max_restarts,
            max_convergence_cost,
            trajectory,
        )
        self.traces.append(self._tr)
        return res

    def _mc_block(  # type: ignore[override]
        self, start_q: np.ndarray, step: float, n_steps: int
    ) -> Tuple[np.ndarray, float, Any]:
        """Copy of MCOptimizer._mc_block that records each trial."""
        tr, lo = self._tr, self._block_start
        opt_q = start_q.copy()
        opt_info = self.cost_fn.evaluate(
            quaternion_to_matrix(opt_q), self.voxel_vertices, self.phase_index
        )
        opt_cost = opt_info.cost
        if math.isnan(self._best_cost):
            self._best_cost = opt_cost
        radius = math.tan(step) / math.sqrt(12.0) if step > 0 else 0.01
        for i in range(n_steps):
            t = lo + i
            tr["step_deg"][t] = math.degrees(step)
            tr["min_erg"][t] = self._n_min
            x = self._rng.uniform(-radius, radius)
            y = self._rng.uniform(-radius, radius)
            z = self._rng.uniform(-radius, radius)
            delta_q = self._grid_gen.get_near_identity_point(x, y, z)
            tr["trial_ang"][t] = math.degrees(2.0 * math.acos(min(1.0, abs(float(delta_q[0])))))
            trial_q = _quat_multiply(delta_q, opt_q)
            trial_info = self.cost_fn.evaluate(
                quaternion_to_matrix(trial_q), self.voxel_vertices, self.phase_index
            )
            tr["trial_cost"][t] = trial_info.cost
            tr["cur_cost"][t] = opt_cost
            ev = H.EV_NONE
            if trial_info.cost < opt_cost:
                opt_q, opt_cost, opt_info = trial_q, trial_info.cost, trial_info
                ev = H.EV_LOCAL
            tr["event"][t] = ev
            tr["best_ang"][t], tr["best_cost"][t] = self._best_ang, self._best_cost
        self._blk = (lo, lo + n_steps)
        return opt_q, opt_cost, opt_info

    def _on_block_end(  # type: ignore[override]
        self,
        event: str,
        block_step: float,
        n_succ_restarts: int,
        global_min_cost: float,
        current_q: np.ndarray,
        best_q: np.ndarray,
    ) -> None:
        tr, (lo, hi) = self._tr, self._blk
        self._best_ang = self._ang(best_q)
        self._best_cost = global_min_cost
        if hi > lo:
            tr["n_since"][lo:hi] = n_succ_restarts
            last = hi - 1
            if event == "mc_accept":
                tr["event"][last] = H.EV_GLOBAL
            elif n_succ_restarts > self._max_restarts:
                tr["event"][last] = H.EV_EXHAUSTED
            else:
                tr["event"][last] = H.EV_RESTART
            tr["best_ang"][last], tr["best_cost"][last] = self._best_ang, self._best_cost


def run_traced_finisher(
    start: np.ndarray, vctx: Any, seed: int
) -> Tuple[np.ndarray, TracedMC, Dict[str, Any]]:
    import contextlib
    import io

    rec, lf = fs._W.rec, fs._W.local_fn
    mc = TracedMC(
        cost_fn=lf,
        voxel_vertices=vctx.vertices,
        phase_index=vctx.voxel.phase,
        rng=np.random.default_rng(seed),
    )
    mc.R_true = np.asarray(vctx.R_true, dtype=np.float64)
    n0 = lf.eval_count
    with contextlib.redirect_stdout(io.StringIO()):
        res = rec.refine_from_candidates(
            [SearchCandidate(orientation=start.astype(np.float32), cost=1.0)],
            vctx.vertices,
            vctx.voxel.phase,
            local_cost_fn=lf,
            mc_optimizer=mc,
        )
    info = dict(evals=lf.eval_count - n0)
    return np.asarray(res.orientation, dtype=np.float64), mc, info


def probability_curve(
    points: np.ndarray, R_true: np.ndarray, vctx: Any, seed: int, grid_gen: Any
) -> Dict[str, np.ndarray]:
    """For each point (n_pt, 3, 3) and step in STEPS_DEG: N_PROP proposals, statistics."""
    lf = fs._W.local_fn
    rng = np.random.default_rng(seed)
    keys = ("p_improve", "cost_prog", "dist_prog", "n_improve")
    out = {k: np.zeros((len(points), len(STEPS_DEG))) for k in keys}
    out["trial_ang_med"] = np.zeros((len(points), len(STEPS_DEG)))
    out["cost0"] = np.zeros(len(points))
    out["dist0"] = np.zeros(len(points))
    out["c0_cast_diff"] = np.zeros(len(points))
    for pi, R in enumerate(points):
        # float64 like the optimizer's trial matrices; the float32 value is kept as a check
        c0 = float(lf.evaluate(R, vctx.vertices, vctx.voxel.phase).cost)
        c32 = float(lf.evaluate(R.astype(np.float32), vctx.vertices, vctx.voxel.phase).cost)
        d0 = D.angle_deg(R, R_true)
        out["cost0"][pi], out["dist0"][pi] = c0, d0
        out["c0_cast_diff"][pi] = c32 - c0
        for si, s in enumerate(STEPS_DEG):
            mats, ang = H.mc_proposals(R, math.radians(s), N_PROP, rng, grid_gen)
            c = np.array(
                [lf.evaluate(m, vctx.vertices, vctx.voxel.phase).cost for m in mats], dtype=float
            )
            d = np.array([D.angle_deg(m, R_true) for m in mats])
            st = H.expected_progress(c0, d0, c, d)
            for k in keys:
                out[k][pi, si] = st[k]
            out["trial_ang_med"][pi, si] = float(np.median(ang))
    return out


def task(item: Tuple[Any, ...]) -> Tuple[int, int, float]:
    vidx, vpos, ri, dirs, pipe, path, ref_nroi, ref_fail, raw_R, raw_start = item
    assert fs._W is not None
    ctx = fs._W.ctx
    a = ctx.args
    t0 = time.time()
    rec = fs._W.rec
    box, step = D.final_box(rec)
    from icenine.orientation_search import QuaternionGrid

    grid_gen = QuaternionGrid()
    res: Dict[Tuple[int, int], Dict[str, np.ndarray]] = {}
    for vb in nnrun.variant_batches(ctx, vidx, vpos, ri, ref_nroi, ref_fail, VARIANTS):
        for jj, j in enumerate(dirs):
            fs.attach_images(nnrun.case_keys(ctx, vb, j))
            R_true = np.asarray(vb.vctx.R_true, dtype=np.float64)
            seed = nnrun.b_seed(a, vpos, ri, j, VI_REAL)  # the T5 seed, for both variants
            R_res, mc, info = run_traced_finisher(raw_start[j], vb.vctx, seed)
            maxabs = float(np.abs(R_res - raw_R[j]).max())
            if vb.vi == VI_REAL and not TAGGED:
                assert maxabs == 0.0, f"T5 result not reproduced: {maxabs}"
            # probability curve points
            prng = np.random.default_rng(seed + 4177)
            ax = prng.normal(size=(len(OFFSETS_DEG) * DIRS_PER_OFFSET, 3))
            ax /= np.linalg.norm(ax, axis=1, keepdims=True)
            offs = np.repeat(OFFSETS_DEG, DIRS_PER_OFFSET)
            from scipy.spatial.transform import Rotation

            near = Rotation.from_rotvec(np.radians(offs)[:, None] * ax).as_matrix() @ R_true
            pts = np.concatenate([R_res[None], R_true[None], near])
            curve = probability_curve(pts, R_true, vb.vctx, seed + 977, grid_gen)
            r: Dict[str, np.ndarray] = dict(
                result=R_res,
                rerun_maxabs=np.array(maxabs),
                rerun_evals=np.array(info["evals"]),
                ang_res_true=np.array(D.angle_deg(R_res, R_true)),
                log=D._pad(mc.mc_logs[0], D.LOG_KEYS),
            )
            for k, v in mc.traces[0].items():
                r["tr_" + k] = v
            for k, v in curve.items():
                r["pc_" + k] = v
            res[(vb.vi, jj)] = r
    out: Dict[str, np.ndarray] = dict(
        dirs=np.array(dirs),
        box_deg=np.array(math.degrees(box)),
        step0_deg=np.array(math.degrees(step)),
        steps_deg=STEPS_DEG,
        offsets_deg=OFFSETS_DEG,
        n_prop=np.array(N_PROP),
        max_mc_steps=np.array(rec.params.max_mc_steps),
        max_restarts=np.array(rec.params.successive_restarts),
        vidx=np.array(vidx),
        ri=np.array(ri),
        pipe=np.array(pipe),
    )
    for k in res[(0, 0)]:
        out[k] = np.stack(
            [np.stack([res[(vi, jj)][k] for jj in range(len(dirs))]) for vi in range(len(VARIANTS))]
        )  # (variant, case, ...)
    np.savez(path, **out)
    return vidx, ri, time.time() - t0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pilot", "run"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--tag", default="", help="rerun into cache/<cmd>_<tag> (C++-faithful MC)")
    args = ap.parse_args()
    import multiprocessing as mp

    if args.tag:
        os.environ["MC_TRACE_TAG"] = args.tag  # inherited by the spawned workers
    cache = CACHE_DIR / (args.cmd + (f"_{args.tag}" if args.tag else ""))
    cache.mkdir(parents=True, exist_ok=True)
    items, wargs = D.build_items(cache, only_first=2 if args.cmd == "pilot" else 0)
    todo = [it for it in items if not Path(it[5]).exists()]
    print(
        f"{len(items)} tasks ({len(items) - len(todo)} cached), {args.workers} workers", flush=True
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
