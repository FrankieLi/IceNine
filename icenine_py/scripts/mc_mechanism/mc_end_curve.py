#!/usr/bin/env python3
"""
Phase B1 follow-up: the improvement-probability curve at the MC stage's own OUTPUT.

mc_trace.py measures the curve at the finisher's final result (after the VarianceMinimizing stage).
The MC stage ends earlier and elsewhere, so this script re-runs only the first stage, the
MCOptimizer.optimize call of refine_from_candidates (same box, step, steps, restarts and seed; the
rng is not used before it), asserts that its log equals the one stored by mc_trace.py, and
estimates P(a proposal of step s lowers the cost) at the MC output over a grid that includes the
steps the runs actually end on, and at the run's own final step.

Usage (from icenine_py/):
  uv run python scripts/mc_mechanism/mc_end_curve.py pilot
  uv run python scripts/mc_mechanism/mc_end_curve.py run --workers 10
"""

import argparse
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, Tuple

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
import mc_trace as T  # noqa: E402
import run as nnrun  # noqa: E402

CACHE_DIR = HERE / "cache"
VARIANTS = T.VARIANTS
VI_REAL = T.VI_REAL
STEPS_DEG = np.array(
    [
        0.00013,
        0.00026,
        0.0005,
        0.001,
        0.002,
        0.003,
        0.005,
        0.0075,
        0.01,
        0.02,
        0.03,
        0.05,
        0.075,
        0.1,
        0.15,
        0.2,
        0.3,
    ]
)
N_PROP = 400
KEYS = ("p_improve", "cost_prog", "dist_prog", "n_improve")


def stats_at(
    R: np.ndarray, step_deg: float, R_true: np.ndarray, vctx: Any, rng: Any, grid: Any
) -> Dict[str, float]:
    lf = fs._W.local_fn
    # float64 like the optimizer's trial matrices; the float32 value is kept as a check
    c0 = float(lf.evaluate(R, vctx.vertices, vctx.voxel.phase).cost)
    c32 = float(lf.evaluate(R.astype(np.float32), vctx.vertices, vctx.voxel.phase).cost)
    d0 = D.angle_deg(R, R_true)
    mats, _ = H.mc_proposals(R, math.radians(step_deg), N_PROP, rng, grid)
    c = np.array([lf.evaluate(m, vctx.vertices, vctx.voxel.phase).cost for m in mats], dtype=float)
    d = np.array([D.angle_deg(m, R_true) for m in mats])
    out = H.expected_progress(c0, d0, c, d)
    out["c0"], out["c0_cast_diff"] = c0, c32 - c0
    return out


def task(item: Tuple[Any, ...]) -> Tuple[int, int, float]:
    vidx, vpos, ri, dirs, pipe, path, ref_nroi, ref_fail, raw_R, raw_start = item
    assert fs._W is not None
    ctx, rec = fs._W.ctx, fs._W.rec
    a = ctx.args
    t0 = time.time()
    box, step = D.final_box(rec)
    p = rec.params
    stored = np.load(
        Path(path).parent.parent / Path(path).parent.name[len("end_") :] / Path(path).name
    )[
        "log"
    ]  # (variant, case, 10): the mc_trace.py cache of the same name without the "end_" prefix
    from icenine.orientation_search import QuaternionGrid

    grid = QuaternionGrid()
    res: Dict[Tuple[int, int], Dict[str, np.ndarray]] = {}
    for vb in nnrun.variant_batches(ctx, vidx, vpos, ri, ref_nroi, ref_fail, VARIANTS):
        for jj, j in enumerate(dirs):
            fs.attach_images(nnrun.case_keys(ctx, vb, j))
            R_true = np.asarray(vb.vctx.R_true, dtype=np.float64)
            seed = nnrun.b_seed(a, vpos, ri, j, VI_REAL)
            lf = fs._W.local_fn
            mc = T.TracedMC(
                cost_fn=lf,
                voxel_vertices=vb.vctx.vertices,
                phase_index=vb.vctx.voxel.phase,
                rng=np.random.default_rng(seed),
            )
            mc.R_true = R_true
            out = mc.optimize(
                raw_start[j].astype(np.float32),
                box,
                step,
                p.max_mc_steps,
                p.successive_restarts,
                p.max_convergence_cost,
            )
            log = D._pad(mc.mc_logs[0], D.LOG_KEYS)
            assert np.array_equal(
                log, stored[vb.vi, jj], equal_nan=True
            ), "MC stage differs from the stored run"
            R_mc = np.asarray(out.orientation, dtype=np.float64)
            prng = np.random.default_rng(seed + 55001)
            cur = {k: np.zeros(len(STEPS_DEG)) for k in KEYS}
            for si, s in enumerate(STEPS_DEG):
                st = stats_at(R_mc, float(s), R_true, vb.vctx, prng, grid)
                for k in KEYS:
                    cur[k][si] = st[k]
            own = stats_at(
                R_mc, float(log[5]), R_true, vb.vctx, prng, grid
            )  # log[5] = final step (deg)
            if log[2] > 0:  # improving runs: the output's cost is the logged cost_end exactly
                assert own["c0"] == log[9], (own["c0"], log[9])
            r: Dict[str, np.ndarray] = dict(
                c0_cast_diff=np.array(own["c0_cast_diff"]),
                R_mc=R_mc,
                dist_mc=np.array(D.angle_deg(R_mc, R_true)),
                own_step=np.array(log[5]),
                n_accept=np.array(log[2]),
            )
            for k in KEYS:
                r["c_" + k] = cur[k]
                r["own_" + k] = np.array(own[k])
            res[(vb.vi, jj)] = r
    out_d: Dict[str, np.ndarray] = dict(
        dirs=np.array(dirs),
        steps_deg=STEPS_DEG,
        n_prop=np.array(N_PROP),
        vidx=np.array(vidx),
        ri=np.array(ri),
        pipe=np.array(pipe),
    )
    for k in res[(0, 0)]:
        out_d[k] = np.stack(
            [np.stack([res[(vi, jj)][k] for jj in range(len(dirs))]) for vi in range(len(VARIANTS))]
        )
    np.savez(path, **out_d)
    return vidx, ri, time.time() - t0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pilot", "run"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--tag", default="", help="rerun into cache/end_<cmd>_<tag>")
    args = ap.parse_args()
    import multiprocessing as mp

    sfx = f"_{args.tag}" if args.tag else ""
    cache = CACHE_DIR / ("end_" + args.cmd + sfx)
    cache.mkdir(parents=True, exist_ok=True)
    items, wargs = D.build_items(cache, only_first=2 if args.cmd == "pilot" else 0)
    if args.cmd == "pilot":  # the stored mc_trace pilot cache has the same file names
        items = [it for it in items if (CACHE_DIR / ("pilot" + sfx) / Path(it[5]).name).exists()]
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
