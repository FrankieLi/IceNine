#!/usr/bin/env python3
# flake8: noqa
"""
Small probe of why the (1+1)-ES of bench.py stalls (4 T5 H3 realistic cases, 2000 evaluations):
strict acceptance (as benchmarked), ties accepted, and a 0.002 deg step floor. Prints, per case and
variant, the error at 250 and 2000 evaluations, the evaluation of the last success and the step
at evaluations 100 and 250. Four cases only: an illustration, not a statistic.

Usage (from icenine_py/):
  uv run python scripts/finisher_bench/es_probe.py > benchmarks/finisher_bench/es_probe.txt
"""
import math
import sys

import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench as B
from bench import D, FO, fs, nnrun
from icenine.orientation_search import QuaternionGrid, matrix_to_quaternion, quaternion_to_matrix
from scipy.spatial.transform import Rotation

scr = Path("scripts/finisher_bench/cache/probe")
(scr / "t5_tmp").mkdir(parents=True, exist_ok=True)
items, wargs = B.build_task_items(scr, True, 0)
fs.init_worker(wargs)
ctx, rec, lf = fs._W.ctx, fs._W.rec, fs._W.local_fn
a = ctx.args
box, step0 = D.final_box(rec)


def ang(R, Rt):
    return float(np.degrees(Rotation.from_matrix(R @ Rt.T).magnitude()))


def es_logged(f, R0, rng, step0_rad, accept_ties=False, floor=math.radians(2e-4)):
    grid = QuaternionGrid()
    q = matrix_to_quaternion(R0)
    cur = f(R0)
    step = step0_rad
    log = []
    last_succ = 0
    while True:
        tq = FO._propose(q, step, rng, grid)
        c = f(quaternion_to_matrix(tq))
        ok = c < cur
        if ok or (accept_ties and c <= cur):
            cur, q = c, tq
        if ok:
            last_succ = f.n
        step *= math.exp((float(ok) - 0.25) / (3 * 0.75))
        step = min(max(step, floor), math.radians(2.0))
        log.append((f.n, math.degrees(step), ok))
        f.log = log
        f.last = last_succ


done = 0
for it in items:
    if it["kind"] != "H3":
        continue
    for vb in nnrun.variant_batches(
        ctx, it["vidx"], it["vpos"], it["ri"], it["ref_nroi"], it["ref_fail"], ["all"]
    ):
        for j in it["dirs"][:2]:
            start = np.asarray(it["raw_start"][j], dtype=np.float64)
            keys = nnrun.case_keys(ctx, vb, j)
            fs.attach_images(keys)
            seed = nnrun.b_seed(a, it["vpos"], it["ri"], j, B.VI_REAL)
            Rt = np.asarray(vb.vctx.R_true, float)
            for label, kw in (
                ("strict", {}),
                ("ties", dict(accept_ties=True)),
                ("floor0.002", dict(floor=math.radians(0.002))),
            ):
                cc = FO.CountingCost(lf, vb.vctx.vertices, vb.vctx.voxel.phase, 2000, (250,))
                FO.run_budgeted(
                    lambda c: es_logged(
                        c, start, np.random.default_rng(seed + B.CFG_SEED["es_box"]), step0, **kw
                    ),
                    cc,
                )
                L = cc.log
                steps_at = {n: s for n, s, _ in L}
                late = sum(ok for n, _, ok in L if n > 250)
                msg = (
                    f"v{it['vidx']} j{j} {label}: start {ang(start, Rt):.4f} "
                    f"err@250 {ang(cc.snap[250][1], Rt):.4f} err@2000 {ang(cc.best_R, Rt):.4f} "
                    f"last success at eval {cc.last} step@100 {steps_at.get(100, 0):.5f} "
                    f"step@250 {steps_at.get(250, 0):.5f} succ in 250-2000 {late}"
                )
                print(msg, flush=True)
            done += 1
    if done >= 4:
        break
