#!/usr/bin/env python3
"""
Phase B3 timing: single-worker wall time of every finisher on a fixed 40-case subset (realistic
variant: 8 H3 + 4 H0 T5 tasks and 8 sweep tasks, 2 cases each), interleaved across methods, gated
by preflight.require_quiet(). The budgeted methods are timed as separate runs at budgets 250 and
2600 (a run to its budget, not a truncation of a longer one); mc_deployed, gn, adam and the default
finisher at their natural stop. If the machine is not quiet nothing is timed: timing.json then
holds status "skipped" and the reasons.

Usage (from icenine_py/; nothing else should be running):
  uv run python scripts/finisher_bench/finisher_timing.py [--out benchmarks/finisher_bench]
"""

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bench as B  # noqa: E402  (sets OMP/MKL threads to 1 and the import paths)
from bench import D, FO, fs, nnrun, ps  # noqa: E402

import preflight  # noqa: E402

N_TASKS = {"H3": 8, "H0": 4, "SW": 8}
DIRS_PER_TASK = 2
TIMED_BUDGETS = [250, 2600]


def select(items: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rng = np.random.default_rng(20261007)
    out = []
    for kind, n in N_TASKS.items():
        pool = [it for it in items if it["kind"] == kind and len(it["dirs"]) >= DIRS_PER_TASK]
        pick = rng.choice(len(pool), size=n, replace=False)
        out += [pool[int(i)] for i in sorted(pick)]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(B.ICENINE_PY / "benchmarks" / "finisher_bench"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    info = preflight.preflight()
    (out / "preflight.json").write_text(json.dumps(info, indent=1))
    try:
        preflight.require_quiet(info=info)
    except preflight.MachineBusyError as e:
        (out / "timing.json").write_text(
            json.dumps(dict(status="skipped", reason=str(e)), indent=1)
        )
        print(f"SKIPPED: {e}")
        return

    cache = B.CACHE_DIR / "timing_tmp"
    (cache / "t5_tmp").mkdir(parents=True, exist_ok=True)
    items, wargs = B.build_task_items(cache, False, 0)
    chosen = select(items)
    fs.init_worker(wargs)
    ctx, rec, lf = fs._W.ctx, fs._W.rec, fs._W.local_fn
    a = ctx.args
    box, step0 = D.final_box(rec)

    rec_t: Dict[str, List[float]] = {}
    rec_n: Dict[str, List[float]] = {}
    case_no = 0
    names = ["mc_deployed"] + [m for m in B.ANYTIME if m != "mc_deployed"]
    for it in chosen:
        vidx, vpos, ri = it["vidx"], it["vpos"], it["ri"]
        r_told = float(ps.RADII_DEG[ri])
        dirs = it["dirs"][:DIRS_PER_TASK]
        for vb in nnrun.variant_batches(
            ctx, vidx, vpos, ri, it["ref_nroi"], it["ref_fail"], ["all"]
        ):
            for j in dirs:
                start = np.asarray(
                    it["raw_start"][j] if it["raw_start"] is not None else vb.R_nom0[j],
                    dtype=np.float64,
                )
                keys = nnrun.case_keys(ctx, vb, j)
                fs.attach_images(keys)
                seed = nnrun.b_seed(a, vpos, ri, j, B.VI_REAL)
                methods = B.make_methods(box, step0, rec.params, math.radians(r_told))
                # the work items of this case, rotated so no method always runs first
                work = [
                    (m, b) for m in names for b in ([None] if m == "mc_deployed" else TIMED_BUDGETS)
                ]
                work += [("gn", None), ("adam", None), ("finisher", None)]
                k = case_no % len(work)
                for m, b in work[k:] + work[:k]:
                    key = f"{m}@{b if b else 'natural'}"
                    n0 = lf.eval_count
                    t0 = time.perf_counter()
                    if m == "gn":
                        B.gn_from_starts(ctx, vidx, vb, {j: start})
                        n = 3  # window renders + solves
                    elif m == "adam":
                        B.run_adam(ctx, vb.vctx, keys, start, rec.params, box, seed + 104729)
                        n = lf.eval_count - n0
                    elif m == "finisher":
                        nnrun.refine_fo(start.astype(np.float32), vb.vctx, seed)
                        n = lf.eval_count - n0
                    else:
                        cc = FO.CountingCost(
                            lf, vb.vctx.vertices, vb.vctx.voxel.phase, b or B.MAX_BUDGET, ()
                        )
                        FO.run_budgeted(
                            lambda c, f=methods[m], mm=m: f(c, start, seed + B.CFG_SEED[mm]), cc
                        )
                        n = cc.n
                    dt = time.perf_counter() - t0
                    rec_t.setdefault(key, []).append(dt)
                    rec_n.setdefault(key, []).append(float(n))
                case_no += 1
                print(f"  timed case {case_no}", flush=True)
    res: Dict[str, Any] = {}
    for key, ts in rec_t.items():
        m, b = key.split("@")
        t = np.array(ts)
        n = np.array(rec_n[key])
        res[m + ("" if b == "natural" else f"@{b}")] = dict(
            budget=b if b != "natural" else "natural",
            evals_median=float(np.median(n)),
            wall_median_s=float(np.median(t)),
            wall_q25_s=float(np.quantile(t, 0.25)),
            wall_q75_s=float(np.quantile(t, 0.75)),
            ms_per_eval=float(1000 * np.median(t / np.maximum(n, 1))),
            n_cases=int(len(t)),
        )
    info_after = preflight.preflight()
    (out / "timing.json").write_text(
        json.dumps(
            dict(
                status="single-worker, preflight passed before the run",
                n_cases=case_no,
                preflight_after=info_after,
                methods=res,
            ),
            indent=1,
        )
    )
    print("done", flush=True)


if __name__ == "__main__":
    main()
