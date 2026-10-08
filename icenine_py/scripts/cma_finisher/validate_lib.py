#!/usr/bin/env python3
"""Phase C validation 1: does the library CMAOptimizer reproduce B3's `cma_02` numbers?

Cases: the first N_CASES (50) T5 H3 cases of Phase B3, realistic variant (same start, same images,
same seed stream as scripts/finisher_bench/bench.py). Per case:
  * lib250 / lib1000: icenine.orientation_search.CMAOptimizer (sigma0 0.2 deg, max_evals 250 /
    1000, seed = the B3 seed for cma_02) from the case's start;
  * b3_250 / b3_1000: B3's stored result of the same case (res_R of cma_02 at the checkpoints);
  * refine1000: refine_from_candidates([start]) with SearchParameters.local_optimizer = "cma"
    (what the switch does inside the reconstructor; one candidate, so one CMA run, then the
    reconstructor's final overlap evaluation).
Error = angle between the result and the truth (no symmetry reduction, local search).

  uv run python scripts/cma_finisher/validate_lib.py [--workers 10] [--n-cases 50]
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]
sys.path.insert(0, str(ICENINE_PY / "scripts" / "finisher_bench"))
import bench as B  # noqa: E402  (also puts nn_hybrid, finisher_diagnosis, common on sys.path)

fs, nnrun = B.fs, B.nnrun
sys.path.insert(0, str(ICENINE_PY / "scripts" / "common"))
from stats import wilson  # noqa: E402

OUT_DIR = ICENINE_PY / "benchmarks" / "cma_finisher"
CACHE = HERE / "cache" / "validate_lib"
B3_CACHE = B.CACHE_DIR / "run"
CMA_IDX = B.ANYTIME.index("cma_02")
K250, K1000 = B.CKPTS.index(250), B.CKPTS.index(1000)
N_CASES = 50


def angle_deg(R: np.ndarray, Rt: np.ndarray) -> float:
    return float(np.degrees(Rotation.from_matrix(R @ Rt.T).magnitude()))


def task(item: Dict[str, Any]) -> str:
    from icenine.orientation_search import CMAOptimizer, SearchCandidate

    ctx = fs._W.ctx
    a = ctx.args
    rec, lf = fs._W.rec, fs._W.local_fn
    vidx, vpos, ri, dirs = item["vidx"], item["vpos"], item["ri"], item["dirs"]
    stored = np.load(item["b3_path"])
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for vb in nnrun.variant_batches(
        ctx, vidx, vpos, ri, item["ref_nroi"], item["ref_fail"], [B.VARIANTS[B.VI_REAL]]
    ):
        assert vb.vi == B.VI_REAL
        for jj, j in enumerate(dirs):
            keys = nnrun.case_keys(ctx, vb, j)
            fs.attach_images(keys)
            seed = nnrun.b_seed(a, vpos, ri, j, B.VI_REAL)
            R0 = np.asarray(item["raw_start"][j], dtype=np.float64)
            Rt = np.asarray(vb.vctx.R_true, dtype=np.float64)
            vv, ph = vb.vctx.vertices, vb.vctx.voxel.phase
            row: Dict[str, Any] = dict(vidx=vidx, j=int(j), start_err=angle_deg(R0, Rt))
            cma_seed = seed + B.CFG_SEED["cma_02"]
            for budget in (250, 1000):
                n0 = lf.eval_count
                res = CMAOptimizer(lf, vv, ph, sigma0_deg=0.2, max_evals=budget).optimize(
                    R0, seed=cma_seed
                )
                assert lf.eval_count - n0 == res.n_evals <= budget
                R_b3 = stored["res_R"][B.VI_REAL, jj, CMA_IDX, K250 if budget == 250 else K1000]
                row[f"lib{budget}_err"] = angle_deg(res.orientation, Rt)
                row[f"lib{budget}_cost"] = float(res.cost)
                row[f"lib{budget}_evals"] = int(res.n_evals)
                row[f"lib{budget}_stop"] = res.stop_reason
                row[f"b3_{budget}_err"] = angle_deg(R_b3, Rt)
                row[f"lib{budget}_vs_b3_maxabs"] = float(np.abs(res.orientation - R_b3).max())
            rec.params.local_optimizer, rec.params.cma_max_evals = "cma", 1000
            try:
                n0 = lf.eval_count
                out = rec.refine_from_candidates(
                    [SearchCandidate(orientation=R0, cost=1.0)], vv, ph,
                    rng=np.random.default_rng(seed),
                )  # fmt: skip
                row["refine_err"] = angle_deg(out.orientation, Rt)
                row["refine_evals"] = int(rec.last_eval_counts[1])
            finally:
                rec.params.local_optimizer = "mc"
            # BFS neighbour refinement call site: local_optimization from the same start
            from icenine.cost_functions import VoxelCostFunction

            orig_eval = VoxelCostFunction.evaluate
            count = [0]

            def counted(self: Any, *args: Any, **kw: Any) -> Any:
                count[0] += 1
                return orig_eval(self, *args, **kw)

            VoxelCostFunction.evaluate = counted  # type: ignore[method-assign]
            try:
                for mode in ("mc", "cma"):
                    rec.params.local_optimizer = mode
                    count[0] = 0
                    t1 = time.perf_counter()
                    lo = rec.local_optimization(vv, ph, R0, rng=np.random.default_rng(seed))
                    row[f"lo_{mode}_secs"] = time.perf_counter() - t1
                    row[f"lo_{mode}_err"] = angle_deg(lo.orientation, Rt)
                    row[f"lo_{mode}_evals"] = int(count[0])
                    row[f"lo_{mode}_moved_deg"] = angle_deg(lo.orientation, R0)
            finally:
                VoxelCostFunction.evaluate = orig_eval  # type: ignore[method-assign]
                rec.params.local_optimizer = "mc"
            rows.append(row)
    (CACHE / f"v{vidx}_r{ri}.json").write_text(json.dumps(rows))
    return f"voxel {vidx} r#{ri} {len(rows)} cases {time.time() - t0:.0f}s"


def summarize(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    n = len(rows)
    out: Dict[str, Any] = dict(n_cases=n)
    for name in ("start", "b3_250", "lib250", "b3_1000", "lib1000", "refine"):
        e = np.array([r[f"{name}_err"] for r in rows])
        k = int((e < 0.02).sum())
        lo, hi = wilson(k, n)
        out[name] = dict(
            median_err_deg=float(np.median(e)), n_below_0p02=k, frac_below_0p02=k / n,
            wilson_lo=lo, wilson_hi=hi, n_wrong_gt1=int((e > 1.0).sum()),
        )  # fmt: skip
    for b in (250, 1000):
        d = np.array([r[f"lib{b}_vs_b3_maxabs"] for r in rows])
        out[f"lib{b}_vs_b3"] = dict(
            n_identical=int((d == 0).sum()), max_abs_diff=float(d.max()),
            median_evals=float(np.median([r[f"lib{b}_evals"] for r in rows])),
        )  # fmt: skip
    out["lib1000_stop_reasons"] = {
        s: sum(1 for r in rows if r["lib1000_stop"] == s)
        for s in sorted({r["lib1000_stop"] for r in rows})
    }
    out["refine_median_evals"] = float(np.median([r["refine_evals"] for r in rows]))
    mv = np.array([r["lo_mc_moved_deg"] for r in rows])
    d_err = np.array([r["lo_mc_err"] - r["start_err"] for r in rows])
    out["local_optimization_mc_outcome"] = dict(
        unchanged_to_1e5_deg=int((mv < 1e-5).sum()),
        moved_and_closer=int(((mv >= 1e-5) & (d_err < 0)).sum()),
        moved_and_farther=int(((mv >= 1e-5) & (d_err >= 0)).sum()),
    )
    for mode in ("mc", "cma"):
        e = np.array([r[f"lo_{mode}_err"] for r in rows])
        k = int((e < 0.02).sum())
        lo_, hi_ = wilson(k, n)
        out[f"local_optimization_{mode}"] = dict(
            median_err_deg=float(np.median(e)), n_below_0p02=k, wilson_lo=lo_, wilson_hi=hi_,
            median_evals=float(np.median([r[f"lo_{mode}_evals"] for r in rows])),
            mean_evals=float(np.mean([r[f"lo_{mode}_evals"] for r in rows])),
            n_wrong_gt1=int((e > 1.0).sum()),
            mean_secs_contended=float(np.mean([r[f"lo_{mode}_secs"] for r in rows])),
        )  # fmt: skip
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--n-cases", type=int, default=N_CASES)
    args = ap.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    items, wargs = B.build_task_items(B3_CACHE, False, 0)
    h3 = sorted((it for it in items if it["kind"] == "H3"), key=lambda it: (it["vidx"], it["ri"]))
    sel: List[Dict[str, Any]] = []
    left = args.n_cases
    for it in h3:
        if left <= 0:
            break
        it = dict(it)
        it["dirs"] = it["dirs"][:left]
        left -= len(it["dirs"])
        it["b3_path"] = it["path"]
        sel.append(it)
    todo = [it for it in sel if not (CACHE / f"v{it['vidx']}_r{it['ri']}.json").exists()]
    print(
        f"{len(sel)} tasks ({len(todo)} to run), {sum(len(i['dirs']) for i in sel)} cases",
        flush=True,
    )
    import multiprocessing as mp

    t0 = time.time()
    with mp.get_context("spawn").Pool(
        min(args.workers, max(len(todo), 1)), initializer=fs.init_worker, initargs=(wargs,)
    ) as pool:
        for k, msg in enumerate(pool.imap_unordered(task, todo), 1):
            print(f"  [{k}/{len(todo)}] {msg} ({time.time() - t0:.0f}s)", flush=True)
    rows: List[Dict[str, Any]] = []
    for it in sel:
        rows += json.loads((CACHE / f"v{it['vidx']}_r{it['ri']}.json").read_text())
    summ = summarize(rows)
    summ["note"] = (
        "T5 H3 realistic, first cases by (voxel, radius index); contended 10-worker run, no timing"
    )
    (OUT_DIR / "validate_lib.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
