"""Task G1 step 2c: whole-voxel wall time on the full sample with the SmoothCost cache.

One process, one torch thread, lenient rule, the probe's settings (``smooth_rank_probe``: clean
full-sample images, Q-max 8 configuration, seed 1, ``min_sin_eta`` 0). Two voxels: A = 15901 (the
voxel Task E timed) and B = 17034. Per carried arm, in this order:

  1. cache cleared, A          cold: the call builds the SmoothCost (the as-shipped cost)
  2. B                         warm: served by the cache built on A (the BFS case)
  3. B, ``smooth_cost_cache`` None   fresh: builds its own SmoothCost (behaviour before the cache)
  4. A                         warm again

with the key off for A and B first. Steps 2 and 3 must give the same search bit for bit (final
orientation, cost, evaluation counts, smooth evaluations); 1 and 4 too. The SmoothCost build
itself is timed once on its own (constructor, lit-pixel trees, observed rays).

Preflight-gated (``preflight.require_quiet``), retried once a minute up to ``--retries`` times;
then it stops, unless ``--contended`` is given (every time is then labelled "contended").

  uv run python scripts/eval_reduction/g1_cache_timing.py --retries 20
Writes benchmarks/eval_reduction/G/g1/timing/{cache_timing.json, preflight.json}.
"""

import argparse
import contextlib
import io
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import preflight  # noqa: E402
import smooth_rank_probe as P  # noqa: E402

OUT = ICE / "benchmarks" / "eval_reduction" / "G" / "g1" / "timing"
VOX_A, VOX_B = 15901, 17034
ARMS = {"E3+E4,tau=0.25": "e34t25", "E1+E2+E3+E4": "e1234a"}


def one(vidx: int, spec: str, cached: bool, step: str) -> Dict[str, Any]:
    from csl import reduced_misorientation_deg
    from icenine.smooth_cost import SmoothCostCache

    ctx = P._W["ctx"]
    rec = ctx.rec
    Rt = ctx.truth(vidx)
    verts, phase = ctx.verts(vidx), int(ctx.mic.voxels[vidx].phase)
    rng = np.random.default_rng(1)
    rec.min_sin_eta = 0.0
    keep = rec.smooth_cost_cache
    if not cached:
        rec.smooth_cost_cache = None
    sp = rec.setup.search_params
    P._set_rule(sp, "old", spec)
    builds0 = None if keep is None else keep.n_builds
    try:
        t0 = time.perf_counter()
        with contextlib.redirect_stdout(io.StringIO()):
            res = rec.reconstruct_voxel(voxel_vertices=verts, phase_index=phase, rng=rng)
        wall = time.perf_counter() - t0
    finally:
        P._set_rule(sp, "old", "off")
        rec.smooth_cost_cache = keep if keep is not None else SmoothCostCache()
    R = np.asarray(res.orientation, dtype=np.float64)
    ranker = rec.last_smooth_ranker
    g, loc, tot = rec.last_eval_counts
    return dict(
        voxel=vidx,
        spec=spec,
        step=step,
        cached=cached,
        built=None if (not cached or keep is None) else int(keep.n_builds - (builds0 or 0)),
        wall_s=wall,
        err_deg=float(reduced_misorientation_deg(R, Rt)),
        R_final=R.tolist(),
        cost=float(res.cost),
        evals_total=int(tot),
        smooth_evals=0 if ranker is None else int(ranker.n_evals),
    )


def build_parts() -> Dict[str, float]:
    """Wall time of a SmoothCost build on its own: constructor, trees (E1-E3), rays (E4)."""
    from icenine.smooth_cost import SmoothCost

    s = P._W["ctx"].setup
    t0 = time.perf_counter()
    sc = SmoothCost(
        s.simulator,
        s.detector_list,
        s.range_map,
        s.exp_data,
        s.sample,
        s.structure_list,
        eta_limit=s.exp_setup.get_eta_limit(),
        min_sin_eta=0.0,
    )
    t1 = time.perf_counter()
    sc._ensure_trees()
    t2 = time.perf_counter()
    sc.observed_rays()
    t3 = time.perf_counter()
    return dict(constructor_s=t1 - t0, trees_s=t2 - t1, rays_s=t3 - t2, total_s=t3 - t0)


def same(a: Dict[str, Any], b: Dict[str, Any]) -> bool:
    return bool(
        a["R_final"] == b["R_final"]
        and a["cost"] == b["cost"]
        and a["evals_total"] == b["evals_total"]
        and a["smooth_evals"] == b["smooth_evals"]
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--retries", type=int, default=0, help="preflight retries, one a minute")
    ap.add_argument("--contended", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    label = "contended"
    pf: Dict[str, Any] = {}
    for i in range(a.retries + 1):
        pf = preflight.preflight()
        try:
            preflight.require_quiet(info=pf)
            label = "single-worker"
            break
        except preflight.MachineBusyError as e:
            print(f"preflight {i + 1}: {str(e)[:300]}", flush=True)
            pf = dict(pf, not_quiet=str(e)[:400])
            if i < a.retries:
                time.sleep(60)
    if label != "single-worker" and not a.contended:
        (OUT / "preflight.json").write_text(json.dumps(pf, indent=1) + "\n")
        print("not quiet; stopping (use --contended to run labelled contended)")
        sys.exit(2)
    pf["attempts"] = i + 1
    (OUT / "preflight.json").write_text(json.dumps(pf, indent=1) + "\n")
    P._init("full")
    t_all = time.time()
    parts = build_parts()
    print(json.dumps(dict(build=parts)), flush=True)
    runs: List[Dict[str, Any]] = []

    def go(v: int, spec: str, cached: bool, step: str) -> Dict[str, Any]:
        x = one(v, spec, cached, step)
        runs.append(x)
        print(json.dumps({k: x[k] for k in x if k != "R_final"}), flush=True)
        return x

    go(VOX_A, "off", True, "off A")
    go(VOX_B, "off", True, "off B")
    checks: Dict[str, Any] = {}
    for spec in ARMS:
        P._W["ctx"].rec.smooth_cost_cache.clear()
        a1 = go(VOX_A, spec, True, "1 cold A")
        b2 = go(VOX_B, spec, True, "2 warm B")
        b3 = go(VOX_B, spec, False, "3 fresh B")
        a4 = go(VOX_A, spec, True, "4 warm A")
        checks[spec] = dict(
            warm_B_equals_fresh_B=same(b2, b3),
            cold_A_equals_warm_A=same(a1, a4),
            builds=[a1["built"], b2["built"], a4["built"]],
        )
    out = dict(
        label=label,
        rule="lenient",
        voxels=[VOX_A, VOX_B],
        note="whole reconstruct_voxel wall time, one process, one torch thread",
        build=parts,
        runs=runs,
        checks=checks,
        wall_total_s=time.time() - t_all,
    )
    (OUT / "cache_timing.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(checks, indent=1))
    print(f"done: {len(runs)} runs, {label}, {time.time() - t_all:.0f}s")


if __name__ == "__main__":
    main()
