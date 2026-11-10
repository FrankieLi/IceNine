"""Task E review fix 1: wall time of whole voxels, SmoothCostRank off against on.

One process, one torch thread. Per voxel, the no-start search (``reconstruct_voxel``, lenient
rule, the probe's seeds and settings, ``smooth_rank_probe.task`` without the probe's counting
wrappers) runs with the key off and with each carried arm, in that order, ``--repeats`` times
(interleaved); the wall time of the call is recorded, so it includes every build of the
``SmoothCost`` the library does today (one per ``reconstruct_voxel`` call). Each run's final
error is checked against the end-to-end caches (``cache/rank``).

Preflight-gated (``preflight.require_quiet``): without a quiet machine the script stops, unless
``--contended`` is given; then every time is labelled "contended" and the preflight JSON records
why (owner decision D-20261018-N16: the quiet repeat is then listed as not done).

  uv run python scripts/eval_reduction/smooth_voxel_timing.py --set e0 --repeats 2
  uv run python scripts/eval_reduction/smooth_voxel_timing.py --set full
Writes benchmarks/eval_reduction/E/timing/{voxel_wall_<set>.json, preflight_voxel_<set>.json}.
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

OUT = ICE / "benchmarks" / "eval_reduction" / "E" / "timing"
E0_VOXELS = [7336, 3108, 2910, 16182]  # the four E0 voxels of the review's contended check
FULL_VOXELS = [15901]
SPECS = {"off": "base", "E3+E4,tau=0.25": "e34t25", "E1+E2+E3+E4": "e1234a"}


def one(kind: str, vidx: int, vpos: int, spec: str) -> Dict[str, Any]:
    """Wall time and result of one whole-voxel search (the probe's setup, lenient rule)."""
    from csl import reduced_misorientation_deg

    if kind == "e0":
        from optimizer_sweep import voxel_context

        C = P._W["C"]
        W = C.get_worker()
        keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_clean.npz")["keys"]
        C.attach(keys)
        vctx = voxel_context(W.ctx, vidx)
        rec = W.rec
        Rt = np.asarray(vctx.R_true, dtype=np.float64)
        verts, phase = vctx.vertices, int(vctx.voxel.phase)
        rng = C.run_seed(vpos, 0)
    else:
        ctx = P._W["ctx"]
        rec = ctx.rec
        Rt = ctx.truth(vidx)
        verts, phase = ctx.verts(vidx), int(ctx.mic.voxels[vidx].phase)
        rng = np.random.default_rng(1)
    rec.min_sin_eta = 0.0
    sp = rec.setup.search_params
    P._set_rule(sp, "old", spec)
    try:
        t0 = time.perf_counter()
        with contextlib.redirect_stdout(io.StringIO()):
            res = rec.reconstruct_voxel(voxel_vertices=verts, phase_index=phase, rng=rng)
        wall = time.perf_counter() - t0
    finally:
        P._set_rule(sp, "old", "off")
        if kind == "e0":
            rec.min_sin_eta = 0.3
    err = float(reduced_misorientation_deg(np.asarray(res.orientation, dtype=np.float64), Rt))
    ranker = rec.last_smooth_ranker
    cached = P.run_dir(kind, "old", SPECS[spec]) / f"v{vidx}.npz"
    err_cached = float(np.load(cached)["err_deg"]) if cached.exists() else float("nan")
    return dict(
        voxel=vidx,
        spec=spec,
        wall_s=wall,
        err_deg=err,
        err_deg_cached=err_cached,
        same_as_cache=bool(err == err_cached),
        smooth_evals=0 if ranker is None else int(ranker.n_evals),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="e0", choices=["e0", "full"])
    ap.add_argument("--repeats", type=int, default=1)
    ap.add_argument("--contended", action="store_true")
    ap.add_argument("--voxels", default="", help="comma-separated voxel ids (default: the set)")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    pf = preflight.preflight()
    try:
        preflight.require_quiet(info=pf)
        label = "single-worker"
    except preflight.MachineBusyError as e:
        if not a.contended:
            print(f"not quiet: {e}")
            sys.exit(2)
        label = "contended"
        pf = dict(pf, not_quiet=str(e)[:400])
    (OUT / f"preflight_voxel_{a.set}.json").write_text(json.dumps(pf, indent=1) + "\n")
    P._init(a.set)
    pos = {v: p for v, p in P.voxel_list(a.set, "30")}
    voxels = E0_VOXELS if a.set == "e0" else FULL_VOXELS
    if a.voxels:
        voxels = [int(v) for v in a.voxels.split(",")]
    runs: List[Dict[str, Any]] = []
    t_all = time.time()
    for v in voxels:
        for r in range(a.repeats):
            for spec in SPECS:
                x = one(a.set, v, pos[v], spec)
                x["repeat"] = r
                runs.append(x)
                print(json.dumps(x), flush=True)
    out = dict(
        set=a.set,
        label=label,
        repeats=a.repeats,
        rule="lenient",
        note="whole reconstruct_voxel wall time, one process, one torch thread; includes the "
        "SmoothCost build the library does on every reconstruct_voxel call",
        runs=runs,
        wall_total_s=time.time() - t_all,
    )
    (OUT / f"voxel_wall_{a.set}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(f"done: {len(runs)} runs, {label}, {time.time() - t_all:.0f}s")


if __name__ == "__main__":
    main()
