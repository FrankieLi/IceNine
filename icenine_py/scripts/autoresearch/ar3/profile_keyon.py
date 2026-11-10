"""AR3 step 1, "where the 213 s go": rerun the key-on pilot (SeedByIndexing, FusedCost, one worker)
under cProfile and split the time by part with in-process timers.

  uv run python scripts/autoresearch/ar3/profile_keyon.py --name keyon_profile

The parts (inclusive wall time; ``eval`` = time inside ``VoxelCostFunction.evaluate`` called while
that part was the innermost running):
  local fit      ``AdaptiveVoxelReconstructor.local_optimization``  (eval = the cost evaluations;
                 the rest is the optimizer's own work)
  indexing       ``ray_indexing.index_voxel`` (voting, clustering, refining)
  finisher       ``refine_from_candidates`` (the seed's local runs; eval inside)
  seed path      ``index_seed`` in total (indexing + strict costs on candidates + finisher)
  pairing        the one-time SmoothCost build and ray pairing (``index_pairing_s``)
  overlap        ``evaluate_overlap`` (the BFS's own hit-ratio evaluations)
The profile slows the run (cProfile overhead), so its wall time is not a timing; the shares are
the reading. The run is also checked against the stored key-on-fused run: evaluations must match.
"""

import argparse
import cProfile
import json
import pstats
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List

import ar3_common as C
import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="keyon_profile")
    ap.add_argument("--no-profile", action="store_true", help="timers only (no cProfile)")
    a = ap.parse_args()
    out = C.OUT / "runs" / a.name
    out.mkdir(parents=True, exist_ok=True)
    label = C.preflight_label(out)
    setup, bfs, t_setup = C.build(C.CONFIG, out)

    from icenine import ray_indexing as RI
    from icenine.cost_functions import VoxelCostFunction
    from icenine.reconstructor import AdaptiveVoxelReconstructor as ARec

    incl: Dict[str, float] = {}
    calls: Dict[str, int] = {}
    ev_t: Dict[str, float] = {}
    ev_n: Dict[str, int] = {}
    stack: List[str] = ["other"]

    def wrap(owner: Any, attr: str, name: str) -> None:
        f = getattr(owner, attr)

        def w(*args: Any, **kw: Any) -> Any:
            stack.append(name)
            t0 = time.perf_counter()
            try:
                return f(*args, **kw)
            finally:
                incl[name] = incl.get(name, 0.0) + time.perf_counter() - t0
                calls[name] = calls.get(name, 0) + 1
                stack.pop()

        setattr(owner, attr, w)

    wrap(ARec, "local_optimization", "local_fit")
    wrap(ARec, "refine_from_candidates", "finisher")
    wrap(ARec, "index_seed", "seed_path")
    wrap(ARec, "evaluate_overlap", "overlap")
    wrap(ARec, "reconstruct_voxel", "grid_search")
    wrap(RI, "index_voxel", "indexing")
    ev = VoxelCostFunction.evaluate

    def timed_eval(self: Any, *args: Any, **kw: Any) -> Any:
        k = stack[-1]
        t0 = time.perf_counter()
        try:
            return ev(self, *args, **kw)
        finally:
            ev_t[k] = ev_t.get(k, 0.0) + time.perf_counter() - t0
            ev_n[k] = ev_n.get(k, 0) + 1

    VoxelCostFunction.evaluate = timed_eval  # type: ignore[method-assign]

    rng = np.random.default_rng(0)
    prof = cProfile.Profile()
    t0 = time.time()
    if not a.no_profile:
        prof.enable()
    bfs.reconstruct_sample(output_mic=str(out / "recon.mic"), rng=rng)
    if not a.no_profile:
        prof.disable()
    wall = time.time() - t0
    VoxelCostFunction.evaluate = ev  # type: ignore[method-assign]
    st = bfs.stats
    evals = {k: int(st[f"n_evals_{k}"]) for k in ("seed", "reseed", "neighbor", "revisit")}
    pairing = float(st.get("index_pairing_s", 0.0))
    parts: Dict[str, Any] = {}
    for k in ("local_fit", "finisher", "indexing", "seed_path", "overlap", "grid_search"):
        parts[k] = dict(
            incl_s=incl.get(k, 0.0),
            calls=calls.get(k, 0),
            eval_s=ev_t.get(k, 0.0),
            evals=ev_n.get(k, 0),
        )
    parts["other_evals"] = dict(eval_s=ev_t.get("other", 0.0), evals=ev_n.get("other", 0))
    top: List[Dict[str, Any]] = []
    if not a.no_profile:
        ps = pstats.Stats(prof)
        rows = sorted(ps.stats.items(), key=lambda kv: -kv[1][2])[:30]  # type: ignore[attr-defined]
        for (fn, line, name), (cc, nc, tt, ct, _) in rows:
            top.append(
                dict(func=f"{Path(fn).name}:{line}:{name}", calls=nc, tottime_s=tt, cumtime_s=ct)
            )
    info = dict(
        name=a.name,
        wall_total_s=wall,
        profiled=not a.no_profile,
        timing_label=label,
        setup_s=t_setup,
        pairing_s=pairing,
        evals=evals,
        evals_total=int(sum(evals.values())),
        parts=parts,
        top_by_tottime=top,
        fused_cost=bool(getattr(setup.search_params, "fused_cost", False)),
        torch_threads=C.torch.get_num_threads(),
        stats={k: v for k, v in st.items() if k != "reseeds"},
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    (out / "records.json").write_text(
        json.dumps({str(k): asdict(v) for k, v in bfs.records.items()}) + "\n"
    )
    print(json.dumps({k: v for k, v in info.items() if k != "stats"}, indent=2, default=float))


if __name__ == "__main__":
    main()
