"""AR3 step 8, part 1.1: replay of the stored seed steps, split into costing and refinements.

  RESEED_SAMPLE=big uv run python scripts/autoresearch/ar3/reseed_replay.py --run S5 \
      --shard 0 --of 4 [--limit N]

The seed step reads only the voxel's observed rays and uses ``default_rng(0)`` (it does not read the
BFS state), so it can be run again at any voxel (as ``seed_replay.py`` of step 6 did). For every
stored seed step (voxel, reseed) of run ``S5`` (20 candidates costed, 5 refined) or ``D4`` (5 and 2)
the library's ``seed_by_indexing`` runs again and one row is kept:

  costed (number of strict cost calls), strict_costs, refined (candidate indices, in order),
  fin_evals, fin_costs, fin_hit (hit ratio of each refinement's end), fin_err (its error to the
  truth, degrees; analysis only), taken, evals (total), n_cand, n_att.

The hit ratio and the error are evaluated after the seed step and add no evaluation to it.
The row's ``evals`` is compared with the stored run's by the summary (``stored_evals``).
Writes ``reseed/cache/replay_<sample>_<run>_<shard>of<n>.json.gz``.
"""

import argparse
import time
from typing import Any, Dict, List

import numpy as np
import reseed_common as RC


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, choices=["S5", "D4"])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    steps = RC.seed_steps(a.run)
    steps = sorted(steps, key=lambda s: (s["voxel"], s["reseed"]))
    steps = [s for i, s in enumerate(steps) if i % a.of == a.shard]
    if a.limit:
        steps = steps[: a.limit]
    env = RC.make_env()
    rec = env.hooks.rec
    n_cost, n_ref = RC.NCOST_REFINE[a.run]
    T = env.T
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for i, s in enumerate(steps):
        v = s["voxel"]
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        r = rec.index_seed(verts, ph, n_cost=n_cost, n_refine=n_ref)
        hits = []
        errs = []
        for R in r.finisher_R:
            info = rec.evaluate_overlap(np.asarray(R, np.float64), verts, ph)
            hits.append(float(env.bfs._ratios(info)[0]))
            errs.append(float(T.err(np.asarray(R, np.float64)[None], [v])[0]))
        rows.append(
            dict(
                voxel=v,
                reseed=s["reseed"],
                stored_evals=s["stored_evals"],
                grid=s["grid"],
                evals=int(r.evals),
                costed=len(r.strict_costs),
                strict_costs=[float(c) for c in r.strict_costs],
                refined=[int(j) for j in r.refined],
                fin_evals=[int(e) for e in r.finisher_evals],
                fin_costs=[float(c) for c in r.finisher_costs],
                fin_hit=hits,
                fin_err=errs,
                taken=int(r.taken),
                n_cand=int(r.n_cand),
                n_att=int(r.n_att),
            )
        )
        if (i + 1) % 10 == 0 or i + 1 == len(steps):
            print(
                f"{a.run} shard {a.shard}: {i + 1}/{len(steps)} {time.time() - t0:.0f}s", flush=True
            )
    name = f"replay_{RC.SAMPLE}_{a.run}_{a.shard}of{a.of}"
    if a.limit:
        name += f"_limit{a.limit}"
    RC.write_gz(rows, RC.CACHE / f"{name}.json.gz")
    (RC.CACHE / f"{name}.done").write_text(f"{time.time() - t0:.1f}\n")


if __name__ == "__main__":
    main()
