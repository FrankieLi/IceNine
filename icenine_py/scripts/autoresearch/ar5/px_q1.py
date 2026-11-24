"""AR5 step 27, gate q1: does the finally accepted orientation stay among the top 20 candidates?

  RESEED_SAMPLE=<sample> OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 uv run python \
      scripts/autoresearch/ar5/px_q1.py [--limit N]

The voxels are those of ``ar5/px/voxels.json`` key ``q1`` (frozen): 300 stored seed steps of the
sample (generator 27). For every voxel the old indexer (``index_rays``, as ``index_voxel`` calls
it) runs with ``max_pairs`` 600,000 (base) and 50,000 (new), the order alternating per voxel.
Rank: position of the first candidate within 1 degree (symmetry reduced) of the orientation the
stored run's seed step finally accepted (the step's best fit by hit ratio, as the sizing scripts
did), -1 when none of the 20 is. No evaluation is made; the truth is not opened. Writes
``<px>/q1_<sample>.json``.
"""

import argparse
import dataclasses
import json
import sys
import time
from typing import Any, Dict, List

import px_common as PC
import px_rules as PR

import numpy as np  # noqa: E402


def final_R(e: Dict[str, Any]) -> Any:
    best = max(e["fits"], key=lambda f: f[3])
    return np.array(best[6]).reshape(3, 3), float(best[3])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    s = PC.SAMPLE
    q = json.loads((PC.RD / "voxels.json").read_text())["q1"][s]
    e1 = PC.load_gz(PC.CACHE / f"J5D_{s}" / "e1_log.json.gz")
    steps = q["steps"][: a.limit] if a.limit else q["steps"]
    S = PC.setup()
    env, rec = S["env"], S["rec"]
    from icenine import ray_indexing as RI

    rows: List[Dict[str, Any]] = []
    t_start = time.time()
    for b, i in enumerate(steps):
        e = e1[i]
        v = int(e["voxel"])
        assert v == q["voxels"][b], "the stored step does not match the frozen voxel list"
        Rf, hit = final_R(e)
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        sc, crystal, params, _ = rec._index_context(ph)
        att = sc.attached_rays(verts)
        g = RI.scattering_vectors(sc, att)
        row: Dict[str, Any] = dict(step=int(i), voxel=v, reseed=bool(e["reseed"]), hit=hit)
        order = ["base", "new"] if b % 2 == 0 else ["new", "base"]
        for name in order:
            mp = PR.OLD_MAX_PAIRS if name == "base" else PR.NEW_MAX_PAIRS
            p = dataclasses.replace(params, max_pairs=mp)
            t0 = time.perf_counter()
            res = RI.index_rays(g, crystal, p, RI.SEED_TOP_K)
            dt = time.perf_counter() - t0
            if len(res.R):
                d = RI.misorientation_deg(res.R, np.broadcast_to(Rf, res.R.shape), crystal.sym)
                hp = np.nonzero(d <= 1.0)[0]
                rank, best_d = (int(hp[0]) if len(hp) else -1), float(d.min())
            else:
                rank, best_d = -1, 999.0
            row[name] = dict(
                rank=rank,
                best_d=best_d,
                n_pairs=int(res.n_pairs),
                n_votes=int(res.n_votes),
                n_cand=int(len(res.R)),
                s=dt,
            )
        row["n_rays"] = int(len(att))
        rows.append(row)
        if b % 50 == 0:
            print(
                s,
                b,
                v,
                row["base"]["rank"],
                row["new"]["rank"],
                f"{time.time() - t_start:.0f}s",
                flush=True,
            )
    base_in = [r for r in rows if r["base"]["rank"] >= 0]
    lost = [r for r in base_in if r["new"]["rank"] < 0]
    gained = [r for r in rows if r["base"]["rank"] < 0 and r["new"]["rank"] >= 0]
    out = dict(
        sample=s,
        n=len(rows),
        was_in_top20=len(base_in),
        lost=len(lost),
        gained=len(gained),
        lost_voxels=[r["voxel"] for r in lost],
        gained_voxels=[r["voxel"] for r in gained],
        rows=rows,
        wall_s=time.time() - t_start,
    )
    name = f"q1_{s}{'_limit' if a.limit else ''}.json"
    (PC.RD / name).write_text(json.dumps(out) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}))


if __name__ == "__main__":
    main()
