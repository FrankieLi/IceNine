"""AR3 step 5, part 1: the pass on one stored map (every try of every boundary voxel).

  uv run python scripts/autoresearch/ar3/boundary_pass.py --map D4 [--limit N]
  uv run python scripts/autoresearch/ar3/boundary_pass.py --map D4 --input after.npy --tag second

The candidates of every boundary voxel are read from the map as it was before the pass, so the
order of the voxels does not matter. One try (``boundary_common.try_fit``) is made from each
candidate: boundary voxels in index order, candidates in neighbour-index order, one worker, one
generator state (``default_rng(0)``). The images and Q-max are those of the run that made the map
(gradient sample: ``grad_clean.npy``; layer: ``clean.npy``; ``key_on_fused.config``, Q-max 8).
Every voxel's current end is also scored once with the strict cost (cost, hit ratio, confidence):
the keep rules and the skip rule read these scores and nothing else about the truth-free side.

Writes ``boundary/runs/<map>[_<tag>]/``: ``tries.npz`` (one row per try: voxel, candidate voxel,
end orientation as float32, strict cost, hit ratio, confidence, evaluations; and the old end's
scores of all voxels), ``pass.json`` (counts, wall time, label), ``progress.json``,
``preflight.json``.
"""

import argparse
import json
import time
from typing import Any, Dict, List

import ar3_common as C
import boundary_common as B
import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True, choices=sorted(B.ALL_MAPS))
    ap.add_argument("--input", default="", help="a .npy map (n,3,3) instead of the stored one")
    ap.add_argument("--tag", default="")
    ap.add_argument("--limit", type=int, default=0, help="pilot: first N boundary voxels")
    ap.add_argument("--side-by-side", type=int, default=1, help="heavy processes running together")
    a = ap.parse_args()
    name = a.map + (f"_{a.tag}" if a.tag else "")
    out = B.RUNS / (f"pilot_{name}" if a.limit else name)
    out.mkdir(parents=True, exist_ok=True)
    label = C.preflight_label(out)
    if a.side_by_side > 1:
        label += f"; {a.side_by_side} passes ran side by side"
    env = B.Env(B.ALL_MAPS[a.map][0])
    R = np.load(a.input).astype(np.float64) if a.input else B.load_map(a.map, env.n)
    isb, cands = env.boundary(R)
    voxels = sorted(cands)
    if a.limit:
        voxels = voxels[: a.limit]
    rng = np.random.default_rng(0)
    t0 = time.time()
    score = np.zeros((env.n, 3))
    for v in range(env.n):
        score[v] = env.score3(R[v], v)
    t_score = time.time() - t0
    cols: Dict[str, List[Any]] = {k: [] for k in ("voxel", "cand", "cost", "hit", "conf", "evals")}
    ends: List[np.ndarray] = []
    t1 = time.time()
    last = t1
    for i, v in enumerate(voxels):
        for c in cands[v]:
            t = B.try_fit(env.hooks, env.bfs, v, R[c], rng)
            cols["voxel"].append(v)
            cols["cand"].append(c)
            for k in ("cost", "hit", "conf", "evals"):
                cols[k].append(t[k])
            ends.append(t["end"])
        if time.time() - last > 60:
            last = time.time()
            (out / "progress.json").write_text(
                json.dumps(dict(boundary_voxels_done=i + 1, of=len(voxels), s=last - t1)) + "\n"
            )
    wall = time.time() - t1
    np.savez_compressed(
        out / "tries.npz",
        voxel=np.array(cols["voxel"], np.int32),
        cand=np.array(cols["cand"], np.int32),
        end=np.array(ends, np.float32).reshape(-1, 3, 3),
        cost=np.array(cols["cost"]),
        hit=np.array(cols["hit"]),
        conf=np.array(cols["conf"]),
        evals=np.array(cols["evals"], np.int32),
        old_score=score,
        old_R=R.astype(np.float32),
    )
    info = dict(
        map=a.map,
        tag=a.tag,
        n_boundary_voxels=len(voxels),
        n_boundary_all=int(isb.sum()),
        n_tries=len(cols["voxel"]),
        evals=int(np.sum(cols["evals"])),
        score_evals=env.n,
        wall_tries_s=wall,
        wall_score_s=t_score,
        timing_label=label,
        limit=a.limit,
        rng="default_rng(0), one generator for the whole pass",
    )
    (out / "pass.json").write_text(json.dumps(info, indent=2) + "\n")
    print("done", json.dumps(info), flush=True)


if __name__ == "__main__":
    main()
