"""AR3 step 3, Reading 0: what an offset costs (uses the truth; for analysis only, in no decision).

For each voxel of the gradient sample: one strict score of its true orientation and one of its
grain's centre orientation (cost and hit ratio). Written to ``gradient/reading0.json.gz``.

  uv run python scripts/autoresearch/ar3/reading0.py [--limit N]   # N: first N voxels (pilot)
"""

import argparse
import json
import time
from pathlib import Path

import ar3_common as C
import gradient_common as G
import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default=str(G.GRAD / "reading0.json.gz"))
    a = ap.parse_args()
    out = G.GRAD / "runs" / "reading0"
    out.mkdir(parents=True, exist_ok=True)
    setup, bfs, _t = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, extra="PaintGrains 0", stack_name=G.STACK
    )
    from icenine.reconstructor import _get_voxel_vertices

    T = G.Truth()
    mic = setup.sample.get_mic()
    assert len(mic.voxels) == T.n
    cf = bfs.reconstructor._make_local_cost_fn()
    n = a.limit or T.n
    rows = []
    t0 = time.time()
    for v in range(n):
        vox = mic.voxels[v]
        vert = _get_voxel_vertices(vox)
        i1 = cf.evaluate(T.R[v], vert, vox.phase)
        i2 = cf.evaluate(T.R_centre[v], vert, vox.phase)
        h1, _ = bfs._ratios(i1)
        h2, _ = bfs._ratios(i2)
        rows.append([v, int(T.gmap[v]), float(T.level[v]), float(i1.cost), h1, float(i2.cost), h2])
    wall = time.time() - t0
    G.dump_gz(
        dict(
            cols=["voxel", "grain", "level", "cost_true", "hit_true", "cost_centre", "hit_centre"],
            rows=rows,
            wall_s=wall,
            n_evals=2 * n,
        ),
        Path(a.out),
    )
    print(json.dumps(dict(n=n, wall_s=wall, per_voxel_s=wall / n)))


if __name__ == "__main__":
    main()
