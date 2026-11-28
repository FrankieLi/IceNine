"""AR3 step 6: facts about the voxels of part 0 that need the sample but not a search.

For every fallback voxel and control: whether it lies at a true grain boundary, its true level and
its true grain. Truth is used for analysis only. The reconstructor's neighbour relation reaches
about 29 voxels, so it is too wide for "at a boundary"; here a voxel is at a true grain boundary
when a voxel within 1.5 times the median nearest-neighbour distance (the nearest ring) belongs to
another true grain. Writes ``seed/cache/voxel_facts.json``.
"""

import json

import ar3_common as C  # noqa: F401
import boundary_common as B
import numpy as np
import seed_common as S
from scipy.spatial import cKDTree


def main() -> None:
    env = B.Env("grad")
    xy = np.array([v.position[:2] for v in env.mic.voxels], dtype=float)
    tree = cKDTree(xy)
    dd = tree.query(xy, k=8)[0]
    d1 = float(np.median([row[row > 1e-9][0] for row in dd]))  # coincident positions skipped
    vox = set()
    for line in (S.CACHE / "part01.jsonl").read_text().splitlines():
        vox.add(int(json.loads(line)["voxel"]))
    for f in S.all_fallbacks():
        vox.add(f["voxel"])
    out = {}
    for v in sorted(vox):
        ring = [j for j in tree.query_ball_point(xy[v], 1.5 * d1) if j != v]
        out[str(v)] = dict(
            boundary=bool(np.any(env.gmap[ring] != env.gmap[v])),
            n_ring=len(ring),
            level=float(env.level[v]),
            grain=int(env.gmap[v]),
            isA=bool(env.isA[v]),
        )
    (S.SEED / "voxel_facts.json").write_text(json.dumps(out) + "\n")
    (S.SEED / "lists.json").write_text(json.dumps(S.counts()) + "\n")
    n_b = sum(x["boundary"] for x in out.values())
    print(len(out), "voxels", n_b, "at a boundary", "d1", d1)


if __name__ == "__main__":
    main()
