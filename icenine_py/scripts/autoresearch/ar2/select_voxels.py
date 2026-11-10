"""AR2 step 1: fix the holdout voxel sets (before any indexing or cost run).

Eligible: voxels of the layer that are not in the pilot region. Boundary voxel: it shares a vertex
(so also an edge) with a voxel of a grain more than 1 degree away (symmetry-reduced) in the true
map; interior: not boundary. Each set takes every s-th eligible voxel in voxel-index order,
s = floor(count / 150), stopping at 150. No random draw.

  uv run python scripts/autoresearch/ar2/select_voxels.py
Writes benchmarks/autoresearch/ar2/holdout/voxels.json.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ar1"))
import ar1_common as AC  # noqa: E402
import misses as M  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

N_PER_SET = 150
GRAIN_TOL_DEG = 1.0
OUT = AC.ICE / "benchmarks" / "autoresearch" / "ar2" / "holdout"


def vertex_neighbours(vxy: np.ndarray, side: np.ndarray) -> List[Set[int]]:
    """Voxels that share at least one vertex with each voxel (itself excluded)."""
    q = float(np.median(side)) * 1e-3
    owners: Dict[Tuple[int, int], List[int]] = {}
    keys = np.rint(vxy / q).astype(np.int64)  # (N, 3, 2)
    for v in range(len(vxy)):
        for k in range(3):
            owners.setdefault((int(keys[v, k, 0]), int(keys[v, k, 1])), []).append(v)
    nb: List[Set[int]] = [set() for _ in range(len(vxy))]
    for vs in owners.values():
        for a in vs:
            nb[a].update(b for b in vs if b != a)
    return nb


def main() -> None:
    gmap = np.load(M.GRAINMAP)
    vxy, side = M.voxel_vertices_xy(M.MIC)
    assert len(vxy) == len(gmap)
    region = np.load(AC.ICE / "benchmarks" / "phase_d_pilot" / "region_index.npy")
    in_region = np.zeros(len(gmap), dtype=bool)
    in_region[region] = True
    pilot_grains = set(int(g) for g in np.unique(gmap[region]))
    AC.init("full")
    ctx = AC.bench().RUN._ENV["ctx"]
    first: Dict[int, int] = {}
    for v, g in enumerate(gmap):
        first.setdefault(int(g), v)
    Rg = {g: ctx.truth(v) for g, v in first.items()}
    e0 = AC.env("full", int(np.nonzero(~in_region)[0][0]), need_sc=False)
    sym = RI.proper_rotations(AC.symmetry_of(e0))
    far_cache: Dict[Tuple[int, int], bool] = {}

    def far(a: int, b: int) -> bool:
        if a == b:
            return False
        k = (min(a, b), max(a, b))
        if k not in far_cache:
            far_cache[k] = bool(
                RI.misorientation_deg(Rg[a][None], Rg[b][None], sym)[0] > GRAIN_TOL_DEG
            )
        return far_cache[k]

    nb = vertex_neighbours(vxy, side)
    boundary = np.array(
        [any(far(int(gmap[v]), int(gmap[u])) for u in nb[v]) for v in range(len(gmap))]
    )
    elig = np.nonzero(~in_region)[0]
    sets: Dict[str, dict] = {}
    for name, flag in (("boundary", True), ("interior", False)):
        pool = [int(v) for v in elig if boundary[v] == flag]
        s = len(pool) // N_PER_SET
        chosen = pool[::s][:N_PER_SET]
        sets[name] = dict(
            eligible=len(pool),
            stride=s,
            n=len(chosen),
            voxels=chosen,
            grain=[int(gmap[v]) for v in chosen],
            grain_in_pilot_region=[int(gmap[v]) in pilot_grains for v in chosen],
        )
    out = dict(
        rule="every s-th eligible voxel in index order, s = floor(count / 150), stop at 150",
        random_draw=False,
        layer_voxels=int(len(gmap)),
        pilot_region_voxels=int(len(region)),
        eligible_voxels=int(len(elig)),
        boundary_test="shares a vertex with a voxel of a grain more than 1 degree away",
        grain_tol_deg=GRAIN_TOL_DEG,
        sets=sets,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "voxels.json").write_text(json.dumps(out, indent=1) + "\n")
    for n, d in sets.items():
        print(
            n,
            "eligible",
            d["eligible"],
            "stride",
            d["stride"],
            "n",
            d["n"],
            "grains",
            len(set(d["grain"])),
            "grain in pilot region",
            sum(d["grain_in_pilot_region"]),
        )


if __name__ == "__main__":
    main()
