"""Select the Phase D pilot region: the N voxels of the new-orientation sample whose triangle
centroids are nearest to a fixed centre (a disc), reproducible from (centre, N).

Writes to benchmarks/phase_d_pilot/: region_grid.mic (the selected lines of the grid-only mic,
orientations zeroed), region_index.npy (index of each region voxel in the full sample file),
region.json (selection and description). The images still contain the whole sample.

Usage (from icenine_py/): uv run python scripts/phase_d/make_region.py [--cx 0.17 --cy 0 --n 2000]
"""

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
EX = HERE.parents[2] / "Examples" / "Example2.ManyGrains" / "SimInput"
OUT = HERE.parents[1] / "benchmarks" / "phase_d_pilot"
NEAR_AXIS_R = 0.12  # mm (.mic coordinates are in mm); "near the rotation axis"


def centroids_from_lines(raw: np.ndarray, header_side: float) -> Tuple[np.ndarray, float]:
    """Triangle centroids (N, 2) from .mic columns (x, y, z, direction, generation, ...)."""
    side = header_side / 2 ** int(raw[0, 4])
    cy = raw[:, 1] + np.where(raw[:, 3] == 1, 1.0, -1.0) * side * np.sqrt(3.0) / 6.0
    return np.stack([raw[:, 0] + side / 2.0, cy], axis=1), side


def boundary_flags(centroids: np.ndarray, grain: np.ndarray, side: float) -> np.ndarray:
    """True for voxels with a different-grain voxel within 1.01 x side (centroid distance), the
    adjacency of make_neworient.py: the 3 edge neighbours and 6 of the 9 vertex-only neighbours (the
    other 3 sit at 2 side / sqrt(3)). The BFS neighbour radius is 2 sides, wider than this."""
    pairs = cKDTree(centroids).query_pairs(1.01 * side, output_type="ndarray")
    flag = np.zeros(len(grain), dtype=bool)
    diff = grain[pairs[:, 0]] != grain[pairs[:, 1]]
    flag[pairs[diff, 0]] = True
    flag[pairs[diff, 1]] = True
    return flag


def select_disc(centroids: np.ndarray, centre: Tuple[float, float], n: int) -> np.ndarray:
    """Indices of the n voxels nearest to `centre` (stable sort, so ties break by file index),
    returned in file order."""
    d = np.hypot(centroids[:, 0] - centre[0], centroids[:, 1] - centre[1])
    return np.sort(np.argsort(d, kind="stable")[:n])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cx", type=float, default=0.17)
    ap.add_argument("--cy", type=float, default=0.0)
    ap.add_argument("--n", type=int, default=2000)
    a = ap.parse_args()
    src = EX / "rand_500grains_1mm_neworient_s0_grid.mic"
    lines = src.read_text().splitlines()
    body = [ln for ln in lines[1:] if ln.strip()]
    raw = np.array([[float(t) for t in ln.split()[:6]] for ln in body])
    cent, side = centroids_from_lines(raw, float(lines[0].split()[0]))
    grain = np.load(EX / "rand_500grains_1mm_neworient_s0_grainmap.npy")
    assert len(grain) == len(body)
    idx = select_disc(cent, (a.cx, a.cy), a.n)
    # contiguous: one connected component under edge adjacency (centroid distance s/sqrt(3))
    sub = cent[idx]
    e = cKDTree(sub).query_pairs(0.6 * side, output_type="ndarray")
    parent = np.arange(len(idx))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in e:
        parent[find(int(i))] = find(int(j))
    n_comp = len({find(i) for i in range(len(idx))})
    bnd = boundary_flags(cent, grain, side)[idx]
    r = np.hypot(sub[:, 0], sub[:, 1])
    g = grain[idx]
    uniq, cnt = np.unique(g, return_counts=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "region_grid.mic").write_text("\n".join([lines[0]] + [body[i] for i in idx]) + "\n")
    np.save(OUT / "region_index.npy", idx.astype(np.int64))
    info: Dict[str, object] = {
        "selection": f"{a.n} voxels with triangle centroids nearest to ({a.cx}, {a.cy}) mm "
        "(stable sort by distance, ties by file index) of "
        "rand_500grains_1mm_neworient_s0_grid.mic",
        "centre_xy_mm": [a.cx, a.cy],
        "n_voxels": int(len(idx)),
        "side_length_mm": side,
        "disc_radius_mm": float(np.hypot(*(sub - np.array([a.cx, a.cy])).T).max()),
        "n_edge_connected_components": int(n_comp),
        "n_grains": int(len(uniq)),
        "voxels_per_grain_min_median_max": [int(cnt.min()), float(np.median(cnt)), int(cnt.max())],
        "n_boundary_voxels": int(bnd.sum()),
        "n_interior_voxels": int((~bnd).sum()),
        "near_axis_r_mm": NEAR_AXIS_R,
        "n_near_axis": int((r < NEAR_AXIS_R).sum()),
        "n_far_axis": int((r >= NEAR_AXIS_R).sum()),
        "r_min_max_mm": [float(r.min()), float(r.max())],
        "boundary_definition": "a different-grain voxel (full sample) within 1.01 x side of the "
        "centroid, as the adjacency of make_neworient.py",
    }
    (OUT / "region.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
