"""AR3 step 7: the sample ``big`` (32 grains on the gradient sample's grid, with gradients).

  AR3_SAMPLE=big uv run python scripts/autoresearch/ar3/big_make_sample.py

Each voxel of ``rand_500grains_1mm_neworient_s0_grid.mic`` (24,570 voxels) takes the orientation of
the voxel of ``random_cubic_32grains.mic`` (6,144 voxels, one size coarser) that holds its centre.
Checks (written to ``big/gradient/sample/big_check.json``): every fine centre lies in exactly one
coarse triangle; the 32 coarse orientations are distinct. The grain map is the connected voxels of
one orientation (voxels whose centres are within 1.2 voxel sides are neighbours). Gradients: the
library's ``gradient_orientations``, generator seed 0, one level per grain drawn uniformly from
LEVELS (0, 0.2, 0.5, 1, 2, 4 degrees), as ``make_gradient_sample.py`` does for its five levels.
Writes ``big/gradient/sample/``: ``gradient_truth.mic``, ``gradient_grid.mic`` (orientations
zeroed), ``grainmap.npy``, ``sample.json``, ``big_check.json``.
"""

import hashlib
import json
import sys

import ar3_common as C
import gradient_common as G
import numpy as np

assert G.BIG, "run with AR3_SAMPLE=big"
SEED = 0


def main() -> None:
    sys.path.insert(0, str(C.ICE / "scripts" / "autoresearch" / "ar1"))
    import misses as M
    from icenine.gradient_sample import gradient_orientations, spread_deg
    from icenine.mic_file import MicFile
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    sim = C.EX / "SimInput"
    coarse_path = sim / "random_cubic_32grains.mic"
    fine_path = sim / "rand_500grains_1mm_neworient_s0_grid.mic"
    G.SAMPLE.mkdir(parents=True, exist_ok=True)
    coarse = MicFile.read(str(coarse_path))
    fine = MicFile.read(str(fine_path))
    cv, _cs = M.voxel_vertices_xy(coarse_path)  # (nc, 3, 2)
    fv, fside = M.voxel_vertices_xy(fine_path)
    fc = fv.mean(axis=1)
    nf, nc = len(fine.voxels), len(coarse.voxels)
    # which coarse triangle holds each fine centre (barycentric test, with a small tolerance)
    tree = cKDTree(cv.mean(axis=1))
    _d, near = tree.query(fc, k=8)
    a, b, c = cv[:, 0], cv[:, 1], cv[:, 2]
    n_in = np.zeros(nf, dtype=int)
    owner = np.full(nf, -1)
    tol = 1e-9
    for k in range(near.shape[1]):
        j = near[:, k]
        v0, v1, v2 = b[j] - a[j], c[j] - a[j], fc - a[j]
        den = v0[:, 0] * v1[:, 1] - v1[:, 0] * v0[:, 1]
        u = (v2[:, 0] * v1[:, 1] - v1[:, 0] * v2[:, 1]) / den
        w = (v0[:, 0] * v2[:, 1] - v2[:, 0] * v0[:, 1]) / den
        inside = (u >= -tol) & (w >= -tol) & (u + w <= 1 + tol)
        n_in += inside
        owner = np.where(inside & (owner < 0), j, owner)
    Rc = np.array([v.orientation for v in coarse.voxels], dtype=np.float64)
    keys = np.round(Rc.reshape(nc, 9), 6)
    uniq, inv = np.unique(keys, axis=0, return_inverse=True)
    R_fine = Rc[owner]
    fkeys = inv.reshape(-1)[owner]
    # grains: connected voxels of one orientation
    pairs = cKDTree(fc).query_pairs(1.2 * float(fside.max()), output_type="ndarray")
    same = fkeys[pairs[:, 0]] == fkeys[pairs[:, 1]]
    p = pairs[same]
    adj = coo_matrix((np.ones(len(p)), (p[:, 0], p[:, 1])), shape=(nf, nf))
    ncomp, comp = connected_components(adj, directed=False)
    gmap = comp.astype(np.int64)
    gids = sorted(int(g) for g in np.unique(gmap))
    size = {g: int((gmap == g).sum()) for g in gids}
    np.save(G.GRAINMAP_BIG, gmap)
    R_old = {g: R_fine[int(np.nonzero(gmap == g)[0][0])].copy() for g in gids}
    for g in gids:  # one stored orientation per grain
        assert np.array_equal(R_fine[gmap == g], np.broadcast_to(R_old[g], (size[g], 3, 3)))
    check = dict(
        n_fine=nf,
        n_coarse=nc,
        fine_with_exactly_one_coarse=int((n_in == 1).sum()),
        fine_with_none=int((n_in == 0).sum()),
        fine_with_two_or_more=int((n_in >= 2).sum()),
        distinct_orientations=int(len(uniq)),
        coarse_voxels_per_orientation=np.bincount(inv.reshape(-1)).tolist(),
        n_grains_connected=len(gids),
        grain_sizes=[size[g] for g in gids],
    )
    print(json.dumps(check))
    assert check["fine_with_exactly_one_coarse"] == nf and check["distinct_orientations"] == 32
    # gradients
    rng = np.random.default_rng(SEED)
    lv_idx = rng.integers(0, len(G.LEVELS), size=len(gids))
    levels = {g: float(G.LEVELS[i]) for g, i in zip(gids, lv_idx)}
    R_new = gradient_orientations(fc, gmap, R_old, levels, rng)
    C.write_mic(fine, {v: R_new[v] for v in range(nf)}, G.TRUTH_MIC)
    back = MicFile.read(str(G.TRUTH_MIC))
    Rb = np.array([v.orientation for v in back.voxels], dtype=np.float64)
    sym = _cubic()
    spread = spread_deg(Rb, gmap, R_old, sym)
    C.write_mic(back, {v: np.zeros((3, 3)) for v in range(nf)}, G.GRID_MIC)
    sha = hashlib.sha256(G.TRUTH_MIC.read_bytes()).hexdigest()
    out = dict(
        seed=SEED,
        levels_deg=list(G.LEVELS),
        n_voxels=nf,
        n_grains=len(gids),
        truth_sha256=sha,
        grains=[
            dict(
                id=g,
                level_deg=levels[g],
                measured_spread_deg=spread[g],
                n_voxels=size[g],
                R_centre=R_old[g].tolist(),
            )
            for g in gids
        ],
        n_by_level={str(L): int(sum(1 for g in gids if levels[g] == L)) for L in G.LEVELS},
    )
    (G.SAMPLE / "sample.json").write_text(json.dumps(out) + "\n")
    err = max(abs(spread[g] - levels[g]) for g in gids if size[g] > 1)
    check["max_abs_spread_minus_level_deg"] = float(err)
    check["n_by_level"] = out["n_by_level"]
    (G.SAMPLE / "big_check.json").write_text(json.dumps(check, indent=1) + "\n")
    print("grains by level", out["n_by_level"], "max |spread - level| (deg)", err)


def _cubic() -> np.ndarray:
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    return proper_rotations(create_cubic_symmetry(3.6))


if __name__ == "__main__":
    main()
