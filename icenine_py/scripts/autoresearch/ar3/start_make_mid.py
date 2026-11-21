"""AR3 step 9: the sample ``mid`` (about 120 grains on the gradient sample's 24,570-voxel grid).

  AR3_SAMPLE=mid uv run python scripts/autoresearch/ar3/start_make_mid.py

Step 7's generator (``big_make_sample.py``) read a 32-grain file that is not available at 120
grains, so the grain map is new; everything after it is step 7's code with its stored settings:
the library's ``gradient_orientations``, generator seed 0, one level per grain drawn uniformly
from LEVELS (0, 0.2, 0.5, 1, 2, 4 degrees), the written truth map read back for the spread.

The grain map: 120 voxel centres drawn with ``default_rng(0)`` are the seeds of a Voronoi
tessellation of the voxel centres (each voxel joins the nearest seed). Grains are the connected
pieces of one cell (centres within 1.2 voxel sides are neighbours), as in step 7. Orientations:
one per grain, uniform on SO(3), reduced to the cubic fundamental zone, drawn with the phase D
helper ``grains.draw_orientation`` (seed 0, grain index, attempt) and redrawn (attempt + 1) for
the later grain of any adjacent pair closer than MIN_SEP_DEG (5 degrees; set here, not tuned) so
that a grain at 4 degrees does not look like its neighbour. The Euler angles are written to the .mic
precision and read back, as ``make_neworient.py`` does.

Writes ``mid/gradient/sample/``: ``gradient_truth.mic``, ``gradient_grid.mic`` (orientations
zeroed), ``grainmap.npy``, ``sample.json``, ``mid_check.json``.
"""

import hashlib
import json
import sys
import tempfile
from pathlib import Path

import ar3_common as C
import gradient_common as G
import numpy as np

assert C.SAMPLE_KIND == "mid", "run with AR3_SAMPLE=mid"
SEED = 0
N_SEEDS = 120
MIN_SEP_DEG = 5.0


def main() -> None:
    sys.path.insert(0, str(C.ICE / "scripts" / "autoresearch" / "ar1"))
    sys.path.insert(0, str(C.ICE / "scripts" / "phase_d"))
    import grains as GR
    import misses as M
    from icenine.gradient_sample import gradient_orientations, spread_deg
    from icenine.mic_file import MicFile
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    fine_path = C.EX / "SimInput" / "rand_500grains_1mm_neworient_s0_grid.mic"
    G.SAMPLE.mkdir(parents=True, exist_ok=True)
    fine = MicFile.read(str(fine_path))
    fv, fside = M.voxel_vertices_xy(fine_path)
    fc = fv.mean(axis=1)
    nf = len(fine.voxels)
    rng = np.random.default_rng(SEED)
    seeds = np.sort(rng.choice(nf, size=N_SEEDS, replace=False))
    cell = cKDTree(fc[seeds]).query(fc, k=1)[1]
    pairs = cKDTree(fc).query_pairs(1.2 * float(fside.max()), output_type="ndarray")
    same = cell[pairs[:, 0]] == cell[pairs[:, 1]]
    p = pairs[same]
    adj = coo_matrix((np.ones(len(p)), (p[:, 0], p[:, 1])), shape=(nf, nf))
    ncomp, comp = connected_components(adj, directed=False)
    gmap = comp.astype(np.int64)
    gids = sorted(int(g) for g in np.unique(gmap))
    size = {g: int((gmap == g).sum()) for g in gids}
    np.save(G.GRAINMAP_BIG, gmap)
    # orientations, one per grain, separated from the adjacent grains
    sym_q = _cubic_quats()
    adjacency = GR.grain_adjacency(fc, gmap, 1.2 * float(fside.max()))
    attempts = np.zeros(len(gids), dtype=np.int64)
    q = np.array([GR.draw_orientation(SEED, g, 0, sym_q) for g in gids])
    for _round in range(200):
        a = np.array([x for x, _ in adjacency])
        b = np.array([y for _, y in adjacency])
        ang = np.array(
            [
                GR.misorientation_matrix_deg(q[i : i + 1], q[j : j + 1], sym_q)[0, 0]
                for i, j in zip(a, b)
            ]
        )
        bad = sorted({int(max(i, j)) for i, j, d in zip(a, b, ang) if d < MIN_SEP_DEG})
        if not bad:
            break
        for g in bad:
            attempts[g] += 1
            q[g] = GR.draw_orientation(SEED, g, int(attempts[g]), sym_q)
    else:
        raise RuntimeError("orientations did not separate")
    eul = GR.quat_to_euler(q)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "mid_cells.mic"
        GR.write_mic_with_euler(str(fine_path), str(tmp), gmap, np.round(eul, GR.EULER_DECIMALS))
        back0 = MicFile.read(str(tmp))
    R_fine = np.array([v.orientation for v in back0.voxels], dtype=np.float64)
    R_old = {g: R_fine[int(np.nonzero(gmap == g)[0][0])].copy() for g in gids}
    for g in gids:
        assert np.array_equal(R_fine[gmap == g], np.broadcast_to(R_old[g], (size[g], 3, 3)))
    sizes = np.array([size[g] for g in gids])
    check = dict(
        n_fine=nf,
        n_voronoi_cells=N_SEEDS,
        n_grains_connected=len(gids),
        grain_size_min_median_max=[int(sizes.min()), float(np.median(sizes)), int(sizes.max())],
        mean_voxels_per_grain=float(nf / len(gids)),
        n_adjacent_pairs=len(adjacency),
        min_sep_deg_rule=MIN_SEP_DEG,
        n_redrawn=int((attempts > 0).sum()),
        grain_sizes=[size[g] for g in gids],
    )
    print(json.dumps(check))
    # gradients (step 7's code and settings)
    grng = np.random.default_rng(SEED)
    lv_idx = grng.integers(0, len(G.LEVELS), size=len(gids))
    levels = {g: float(G.LEVELS[i]) for g, i in zip(gids, lv_idx)}
    R_new = gradient_orientations(fc, gmap, R_old, levels, grng)
    C.write_mic(fine, {v: R_new[v] for v in range(nf)}, G.TRUTH_MIC)
    back = MicFile.read(str(G.TRUTH_MIC))
    Rb = np.array([v.orientation for v in back.voxels], dtype=np.float64)
    sym = proper_rotations(create_cubic_symmetry(3.6))
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
    (G.SAMPLE / "mid_check.json").write_text(json.dumps(check, indent=1) + "\n")
    print("grains by level", out["n_by_level"], "max |spread - level| (deg)", err)


def _cubic_quats() -> np.ndarray:
    """The cubic symmetry quaternions of the simulation config, the form ``grains.py`` takes
    (as ``make_neworient.py`` reads them)."""
    from icenine.config_file import ConfigFile
    from icenine.experiment_setup import XDMExperimentSetup
    from icenine.orientation_search import get_symmetry_quaternions

    cfg = C.EX / "ConfigFiles" / "Example2.Simulation.config"
    return get_symmetry_quaternions(
        XDMExperimentSetup(ConfigFile.from_file(str(cfg))).get_sample_symmetry()
    )


if __name__ == "__main__":
    main()
