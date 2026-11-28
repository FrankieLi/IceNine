"""fine1: the sample ``fine7`` (250 grains on the 98,280-voxel generation-7 grid, generator seed 7).

  uv run python scripts/autoresearch/fine1/fine1_make_sample.py

``held_make_sample.py``'s steps with its helpers imported unchanged (``held_sample_lib``,
``scripts/phase_d/grains.py``, ``scripts/autoresearch/ar1/misses.py``): 250 Voronoi cells, grains
as the connected pieces of the cells, one orientation per grain drawn from ``grains.py`` and
redrawn if under 5 degrees from a neighbour grain, gradient levels 0, 0.2, 0.5, 1, 2 and 4
degrees drawn per grain, the linear rotation-vector model. One addition, the gradient cap
(``fine1_lib.capped_gradient_orientations``): a grain whose largest neighbour-pair
misorientation, after its gradient is drawn, is 0.45 degree or more is scaled down to 0.45
degree. A new generator seed is a new sample (as ``held3`` to ``held6``), not an extra seed for a
comparison. Needs ``fine1_make_grid.py`` first.

Writes ``<cache>/sample/gradient_truth.mic`` (not tracked; its sha256 is) and the tracked
``benchmarks/autoresearch/fine1/sample/{sample.json, grainmap.npz, caps.json}``.
"""

import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fine1_lib as F  # noqa: E402
import fine1_paths as P  # noqa: E402

ICE = P.ICE
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar3"))
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar1"))
sys.path.insert(0, str(ICE / "scripts" / "phase_d"))
import grains as GR  # noqa: E402
import held_sample_lib as H  # noqa: E402
import misses as M  # noqa: E402


def cubic_quats() -> np.ndarray:
    """The cubic symmetry quaternions of the simulation config, the form ``grains.py`` takes
    (as ``held_make_sample._cubic_quats``)."""
    from icenine.config_file import ConfigFile
    from icenine.experiment_setup import XDMExperimentSetup
    from icenine.orientation_search import get_symmetry_quaternions

    cfg = P.EX / "ConfigFiles" / "Example2.Simulation.config"
    return get_symmetry_quaternions(
        XDMExperimentSetup(ConfigFile.from_file(str(cfg))).get_sample_symmetry()
    )


def main() -> None:
    from icenine.gradient_sample import spread_deg
    from icenine.mic_file import MicFile
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    seed = P.SEED
    P.SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    P.SAMPLE_JSON.parent.mkdir(parents=True, exist_ok=True)
    fine = MicFile.read(str(P.GRID))
    fv, fside = M.voxel_vertices_xy(P.GRID)
    fc = fv.mean(axis=1)
    nf = len(fine.voxels)
    assert nf == P.N_FINE, nf
    side = float(fside.max())
    seeds = H.draw_cell_seeds(nf, P.N_CELLS, seed)
    cell = H.voronoi_cells(fc, seeds)
    gmap = H.connected_grains(fc, cell, 1.2 * side)
    gids = sorted(int(g) for g in np.unique(gmap))
    size = {g: int((gmap == g).sum()) for g in gids}
    np.savez_compressed(P.GRAINMAP, grain=gmap.astype(np.int32))
    # orientations, one per grain, separated from the adjacent grains
    sym_q = cubic_quats()
    adjacency = GR.grain_adjacency(fc, gmap, 1.2 * side)
    attempts = np.zeros(len(gids), dtype=np.int64)
    q = np.array([GR.draw_orientation(seed, g, 0, sym_q) for g in gids])
    for _round in range(200):
        a = np.array([x for x, _ in adjacency])
        b = np.array([y for _, y in adjacency])
        ang = np.array(
            [
                GR.misorientation_matrix_deg(q[i : i + 1], q[j : j + 1], sym_q)[0, 0]
                for i, j in zip(a, b)
            ]
        )
        bad = H.grains_to_redraw(a, b, ang, P.MIN_SEP_DEG)
        if not bad:
            break
        for g in bad:
            attempts[g] += 1
            q[g] = GR.draw_orientation(seed, g, int(attempts[g]), sym_q)
    else:
        raise RuntimeError("orientations did not separate")
    eul = GR.quat_to_euler(q)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "fine7_cells.mic"
        GR.write_mic_with_euler(str(P.GRID), str(tmp), gmap, np.round(eul, GR.EULER_DECIMALS))
        back0 = MicFile.read(str(tmp))
    R_fine = np.array([v.orientation for v in back0.voxels], dtype=np.float64)
    R_old = {g: R_fine[int(np.nonzero(gmap == g)[0][0])].copy() for g in gids}
    for g in gids:
        assert np.array_equal(R_fine[gmap == g], np.broadcast_to(R_old[g], (size[g], 3, 3)))
    sizes = np.array([size[g] for g in gids])
    check: Dict[str, Any] = dict(
        n_voxels=nf,
        n_voronoi_cells=P.N_CELLS,
        n_grains_connected=len(gids),
        generator_seed=seed,
        grain_size_min_median_max=[int(sizes.min()), float(np.median(sizes)), int(sizes.max())],
        mean_voxels_per_grain=float(nf / len(gids)),
        n_adjacent_pairs=len(adjacency),
        min_sep_deg_rule=P.MIN_SEP_DEG,
        n_redrawn=int((attempts > 0).sum()),
    )
    print(json.dumps(check))
    # gradients: step 7's code and settings, with the cap
    grng, lv = H.draw_levels(seed, len(gids), P.LEVELS)
    levels = dict(zip(gids, lv))
    pairs = F.edge_pairs(fc, side)
    R_new, capped = F.capped_gradient_orientations(fc, gmap, R_old, levels, grng, pairs)
    for v in range(nf):
        fine.voxels[v].orientation = np.asarray(R_new[v], dtype=np.float64).copy()
    fine.write(str(P.TRUTH))
    back = MicFile.read(str(P.TRUTH))
    Rb = np.array([v.orientation for v in back.voxels], dtype=np.float64)
    sym = proper_rotations(create_cubic_symmetry(3.6))
    spread = spread_deg(Rb, gmap, R_old, sym)
    sha = P.sha256(P.TRUTH)
    out = dict(
        seed=seed,
        levels_deg=list(P.LEVELS),
        n_voxels=nf,
        n_grains=len(gids),
        truth_sha256=sha,
        gradient_cap_deg=F.GRADIENT_CAP_DEG,
        n_capped=len(capped),
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
        n_by_level={str(L): int(sum(1 for g in gids if levels[g] == L)) for L in P.LEVELS},
        check=check,
    )
    P.SAMPLE_JSON.write_text(json.dumps(out) + "\n")
    P.CAPS_JSON.write_text(
        json.dumps(
            dict(
                cap_deg=F.GRADIENT_CAP_DEG,
                n_edge_pairs=int(len(pairs)),
                n_grains=len(gids),
                n_capped=len(capped),
                capped=capped,
            ),
            indent=1,
        )
        + "\n"
    )
    print("voxels", nf, "grains", len(gids), "capped", len(capped), "truth sha256", sha)
    print("grains by level", out["n_by_level"])


if __name__ == "__main__":
    main()
