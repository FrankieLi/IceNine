"""AR3 step 11: the held-out sample ``held`` (250 grains on the 24,570-voxel gradient grid).

  AR3_SAMPLE=held uv run python scripts/autoresearch/ar3/held_make_sample.py

``start_make_mid.py`` (the builder of ``mid``) with three differences, each a constant here:
250 Voronoi cells (``mid``: 120), generator seed 1 (``mid``: 0) for the centres, the orientations
and the gradient levels, and the tag ``held``. Same 24,570 voxel centres, same six gradient levels
(0, 0.2, 0.5, 1, 2 and 4 degrees), same rule that redraws an orientation closer than 5 degrees to a
neighbour grain's. The seed is new because ``held`` is a different sample for a carry-over test
(every other sample has been read for a choice); it is not a repeat to firm up a comparison.

Writes ``held/gradient/sample/``: ``gradient_truth.mic``, ``gradient_grid.mic`` (orientations
zeroed), ``grainmap.npy``, ``sample.json``, ``held_check.json`` (grain sizes, grains per level, the
measured spread per level; the lit pixels and spots are in ``held_stack.json``).
"""

import hashlib
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict

import ar3_common as C
import gradient_common as G
import held_sample_lib as H
import numpy as np

assert C.SAMPLE_KIND in ("held", "held2"), "run with AR3_SAMPLE=held or held2"
# AR5 step 12: ``held2`` is this builder with generator seed 2 (``held``: seed 1); no other change
SEED = {"held": 1, "held2": 2}[C.SAMPLE_KIND]
N_CELLS = 250
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

    fine_path = C.EX / "SimInput" / "rand_500grains_1mm_neworient_s0_grid.mic"
    G.SAMPLE.mkdir(parents=True, exist_ok=True)
    fine = MicFile.read(str(fine_path))
    fv, fside = M.voxel_vertices_xy(fine_path)
    fc = fv.mean(axis=1)
    nf = len(fine.voxels)
    seeds = H.draw_cell_seeds(nf, N_CELLS, SEED)
    cell = H.voronoi_cells(fc, seeds)
    gmap = H.connected_grains(fc, cell, 1.2 * float(fside.max()))
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
        bad = H.grains_to_redraw(a, b, ang, MIN_SEP_DEG)
        if not bad:
            break
        for g in bad:
            attempts[g] += 1
            q[g] = GR.draw_orientation(SEED, g, int(attempts[g]), sym_q)
    else:
        raise RuntimeError("orientations did not separate")
    eul = GR.quat_to_euler(q)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "held_cells.mic"
        GR.write_mic_with_euler(str(fine_path), str(tmp), gmap, np.round(eul, GR.EULER_DECIMALS))
        back0 = MicFile.read(str(tmp))
    R_fine = np.array([v.orientation for v in back0.voxels], dtype=np.float64)
    R_old = {g: R_fine[int(np.nonzero(gmap == g)[0][0])].copy() for g in gids}
    for g in gids:
        assert np.array_equal(R_fine[gmap == g], np.broadcast_to(R_old[g], (size[g], 3, 3)))
    sizes = np.array([size[g] for g in gids])
    check: Dict[str, Any] = dict(
        n_fine=nf,
        n_voronoi_cells=N_CELLS,
        n_grains_connected=len(gids),
        generator_seed=SEED,
        grain_size_min_median_max=[int(sizes.min()), float(np.median(sizes)), int(sizes.max())],
        mean_voxels_per_grain=float(nf / len(gids)),
        n_adjacent_pairs=len(adjacency),
        min_sep_deg_rule=MIN_SEP_DEG,
        n_redrawn=int((attempts > 0).sum()),
        grain_sizes=[size[g] for g in gids],
    )
    print(json.dumps(check))
    # gradients (step 7's code and settings)
    grng, lv = H.draw_levels(SEED, len(gids), G.LEVELS)
    levels = dict(zip(gids, lv))
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
    by_level: Dict[str, Any] = {}
    for L in G.LEVELS:
        sp = [spread[g] for g in gids if levels[g] == L and size[g] > 1]
        by_level[str(L)] = dict(
            n_grains=int(sum(1 for g in gids if levels[g] == L)),
            n_voxels=int(sum(size[g] for g in gids if levels[g] == L)),
            measured_spread_deg_min_median_max=(
                [float(np.min(sp)), float(np.median(sp)), float(np.max(sp))] if sp else None
            ),
        )
    check["by_level"] = by_level
    (G.SAMPLE / "held_check.json").write_text(json.dumps(check, indent=1) + "\n")
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
