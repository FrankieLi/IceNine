"""tex1: the sample ``tex24`` (250 grains on fine7's 98,280-voxel grid, generator seed 24).

  uv run python scripts/autoresearch/tex1/tex1_make_sample.py

``fine1_make_sample.py``'s grain steps with its helpers imported unchanged (``held_sample_lib``,
``scripts/phase_d/grains.py``, ``scripts/autoresearch/ar1/misses.py``): 250 Voronoi cells, grains
as the connected pieces of the cells, one orientation ``R_g`` per grain drawn from ``grains.py``
and redrawn if under 5 degrees from a neighbour grain. The texture replaces fine1's gradient
levels (``tex1_lib``): 20% of the grains get voxel-to-voxel noise (rotation vector uniform in a
1-degree ball, so any two voxels differ by at most 2 degrees), 50% get a smooth gradient field
(largest edge-pair misorientation ``t_g`` uniform on [0.1, 0.5) degree); the two sets are drawn
independently. A new generator seed is a new sample, not an extra seed for a comparison.
Needs ``fine7_grid.mic`` in the cache (a copy of fine1's, sha256 checked against
``fine1/frozen.json``).

Writes ``<cache>/sample/tex24_truth.mic`` (not tracked; its sha256 is) and the tracked
``benchmarks/autoresearch/tex1/sample/{sample.json, grainmap.npz}``.
"""

import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tex1_lib as F  # noqa: E402
import tex1_paths as P  # noqa: E402

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
    from icenine.mic_file import MicFile

    seed = P.SEED
    frozen_grid = json.loads(P.FINE1_FROZEN.read_text())["sha256"]["grid"]
    assert P.sha256(P.GRID) == frozen_grid, "fine7_grid.mic differs from fine1's frozen grid"
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
    # the texture
    noise_g, grad_g = F.draw_classes(gids, seed, P.NOISE_FRACTION, P.GRADIENT_FRACTION)
    pairs = F.edge_pairs(fc, side)
    R_new, R_smooth, grecs = F.textured_orientations(fc, gmap, R_old, noise_g, grad_g, pairs, seed)
    for v in range(nf):
        fine.voxels[v].orientation = np.asarray(R_new[v], dtype=np.float64).copy()
    fine.write(str(P.TRUTH))
    sha = P.sha256(P.TRUTH)
    nset, gset = set(noise_g), set(grad_g)
    classes = {g: F.class_of(g in nset, g in gset) for g in gids}
    out = dict(
        seed=seed,
        n_voxels=nf,
        n_grains=len(gids),
        truth_sha256=sha,
        noise_ball_deg=F.NOISE_BALL_DEG,
        gradient_cap_deg=F.GRADIENT_CAP_DEG,
        gradient_min_deg=F.GRADIENT_MIN_DEG,
        curvature=F.CURVATURE,
        n_noise_grains=len(noise_g),
        n_gradient_grains=len(grad_g),
        n_by_class={
            F.CLASS_NAMES[c]: int(sum(1 for g in gids if classes[g] == c)) for c in range(4)
        },
        grains=[
            dict(
                id=g,
                cls=classes[g],
                n_voxels=size[g],
                R_centre=R_old[g].tolist(),
            )
            for g in gids
        ],
        gradient_grains=grecs,
        check=check,
    )
    P.SAMPLE_JSON.write_text(json.dumps(out) + "\n")
    print("voxels", nf, "grains", len(gids), "truth sha256", sha)
    print("grains by class", out["n_by_class"])


if __name__ == "__main__":
    main()
