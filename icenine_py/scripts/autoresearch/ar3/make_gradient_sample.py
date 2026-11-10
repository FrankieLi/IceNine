"""AR3 step 3: write the gradient truth map, its grid-only copy and the record JSON.

  uv run python scripts/autoresearch/ar3/make_gradient_sample.py

Reads the layer's truth (``rand_500grains_1mm_neworient_s0.mic``) and its voxel-to-grain map.
Levels: one draw per grain, uniformly from the five values in grain-id order, generator seed 0 (the
same generator then draws the gradient directions; see ``icenine.gradient_sample``). Writes to
``benchmarks/autoresearch/ar3/gradient/sample/``: ``gradient_truth.mic``, ``gradient_grid.mic``
(orientation columns zeroed), ``sample.json`` (seed, level and measured spread of each grain,
sha256 of the truth map). The truth map is read back from the written file for the spread, so the
record states what the images were rendered from (the .mic stores Euler angles to 1e-4 degrees).
"""

import hashlib
import json
import sys
from pathlib import Path

import ar3_common as C
import numpy as np

SAMPLE = C.OUT / "gradient" / "sample"
SEED = 0


def main() -> None:
    sys.path.insert(0, str(C.ICE / "scripts" / "autoresearch" / "ar1"))
    import misses as M
    from icenine.gradient_sample import LEVELS_DEG, gradient_orientations, spread_deg
    from icenine.mic_file import MicFile

    SAMPLE.mkdir(parents=True, exist_ok=True)
    truth = MicFile.read(str(M.MIC))
    gmap = np.load(M.GRAINMAP)
    n = len(truth.voxels)
    assert len(gmap) == n
    vxy, _side = M.voxel_vertices_xy(M.MIC)
    centres = vxy.mean(axis=1)
    gids = sorted(int(g) for g in np.unique(gmap))
    R_old = {
        g: np.asarray(truth.voxels[int(np.nonzero(gmap == g)[0][0])].orientation, np.float64)
        for g in gids
    }
    rng = np.random.default_rng(SEED)
    lv_idx = rng.integers(0, len(LEVELS_DEG), size=len(gids))
    levels = {g: float(LEVELS_DEG[i]) for g, i in zip(gids, lv_idx)}
    R_new = gradient_orientations(centres, gmap, R_old, levels, rng)
    C.write_mic(truth, {v: R_new[v] for v in range(n)}, SAMPLE / "gradient_truth.mic")

    # read back what was written (the rendering uses this file)
    back = MicFile.read(str(SAMPLE / "gradient_truth.mic"))
    Rb = np.array([v.orientation for v in back.voxels], dtype=np.float64)
    sym = _cubic()
    spread = spread_deg(Rb, gmap, R_old, sym)
    grid = MicFile.read(str(M.MIC.parent / "rand_500grains_1mm_neworient_s0_grid.mic"))
    gxy = np.array([v.position[:2] for v in grid.voxels])
    assert np.abs(gxy - np.array([v.position[:2] for v in back.voxels])).max() < 1e-6
    C.write_mic(
        back,
        {v: np.zeros((3, 3)) for v in range(n)},
        SAMPLE / "gradient_grid.mic",
    )
    sha = hashlib.sha256((SAMPLE / "gradient_truth.mic").read_bytes()).hexdigest()
    size = np.array([int((gmap == g).sum()) for g in gids])
    out = dict(
        seed=SEED,
        levels_deg=list(LEVELS_DEG),
        n_voxels=n,
        n_grains=len(gids),
        truth_sha256=sha,
        grains=[
            dict(
                id=g,
                level_deg=levels[g],
                measured_spread_deg=spread[g],
                n_voxels=int(size[i]),
                R_centre=R_old[g].tolist(),
            )
            for i, g in enumerate(gids)
        ],
        n_by_level={str(L): int(sum(1 for g in gids if levels[g] == L)) for L in LEVELS_DEG},
    )
    (SAMPLE / "sample.json").write_text(json.dumps(out) + "\n")
    err = max(abs(spread[g] - levels[g]) for g in gids if size[gids.index(g)] > 1)
    print("grains by level", out["n_by_level"], "max |spread - level| (deg)", err)


def _cubic() -> np.ndarray:
    from icenine.ray_indexing import proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    return proper_rotations(create_cubic_symmetry(3.6))


if __name__ == "__main__":
    main()
