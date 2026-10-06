"""AdaptiveVoxelReconstructor.refine_from_candidates (the final stage factored out of
reconstruct_voxel) leaves reconstruct_voxel unchanged, and runs standalone.

The golden numbers (GOLDEN) were recorded with reconstruct_voxel BEFORE the refactor, on a small
deterministic problem (ThreeVoxels voxel 0, a 3-orientation FZ set, one level, 150 MC steps,
seed 7).
Skipped without the ThreeVoxels Python-simulated data.
"""

import math
import os
from pathlib import Path
from typing import Tuple

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from icenine.config_file import ConfigFile
from icenine.experimental_data import ExperimentalData
from icenine.mic_file import MicFile
from icenine.orientation_search import SearchCandidate
from icenine.reconstructor import (
    AdaptiveVoxelReconstructor,
    _get_voxel_vertices,
    setup_reconstruction,
)

EXAMPLE = Path(__file__).parent.parent.parent / "Examples" / "Example2.ThreeVoxels"

# (final orientation as rotvec of R_final R_true^T in degrees, cost) recorded before the refactor
GOLDEN: Tuple[np.ndarray, float] = (
    np.array([-0.12518150458946703, 0.16879997485875445, 0.02470329742498792]),
    0.8180327868852459,
)


def _build(min_sin_eta: float = 0.0):
    data_dir = EXAMPLE / "ScatteringData_Python"
    if not data_dir.exists() or len(list(data_dir.glob("*.d*"))) < 360:
        pytest.skip("ThreeVoxels Python-simulated data not available")
    cwd = os.getcwd()
    os.chdir(EXAMPLE)
    try:
        config = ConfigFile.from_file("ConfigFiles/ReconstructQ8.config")
        config.out_file_basename = "3Grains.sim"
        exp_data = ExperimentalData.from_image_directory(
            directory=str(data_dir),
            basename="3Grains.sim",
            ext="d",
            serial_length=5,
            n_omega=180,
            n_detectors=2,
            num_rows=2048,
            num_cols=2048,
        )
        mic = MicFile.read("SimInput/three_voxels.mic")
        voxel = mic.voxels[0]
        R_true = np.asarray(voxel.orientation, dtype=np.float64)
        fz = np.stack(
            [
                Rotation.from_rotvec(np.radians([1.5, -1.0, 0.5])).as_matrix() @ R_true,
                Rotation.random(random_state=1).as_matrix(),
                Rotation.random(random_state=2).as_matrix(),
            ]
        )
        setup = setup_reconstruction(config, exp_data=exp_data, fz_orientations=fz)
        setup.search_params.max_local_resolution = 0
        setup.search_params.max_mc_steps = 150
        rec = AdaptiveVoxelReconstructor(setup, min_sin_eta=min_sin_eta)
    finally:
        os.chdir(cwd)
    return rec, voxel, R_true


def _rotvec_deg(R: np.ndarray, R_true: np.ndarray) -> np.ndarray:
    return np.degrees(Rotation.from_matrix(R @ R_true.T).as_rotvec())


def test_reconstruct_voxel_unchanged_by_refactor():
    rec, voxel, R_true = _build()
    res = rec.reconstruct_voxel(
        _get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7)
    )
    err = _rotvec_deg(res.orientation, R_true)
    np.testing.assert_allclose(err, GOLDEN[0], atol=1e-9)
    assert res.cost == pytest.approx(GOLDEN[1], abs=1e-12)
    assert rec.last_find_optimal["n_candidates"] >= 1
    assert len(rec.last_level_best) == 1


def test_refine_from_candidates_standalone_converges():
    rec, voxel, R_true = _build()
    start = Rotation.from_rotvec(np.radians([0.2, -0.1, 0.15])).as_matrix() @ R_true
    res = rec.refine_from_candidates(
        [SearchCandidate(orientation=start, cost=1.0)],
        _get_voxel_vertices(voxel),
        voxel.phase,
        rng=np.random.default_rng(3),
    )
    assert np.linalg.norm(_rotvec_deg(res.orientation, R_true)) < 0.2
    g, loc, tot = rec.last_eval_counts
    assert g == 0 and loc > 0 and tot == loc
    assert math.isfinite(res.cost) and res.overlap_info is not None


def test_refine_from_candidates_empty_returns_identity():
    rec, voxel, _ = _build()
    res = rec.refine_from_candidates([], _get_voxel_vertices(voxel), voxel.phase)
    assert res.cost == 1.0 and np.array_equal(res.orientation, np.eye(3))


def test_recorder_hook_leaves_reconstruct_voxel_bit_identical():
    """With a recorder attached (and the knobs at their defaults) reconstruct_voxel returns exactly
    the GOLDEN numbers, and the recorder sees every stage."""
    rec, voxel, R_true = _build()
    events = []
    rec.recorder = lambda name, data: events.append((name, data))
    res = rec.reconstruct_voxel(
        _get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7)
    )
    err = _rotvec_deg(res.orientation, R_true)
    np.testing.assert_allclose(err, GOLDEN[0], atol=1e-9)
    assert res.cost == pytest.approx(GOLDEN[1], abs=1e-12)
    names = [n for n, _ in events]
    assert names[0] == "discrete" and names[1] == "quick_mc"
    assert names.count("find_candidate") == rec.last_find_optimal["n_evaluated"]
    assert names[-2:] == ["variance", "final"]
    quick = dict(events[1][1])
    assert sorted(quick["perm"].tolist()) == list(range(len(quick["perm"])))
    assert np.all(np.diff(quick["cost"]) >= 0)  # sorted best first
    np.testing.assert_allclose(events[-1][1]["R"], res.orientation)
