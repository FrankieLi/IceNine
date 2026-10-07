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
import platform
from typing import Dict, Tuple

import numpy as np
import pytest
import scipy
import torch
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

# Environment the golden values were recorded in (the one `uv sync --extra dev` builds from
# uv.lock with .python-version). The orientation differs at ~4e-7 deg in other environments
# (see MIGRATION_HISTORY "Follow-ups (2026-10-06)"), so the golden tests refuse to run elsewhere.
# Python is compared on major.minor; numpy, torch and scipy on their release (local build tag
# such as "+cpu" ignored).
GOLDEN_ENV: Dict[str, str] = {
    "python": "3.9",
    "numpy": "2.0.2",
    "torch": "2.8.0",
    "scipy": "1.13.1",
}


def _current_env() -> Dict[str, str]:
    return {
        "python": ".".join(platform.python_version_tuple()[:2]),
        "numpy": np.__version__.split("+")[0],
        "torch": torch.__version__.split("+")[0],
        "scipy": scipy.__version__.split("+")[0],
    }


def _require_data() -> Path:
    data_dir = EXAMPLE / "ScatteringData_Python"
    if not data_dir.exists() or len(list(data_dir.glob("*.d*"))) < 360:
        pytest.skip("ThreeVoxels Python-simulated data not available")
    return data_dir


@pytest.fixture
def golden_env() -> None:
    """Fail (never skip, never loosen the tolerance) if the golden values were recorded with
    different library versions than the ones running."""
    _require_data()  # skip (as _build does) when the data is absent; only then compare versions
    now = _current_env()
    if now != GOLDEN_ENV:
        fmt = lambda d: ", ".join(f"{k} {v}" for k, v in d.items())  # noqa: E731
        pytest.fail(
            f"golden values recorded with {fmt(GOLDEN_ENV)}; running {fmt(now)}. "
            "Run `uv sync --extra dev` in icenine_py (locked environment) and use `uv run`.",
            pytrace=False,
        )


def _build(min_sin_eta: float = 0.0):
    data_dir = _require_data()
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


def test_reconstruct_voxel_unchanged_by_refactor(golden_env):
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


def test_recorder_hook_leaves_reconstruct_voxel_bit_identical(golden_env):
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


def _run_golden_problem(configure=None):
    """reconstruct_voxel on the GOLDEN problem with the recorder attached; `configure(rec, R_true)`
    may set knobs. Returns (result, events, R_true)."""
    rec, voxel, R_true = _build()
    events = []
    rec.recorder = lambda name, data: events.append((name, data))
    if configure is not None:
        configure(rec, R_true)
    res = rec.reconstruct_voxel(
        _get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7)
    )
    return res, events, R_true


def _assert_golden(res, R_true):
    np.testing.assert_allclose(_rotvec_deg(res.orientation, R_true), GOLDEN[0], atol=1e-9)
    assert res.cost == pytest.approx(GOLDEN[1], abs=1e-12)


def test_rank_key_with_post_mc_costs_reproduces_golden(golden_env):
    def configure(rec, R_true):
        rec.rank_key = lambda level, cands: np.array([c.cost for c in cands])

    res, _, R_true = _run_golden_problem(configure)
    _assert_golden(res, R_true)


def test_extra_candidates_returning_nothing_reproduces_golden(golden_env):
    def configure(rec, R_true):
        rec.extra_candidates = lambda level, cands: []

    res, _, R_true = _run_golden_problem(configure)
    _assert_golden(res, R_true)


def test_keep_fraction_changes_n_keep():
    """The GOLDEN problem has one candidate per level, so three extra candidates are added (which
    also exercises extra_candidates): 4 candidates, n_keep = int(4 * f)."""

    def with_fraction(f):
        def configure(rec, R_true):
            rec.keep_fraction = f
            rec.extra_candidates = lambda level, cands: [
                SearchCandidate(
                    orientation=Rotation.from_rotvec(np.radians(v)).as_matrix() @ R_true, cost=1.0
                )
                for v in ([0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5])
            ]

        return configure

    n_keep = {}
    for f in (0.25, 0.5):
        _, events, _ = _run_golden_problem(with_fraction(f))
        quick = [d for n, d in events if n == "quick_mc"][0]
        assert len(quick["cost"]) == 4
        n_keep[f] = quick["n_keep"]
    assert n_keep == {0.25: 1, 0.5: 2}
