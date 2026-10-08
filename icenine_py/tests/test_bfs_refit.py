"""BFS refit pass, CMA neighbour budget, wider-sigma retry and per-voxel provenance.

The scripted tests replace the reconstructor with a fake whose local fit succeeds only within a
capture angle of the voxel's true orientation (1 deg, 1 deg + sigma0 for the retry), on a three
voxel fan where every voxel neighbours every other. They check the BFS bookkeeping exactly. The
ThreeVoxels tests (skipped without the Python-simulated data) check that the default path is
bit-identical to numbers recorded from the code BEFORE these options existed, and that the refit
pass runs on the real cost function.
"""

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

import icenine.reconstructor as reconstructor_module
from icenine.config_file import ConfigFile
from icenine.mic_file import MicFile, ReconstructionState, Voxel
from icenine.orientation_search import SearchCandidate, SearchParameters
from icenine.reconstructor import BFSReconstruction

# tests/ is on sys.path under pytest
from test_findoptimal_refactor import EXAMPLE, _build, golden_env  # noqa: F401

# ---------------------------------------------------------------------------
# Scripted world
# ---------------------------------------------------------------------------

POSITIONS = [(0.0, 0.0), (1.5, 0.0), (0.75, 1.3)]  # side 1, radius 2: all mutual neighbours
GOOD, BAD = 0.95, 0.3  # hit ratio of a fit that found / missed the truth
FITTED, REFIT = ReconstructionState.FITTED, ReconstructionState.REFIT


def _rz(deg: float) -> np.ndarray:
    return Rotation.from_rotvec(np.radians([0.0, 0.0, deg])).as_matrix().astype(np.float32)


def _angle(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.degrees(np.linalg.norm(Rotation.from_matrix(a @ b.T).as_rotvec())))


def _info(hit: float) -> Any:
    return SimpleNamespace(
        pixel_overlap=int(hit * 1000),
        pixel_on_detector=1000,
        peak_overlap=int(hit * 100),
        peak_on_detector=100,
        cost=1.0 - hit,
    )


class FakeReconstructor:
    """Stands in for AdaptiveVoxelReconstructor: a full search returns the truth; a local fit
    returns the truth if the start is within the capture angle, else the start unchanged."""

    truth: List[np.ndarray] = []
    local_calls: List[Dict[str, Any]] = []
    full_calls: List[int] = []

    def __init__(self, setup: Any) -> None:
        self.last_local_optimization_evals = 0
        self._counts = (0, 0, 0)

    @property
    def last_eval_counts(self) -> Tuple[int, int, int]:
        return self._counts

    @staticmethod
    def _idx(vertices: Any) -> int:
        x, y = float(vertices[0][0]), float(vertices[0][1])
        return min(range(3), key=lambda i: abs(POSITIONS[i][0] - x) + abs(POSITIONS[i][1] - y))

    def reconstruct_voxel(self, voxel_vertices: Any, phase_index: int, rng: Any) -> SearchCandidate:
        i = self._idx(voxel_vertices)
        FakeReconstructor.full_calls.append(i)
        self._counts = (900, 100, 1000)
        return SearchCandidate(self.truth[i].astype(np.float64), 0.05, _info(GOOD))

    def evaluate_overlap(self, R: np.ndarray, vertices: Any, phase: int) -> Any:
        return _info(GOOD if _angle(R, self.truth[self._idx(vertices)]) < 0.3 else BAD)

    def local_optimization(
        self,
        voxel_vertices: Any,
        phase_index: int,
        initial_orientation: np.ndarray,
        rng: Any = None,
        cma_sigma0_deg: Optional[float] = None,
        cma_max_evals: Optional[int] = None,
    ) -> SearchCandidate:
        i = self._idx(voxel_vertices)
        FakeReconstructor.local_calls.append(
            dict(idx=i, sigma0=cma_sigma0_deg, max_evals=cma_max_evals)
        )
        capture = 1.0 + (cma_sigma0_deg or 0.0)
        hit_start = _angle(initial_orientation, self.truth[i])
        if hit_start < capture:
            R, hit = self.truth[i].astype(np.float64), GOOD
        else:
            R, hit = np.asarray(initial_orientation, dtype=np.float64), BAD
        self.last_local_optimization_evals = cma_max_evals or 100
        return SearchCandidate(R, 1.0 - hit, _info(hit))


def _run(
    deltas: Tuple[float, float],
    monkeypatch: pytest.MonkeyPatch,
    gate: float = 0.0,
    **params: Any,
) -> Tuple[BFSReconstruction, MicFile]:
    """BFS over the fan: voxel 0 is the seed (truth 0 deg); voxels 1, 2 are `deltas` degrees
    away from it about z."""
    truth = [_rz(0.0), _rz(deltas[0]), _rz(deltas[1])]
    FakeReconstructor.truth = truth
    FakeReconstructor.local_calls = []
    FakeReconstructor.full_calls = []
    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", FakeReconstructor)
    voxels = [
        Voxel(position=np.array([x, y, 0.0]), orientation=np.eye(3), side_length=1.0)
        for x, y in POSITIONS
    ]
    mic = MicFile(voxels)
    setup = SimpleNamespace(
        sample=SimpleNamespace(get_mic=lambda: mic),
        config=SimpleNamespace(min_acceleration_threshold=0.8, partial_result_acceptance_conf=gate),
        search_params=SearchParameters(**params),
    )
    bfs = BFSReconstruction(setup)  # type: ignore[arg-type]
    order = SimpleNamespace(shuffle=lambda x: None)  # seed order 0, 1, 2
    done = bfs.reconstruct_sample(rng=order)  # type: ignore[arg-type]
    assert sorted(done) == [0, 1, 2]
    return bfs, mic


def _states(mic: MicFile) -> List[int]:
    return [v.reconstruction_id for v in mic.voxels]


# ---------------------------------------------------------------------------
# Default path: no retry, no refit, mc call signature unchanged
# ---------------------------------------------------------------------------


def test_default_leaves_rejected_voxels_refit(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run((3.0, 3.5), monkeypatch)
    assert _states(mic) == [FITTED, REFIT, REFIT]
    s = bfs.stats
    assert s["n_seeds"] == 1 and s["n_neighbor_fits"] == 2 and s["n_unresolved"] == 2
    assert s["n_refit_attempted"] == 0 and s["n_retry_attempted"] == 0
    # mc mode: no CMA override reaches local_optimization
    assert all(
        c["sigma0"] is None and c["max_evals"] is None for c in FakeReconstructor.local_calls
    )
    assert [bfs.records[i].source for i in range(3)] == ["seed", "unresolved", "unresolved"]
    assert s["n_evals_seed"] == 1000 and s["n_evals_neighbor"] == 200 and s["n_evals_refit"] == 0
    assert bfs.records[1].n_evals == 100 and bfs.records[0].n_evals == 1000


def test_provenance_neighbor_accepted(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run((0.5, -0.5), monkeypatch)
    assert _states(mic) == [FITTED, FITTED, FITTED]
    assert [bfs.records[i].source for i in range(3)] == ["seed", "neighbor", "neighbor"]
    assert bfs.stats["n_neighbor_accepted_first"] == 2 and bfs.stats["n_unresolved"] == 0
    assert all(r.wall_s >= 0 for r in bfs.records.values())


# ---------------------------------------------------------------------------
# Refit pass
# ---------------------------------------------------------------------------


def test_refit_full_search_and_expansion(monkeypatch: pytest.MonkeyPatch) -> None:
    """Gate above the local confidence: voxel 1 gets a full search, its expansion then fixes
    voxel 2 (0.5 deg from voxel 1) with a local fit and no search of its own."""
    bfs, mic = _run((3.0, 3.5), monkeypatch, bfs_refit=True, bfs_refit_conf=0.99)
    assert _states(mic) == [FITTED, FITTED, FITTED]
    s = bfs.stats
    assert s["n_refit_candidates"] == 2 and s["n_refit_attempted"] == 1
    assert s["n_refit_full"] == 1 and s["n_refit_local"] == 0
    assert s["n_refit_resolved"] == 2 and s["n_unresolved"] == 0
    assert FakeReconstructor.full_calls == [0, 1]  # seed, then the refit of voxel 1 only
    r1, r2 = bfs.records[1], bfs.records[2]
    assert r1.source == "refit" and r1.refit_tried and r1.refit_mode == "full"
    assert r2.source == "refit" and not r2.refit_tried
    assert r1.n_evals == 100 + 100 + 1000  # first local fit, refit local fit, full search
    assert s["n_evals_refit"] == 100 + 1000 + 100  # refit local + full for 1, expansion fit of 2
    assert _angle(mic.voxels[1].orientation, _rz(3.0)) < 0.01
    assert _angle(mic.voxels[2].orientation, _rz(3.5)) < 0.01


def test_refit_gate_from_config_when_not_set(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run((3.0, 3.5), monkeypatch, gate=0.99, bfs_refit=True)
    assert bfs.records[1].refit_mode == "full" and _states(mic) == [FITTED] * 3


def test_refit_local_start_keeps_failed_fit_state(monkeypatch: pytest.MonkeyPatch) -> None:
    """Gate 0: the refit keeps the local fit as the centre (skip-discrete); it misses the truth
    (3 deg away), fails the hit-ratio test and the voxels stay REFIT with their old fit."""
    bfs, mic = _run((3.0, 3.5), monkeypatch, bfs_refit=True, bfs_refit_conf=0.0)
    assert _states(mic) == [FITTED, REFIT, REFIT]
    s = bfs.stats
    assert s["n_refit_attempted"] == 2 and s["n_refit_local"] == 2 and s["n_refit_full"] == 0
    assert s["n_refit_resolved"] == 0 and s["n_unresolved"] == 2
    assert FakeReconstructor.full_calls == [0]
    assert all(bfs.records[i].source == "unresolved" for i in (1, 2))
    assert _angle(mic.voxels[1].orientation, _rz(0.0)) < 0.01  # stored fit unchanged


def test_refit_off_by_default_runs_nothing_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    on, _ = _run((3.0, 3.5), monkeypatch, bfs_refit=True, bfs_refit_conf=0.99)
    off, mic = _run((3.0, 3.5), monkeypatch, bfs_refit=False, bfs_refit_conf=0.99)
    assert on.stats["n_unresolved"] == 0 and off.stats["n_unresolved"] == 2
    assert _states(mic) == [FITTED, REFIT, REFIT]


# ---------------------------------------------------------------------------
# CMA neighbour budget and wider-sigma retry
# ---------------------------------------------------------------------------


def test_cma_neighbor_budget_and_retry(monkeypatch: pytest.MonkeyPatch) -> None:
    """Voxel 1 is 1.8 deg away: the first fit (capture 1 deg) fails, the retry (sigma0 1.5, capture
    2.5 deg) succeeds. Voxel 2 is 4 deg away: both fail and it stays REFIT."""
    bfs, mic = _run((1.8, 4.0), monkeypatch, local_optimizer="cma", cma_neighbor_max_evals=250)
    assert _states(mic) == [FITTED, FITTED, REFIT]
    calls = FakeReconstructor.local_calls
    assert all(c["max_evals"] == 250 for c in calls)  # neighbours: the neighbour budget
    assert sum(c["sigma0"] == 1.5 for c in calls) == 2 and len(calls) == 4
    s = bfs.stats
    assert s["n_retry_attempted"] == 2 and s["n_retry_accepted"] == 1
    assert s["n_neighbor_accepted_first"] == 0
    r1, r2 = bfs.records[1], bfs.records[2]
    assert r1.source == "neighbor_retry" and r1.retried and r1.retry_accepted and r1.local_rejected
    assert r2.source == "unresolved" and r2.retried and not r2.retry_accepted
    assert r1.n_evals == 500 and r2.n_evals == 500
    assert s["n_evals_neighbor"] == 1000
    assert _angle(mic.voxels[1].orientation, _rz(1.8)) < 0.01


def test_retry_disabled_with_zero_sigma(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run((1.8, 4.0), monkeypatch, local_optimizer="cma", cma_retry_sigma0_deg=0.0)
    assert _states(mic) == [FITTED, REFIT, REFIT]
    assert bfs.stats["n_retry_attempted"] == 0
    assert all(c["sigma0"] is None for c in FakeReconstructor.local_calls)


def test_failed_retry_keeps_first_fit(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run((4.0, 5.0), monkeypatch, local_optimizer="cma")
    assert _states(mic) == [FITTED, REFIT, REFIT]
    assert _angle(mic.voxels[1].orientation, _rz(0.0)) < 0.01  # the inherited start, as before


def test_cma_retry_then_refit(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run(
        (4.0, 4.5), monkeypatch, local_optimizer="cma", bfs_refit=True, bfs_refit_conf=0.99
    )
    assert _states(mic) == [FITTED] * 3
    assert bfs.records[1].retried and bfs.records[1].source == "refit"
    # the refit fit (local, budget 250) carries no wider sigma; the retry does not run again
    assert bfs.stats["n_retry_attempted"] == 2


# ---------------------------------------------------------------------------
# Parameters and config keys
# ---------------------------------------------------------------------------


def test_defaults_and_validation() -> None:
    sp = SearchParameters()
    assert sp.bfs_refit is False and sp.bfs_refit_conf is None
    assert sp.cma_neighbor_max_evals == 250 and sp.cma_retry_sigma0_deg == 1.5
    SearchParameters(cma_retry_sigma0_deg=0.0)
    for bad in (
        dict(cma_neighbor_max_evals=1),
        dict(cma_retry_sigma0_deg=-0.1),
        dict(bfs_refit_conf=1.5),
        dict(bfs_refit_conf=-0.1),
    ):
        with pytest.raises(ValueError):
            SearchParameters(**bad)  # type: ignore[arg-type]


_CFG = EXAMPLE / "ConfigFiles" / "ReconstructQ8.config"


def _config_with(tmp_path: Path, extra: str) -> ConfigFile:
    if not _CFG.exists():
        pytest.skip("ThreeVoxels config not available")
    p = tmp_path / "x.config"
    p.write_text(_CFG.read_text() + "\n" + extra + "\n")
    return ConfigFile.from_file(str(p))


def test_config_absent_keys_are_defaults(tmp_path: Path) -> None:
    sp = SearchParameters.from_config(_config_with(tmp_path, ""))
    assert sp.bfs_refit is False and sp.bfs_refit_conf is None
    assert sp.cma_neighbor_max_evals == 250 and sp.cma_retry_sigma0_deg == 1.5


def test_config_keys(tmp_path: Path) -> None:
    cfg = _config_with(
        tmp_path, "BFSRefit 1\nBFSRefitConf 0.6\nCMANeighborMaxEvals 120\nCMARetrySigma0 0"
    )
    sp = SearchParameters.from_config(cfg)
    assert sp.bfs_refit is True and sp.bfs_refit_conf == pytest.approx(0.6)
    assert sp.cma_neighbor_max_evals == 120 and sp.cma_retry_sigma0_deg == 0.0


@pytest.mark.parametrize(
    "line",
    ["BFSRefit 2", "BFSRefitConf 1.5", "CMANeighborMaxEvals 1", "CMARetrySigma0 -1", "BFSRefit"],
)
def test_config_bad_values(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


# ---------------------------------------------------------------------------
# ThreeVoxels (real cost function): default parity and a real refit pass
# ---------------------------------------------------------------------------

# Recorded with the code BEFORE these options (BFS seed 1, max_local_resolution 0,
# max_mc_steps 150, the fz set of test_cma_optimizer's BFS test): voxel 0 is rejected as a seed
# (hit 0.797 < 0.8); voxels 1 and 2 are seeds of their own (the voxels are far apart).
GOLDEN_STATES = [REFIT, FITTED, FITTED]
GOLDEN_HIT = [0.7970479704797048, 0.9774436090225563, 1.0]


def _three_voxel_setup() -> Any:
    rec, _, _ = _build()
    setup = rec.setup
    mic = setup.sample.get_mic()
    truth = [np.asarray(v.orientation, dtype=np.float64).copy() for v in mic.voxels]
    setup.fz_orientations = np.stack(
        [
            Rotation.from_rotvec(np.radians([1.5, -1.0, 0.5])).as_matrix() @ truth[0],
            truth[1],
            truth[2],
        ]
    )
    return setup


def test_default_bfs_matches_recorded_old_result(golden_env: None) -> None:  # noqa: F811
    setup = _three_voxel_setup()
    bfs = BFSReconstruction(setup)
    done = bfs.reconstruct_sample(rng=np.random.default_rng(1))
    mic = setup.sample.get_mic()
    assert done == [0, 1, 2] and _states(mic) == GOLDEN_STATES
    assert [float(v.overlap_ratio) for v in mic.voxels] == GOLDEN_HIT
    assert bfs.stats["n_refit_attempted"] == 0 and bfs.stats["n_unresolved"] == 1
    assert [bfs.records[i].source for i in range(3)] == ["unresolved", "seed", "seed"]
    assert all(bfs.records[i].n_evals > 0 for i in range(3))


@pytest.mark.parametrize("gate", [0.0, 0.99])
def test_refit_runs_on_real_data(golden_env: None, gate: float) -> None:  # noqa: F811
    """The rejected seed (hit 0.797) is revisited once: a local start (gate 0) or a full search
    (gate 0.99). Its outcome is whatever the cost function gives; the bookkeeping is checked."""
    setup = _three_voxel_setup()
    setup.search_params.bfs_refit = True
    setup.search_params.bfs_refit_conf = gate
    bfs = BFSReconstruction(setup)
    bfs.reconstruct_sample(rng=np.random.default_rng(1))
    s, r0 = bfs.stats, bfs.records[0]
    assert s["n_refit_candidates"] == 1 and s["n_refit_attempted"] == 1 and r0.refit_tried
    assert r0.refit_mode == ("local" if gate == 0.0 else "full")
    assert s["n_refit_resolved"] + s["n_unresolved"] == 1
    assert r0.n_evals > 0
    assert s["n_evals_refit"] > 0
