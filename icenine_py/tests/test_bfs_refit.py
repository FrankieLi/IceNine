"""BFS revisit of REFIT voxels, restart pass, CMA neighbour budget / retry, per-voxel provenance.

Three kinds of test:
- A scripted 14 x 14 grid (tests/bfs_scripted.py) whose golden (tests/data/bfs_grid_golden.npz)
  was recorded from the commit BEFORE these options existed: with every option off, the new code
  reproduces the done order, states, orientations, hit ratios and the whole fit sequence (an
  rng-drawing fake makes the final generator state depend on every fit), in mc mode and in cma
  mode with the retry off.
- Small scripted voxel groups whose fake reconstructor succeeds only within a capture angle of
  the truth (1 deg, 1 deg + sigma0 for the retry), to check the revisit bookkeeping exactly.
- ThreeVoxels (skipped without the Python-simulated data): default parity with recorded numbers.
"""

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

import icenine.mic_file as mic_module
import icenine.reconstructor as reconstructor_module
from icenine.config_file import ConfigFile
from icenine.mic_file import MicFile, ReconstructionState, Voxel
from icenine.orientation_search import SearchCandidate, SearchParameters
from icenine.reconstructor import BFSReconstruction

# tests/ is on sys.path under pytest
import bfs_scripted
from test_findoptimal_refactor import EXAMPLE, _build, golden_env  # noqa: F401

FITTED, REFIT = ReconstructionState.FITTED, ReconstructionState.REFIT

# ---------------------------------------------------------------------------
# Default-path parity on the scripted multi-grain grid (golden from the parent commit)
# ---------------------------------------------------------------------------

GOLDEN = Path(__file__).parent / "data" / "bfs_grid_golden.npz"


@pytest.mark.parametrize(
    "label,params",
    [("mc", {}), ("cma", {"local_optimizer": "cma", "cma_retry_sigma0_deg": 0.0})],
)
def test_default_options_reproduce_parent_commit(label: str, params: Dict[str, Any]) -> None:
    gold = np.load(GOLDEN)
    out = bfs_scripted.run_grid(
        reconstructor_module, mic_module, SearchCandidate, SearchParameters, params
    )
    assert np.array_equal(out["done"], gold[f"{label}_done"])
    assert np.array_equal(out["states"], gold[f"{label}_states"])
    assert np.array_equal(out["log_kind"], gold[f"{label}_log_kind"])
    assert np.array_equal(out["log_idx"], gold[f"{label}_log_idx"])
    np.testing.assert_allclose(out["log_val"], gold[f"{label}_log_val"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(out["R"], gold[f"{label}_R"], rtol=0, atol=1e-6)  # float32 store
    np.testing.assert_allclose(out["hit"], gold[f"{label}_hit"], rtol=0, atol=1e-12)
    assert float(out["rng_after"]) == float(gold[f"{label}_rng_after"])
    # the grid exercises the interesting case: many voxels are left REFIT
    assert int((out["states"] == REFIT).sum()) > 20
    s = out["bfs"].stats
    assert s["n_revisit_attempted"] == 0 and s["n_restart_attempted"] == 0


def test_revisit_improves_the_scripted_grid() -> None:
    """With the revisit on, the grid ends with fewer REFIT and fewer wrong (> 1 deg) voxels.
    The counts are of the scripted fake (acceptance and error against its truth), a mechanism
    check, not a statement about real data."""
    base = bfs_scripted.run_grid(
        reconstructor_module, mic_module, SearchCandidate, SearchParameters, {}
    )
    rev = bfs_scripted.run_grid(
        reconstructor_module,
        mic_module,
        SearchCandidate,
        SearchParameters,
        {"bfs_revisit_refit": True},
    )
    s = rev["bfs"].stats
    assert s["n_revisit_attempted"] > 0 and s["n_revisit_accepted"] > 0
    assert s["n_revisit_accepted"] + s["n_revisit_rejected"] == s["n_revisit_attempted"]
    assert int((rev["states"] == REFIT).sum()) < int((base["states"] == REFIT).sum())
    assert bfs_scripted.n_wrong(rev) < bfs_scripted.n_wrong(base)
    assert sorted(rev["done"]) == sorted(base["done"])  # same voxels processed, none twice
    assert len(set(rev["done"].tolist())) == len(rev["done"])
    assert all(r.n_revisits <= 3 for r in rev["bfs"].records.values())


# ---------------------------------------------------------------------------
# Small scripted groups
# ---------------------------------------------------------------------------

# voxel 0 is the first seed. Voxels 1 and 2 touch 0 and each other; voxel 3 touches 1 and 2 only.
POSITIONS = [(0.0, 0.0), (1.5, 0.0), (0.75, 1.3), (2.25, 1.3)]
GOOD = 0.95


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
    """A full search returns the truth. A local fit returns the truth if the start is within the
    capture angle, else the start with a hit ratio that falls with the distance to the truth."""

    truth: List[np.ndarray] = []
    local_calls: List[Dict[str, Any]] = []

    def __init__(self, setup: Any) -> None:
        self.last_local_optimization_evals = 0
        self._counts = (0, 0, 0)

    @property
    def last_eval_counts(self) -> Tuple[int, int, int]:
        return self._counts

    @staticmethod
    def _idx(vertices: Any) -> int:
        x, y = float(vertices[0][0]), float(vertices[0][1])
        return min(range(4), key=lambda i: abs(POSITIONS[i][0] - x) + abs(POSITIONS[i][1] - y))

    def reconstruct_voxel(self, voxel_vertices: Any, phase_index: int, rng: Any) -> SearchCandidate:
        self._counts = (900, 100, 1000)
        return SearchCandidate(
            self.truth[self._idx(voxel_vertices)].astype(np.float64), 0.05, _info(GOOD)
        )

    def evaluate_overlap(self, R: np.ndarray, vertices: Any, phase: int) -> Any:
        return _info(GOOD)

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
        a = _angle(initial_orientation, self.truth[i])
        if a < 1.0 + (cma_sigma0_deg or 0.0):
            R, hit = self.truth[i].astype(np.float64), GOOD
        else:
            R, hit = np.asarray(initial_orientation, dtype=np.float64), max(0.05, 0.6 - 0.1 * a)
        self.last_local_optimization_evals = cma_max_evals or 100
        return SearchCandidate(R, 1.0 - hit, _info(hit))


def _run(
    deltas: Tuple[float, ...],
    monkeypatch: pytest.MonkeyPatch,
    **params: Any,
) -> Tuple[BFSReconstruction, MicFile]:
    """BFS over the first len(deltas)+1 positions, seeds in index order: voxel 0 has truth 0 deg,
    voxel k has truth deltas[k-1] degrees (about z)."""
    n = len(deltas) + 1
    FakeReconstructor.truth = [_rz(0.0)] + [_rz(d) for d in deltas] + [_rz(0.0)] * (4 - n)
    FakeReconstructor.local_calls = []
    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", FakeReconstructor)
    voxels = [
        Voxel(position=np.array([x, y, 0.0]), orientation=np.eye(3), side_length=1.0)
        for x, y in POSITIONS[:n]
    ]
    mic = MicFile(voxels)
    setup = SimpleNamespace(
        sample=SimpleNamespace(get_mic=lambda: mic),
        config=SimpleNamespace(min_acceleration_threshold=0.8),
        search_params=SearchParameters(**params),
    )
    bfs = BFSReconstruction(setup)  # type: ignore[arg-type]
    done = bfs.reconstruct_sample(rng=SimpleNamespace(shuffle=lambda x: None))  # type: ignore
    assert sorted(done) == list(range(n))
    return bfs, mic


def _states(mic: MicFile) -> List[int]:
    return [v.reconstruction_id for v in mic.voxels]


def test_default_leaves_rejected_voxels_refit(monkeypatch: pytest.MonkeyPatch) -> None:
    """Voxels 1, 2 (3 and 3.5 deg away) fail from seed 0; voxel 3 is a seed of its own and its
    expansion does not touch them: they stay REFIT, nothing is revisited."""
    bfs, mic = _run((3.0, 3.5, 3.2), monkeypatch)
    assert _states(mic) == [FITTED, REFIT, REFIT, FITTED]
    s = bfs.stats
    assert s["n_seeds"] == 2 and s["n_neighbor_fits"] == 2 and s["n_unresolved"] == 2
    assert s["n_revisit_attempted"] == 0 and len(FakeReconstructor.local_calls) == 2
    assert all(
        c["sigma0"] is None and c["max_evals"] is None for c in FakeReconstructor.local_calls
    )
    assert [bfs.records[i].source for i in range(4)] == ["seed", "unresolved", "unresolved", "seed"]
    assert s["n_evals_seed"] == 2000 and s["n_evals_neighbor"] == 200
    assert bfs.records[1].n_evals == 100 and bfs.records[0].n_evals == 1000


def test_provenance_neighbor_accepted(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run((0.5, -0.5), monkeypatch)
    assert _states(mic) == [FITTED, FITTED, FITTED]
    assert [bfs.records[i].source for i in range(3)] == ["seed", "neighbor", "neighbor"]
    assert bfs.stats["n_neighbor_accepted_first"] == 2 and bfs.stats["n_unresolved"] == 0


# ---------------------------------------------------------------------------
# Revisit
# ---------------------------------------------------------------------------


def test_revisit_from_the_true_grain_is_accepted(monkeypatch: pytest.MonkeyPatch) -> None:
    """Voxel 3 (seed of the grain that voxels 1, 2 belong to) reaches the REFIT voxels, which
    inherit its orientation, fit within the capture angle and are accepted; each is tried once
    even though both voxel 3 and the freshly fitted voxel 1 reach voxel 2."""
    bfs, mic = _run((3.0, 3.5, 3.2), monkeypatch, bfs_revisit_refit=True)
    assert _states(mic) == [FITTED] * 4
    s = bfs.stats
    assert s["n_revisit_attempted"] == 2 and s["n_revisit_accepted"] == 2
    assert s["n_revisit_rejected"] == 0 and s["n_unresolved"] == 0
    assert len(FakeReconstructor.local_calls) == 4  # 2 neighbour fits + 2 revisits, no repeats
    for i in (1, 2):
        r = bfs.records[i]
        assert r.source == "revisit" and r.n_revisits == 1 and r.local_rejected
        assert r.n_evals == 200  # the failed first fit and the accepted revisit
    assert _angle(mic.voxels[1].orientation, _rz(3.0)) < 0.01
    assert _angle(mic.voxels[2].orientation, _rz(3.5)) < 0.01
    assert s["n_evals_revisit"] == 200 and s["n_evals_neighbor"] == 200


def test_revisit_rejected_keeps_the_better_old_fit(monkeypatch: pytest.MonkeyPatch) -> None:
    """Voxel 3 is 8 deg from the seed's grain: inheriting from it gives a worse fit of voxels 1
    and 2 than the old one (3 and 3.5 deg off), so the snapshot is restored (Push keep-best)."""
    bfs, mic = _run((3.0, 3.5, 8.0), monkeypatch, bfs_revisit_refit=True)
    assert _states(mic) == [FITTED, REFIT, REFIT, FITTED]
    s = bfs.stats
    assert s["n_revisit_attempted"] == 2 and s["n_revisit_rejected"] == 2
    for i, d in ((1, 3.0), (2, 3.5)):
        v = mic.voxels[i]
        assert _angle(v.orientation, _rz(0.0)) < 0.01  # still the first fit (the seed's)
        assert v.confidence == pytest.approx(max(0.05, 0.6 - 0.1 * d), abs=0.011)
        assert bfs.records[i].n_revisits == 1 and bfs.records[i].source == "unresolved"
        assert bfs.records[i].hit_ratio == pytest.approx(v.overlap_ratio)


def test_rejected_again_with_a_better_new_fit_replaces_the_old(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Voxel 3 is 5.5 deg from the seed's grain: still outside the capture angle of voxels 1 and
    2, but closer to them (2.5 / 2 deg) than the seed's orientation was (3 / 3.5 deg), so the
    rejected revisit has the higher confidence and replaces the stored fit."""
    bfs, mic = _run((3.0, 3.5, 5.5), monkeypatch, bfs_revisit_refit=True)
    assert _states(mic) == [FITTED, REFIT, REFIT, FITTED]
    for i, d in ((1, 3.0), (2, 3.5)):
        v = mic.voxels[i]
        assert _angle(v.orientation, _rz(5.5)) < 0.01  # the revisit's (inherited) fit
        assert v.confidence > 0.6 - 0.1 * d + 0.01  # higher than the old fit's
    assert bfs.stats["n_revisit_rejected"] == 2


def test_revisit_cap_is_respected(monkeypatch: pytest.MonkeyPatch) -> None:
    """With the restart pass on top, FITTED voxel 3 reaches voxels 1 and 2 a second time (it did
    in its own expansion). The cap decides whether that second visit happens."""
    capped, _ = _run(
        (3.0, 3.5, 8.0),
        monkeypatch,
        bfs_revisit_refit=True,
        bfs_restart_pass=True,
        bfs_revisit_max=1,
    )
    assert capped.stats["n_revisit_capped"] >= 2
    assert all(capped.records[i].n_revisits == 1 for i in (1, 2))
    assert capped.stats["n_restart_attempted"] == 0
    loose, _ = _run(
        (3.0, 3.5, 8.0),
        monkeypatch,
        bfs_revisit_refit=True,
        bfs_restart_pass=True,
        bfs_revisit_max=3,
    )
    assert loose.stats["n_restart_attempted"] > 0
    assert all(1 <= loose.records[i].n_revisits <= 3 for i in (1, 2))


def test_revisit_off_changes_nothing(monkeypatch: pytest.MonkeyPatch) -> None:
    on, _ = _run((3.0, 3.5, 3.2), monkeypatch, bfs_revisit_refit=True)
    off, mic = _run((3.0, 3.5, 3.2), monkeypatch, bfs_revisit_refit=False)
    assert on.stats["n_unresolved"] == 0 and off.stats["n_unresolved"] == 2
    assert _states(mic) == [FITTED, REFIT, REFIT, FITTED]


def test_restart_pass_pushes_fitted_borders_onto_refit(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run((3.0, 3.5, 3.2), monkeypatch, bfs_restart_pass=True)
    assert _states(mic) == [FITTED] * 4
    s = bfs.stats
    assert s["n_restart_candidates"] == 2 and s["n_restart_accepted"] == 2
    assert s["n_revisit_attempted"] == 0
    assert [bfs.records[i].source for i in (1, 2)] == ["restart", "restart"]


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
    assert s["n_retry_attempted_neighbor"] == 2 and s["n_retry_accepted_neighbor"] == 1
    assert s["n_retry_attempted_revisit"] == 0 and s["n_neighbor_accepted_first"] == 0
    r1, r2 = bfs.records[1], bfs.records[2]
    assert r1.source == "neighbor_retry" and r1.retried and r1.retry_accepted and r1.local_rejected
    assert r2.source == "unresolved" and r2.retried and not r2.retry_accepted
    assert r1.n_evals == 500 and r2.n_evals == 500 and s["n_evals_neighbor"] == 1000
    assert _angle(mic.voxels[1].orientation, _rz(1.8)) < 0.01


def test_retry_disabled_with_zero_sigma(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run((1.8, 4.0), monkeypatch, local_optimizer="cma", cma_retry_sigma0_deg=0.0)
    assert _states(mic) == [FITTED, REFIT, REFIT]
    assert bfs.stats["n_retry_attempted_neighbor"] == 0
    assert all(c["sigma0"] is None for c in FakeReconstructor.local_calls)


def test_failed_retry_keeps_first_fit(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run((4.0, 5.0), monkeypatch, local_optimizer="cma")
    assert _states(mic) == [FITTED, REFIT, REFIT]
    assert _angle(mic.voxels[1].orientation, _rz(0.0)) < 0.01  # the inherited start, as before


def test_cma_revisit_uses_the_retry_and_counts_it_separately(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Voxels 1 and 2 are 4 deg from seed 0 (plain fit and retry both fail) but 1.5 deg from
    voxel 3: the plain revisit fit (capture 1 deg) fails, the retry (2.5 deg) wins."""
    bfs, mic = _run((4.0, 4.0, 5.5), monkeypatch, local_optimizer="cma", bfs_revisit_refit=True)
    s = bfs.stats
    assert s["n_retry_attempted_neighbor"] == 2 and s["n_retry_accepted_neighbor"] == 0
    assert s["n_retry_attempted_revisit"] == 2 and s["n_retry_accepted_revisit"] == 2
    assert _states(mic) == [FITTED] * 4
    assert all(bfs.records[i].source == "revisit" and bfs.records[i].retry_accepted for i in (1, 2))


# ---------------------------------------------------------------------------
# Parameters and config keys
# ---------------------------------------------------------------------------


def test_defaults_and_validation() -> None:
    sp = SearchParameters()
    assert sp.bfs_revisit_refit is False and sp.bfs_restart_pass is False
    assert sp.bfs_revisit_max == 3
    assert sp.cma_neighbor_max_evals == 250 and sp.cma_retry_sigma0_deg == 1.5
    SearchParameters(cma_retry_sigma0_deg=0.0)
    for bad in (
        dict(cma_neighbor_max_evals=1),
        dict(cma_retry_sigma0_deg=-0.1),
        dict(bfs_revisit_max=0),
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
    assert sp.bfs_revisit_refit is False and sp.bfs_restart_pass is False
    assert sp.bfs_revisit_max == 3
    assert sp.cma_neighbor_max_evals == 250 and sp.cma_retry_sigma0_deg == 1.5


def test_config_keys(tmp_path: Path) -> None:
    cfg = _config_with(
        tmp_path,
        "BFSRevisitRefit 1\nBFSRevisitMax 2\nBFSRestartPass 1\n"
        "CMANeighborMaxEvals 120\nCMARetrySigma0 0",
    )
    sp = SearchParameters.from_config(cfg)
    assert sp.bfs_revisit_refit is True and sp.bfs_restart_pass is True
    assert sp.bfs_revisit_max == 2
    assert sp.cma_neighbor_max_evals == 120 and sp.cma_retry_sigma0_deg == 0.0


@pytest.mark.parametrize(
    "line",
    [
        "BFSRevisitRefit 2",
        "BFSRestartPass 3",
        "BFSRevisitMax 0",
        "CMANeighborMaxEvals 1",
        "CMARetrySigma0 -1",
        "BFSRevisitRefit",
    ],
)
def test_config_bad_values(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


# ---------------------------------------------------------------------------
# ThreeVoxels (real cost function): default parity
# ---------------------------------------------------------------------------

# Recorded with the code BEFORE these options (BFS seed 1, max_local_resolution 0,
# max_mc_steps 150, the fz set of test_cma_optimizer's BFS test): voxel 0 is rejected as a seed
# (hit 0.797 < 0.8); voxels 1 and 2 are seeds of their own (the voxels are far apart, so no BFS
# neighbour is fitted: the neighbour path is covered by the scripted tests above).
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
    assert done == [0, 1, 2] and [v.reconstruction_id for v in mic.voxels] == GOLDEN_STATES
    assert [float(v.overlap_ratio) for v in mic.voxels] == GOLDEN_HIT
    assert bfs.stats["n_unresolved"] == 1
    assert [bfs.records[i].source for i in range(3)] == ["unresolved", "seed", "seed"]
    assert all(bfs.records[i].n_evals > 0 for i in range(3))
