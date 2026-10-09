"""Opt-in CMA-ES local refinement (CMAOptimizer, SearchParameters.local_optimizer).

Pure tests use a synthetic SO(3) quadratic. The reconstruction tests use ThreeVoxels voxel 0 (as
test_findoptimal_refactor) and skip without its Python-simulated data. The default ("mc") path is
checked against golden numbers recorded BEFORE the switch was added (same environment guard as
test_findoptimal_refactor).
"""

import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, List

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from icenine.config_file import ConfigFile
from icenine.mic_file import ReconstructionState
from icenine.orientation_search import (
    CMAOptimizer,
    CMAResult,
    SearchCandidate,
    SearchParameters,
)
from icenine.reconstructor import BFSReconstruction, _get_voxel_vertices

# tests/ is on sys.path under pytest
from test_findoptimal_refactor import (  # noqa: F401  (golden_env is a fixture)
    EXAMPLE,
    _build,
    _rotvec_deg,
    golden_env,
)

# ---------------------------------------------------------------------------
# Synthetic SO(3) quadratic
# ---------------------------------------------------------------------------


class QuadCost:
    """cost = (misorientation angle to R_true in degrees)^2, counting evaluations."""

    def __init__(self, R_true: np.ndarray) -> None:
        self.R_true = R_true
        self.n = 0

    def evaluate(self, R: np.ndarray, vertices: Any = None, phase: int = 0) -> Any:
        self.n += 1
        ang = np.degrees(np.linalg.norm(Rotation.from_matrix(R @ self.R_true.T).as_rotvec()))
        return SimpleNamespace(cost=float(ang) ** 2)


def _quad(offset_deg: List[float] = (0.3, -0.2, 0.25)):  # type: ignore[assignment]
    R_true = Rotation.from_rotvec(np.radians([10.0, 20.0, -5.0])).as_matrix()
    start = Rotation.from_rotvec(np.radians(offset_deg)).as_matrix() @ R_true
    return QuadCost(R_true), start, R_true


def test_converges_on_so3_quadratic() -> None:
    f, start, R_true = _quad()
    res = CMAOptimizer(f, None, sigma0_deg=0.2, max_evals=1000).optimize(start, seed=0)
    assert isinstance(res, CMAResult) and isinstance(res, SearchCandidate)
    assert (
        np.degrees(np.linalg.norm(Rotation.from_matrix(res.orientation @ R_true.T).as_rotvec()))
        < 1e-3
    )
    assert res.cost < 1e-6 < 0.3**2
    assert res.n_evals == f.n <= 1000
    assert res.stop_reason  # tolx here, otherwise the budget


def test_never_worse_than_start_and_returns_lowest_evaluated() -> None:
    f, start, _ = _quad()
    seen: List[float] = []
    inner = f.evaluate

    def spy(R: np.ndarray, v: Any = None, p: int = 0) -> Any:
        info = inner(R, v, p)
        seen.append(info.cost)
        return info

    f.evaluate = spy  # type: ignore[method-assign]
    res = CMAOptimizer(f, None, max_evals=60).optimize(start, seed=3)
    assert res.cost == min(seen) and res.cost <= seen[0]
    assert res.overlap_info.cost == res.cost  # info of the best evaluation, no extra call


@pytest.mark.parametrize("budget", [1, 2, 7, 8, 9, 50, 123])
def test_respects_max_evals_exactly(budget: int) -> None:
    f, start, _ = _quad()
    res = CMAOptimizer(f, None, sigma0_deg=0.2, max_evals=budget, tolx=0.0).optimize(start, seed=1)
    assert f.n == budget and res.n_evals == budget
    assert res.stop_reason == "max_evals"


def test_deterministic_given_seed_and_given_generator() -> None:
    f, start, _ = _quad()
    a = CMAOptimizer(f, None, max_evals=80).optimize(start, seed=5)
    b = CMAOptimizer(f, None, max_evals=80).optimize(start, seed=5)
    c = CMAOptimizer(f, None, max_evals=80).optimize(start, seed=6)
    assert np.array_equal(a.orientation, b.orientation) and a.cost == b.cost
    assert not np.array_equal(a.orientation, c.orientation)
    # seeds drawn from a generator: same generator state -> same run, successive runs differ
    g1, g2 = np.random.default_rng(11), np.random.default_rng(11)
    r1 = [CMAOptimizer(f, None, rng=g1, max_evals=80).optimize(start) for _ in range(2)]
    r2 = [CMAOptimizer(f, None, rng=g2, max_evals=80).optimize(start) for _ in range(2)]
    assert np.array_equal(r1[0].orientation, r2[0].orientation)
    assert np.array_equal(r1[1].orientation, r2[1].orientation)
    assert not np.array_equal(r1[0].orientation, r1[1].orientation)


def test_does_not_touch_global_numpy_rng() -> None:
    f, start, _ = _quad()
    np.random.seed(123)
    before = np.random.get_state()
    CMAOptimizer(f, None, max_evals=60).optimize(start, seed=2)
    after = np.random.get_state()
    assert before[0] == after[0] and np.array_equal(before[1], after[1])
    assert before[2:] == after[2:]
    with pytest.raises(ValueError):
        CMAOptimizer(f, None, max_evals=60).optimize(start, seed=-1)


def test_converged_cost_stops_early() -> None:
    f, start, _ = _quad()
    res = CMAOptimizer(f, None, max_evals=1000, max_convergence_cost=0.05).optimize(start, seed=0)
    assert res.stop_reason == "converged_cost" and res.cost < 0.05 and res.n_evals < 1000


def test_invalid_options() -> None:
    f, _, _ = _quad()
    with pytest.raises(ValueError):
        CMAOptimizer(f, None, sigma0_deg=0.0)
    with pytest.raises(ValueError):
        CMAOptimizer(f, None, max_evals=0)
    with pytest.raises(ValueError):
        SearchParameters(local_optimizer="adam")
    with pytest.raises(ValueError):
        SearchParameters(cma_max_evals=1)


# ---------------------------------------------------------------------------
# Config parsing
# ---------------------------------------------------------------------------

_CFG = EXAMPLE / "ConfigFiles" / "ReconstructQ8.config"


def _config_with(tmp_path: Path, extra: str) -> ConfigFile:
    if not _CFG.exists():
        pytest.skip("ThreeVoxels config not available")
    p = tmp_path / "x.config"
    p.write_text(_CFG.read_text() + "\n" + extra + "\n")
    return ConfigFile.from_file(str(p))


def test_config_absent_keys_mean_mc(tmp_path: Path) -> None:
    sp = SearchParameters.from_config(_config_with(tmp_path, ""))
    assert sp.local_optimizer == "mc" and sp.cma_sigma0_deg == 0.2
    assert sp.cma_max_evals == 1000 and sp.cma_popsize is None


def test_config_cma_keys(tmp_path: Path) -> None:
    cfg = _config_with(tmp_path, "LocalOptimizer CMA\nCMASigma0 0.1\nCMAMaxEvals 400\nCMAPopSize 9")
    sp = SearchParameters.from_config(cfg)
    assert sp.local_optimizer == "cma"
    assert sp.cma_sigma0_deg == pytest.approx(0.1) and sp.cma_max_evals == 400
    assert sp.cma_popsize == 9


@pytest.mark.parametrize(
    "line",
    [
        "LocalOptimizer adam",
        "CMAPopSize 1",
        "CMASigma0 0",
        "CMASigma0 -1",
        "CMAMaxEvals 1",
        "LocalOptimizer",
    ],
)
def test_config_bad_values(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


# ---------------------------------------------------------------------------
# Default path unchanged (golden recorded before the switch existed)
# ---------------------------------------------------------------------------

LOCAL_START = [0.3, -0.2, 0.25]  # deg
GOLDEN_LOCAL = (
    np.array([0.006113616625963422, -0.008326472025432697, 0.4981822883419766]),
    0.4524590163934423,
)
# Re-recorded twice. (1) After the VarianceMinimizing restart fix (C++ restart semantics,
# MIGRATION_HISTORY "C++ vs Python variance stage"); before: [-0.058482767583009056,
# -0.009258290608073556, 0.08004536647020442], cost 0.5508196721311475, 470 evaluations. Over 30
# seeds the old and new stage give the same error (median 0.13 vs 0.14 deg) and cost (mean 0.303 vs
# 0.294). (2) After the port of MCOptimizer.optimize to C++ RandomRestartZeroTemp ("C++-faithful MC:
# reruns"): the FindOptimal MC draws differ; before: [0.05764184647751444, 0.0032931430107475605,
# -0.13882759611192275], cost 0.5901639344262296, 392 evaluations. One seed of a toy problem.
GOLDEN_REFINE = (
    np.array([-0.00795524046660142, -0.04733991639366214, 0.19163534607486127]),
    0.47322404371584703,
    304,  # local cost evaluations
)


def _start(R_true: np.ndarray, deg: List[float] = LOCAL_START) -> np.ndarray:
    return Rotation.from_rotvec(np.radians(deg)).as_matrix() @ R_true


def test_default_local_optimization_unchanged(golden_env) -> None:  # noqa: F811
    rec, voxel, R_true = _build()
    assert rec.params.local_optimizer == "mc"
    res = rec.local_optimization(
        _get_voxel_vertices(voxel), voxel.phase, _start(R_true), rng=np.random.default_rng(5)
    )
    np.testing.assert_allclose(_rotvec_deg(res.orientation, R_true), GOLDEN_LOCAL[0], atol=1e-9)
    assert res.cost == pytest.approx(GOLDEN_LOCAL[1], abs=1e-12)


def test_default_refine_from_candidates_unchanged(golden_env) -> None:  # noqa: F811
    rec, voxel, R_true = _build()
    res = rec.refine_from_candidates(
        [SearchCandidate(orientation=_start(R_true), cost=1.0)],
        _get_voxel_vertices(voxel),
        voxel.phase,
        rng=np.random.default_rng(3),
    )
    np.testing.assert_allclose(_rotvec_deg(res.orientation, R_true), GOLDEN_REFINE[0], atol=1e-9)
    assert res.cost == pytest.approx(GOLDEN_REFINE[1], abs=1e-12)
    assert rec.last_eval_counts == (0, GOLDEN_REFINE[2], GOLDEN_REFINE[2])


def test_explicit_mc_equals_default() -> None:
    a, voxel, R_true = _build()
    b, _, _ = _build()
    b.params.local_optimizer = "mc"
    args = (_get_voxel_vertices(voxel), voxel.phase, _start(R_true))
    ra = a.local_optimization(*args, rng=np.random.default_rng(5))
    rb = b.local_optimization(*args, rng=np.random.default_rng(5))
    assert np.array_equal(ra.orientation, rb.orientation) and ra.cost == rb.cost


# ---------------------------------------------------------------------------
# CMA on a real small case
# ---------------------------------------------------------------------------


def test_refine_from_candidates_cma_reaches_default_quality() -> None:
    rec, voxel, R_true = _build()

    def cands() -> List[SearchCandidate]:
        return [SearchCandidate(orientation=_start(R_true, [0.2, -0.1, 0.15]), cost=1.0)]

    vv = _get_voxel_vertices(voxel)
    res_mc = rec.refine_from_candidates(cands(), vv, voxel.phase, rng=np.random.default_rng(3))
    rec.params.local_optimizer = "cma"
    rec.params.cma_max_evals = 400
    res_cma = rec.refine_from_candidates(cands(), vv, voxel.phase, rng=np.random.default_rng(3))
    again = rec.refine_from_candidates(cands(), vv, voxel.phase, rng=np.random.default_rng(3))
    assert res_cma.cost <= res_mc.cost
    assert np.linalg.norm(_rotvec_deg(res_cma.orientation, R_true)) < 0.1
    g, loc, tot = rec.last_eval_counts
    assert g == 0 and 1 <= loc <= 401  # at most one CMA run (budget 400) plus the final overlap
    assert np.array_equal(res_cma.orientation, again.orientation)  # deterministic given the rng
    assert res_cma.overlap_info is not None and math.isfinite(res_cma.cost)


def test_cma_and_hybrid_are_exclusive() -> None:
    """Rejected up front from the parameters alone: before any search, and with no diff cost
    function (setup.diff_cost_fn is None here)."""
    rec, voxel, R_true = _build()
    assert rec.setup.diff_cost_fn is None
    rec.params.local_optimizer = "cma"
    rec.params.use_hybrid_optimizer = True
    vv = _get_voxel_vertices(voxel)
    with pytest.raises(ValueError, match="mutually exclusive"):
        rec.reconstruct_voxel(vv, voxel.phase, rng=np.random.default_rng(0))
    with pytest.raises(ValueError, match="mutually exclusive"):
        rec.refine_from_candidates(
            [SearchCandidate(orientation=_start(R_true), cost=1.0)], vv, voxel.phase
        )
    with pytest.raises(ValueError, match="mutually exclusive"):
        rec.local_optimization(vv, voxel.phase, _start(R_true))


def test_local_optimization_cma_uses_cma(monkeypatch: pytest.MonkeyPatch) -> None:
    rec, voxel, R_true = _build()
    rec.params.local_optimizer = "cma"
    rec.params.cma_max_evals = 300
    calls: List[int] = []
    orig = CMAOptimizer.optimize

    def spy(self: CMAOptimizer, R0: np.ndarray, seed: Any = None) -> CMAResult:
        out = orig(self, R0, seed)
        calls.append(out.n_evals)
        return out

    monkeypatch.setattr(CMAOptimizer, "optimize", spy)
    res = rec.local_optimization(
        _get_voxel_vertices(voxel), voxel.phase, _start(R_true), rng=np.random.default_rng(5)
    )
    assert len(calls) == 1 and calls[0] <= 300
    assert np.linalg.norm(_rotvec_deg(res.orientation, R_true)) < 0.1
    assert res.overlap_info is not None and res.cost == res.overlap_info.cost


def test_small_bfs_runs_end_to_end_with_cma(monkeypatch: pytest.MonkeyPatch) -> None:
    rec, voxel, R_true = _build()
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
    setup.search_params.local_optimizer = "cma"
    setup.search_params.cma_max_evals = 300
    calls: List[int] = []
    orig = CMAOptimizer.optimize

    def spy(self: CMAOptimizer, R0: np.ndarray, seed: Any = None) -> CMAResult:
        out = orig(self, R0, seed)
        calls.append(out.n_evals)
        return out

    monkeypatch.setattr(CMAOptimizer, "optimize", spy)
    for v in mic.voxels:
        v.reconstruction_id = ReconstructionState.NOT_VISITED
    done = BFSReconstruction(setup).reconstruct_sample(rng=np.random.default_rng(1))
    assert sorted(done) == [0, 1, 2]
    states = {v.reconstruction_id for v in mic.voxels}
    assert states <= {ReconstructionState.FITTED, ReconstructionState.REFIT}
    assert ReconstructionState.FITTED in states
    # seed refinement (FindOptimal) and neighbour refinement both used CMA, within the budget
    assert len(calls) >= 3 and max(calls) <= 300
    # voxel 0's reconstruction lands within 0.1 deg of its truth
    assert np.linalg.norm(_rotvec_deg(mic.voxels[0].orientation, truth[0])) < 0.1
