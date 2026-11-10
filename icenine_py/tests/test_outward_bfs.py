"""OutwardSearch (Task G2, opt-in, default off): BFS neighbours and reseeds start with the Task F
fallback policy (F4; where F4 rejects, F5 then F4; 1-frame window).

Two kinds of test:
- The scripted voxel groups of ``test_bfs_refit`` (fake reconstructor, no images) with a fake
  policy that captures a voxel within 5 deg of its truth: the key's bookkeeping (where the local
  fit starts, which searches run, the counts), and that the key off never calls the policy.
- ``AdaptiveVoxelReconstructor.outward_search`` on the ThreeVoxels toy problem (real images): the
  policy equals the composition of the library's ``RayGaussNewton`` and ``ShellSearch`` with Task
  F's settings, its counts match the cost functions' evaluation counters, and it is deterministic.
"""

import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest

import icenine.mic_file as mic_module
import icenine.reconstructor as reconstructor_module
from icenine import cost_functions as CF
from icenine.config_file import ConfigFile
from icenine.mic_file import ReconstructionState
from icenine.orientation_search import SearchCandidate, SearchParameters
from icenine.reconstructor import OutwardPolicyResult, _get_voxel_vertices

# tests/ is on sys.path under pytest
import bfs_scripted
from test_bfs_refit import FakeReconstructor, RejectingFull, _angle, _info, _run, _rz, _states

FITTED, REFIT = ReconstructionState.FITTED, ReconstructionState.REFIT
CAPTURE_DEG = 5.0
needs_c = pytest.mark.skipif(not CF._HAS_C_RASTERIZE, reason="C extension missing")


class OutwardFake(FakeReconstructor):
    """The fake of ``test_bfs_refit`` plus a policy: within ``CAPTURE_DEG`` of the voxel's truth it
    returns the truth (accepted), else the start (rejected, after the F5 fallback). Every policy
    call and every local fit's (voxel, start) is logged."""

    outward_calls: List[Dict[str, Any]] = []
    local_starts: List[Tuple[int, np.ndarray]] = []

    def outward_search(self, voxel_vertices: Any, phase_index: int, start: np.ndarray) -> Any:
        i = self._idx(voxel_vertices)
        a = _angle(start, self.truth[i])
        ok = a < CAPTURE_DEG
        OutwardFake.outward_calls.append(dict(idx=i, angle=a, accepted=ok))
        R = self.truth[i].astype(np.float64) if ok else np.asarray(start, dtype=np.float64)
        return OutwardPolicyResult(
            orientation=R,
            accepted=ok,
            strict_cost=0.1 if ok else 0.97,
            fallback=not ok,
            n_strict=2 if ok else 10,
            n_screen=0 if ok else 500,
            n_smooth=20 if ok else 80,
            n_jac=5,
            wall_s=0.0,
        )

    def local_optimization(self, voxel_vertices: Any, phase_index: int, initial_orientation: Any,
                           rng: Any = None, cma_sigma0_deg: Any = None,
                           cma_max_evals: Any = None) -> SearchCandidate:  # fmt: skip
        R0 = np.asarray(initial_orientation, dtype=np.float64).copy()
        OutwardFake.local_starts.append((self._idx(voxel_vertices), R0))
        return super().local_optimization(
            voxel_vertices, phase_index, initial_orientation, rng, cma_sigma0_deg, cma_max_evals
        )


class OutwardRejectingFull(OutwardFake, RejectingFull):
    """Policy of ``OutwardFake``; full searches logged (and rejected for ``bad``) as in
    ``RejectingFull``."""


class Shifted(OutwardRejectingFull):
    """A local fit of voxel 1 or 2 from more than 5 deg out ends 3 deg short of the truth (and is
    rejected), so the voxel's stored fit, the start of a later reseed, is within the policy's
    reach."""

    def local_optimization(self, voxel_vertices: Any, phase_index: int, initial_orientation: Any,
                           rng: Any = None, cma_sigma0_deg: Any = None,
                           cma_max_evals: Any = None) -> SearchCandidate:  # fmt: skip
        i = self._idx(voxel_vertices)
        out = super().local_optimization(voxel_vertices, phase_index, initial_orientation, rng,
                                         cma_sigma0_deg, cma_max_evals)  # fmt: skip
        if i in (1, 2) and _angle(initial_orientation, self.truth[i]) > CAPTURE_DEG:
            R = self.truth[i] @ _rz(3.0).T
            return SearchCandidate(R.astype(np.float64), 0.7, _info(0.3))
        return out


class NoPolicy(FakeReconstructor):
    def outward_search(self, *a: Any, **k: Any) -> Any:  # pragma: no cover - must not run
        raise AssertionError("the policy ran with OutwardSearch off")


def _reset(bad: Tuple[int, ...] = ()) -> None:
    OutwardFake.outward_calls, OutwardFake.local_starts = [], []
    RejectingFull.bad, RejectingFull.full_calls = bad, []


# ------------------------------------------------------------------------------ config key
_CFG = Path(__file__).parent.parent.parent / "Examples" / "Example2.ThreeVoxels"
_CFG = _CFG / "ConfigFiles" / "ReconstructQ8.config"


def _config_with(tmp_path: Path, extra: str) -> ConfigFile:
    if not _CFG.exists():
        pytest.skip("ThreeVoxels config not available")
    p = tmp_path / "x.config"
    p.write_text(_CFG.read_text() + "\n" + extra + "\n")
    return ConfigFile.from_file(str(p))


def test_outward_search_key_default_off_and_parsed(tmp_path: Path) -> None:
    assert SearchParameters().outward_search is False
    assert SearchParameters.from_config(_config_with(tmp_path, "")).outward_search is False
    on = SearchParameters.from_config(_config_with(tmp_path, "OutwardSearch 1"))
    off = SearchParameters.from_config(_config_with(tmp_path, "OutwardSearch 0"))
    assert on.outward_search is True and off.outward_search is False


@pytest.mark.parametrize("line", ["OutwardSearch 2", "OutwardSearch"])
def test_outward_search_key_bad_values(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


# ------------------------------------------------------------------------------ key off
def test_key_off_never_calls_the_policy(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _run((3.0, 3.5, 3.2), monkeypatch, fake=NoPolicy, bfs_reseed=True)
    assert _states(mic) == [FITTED] * 4
    assert bfs.outward_log == []
    assert not any(k.startswith("outward_") for k in bfs.stats)


def test_key_off_scripted_grid_still_matches_the_golden() -> None:
    """The scripted 14 x 14 grid (golden recorded before the Phase D options) with the key
    written explicitly off: the whole fit sequence is unchanged."""
    gold = np.load(Path(__file__).parent / "data" / "bfs_grid_golden.npz")
    out = bfs_scripted.run_grid(
        reconstructor_module, mic_module, SearchCandidate, SearchParameters,
        {"outward_search": False},
    )  # fmt: skip
    assert np.array_equal(out["done"], gold["mc_done"])
    assert np.array_equal(out["log_idx"], gold["mc_log_idx"])
    assert float(out["rng_after"]) == float(gold["mc_rng_after"])


# ------------------------------------------------------------------------------ neighbours
def test_accepted_policy_moves_the_local_fit_start(monkeypatch: pytest.MonkeyPatch) -> None:
    """Voxels 1, 2 are 3 and 3.5 deg from seed 0: outside the fake local fit's 1 deg capture, but
    within the policy's 5 deg. With the key on the policy runs first on each neighbour, its end
    point (the truth) is the start of the local fit, and every voxel is fitted from seed 0."""
    _reset()
    bfs, mic = _run((3.0, 3.5, 3.2), monkeypatch, fake=OutwardFake, outward_search=True)
    assert _states(mic) == [FITTED] * 4
    assert bfs.stats["n_seeds"] == 1 and bfs.stats["n_neighbor_accepted_first"] == 3
    assert sorted(c["idx"] for c in OutwardFake.outward_calls) == [1, 2, 3]
    assert all(c["accepted"] for c in OutwardFake.outward_calls)
    for i, start in OutwardFake.local_starts:
        assert _angle(start, FakeReconstructor.truth[i]) < 1e-6
    s = bfs.stats
    assert s["outward_neighbor_runs"] == 3 and s["outward_neighbor_accepted"] == 3
    assert s["outward_neighbor_fallback"] == 0 and s["outward_neighbor_smooth"] == 60
    # the policy's strict + screen evaluations are counted with the local fit's
    assert s["n_evals_neighbor"] == 3 * (100 + 2) and bfs.records[1].n_evals == 102
    assert sorted(e["voxel"] for e in bfs.outward_log) == [1, 2, 3]
    assert all(e["phase"] == "neighbor" for e in bfs.outward_log)


def test_rejected_policy_falls_through_to_the_inherited_start(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Voxels 1, 2 are 8 and 8.5 deg from seed 0 (outside the policy too): the policy is rejected
    and the local fit starts from the inherited orientation, as with the key off."""
    _reset()
    bfs, mic = _run((8.0, 8.5, 0.5), monkeypatch, fake=OutwardFake, outward_search=True)
    assert _states(mic) == [FITTED, REFIT, REFIT, FITTED]
    rejected = [c["idx"] for c in OutwardFake.outward_calls if not c["accepted"]]
    assert sorted(rejected) == [1, 2]
    for i, start in OutwardFake.local_starts:
        assert _angle(start, _rz(0.0)) < 1e-6  # the seed's orientation, inherited
    s = bfs.stats
    assert s["outward_neighbor_fallback"] == 2 and s["outward_neighbor_screen"] == 1000
    assert s["n_evals_neighbor"] == 2 * (100 + 10 + 500)


def test_retry_does_not_rerun_the_policy(monkeypatch: pytest.MonkeyPatch) -> None:
    """cma mode with the retry: the wider-sigma retry starts from the inherited orientation and
    does not call the policy again (one policy run per fit)."""
    _reset()
    bfs, _ = _run((8.0,), monkeypatch, fake=OutwardFake, outward_search=True, local_optimizer="cma")
    assert len(OutwardFake.outward_calls) == 1
    assert bfs.stats["n_retry_attempted_neighbor"] == 1


# ------------------------------------------------------------------------------ reseeds
def test_reseed_policy_rejected_runs_the_full_search(monkeypatch: pytest.MonkeyPatch) -> None:
    """Voxels 1, 2 are 6 and 6.5 deg from seed 0, beyond the policy's 5 deg: both neighbour fits
    are rejected and keep seed 0's orientation. Voxel 3 is a seed of its own. The cluster {1, 2}
    is reseeded on voxel 1, whose start (6 deg out) the policy rejects: the full search runs. Its
    forced revisit then fits voxel 2 from 0.5 deg through the policy."""
    _reset()
    bfs, mic = _run(
        (6.0, 6.5, 0.2), monkeypatch, fake=OutwardRejectingFull, outward_search=True,
        bfs_reseed=True,
    )  # fmt: skip
    assert _states(mic) == [FITTED] * 4
    assert RejectingFull.full_calls == [0, 3, 1]  # two seeds, then the reseed's full search
    log = [e for e in bfs.outward_log if e["phase"] == "reseed"]
    assert len(log) == 1 and log[0]["voxel"] == 1
    assert not log[0]["accepted"] and log[0]["seed_test"] is None
    s = bfs.stats
    assert s["outward_reseed_runs"] == 1 and "outward_reseed_used" not in s
    assert s["n_evals_reseed"] == 10 + 500 + 1000  # the policy, then the full search
    assert bfs.records[2].source == "revisit"


def test_reseed_policy_accepted_replaces_the_full_search(monkeypatch: pytest.MonkeyPatch) -> None:
    """As above, but the failed neighbour fits end 3 deg short of the truth, so the reseed's start
    is within the policy's reach: the policy and one local fit replace the full search."""
    _reset()
    bfs, mic = _run(
        (6.0, 6.5, 0.2), monkeypatch, fake=Shifted, outward_search=True, bfs_reseed=True
    )
    assert _states(mic) == [FITTED] * 4
    assert RejectingFull.full_calls == [0, 3]  # the reseed of voxel 1 made no full search
    s = bfs.stats
    assert s["n_reseeds"] == 1 and s["n_reseed_accepted"] == 1 and s["outward_reseed_used"] == 1
    (entry,) = bfs.reseed_log
    assert entry["voxel"] == 1 and entry["accepted"]
    # the reseed's cost: the policy's strict + screen evaluations and one local fit
    assert entry["full_search_evaluations"] == 2 + 100
    assert s["n_evals_reseed"] == 2 + 100 and s["n_evals_seed"] == 2000
    assert bfs.records[1].reseeded and bfs.records[1].source == "seed"
    log = [e for e in bfs.outward_log if e["phase"] == "reseed"]
    assert len(log) == 1 and log[0]["accepted"] and log[0]["seed_test"] is True
    assert log[0]["finisher_evals"] == 100


def test_reseed_seed_test_rejection_falls_back_to_the_full_search(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The policy accepts voxel 1's reseed but the overlap of voxel 1 is always below the seed
    test: the full search runs too (a voxel is final only after a full search), is rejected, and
    voxel 1 is failed. The cluster's other voxel 2 is reseeded through the policy, and its forced
    revisit fits voxel 1 (local fits use their own overlap, not the seed test's)."""

    class Bad(Shifted):
        def evaluate_overlap(self, R: np.ndarray, vertices: Any, phase: int) -> Any:
            return _info(0.3 if self._idx(vertices) == 1 else 0.95)

    _reset()
    bfs, mic = _run((6.0, 6.5, 0.2), monkeypatch, fake=Bad, outward_search=True, bfs_reseed=True)
    assert RejectingFull.full_calls == [0, 3, 1]
    log = [e for e in bfs.outward_log if e["phase"] == "reseed"]
    assert [e["voxel"] for e in log] == [1, 2]
    assert [e["seed_test"] for e in log] == [False, True]
    s = bfs.stats
    assert s["outward_reseed_seed_test_rejected"] == 1 and s["outward_reseed_used"] == 1
    assert s["n_reseed_rejected"] == 1 and s["n_reseed_accepted"] == 1
    assert _states(mic) == [FITTED] * 4


# ------------------------------------------------------------------------------ real policy
def _toy() -> Any:
    from test_findoptimal_refactor import _build

    return _build()


@needs_c
def test_policy_on_the_toy_matches_the_library_composition(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """From 2 deg the policy is F4 alone (accepted, two strict evaluations); from 30 deg F4 is
    rejected and F5 then F4 run (outer shell cut to 0.9 deg to keep the test short). Both equal
    the same composition built by hand from ``RayGaussNewton`` and ``ShellSearch`` with Task F's
    settings, and the counts equal the cost functions' counters."""
    from icenine.outward_search import (
        RayGaussNewton,
        ShellSearch,
        angle_deg,
        rotate_left,
        screen_fn_from,
        smooth_fn_from,
    )

    monkeypatch.setattr(reconstructor_module, "OUTWARD_R_MAX_DEG", 0.9)
    rec, voxel, R_true = _toy()
    verts = _get_voxel_vertices(voxel)
    d = np.array([0.3, -0.8, 0.52])
    d /= np.linalg.norm(d)

    near = rotate_left(R_true, math.radians(2.0) * d)
    res = rec.outward_search(verts, voxel.phase, near)
    strict_cf, screen_cf, sc = rec._outward_costs()
    assert res.accepted and not res.fallback and res.n_strict == 2 and res.n_screen == 0
    assert angle_deg(res.orientation, R_true) < 0.1
    assert res.n_jac >= 1 and res.n_smooth > res.n_jac

    def strict(R: np.ndarray) -> float:
        return float(strict_cf.evaluate(np.asarray(R, np.float64), verts, voxel.phase).cost)

    gn = RayGaussNewton(sc, verts, strict, voxel.phase, tau_schedule_deg=(0.25,),
                        accept_cost=0.8, frame_window=1)  # fmt: skip
    ref = gn.solve(near)
    assert np.array_equal(ref.R, res.orientation) and ref.strict_cost == res.strict_cost

    far = rotate_left(R_true, math.radians(30.0) * d)
    n0_strict, n0_screen = strict_cf.eval_count, screen_cf.eval_count
    res2 = rec.outward_search(verts, voxel.phase, far)
    assert rec._outward_costs()[0] is strict_cf  # cached per image set
    assert not res2.accepted and res2.fallback and res2.n_screen > 0
    assert strict_cf.eval_count - n0_strict == res2.n_strict
    assert screen_cf.eval_count - n0_screen == res2.n_screen
    shells = ShellSearch(
        screen_fn_from(screen_cf, verts, voxel.phase), strict,
        smooth_fn_from(sc, verts, voxel.phase, "E4,tau=0.25"), spacing_deg=0.3, r_max_deg=0.9,
        accept_cost=0.8,
    )  # fmt: skip
    r5 = shells.search(far)
    r54 = gn.solve(r5.R)
    assert np.array_equal(r54.R, res2.orientation) and r54.accepted == res2.accepted
    assert res2.n_strict == 2 + r5.n_strict + 2 and res2.n_screen == r5.n_screen
    # deterministic: a second call gives the same end point
    res3 = rec.outward_search(verts, voxel.phase, far)
    assert np.array_equal(res3.orientation, res2.orientation)
