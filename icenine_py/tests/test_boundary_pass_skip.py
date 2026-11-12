"""AR5 step 7: the opt-in key ``BoundaryPassSkipHitRatio h`` (off at 0).

A boundary voxel whose own hit ratio (the hit ratio of its end before the pass) is at least ``h``
gets no try and keeps its end; its score is still recorded and counted. Scripted groups with the
fake reconstructor of ``test_bfs_ar3_keys`` (no images): defaults and parsing, bad values and
combinations, key off equals key absent (same tries, evaluations, map and stats), key on skips
exactly the voxels at or above ``h``, the counts, and the call from the BFS.
"""

from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List

import numpy as np
import pytest

import icenine.reconstructor as reconstructor_module
from icenine import boundary_pass as BP
from icenine.mic_file import MicFile
from icenine.orientation_search import SearchParameters
from icenine.reconstructor import BFSReconstruction

# tests/ is on sys.path under pytest
from test_bfs_ar3_keys import PassFake, _pass_bfs
from test_bfs_refit import FakeReconstructor, _angle, _run, _rz
from test_seed_by_indexing import IndexFake, _config_with, _reset

KEY = "BoundaryPassSkipHitRatio"


# ------------------------------------------------------------------------------ 1. parsing
def test_default_off_and_parsed(tmp_path: Path) -> None:
    assert SearchParameters().boundary_pass_skip_hit == 0.0
    assert SearchParameters.from_config(_config_with(tmp_path, "")).boundary_pass_skip_hit == 0.0
    zero = SearchParameters.from_config(_config_with(tmp_path, f"BoundaryPass 1\n{KEY} 0"))
    assert zero.boundary_pass_skip_hit == 0.0 and zero.boundary_pass is True
    on = SearchParameters.from_config(_config_with(tmp_path, f"BoundaryPass 1\n{KEY} 0.95"))
    assert on.boundary_pass_skip_hit == 0.95 and on.boundary_pass is True


@pytest.mark.parametrize("line", [f"{KEY} 1.5", f"{KEY} -0.1", KEY])
def test_bad_values_raise_at_parse_time(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


def test_bad_value_raises_in_search_parameters() -> None:
    with pytest.raises(ValueError):
        SearchParameters(boundary_pass=True, boundary_pass_skip_hit=1.5)
    with pytest.raises(ValueError):
        SearchParameters(boundary_pass=True, boundary_pass_skip_hit=-0.1)


@pytest.mark.parametrize(
    "kw",
    [
        dict(boundary_pass_skip_hit=0.95),  # needs BoundaryPass 1
        dict(
            seed_by_indexing=True, paint_grains=True, boundary_pass=True, boundary_pass_skip_hit=0.9
        ),
    ],
)
def test_invalid_combinations_raise_before_any_search(
    monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any]
) -> None:
    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", FakeReconstructor)
    setup = SimpleNamespace(
        sample=SimpleNamespace(get_mic=lambda: None),
        config=SimpleNamespace(min_acceleration_threshold=0.8),
        search_params=SearchParameters(**kw),
    )
    with pytest.raises(ValueError):
        BFSReconstruction(setup)  # type: ignore[arg-type]


def test_key_at_default_never_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    _reset()
    for extra in ({}, {"seed_by_indexing": True}):
        _run((0.5,), monkeypatch, fake=IndexFake, boundary_pass_skip_hit=0.0, **extra)


# ------------------------------------------------------------------------------ 2. the pass
def _summary(res: BP.BoundaryPassResult, mic: MicFile, bfs: Any) -> Dict[str, Any]:
    return dict(
        boundary_voxels=res.boundary_voxels,
        tries=[(t["voxel"], t["cand"], t["evals"], t["cost"], t["hit"]) for t in res.tries],
        ends=[t["end"].tobytes() for t in res.tries],
        old_score=res.old_score,
        changed=res.changed,
        n_evals_tries=res.n_evals_tries,
        n_evals_score=res.n_evals_score,
        orientations=[np.asarray(v.orientation).tobytes() for v in mic.voxels],
        stats={k: v for k, v in bfs.stats.items() if not k.startswith("wall")},
    )


def test_key_off_is_identical_to_key_absent(monkeypatch: pytest.MonkeyPatch) -> None:
    """``skip_hit`` left out, ``0.0`` given: the same tries, evaluations, scores, map and stats,
    and no ``boundary_skipped`` entry in the stats."""
    outs = []
    for kw in ({}, {"skip_hit": 0.0}):
        bfs, mic = _pass_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
        res = BP.run_boundary_pass(bfs, mic, **kw)
        outs.append(_summary(res, mic, bfs))
        assert "boundary_skipped" not in bfs.stats and res.skipped == []
        assert PassFake.fits  # the fake ran
    assert outs[0] == outs[1]


def test_key_on_skips_exactly_the_voxels_at_or_above_h(monkeypatch: pytest.MonkeyPatch) -> None:
    """Truth of every voxel is 2 degrees about z. Voxel 0 sits at 0 degrees (own hit 0.85); voxels
    1 to 3 sit at truth (own hit 0.95, a hit ratio computed by the same function as the pass's).
    Boundary voxels: 0, 1, 2. With ``h`` at the hit ratio of voxels 1 and 2 they are skipped (at
    least ``h``); with ``h`` a hair above they are not."""
    bfs, mic = _pass_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    own = PassFake(None)._hit(_rz(2.0).astype(np.float64), 1)
    res = BP.run_boundary_pass(bfs, mic, skip_hit=own)
    assert res.boundary_voxels == [0, 1, 2]
    assert res.skipped == [1, 2]
    assert [(t["voxel"], t["cand"]) for t in res.tries] == [(0, 1)]
    assert [f[0] for f in PassFake.fits] == [0]  # no fit was started for a skipped voxel
    assert sorted(res.old_score) == [0, 1, 2]  # scores are still recorded for all of them
    assert res.old_score[1][1] == pytest.approx(own)
    for i in (1, 2, 3):
        assert np.array_equal(mic.voxels[i].orientation, _rz(2.0).astype(np.float64))
    assert res.changed == [0] and _angle(mic.voxels[0].orientation, _rz(2.0)) < 1e-4
    s = bfs.stats
    assert s["boundary_voxels"] == 3 and s["boundary_skipped"] == 2
    assert s["boundary_tries"] == 1 and s["boundary_changed"] == 1
    assert s["n_evals_boundary_pass"] == 61 and s["n_evals_boundary_score"] == 3

    bfs2, mic2 = _pass_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    res2 = BP.run_boundary_pass(bfs2, mic2, skip_hit=float(np.nextafter(own, 1.0)))
    assert res2.skipped == [] and len(res2.tries) == 3
    assert bfs2.stats["boundary_skipped"] == 0  # key on: the entry is present, with 0


def test_score_evaluation_is_counted_the_same_with_and_without_the_key(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The own-hit-ratio evaluation of a voxel is the pass's old-end score: one evaluation per
    boundary voxel in ``n_evals_boundary_score``, skipped or not, key on or off."""
    counts = []
    for kw in ({}, {"skip_hit": 0.9}, {"skip_hit": 0.5}, {"skip_hit": 1.0}):
        bfs, mic = _pass_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
        res = BP.run_boundary_pass(bfs, mic, **kw)
        counts.append((res.n_evals_score, bfs.stats["n_evals_boundary_score"]))
    assert counts == [(3, 3)] * 4


def test_tries_of_kept_voxels_are_the_tries_of_the_full_pass(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The kept voxels' tries (candidate, end, cost, hit, evaluations) equal those of the full
    pass; the kept tries' evaluations are 61 each and the skipped voxels add none."""
    bfs, mic = _pass_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    full = BP.run_boundary_pass(bfs, mic)
    bfs2, mic2 = _pass_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    part = BP.run_boundary_pass(bfs2, mic2, skip_hit=0.9)
    keep = [t for t in full.tries if t["voxel"] not in part.skipped]
    assert len(keep) == len(part.tries) == 1
    for a, b in zip(keep, part.tries):
        assert (a["voxel"], a["cand"], a["evals"], a["cost"], a["hit"]) == (
            b["voxel"],
            b["cand"],
            b["evals"],
            b["cost"],
            b["hit"],
        )
        assert np.array_equal(a["end"], b["end"])
    assert part.n_evals_tries == 61 * len(part.tries)
    assert full.n_evals_tries - part.n_evals_tries == 61 * 2


def test_everything_skipped_leaves_the_map_alone(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _pass_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    before = [np.asarray(v.orientation).copy() for v in mic.voxels]
    res = BP.run_boundary_pass(bfs, mic, skip_hit=0.05)
    assert res.tries == [] and res.changed == [] and res.skipped == res.boundary_voxels
    assert all(np.array_equal(b, v.orientation) for b, v in zip(before, mic.voxels))
    assert bfs.stats["boundary_tries"] == 0 and bfs.stats["n_evals_boundary_pass"] == 0
    assert bfs.stats["boundary_skipped"] == 3 and bfs.boundary_log == []


# ------------------------------------------------------------------------------ 3. in the BFS
def test_the_bfs_passes_the_key_only_when_it_is_on(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: List[Dict[str, Any]] = []

    def spy(bfs: Any, mic: Any, **kw: Any) -> None:
        calls.append(kw)

    monkeypatch.setattr(BP, "run_boundary_pass", spy)
    _run((0.5, -0.5), monkeypatch, boundary_pass=True)
    _run((0.5, -0.5), monkeypatch, boundary_pass=True, boundary_pass_skip_hit=0.0)
    _run((0.5, -0.5), monkeypatch, boundary_pass=True, boundary_pass_skip_hit=0.95)
    _run((0.5, -0.5), monkeypatch, boundary_pass=False)
    assert calls == [{}, {}, {"skip_hit": 0.95}]


# ------------------------------------------------------------------------------ 4. aligned draws
# Variant 2 (decided after the first live result, D-20261112-N5): a skipped try consumes the one
# CMA seed it would have drawn, so every kept try gets the draw it gets in the full pass.
class DrawFake(PassFake):
    """``PassFake`` whose fit draws a CMA seed as the real ``CMAOptimizer`` does and records it."""

    draws: List[Any] = []

    def local_optimization(self, *a: Any, **k: Any) -> Any:
        rng = k["rng"]
        i = self._idx(k["voxel_vertices"])
        DrawFake.draws.append((i, int(rng.integers(0, BP.CMA_SEED_BOUND))))
        k = dict(k, rng=np.random.default_rng(0))  # PassFake's own rng.integers(0, 10) is private
        return super().local_optimization(*a, **k)


def _draw_run(monkeypatch: pytest.MonkeyPatch, **kw: Any) -> Any:
    bfs, mic = _pass_bfs(monkeypatch, [0.0, 2.0, 2.0, 2.0])
    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", DrawFake)
    bfs.reconstructor = DrawFake(None)
    DrawFake.draws = []
    res = BP.run_boundary_pass(bfs, mic, **kw)
    return res, list(DrawFake.draws)


def test_a_cut_that_skips_nothing_equals_key_off_try_by_try(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    off, d_off = _draw_run(monkeypatch)
    on, d_on = _draw_run(monkeypatch, skip_hit=1.0)  # no own hit ratio reaches 1.0 here
    assert on.skipped == [] and len(d_off) == 3
    assert d_on == d_off
    assert [(t["voxel"], t["cand"], t["cost"], t["hit"]) for t in on.tries] == [
        (t["voxel"], t["cand"], t["cost"], t["hit"]) for t in off.tries
    ]


def test_every_kept_try_gets_the_draw_it_gets_in_the_full_pass(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    full, d_full = _draw_run(monkeypatch)
    # case 1: voxels 1 and 2 are skipped, voxel 0 (own hit 0.85) is kept and comes first
    part, d_part = _draw_run(monkeypatch, skip_hit=0.9)
    assert part.skipped == [1, 2]
    assert d_part == [d for d in d_full if d[0] == 0]  # voxel 0 first in the pass: same draw
    # case 2: voxel 0 now fits its start and is skipped; the kept voxels come AFTER it
    PassFake.TRUTH_DEG = [0.0, 2.0, 2.0, 2.0]
    try:
        full2, d_full2 = _draw_run(monkeypatch)
        part2, d_part2 = _draw_run(monkeypatch, skip_hit=0.9)
    finally:
        PassFake.TRUTH_DEG = [2.0, 2.0, 2.0, 2.0]
    assert part2.skipped  # voxel 0 now fits its start: skipped, voxels after it are kept
    kept_full = [d for d in d_full2 if d[0] not in part2.skipped]
    assert d_part2 == kept_full and len(d_part2) < len(d_full2)


def test_key_off_draws_are_unchanged(monkeypatch: pytest.MonkeyPatch) -> None:
    a, d_a = _draw_run(monkeypatch)
    b, d_b = _draw_run(monkeypatch, skip_hit=0.0)
    assert d_a == d_b and len(d_a) == len(a.tries) == 3
    assert [t["end"].tobytes() for t in a.tries] == [t["end"].tobytes() for t in b.tries]


# ------------------------------------------------------------------------------ 5. draw premise
# The aligned draws rest on two facts about the real ``CMAOptimizer.optimize``: one try takes
# exactly one draw from the generator, and that draw is ``integers(0, BP.CMA_SEED_BOUND)``. If
# either changes, the skipped tries would consume a different stream and the key would silently
# fall back to unaligned draws (variant 1). This test ties both to the real optimizer.
class _StubCost:
    """A smooth cost on rotations (angle to a fixed target, radians); no images."""

    def __init__(self) -> None:
        self.target = _rz(0.3)
        self.calls = 0

    def evaluate(self, R: np.ndarray, verts: Any, phase: int) -> Any:
        self.calls += 1
        c = float(np.arccos(np.clip((np.trace(self.target.T @ R) - 1.0) / 2.0, -1.0, 1.0)))
        return SimpleNamespace(cost=c)


@pytest.mark.parametrize(
    "kw",
    [
        dict(sigma0_deg=BP.PASS_SIGMA0_DEG, max_evals=BP.PASS_MAX_EVALS),  # the pass's settings
        dict(sigma0_deg=0.2, max_evals=1),  # early return at the budget
        dict(sigma0_deg=0.2, max_evals=50, max_convergence_cost=10.0),  # early return, converged
    ],
)
@pytest.mark.parametrize("seed", [0, 7, 12345])
def test_one_cma_try_takes_exactly_one_seed_draw_of_the_pass_bound(
    kw: Dict[str, Any], seed: int
) -> None:
    from icenine.orientation_search import CMAOptimizer

    rng = np.random.default_rng(seed)
    rng.integers(0, 10, size=seed % 5)  # start from a generator that has already drawn
    expected = np.random.default_rng(seed)
    expected.integers(0, 10, size=seed % 5)
    cost = _StubCost()
    res = CMAOptimizer(cost_fn=cost, voxel_vertices=None, rng=rng, **kw).optimize(_rz(0.0))
    assert cost.calls == res.n_evals >= 1
    expected.integers(0, BP.CMA_SEED_BOUND)
    assert rng.bit_generator.state == expected.bit_generator.state
    # the next draw matches too (same position in the stream)
    assert rng.integers(0, 2**63 - 1) == expected.integers(0, 2**63 - 1)
