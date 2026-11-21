"""AR5 step 21, phase B: the opt-in key ``BoundaryPassScreen`` (off by default).

With the key on every boundary-pass try has a budget of one evaluation (a screen try: the
candidate's start is scored, the final evaluation makes 2) and the repair of ``BoundaryPassRepair``
takes the wide flag: contested voxels, plus voxels whose hit ratio over the median hit ratio of
their neighbours is under 0.95 and that have a pass candidate more than 1 degree from the kept end.
Scripted groups with the fake reconstructor of ``test_bfs_ar3_keys`` (no images): defaults and
parsing, invalid combinations, the wide flag equals the script rule's, key off equals key absent
(same calls, same draws), key on counts every evaluation, and the call from the BFS.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List

import numpy as np
import pytest

import icenine.reconstructor as reconstructor_module
from icenine import boundary_pass as BP
from icenine import boundary_repair as BR
from icenine.orientation_search import SearchCandidate, SearchParameters
from icenine.reconstructor import BFSReconstruction

# tests/ is on sys.path under pytest
from test_bfs_ar3_keys import PassFake, _info, _pass_bfs
from test_bfs_refit import FakeReconstructor, _angle, _run, _rz
from test_boundary_repair import RepairFake
from test_seed_by_indexing import _config_with

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts" / "autoresearch" / "ar5"))
import screen_rules as SR  # noqa: E402

KEY = "BoundaryPassScreen"
I3 = np.eye(3)[None]


# ------------------------------------------------------------------------------ 1. parsing
def test_default_off_and_parsed(tmp_path: Path) -> None:
    assert SearchParameters().boundary_pass_screen is False
    assert SearchParameters.from_config(_config_with(tmp_path, "")).boundary_pass_screen is False
    on = SearchParameters.from_config(
        _config_with(tmp_path, f"BoundaryPass 1\nBoundaryPassRepair 1\n{KEY} 1")
    )
    assert on.boundary_pass_screen is True and on.boundary_pass_repair is True
    off = SearchParameters.from_config(_config_with(tmp_path, f"{KEY} 0"))
    assert off.boundary_pass_screen is False


def test_bad_value_raises_at_parse_time(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, f"{KEY} maybe")


@pytest.mark.parametrize(
    "kw",
    [
        dict(boundary_pass=True, boundary_pass_screen=True),  # needs BoundaryPassRepair 1
        dict(boundary_pass_screen=True, boundary_pass_repair=True),  # needs BoundaryPass 1
        dict(
            boundary_pass=True,
            boundary_pass_repair=True,
            boundary_pass_screen=True,
            boundary_pass_skip_hit=0.9,  # no voxel is skipped by its own hit ratio
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


def test_constants_are_the_plans_and_the_scripts() -> None:
    assert BP.SCREEN_MAX_EVALS == 1 == SR.SCREEN_MAX_EVALS
    assert BR.WIDE_DEG == 1.0 == SR.WIDE_DEG and BR.RATIO_MAX == SR.RATIO_MAX == 0.95


def _real_local_optimization(monkeypatch: pytest.MonkeyPatch) -> Any:
    """The library's ``AdaptiveVoxelReconstructor.local_optimization`` on a stub ``self`` (real
    ``_make_cma_optimizer`` and ``_check_local_optimizer``, a counting fake cost function)."""
    from icenine.reconstructor import AdaptiveVoxelReconstructor as A

    state = SimpleNamespace(n=0)

    class CF:
        def __init__(self, **kw: Any) -> None:
            self.eval_count = 0

        def evaluate(self, orientation: Any, voxel_vertices: Any, phase_index: int = 0) -> Any:
            self.eval_count += 1
            return _info(0.5)

    monkeypatch.setattr(reconstructor_module, "VoxelCostFunction", CF)
    me = SimpleNamespace(
        params=SimpleNamespace(
            local_optimizer="mc",
            use_hybrid_optimizer=False,
            cma_sigma0_deg=0.2,
            cma_max_evals=1000,
            cma_popsize=None,
        ),  # fmt: skip
        setup=SimpleNamespace(
            exp_setup=SimpleNamespace(get_eta_limit=lambda: 0.0),
            simulator=None,
            detector_list=None,
            range_map=None,
            exp_data=None,
            sample=None,
            structure_list=None,
        ),  # fmt: skip
        _cost_flags=lambda: {},
        last_local_optimization_evals=0,
    )
    me._check_local_optimizer = lambda: A._check_local_optimizer(me)
    me._make_cma_optimizer = lambda *a, **k: A._make_cma_optimizer(me, *a, **k)

    def call(budget: Any, **kw: Any) -> Any:
        rng = np.random.default_rng(3)
        r = A.local_optimization(
            me, np.zeros((3, 3)), 0, np.eye(3), rng=rng, cma_sigma0_deg=0.05,
            cma_max_evals=budget, optimizer="cma", **kw,
        )  # fmt: skip
        return r, me.last_local_optimization_evals, rng

    return call


def test_a_budget_of_one_raises_without_cma_screen_and_runs_with_it(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    call = _real_local_optimization(monkeypatch)
    with pytest.raises(ValueError, match="cma_max_evals must be >= 2"):
        call(1)
    with pytest.raises(ValueError, match="cma_max_evals must be >= 2"):
        call(1, cma_screen=False)
    with pytest.raises(ValueError):
        call(0, cma_screen=True)  # only a budget of exactly 1 is the screen try
    r, evals, _rng = call(1, cma_screen=True)
    assert evals == 2  # the start's own evaluation plus the final one
    assert np.array_equal(r.orientation, np.eye(3))  # the end is the start
    r2, evals2, _ = call(2)  # key-off budgets are unchanged
    assert evals2 == 3


def test_cma_optimizer_with_budget_one_returns_the_start_and_draws_one_seed() -> None:
    from icenine.orientation_search import CMAOptimizer

    n = [0]

    class CF:
        def evaluate(self, R: np.ndarray, vertices: Any, phase: int) -> Any:
            n[0] += 1
            return _info(0.25)

    rng = np.random.default_rng(11)
    res = CMAOptimizer(CF(), None, 0, rng, sigma0_deg=0.05, max_evals=1).optimize(_rz(3.0))
    assert res.n_evals == 1 and n[0] == 1
    assert np.array_equal(res.orientation, _rz(3.0))
    ref = np.random.default_rng(11)
    ref.integers(0, 2**31 - 2)  # exactly one seed draw
    assert rng.bit_generator.state == ref.bit_generator.state


# ------------------------------------------------------------------------------ 2. wide flag
def test_the_wide_flag_equals_the_script_rule() -> None:
    rng = np.random.default_rng(7)
    n = 80
    R = np.stack([_rz(d).astype(np.float64) for d in rng.uniform(0, 4, n)])
    ratio = rng.uniform(0.6, 1.1, n)
    flag = rng.random(n) < 0.15
    tv = rng.integers(0, n, 120)
    te = np.stack([_rz(d).astype(np.float64) for d in rng.uniform(0, 4, 120)])
    for a, b in zip(
        BR.wide_flag(flag, ratio, tv, te, R, I3), SR.wide_flag(flag, ratio, tv, te, R, I3)
    ):
        assert np.array_equal(a, b)
    w, extra = BR.wide_flag(flag, ratio, tv[:0], te[:0], R, I3)
    assert np.array_equal(w, flag) and not extra.any()


def test_the_wide_ratio_is_relative_to_the_neighbours() -> None:
    """A voxel with a high own hit ratio but a higher neighbour median is wide; one with a low own
    hit ratio among equally low neighbours is not."""
    hit = np.array([0.9, 0.99, 0.99, 0.99, 0.5, 0.5])
    nbrs = [[1, 2, 3], [0], [0], [0], [5], [4]]
    ratio = BR.ratio_of(hit, nbrs)
    assert ratio[0] < 0.95 and ratio[4] >= 0.95  # 0.9 / 0.99 against 0.5 / 0.5
    R = np.stack([_rz(0.0)] * 6).astype(np.float64)
    tv = np.array([0, 4])
    te = np.stack([_rz(5.0), _rz(5.0)]).astype(np.float64)  # both have a far candidate
    wide, extra = BR.wide_flag(np.zeros(6, bool), ratio, tv, te, R, I3)
    assert wide.tolist() == [True, False, False, False, False, False] and extra[0]


# ------------------------------------------------------------------------------ 3. the pass
class ScreenFake(RepairFake):
    """``RepairFake`` that also takes the screen try: budget 1 needs ``cma_screen=True`` (as the
    library), draws one seed from the generator, scores the start, returns it, 2 evaluations."""

    TRUTH_DEG = [10.0, 2.0, 8.0, 8.0]
    calls: List[Dict[str, Any]] = []

    def local_optimization(  # type: ignore[override]
        self, voxel_vertices: Any, phase_index: int, initial_orientation: np.ndarray,
        rng: Any = None, cma_sigma0_deg: Any = None, cma_max_evals: Any = None,
        optimizer: Any = None, cma_screen: bool = False,
    ) -> SearchCandidate:  # fmt: skip
        ScreenFake.calls.append(dict(budget=cma_max_evals, screen=cma_screen))
        if cma_max_evals == 1:
            assert cma_screen is True and optimizer == "cma"
            rng.integers(0, 10)
            i = self._idx(voxel_vertices)
            R0 = np.asarray(initial_orientation, np.float64)
            self.last_local_optimization_evals = 2
            return SearchCandidate(R0, 1.0 - self._hit(R0, i), _info(self._hit(R0, i)))
        assert not cma_screen
        return super().local_optimization(
            voxel_vertices, phase_index, initial_orientation, rng, cma_sigma0_deg,
            cma_max_evals, optimizer,
        )  # fmt: skip


def _screen_bfs(monkeypatch: pytest.MonkeyPatch, deg: List[float]) -> Any:
    bfs, mic = _pass_bfs(monkeypatch, deg)
    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", ScreenFake)
    bfs.reconstructor = ScreenFake(None)
    ScreenFake.states = []
    ScreenFake.calls = []
    return bfs, mic


def test_key_off_makes_the_same_calls_and_draws(monkeypatch: pytest.MonkeyPatch) -> None:
    """Key absent and ``screen=False``: full tries, no ``cma_screen``, the stats of today."""
    outs = []
    for kw in ({}, {"screen": False}):
        bfs, mic = _screen_bfs(monkeypatch, [0.0, 2.0, 8.0, 8.0])
        res = BP.run_boundary_pass(bfs, mic, repair=True, **kw)
        outs.append(
            (
                [(t["voxel"], t["cand"], t["evals"], t["hit"]) for t in res.tries],
                res.repair.flagged,
                [np.asarray(v.orientation).tobytes() for v in mic.voxels],
                dict(ScreenFake.calls[0]),
                sorted(bfs.stats),
            )
        )
        assert all(c == dict(budget=60, screen=False) for c in ScreenFake.calls[: len(res.tries)])
        assert not any("screen" in k or "extra" in k for k in bfs.stats)
        assert res.repair.extra == []
    assert outs[0] == outs[1]


def test_key_on_screen_tries_and_the_wide_repair(monkeypatch: pytest.MonkeyPatch) -> None:
    """Voxel 0 starts at 0 degrees (truth 10); its candidates are 2 and 8 degrees (neighbours 1 and
    2 and 3; 2 and 3 are merged). Screen tries score the starts: 2 degrees hit 0.55, 8 degrees 0.85;
    the pass keeps 8 degrees. Margin 0.85 - 0.55 = 0.30: not contested. Its ratio 0.85 / 0.95 is
    under 0.95 and the 2-degree candidate is more than 1 degree from the kept end: wide."""
    bfs, mic = _screen_bfs(monkeypatch, [0.0, 2.0, 8.0, 8.0])
    res = BP.run_boundary_pass(bfs, mic, repair=True, screen=True)
    assert all(t["evals"] == 2 for t in res.tries) and res.n_evals_tries == 2 * len(res.tries)
    assert [c for c in ScreenFake.calls[: len(res.tries)]] == [dict(budget=1, screen=True)] * len(
        res.tries
    )
    for t in res.tries:  # the end of a screen try is its start, bit for bit
        assert np.array_equal(t["end"], _rz(float([0.0, 2.0, 8.0, 8.0][t["cand"]])).astype(float))
    rep = res.repair
    assert rep.flagged == [0] and rep.extra == [0]
    assert [(f["voxel"], f["start"]) for f in rep.fits][0] == (0, None)
    assert all(f["evals"] == 121 for f in rep.fits)
    assert rep.n_evals_fits == 121 * len(rep.fits) and rep.n_evals_score == 0
    assert _angle(mic.voxels[0].orientation, _rz(10.0)) < 1e-4  # the repair's fit reached it
    s = bfs.stats
    assert s["repair_extra"] == 1 and s["repair_flagged"] == 1
    assert s["n_evals_boundary_pass"] == 2 * len(res.tries)
    assert s["n_evals_boundary_repair"] == rep.n_evals


def test_key_on_without_a_wide_voxel_adds_none(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _screen_bfs(monkeypatch, [2.0, 2.0, 2.0, 2.0])
    ScreenFake.TRUTH_DEG = [2.0, 2.0, 2.0, 2.0]
    try:
        res = BP.run_boundary_pass(bfs, mic, repair=True, screen=True)
    finally:
        ScreenFake.TRUTH_DEG = [10.0, 2.0, 8.0, 8.0]
    assert res.repair.flagged == [] and res.repair.extra == [] and res.repair.n_evals == 0


def test_screen_needs_the_repair_and_no_skip(monkeypatch: pytest.MonkeyPatch) -> None:
    bfs, mic = _screen_bfs(monkeypatch, [0.0, 2.0, 8.0, 8.0])
    with pytest.raises(ValueError):
        BP.run_boundary_pass(bfs, mic, screen=True)
    with pytest.raises(ValueError):
        BP.run_boundary_pass(bfs, mic, repair=True, screen=True, skip_hit=0.9)
    assert PassFake is not None


# ------------------------------------------------------------------------------ 4. the BFS call
def test_the_bfs_passes_the_key_only_when_it_is_on(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: List[Dict[str, Any]] = []

    def spy(bfs: Any, mic: Any, **kw: Any) -> None:
        calls.append(kw)

    monkeypatch.setattr(BP, "run_boundary_pass", spy)
    _run((0.5, -0.5), monkeypatch, boundary_pass=True, boundary_pass_repair=True)
    _run(
        (0.5, -0.5), monkeypatch, boundary_pass=True, boundary_pass_repair=True,
        boundary_pass_screen=True,
    )  # fmt: skip
    assert calls == [
        {"skip_hit": 0.0, "repair": True},  # today's call, unchanged
        {"skip_hit": 0.0, "repair": True, "screen": True},
    ]
