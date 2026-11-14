"""AR5 step 12, phase B: the opt-in key ``SeedIndexMCStallHit k`` (off by default).

In the indexed seed step's MC refinement of a candidate (MC stage, variance stage and final
evaluation, one counter), after ``k`` strict-cost evaluations with no strict gain the hit ratio of
the kept best (the first lowest-cost orientation) is evaluated once and counted as one evaluation;
at ``SeedIndexStopHitRatio`` or more the refinement ends there, else it goes on and checks again
``k`` evaluations after the next strict gain. Rule ``PH(k)`` of
``scripts/autoresearch/ar5/stallhit_rules.py`` (frozen in
``benchmarks/autoresearch/ar5/stallhit/frozen.json``), moved into the library.

Scripted tests (fakes, no images): defaults, parsing, bad values, refused combinations; the
counting wrapper (stop at the ratio, go on under it, re-check after a gain, a gap of exactly k,
first lowest cost kept and copied, one evaluation per check); that it equals the script's
``patch_evaluate`` and the counting of ``stallhit_run.StallHitPick._fit`` on random cost
sequences; ``seed_mc_refine`` against the stock refinement; the wiring in ``index_seed``; the BFS
with the key absent, written 0 and on, with ``SeedRescueCMA`` off and on.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pytest

import icenine.reconstructor as reconstructor_module
from icenine import ray_indexing as RI
from icenine.orientation_search import (
    SearchCandidate,
    SearchParameters,
    check_seed_index_mc_stall_hit,
)
from icenine.reconstructor import (
    AdaptiveVoxelReconstructor,
    SeedMCStallHit,
    SeedStallHitStop,
)

# tests/ is on sys.path under pytest
from test_bfs_refit import _run, _rz
from test_bfs_seed_rescue import RescueFake, _reset_rescue
from test_seed_by_indexing import _config_with, _fingerprint, _reset

BASE: Dict[str, Any] = dict(seed_by_indexing=True, seed_index_stop_hit=0.95)
ON: Dict[str, Any] = dict(BASE, seed_index_mc_stall_hit=3000)
S11 = dict(
    bfs_neighbor_cma_first=True,
    seed_index_cma_first=200,
    bfs_neighbor_cma_no_fallback=True,
    seed_index_n_cost=20,
    seed_index_n_refine=5,
    boundary_pass=True,
)
RESCUE = dict(seed_rescue_sigma_deg=0.25, seed_rescue_evals=1000)


# ------------------------------------------------------------------------------ 1. defaults
def test_defaults_and_parsing(tmp_path: Path) -> None:
    assert SearchParameters().seed_index_mc_stall_hit == 0
    for extra in ("", "SeedIndexMCStallHit 0"):
        p = SearchParameters.from_config(_config_with(tmp_path, extra))
        assert p.seed_index_mc_stall_hit == 0, extra
    on = SearchParameters.from_config(
        _config_with(
            tmp_path, "SeedByIndexing 1\nSeedIndexStopHitRatio 0.95\nSeedIndexMCStallHit 3000"
        )
    )
    assert on.seed_index_mc_stall_hit == 3000 and isinstance(on.seed_index_mc_stall_hit, int)


@pytest.mark.parametrize(
    "line",
    [
        "SeedIndexMCStallHit",
        "SeedIndexMCStallHit -1",
        "SeedIndexMCStallHit 2.5",
        "SeedIndexMCStallHit x",
        "SeedIndexMCStallHit 3000 5",
    ],
)
def test_bad_values_raise_at_parse_time(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


@pytest.mark.parametrize("k", [-1, 2.5, True, "3000", None])
def test_bad_values_raise_in_search_parameters(k: Any) -> None:
    with pytest.raises(ValueError):
        SearchParameters(**dict(BASE, seed_index_mc_stall_hit=k))
    with pytest.raises(ValueError):
        check_seed_index_mc_stall_hit(k)


@pytest.mark.parametrize("k", [0, 1, 3000, np.int64(5)])
def test_good_values_pass_the_check(k: Any) -> None:
    check_seed_index_mc_stall_hit(k)


# ------------------------------------------------------------------------------ 2. combinations
def _bfs_with(monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any]) -> None:
    from test_bfs_refit import FakeReconstructor

    monkeypatch.setattr(reconstructor_module, "AdaptiveVoxelReconstructor", FakeReconstructor)
    setup = SimpleNamespace(
        sample=SimpleNamespace(get_mic=lambda: None),
        config=SimpleNamespace(min_acceleration_threshold=0.8),
        search_params=SearchParameters(**kw),
    )
    reconstructor_module.BFSReconstruction(setup)  # type: ignore[arg-type]


def _on(*drop: str, **change: Any) -> Dict[str, Any]:
    d = {k: v for k, v in ON.items() if k not in drop}
    d.update(change)
    return d


@pytest.mark.parametrize(
    "kw, message",
    [
        (
            _on("seed_by_indexing", "seed_index_stop_hit"),
            "SeedIndexMCStallHit needs SeedByIndexing 1",
        ),
        (_on("seed_index_stop_hit"), "SeedIndexMCStallHit needs SeedIndexStopHitRatio > 0"),
        (_on(local_optimizer="cma"), "SeedIndexMCStallHit needs LocalOptimizer mc"),
        (_on(use_hybrid_optimizer=True), "SeedIndexMCStallHit cannot be combined with"),
    ],
)
def test_refused_combinations_raise_before_any_search(
    monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any], message: str
) -> None:
    with pytest.raises(ValueError) as e:
        _bfs_with(monkeypatch, kw)
    assert str(e.value).startswith(message), str(e.value)


@pytest.mark.parametrize("kw", [ON, dict(ON, **S11), dict(ON, **S11, **RESCUE), dict(ON, **RESCUE)])
def test_valid_combinations_build(monkeypatch: pytest.MonkeyPatch, kw: Dict[str, Any]) -> None:
    _bfs_with(monkeypatch, kw)


def test_key_off_never_needs_the_seed_keys(monkeypatch: pytest.MonkeyPatch) -> None:
    _bfs_with(monkeypatch, dict(seed_index_mc_stall_hit=0))


# ------------------------------------------------------------------------------ 3. the wrapper
class _Hits:
    """Hit ratio per call (in order), recording the orientations it is asked about."""

    def __init__(self, hits: List[float]) -> None:
        self.hits = list(hits)
        self.asked: List[np.ndarray] = []

    def __call__(self, R: np.ndarray) -> float:
        self.asked.append(np.array(R, copy=True))
        return float(self.hits.pop(0))


def _feed(m: SeedMCStallHit, costs: List[float]) -> Optional[int]:
    """Feed ``costs`` (orientation i = a 3x3 filled with i); the evaluation number of a stop."""
    for i, c in enumerate(costs):
        try:
            m.feed(np.full((3, 3), float(i)), c)
        except SeedStallHitStop:
            return i + 1
    return None


def test_wrapper_stops_at_the_ratio_after_k_non_gains() -> None:
    m = SeedMCStallHit(3, 0.95, _Hits([0.95]))  # exactly at the ratio stops
    stop = _feed(m, [0.5, 0.4, 0.4, 0.6, 0.4, 0.3])
    assert stop == 5 and m.stopped and m.n == 5 and m.continued == 0
    assert m.checks == [[5, 0.95, True]]
    assert m.best == 0.4 and np.all(m.best_R == 1.0)  # the first 0.4 (evaluation 2) is kept


def test_wrapper_goes_on_under_the_ratio_and_checks_again_after_the_next_gain() -> None:
    hits = _Hits([0.9, 0.97])
    m = SeedMCStallHit(2, 0.95, hits)
    # gains at 1, 2; non-gains 3, 4 -> check at 4 (0.9, go on); 5, 6 no check (disarmed);
    # gain at 7; non-gains 8, 9 -> check at 9 (0.97, stop)
    stop = _feed(m, [0.5, 0.4, 0.4, 0.5, 0.6, 0.4, 0.3, 0.3, 0.9])
    assert stop == 9 and m.continued == 1
    assert m.checks == [[4, 0.9, False], [9, 0.97, True]]
    assert np.all(hits.asked[0] == 1.0) and np.all(hits.asked[1] == 6.0)


def test_wrapper_without_a_gain_makes_no_second_check() -> None:
    m = SeedMCStallHit(2, 0.95, _Hits([0.5]))
    assert _feed(m, [0.5] + [0.6] * 20) is None
    assert m.checks == [[3, 0.5, False]] and m.continued == 1 and not m.stopped


def test_wrapper_a_gap_of_exactly_k_makes_no_check() -> None:
    # gains at 1 and 4 (k = 3): 2 non-gains between them, no check; then 3 non-gains: a check
    m = SeedMCStallHit(3, 0.95, _Hits([0.99]))
    stop = _feed(m, [0.5, 0.6, 0.6, 0.4, 0.6, 0.6, 0.6])
    assert stop == 7 and m.checks == [[7, 0.99, True]]


def test_wrapper_keeps_a_copy_of_the_first_lowest_cost_orientation() -> None:
    hits = _Hits([0.99])
    m = SeedMCStallHit(2, 0.95, hits)
    R = np.eye(3)
    m.feed(R, 0.2)
    R[0, 0] = 7.0  # the caller reuses its array: the kept best must not move
    m.feed(np.zeros((3, 3)), 0.2)  # a tie does not replace it
    with pytest.raises(SeedStallHitStop):
        m.feed(np.ones((3, 3)), 0.3)
    assert np.array_equal(m.best_R, np.eye(3)) and np.array_equal(hits.asked[0], np.eye(3))
    assert m.best_R is not None and m.best_R.dtype == np.float64


def test_wrapper_wrap_counts_every_call_and_returns_the_value() -> None:
    m = SeedMCStallHit(5, 0.95, _Hits([]))
    seen: List[Tuple[Any, Any]] = []

    def evaluate(orientation: Any, voxel_vertices: Any = None, phase_index: int = 0) -> Any:
        seen.append((voxel_vertices, phase_index))
        return SimpleNamespace(cost=0.1 * len(seen))

    w = m.wrap(evaluate)
    out = w(np.eye(3), "V", 2)
    out2 = w(orientation=np.eye(3), voxel_vertices="W", phase_index=1)
    assert out.cost == pytest.approx(0.1) and out2.cost == pytest.approx(0.2)
    assert seen == [("V", 2), ("W", 1)] and m.n == 2 and m.best == pytest.approx(0.1)


def test_wrapper_refuses_k_0() -> None:
    with pytest.raises(ValueError):
        SeedMCStallHit(0, 0.95, _Hits([]))


# ------------------------------------------------------------------------------ 4. seed_mc_refine
class FakeCF:
    """A local cost function: ``evaluate`` returns the next scripted cost (the final evaluation of
    a refinement is the fake refinement's own call). Counts calls in ``eval_count``."""

    def __init__(self, costs: List[float]) -> None:
        self.costs = list(costs)
        self.eval_count = 0
        self.final_cost: Optional[float] = None

    def evaluate(self, orientation: Any, voxel_vertices: Any = None, phase_index: int = 0) -> Any:
        self.eval_count += 1
        if self.final_cost is not None:
            c = self.final_cost
        else:
            c = self.costs.pop(0)
        return SimpleNamespace(cost=float(c))


class FakeRec:
    """Stand-in for ``AdaptiveVoxelReconstructor`` in ``seed_mc_refine``: ``refine_from_candidates``
    walks the scripted costs (orientation i: the start rotated by i degrees about z, drawn in order;
    a draw from the generator per step) and ends with the final evaluation of its best. With no
    ``local_cost_fn`` it is the standalone call (sets ``_last_local_evals``)."""

    seed_mc_refine = AdaptiveVoxelReconstructor.seed_mc_refine
    last_eval_counts = AdaptiveVoxelReconstructor.last_eval_counts

    def __init__(self, costs: List[float]) -> None:
        self.script = list(costs)
        self.cfs: List[FakeCF] = []
        self.last_seed_stall_hit: List[Dict[str, Any]] = []
        self._last_global_evals = 0
        self._last_local_evals = 0
        self.draws: List[float] = []

    def _make_local_cost_fn(self) -> FakeCF:
        cf = FakeCF(self.script)
        self.cfs.append(cf)
        return cf

    def refine_from_candidates(
        self,
        candidates: List[SearchCandidate],
        voxel_vertices: Any,
        phase_index: int = 0,
        rng: Any = None,
        local_cost_fn: Any = None,
    ) -> SearchCandidate:
        standalone = local_cost_fn is None
        cf = self._make_local_cost_fn() if standalone else local_cost_fn
        R0 = np.asarray(candidates[0].orientation, np.float64)
        best_R, best = R0, float("inf")
        for i in range(len(self.script)):
            self.draws.append(float(rng.random()))
            R = R0 @ _rz(float(i)).astype(np.float64)
            c = cf.evaluate(R, voxel_vertices, phase_index).cost
            if c < best:
                best_R, best = R, c
        cf.final_cost = best
        final = cf.evaluate(
            orientation=best_R, voxel_vertices=voxel_vertices, phase_index=phase_index
        )
        if standalone:
            self._last_global_evals, self._last_local_evals = 0, cf.eval_count
        return SearchCandidate(best_R, float(final.cost))


def _hit_by_angle(table: Dict[int, float]) -> Any:
    """Hit ratio of an orientation from its angle (degrees, about z) to the identity start."""

    def hit(R: np.ndarray) -> float:
        deg = int(round(np.degrees(np.arctan2(R[1, 0], R[0, 0]))))
        return float(table.get(deg, 0.5))

    return hit


def test_seed_mc_refine_uncut_equals_the_stock_refinement() -> None:
    # gains at 1, 2, 4; then 3 non-gains with the final evaluation: under k = 4, no check
    costs = [0.5, 0.4, 0.45, 0.3, 0.35, 0.36]
    stock = FakeRec(costs)
    r = stock.refine_from_candidates(
        [SearchCandidate(np.eye(3), 0.6)], "V", 0, rng=np.random.default_rng(0)
    )
    ev_stock = stock.last_eval_counts[1]
    lib = FakeRec(costs)
    end, c, ev, log = lib.seed_mc_refine(np.eye(3), 0.6, "V", 0, 4, 0.95, _hit_by_angle({}))
    assert np.array_equal(end, r.orientation) and c == r.cost and ev == ev_stock == 7
    assert log == dict(n=7, evals=7, checks=[], continued=0, stopped=False)
    assert lib.draws == stock.draws  # the refinement's own default_rng(0)
    assert lib.last_seed_stall_hit == [log] and lib.last_eval_counts[1] == 7


def test_seed_mc_refine_a_check_that_goes_on_is_counted() -> None:
    costs = [0.5, 0.4, 0.45, 0.46, 0.47, 0.3, 0.31]  # check at 5 (best at 1 deg), then a gain
    lib = FakeRec(costs)
    end, c, ev, log = lib.seed_mc_refine(np.eye(3), 0.6, "V", 0, 3, 0.95, _hit_by_angle({1: 0.9}))
    assert c == 0.3 and ev == 8 + 1  # 7 loop + final evaluation + one counted check
    assert log["checks"] == [[5, 0.9, False]] and log["continued"] == 1 and not log["stopped"]


def test_seed_mc_refine_a_stop_keeps_the_best_and_counts_the_check() -> None:
    costs = [0.5, 0.4, 0.45, 0.46, 0.47, 0.3, 0.31, 0.32, 0.33, 0.1]
    lib = FakeRec(costs)
    hit = _hit_by_angle({1: 0.9, 5: 0.96})
    end, c, ev, log = lib.seed_mc_refine(np.eye(3), 0.6, "V", 0, 3, 0.95, hit)
    # check at 5 goes on (0.9); gain at 6 (5 deg); check at 9 stops (0.96): the later 0.1 is lost
    assert c == 0.3 and np.allclose(end, _rz(5.0))
    assert ev == 9 + 1 + 1 and log["stopped"] and log["checks"][-1] == [9, 0.96, True]
    assert lib.cfs[0].eval_count == 9  # nothing evaluated after the stop


# ------------------------------------------------------------------------------ 5. vs the script
def _script_fit(costs: List[float], hits: Dict[int, float], k: int) -> Tuple[Any, float, int, Any]:
    """The script's arm H fit: ``stallhit_rules.patch_evaluate`` on the cost function class and
    the counting of ``stallhit_run.StallHitPick._fit`` (cut: kept best, re-evaluated cost,
    ``n + 1 + continued``; uncut: the stock end, ``evals + continued``)."""
    sys.path.insert(0, str(Path(__file__).parent.parent / "scripts" / "autoresearch" / "ar5"))
    import stallhit_rules as SR

    class CF(FakeCF):
        pass

    meter = SR.Meter()
    CF.evaluate = SR.patch_evaluate(FakeCF.evaluate, meter)  # type: ignore[method-assign]
    rec = FakeRec(costs)
    rec._make_local_cost_fn = lambda: CF(rec.script)  # type: ignore[method-assign]
    meter.begin(k)
    meter.hit_fn = _hit_by_angle(hits)
    try:
        r = rec.refine_from_candidates(
            [SearchCandidate(np.eye(3), 0.6)], "V", 0, rng=np.random.default_rng(0)
        )
        R, c, ev = r.orientation, r.cost, rec.last_eval_counts[1] + meter.continued
    except SR.StallHitStop:
        R = meter.best_R
        cf = FakeCF([])
        cf.final_cost = meter.best  # the re-evaluation of the kept best gives its cost
        c = float(cf.evaluate(R).cost)
        ev = meter.n + 1 + meter.continued
    meter.on = False
    return R, c, ev, [list(x) for x in meter.checks]


@pytest.mark.parametrize("seed", range(12))
def test_seed_mc_refine_equals_the_script_counting(seed: int) -> None:
    g = np.random.default_rng(100 + seed)
    n = int(g.integers(5, 60))
    costs = list(np.round(g.random(n), 1))  # many ties
    hits = {d: float(g.choice([0.5, 0.94999, 0.95, 0.99])) for d in range(n)}
    k = int(g.integers(1, 6))
    R0, c0, ev0, checks0 = _script_fit(costs, hits, k)
    lib = FakeRec(costs)
    R1, c1, ev1, log = lib.seed_mc_refine(np.eye(3), 0.6, "V", 0, k, 0.95, _hit_by_angle(hits))
    assert np.array_equal(np.asarray(R0), R1) and c0 == c1 and ev0 == ev1
    assert checks0 == log["checks"]


# ------------------------------------------------------------------------------ 6. index_seed
class _IndexCtx:
    """``index_seed`` on a stand-in: the index context and ``seed_by_indexing`` are faked so the
    test sees which refinement each candidate gets."""

    index_seed = AdaptiveVoxelReconstructor.index_seed

    def __init__(self) -> None:
        self.calls: List[Tuple[str, Dict[str, Any]]] = []
        self.last_seed_stall_hit: List[Dict[str, Any]] = [dict(stale=True)]
        self.last_eval_counts = (0, 40, 40)
        self.last_local_optimization_evals = 20

    def _index_context(self, phase: int) -> Any:
        strict_cf = SimpleNamespace(evaluate=lambda R, v, p: SimpleNamespace(cost=0.5))
        return "SC", "CRYSTAL", SimpleNamespace(), strict_cf

    def evaluate_overlap(self, R: Any, v: Any, p: int) -> Any:
        return SimpleNamespace(
            pixel_overlap=10, pixel_on_detector=100, peak_overlap=1, peak_on_detector=10
        )

    def seed_mc_refine(
        self, R: Any, cost: float, v: Any, p: int, k: int, tau: float, hit: Any
    ) -> Any:
        self.calls.append(("stall", dict(k=k, tau=tau, hit=hit(R))))
        return np.asarray(R), 0.2, 77, dict(stopped=True)

    def refine_from_candidates(self, cands: Any, v: Any, p: int, rng: Any = None) -> Any:
        self.calls.append(("stock", {}))
        return SearchCandidate(np.asarray(cands[0].orientation), 0.3)

    def local_optimization(self, **kw: Any) -> Any:
        self.calls.append(("cma", {}))
        return SearchCandidate(np.asarray(kw["initial_orientation"]), 0.4)


def _fake_seed_by_indexing(monkeypatch: pytest.MonkeyPatch) -> List[Any]:
    got: List[Any] = []

    def fake(rec: Any, sc: Any, crystal: Any, verts: Any, phase: int, params: Any, strict: Any,
             n_cost: int, n_refine: int, stop_hit: float = 0.0, hit_fn: Any = None,
             refine_fn: Any = None) -> Any:  # fmt: skip
        got.append(dict(refine_fn=refine_fn, stop_hit=stop_hit))
        if refine_fn is not None:
            got.append(refine_fn(np.eye(3), 0.5))
        return "RESULT"

    monkeypatch.setattr(RI, "seed_by_indexing", fake)
    return got


def test_index_seed_key_on_routes_the_finisher_through_the_stall_check(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    got = _fake_seed_by_indexing(monkeypatch)
    rec = _IndexCtx()
    out = rec.index_seed("V", 0, stop_hit=0.95, mc_stall_hit=3000)
    assert out == "RESULT" and rec.last_seed_stall_hit == []  # reset for this call
    assert [c[0] for c in rec.calls] == ["stall"]
    assert rec.calls[0][1] == dict(k=3000, tau=0.95, hit=pytest.approx(0.1))
    assert got[1][:3][2] == 77


def test_index_seed_key_on_with_cma_first_puts_the_check_in_the_mc_refinement_only(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    got = _fake_seed_by_indexing(monkeypatch)
    rec = _IndexCtx()
    rec.index_seed("V", 0, stop_hit=0.95, cma_first_evals=200, mc_stall_hit=3000)
    assert [c[0] for c in rec.calls] == ["cma", "stall"]  # CMA end under 0.95: MC with check
    end, c, ev, info = got[1]
    assert c == 0.2 and ev == 20 + 77 and info == dict(cma_hit=pytest.approx(0.1), mc_ran=True)


def test_index_seed_key_off_is_the_stock_path(monkeypatch: pytest.MonkeyPatch) -> None:
    got = _fake_seed_by_indexing(monkeypatch)
    rec = _IndexCtx()
    rec.index_seed("V", 0, stop_hit=0.95)
    assert got[0]["refine_fn"] is None and rec.calls == []
    rec.index_seed("V", 0, stop_hit=0.95, cma_first_evals=200)
    assert [c[0] for c in rec.calls] == ["cma", "stock"]
    assert rec.last_seed_stall_hit == [dict(stale=True)]  # untouched with the key off


def test_index_seed_key_on_needs_stop_hit() -> None:
    with pytest.raises(ValueError):
        _IndexCtx().index_seed("V", 0, stop_hit=0.0, mc_stall_hit=3000)


# ------------------------------------------------------------------------------ 7. BFS
class StallFake(RescueFake):
    """``RescueFake`` whose ``index_seed`` records its keywords and, with the key on, leaves one
    stall-check log per refined candidate (the first stopped, the second not)."""

    kws: List[Dict[str, Any]] = []

    def index_seed(self, voxel_vertices: Any, phase_index: int, **kw: Any) -> RI.SeedResult:
        StallFake.kws.append(dict(kw))
        out = super().index_seed(voxel_vertices, phase_index, **kw)
        if kw.get("mc_stall_hit"):
            self.last_seed_stall_hit = [
                dict(n=3000, evals=3002, checks=[[2999, 0.9, False], [3000, 0.97, True]],
                     continued=1, stopped=True),
                dict(n=600, evals=600, checks=[], continued=0, stopped=False),
            ]  # fmt: skip
        return out


@pytest.mark.parametrize("rescue", [False, True])
def test_key_absent_and_written_0_are_identical(
    monkeypatch: pytest.MonkeyPatch, rescue: bool
) -> None:
    runs = []
    extra = RESCUE if rescue else {}
    for params in ({}, dict(seed_index_mc_stall_hit=0)):
        _reset_rescue(failing=(0,), rescuable=(0,))
        StallFake.kws = []
        bfs, mic = _run((3.0, 3.5, 3.2), monkeypatch, fake=StallFake, **BASE, **extra, **params)
        runs.append((_fingerprint(bfs, mic), [tuple(sorted(e)) for e in bfs.index_log]))
        assert all("mc_stall_hit" not in kw for kw in StallFake.kws)
    assert runs[0] == runs[1]
    assert not any("stallhit" in k for k in runs[1][0]["stats"])
    assert all("stall_hit" not in e for e in runs[1][1])


@pytest.mark.parametrize("rescue", [False, True])
def test_key_on_logs_the_checks_and_keeps_the_rest(
    monkeypatch: pytest.MonkeyPatch, rescue: bool
) -> None:
    extra = RESCUE if rescue else {}
    _reset_rescue(failing=(0,), rescuable=(0,))
    StallFake.kws = []
    off, mic_off = _run((0.5,), monkeypatch, fake=StallFake, **BASE, **extra)
    _reset_rescue(failing=(0,), rescuable=(0,))
    StallFake.kws = []
    on, mic_on = _run((0.5,), monkeypatch, fake=StallFake, **ON, **extra)
    assert StallFake.kws and all(kw["mc_stall_hit"] == 3000 for kw in StallFake.kws)
    s = on.stats
    n = len(on.index_log)
    assert s["index_stallhit_refinements"] == 2 * n and s["index_stallhit_checks"] == 2 * n
    assert s["index_stallhit_stops"] == n and s["index_stallhit_went_on"] == n
    assert s["index_stallhit_evals_in_stopped"] == 3002 * n
    assert all(len(e["stall_hit"]) == 2 for e in on.index_log)
    # the fake's index results do not depend on the key: everything else is unchanged
    fo, fn = _fingerprint(off, mic_off), _fingerprint(on, mic_on)
    fn["stats"] = {k: v for k, v in fn["stats"].items() if "stallhit" not in k}
    assert fo == fn
    assert ("index_rescue_runs" in s) == rescue
