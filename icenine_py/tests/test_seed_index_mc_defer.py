"""AR5 step 22, phase B: the opt-in key ``SeedIndexMCDefer 1`` (off by default).

With ``SeedIndexMCStallHit k`` on, a seed-step MC refinement whose stall check goes on (the kept
best's hit ratio under ``SeedIndexStopHitRatio``) is PAUSED while the step has an untried candidate:
its state and generator are kept, the next candidate is tried, and the step stops at the first end
at the ratio or more. If no later candidate reaches it, the paused refinements are resumed in
candidate order, exactly where they stopped (no re-run), until one reaches it; else the lowest-cost
end is taken, as before. The last candidate is never paused. Rule ``DEF`` of
``scripts/autoresearch/ar5/defer_rules.py`` (frozen in ``benchmarks/autoresearch/ar5/defer/``),
moved into the library.

Behaviour tests with scripted refinements (no images): parsing and refused combinations; a paused
refinement resumes to exactly the end, cost, evaluation count and draws of the uninterrupted one
(the generator is continued, not restarted); a pause only at a check that goes on and never when it
is not allowed; an abandoned refinement is torn down; the seed step's order of fits, which ends are
taken and what is counted, in step types A, N and C; the keyword wiring of ``index_seed``; the BFS
with the key absent and written 0 identical.
"""

import threading
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pytest

from icenine import ray_indexing as RI
from icenine.orientation_search import SearchParameters, check_seed_index_mc_defer
from icenine.reconstructor import AdaptiveVoxelReconstructor, SeedMCRun

# tests/ is on sys.path under pytest
from test_bfs_refit import _run as _bfs_run
from test_bfs_seed_rescue import _reset_rescue
from test_seed_by_indexing import _config_with, _fingerprint
from test_seed_index_stall_hit import (
    BASE,
    ON,
    RESCUE,
    FakeRec,
    StallFake,
    _bfs_with,
    _hit_by_angle,
    _on,
)

DEFER = dict(ON, seed_index_mc_defer=1)
STOCK_SCRIPT = [0.5, 0.9, 0.9, 0.9, 0.4, 0.9, 0.9, 0.9, 0.3, 0.9]  # k = 3: checks follow gains


# ------------------------------------------------------------------------------ 1. the key
def test_default_is_off_and_the_key_parses(tmp_path: Path) -> None:
    assert SearchParameters().seed_index_mc_defer == 0
    for extra in ("", "SeedIndexMCDefer 0"):
        assert SearchParameters.from_config(_config_with(tmp_path, extra)).seed_index_mc_defer == 0
    on = SearchParameters.from_config(
        _config_with(
            tmp_path,
            "SeedByIndexing 1\nSeedIndexStopHitRatio 0.95\nSeedIndexMCStallHit 500\n"
            "SeedIndexMCDefer 1",
        )
    )
    assert on.seed_index_mc_defer == 1 and isinstance(on.seed_index_mc_defer, int)


@pytest.mark.parametrize(
    "line",
    ["SeedIndexMCDefer", "SeedIndexMCDefer 2", "SeedIndexMCDefer -1", "SeedIndexMCDefer x",
     "SeedIndexMCDefer 1.5", "SeedIndexMCDefer 1 1"],
)  # fmt: skip
def test_bad_values_raise_at_parse_time(tmp_path: Path, line: str) -> None:
    with pytest.raises(ValueError):
        _config_with(tmp_path, line)


@pytest.mark.parametrize("v", [2, -1, 1.0, True, "1", None])
def test_bad_values_raise_in_search_parameters(v: Any) -> None:
    with pytest.raises(ValueError):
        check_seed_index_mc_defer(v)
    with pytest.raises(ValueError):
        SearchParameters(**dict(ON, seed_index_mc_defer=v))


def test_the_key_needs_the_stall_check_and_says_so(monkeypatch: pytest.MonkeyPatch) -> None:
    with pytest.raises(ValueError) as e:
        _bfs_with(monkeypatch, dict(BASE, seed_index_mc_defer=1))
    msg = str(e.value)
    assert msg.startswith("SeedIndexMCDefer needs SeedIndexMCStallHit")
    assert "k >= 1" in msg
    with pytest.raises(ValueError):
        _bfs_with(monkeypatch, dict(seed_index_mc_defer=1))  # nothing else on
    _bfs_with(monkeypatch, DEFER)  # valid
    _bfs_with(monkeypatch, _on(seed_index_mc_defer=0))  # off needs nothing extra
    _bfs_with(monkeypatch, dict(seed_index_mc_defer=0))


# ------------------------------------------------------------------------------ 2. one refinement
def _run(
    costs: List[float], hits: Dict[int, float], k: int = 3, allow: bool = True
) -> Tuple[FakeRec, SeedMCRun]:
    rec = FakeRec(costs)
    run = SeedMCRun(rec, np.eye(3), 0.6, "V", 0, k, 0.95, _hit_by_angle(hits), allow)
    run.start()
    return rec, run


def _stock(costs: List[float], hits: Dict[int, float], k: int = 3) -> Tuple[FakeRec, Any]:
    rec = FakeRec(costs)
    return rec, rec.seed_mc_refine(np.eye(3), 0.6, "V", 0, k, 0.95, _hit_by_angle(hits))


def test_a_check_that_goes_on_pauses_the_refinement_with_its_state() -> None:
    # gain at 1 (0.5); non-gains 2, 3, 4 -> the check at evaluation 4 (hit 0.6: goes on): pause
    rec, run = _run(STOCK_SCRIPT, {0: 0.6})
    assert run.paused and not run.done
    end, best, evals, hit = run.partial()
    assert evals == 4 + 1 and best == 0.5 and hit == pytest.approx(0.6)  # the check is counted
    assert np.allclose(end, np.eye(3))  # the first lowest-cost orientation
    assert len(rec.draws) == 4  # four evaluations were made, then it waited
    run.abort()


def test_resume_equals_the_uninterrupted_refinement_and_does_not_rerun() -> None:
    hits = {0: 0.6, 4: 0.6, 8: 0.6}
    srec, (s_end, s_cost, s_evals, s_log) = _stock(STOCK_SCRIPT, hits)
    rec, run = _run(STOCK_SCRIPT, hits)
    n_before = len(rec.draws)
    assert run.paused and n_before < len(srec.draws)
    end, cost, evals = run.resume()
    assert run.done and not run.paused
    assert np.array_equal(end, s_end) and cost == s_cost and evals == s_evals
    # the generator was continued, not restarted: the draws are the stock refinement's, once each
    assert rec.draws == srec.draws and len(rec.draws) == len(STOCK_SCRIPT)
    assert rec.last_seed_stall_hit == [s_log] and run.log == s_log
    assert rec.last_eval_counts[1] == srec.last_eval_counts[1]


def test_a_check_that_stops_is_not_paused() -> None:
    srec, (s_end, s_cost, s_evals, s_log) = _stock(STOCK_SCRIPT, {0: 0.97})
    rec, run = _run(STOCK_SCRIPT, {0: 0.97})
    assert run.done and not run.paused and s_log["stopped"]
    assert np.array_equal(run.result[0], s_end) and run.result[1:] == (s_cost, s_evals)


def test_a_pause_is_not_made_when_it_is_not_allowed() -> None:
    srec, (s_end, s_cost, s_evals, s_log) = _stock(STOCK_SCRIPT, {0: 0.6})
    rec, run = _run(STOCK_SCRIPT, {0: 0.6}, allow=False)  # e.g. the last candidate
    assert run.done and not run.paused
    assert np.array_equal(run.result[0], s_end) and run.result[1:] == (s_cost, s_evals)
    assert rec.draws == srec.draws


def test_a_refinement_with_no_check_that_goes_on_is_not_paused() -> None:
    rec, run = _run([0.5, 0.4, 0.3, 0.2], {}, k=10)
    assert run.done and not run.paused


def test_an_abandoned_refinement_is_torn_down_and_its_log_says_paused() -> None:
    before = threading.active_count()
    rec, run = _run(STOCK_SCRIPT, {0: 0.6})
    assert threading.active_count() == before + 1  # the paused refinement waits
    run.abort()
    assert threading.active_count() == before
    assert not run.paused and run.done is False and run.aborted
    (log,) = rec.last_seed_stall_hit
    assert log["checks"] == [[4, 0.6, False]] and log["continued"] == 1
    assert log["evals"] == 5 and log["stopped"] is False  # what the step was charged


def test_an_error_in_the_refinement_reaches_the_caller() -> None:
    class Boom(FakeRec):
        def refine_from_candidates(self, *a: Any, **k: Any) -> Any:
            raise RuntimeError("scripted failure")

    with pytest.raises(RuntimeError, match="scripted failure"):
        SeedMCRun(Boom([0.5]), np.eye(3), 0.6, "V", 0, 3, 0.95, lambda R: 0.1, True).start()


# ------------------------------------------------------------------------------ 3. the seed step
class Step:
    """A scripted seed step for ``seed_by_indexing``: candidate j is the rotation by j degrees; its
    CMA end hits ``cma_hit[j]`` at cost ``cma_cost[j]``; its MC refinement is a ``SeedMCRun``-like
    object scripted by ``mc[j]``: ``(pause_hit or None, final_hit, final_cost, evals)``."""

    def __init__(self, n: int, mc: List[Tuple[Optional[float], float, float, int]]) -> None:
        self.n = n
        self.mc = mc
        self.events: List[Tuple[str, int]] = []
        self.hit_of: Dict[int, float] = {}

    def run(self, defer: bool, stop_hit: float = 0.95) -> RI.SeedResult:
        rs = [np.eye(3) @ _rz(float(j)) for j in range(self.n)]
        res = SimpleNamespace(R=np.array(rs))
        orig = RI.index_voxel
        RI.index_voxel = lambda *a, **k: (res, [])  # type: ignore[assignment]
        try:
            return RI.seed_by_indexing(
                None, None, None, "V", 0, SimpleNamespace(), lambda R: 0.5, self.n, self.n,
                stop_hit=stop_hit, hit_fn=self._hit, refine_fn=self._refine, defer=defer,
            )  # fmt: skip
        finally:
            RI.index_voxel = orig  # type: ignore[assignment]

    def _hit(self, R: np.ndarray) -> float:
        return self.hit_of[int(round(np.degrees(np.arctan2(R[1, 0], R[0, 0]))))]

    def _refine(self, R: np.ndarray, cost: float, allow_pause: bool = False) -> Any:
        j = int(round(np.degrees(np.arctan2(R[1, 0], R[0, 0]))))
        pause, h, c, ev = self.mc[j]
        end = R @ _rz(0.25)  # the MC end is a small step from the candidate
        self.hit_of[int(round(np.degrees(np.arctan2(end[1, 0], end[0, 0]) - 0.0)))] = h
        self.events.append(("mc", j))
        info: Dict[str, Any] = dict(cma_hit=0.1, mc_ran=True)
        if allow_pause and pause is not None:
            info["paused"] = True

            def resume() -> Any:
                self.events.append(("resume", j))
                return end, c, ev, dict(cma_hit=0.1, mc_ran=True)

            def abort() -> None:
                self.events.append(("abort", j))

            info["resume"], info["abort"] = resume, abort
            return R, 0.45, 100, info
        return end, c, ev, info


def _rz(deg: float) -> np.ndarray:
    a = np.radians(deg)
    return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])


def test_off_never_asks_for_a_pause_and_is_the_stock_loop() -> None:
    s = Step(3, [(0.6, 0.5, 0.4, 900), (0.6, 0.99, 0.2, 800), (None, 0.5, 0.3, 700)])
    out = s.run(defer=False)
    assert out.stopped_at == 1 and s.events == [("mc", 0), ("mc", 1)]
    assert out.finisher_evals == [900, 800]


def test_type_a_the_later_candidate_reaches_the_ratio_and_the_paused_one_is_never_resumed() -> None:
    s = Step(3, [(0.6, 0.5, 0.4, 900), (None, 0.99, 0.2, 800), (None, 0.5, 0.3, 700)])
    out = s.run(defer=True)
    assert out.stopped_at == 1 and out.taken == 1 and out.paused == [0] and out.resumed == []
    assert s.events == [("mc", 0), ("mc", 1), ("abort", 0)]
    assert out.finisher_evals == [100, 800]  # the paused part only, then the found refinement
    assert out.evals == 3 + 100 + 800  # the strict costs, then the fits


def test_type_c_none_later_reaches_so_the_paused_refinement_is_resumed_and_found() -> None:
    s = Step(3, [(0.6, 0.99, 0.2, 900), (None, 0.5, 0.4, 800), (None, 0.5, 0.3, 700)])
    out = s.run(defer=True)
    assert out.stopped_at == 0 and out.taken == 0 and out.resumed == [0]
    assert s.events == [("mc", 0), ("mc", 1), ("mc", 2), ("resume", 0)]  # later ones run first
    assert out.finisher_evals == [900, 800, 700]  # the resumed refinement counts as one run


def test_type_n_every_paused_refinement_is_resumed_in_order_and_the_step_is_the_stock_one() -> None:
    mc = [(0.6, 0.5, 0.40, 900), (0.6, 0.5, 0.30, 800), (None, 0.5, 0.35, 700)]
    stock = Step(3, mc).run(defer=False)
    s = Step(3, mc)
    out = s.run(defer=True)
    assert out.stopped_at == -1 and out.resumed == [0, 1] and out.paused == [0, 1]
    assert [e for e in s.events if e[0] == "resume"] == [("resume", 0), ("resume", 1)]
    assert out.finisher_evals == stock.finisher_evals and out.evals == stock.evals
    assert out.finisher_costs == stock.finisher_costs and out.taken == stock.taken
    for a, b in zip(out.finisher_R, stock.finisher_R):
        assert np.array_equal(a, b)


def test_the_last_candidate_is_never_asked_to_pause() -> None:
    asked: List[bool] = []

    class S2(Step):
        def _refine(self, R: np.ndarray, cost: float, allow_pause: bool = False) -> Any:
            asked.append(allow_pause)
            return super()._refine(R, cost, allow_pause)

    S2(3, [(0.6, 0.5, 0.4, 900)] * 3).run(defer=True)
    assert asked == [True, True, False]


def test_resumes_stop_at_the_first_end_at_the_ratio_and_later_ones_are_torn_down() -> None:
    mc = [(0.6, 0.5, 0.4, 900), (0.6, 0.99, 0.3, 800), (0.6, 0.99, 0.2, 700), (None, 0.5, 0.5, 600)]
    s = Step(4, mc)
    out = s.run(defer=True)
    assert out.stopped_at == 1 and out.resumed == [0, 1]
    assert ("resume", 2) not in s.events and ("abort", 2) in s.events
    assert out.taken == 1


# ------------------------------------------------------------------------------ 4. wiring
class _IndexCtx2:
    index_seed = AdaptiveVoxelReconstructor.index_seed

    def __init__(self) -> None:
        self.last_seed_stall_hit: List[Dict[str, Any]] = []

    def _index_context(self, phase: int) -> Any:
        return (
            "SC",
            "CRYSTAL",
            SimpleNamespace(),
            SimpleNamespace(evaluate=lambda R, v, p: SimpleNamespace(cost=0.5)),
        )

    def evaluate_overlap(self, R: Any, v: Any, p: int) -> Any:
        return SimpleNamespace(
            pixel_overlap=10, pixel_on_detector=100, peak_overlap=1, peak_on_detector=10
        )


def test_index_seed_passes_defer_to_the_seed_step_only_when_asked(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    got: List[Dict[str, Any]] = []

    def fake(rec: Any, sc: Any, crystal: Any, verts: Any, phase: int, params: Any, strict: Any,
             n_cost: int, n_refine: int, **kw: Any) -> Any:  # fmt: skip
        got.append(kw)
        return "RESULT"

    monkeypatch.setattr(RI, "seed_by_indexing", fake)
    rec = _IndexCtx2()
    rec.index_seed("V", 0, stop_hit=0.95, mc_stall_hit=500)
    rec.index_seed("V", 0, stop_hit=0.95, mc_stall_hit=500, mc_defer=True)
    rec.index_seed("V", 0, stop_hit=0.95, cma_first_evals=200, mc_stall_hit=500, mc_defer=True)
    assert not got[0].get("defer") and got[1]["defer"] is True and got[2]["defer"] is True
    with pytest.raises(ValueError):
        rec.index_seed("V", 0, stop_hit=0.95, mc_defer=True)  # no stall check
    with pytest.raises(ValueError):
        rec.index_seed("V", 0, stop_hit=0.0, mc_stall_hit=500, mc_defer=True)


# ------------------------------------------------------------------------------ 5. the BFS
class DeferFake(StallFake):
    """``StallFake`` whose ``index_seed`` also leaves the pause lists of a deferring step."""

    def index_seed(self, voxel_vertices: Any, phase_index: int, **kw: Any) -> RI.SeedResult:
        out = super().index_seed(voxel_vertices, phase_index, **kw)
        if kw.get("mc_defer"):
            out.paused, out.resumed = [0, 2], [0]
        return out


@pytest.mark.parametrize("rescue", [False, True])
def test_key_absent_and_written_0_are_identical_and_send_nothing(
    monkeypatch: pytest.MonkeyPatch, rescue: bool
) -> None:
    runs = []
    extra = RESCUE if rescue else {}
    for params in ({}, dict(seed_index_mc_defer=0)):
        _reset_rescue(failing=(0,), rescuable=(0,))
        StallFake.kws = []
        bfs, mic = _bfs_run((3.0, 3.5, 3.2), monkeypatch, fake=StallFake, **ON, **extra, **params)
        runs.append((_fingerprint(bfs, mic), [tuple(sorted(e)) for e in bfs.index_log]))
        assert all("mc_defer" not in kw for kw in StallFake.kws)
    assert runs[0] == runs[1]
    assert not any("defer" in k for k in runs[1][0]["stats"])
    assert all("defer" not in e for e in runs[1][1])


def test_key_on_counts_pauses_and_resumes_and_keeps_the_rest(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _reset_rescue(failing=(0,), rescuable=(0,))
    off, mic_off = _bfs_run((0.5,), monkeypatch, fake=DeferFake, **ON)
    _reset_rescue(failing=(0,), rescuable=(0,))
    StallFake.kws = []
    on, mic_on = _bfs_run((0.5,), monkeypatch, fake=DeferFake, **DEFER)
    assert StallFake.kws and all(kw["mc_defer"] is True for kw in StallFake.kws)
    n = len(on.index_log)
    s = on.stats
    assert s["index_defer_steps"] == n and s["index_defer_pauses"] == 2 * n
    assert s["index_defer_resumes"] == n
    assert all(e["defer"] == dict(paused=[0, 2], resumed=[0]) for e in on.index_log)
    fo, fn = _fingerprint(off, mic_off), _fingerprint(on, mic_on)
    fn["stats"] = {k: v for k, v in fn["stats"].items() if "defer" not in k}
    assert fo == fn
