"""AR5 step 15 (task ar5o): rule PB(k), a two-sided stall check in the seed step's MC refinement, on
synthetic sequences. Plan "AR5 step 15", phase A step 2: stop at 0.95 and under 0.80, go on between,
re-check after a gain, the first lowest cost kept, one evaluation per check, the stopped end kept,
the walk of a whole seed step, other draws untouched, guards on ``held2``.
"""

import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import band_guard as BG  # noqa: E402
import band_rules as BR  # noqa: E402
import stallhit_rules as SR  # noqa: E402


class _Out:
    def __init__(self, cost: float) -> None:
        self.cost = cost


class _FakeCost:
    """A cost function whose cost is the first coordinate of the orientation."""

    def evaluate(self, orientation: Any, *args: Any, **kw: Any) -> _Out:
        return _Out(float(np.asarray(orientation).ravel()[0]))


def _drive(
    costs: List[float],
    k: int,
    hit: Callable[[float], float],
    low: Optional[float] = BR.LOW,
    aux: Sequence[int] = (),
) -> Dict[str, Any]:
    """Feed ``costs`` through the wrapped cost function; ``hit(best cost)`` plays the hit ratio."""
    meter = BR.BandMeter()
    calls: List[float] = []

    def hit_fn(R: Any) -> float:
        calls.append(float(np.asarray(R).ravel()[0]))
        return hit(float(np.asarray(R).ravel()[0]))

    meter.hit_fn = hit_fn
    patched = BR.patch_evaluate(_FakeCost.evaluate, meter)
    fake = _FakeCost()
    meter.begin(k, low, aux)
    stopped = False
    try:
        for c in costs:
            patched(fake, np.full((3, 3), c))
    except SR.StallHitStop:
        stopped = True
    return dict(stopped=stopped, meter=meter, calls=calls)


def test_constants_are_the_frozen_values() -> None:
    assert BR.TAU == 0.95 and BR.LOW == 0.80 and BR.PH_K == 3000
    assert BR.FORMS == {"PB1000": 1000, "PB500": 500} and BR.RECORD_KS == (500, 1000)
    assert "PB(k)" in BR.RULE and "0.95" in BR.RULE and "0.80" in BR.RULE


def test_decide_is_two_sided_with_closed_high_end_and_open_low_end() -> None:
    assert BR.decide(0.95) == "high" and BR.decide(0.9499999) == "go"
    assert BR.decide(0.80) == "go" and BR.decide(0.7999999) == "low"
    assert BR.decide(0.1, low=None) == "go"  # PH(3000): no low side


def test_stops_high_at_0_95_and_the_check_is_the_final_evaluation() -> None:
    costs = [0.9, 0.5] + [0.7] * 10 + [0.4]
    r = _drive(costs, 5, lambda c: 0.96)
    m = r["meter"]
    assert r["stopped"] and not m.on and m.n == 7 and m.stop_kind == "high"
    assert m.checks == [[7, 0.96, True]] and m.continued == 0
    assert float(m.best_R[0, 0]) == 0.5  # the stopped end is the running best, not the later 0.4


def test_stops_low_under_0_80_and_keeps_the_best_at_the_check() -> None:
    costs = [0.9, 0.5] + [0.7] * 10 + [0.4]
    r = _drive(costs, 5, lambda c: 0.5)
    m = r["meter"]
    assert r["stopped"] and m.stop_kind == "low" and m.n == 7
    assert m.checks == [[7, 0.5, True]] and float(m.best_R[0, 0]) == 0.5


def test_stops_at_exactly_0_95_not_just_under_and_not_at_exactly_0_80() -> None:
    costs = [0.9, 0.5] + [0.7] * 10
    assert _drive(costs, 5, lambda c: 0.95)["stopped"]
    assert not _drive(costs, 5, lambda c: 0.9499999)["stopped"]
    assert not _drive(costs, 5, lambda c: 0.80)["stopped"]  # 0.80 is the seed test: passes
    assert _drive(costs, 5, lambda c: 0.7999999)["stopped"]


def test_goes_on_between_and_counts_the_check_as_an_evaluation() -> None:
    costs = [0.9, 0.5] + [0.7] * 10
    r = _drive(costs, 5, lambda c: 0.88)
    m = r["meter"]
    assert not r["stopped"] and m.on and m.n == len(costs)
    assert m.checks == [[7, 0.88, False]] and m.continued == 1
    assert len(r["calls"]) == 1  # one hit evaluation for the check; no second without a gain


def test_rechecks_after_a_gain_and_stops_when_the_new_best_leaves_the_band() -> None:
    costs = [0.9, 0.6] + [0.7] * 6 + [0.4] + [0.8] * 10  # stall at 7, gain at 9, stall at 14

    def hit(c: float) -> float:
        return 0.88 if c == 0.6 else 0.97

    r = _drive(costs, 5, hit)
    m = r["meter"]
    assert r["stopped"] and m.n == 14 and m.stop_kind == "high"
    assert m.checks == [[7, 0.88, False], [14, 0.97, True]] and m.continued == 1
    assert float(m.best_R[0, 0]) == 0.4 and r["calls"] == [0.6, 0.4]


def test_first_lowest_cost_is_kept_on_a_tie() -> None:
    costs = [0.9, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3]
    r = _drive(costs, 3, lambda c: 1.0)
    m = r["meter"]
    assert r["stopped"] and m.best == 0.3 and m.trace == [[1, 0.9], [2, 0.3]]
    assert m.counter.best_n == 2


def test_exactly_one_hit_evaluation_per_stall_check() -> None:
    costs = [0.9, 0.5] + [0.7] * 40
    r = _drive(costs, 5, lambda c: 0.9)
    assert len(r["calls"]) == len(r["meter"].checks) == 1
    costs = [0.9] + [0.8, 0.7, 0.6, 0.5, 0.4, 0.3] + [0.9] * 5 + [0.2] + [0.9] * 5
    r = _drive(costs, 5, lambda c: 0.9)
    assert len(r["calls"]) == len(r["meter"].checks) == 2


def test_a_gap_of_exactly_k_does_not_check() -> None:
    costs = [0.9, 0.5, 0.7, 0.7, 0.7, 0.7, 0.4]
    r = _drive(costs, 5, lambda c: 1.0)
    assert not r["stopped"] and r["calls"] == []


def test_ph_mode_without_a_low_side_never_stops_low() -> None:
    costs = [0.9, 0.5] + [0.7] * 10
    assert not _drive(costs, 5, lambda c: 0.1, low=None)["stopped"]
    assert _drive(costs, 5, lambda c: 0.99, low=None)["stopped"]


def test_off_only_records_and_calls_nothing() -> None:
    r = _drive([0.9, 0.5, 0.7, 0.7, 0.4], 0, lambda c: 0.1)
    m = r["meter"]
    assert not r["stopped"] and m.n == 5 and r["calls"] == []
    assert m.trace == [[1, 0.9], [2, 0.5], [5, 0.4]]


def test_wrapped_evaluate_ignores_calls_while_off_and_copies_the_orientation() -> None:
    meter = BR.BandMeter()
    patched = BR.patch_evaluate(_FakeCost.evaluate, meter)
    out = patched(_FakeCost(), np.full((3, 3), 0.7))
    assert out.cost == 0.7 and meter.n == 0
    meter.begin(5, BR.LOW)
    R = np.full((3, 3), 0.3)
    patched(_FakeCost(), R)
    R[:] = 9.0
    assert float(meter.best_R[0, 0]) == 0.3


def test_the_check_is_not_counted_by_the_wrapper_itself() -> None:
    meter = BR.BandMeter()
    patched = BR.patch_evaluate(_FakeCost.evaluate, meter)
    fake = _FakeCost()

    def hit_fn(R: Any) -> float:
        patched(fake, np.full((3, 3), 0.1))  # an overlap evaluation that reached the wrapper
        return 0.9

    meter.hit_fn = hit_fn
    meter.begin(2, BR.LOW)
    for c in (0.5, 0.6, 0.6):
        patched(fake, np.full((3, 3), c))
    assert meter.n == 3 and meter.continued == 1 and meter.on


def test_recording_counters_record_every_500_style_stall_and_never_act() -> None:
    # k = 0 (record only) with auxiliary stall sizes 3 and 5
    costs = [0.9, 0.5] + [0.7] * 6 + [0.4] + [0.8] * 6
    r = _drive(costs, 0, lambda c: 0.1, low=None, aux=(3, 5))
    m = r["meter"]
    assert not r["stopped"] and m.n == len(costs)
    assert m.events[3] == [[5, 0.5, 0.1], [12, 0.4, 0.1]]
    assert m.events[5] == [[7, 0.5, 0.1], [14, 0.4, 0.1]]
    assert all(e[2] == 0.1 for e in m.events[3])


def test_recording_does_not_change_the_acting_rule() -> None:
    """PH(3000)-style acting counter plus recording: same stop and checks as without recording."""
    costs = [0.9, 0.5] + [0.7] * 10
    plain = _drive(costs, 6, lambda c: 0.96, low=None)
    rec = _drive(costs, 6, lambda c: 0.96, low=None, aux=(2, 4))
    assert plain["stopped"] == rec["stopped"] is True
    assert plain["meter"].checks == rec["meter"].checks
    assert plain["meter"].n == rec["meter"].n
    assert [e[0] for e in rec["meter"].events[2]] == [4]  # fired once, then disarmed
    assert [e[0] for e in rec["meter"].events[4]] == [6]


def test_walk_events_stops_high_low_or_never() -> None:
    ev = [[1000, 0.5, 0.88], [2500, 0.3, 0.97]]
    w = BR.walk_events(ev, 4000)
    assert w["kind"] == "high" and w["evals"] == 2500 + 1 + 1 and w["continued"] == 1
    assert w["cost"] == 0.3 and w["i"] == 1
    ev = [[1000, 0.5, 0.4]]
    w = BR.walk_events(ev, 4000)
    assert w["kind"] == "low" and w["evals"] == 1001 and w["hit"] == 0.4
    w = BR.walk_events([[1000, 0.5, 0.9], [2200, 0.4, 0.85]], 3000)
    assert w["kind"] == "none" and w["evals"] == 3000 + 2 and w["continued"] == 2
    w = BR.walk_events([], 3000)
    assert w["kind"] == "none" and w["evals"] == 3000


def test_cut_entry_marks_a_stop_before_a_later_gain_as_moved() -> None:
    w = BR.walk_events([[1000, 0.5, 0.4]], 4000)
    assert BR.cut_entry(w, 0.5)["kind"] == "after_last"  # J's end is this best
    assert BR.cut_entry(w, 0.3)["kind"] == "moved"  # a later gain existed
    c = BR.cut_entry(w, 0.3)
    assert c["evals"] == 1001 and c["cost"] == 0.5 and c["hit"] == 0.4
    n = BR.cut_entry(BR.walk_events([], 3000), 0.3)
    assert n["kind"] == "none" and n["evals"] == 3000


def _fits() -> List[Dict[str, Any]]:
    return [
        dict(pos=0, kind="C", hit=0.60, cost=0.40, evals=100),
        dict(pos=0, kind="M", hit=0.97, cost=0.30, evals=5000),
        dict(pos=1, kind="C", hit=0.50, cost=0.45, evals=100),
        dict(pos=1, kind="M", hit=0.70, cost=0.35, evals=4000),
    ]


def test_whole_step_walk_with_a_stop_after_the_last_gain_keeps_the_path() -> None:
    w = BR.walk_events([[1000, 0.30, 0.97]], 4999)  # J's end (0.30) is the stopped best
    cut = BR.cut_entry(w, 0.30)
    r = SR.walk_step(_fits(), {0: cut})
    assert r["stopped"] and r["taken"] == 0 and not r["lost"]
    assert r["evals"] == 100 + 1001  # stopped at the first candidate's MC end


def test_whole_step_walk_with_a_low_stop_before_a_later_gain_moves_the_end() -> None:
    w = BR.walk_events([[1000, 0.38, 0.5]], 4999)  # stall best 0.38 at 0.5; J's end was 0.30
    cut = BR.cut_entry(w, 0.30)
    assert cut["kind"] == "moved"
    r = SR.walk_step(_fits(), {0: cut})
    assert r["kept"][0]["hit"] == 0.5 and r["kept"][0]["cost"] == 0.38
    assert r["lost"] and not r["stopped"]  # J's stop at 0.95 is lost: the step goes on
    assert r["taken"] == 1  # the lowest-cost kept end: candidate 1's 0.35 against 0.38


def test_guard_refuses_held2_as_a_sample() -> None:
    for name in ("held2", "held2x"):
        with pytest.raises(SystemExit):
            BG.refuse_sample(name)
    with pytest.raises(SystemExit):
        BG.refuse_sample("other")
    for ok in ("big", "grad", "mid", "held", "held3", ""):
        BG.refuse_sample(ok)


def test_guard_refuses_any_path_that_mentions_held2(tmp_path: Path) -> None:
    for p in (
        tmp_path / "J_held2" / "run.json",
        tmp_path / "passes_held2" / "D4_J" / "pass.json",
        Path("a/b/pass_J_held2.log"),
        tmp_path / "stats_held2.json",
        tmp_path / "x" / "held2_clean.npy",
    ):
        with pytest.raises(SystemExit):
            BG.guarded(p)
    real = tmp_path / "passes_held2"
    real.mkdir()
    link = tmp_path / "link"
    link.symlink_to(real)
    with pytest.raises(SystemExit):
        BG.guarded(link / "x.json")  # a link into held2 is refused as well
    assert BG.guarded(tmp_path / "J_held3" / "run.json") == tmp_path / "J_held3" / "run.json"
    assert BG.guarded(tmp_path / "J_held" / "run.json").name == "run.json"


def test_other_random_draws_are_untouched_by_a_stop() -> None:
    """The meter draws no random number: a global generator is in the same state after a run
    with a stop as after a run without one."""
    state = np.random.get_state()
    _drive([0.9, 0.5] + [0.7] * 10, 5, lambda c: 0.96)
    a = np.random.random()
    np.random.set_state(state)
    _drive([0.9, 0.5] + [0.7] * 10, 0, lambda c: 0.96)
    b = np.random.random()
    assert a == b
    src = open(BR.__file__).read()
    assert "np.random" not in src and "default_rng" not in src
