"""AR5 step 11 (task ar5k): rule P(k), patience of the seed step's MC refinement, on synthetic
traces. Written before the scripts ``patience_*.py`` (plan "AR5 step 11", phase A step 2): the stop
point, the first lowest cost kept, the final evaluation, the walk of a whole seed step, and that no
other random draw moves.
"""

import inspect
import random
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import patience_rules as PR  # noqa: E402


def _trace_from_costs(costs: List[float]) -> List[List[float]]:
    out: List[List[float]] = []
    best = float("inf")
    for i, c in enumerate(costs):
        if c < best:
            best = c
            out.append([i + 1, c])
    return out


def test_constants_are_the_frozen_values() -> None:
    assert PR.K == 3000 and PR.TAU == 0.95 and PR.SEED_TEST == 0.80


def test_counter_fires_after_k_evaluations_without_a_gain() -> None:
    c = PR.PatienceCounter(3)
    assert [c.feed(x) for x in (0.5, 0.4, 0.4, 0.45, 0.41)] == [False, False, False, False, True]
    assert c.n == 5 and c.best == 0.4 and c.best_n == 2  # the tie at 3 does not move it


def test_counter_off_never_fires() -> None:
    c = PR.PatienceCounter(0)
    assert not any(c.feed(0.5) for _ in range(10_000))


def test_first_lowest_cost_is_kept_on_a_tie() -> None:
    costs = [0.9, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3]  # final evaluation: 0.3
    r = PR.run_costs(costs, 3)
    assert r["stopped"] and r["best_n"] == 2 and r["cost"] == 0.3


def test_stop_point_and_final_evaluation() -> None:
    costs = [0.9, 0.5] + [0.7] * 20 + [0.5]  # the last one is the final evaluation of 0.5
    r = PR.run_costs(costs, 5)
    # last gain at evaluation 2; the counter fires at evaluation 7; one final evaluation: 8
    assert r["stopped"] and r["best_n"] == 2 and r["evals"] == 8 and r["cost"] == 0.5


def test_a_run_that_ends_within_k_is_not_stopped() -> None:
    costs = [0.9, 0.5, 0.7, 0.7, 0.5]
    r = PR.run_costs(costs, 5)
    assert not r["stopped"] and r["evals"] == len(costs)


def test_a_gain_resets_the_count() -> None:
    costs = [0.9] + [0.95] * 4 + [0.8] + [0.95] * 4 + [0.7] + [0.95] * 4 + [0.7]
    r = PR.run_costs(costs, 5)
    assert not r["stopped"] and r["cost"] == 0.7


def test_cut_point_after_the_last_improvement_keeps_the_end() -> None:
    trace = [[1, 0.5], [10, 0.3], [40, 0.2]]
    cp = PR.cut_point(trace, 5000, 1000)
    assert cp["kind"] == "after_last" and cp["cost"] == 0.2 and cp["evals"] == 40 + 1000 + 1


def test_cut_point_before_a_later_improvement_moves_the_end() -> None:
    trace = [[1, 0.5], [10, 0.3], [2000, 0.2]]
    cp = PR.cut_point(trace, 2100, 1000)
    assert cp["kind"] == "moved" and cp["cost"] == 0.3 and cp["evals"] == 10 + 1000 + 1


def test_cut_point_gap_of_exactly_k_does_not_fire() -> None:
    # the gain at evaluation 1010 comes at the k-th evaluation after 10: the counter has not fired
    trace = [[1, 0.5], [10, 0.3], [1010, 0.2]]
    cp = PR.cut_point(trace, 1011, 1000)
    assert cp["kind"] == "none" and cp["evals"] == 1011
    trace2 = [[1, 0.5], [10, 0.3], [1011, 0.2]]
    assert PR.cut_point(trace2, 1012, 1000)["kind"] == "moved"


def test_cut_point_never_costs_more_than_the_uncut_run() -> None:
    trace = [[1, 0.5], [10, 0.3]]
    # last gain 10, final evaluation 1011: the counter fires at 1010, the final one is 1011
    assert PR.cut_point(trace, 1011, 1000)["kind"] == "none"
    assert PR.cut_point(trace, 1012, 1000) == dict(kind="after_last", evals=1011, cost=0.3)
    cp = PR.cut_point(trace, 1500, 1000)
    assert cp["evals"] == 1011 and cp["kind"] == "after_last"


def test_cut_point_gate_stops_only_while_worse_than_the_reference() -> None:
    trace = [[1, 0.5], [10, 0.3], [2000, 0.2]]
    assert PR.cut_point(trace, 2100, 1000, worse_than=0.35)["kind"] == "none"
    assert PR.cut_point(trace, 2100, 1000, worse_than=0.25)["kind"] == "moved"


def test_cut_point_off_or_short_is_untouched() -> None:
    trace = [[1, 0.5], [10, 0.3]]
    assert PR.cut_point(trace, 500, 0)["kind"] == "none"
    assert PR.cut_point(trace, 500, 1000)["evals"] == 500


def test_cut_point_agrees_with_the_counter_on_random_sequences() -> None:
    rng = np.random.default_rng(1)
    for _ in range(300):
        n = int(rng.integers(5, 400))
        costs = [float(x) for x in np.round(np.cumsum(rng.normal(0.0, 0.05, n)) + 1.0, 2)]
        costs.append(min(costs))  # the final evaluation of the kept orientation
        k = int(rng.integers(2, 60))
        r = PR.run_costs(costs, k)
        cp = PR.cut_point(_trace_from_costs(costs), len(costs), k)
        assert cp["evals"] == r["evals"]
        assert (cp["kind"] != "none") == r["stopped"]
        assert abs(cp["cost"] - r["cost"]) < 1e-12


def _fits(cma_hit: List[float], mc: Dict[int, Dict[str, float]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for p, h in enumerate(cma_hit):
        out.append(dict(pos=p, kind="C", hit=h, cost=0.30 - 0.01 * p, evals=201))
        if p in mc:
            out.append(
                dict(
                    pos=p,
                    kind="M",
                    evals=int(mc[p]["evals"]),
                    **{k: mc[p][k] for k in ("hit", "cost")},
                )
            )
    return out


def test_walk_without_cuts_is_the_stored_step() -> None:
    fits = _fits([0.5, 0.97], {0: dict(hit=0.6, cost=0.2, evals=4000)})
    w = PR.walk_step(fits, {})
    assert w["evals"] == 201 + 4000 + 201 and w["taken"] == 1 and w["stopped"] and not w["lost"]


def test_walk_cut_after_last_improvement_changes_only_the_count() -> None:
    fits = _fits([0.5, 0.97], {0: dict(hit=0.6, cost=0.2, evals=4000)})
    cut = dict(kind="after_last", evals=1500, cost=0.2)
    w = PR.walk_step(fits, {0: cut})
    assert w["evals"] == 201 + 1500 + 201 and w["taken"] == 1 and w["kept"][0]["cost"] == 0.2


def test_walk_moved_end_that_reached_tau_is_lost_when_it_turns_dearer_than_the_cma_end() -> None:
    fits = _fits([0.5, 0.97, 0.9], {0: dict(hit=0.96, cost=0.1, evals=4000)})
    # the stored step stops at candidate 0 (MC end hit 0.96); the cut end is dearer than the CMA end
    assert PR.walk_step(fits, {})["taken"] == 0
    cut = dict(kind="moved", evals=1500, cost=0.35)  # CMA end cost is 0.30
    w = PR.walk_step(fits, {0: cut})
    assert w["kept"][0]["which"] == "C" and w["kept"][0]["hit"] == 0.5
    assert w["taken"] == 1 and w["stopped"]


def test_walk_pessimistic_loses_a_moved_end_that_reached_tau() -> None:
    fits = _fits([0.5], {0: dict(hit=0.96, cost=0.1, evals=4000)})
    cut = dict(kind="moved", evals=1500, cost=0.15)  # still below the CMA end (0.30)
    opt = PR.walk_step(fits, {0: cut})
    assert opt["stopped"] and opt["taken"] == 0 and not opt["lost"]
    pes = PR.walk_step(fits, {0: cut}, pessimistic=True, n_avail=5)
    assert not pes["stopped"] and pes["lost"] and pes["extra_candidates"] == 4


def test_walk_without_a_tau_end_takes_the_lowest_cost() -> None:
    fits = _fits(
        [0.5, 0.6],
        {0: dict(hit=0.7, cost=0.25, evals=3000), 1: dict(hit=0.7, cost=0.22, evals=3000)},
    )
    w = PR.walk_step(fits, {})
    assert not w["stopped"] and w["taken"] == 1


def test_the_rule_draws_no_random_number() -> None:
    src = inspect.getsource(PR)
    assert "random" not in src.replace("no random number", "").replace("random draw", "")
    state = random.getstate()
    nstate = np.random.get_state()[1].copy()
    c = PR.PatienceCounter(4)
    for x in (0.5, 0.6, 0.6, 0.6, 0.6):
        c.feed(x)
    PR.cut_point([[1, 0.5]], 100, 4)
    assert random.getstate() == state
    assert (np.random.get_state()[1] == nstate).all()
