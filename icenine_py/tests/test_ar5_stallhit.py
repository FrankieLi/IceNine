"""AR5 step 12 (task ar5l): rule PH(k), stop a stalled seed-step MC refinement only when its best
is already at hit ratio 0.95, on synthetic sequences. Plan "AR5 step 12", phase A step 2: stop at
0.95, go on under it, re-check after a gain, the first lowest cost kept, one evaluation per check,
the walk of a whole seed step, and that no other random draw moves.
"""

import inspect
import random
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import stallhit_rules as SR  # noqa: E402


class _Out:
    def __init__(self, cost: float) -> None:
        self.cost = cost


class _FakeCost:
    """A cost function whose cost is the first coordinate of the orientation."""

    def evaluate(self, orientation: Any, *args: Any, **kw: Any) -> _Out:
        return _Out(float(np.asarray(orientation).ravel()[0]))


def _drive(
    costs: List[float], k: int, hit: Callable[[float], float], tau: float = 0.95
) -> Dict[str, Any]:
    """Feed ``costs`` through the wrapped cost function; ``hit(best cost)`` plays the hit ratio of
    the running best. Returns the meter, whether the refinement was stopped and the hit calls."""
    meter = SR.Meter()
    calls: List[float] = []

    def hit_fn(R: Any) -> float:
        calls.append(float(np.asarray(R).ravel()[0]))
        return hit(float(np.asarray(R).ravel()[0]))

    meter.hit_fn = hit_fn
    patched = SR.patch_evaluate(_FakeCost.evaluate, meter, tau)
    fake = _FakeCost()
    meter.begin(k)
    stopped = False
    try:
        for c in costs:
            patched(fake, np.full((3, 3), c))
    except SR.StallHitStop:
        stopped = True
    return dict(stopped=stopped, meter=meter, calls=calls)


def test_constants_are_the_frozen_values() -> None:
    assert SR.K == 3000 and SR.TAU == 0.95 and SR.SEED_TEST == 0.80
    assert "PH(3000)" in SR.RULE and "0.95" in SR.RULE


def test_counter_fires_at_a_stall_and_is_disarmed_until_the_next_gain() -> None:
    c = SR.StallHitCounter(3)
    assert [c.feed(x) for x in (0.5, 0.6, 0.6, 0.6)] == [False, False, False, True]
    c.go_on()
    assert not any(c.feed(0.6) for _ in range(20))  # no gain, no second check
    assert c.feed(0.4) is False  # a gain re-arms it
    assert [c.feed(0.5) for _ in range(3)] == [False, False, True]


def test_counter_off_never_fires() -> None:
    c = SR.StallHitCounter(0)
    assert not any(c.feed(0.5) for _ in range(10_000))


def test_stops_at_a_stall_when_the_best_is_at_0_95() -> None:
    costs = [0.9, 0.5] + [0.7] * 10 + [0.4]
    r = _drive(costs, 5, lambda c: 0.96)
    m = r["meter"]
    # the last gain is evaluation 2; the check is at evaluation 7
    assert r["stopped"] and not m.on and m.n == 7
    assert m.checks == [[7, 0.96, True]] and m.continued == 0
    assert float(m.best_R[0, 0]) == 0.5


def test_stops_exactly_at_0_95_and_not_just_under() -> None:
    costs = [0.9, 0.5] + [0.7] * 10
    assert _drive(costs, 5, lambda c: 0.95)["stopped"]
    assert not _drive(costs, 5, lambda c: 0.9499999)["stopped"]


def test_goes_on_under_0_95_and_counts_the_check_as_an_evaluation() -> None:
    costs = [0.9, 0.5] + [0.7] * 10
    r = _drive(costs, 5, lambda c: 0.5)
    m = r["meter"]
    assert not r["stopped"] and m.on and m.n == len(costs)
    assert m.checks == [[7, 0.5, False]] and m.continued == 1
    assert len(r["calls"]) == 1  # one hit evaluation for the check; no second without a gain


def test_rechecks_after_a_gain_and_stops_when_the_new_best_is_at_0_95() -> None:
    costs = [0.9, 0.6] + [0.7] * 6 + [0.4] + [0.8] * 10  # stall at 7, gain at 9, stall at 14

    def hit(c: float) -> float:
        return 0.5 if c == 0.6 else 0.97

    r = _drive(costs, 5, hit)
    m = r["meter"]
    assert r["stopped"] and m.n == 14
    assert m.checks == [[7, 0.5, False], [14, 0.97, True]] and m.continued == 1
    assert float(m.best_R[0, 0]) == 0.4 and r["calls"] == [0.6, 0.4]


def test_first_lowest_cost_is_kept_on_a_tie() -> None:
    costs = [0.9, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5, 0.3, 0.3]
    r = _drive(costs, 3, lambda c: 1.0)
    m = r["meter"]
    assert r["stopped"] and m.best == 0.3 and m.trace == [[1, 0.9], [2, 0.3]]
    assert m.counter.best_n == 2  # the tie at evaluation 3 does not move it


def test_exactly_one_hit_evaluation_per_stall_check() -> None:
    costs = [0.9, 0.5] + [0.7] * 40
    r = _drive(costs, 5, lambda c: 0.1)
    assert len(r["calls"]) == len(r["meter"].checks) == 1
    costs = [0.9] + [0.8, 0.7, 0.6, 0.5, 0.4, 0.3] * 1 + [0.9] * 5 + [0.2] + [0.9] * 5
    r = _drive(costs, 5, lambda c: 0.1)
    assert len(r["calls"]) == len(r["meter"].checks) == 2


def test_a_gap_of_exactly_k_does_not_check() -> None:
    costs = [0.9, 0.5, 0.7, 0.7, 0.7, 0.7, 0.4]  # gain at 2, next gain 4 evaluations later
    r = _drive(costs, 5, lambda c: 1.0)
    assert not r["stopped"] and r["calls"] == []


def test_off_only_records_and_calls_nothing() -> None:
    r = _drive([0.9, 0.5, 0.7, 0.7, 0.4], 0, lambda c: 1.0)
    m = r["meter"]
    assert not r["stopped"] and m.n == 5 and r["calls"] == []
    assert m.trace == [[1, 0.9], [2, 0.5], [5, 0.4]]


def test_wrapped_evaluate_ignores_calls_while_off_and_copies_the_orientation() -> None:
    meter = SR.Meter()
    patched = SR.patch_evaluate(_FakeCost.evaluate, meter)
    out = patched(_FakeCost(), np.full((3, 3), 0.7))
    assert out.cost == 0.7 and meter.n == 0
    meter.begin(5)
    R = np.full((3, 3), 0.3)
    patched(_FakeCost(), R)
    R[:] = 9.0
    assert float(meter.best_R[0, 0]) == 0.3


def test_the_check_is_not_counted_by_the_wrapper_itself() -> None:
    """The hit evaluation runs with the wrapper off: a cost evaluation inside ``hit_fn`` is not
    counted and cannot recurse into another check."""
    meter = SR.Meter()
    patched = SR.patch_evaluate(_FakeCost.evaluate, meter)
    fake = _FakeCost()

    def hit_fn(R: Any) -> float:
        patched(fake, np.full((3, 3), 0.1))  # as an overlap evaluation that reached the wrapper
        return 0.2

    meter.hit_fn = hit_fn
    meter.begin(2)
    for c in (0.5, 0.6, 0.6):
        patched(fake, np.full((3, 3), c))
    assert meter.n == 3 and meter.continued == 1 and meter.on


def _trace_from_costs(costs: List[float]) -> List[List[float]]:
    out: List[List[float]] = []
    best = float("inf")
    for i, c in enumerate(costs):
        if c < best:
            best = c
            out.append([i + 1, c])
    return out


def test_stall_events_gap_rules() -> None:
    tr = [[1, 0.9], [4, 0.8], [20, 0.7]]
    ev = SR.stall_events(tr, 40, 5)
    assert [(e["i"], e["at"], e["last"]) for e in ev] == [(1, 9, False), (2, 25, True)]
    # the last check falls inside the refinement only if n_total > n + k
    assert [e["i"] for e in SR.stall_events(tr, 25, 5)] == [1]
    assert [e["i"] for e in SR.stall_events(tr, 26, 5)] == [1, 2]
    # a gap of exactly k does not stall; off or empty: no events
    assert SR.stall_events([[1, 0.9], [6, 0.8]], 6, 5) == []
    assert SR.stall_events(tr, 40, 0) == [] and SR.stall_events([], 40, 5) == []


def test_refine_walk_after_the_last_gain_keeps_the_end() -> None:
    tr = [[1, 0.9], [2, 0.5]]
    r = SR.refine_walk(tr, 100, 5, lambda i: 0.96)
    assert r["kind"] == "after_last" and r["evals"] == 8 and r["cost"] == 0.5 and r["stop_i"] == 1


def test_refine_walk_stops_before_a_later_gain_and_moves_the_end() -> None:
    tr = [[1, 0.9], [2, 0.6], [30, 0.4]]
    r = SR.refine_walk(tr, 40, 5, lambda i: 0.97 if i == 1 else 0.5)
    assert r["kind"] == "moved" and r["evals"] == 8 and r["cost"] == 0.6 and r["continued"] == 0


def test_refine_walk_goes_on_under_0_95_and_counts_the_checks() -> None:
    tr = [[1, 0.9], [2, 0.6], [30, 0.4]]
    r = SR.refine_walk(tr, 100, 5, lambda i: 0.5)
    assert r["kind"] == "none" and r["evals"] == 102 and r["cost"] == 0.4 and r["continued"] == 2
    # first stall under 0.95, the final one at 0.95: stops at the second
    r = SR.refine_walk(tr, 100, 5, lambda i: 0.96 if i == 2 else 0.5)
    assert r["kind"] == "after_last" and r["evals"] == 36 + 1 and r["continued"] == 1
    assert r["cost"] == 0.4


def test_refine_walk_unknown_stall_is_listed_and_read_by_argument() -> None:
    tr = [[1, 0.9], [2, 0.6], [30, 0.4]]
    h: Callable[[int], Optional[float]] = lambda i: None if i == 1 else 0.5  # noqa: E731
    go = SR.refine_walk(tr, 100, 5, h, unknown="go")
    stop = SR.refine_walk(tr, 100, 5, h, unknown="stop")
    assert go["unknown_events"] == [1] and go["kind"] == "none"
    assert stop["kind"] == "moved" and stop["evals"] == 8


def test_refine_walk_agrees_with_the_live_wrapper_on_random_sequences() -> None:
    rng = random.Random(7)
    for _ in range(300):
        n = rng.randint(2, 60)
        k = rng.randint(1, 6)
        costs = [rng.choice([0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]) for _ in range(n)]
        # the live refinement ends with one final evaluation of the kept orientation
        tr = _trace_from_costs(costs)
        hits = {c: rng.choice([0.2, 0.96]) for _, c in tr}
        live = _drive(costs, k, lambda c: hits[c])
        m = live["meter"]
        walk = SR.refine_walk(tr, n + 1, k, lambda i: hits[tr[i][1]])
        assert live["stopped"] == (walk["kind"] != "none")
        assert m.continued == walk["continued"]
        if live["stopped"]:
            # evaluations seen, the final evaluation of the kept orientation, the checks gone on
            assert m.n + 1 + m.continued == walk["evals"]
            assert float(m.best_R[0, 0]) == walk["cost"]
        else:
            assert n + 1 + m.continued == walk["evals"]


def _fits(*rows: Any) -> List[Dict[str, Any]]:
    return [dict(pos=p, j=p, kind=k, hit=h, cost=c, evals=e) for p, k, h, c, e in rows]


def test_walk_without_cuts_is_the_stored_step() -> None:
    fits = _fits((0, "C", 0.5, 0.3, 200), (0, "M", 0.96, 0.2, 10000))
    w = SR.walk_step(fits, {})
    assert w["evals"] == 10200 and w["taken"] == 0 and w["stopped"] and not w["lost"]


def test_walk_cut_after_last_gain_changes_only_the_count() -> None:
    fits = _fits((0, "C", 0.5, 0.3, 200), (0, "M", 0.96, 0.2, 10000))
    cuts = {0: dict(kind="after_last", evals=4000, cost=0.2, continued=0)}
    w = SR.walk_step(fits, cuts)
    assert w["evals"] == 4200 and w["taken"] == 0 and w["stopped"]
    assert w["kept"][0]["which"] == "M" and not w["kept"][0]["moved"]


def test_walk_moved_end_at_0_95_still_stops_the_step_at_the_same_candidate() -> None:
    fits = _fits(
        (0, "C", 0.5, 0.3, 200),
        (0, "M", 0.96, 0.2, 10000),
    )
    cuts = {0: dict(kind="moved", evals=3000, cost=0.22)}
    w = SR.walk_step(fits, cuts)
    assert w["stopped"] and w["taken"] == 0 and w["kept"][0]["moved"]
    assert w["kept"][0]["hit"] >= 0.95 and w["kept"][0]["cost"] == 0.22
    cuts = {0: dict(kind="moved", evals=3000, cost=0.22, hit=0.97)}
    assert SR.walk_step(fits, cuts)["kept"][0]["hit"] == 0.97


def test_walk_stops_at_the_first_kept_end_at_tau_and_skips_later_candidates() -> None:
    fits = _fits(
        (0, "C", 0.5, 0.3, 200),
        (0, "M", 0.6, 0.25, 5000),
        (1, "C", 0.5, 0.3, 200),
        (1, "M", 0.97, 0.1, 9000),
        (2, "C", 0.5, 0.3, 200),
    )
    w = SR.walk_step(fits, {})
    assert w["taken"] == 1 and w["stopped"] and w["evals"] == 200 + 5000 + 200 + 9000


def test_walk_without_a_tau_end_takes_the_lowest_cost() -> None:
    fits = _fits((0, "C", 0.5, 0.3, 200), (1, "C", 0.4, 0.2, 200), (2, "C", 0.3, 0.2, 200))
    w = SR.walk_step(fits, {})
    assert w["taken"] == 1 and not w["stopped"]


def test_walk_counts_the_checks_that_went_on_in_the_evaluations() -> None:
    fits = _fits((0, "C", 0.5, 0.3, 200), (0, "M", 0.96, 0.2, 10000))
    cuts = {0: dict(kind="none", evals=10002, cost=0.2, continued=2)}
    assert SR.walk_step(fits, cuts)["evals"] == 10202


def test_walk_a_lost_stop_runs_the_remaining_candidates() -> None:
    fits = _fits((0, "C", 0.5, 0.3, 200), (0, "M", 0.96, 0.2, 10000))
    # a (hypothetical) moved end that turns out dearer than the CMA end loses the stop
    cuts = {0: dict(kind="moved", evals=3000, cost=0.35, hit=0.97)}
    w = SR.walk_step(fits, cuts)
    assert w["lost"] and w["extra_candidates"] == 4 and w["extra_evals"] > 0


def test_the_rule_draws_no_random_number() -> None:
    src = inspect.getsource(SR)
    assert "random" not in src.replace("no random number", "").replace("random draw", "")
    state = random.getstate()
    nstate = np.random.get_state()[1].copy()
    c = SR.StallHitCounter(4)
    for x in (0.5, 0.6, 0.6, 0.6, 0.6):
        c.feed(x)
    c.go_on()
    SR.refine_walk([[1, 0.5]], 100, 4, lambda i: 1.0)
    _drive([0.9, 0.5] + [0.7] * 10, 5, lambda c: 0.96)
    assert random.getstate() == state
    assert (np.random.get_state()[1] == nstate).all()


def test_the_rule_agrees_with_step_11_cut_point_when_every_check_stops() -> None:
    """With every stall at 0.95 or more the rule stops where step 11's P(k) does."""
    import patience_rules as PR

    rng = random.Random(3)
    for _ in range(200):
        costs = [rng.choice([0.9, 0.8, 0.7, 0.6, 0.5, 0.4]) for _ in range(rng.randint(3, 60))]
        k = rng.randint(1, 6)
        tr = _trace_from_costs(costs)
        n_total = len(costs) + 1
        a = PR.cut_point(tr, n_total, k)
        b = SR.refine_walk(tr, n_total, k, lambda i: 1.0)
        if a["kind"] == "none":
            # P(k) never fires; PH can still have a check at the very last evaluation
            assert b["kind"] in ("none", "after_last") and b["evals"] <= n_total
        else:
            assert (a["kind"], a["cost"]) == (b["kind"], b["cost"]) and a["evals"] == b["evals"]


def test_freeze_lists_the_rule_scripts_and_the_library() -> None:
    import stallhit_freeze as SF

    assert SF.SCRIPTS[:3] == ["stallhit_rules.py", "stallhit_common.py", "stallhit_run.py"]
    assert SF.CONSTANTS["k"] == 3000 and SF.CONSTANTS["tau"] == 0.95
    cur = SF.current()
    assert "scripts/autoresearch/ar5/stallhit_run.py" in cur["files"]
    assert "scripts/autoresearch/ar3/held_make_sample.py" in cur["files"]
    assert "icenine/reconstructor.py" in cur["library"] and "uv.lock" in cur["library"]
    assert not any("held2" in k for k in list(cur["library"]))


def test_whatif_md_has_reading_lines(tmp_path: Path, monkeypatch: Any) -> None:
    import stallhit_whatif as SW

    run = dict(saved_share=0.07, saved_share_without_extra=0.075, taken_changed=7)
    stored = {
        "mid": dict(
            step11_measured=dict(total_cut_pct=10.6, taken_changed_n=20),
            high_go=run,
            high_stop=dict(saved_share=0.08, saved_share_without_extra=0.083, taken_changed=11),
            low=dict(saved_share=0.03, saved_share_without_extra=0.033, taken_changed=13),
            cases={c: dict(n=1, saved=10, share=0.001) for c in "ABCDE"},
            k_table={str(k): dict(exact_share=0.05, bound_changed=3) for k in (2000, 3000, 5000)},
        )
    }
    monkeypatch.setattr(SW, "SH", tmp_path)
    SW.write_md(stored, ["mid"])
    text = (tmp_path / "whatif.md").read_text()
    assert "Reading line, `stallhit_whatif_rules`" in text
    assert "| gated high: cut only if the best is at 0.95 or more | 7.5% to 8.3%; 7 to 11 |" in text
    assert "table:stallhit_whatif_cases" in text and "table:stallhit_whatif_k" in text


def test_gate_passes_a_clean_replay_and_fails_each_item() -> None:
    import stallhit_gate as SG

    def rows(n: int, changed: int, wrong: bool = False) -> List[Dict[str, Any]]:
        out = []
        for i in range(n):
            ch = i < changed
            out.append(
                dict(
                    step=i,
                    classes=["same", "moved"] if ch else ["same", "after_last"],
                    pred_mismatch=0,
                    bad=False,
                    path="seed" if ch else "same",
                    path_same=not ch,
                    g_total=3000,
                    new_total=1900,
                    g_stopped=True,
                    new_hit=0.97,
                    g_hit=0.97,
                    g_failed=False,
                    new_failed=False,
                    err_g=0.02,
                    err_new=2.0 if (ch and wrong) else 0.02,
                    taken_same=not ch,
                    mis_to_g_deg=0.05 if ch else 0.0,
                )
            )
        return out

    n = SG.STEPS_WITH_MC["mid"]
    ok = SG.sample_gate("mid", rows(n, 8))
    assert ok["g1_ok"] and ok["g2_ok"] and ok["g3_ok"] and ok["g4_ok"]
    assert not SG.sample_gate("mid", rows(n, 13))["g2_ok"]  # more than 12 changed ends
    assert not SG.sample_gate("mid", rows(n, 8, wrong=True))["g3_ok"]
    assert not SG.sample_gate("mid", rows(n - 1, 8))["g1_ok"]  # a step is missing
    low = rows(n, 8)
    for r in low:
        r["new_total"] = 3000
    assert not SG.sample_gate("mid", low)["g4_ok"]  # saves under 6%
    bad = rows(n, 8)
    bad[0]["pred_mismatch"] = 1
    bad[0]["classes"] = ["same", "moved"]
    assert not SG.sample_gate("mid", bad)["g1_ok"]
    assert SG.MAX_CHANGED == {"mid": 12, "grad": 1, "big": 3} and SG.G4_MID == 0.06


def test_run_all_dry_run_arguments_and_job_order() -> None:
    text = (ROOT / "scripts" / "autoresearch" / "ar5" / "stallhit_run_all.sh").read_text()
    for needle in ("--resume", "--dry-run", "stallhit_run_all.done", "exit=$code"):
        assert needle in text
    assert text.index("held2 build") < text.index("S8 held2")
