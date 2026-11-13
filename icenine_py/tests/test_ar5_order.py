"""AR5 step 8 (task ar5h): the sweep rule on synthetic steps and the aligned-draw seeding.
Written before the scripts ``order_rules.py`` and ``order_draws.py`` (plan "AR5 step 8")."""

import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import mix_rules as MR  # noqa: E402
import order_draws as OD  # noqa: E402
import order_logs as OL  # noqa: E402
import order_rules as OR  # noqa: E402
import order_summary as OS  # noqa: E402


def _arm(hit: List[float], cost: List[float], evals: List[int], err: Any = None) -> Dict[str, Any]:
    err = err if err is not None else [0.1] * len(hit)
    return dict(fin_hit=hit, fin_cost=cost, fin_evals=evals, fin_err=err)


# ---- the sweep rule ------------------------------------------------------------------------
def test_first_cma_end_at_tau_stops_the_sweep_and_mc_never_runs() -> None:
    c = _arm([0.5, 0.95], [0.4, 0.2], [201, 201])
    m = _arm([0.99, 0.99], [0.1, 0.1], [1000, 1000])
    s = OR.sweep_step(20, c, m)
    assert s["taken"] == 1 and s["mc_ran"] == 0 and s["stopped"]
    assert s["evals"] == 20 + 402 and s["hit"] == 0.95 and s["end"] == "C"


def test_interleaved_order_pays_for_mc_on_the_discarded_candidate() -> None:
    c = _arm([0.5, 0.96], [0.4, 0.2], [201, 201])
    m = _arm([0.7], [0.3], [1000])
    a = MR.mixed_step(20, c, m)
    y = OR.sweep_step(20, c, m)
    assert a["evals"] == 20 + 201 + 1000 + 201 and a["taken"] == 1
    assert y["evals"] == 20 + 402 and y["taken"] == 1  # the same end, 1000 evaluations fewer


def test_no_cma_end_reaches_tau_then_mc_in_order_with_lower_cost_end() -> None:
    c = _arm([0.5, 0.6, 0.7], [0.4, 0.3, 0.5], [201, 201, 201])
    m = _arm([0.7, 0.96, 0.99], [0.35, 0.2, 0.1], [500, 600, 700])
    s = OR.sweep_step(20, c, m)
    assert s["taken"] == 1 and s["end"] == "M" and s["hit"] == 0.96 and s["stopped"]
    assert s["evals"] == 20 + 3 * 201 + 500 + 600 and s["mc_ran"] == 2


def test_stop_is_tested_on_the_kept_lower_cost_end() -> None:
    # MC end reaches tau but the CMA end has the lower cost: the kept end is the CMA end (hit 0.6),
    # the sweep goes on to the next candidate
    c = _arm([0.5, 0.6], [0.1, 0.3], [201, 201])
    m = _arm([0.97, 0.98], [0.2, 0.25], [500, 600])
    s = OR.sweep_step(0, c, m)
    assert s["taken"] == 1 and s["end"] == "M" and s["mc_ran"] == 2


def test_cost_tie_keeps_the_cma_end() -> None:
    c = _arm([0.5], [0.3], [201])
    m = _arm([0.97], [0.3], [500])
    s = OR.sweep_step(0, c, m)
    assert s["end"] == "C" and s["hit"] == 0.5 and not s["stopped"]


def test_nothing_reaches_tau_takes_the_lowest_cost_first_on_a_tie() -> None:
    c = _arm([0.5, 0.6, 0.4], [0.4, 0.3, 0.3], [201, 201, 201])
    m = _arm([0.6, 0.7, 0.8], [0.45, 0.31, 0.5], [500, 500, 500])
    s = OR.sweep_step(0, c, m)
    assert s["taken"] == 1 and s["end"] == "C" and not s["stopped"]  # costs .4 .3 .3 -> first .3


def test_missing_mc_end_is_reported() -> None:
    c = _arm([0.5, 0.6], [0.4, 0.3], [201, 201])
    m = _arm([0.97], [0.5], [500])  # kept end of candidate 0 is the CMA end; candidate 1 absent
    s = OR.sweep_step(5, c, m)
    assert s["missing"] and s["evals"] == 5


def test_step_without_candidates() -> None:
    s = OR.sweep_step(0, None, None)
    assert s["evals"] == 0 and s["taken"] is None and not s["passed"]


def test_seed_test_and_right_flags_follow_the_taken_end() -> None:
    c = _arm([0.5], [0.4], [201], err=[3.0])
    m = _arm([0.85], [0.1], [500], err=[0.2])
    s = OR.sweep_step(0, c, m)
    assert s["end"] == "M" and s["passed"] and s["right"] and s["err"] == 0.2


def test_sweep_choose_never_calls_mc_when_a_cma_end_passes() -> None:
    calls: List[str] = []

    def cma(p: int) -> Any:
        calls.append(f"c{p}")
        return (0.96 if p == 2 else 0.1, 1.0, 10)

    def mc(p: int) -> Any:
        calls.append(f"m{p}")
        return (0.0, 1.0, 100)

    r = OR.sweep_choose(5, cma, mc)
    assert calls == ["c0", "c1", "c2"] and r["taken"] == 2 and r["evals"] == 30


# ---- the aligned draws ---------------------------------------------------------------------
def test_same_voxel_kind_count_gives_the_same_draws_whatever_ran_before() -> None:
    b1, b2 = OD.DrawBook(), OD.DrawBook()
    for k in (5, 6, 7, 5):  # b1 has other history
        b1.rng(k, "grid")
        b1.rng(k, "_local_fit")
    x1 = [b1.rng(42, "_local_fit").random(4) for _ in range(3)]
    x2 = [b2.rng(42, "_local_fit").random(4) for _ in range(3)]
    for a, b in zip(x1, x2):
        assert np.array_equal(a, b)


def test_count_kind_and_voxel_change_the_draws() -> None:
    b = OD.DrawBook()
    first = b.rng(1, "k").random(3)
    second = b.rng(1, "k").random(3)  # count 1
    other_kind = b.rng(1, "j").random(3)
    other_voxel = b.rng(2, "k").random(3)
    assert not np.array_equal(first, second)
    assert not np.array_equal(first, other_kind)
    assert not np.array_equal(first, other_voxel)
    assert np.array_equal(first, OD.make_rng(1, "k", 0).random(3))
    assert np.array_equal(second, OD.make_rng(1, "k", 1).random(3))


def test_voxel_key_depends_on_values_only() -> None:
    a = np.arange(24, dtype=np.float64).reshape(8, 3)
    assert OD.voxel_key(a) == OD.voxel_key(a.copy())
    assert OD.voxel_key(a) != OD.voxel_key(a + 1e-9)
    assert OD.voxel_key(7) == 7


def test_wrappers_replace_only_the_generator_and_name_the_calling_function() -> None:
    seen: List[Dict[str, Any]] = []

    def inner(**kw: Any) -> str:
        seen.append(kw)
        return "ok"

    book = OD.DrawBook()
    lo = OD.aligned_local_optimization(book, inner)
    verts = np.ones((8, 3))

    def _local_fit() -> Any:
        return lo(voxel_vertices=verts, phase_index=0, initial_orientation=np.eye(3), rng="shared")

    assert _local_fit() == "ok" and _local_fit() == "ok"
    assert seen[0]["rng"] is not seen[1]["rng"] and seen[0]["rng"] != "shared"
    expect = OD.make_rng(OD.voxel_key(verts), "_local_fit", 0).random(3)
    assert np.array_equal(seen[0]["rng"].random(3), expect)
    assert seen[0]["phase_index"] == 0

    grid = OD.aligned_reconstruct_voxel(book, inner)
    grid(voxel_vertices=verts, phase_index=1, rng="shared")
    assert seen[-1]["rng"] != "shared" and seen[-1]["phase_index"] == 1
    assert book.counts[(OD.voxel_key(verts), "grid")] == 1

    calls: List[Any] = []

    def try_fit(hooks: Any, bfs: Any, v: int, R0: Any, rng: Any) -> Any:
        calls.append(rng)
        return v

    tf = OD.aligned_try_fit(book, try_fit)
    assert tf("h", "b", 9, None, "shared") == 9 and calls[0] != "shared"
    assert np.array_equal(calls[0].random(2), OD.make_rng(9, "pass", 0).random(2))


# ---- the log readers and the verdict (review fix of 2026-11-12; written after the runs) --------
def _fit(pos: int, kind: str, hit: float, cost: float, ev: int = 201) -> List[Any]:
    return [pos, pos, kind, hit, cost, ev, [float(pos)] * 9]


def _step(voxel: int, stop: int, fits: List[List[Any]]) -> Dict[str, Any]:
    return dict(voxel=voxel, stop=stop, fits=fits)


def test_kept_end_is_the_lower_cost_fit_and_a_tie_keeps_the_cma_end() -> None:
    e = _step(1, -1, [_fit(0, "C", 0.6, 0.30), _fit(0, "M", 0.7, 0.30, 900)])
    assert OL.kept_ends(e)[0]["kind"] == "C"  # tie: the CMA end stays
    e = _step(1, -1, [_fit(0, "C", 0.6, 0.30), _fit(0, "M", 0.7, 0.29, 900)])
    assert OL.kept_ends(e)[0]["kind"] == "M"


def test_taken_end_follows_stop_else_lowest_cost_first_on_a_tie() -> None:
    fits = [_fit(0, "C", 0.5, 0.4), _fit(1, "C", 0.96, 0.5), _fit(2, "C", 0.4, 0.2)]
    assert OL.taken_end(_step(1, 1, fits))["pos"] == 1  # stop index, not the lowest cost
    assert OL.taken_end(_step(1, -1, fits))["pos"] == 2
    tie = [_fit(0, "C", 0.5, 0.2), _fit(1, "C", 0.5, 0.2)]
    assert OL.taken_end(_step(1, -1, tie))["pos"] == 0


def test_log_counts_flags_mc_after_a_passing_cma_end_in_y_only() -> None:
    ok = _step(1, 1, [_fit(0, "C", 0.5, 0.4), _fit(1, "C", 0.96, 0.2)])
    bad = _step(2, 0, [_fit(0, "C", 0.96, 0.2), _fit(0, "M", 0.97, 0.3, 500)])
    y = OL.log_counts([ok, bad], "Y")
    assert y["bad_mc_after_cma_pass"] == 1 and y["mc_evals"] == 500 and y["steps_with_mc"] == 1
    a_ok = _step(
        3, 1, [_fit(0, "C", 0.5, 0.4), _fit(0, "M", 0.6, 0.3, 500), _fit(1, "C", 0.96, 0.2)]
    )
    assert OL.log_counts([a_ok], "A")["bad_mc_after_cma_pass"] == 0


def test_compare_logs_counts_matched_one_arm_and_differing_ends() -> None:
    la = [_step(5, 0, [_fit(0, "C", 0.96, 0.1)]), _step(6, 0, [_fit(0, "C", 0.96, 0.1)])]
    lb = [
        _step(5, 1, [_fit(0, "C", 0.5, 0.3), _fit(1, "C", 0.96, 0.2)]),
        _step(7, 0, [_fit(0, "C", 0.96, 0.1)]),
    ]
    c = OL.compare_logs(la, lb, lambda v, R: R[0] > 0.5)  # candidate 1's end is "wrong"
    assert (c["matched"], c["only_a"], c["only_b"], c["taken_differs"]) == (1, 1, 1, 1)
    assert (c["differing_wrong_first"], c["differing_wrong_second"]) == (0, 1)
    assert c["shared_fits_compared"] == 1 and c["shared_fits_different"] == 1
    assert OL.same_voxel_prefix(la, lb) == 1 and OL.first_taken_difference(la, lb) == 0


def test_flip_origin() -> None:
    assert OL.flip_origin(False, True, False, True) == "before the pass"
    assert OL.flip_origin(False, False, False, True) == "in the pass"
    assert OL.flip_origin(True, False, False, True) == "in the pass"  # reversed by the pass


def test_verdict_clause_judges_y_without_rule_5_only_where_a_fails_it() -> None:
    r = dict(r1=True, r2=True, r3=True, r4=True, r5=False)
    assert not OS.y_passes(r, a_fails_5=False)  # this step: A passes (5), Y fails it
    assert OS.y_passes(r, a_fails_5=True)
    assert not OS.y_passes(dict(r, r2=False, r5=True), a_fails_5=True)
    per = dict(held=True, mid=True, grad=True, big=True)
    assert OS.step_passes(per) and not OS.step_passes(dict(per, mid=False))


def test_table_headers_carry_no_bare_abbreviation_and_label_the_wall() -> None:
    t = OS.labelled([dict(sample="big", wall="9")], ["sample", "wall"])
    assert "contended, not a timing result" in t.splitlines()[0]
    t = OS.labelled([dict(grid="0 to 0")], ["grid"], table="y_vs_a")
    assert t.splitlines()[0] == "| grid fallbacks, A to Y |"
