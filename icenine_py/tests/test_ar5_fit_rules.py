"""Pure rule functions of AR5 step 2 (``scripts/autoresearch/ar5/fit_rules.py``): the early stop of
the seed step, the paired split, the pick rule and mid's pass rule of Part 1, and the pass rule of
run S10. Written before the module (plan "AR5 step 2")."""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import fit_rules as FR  # noqa: E402


# ---- the early stop of the seed step --------------------------------------------------------
def test_early_stop_takes_the_first_end_at_or_above_tau() -> None:
    hits = [0.50, 0.95, 0.99]
    costs = [0.30, 0.20, 0.10]
    # the third end is better but is never refined: the loop stopped at the second
    assert FR.early_stop(hits, costs, 0.95) == (1, True)


def test_early_stop_without_a_pass_takes_the_lowest_cost() -> None:
    assert FR.early_stop([0.5, 0.6, 0.7], [0.30, 0.10, 0.20], 0.95) == (1, False)


def test_early_stop_cost_tie_takes_the_earlier_end() -> None:
    assert FR.early_stop([0.5, 0.6], [0.2, 0.2], 0.95) == (0, False)


def test_stop_position_is_none_without_a_pass() -> None:
    assert FR.stop_position([0.1, 0.94], 0.95) is None
    assert FR.stop_position([0.1, 0.95], 0.95) == 1


def test_seed_test_edge() -> None:
    assert FR.seed_test(0.80)
    assert not FR.seed_test(0.7999)


# ---- paired split ---------------------------------------------------------------------------
def test_paired_split_counts() -> None:
    a = np.array([True, True, False, False, True])
    b = np.array([True, False, True, False, False])
    assert FR.paired_split(a, b) == dict(both=1, only_a=2, only_b=1, neither=1)


# ---- Part 1: pick rule and mid's pass rule --------------------------------------------------
def _counts(n_pass: int, n_right: int, evals: int) -> dict:
    return dict(n_pass=n_pass, n_right=n_right, evals=evals)


def test_pick_takes_the_cheapest_arm_that_loses_at_most_one_percent() -> None:
    n = 1000  # 1% = 10 steps
    arms = {
        "M": _counts(990, 985, 1_000_000),
        "C60": _counts(985, 980, 50_000),  # 5 fewer each: ok
        "C200": _counts(989, 985, 120_000),
        "C0": _counts(900, 900, 0),  # 90 fewer: fails (a)
    }
    assert FR.pick_arm(arms, n) == "C60"


def test_pick_exactly_one_percent_fewer_passes() -> None:
    n = 1000
    arms = {"M": _counts(990, 985, 100), "C60": _counts(980, 985, 10)}
    assert FR.pick_arm(arms, n) == "C60"  # exactly 10 fewer is "at most 1%"
    arms["C60"] = _counts(979, 985, 10)
    assert FR.pick_arm(arms, n) == "M"  # 11 fewer: no cheaper arm qualifies


def test_pick_fails_on_rule_b_alone() -> None:
    n = 1000
    arms = {"M": _counts(990, 985, 100), "C60": _counts(990, 960, 10)}
    assert FR.pick_arm(arms, n) == "M"


def test_pick_more_passes_than_m_is_fine() -> None:
    n = 100
    arms = {"M": _counts(90, 90, 100), "C60": _counts(95, 94, 10)}
    assert FR.pick_arm(arms, n) == "C60"


def test_mid_pass_needs_a_b_and_an_eighty_percent_fall() -> None:
    ref = _counts(990, 985, 1_000_000)
    ok = FR.mid_pass(_counts(985, 980, 200_000), ref, 1000)
    assert ok["a"] and ok["b"] and ok["c"] and ok["passes"]
    assert ok["evals_fall"] == 0.8
    short = FR.mid_pass(_counts(985, 980, 200_001), ref, 1000)
    assert not short["c"] and not short["passes"]
    worse = FR.mid_pass(_counts(975, 980, 1), ref, 1000)
    assert not worse["a"] and not worse["passes"]


# ---- Part 2: the pass rule of S10 -----------------------------------------------------------
def test_s10_pass_rule_edges() -> None:
    r = FR.s10_pass(
        evals=8_000_000, evals_s8=10_000_000, broken=7, fixed=2, med=0.0300, med_s8=0.0250
    )
    assert r["fewer"] == 0.2 and r["a"] and r["b"] and r["c"] and r["passes"]
    r = FR.s10_pass(8_000_001, 10_000_000, 7, 2, 0.03, 0.025)
    assert not r["a"] and not r["passes"]
    r = FR.s10_pass(5_000_000, 10_000_000, 8, 2, 0.03, 0.025)  # broken - fixed = 6
    assert not r["b"] and not r["passes"]
    r = FR.s10_pass(5_000_000, 10_000_000, 2, 7, 0.03, 0.025)  # net fixed: fine
    assert r["b"] and r["passes"]
    r = FR.s10_pass(5_000_000, 10_000_000, 0, 0, 0.0301, 0.0250)  # 0.0051 above
    assert not r["c"] and not r["passes"]


def test_share_and_fall() -> None:
    assert FR.fall(25, 100) == 0.75
    assert FR.fall(0, 0) == 0.0
