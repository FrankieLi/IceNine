"""Pure rule functions of AR5 step 3 (``scripts/autoresearch/ar5/mix_rules.py``): the mixed seed
refinement of one candidate and of one seed step, built from the stored CMA and MC ends of the
replay. Written before the module (plan "AR5 step 3")."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import fit_rules as FR  # noqa: E402
import mix_rules as MR  # noqa: E402


def _arm(hit, cost, evals, err=None):
    err = err if err is not None else [0.1] * len(hit)
    return dict(fin_hit=hit, fin_cost=cost, fin_evals=evals, fin_err=err)


# ---- one candidate --------------------------------------------------------------------------
def test_cma_end_at_tau_is_kept_and_mc_does_not_run() -> None:
    r = MR.mixed_candidate(0.95, 0.20, 100, 0.99, 0.10, 500)
    assert r["end"] == "C" and r["mc_ran"] is False and r["evals"] == 100


def test_cma_end_below_tau_runs_mc_and_counts_both() -> None:
    r = MR.mixed_candidate(0.90, 0.20, 100, 0.99, 0.30, 500)
    assert r["mc_ran"] is True and r["evals"] == 600


def test_lower_cost_end_is_kept_mc_when_strictly_lower() -> None:
    r = MR.mixed_candidate(0.90, 0.20, 100, 0.96, 0.10, 500)
    assert r["end"] == "M" and r["hit"] == 0.96 and r["cost"] == 0.10


def test_cost_tie_keeps_the_cma_end() -> None:
    r = MR.mixed_candidate(0.90, 0.20, 100, 0.96, 0.20, 500)
    assert r["end"] == "C" and r["hit"] == 0.90


def test_the_stop_test_is_on_the_lower_cost_end() -> None:
    # MC has hit 0.96 but the CMA end has the lower cost and hit 0.94: the step does not stop
    r = MR.mixed_candidate(0.94, 0.10, 100, 0.96, 0.20, 500)
    assert r["end"] == "C" and r["hit"] < FR.TAU


# ---- one seed step --------------------------------------------------------------------------
def test_step_stops_at_first_cma_end_with_tau() -> None:
    c = _arm([0.97, 0.99], [0.1, 0.1], [50, 50])
    m = _arm([0.5, 0.5], [0.2, 0.2], [400, 400])
    s = MR.mixed_step(20, c, m)
    assert s["evals"] == 20 + 50 and s["taken"] == 0 and s["stopped"] and s["mc_ran"] == 0


def test_step_falls_to_mc_then_stops_on_its_end() -> None:
    c = _arm([0.50, 0.99], [0.30, 0.1], [50, 50])
    m = _arm([0.96, 0.5], [0.20, 0.2], [400, 400])
    s = MR.mixed_step(20, c, m)
    assert s["evals"] == 20 + 50 + 400 and s["taken"] == 0 and s["hit"] == 0.96
    assert s["mc_ran"] == 1 and s["stopped"]


def test_step_without_a_pass_takes_the_lowest_cost_end_of_all() -> None:
    c = _arm([0.5, 0.6, 0.7], [0.30, 0.25, 0.40], [10, 10, 10], [1.0, 2.0, 3.0])
    m = _arm([0.4, 0.8, 0.6], [0.35, 0.22, 0.45], [100, 100, 100], [4.0, 5.0, 6.0])
    s = MR.mixed_step(20, c, m)
    assert not s["stopped"] and s["taken"] == 1 and s["err"] == 5.0
    assert s["evals"] == 20 + 3 * 110 and s["mc_ran"] == 3


def test_step_missing_mc_refinement_is_flagged_not_invented() -> None:
    c = _arm([0.5, 0.6], [0.3, 0.2], [10, 10])
    m = _arm([0.4], [0.35], [100])  # MC stopped after the first candidate in the replay
    s = MR.mixed_step(20, c, m)
    assert s["missing"] is True and s["evals"] == 20 + 110 + 10


def test_step_without_candidates_costs_nothing() -> None:
    s = MR.mixed_step(0, None, None)
    assert s["evals"] == 0 and s["taken"] is None and not s["passed"] and not s["right"]


def test_passed_and_right_flags_follow_the_seed_test_and_one_degree() -> None:
    c = _arm([0.85], [0.1], [10], [0.5])
    s = MR.mixed_step(5, c, _arm([0.0], [9.0], [1]))
    assert s["passed"] and s["right"]
    c = _arm([0.79], [0.1], [10], [1.5])
    s = MR.mixed_step(5, c, _arm([0.0], [9.0], [1]))
    assert not s["passed"] and not s["right"]


# ---- the pick uses step 2's rule ------------------------------------------------------------
def test_pick_takes_the_cheaper_arm_that_meets_both_rules() -> None:
    arms = {
        "M": dict(n_pass=1592, n_right=1583, evals=5_865_469),
        "X60": dict(n_pass=1594, n_right=1584, evals=4_390_198),
        "X200": dict(n_pass=1596, n_right=1589, evals=2_501_991),
    }
    assert FR.pick_arm(arms, 1608) == "X200"
    arms["X200"]["n_pass"] = 1570  # 22 below M: more than 1% of 1,608 steps
    assert FR.pick_arm(arms, 1608) == "X60"
