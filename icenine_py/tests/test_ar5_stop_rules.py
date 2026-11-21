"""Pure rule functions of AR5 step 1 (``scripts/autoresearch/ar5/stop_rules.py``): the try test of
the boundary pass (P1, P2, P3), the choice of a try's end, the pick, the hit-ratio flag, the ladder
accounting and the cut of the vote count. Written before the module (plan "AR5 step 1")."""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import stop_rules as SR  # noqa: E402


# ---- Part 1: the try test -----------------------------------------------------------------
def test_p1_keeps_equal_and_above_skips_lower() -> None:
    start = np.array([0.90, 0.95, 0.89, 0.50])
    old = np.array([0.90, 0.90, 0.90, 0.90])
    assert SR.p1_keep(start, old, 0.0).tolist() == [True, True, False, False]


def test_p1_tie_at_the_edge_is_kept() -> None:
    # 0.9 - 0.02 is 0.88 up to rounding; a start at exactly 0.88 is not "more than d under"
    assert SR.p1_keep(np.array([0.88]), np.array([0.90]), 0.02).tolist() == [True]
    assert SR.p1_keep(np.array([0.87]), np.array([0.90]), 0.02).tolist() == [False]
    assert SR.p1_keep(np.array([0.70]), np.array([0.90]), 0.20).tolist() == [True]


def test_p2_keeps_m_highest_start_hit_per_voxel() -> None:
    voxel = np.array([5, 5, 5, 9, 9, 12])
    hit = np.array([0.2, 0.9, 0.5, 0.1, 0.1, 0.0])
    cost = np.array([0.5, 0.5, 0.5, 0.7, 0.6, 0.9])
    assert SR.p2_keep(voxel, hit, cost, 1).tolist() == [False, True, False, False, True, True]
    # m = 2: voxel 9 has a tie in hit ratio; the lower start cost goes first, both are kept anyway
    assert SR.p2_keep(voxel, hit, cost, 2).tolist() == [False, True, True, True, True, True]


def test_p2_tie_goes_to_lower_cost_then_first_in_order() -> None:
    voxel = np.array([1, 1, 1])
    hit = np.array([0.5, 0.5, 0.5])
    cost = np.array([0.4, 0.3, 0.3])
    assert SR.p2_keep(voxel, hit, cost, 1).tolist() == [False, True, False]


def test_p3_is_p1_and_p2() -> None:
    voxel = np.array([1, 1, 1, 2, 2])
    start = np.array([0.9, 0.8, 0.1, 0.3, 0.2])
    cost = np.zeros(5)
    old = np.array([0.9, 0.9, 0.9, 0.9, 0.9])
    k = SR.rule_keep("P3", voxel, start, cost, old, d=0.05, m=2)
    assert k.tolist() == [True, False, False, False, False]
    assert SR.rule_keep("P1", voxel, start, cost, old, d=0.20, m=None).tolist() == [
        True,
        True,
        False,
        False,
        False,
    ]


def test_choose_try_rule_hit_old_wins_ties() -> None:
    # voxel 3 has old hit 0.9 / cost 0.1; voxel 7 has old hit 0.5 / cost 0.4
    old_hit = np.zeros(10)
    old_cost = np.ones(10)
    old_hit[3], old_cost[3] = 0.9, 0.1
    old_hit[7], old_cost[7] = 0.5, 0.4
    voxel = np.array([3, 3, 7, 7, 7])
    t_hit = np.array([0.9, 0.8, 0.5, 0.7, 0.7])
    t_cost = np.array([0.1, 0.1, 0.4, 0.3, 0.2])
    keep = np.ones(5, bool)
    ch = SR.choose_try(voxel, old_hit, old_cost, t_hit, t_cost, keep)
    assert ch[3] == -1  # a tie with the old end: the old end stays
    assert ch[7] == 4  # equal hit ratio: the lower cost; cost ties would go to the first
    assert ch[0] == -1  # a voxel with no try


def test_choose_try_ignores_skipped_tries() -> None:
    old_hit, old_cost = np.zeros(4), np.ones(4)
    voxel = np.array([2, 2])
    ch = SR.choose_try(
        voxel,
        old_hit,
        old_cost,
        np.array([0.9, 0.6]),
        np.array([0.1, 0.2]),
        np.array([False, True]),
    )
    assert ch[2] == 1


def test_pass_evals_counts_every_start_and_the_kept_tries() -> None:
    evals = np.array([61, 40, 61, 61])
    keep = np.array([True, False, False, True])
    assert SR.pass_evals(evals, keep) == 4 + 61 + 61
    assert SR.pass_evals(evals, np.ones(4, bool)) == 4 + evals.sum()


def test_counts_paired_against_full_pass() -> None:
    wrong_before = np.array([1, 1, 1, 0, 0, 1], bool)
    wrong_full = np.array([0, 0, 1, 0, 1, 1], bool)
    wrong_rule = np.array([0, 1, 1, 0, 0, 1], bool)
    c = SR.paired_vs_full(wrong_before, wrong_full, wrong_rule)
    assert (c["fixed"], c["broken"], c["both_wrong"]) == (1, 1, 2)
    assert SR.net_gain(wrong_before, wrong_full) == 4 - 3
    assert SR.net_gain(wrong_before, wrong_rule) == 4 - 3


def test_pick_setting_fewest_evals_among_qualifying() -> None:
    rows = [
        dict(name="a", evals=100, net_gain=50),
        dict(name="b", evals=40, net_gain=48),  # loses 4%: qualifies
        dict(name="c", evals=20, net_gain=46),  # loses 8%: does not
        dict(name="d", evals=40, net_gain=50),  # ties b on evals: b is earlier
    ]
    assert SR.pick_setting(rows, full_net_gain=50, max_loss=0.05)["name"] == "b"
    assert SR.pick_setting(rows[2:3], full_net_gain=50, max_loss=0.05) is None


def test_pick_d_falls_back_to_least_loss() -> None:
    rows = [
        dict(name="d=0", d=0.0, evals=10, net_gain=20),
        dict(name="d=0.2", d=0.2, evals=90, net_gain=40),
    ]
    assert SR.pick_d(rows, full_net_gain=50, max_loss=0.05)["d"] == 0.2
    rows.append(dict(name="d=0.1", d=0.1, evals=30, net_gain=48))
    assert SR.pick_d(rows, full_net_gain=50, max_loss=0.05)["d"] == 0.1


def test_pass_rule_on_mid() -> None:
    assert SR.try_test_passes(evals_rule=49, evals_full=100, gain_rule=95, gain_full=100)
    assert not SR.try_test_passes(evals_rule=51, evals_full=100, gain_rule=100, gain_full=100)
    assert not SR.try_test_passes(evals_rule=10, evals_full=100, gain_rule=94, gain_full=100)


# ---- Part 2: the ladder -------------------------------------------------------------------
def test_ladder_total_stops_at_first_rung_that_passes() -> None:
    # R1 passes: only R1 is paid
    t = SR.ladder_case(
        r1=dict(evals=400, hit=0.97), r3=dict(evals=30000, hit=0.9), grid_evals=333000
    )
    assert t["evals"] == 400 and t["rung"] == "R1"
    # R1 fails (hit < 0.80), R3 passes
    t = SR.ladder_case(
        r1=dict(evals=400, hit=0.5), r3=dict(evals=30000, hit=0.9), grid_evals=333000
    )
    assert t["evals"] == 30400 and t["rung"] == "R3"
    # both fail: the full grid search is paid on top, its level 0 counted again
    t = SR.ladder_case(
        r1=dict(evals=400, hit=0.5), r3=dict(evals=30000, hit=0.7), grid_evals=333000
    )
    assert t["evals"] == 400 + 30000 + 333000 and t["rung"] == "grid"


def test_ladder_verdict() -> None:
    # at most 2/3 of the evaluations and at most one case fewer right than the grid search
    v = SR.ladder_verdict(ladder_evals=200, grid_evals=300, right_ladder=14, right_grid=15)
    assert v["passes"] and v["evals_ok"] and v["right_ok"]
    v = SR.ladder_verdict(ladder_evals=201, grid_evals=300, right_ladder=15, right_grid=15)
    assert not v["passes"] and not v["evals_ok"]
    v = SR.ladder_verdict(ladder_evals=100, grid_evals=300, right_ladder=13, right_grid=15)
    assert not v["passes"] and not v["right_ok"]


# ---- Part 3: the flag ---------------------------------------------------------------------
def test_flag_counts() -> None:
    hit = np.array([0.40, 0.60, 0.90, 0.99, 0.50])
    wrong = np.array([True, True, False, True, False])
    c = SR.flag_counts(hit, wrong, 0.60)
    # under 0.60: 0.40 and 0.50 (0.60 itself is not under)
    assert (c["flagged"], c["wrong_flagged"], c["wrong_unflagged"]) == (2, 1, 2)


def test_pick_threshold_lowest_with_enough_recall() -> None:
    hit = [np.array([0.5, 0.7, 0.9, 0.99]), np.array([0.6, 0.95])]
    wrong = [np.array([True, True, True, False]), np.array([True, False])]
    thr = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    # 4 wrong in all; recall 0.8 needs 4 of 4 -> 3.2 -> at least 4: needs hit < thr for 0.9 -> 0.95
    assert SR.pick_threshold(hit, wrong, thr, recall=0.8) == 0.95
    assert SR.pick_threshold(hit, wrong, thr, recall=0.5) == 0.7
    assert SR.pick_threshold(hit, wrong, [0.5], recall=0.8) is None


def test_hit_after_uses_the_taken_try() -> None:
    old_hit = np.array([0.1, 0.2, 0.3])
    voxel = np.array([1, 1, 2])
    t_hit = np.array([0.5, 0.6, 0.1])
    choice = np.array([-1, 1, -1])
    out = SR.hit_after(old_hit, voxel, t_hit, choice)
    assert out.tolist() == [0.1, 0.6, 0.3]


def test_every_20th_unflagged_boundary_voxel() -> None:
    idx = np.arange(100)
    flagged = np.zeros(100, bool)
    flagged[[0, 20]] = True
    boundary = np.ones(100, bool)
    boundary[40] = False
    sel = SR.sample_unflagged(flagged, boundary, every=20)
    # unflagged boundary voxels in index order: the 20th, 40th ... of that list
    un = [i for i in idx if not flagged[i] and boundary[i]]
    assert sel.tolist() == un[19::20]


# ---- Part 3: the vote count cut -------------------------------------------------------------
def test_best_cut_and_overlap() -> None:
    wrong = np.array([5, 6, 7, 2])
    right = np.array([1, 2, 3, 4, 1, 2])
    cut = SR.best_cut(wrong, right)
    # cut at 5: wrong >= 5 is 3 of 4; right < 5 is 6 of 6 -> J = 0.75 + 1 - 1
    assert cut["cut"] == 5
    assert (cut["wrong_ge"], cut["right_ge"]) == (3, 0)
    assert cut["j"] == 0.75


def test_intervals_overlap() -> None:
    assert SR.separates(wrong_ge=40, n_wrong=50, right_ge=2, n_right=50)
    assert not SR.separates(wrong_ge=10, n_wrong=20, right_ge=8, n_right=20)
