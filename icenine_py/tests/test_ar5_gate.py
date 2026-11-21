"""AR5 step 7 (``scripts/autoresearch/ar5/gate_rules.py``): the replay of a boundary pass with the
skip, the fixed / broken count and the decision rule. Pure numpy; written before the module."""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import gate_rules as GR  # noqa: E402


def _pass():
    """5 voxels; voxel 0 has own hit 0.99 (2 tries), voxel 1 has 0.90 (3 tries), voxel 3 has 0.95
    (1 try), voxel 4 has 0.5 (1 try). Voxel 2 has no try."""
    old_hit = np.array([0.99, 0.90, 1.0, 0.95, 0.5])
    old_cost = np.array([0.1, 0.4, 0.0, 0.2, 0.9])
    t_voxel = np.array([0, 0, 1, 1, 1, 3, 4])
    t_hit = np.array([0.995, 0.97, 0.93, 0.99, 0.80, 0.97, 0.6])
    t_cost = np.array([0.09, 0.2, 0.3, 0.1, 0.5, 0.1, 0.8])
    t_evals = np.array([61, 61, 61, 61, 60, 61, 61])
    return old_hit, old_cost, t_voxel, t_hit, t_cost, t_evals


def test_skip_mask_is_at_least_h() -> None:
    m = GR.skip_mask(np.array([0.9499, 0.95, 0.951, 1.0]), 0.95)
    assert m.tolist() == [False, True, True, True]


def test_h_zero_or_above_one_skips_nothing_or_everything() -> None:
    x = np.array([0.0, 0.5, 1.0])
    assert GR.skip_mask(x, 1.01).tolist() == [False, False, False]
    assert GR.skip_mask(x, 0.0).tolist() == [True, True, True]


def test_replay_counts() -> None:
    old_hit, old_cost, tv, th, tc, te = _pass()
    r = GR.replay(old_hit, old_cost, tv, th, tc, te, 0.95)
    assert r["kept"].tolist() == [False, False, True, True, True, False, True]
    assert r["evals_full"] == 61 * 6 + 60
    assert r["evals"] == 61 + 61 + 60 + 61
    assert r["tries_skipped"] == 3 and r["tries_kept"] == 4
    assert r["voxels_skipped"] == 2  # voxels 0 and 3 (voxel 2 has no try, so it is not counted)
    assert r["boundary_voxels"] == 4
    assert abs(r["cut"] - (1 - 243 / 426)) < 1e-12


def test_replay_winner_rule_hit_then_cost_then_earliest() -> None:
    old_hit, old_cost, tv, th, tc, te = _pass()
    r = GR.replay(old_hit, old_cost, tv, th, tc, te, 2.0)  # nothing skipped
    # voxel 0: try 0 (0.995); voxel 1: try 3 (0.99); voxel 3: try 5 (0.97); voxel 4: try 6 (0.6)
    assert r["winner"][0] == 0 and r["winner"][1] == 3 and r["winner"][3] == 5
    assert r["winner"][4] == 6 and r["winner"][2] == -1
    r2 = GR.replay(old_hit, old_cost, tv, th, tc, te, 0.95)
    assert r2["winner"][0] == -1 and r2["winner"][3] == -1  # skipped voxels keep the old end
    assert r2["winner"][1] == 3 and r2["winner"][4] == 6


def test_winner_tie_prefers_old_then_lower_cost_then_earlier() -> None:
    old_hit = np.array([0.5, 0.5])
    old_cost = np.array([0.5, 0.5])
    tv = np.array([0, 0, 1, 1])
    th = np.array([0.5, 0.5, 0.5, 0.5])
    tc = np.array([0.5, 0.4, 0.4, 0.4])
    r = GR.replay(old_hit, old_cost, tv, th, tc, np.full(4, 61), 2.0)
    assert r["winner"][0] == 1  # equal hit; lower cost than the old end wins
    assert r["winner"][1] == 2  # equal everything between tries: earliest; beats old on cost


def test_fixed_broken_pairs() -> None:
    a = np.array([True, True, False, False, True])  # wrong under the reference
    b = np.array([False, True, True, False, True])  # wrong under the new
    fixed, broken = GR.fixed_broken(a, b)
    assert (fixed, broken) == (1, 1)


def test_net_gain_and_allowed_loss_rounds_up_min_one() -> None:
    assert GR.net_gain(75, 35) == 40
    assert GR.allowed_loss(40) == 2  # 5% of 40 = 2.0 exactly (integer arithmetic)
    assert GR.allowed_loss(41) == 3  # 2.05 rounds up
    assert GR.allowed_loss(3) == 1  # at least 1 voxel
    assert GR.allowed_loss(0) == 1
    assert GR.allowed_loss(-5) == 1


def _args(**kw):
    d = dict(
        cut=0.689,
        wrong_before=75,
        wrong_full=35,
        fixed=0,
        broken=0,
        med_full=0.0291,
        med_skip=0.0291,
        total_skip=12_000_000,
        total_s8=26_000_000,
        fixed_s8=13,
        broken_s8=17,
        med_s8=0.0297,
    )
    d.update(kw)
    return d


def test_gate_rule_passes_and_each_part_can_fail() -> None:
    assert GR.gate_rule(**_args())["passes"] is True
    assert GR.gate_rule(**_args(cut=0.49))["passes"] is False
    assert GR.gate_rule(**_args(broken=3))["pass_loss"] is False  # 3 > 2
    assert GR.gate_rule(**_args(broken=3, fixed=1))["pass_loss"] is True  # 2 <= 2
    assert GR.gate_rule(**_args(med_skip=0.0302))["pass_median"] is False  # 0.0011 higher
    assert GR.gate_rule(**_args(total_skip=21_000_000))["pass_s8_evals"] is False  # 19.2% fewer
    assert GR.gate_rule(**_args(broken_s8=19))["pass_s8_wrong"] is False  # 19 - 13 = 6 > 5
    assert GR.gate_rule(**_args(med_s8=0.0240))["pass_s8_median"] is False  # 0.0051 higher


def test_gate_rule_boundaries_inclusive() -> None:
    assert GR.gate_rule(**_args(cut=0.5))["pass_cut"] is True
    assert GR.gate_rule(**_args(med_skip=0.0301, med_full=0.0291))["pass_median"] is True
    assert GR.gate_rule(**_args(total_skip=20_800_000))["pass_s8_evals"] is True  # exactly 20%
    assert GR.gate_rule(**_args(broken_s8=18))["pass_s8_wrong"] is True  # 18 - 13 = 5
