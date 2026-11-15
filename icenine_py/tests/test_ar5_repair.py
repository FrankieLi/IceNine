"""AR5 step 19 (task ar5s): the repair of contested voxels after the boundary pass, on hand-made
maps and a fake fit. Plan "AR5 step 19", phase A step 2: ratio and margin; a voxel with no other
pass end is never flagged; candidates as the pass; rule ``hit`` with the current end first on ties;
one MC fit per start; every evaluation counted; its own generator leaves other draws alone; guards
on ``held2`` and ``held4``.
"""

import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "autoresearch" / "ar5"))

import repair_guard as RG  # noqa: E402
import repair_rules as RR  # noqa: E402

from icenine.boundary_pass import boundary_candidates  # noqa: E402
from icenine.ray_indexing import proper_rotations  # noqa: E402
from icenine.symmetry import create_cubic_symmetry  # noqa: E402

SYM = proper_rotations(create_cubic_symmetry(3.61))


def rot_z(deg: float) -> np.ndarray:
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def line_nbrs(n: int) -> List[List[int]]:
    """A row of ``n`` voxels, each next to the one or two beside it (radius 1 on a line)."""
    return [[j for j in (i - 1, i + 1) if 0 <= j < n] for i in range(n)]


# ---- constants --------------------------------------------------------------------------------


def test_the_frozen_constants_are_the_plans() -> None:
    assert RR.RATIO_MAX == 0.95 and RR.MARGIN_MAX == 0.15
    assert RR.FAR_DEG == 1.0 and RR.CAND_DEG == 0.5 and RR.WRONG_DEG == 1.0
    assert RR.REPAIR_SEED == 0 and RR.REPAIR_OPTIMIZER == "mc"
    assert RR.J5_TOTAL == dict(
        held3=23_503_385, held=22_939_812, mid=10_589_493, grad=9_041_686, big=3_567_098
    )
    assert RR.FLAGGED == dict(held3=193, held=183, mid=143, grad=238, big=56)
    assert RR.FITS == dict(held3=1_217, held=1_035, mid=643, grad=967, big=161)
    assert RR.J5_WRONG == dict(held3=54, held=67, mid=35, grad=56, big=21)
    assert RR.max_share("grad") == 3.5 and RR.max_share("mid") == 2.5


# ---- ratio and margin on hand-made maps -------------------------------------------------------


def test_ratio_is_the_hit_over_the_median_of_the_neighbours() -> None:
    hit = np.array([0.9, 1.0, 0.5, 1.0, 1.0])
    nb = line_nbrs(5)
    r = RR.ratio_of(hit, nb)
    # voxel 0: neighbour {1} -> 0.9 / 1.0; voxel 2: neighbours {1, 3} -> median 1.0 -> 0.5
    assert r[0] == pytest.approx(0.9)
    assert r[2] == pytest.approx(0.5)
    assert r[1] == pytest.approx(1.0 / np.median([0.9, 0.5]))
    # a voxel with no neighbour is never under the ratio line
    assert RR.ratio_of(np.array([0.1]), [[]])[0] == np.inf


def test_margin_is_the_hit_over_the_best_other_end_and_a_voxel_with_none_is_never_flagged() -> None:
    hit = np.array([0.90, 0.90, 0.90])
    other = np.array([0.88, 0.70, -1.0])  # voxel 2: no other pass end
    m = RR.margin_of(hit, other)
    assert m[0] == pytest.approx(0.02) and m[1] == pytest.approx(0.20)
    assert m[2] == RR.NO_OTHER >= RR.MARGIN_MAX
    # all three have a low ratio; only voxel 0 is contested (margin under 0.15); voxel 2 never
    nb = [[1], [0, 2], [1]]
    hit_n = np.array([0.90, 1.0, 0.90])
    flag, ratio, margin = RR.contested(hit_n, other, nb)
    assert bool(ratio[0] < RR.RATIO_MAX) and bool(ratio[2] < RR.RATIO_MAX)
    assert list(flag) == [True, False, False]


def test_the_two_tests_are_strict() -> None:
    nb = [[1], [0]]
    # ratio exactly 0.95 (0.95 over a neighbour at 1.0) is not flagged, even with a margin of 0
    flag, ratio, _ = RR.contested(np.array([0.95, 1.0]), np.array([0.95, -1.0]), nb)
    assert ratio[0] == 0.95 and not flag[0]
    # just under 0.95 with a margin under 0.15 is flagged
    flag, _, _ = RR.contested(np.array([0.9499, 1.0]), np.array([0.9, -1.0]), nb)
    assert flag[0]
    # a margin just under 0.15 is flagged, just over is not (hit 0.5: ratio 0.5)
    hit = np.array([0.5, 1.0])
    flag, _, mg = RR.contested(hit, np.array([0.5 - 0.1499, -1.0]), nb)
    assert flag[0] and mg[0] < 0.15
    flag, _, mg = RR.contested(hit, np.array([0.5 - 0.1501, -1.0]), nb)
    assert not flag[0] and mg[0] > 0.15


def test_other_end_hit_takes_far_try_ends_and_a_far_old_end_only() -> None:
    R_before = np.stack([rot_z(0.0), rot_z(0.0), rot_z(30.0)])
    R_after = np.stack([rot_z(0.0), rot_z(30.0), rot_z(30.0)])  # voxel 1 took a try
    old_hit = np.array([0.80, 0.70, 0.90])
    # voxel 0: a try that ended 0.4 degree from the kept end (not "other") and one far (0.6 hit)
    # voxel 1: its kept end is the try; the old end is 30 degrees away: other = old hit 0.70
    t_voxel = np.array([0, 0, 1])
    t_end = np.stack([rot_z(0.4), rot_z(20.0), rot_z(30.0)])
    t_hit = np.array([0.85, 0.60, 0.95])
    other = RR.other_end_hit(R_after, R_before, old_hit, t_voxel, t_end, t_hit, SYM)
    assert other[0] == pytest.approx(0.60)  # the near try (0.85) is the kept end's twin: excluded
    assert other[1] == pytest.approx(0.70)  # the try at voxel 1 IS the kept end: excluded
    assert other[2] == -1.0  # untried, old end kept: no other end


def test_after_pass_state_follows_rule_hit_with_the_old_end_first_on_ties() -> None:
    R_before = np.stack([rot_z(0.0), rot_z(0.0), rot_z(0.0)])
    old_score = np.array([[2.0, 0.90, 0.5], [2.0, 0.90, 0.5], [2.0, 0.90, 0.5]])
    t_voxel = np.array([0, 1, 1, 2])
    t_end = np.stack([rot_z(10.0), rot_z(11.0), rot_z(12.0), rot_z(13.0)])
    t_cost = np.array([1.0, 1.5, 1.4, 2.0])
    t_hit = np.array([0.90, 0.95, 0.95, 0.90])  # voxel 0: tie in hit, lower cost -> the try
    R, cost, hit = RR.after_pass_state(R_before, old_score, t_voxel, t_end, t_cost, t_hit)
    assert np.allclose(R[0], rot_z(10.0))  # same hit, lower cost
    assert np.allclose(R[1], rot_z(12.0)) and cost[1] == 1.4  # same hit, lower cost between tries
    assert np.allclose(R[2], rot_z(0.0)) and cost[2] == 2.0  # all equal: the old end first


# ---- candidates as the pass -------------------------------------------------------------------


def test_candidates_are_the_passs_and_a_non_boundary_voxel_has_none() -> None:
    R = np.stack([rot_z(0.0), rot_z(0.1), rot_z(20.0), rot_z(20.2), rot_z(40.0)])
    nb = [[1, 2, 3], [0, 2, 3], [0, 1, 3, 4], [0, 1, 2], [2]]
    _isb, ref = boundary_candidates(R, nb, SYM)
    got = RR.repair_candidates(R, nb, SYM, [0, 2, 4])
    assert got[0] == ref[0] == [2]  # 20 and 20.2 merged: the first one kept
    assert got[2] == ref[2]
    assert got[4] == ref[4]
    R2 = np.stack([rot_z(0.0), rot_z(0.1)])
    assert RR.repair_candidates(R2, [[1], [0]], SYM, [0, 1]) == {0: [], 1: []}


# ---- the repair ---------------------------------------------------------------------------------


class FakeFit:
    """A fit that records its calls and returns scripted ends."""

    def __init__(self, script: dict) -> None:
        self.calls: List[Tuple[int, np.ndarray]] = []
        self.draws: List[float] = []
        self.script = script  # (voxel, rounded start angle) -> (end angle, cost, hit, evals)

    def __call__(self, v: int, R0: np.ndarray, rng: np.random.Generator):
        self.calls.append((v, R0.copy()))
        self.draws.append(float(rng.random()))
        ang = round(float(np.degrees(np.arctan2(R0[1, 0], R0[0, 0]))), 3)
        end, c, h, ev = self.script[(v, ang)]
        return rot_z(end), c, h, ev


def small_case():
    R = np.stack([rot_z(0.0), rot_z(20.0), rot_z(40.0)])
    cost = np.array([2.0, 2.0, 2.0])
    hit = np.array([0.80, 0.90, 0.90])
    nb = [[1, 2], [0, 2], [0, 1]]
    cands = RR.repair_candidates(R, nb, SYM, [0])
    return R, cost, hit, cands


def test_one_fit_per_start_every_evaluation_counted_and_rule_hit_keeps_the_best() -> None:
    R, cost, hit, cands = small_case()
    assert cands[0] == [1, 2]
    fake = FakeFit(
        {
            (0, 0.0): (0.5, 1.9, 0.82, 101),  # own end: a little better, still low
            (0, 20.0): (20.0, 1.0, 0.97, 130),  # from candidate 1: much better
            (0, 40.0): (40.0, 1.2, 0.97, 117),  # from candidate 2: same hit, higher cost
        }
    )
    out = RR.repair(R, cost, hit, [0], cands, fake)
    assert [c[0] for c in fake.calls] == [0, 0, 0]  # one fit per start, own end first
    assert np.allclose(fake.calls[0][1], rot_z(0.0))
    assert np.allclose(fake.calls[1][1], R[1]) and np.allclose(fake.calls[2][1], R[2])
    assert out.n_evals_fits == 101 + 130 + 117 == out.n_evals and out.n_evals_score == 0
    assert out.changed == [0] and np.allclose(out.R[0], rot_z(20.0))
    assert out.hit[0] == 0.97 and out.cost[0] == 1.0
    assert np.allclose(R[0], rot_z(0.0)) and hit[0] == 0.80  # inputs are not changed


def test_a_tie_keeps_the_current_end() -> None:
    R, cost, hit, cands = small_case()
    fake = FakeFit(
        {
            (0, 0.0): (0.3, 2.0, 0.80, 100),  # the same hit and cost as the current end
            (0, 20.0): (20.0, 2.0, 0.80, 100),
            (0, 40.0): (40.0, 2.5, 0.80, 100),
        }
    )
    out = RR.repair(R, cost, hit, [0], cands, fake)
    assert out.changed == [] and np.allclose(out.R[0], rot_z(0.0))
    assert out.n_evals == 300


def test_a_missing_stored_hit_costs_one_counted_evaluation() -> None:
    R, cost, hit, cands = small_case()
    hit = hit.copy()
    hit[0] = np.nan
    fake = FakeFit(
        {
            (0, 0.0): (0.0, 2.0, 0.5, 90),
            (0, 20.0): (20.0, 2.0, 0.5, 90),
            (0, 40.0): (40.0, 2.0, 0.5, 90),
        }
    )
    with pytest.raises(ValueError):
        RR.repair(R, cost, hit, [0], cands, fake)
    out = RR.repair(R, cost, hit, [0], cands, fake, score=lambda v, Rv: (3.0, 0.4))
    assert out.n_evals_score == 1 and out.n_evals == 271
    # the scored current end (cost 3.0, hit 0.4) loses to every fit (hit 0.5); ties among the
    # fits go to the lowest cost, then the first (the own end's fit)
    assert out.hit[0] == 0.5 and out.changed == [0] and np.allclose(out.R[0], rot_z(0.0))


def test_voxels_are_taken_in_index_order_and_candidates_are_the_snapshots() -> None:
    R = np.stack([rot_z(0.0), rot_z(20.0), rot_z(40.0)])
    cost = np.full(3, 2.0)
    hit = np.array([0.8, 0.8, 0.8])
    nb = [[1, 2], [0, 2], [0, 1]]
    cands = RR.repair_candidates(R, nb, SYM, [0, 1])
    seen: List[Tuple[int, float]] = []

    def fit(v, R0, rng):
        seen.append((v, round(float(np.degrees(np.arctan2(R0[1, 0], R0[0, 0]))), 3)))
        # voxel 0 jumps to 20 degrees; voxel 1's candidate 0 must still start from 0 degrees
        return rot_z(20.0 if v == 0 else 0.0), 1.0, 0.99 if v == 0 else 0.5, 10

    out = RR.repair(R, cost, hit, [1, 0], cands, fit)  # given out of order
    assert [v for v, _ in seen] == [0, 0, 0, 1, 1, 1]
    assert (1, 0.0) in seen  # the snapshot's orientation of voxel 0, not the repaired one
    assert out.flagged == [0, 1] and 0 in out.changed


def test_its_own_generator_default_rng_0_and_other_draws_do_not_move() -> None:
    R, cost, hit, cands = small_case()
    np.random.seed(123)
    before = np.random.get_state()[1].copy()
    g = np.random.default_rng(5)
    g_state = g.bit_generator.state
    fake = FakeFit(
        {
            (0, 0.0): (0.0, 2.0, 0.8, 1),
            (0, 20.0): (20.0, 2.0, 0.8, 1),
            (0, 40.0): (40.0, 2.0, 0.8, 1),
        }
    )
    RR.repair(R, cost, hit, [0], cands, fake)
    ref = np.random.default_rng(0)
    assert fake.draws == [float(ref.random()) for _ in range(3)]  # one generator, seed 0
    assert np.array_equal(np.random.get_state()[1], before)  # the global state did not move
    assert g.bit_generator.state == g_state  # another generator did not move


def test_no_flagged_voxel_no_fit_and_no_evaluation() -> None:
    R, cost, hit, cands = small_case()
    fake = FakeFit({})
    out = RR.repair(R, cost, hit, [], {}, fake)
    assert fake.calls == [] and out.n_evals == 0 and out.changed == []


# ---- scoring ------------------------------------------------------------------------------------


def test_fixed_broken_and_the_gate_tests() -> None:
    old = np.array([True, True, False, False, True])
    new = np.array([False, True, True, False, False])
    assert RR.fixed_broken(old, new) == (2, 1)
    assert RR.gate_g12("mid", 5, 1, 1.0)["ok"]  # net 4 fixed, cost fine
    assert not RR.gate_g12("mid", 3, 1, 1.0)["ok_net"]  # net 2 < 3
    assert not RR.gate_g12("mid", 9, 0, 2.6)["ok_cost"]
    assert RR.gate_g12("held", 0, 0, 1.0)["ok"]  # g2: broken minus fixed at most 0
    assert not RR.gate_g12("held", 3, 4, 1.0)["ok_net"]
    assert RR.gate_g12("grad", 5, 0, 3.4)["ok"] and not RR.gate_g12("big", 5, 0, 3.4)["ok_cost"]
    w = RR.wrong_mask(np.stack([rot_z(0.0), rot_z(1.5)]), np.stack([rot_z(0.5)] * 2), SYM)
    assert list(w) == [False, True]


# ---- guards --------------------------------------------------------------------------------------


def test_held2_is_refused_everywhere(tmp_path: Path) -> None:
    for name in ("held2", "xheld2"):
        with pytest.raises(SystemExit) as e:
            RG.refuse_sample(name)
        assert e.value.code == 2
    for p in (tmp_path / "J_held2" / "run.json", tmp_path / "passes_held2", "x/held2_clean.npy"):
        with pytest.raises(SystemExit) as e:
            RG.guarded(p, cache=tmp_path)
        assert e.value.code == 2
    with pytest.raises(SystemExit):
        RG.refuse_sample("nosuch")


def test_held4_is_refused_before_the_gate_and_its_files_before_its_runs_end(
    tmp_path: Path,
) -> None:
    with pytest.raises(SystemExit) as e:
        RG.refuse_sample("held4", gate=tmp_path / "gate.json")
    assert e.value.code == 2
    (tmp_path / "gate.json").write_text("{}")
    RG.refuse_sample("held4", gate=tmp_path / "gate.json")  # now allowed
    RG.refuse_sample("mid", gate=tmp_path / "none.json")  # disclosed samples need no gate
    cache = tmp_path / "cache"
    cache.mkdir()
    p = cache / "J5_held4" / "run.json"
    with pytest.raises(SystemExit) as e:
        RG.guarded(p, cache=cache)
    assert e.value.code == 2
    RG.guarded(cache / "J5_mid" / "run.json", cache=cache)  # other samples are open
    (cache / RG.H4_MARKER).write_text("{}")
    assert RG.guarded(p, cache=cache) == p  # marker: held4's runs ended


def test_the_freeze_hashes_this_steps_files_and_the_library() -> None:
    import repair_freeze as RF

    cur = RF.current()
    for name in RF.OWN:
        assert f"scripts/autoresearch/ar5/{name}" in cur["files"]
    assert "icenine/boundary_pass.py" in cur["library"]
    assert "uv.lock" in cur["library"]
    # the scripts that only apply the frozen rule are not in the list
    for name in ("repair_run.py", "repair_gate.py", "repair_flags.py", "repair_run_all.sh"):
        assert f"scripts/autoresearch/ar5/{name}" not in cur["files"]
    consts = RF.CONSTANTS
    assert consts["ratio_max"] == 0.95 and consts["margin_max"] == 0.15
    assert consts["repair_optimizer"] == "mc" and consts["repair_seed"] == 0


def test_the_link_script_refuses_held2_and_held4_and_unknown_names(tmp_path: Path) -> None:
    import subprocess

    script = ROOT / "scripts" / "autoresearch" / "ar5" / "repair_link_cache.sh"
    for bad in ("held2", "held4", "nosuch"):
        r = subprocess.run(
            ["bash", str(script), "--dry-run", "--cache", str(tmp_path), bad],
            capture_output=True,
            text=True,
        )
        assert r.returncode == 2 and "refused argument" in r.stderr
    assert not list(tmp_path.iterdir())


def test_the_launcher_only_accepts_four_scripts() -> None:
    src = (ROOT / "scripts" / "autoresearch" / "ar5" / "repair_launch.py").read_text()
    for name in ("fit_run.py", "fit_pass.py", "high500_run.py", "high500_pass.py"):
        assert f'"{name}",' in src
    allowed = src.split("ALLOWED = (")[1].split(")")[0]
    assert allowed.count(".py") == 4


def test_the_guard_module_names_step_19s_gate_and_marker() -> None:
    assert RG.gate_file().parts[-3:] == ("ar5", "repair", "gate.json")
    assert RG.H4_MARKER == "held4_runs_done.json"
    assert RG.HELD_SAMPLES == ("held", "held3", "held4")
