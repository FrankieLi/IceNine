"""scripts/sweep_audit + scripts/mc_mechanism pure helpers."""

import math
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).parent.parent
for sub in ("scripts/sweep_audit", "scripts/mc_mechanism", "scripts/common"):
    sys.path.insert(0, str(ROOT / sub))

import audit_helpers as A  # noqa: E402
import mc_helpers as H  # noqa: E402

from icenine.orientation_search import MCOptimizer, QuaternionGrid  # noqa: E402


def test_parse_log_line() -> None:
    ln = (
        "  riemannian_adam_geoopt       hp=  2  vox=  414  pert=5°  misori=  4.766°  q=0.2611  "
        "evals=  201  0.8s"
    )
    r = A.parse_log_line(ln)
    assert r is not None
    assert (r["opt"], r["hp"], r["vox"], r["pert"], r["ev"]) == (
        "riemannian_adam_geoopt",
        2,
        414,
        5,
        201,
    )
    assert r["mis"] == 4.766 and r["t"] == 0.8
    assert A.parse_log_line("    idx=   34  φ1= 135.07°") is None
    assert A.parse_log_line("Total: 58500 runs in 806.6 min") is None


def test_success_strict() -> None:
    assert A.success([0.4999, 0.5, 0.6], 0.5).tolist() == [True, False, False]


def test_split_complete_blocks() -> None:
    keys = [(0, 1.0), (0, 2.0), (1, 1.0), (0, 1.0), (0, 2.0), (0, 0.5), (0, 1.0), (1, 0.5)]
    # (0, 1.0) after (0, 2.0) restarts; a repeated key in a later block does not split wrongly
    assert A.split_complete_blocks(keys) == [(0, 3), (3, 5), (5, 8)]
    assert A.split_complete_blocks([]) == []


def test_pick_best_is_order_independent() -> None:
    rows = [
        dict(k=5, n_evals=201, hp=3),
        dict(k=5, n_evals=101, hp=4),
        dict(k=5, n_evals=101, hp=2),
        dict(k=4, n_evals=1, hp=0),
    ]
    assert A.pick_best(rows)["hp"] == 2
    assert A.pick_best(rows[::-1])["hp"] == 2


def test_quantiles_ignore_nan() -> None:
    q = A.quantiles([1.0, 2.0, 3.0, float("nan")], (0.5,))
    assert q == [2.0]
    assert math.isnan(A.quantiles([], (0.5,))[0])


def test_min_ergodic() -> None:
    assert H.min_ergodic(2.5, 1.0, 200) == 31
    assert H.min_ergodic(2.5, 0.5, 200) == 250  # 8x per halving
    assert H.min_ergodic(1.0, 0.0, 77) == 77


def test_expected_progress_and_argmax() -> None:
    st = H.expected_progress(1.0, 0.5, np.array([0.5, 2.0, 1.0]), np.array([0.2, 0.9, 0.5]))
    assert st["p_improve"] == 1 / 3  # strict: the tie at 1.0 is not an improvement
    assert st["cost_prog"] == 0.5 / 3
    assert abs(st["dist_prog"] - 0.3 / 3) < 1e-12
    grid = np.array([0.1, 0.2, 0.3])
    assert H.argmax_step(grid, np.array([0.5, 0.5, 0.1])) == 0.1  # tie -> smallest step
    assert H.argmax_step(grid, np.array([0.0, 0.0, 0.0])) is None
    assert H.argmax_step(grid, np.array([np.nan] * 3)) is None


class _Recorder:
    """Cost function that records the orientations MCOptimizer evaluates and never improves."""

    def __init__(self) -> None:
        self.mats: list = []

    def evaluate(self, R, vertices, phase):
        self.mats.append(np.array(R, dtype=float))
        return SimpleNamespace(cost=1.0 + len(self.mats))


def test_mc_proposals_match_mcoptimizer_draws() -> None:
    """With no improvement and no restart, optimize() evaluates the start, then (the C++ block
    evaluates its start again) n proposals around the start; mc_proposals must reproduce them for
    the same RNG stream."""
    R0 = np.eye(3)
    step = math.radians(0.13)
    rec = _Recorder()
    mc = MCOptimizer(rec, None, 0, rng=np.random.default_rng(5))
    mc.optimize(R0, math.radians(0.33), step, 6, 0, 0.0)
    mats, ang = H.mc_proposals(R0, step, 6, np.random.default_rng(5), QuaternionGrid())
    assert len(rec.mats) == 8
    np.testing.assert_allclose(np.array(rec.mats[2:]), mats, atol=0, rtol=0)
    assert np.all(ang > 0) and np.all(ang < 0.5)


def test_corner_angle_ratio_bounds_sampled_proposals() -> None:
    ratio = H.corner_angle_ratio()
    assert 1.7 < ratio < 1.73
    step = math.radians(0.13)
    _, ang = H.mc_proposals(np.eye(3), step, 500, np.random.default_rng(1), QuaternionGrid())
    assert ang.max() <= ratio * math.degrees(step) * (1 + 1e-9)
    assert 0.9 * math.degrees(step) < np.median(ang) < 1.0 * math.degrees(step)


def test_restart_jump_matches_mcoptimizer_restart() -> None:
    """A run with a cost that never improves restarts after its first block; the restart
    orientation (the next block's start) must equal q(delta) * q(initial) with delta drawn by
    restart_jump_angles' convention (U(-r, r), r = tan(box)/sqrt(48), about the initial one)."""
    box = math.radians(0.33)
    rec = _Recorder()
    mc = MCOptimizer(rec, None, 0, rng=np.random.default_rng(9))
    n_ergodic = H.min_ergodic(box, math.radians(0.13), 10**6)
    mc.optimize(np.eye(3), box, math.radians(0.13), n_ergodic + 5, 1, 0.0)
    # evaluations: initial, block start, n_ergodic proposals, then the next block's start = the
    # restart point
    restart_mat = rec.mats[2 + n_ergodic]
    ang_opt = math.degrees(np.arccos(np.clip((np.trace(restart_mat) - 1) / 2, -1, 1)))
    # replay the rng: n_ergodic proposals (3 draws each), then the 3 restart draws
    rng = np.random.default_rng(9)
    rng.uniform(size=3 * n_ergodic)
    ang = H.restart_jump_angles(box, 1, rng, QuaternionGrid())[0]
    assert abs(ang_opt - ang) < 1e-6
    # the jump is smaller than the box (median 0.49 box, corner 0.85 box), unlike the old +-box/2
    big = H.restart_jump_angles(box, 2000, np.random.default_rng(0), QuaternionGrid())
    assert 0.4 * math.degrees(box) < np.median(big) < 0.6 * math.degrees(box)
    assert big.max() < 0.9 * math.degrees(box)
