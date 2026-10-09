"""MCOptimizer.optimize against C++ COrientationMC::RandomRestartZeroTemp
(Src/OrientationSearch.cpp:297-363, inner block ZeroTemperatureOptimization:100-135).

Scripted blocks pin the block structure (block length fixed at nMinErgodicSteps =
int(2 (box/step)^3) from the initial step, last block cut to the remaining budget), the comparison
of a block result with the global best (failure when not strictly lower), the restart
(x, y, z ~ U(-r, r), r = tan(box)/sqrt(48), applied to the INITIAL orientation, step reset), the
update order (success: global and current := block result, step halved, failure counter zeroed),
and the stop rule
(MaxConvergenceCost first, then successive failures > SuccessiveRestarts). The step budget counts
trial steps; a block also evaluates its start. Cases follow the C++ debug trace kept on the
branch feature/finisher-mc-study (benchmarks/phase_d_seed_diag/mc_parity/; not on develop).
"""

import math
from types import SimpleNamespace
from typing import Any, List, Tuple

import numpy as np

from icenine.orientation_search import MCOptimizer
from icenine.sampling import _quat_multiply, matrix_to_quaternion, quaternion_to_matrix

BOX = math.radians(5.0 / 3.0)
STEP = math.radians(2.0 / 3.0)  # 2 (box/step)^3 = 31.25 -> 31 steps per block
N_MIN = 31


class _Cost:
    def __init__(self) -> None:
        self.n = 0

    def evaluate(self, orientation: Any, vertices: Any, phase: int) -> Any:
        self.n += 1
        return SimpleNamespace(cost=1.0)


class _Scripted:
    """Replaces _mc_block: returns the scripted block costs, records (start_q, step, n)."""

    def __init__(self, mc: MCOptimizer, costs: List[float], marker: float = 1e-3) -> None:
        self.costs = list(costs)
        self.calls: List[Tuple[np.ndarray, float, int]] = []
        self.results: List[np.ndarray] = []
        self.marker = marker
        grid = mc._grid_gen
        orig = grid.get_near_identity_point
        self.offsets: List[Tuple[float, float, float]] = []

        def rec_grid(x: float, y: float, z: float) -> np.ndarray:
            self.offsets.append((x, y, z))
            return orig(x, y, z)

        grid.get_near_identity_point = rec_grid  # type: ignore[method-assign]
        self.orig_grid = orig

        def fake(start_q: np.ndarray, step: float, n: int) -> Any:
            k = len(self.calls)
            self.calls.append((np.array(start_q), step, n))
            cost = self.costs[k] if k < len(self.costs) else 2.0  # past the script: failures
            # a distinguishable block result: start rotated by a k-dependent small rotation
            res = _quat_multiply(orig(self.marker * (k + 1), 0.0, 0.0), start_q)
            self.results.append(res)
            return res, cost, SimpleNamespace(cost=cost, hit_ratio=1.0)

        mc._mc_block = fake  # type: ignore[method-assign]


def _mc(seed: int = 0) -> MCOptimizer:
    return MCOptimizer(_Cost(), None, 0, np.random.default_rng(seed))  # type: ignore[arg-type]


def _run(mc: MCOptimizer, steps: int, restarts: int, conv: float = 0.0, **kw: Any) -> Any:
    return mc.optimize(np.eye(3), BOX, STEP, steps, restarts, conv, **kw)


def test_block_length_is_fixed_and_last_block_is_truncated() -> None:
    mc = _mc()
    sc = _Scripted(mc, [0.9, 0.8, 0.7, 0.6])  # four successes, budget 100
    _run(mc, 100, 2)
    ns = [c[2] for c in sc.calls]
    assert ns == [N_MIN, N_MIN, N_MIN, 100 - 3 * N_MIN]  # 31 31 31 7: n not recomputed on halving
    steps = [c[1] for c in sc.calls]
    assert steps == [STEP, STEP / 2, STEP / 4, STEP / 8]  # halved after each success only
    assert mc.last_run["steps_run"] == 100 and mc.last_run["stop"] == 0
    assert mc.last_run["n_blocks"] == 4


def test_failure_is_not_strictly_lower_and_resets_step_and_block_length() -> None:
    mc = _mc()
    # initial cost is 1.0: a block ending at exactly 1.0 fails (>=); then a success; then a
    # block equal to the new global best fails
    sc = _Scripted(mc, [1.0, 0.5, 0.5, 0.4], marker=1e-3)
    _run(mc, 200, 5)
    steps = [c[1] for c in sc.calls]
    assert steps[0] == STEP and steps[1] == STEP  # failure: step stays the initial one
    assert steps[2] == STEP / 2  # success halved it
    assert steps[3] == STEP  # failure at equal cost reset it
    assert [c[2] for c in sc.calls[:4]] == [N_MIN] * 4


def test_restart_base_range_and_draw_order() -> None:
    """A failed block restarts at delta * INITIAL, delta from x, y, z ~ U(-r, r) drawn in that
    order with r = tan(box)/sqrt(48); the block after a success continues from the block result,
    and a later failure again restarts from the initial orientation (not the improved best)."""
    mc = _mc(7)
    sc = _Scripted(mc, [0.8, 2.0, 0.6, 2.0, 0.5])  # success, fail, success, fail, success
    _run(mc, 5 * N_MIN, 5)  # the fifth block starts at the second restart point
    q0 = matrix_to_quaternion(np.eye(3))
    r = math.tan(BOX) / math.sqrt(48.0)
    twin = np.random.default_rng(7)  # the fake blocks consume no random numbers
    off1 = tuple(twin.uniform(-r, r) for _ in range(3))
    off2 = tuple(twin.uniform(-r, r) for _ in range(3))
    assert sc.offsets == [off1, off2]
    # block 0 starts at the initial orientation
    np.testing.assert_allclose(sc.calls[0][0], q0, atol=1e-15)
    # block 1 continues from block 0's result (success)
    np.testing.assert_array_equal(sc.calls[1][0], sc.results[0])
    # block 2 (after the failure of block 1) starts at delta1 * initial
    np.testing.assert_allclose(sc.calls[2][0], _quat_multiply(sc.orig_grid(*off1), q0), atol=1e-15)
    # block 3 continues from block 2's result; block 3 fails -> block 4 starts at delta2 * initial
    np.testing.assert_array_equal(sc.calls[3][0], sc.results[2])
    # the second restart base is the initial orientation, not the improved best (block 2's result)
    np.testing.assert_allclose(sc.calls[4][0], _quat_multiply(sc.orig_grid(*off2), q0), atol=1e-15)
    assert mc.last_run["n_restarts"] == 2


def test_successive_failure_counter_and_stop_rule() -> None:
    # SuccessiveRestarts = 2: the third consecutive failure stops (count 3 > 2)
    mc = _mc()
    sc = _Scripted(mc, [2.0, 2.0, 2.0, 2.0])
    _run(mc, 1000, 2)
    assert len(sc.calls) == 3 and mc.last_run["stop"] == 1
    assert mc.last_run["since_improve"] == 3 and mc.last_run["steps_run"] == 3 * N_MIN
    # SuccessiveRestarts = 0: the first failure stops
    mc = _mc()
    sc = _Scripted(mc, [2.0])
    _run(mc, 1000, 0)
    assert len(sc.calls) == 1 and mc.last_run["stop"] == 1
    # a success resets the counter: F F S F F F stops after the sixth block, not the third
    mc = _mc()
    sc = _Scripted(mc, [2.0, 2.0, 0.9, 2.0, 2.0, 2.0])
    _run(mc, 1000, 2)
    assert len(sc.calls) == 6 and mc.last_run["since_improve"] == 3
    assert mc.last_run["n_accept"] == 1 and mc.last_run["n_restarts"] == 5


def test_convergence_checked_after_every_block_and_before_the_restart_rule() -> None:
    mc = _mc()
    sc = _Scripted(mc, [0.5, 0.00005])
    _run(mc, 1000, 2, conv=1e-4)
    assert len(sc.calls) == 2 and mc.last_run["stop"] == 2

    # C++ quirk: an initial cost already below the threshold is only noticed after the first block
    # (the check follows the block, success or not)
    class _Low(_Cost):
        def evaluate(self, orientation: Any, vertices: Any, phase: int) -> Any:
            return SimpleNamespace(cost=0.00001)

    mc = MCOptimizer(_Low(), None, 0, np.random.default_rng(0))  # type: ignore[arg-type]
    sc = _Scripted(mc, [2.0])
    _run(mc, 1000, 2, conv=1e-4)
    assert len(sc.calls) == 1 and mc.last_run["stop"] == 2

    # convergence is checked before the restart rule: a failed block with max_restarts = 0 would
    # stop with 1 (count 1 > 0), but the global cost is below the threshold, so it stops with 2
    mc = MCOptimizer(_Low(), None, 0, np.random.default_rng(0))  # type: ignore[arg-type]
    sc = _Scripted(mc, [2.0])
    _run(mc, 1000, 0, conv=1e-4)
    assert len(sc.calls) == 1 and mc.last_run["since_improve"] == 1 and mc.last_run["stop"] == 2


def test_trajectory_has_one_record_per_block() -> None:
    mc = _mc()
    _Scripted(mc, [0.5, 2.0, 0.4])
    traj: List[dict] = []
    _run(mc, 3 * N_MIN, 5, trajectory=traj)
    assert [t["event_type"] for t in traj] == ["mc_accept", "mc_restart", "mc_accept"]
    assert [t["step"] for t in traj] == [N_MIN, 2 * N_MIN, 3 * N_MIN]
    assert math.isclose(traj[0]["cur_step_rad"], STEP / 2)
    assert math.isclose(traj[1]["cur_step_rad"], STEP)  # reset on the failure
    assert math.isclose(traj[2]["cur_step_rad"], STEP / 2)


def test_evaluation_count_is_steps_plus_block_starts() -> None:
    """Budget counts trial steps; evaluations = 1 (initial) + sum over blocks of (1 + steps)."""
    cost = _Cost()
    mc = MCOptimizer(cost, None, 0, np.random.default_rng(0))  # type: ignore[arg-type]
    _run(mc, 100, 2)  # constant cost 1.0: every block fails; stops after 3 blocks
    assert mc.last_run["n_blocks"] == 3 and mc.last_run["steps_run"] == 93
    assert cost.n == 1 + 3 * (1 + N_MIN)
    # the quick MC of the reconstruction (10 steps, 5 restarts): one block, 12 evaluations
    cost = _Cost()
    mc = MCOptimizer(cost, None, 0, np.random.default_rng(0))  # type: ignore[arg-type]
    mc.optimize(np.eye(3), math.radians(5 / 3), math.radians(2 / 3), 10, 5, 1e-4)
    assert cost.n == 12 and mc.last_run["n_blocks"] == 1


def test_block_proposals_follow_the_c_inner_loop() -> None:
    """_mc_block: x, y, z ~ U(-r, r), r = tan(step)/sqrt(12); trial = delta * (block best);
    strict improvement only; the start is evaluated first."""
    seen: List[np.ndarray] = []

    class _Rec:
        def evaluate(self, R: Any, vertices: Any, phase: int) -> Any:
            seen.append(np.array(R))
            return SimpleNamespace(cost=1.0 + 0.0 * len(seen))  # never improves

    mc = MCOptimizer(_Rec(), None, 0, np.random.default_rng(3))  # type: ignore[arg-type]
    q0 = matrix_to_quaternion(np.eye(3))
    q, cost, _ = mc._mc_block(q0, STEP, 5)
    assert len(seen) == 6 and cost == 1.0
    np.testing.assert_array_equal(q, q0)
    r = math.tan(STEP) / math.sqrt(12.0)
    twin = np.random.default_rng(3)
    grid = mc._grid_gen
    for k in range(5):
        x, y, z = (twin.uniform(-r, r) for _ in range(3))
        expect = quaternion_to_matrix(_quat_multiply(grid.get_near_identity_point(x, y, z), q0))
        np.testing.assert_allclose(seen[1 + k], expect, atol=1e-15)
