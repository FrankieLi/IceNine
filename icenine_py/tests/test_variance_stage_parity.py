"""VarianceMinimizing stage (MCOptimizer.variance_minimizing_optimize) restart and termination
semantics against C++ AdaptiveSamplingZeroTemp (OrientationSearch.cpp:206-287).

Regression for the Phase D finding: Python restarted from the current global best with offsets
uniform in +-box/2, C++ restarts from the INITIAL orientation with offsets uniform in
+-SubregionRadius (the radius after the doubling), so on a dense image Python restarts landed in the
steep part of the cost basin (cost variance above 0.02^2, budget extended without end) where C++
restarts land on the cost plateau.
"""

import math
from types import SimpleNamespace
from typing import Any, List, Tuple

import numpy as np

from icenine.orientation_search import MCOptimizer
from icenine.sampling import _quat_multiply, matrix_to_quaternion, quaternion_to_matrix


class _Cost:
    def evaluate(self, orientation: Any, vertices: Any, phase: int) -> Any:
        return SimpleNamespace(cost=1.0)


class _Recorder:
    """Scripted subregion runs; records start matrices, radii, steps and restart offsets."""

    def __init__(self, mc: MCOptimizer, script: List[Tuple[float, float]]):
        self.script = list(script)  # (cost, variance) per subregion run; then (1.0, 0.0)
        self.starts: List[np.ndarray] = []
        self.radii: List[float] = []
        self.steps: List[int] = []
        self.offsets: List[Tuple[float, float, float]] = []
        self.mc = mc
        grid = mc._grid_gen
        orig = grid.get_near_identity_point

        def rec_grid(x: float, y: float, z: float) -> np.ndarray:
            self.offsets.append((x, y, z))
            return orig(x, y, z)

        grid.get_near_identity_point = rec_grid  # type: ignore[method-assign]

        def fake(init: np.ndarray, radius: float, n: int) -> Any:
            self.starts.append(np.array(init))
            self.radii.append(radius)
            self.steps.append(n)
            cost, var = self.script.pop(0) if self.script else (1.0, 0.0)
            return init, cost, var, SimpleNamespace(cost=cost)

        mc._zero_temp_with_variance = fake  # type: ignore[method-assign]


def _mc(seed: int = 0) -> MCOptimizer:
    return MCOptimizer(_Cost(), None, 0, np.random.default_rng(seed))  # type: ignore[arg-type]


def test_restart_is_about_initial_orientation_with_subregion_radius() -> None:
    box = math.radians(0.33)
    R0 = quaternion_to_matrix(matrix_to_quaternion(np.eye(3)))
    mc = _mc(1)
    # run 1 improves (cost 0.5 < 1.0): radius halves, continue from the improved state;
    # run 2 fails: radius doubles (capped by box), restart
    rec = _Recorder(mc, [(0.5, 1.0), (0.9, 1.0), (0.9, 1.0)])
    mc.variance_minimizing_optimize(R0, box, 30, 2, 0.0, 0.02**2)
    r0 = math.tan(box) / math.sqrt(48.0)
    assert math.isclose(rec.radii[0], r0)
    assert math.isclose(rec.radii[1], 0.5 * r0)  # shrunk after the improvement
    assert math.isclose(rec.radii[2], r0)  # doubled after the failure (no cap reached)
    # the restart offsets are uniform in +-(doubled radius), not +-box/2
    assert len(rec.offsets) >= 1
    for off in rec.offsets[:1]:
        assert max(abs(c) for c in off) <= r0 + 1e-15
    # the restart is applied to the INITIAL orientation (run 1 improved but kept the initial matrix
    # here, so also check against a global best that differs: see the next test)
    expect = quaternion_to_matrix(
        _quat_multiply(
            mc._grid_gen.get_near_identity_point(*rec.offsets[0]), matrix_to_quaternion(R0)
        )
    )
    np.testing.assert_allclose(rec.starts[2], expect, atol=1e-12)


def test_restart_ignores_improved_global_best() -> None:
    """After an improvement the next failure restarts about the initial orientation (C++), so the
    restart start must differ from the improved state by more than the offset alone."""
    box = math.radians(0.33)
    R0 = np.eye(3)
    mc = _mc(2)
    improved = quaternion_to_matrix(
        _quat_multiply(
            mc._grid_gen.get_near_identity_point(0.003, 0.0, 0.0), matrix_to_quaternion(R0)
        )
    )
    rec = _Recorder(mc, [])

    calls = {"n": 0}

    def fake(init: np.ndarray, radius: float, n: int) -> Any:
        rec.starts.append(np.array(init))
        calls["n"] += 1
        if calls["n"] == 1:
            return improved, 0.5, 1.0, SimpleNamespace(cost=0.5)  # improvement
        return init, 0.9, (1.0 if calls["n"] < 4 else 0.0), SimpleNamespace(cost=0.9)  # failure

    mc._zero_temp_with_variance = fake  # type: ignore[method-assign]
    mc.variance_minimizing_optimize(R0, box, 20, 2, 0.0, 0.02**2)
    assert calls["n"] >= 3
    off = rec.offsets[0]
    q0 = matrix_to_quaternion(R0)
    from_initial = quaternion_to_matrix(
        _quat_multiply(mc._grid_gen.get_near_identity_point(*off), q0)
    )
    np.testing.assert_allclose(rec.starts[2], from_initial, atol=1e-12)


def test_budget_extension_and_termination() -> None:
    """A run with variance above the threshold adds its steps to the budget; the loop ends after
    max_mc_steps worth of runs with variance below it (cost threshold 0 never fires)."""
    box = math.radians(0.33)
    mc = _mc(3)
    rec = _Recorder(mc, [])
    # all runs fail with low variance: ends after max_mc_steps / 10 = 20 runs (10 steps each)
    mc.variance_minimizing_optimize(np.eye(3), box, 200, 2, 0.0, 0.02**2)
    assert sum(rec.steps) == 200 and set(rec.steps) == {10}

    mc2 = _mc(3)
    rec2 = _Recorder(mc2, [(1.0, 1.0)] * 15)  # 15 high-variance runs first: budget +150
    mc2.variance_minimizing_optimize(np.eye(3), box, 200, 2, 0.0, 0.02**2)
    assert sum(rec2.steps) == 350
