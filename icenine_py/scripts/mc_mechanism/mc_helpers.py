"""Pure helpers for the Phase B1 MC mechanism study (no pixel data needed).

mc_proposals draws trial rotations with exactly the proposal MCOptimizer.optimize uses; the other
helpers summarise a per-step trace and the improvement-probability curve.
"""

import math
from typing import Any, Dict, Optional, Tuple

import numpy as np

from icenine.orientation_search import (
    _quat_multiply,
    matrix_to_quaternion,
    quaternion_to_matrix,
)

# Event codes of a traced step
EV_NONE, EV_LOCAL, EV_GLOBAL, EV_RESTART, EV_EXHAUSTED, EV_NOT_RUN = 0, 1, 2, 3, 4, -1


def mc_proposals(
    R: np.ndarray, step_rad: float, n: int, rng: np.random.Generator, grid_gen: Any
) -> Tuple[np.ndarray, np.ndarray]:
    """n trial rotations around R, drawn exactly like MCOptimizer.optimize at step step_rad:
    radius = tan(step)/sqrt(12); x, y, z ~ U(-radius, radius) in that order;
    delta = grid_gen.get_near_identity_point(x, y, z); trial = delta * q(R).
    Returns (matrices (n, 3, 3), rotation angle of each delta in degrees (n,))."""
    q = matrix_to_quaternion(np.asarray(R))
    radius = math.tan(step_rad) / math.sqrt(12.0) if step_rad > 0 else 0.01
    mats = np.empty((n, 3, 3))
    ang = np.empty(n)
    for i in range(n):
        x = rng.uniform(-radius, radius)
        y = rng.uniform(-radius, radius)
        z = rng.uniform(-radius, radius)
        dq = grid_gen.get_near_identity_point(x, y, z)
        mats[i] = quaternion_to_matrix(_quat_multiply(dq, q))
        ang[i] = math.degrees(2.0 * math.acos(min(1.0, abs(float(dq[0])))))
    return mats, ang


def min_ergodic(box: float, step: float, max_mc_steps: int) -> int:
    """MCOptimizer's restart threshold 2 (box/step)^3 (at least 1)."""
    return max(1, int(2.0 * (box / step) ** 3)) if step > 0 else max_mc_steps


def expected_progress(
    c0: float, d0: float, c_new: np.ndarray, d_new: np.ndarray
) -> Dict[str, float]:
    """Per-proposal statistics at one point and step: p_improve = P(c_new < c0) (strict, as the
    optimizer's acceptance), cost_prog = E[max(0, c0 - c_new)], dist_prog = E[(d0 - d_new) 1{c_new
    < c0}] (signed: an accepted move can go away from the truth), n_improve."""
    c_new = np.asarray(c_new, float)
    d_new = np.asarray(d_new, float)
    acc = c_new < c0
    return dict(
        p_improve=float(acc.mean()),
        cost_prog=float(np.maximum(0.0, c0 - c_new).mean()),
        dist_prog=float(((d0 - d_new) * acc).mean()),
        n_improve=float(acc.sum()),
    )


def argmax_step(s_grid: np.ndarray, y: np.ndarray) -> Optional[float]:
    """The grid step whose value is largest; ties go to the smallest step (deterministic).
    None if all values are NaN or the maximum is not positive."""
    y = np.asarray(y, float)
    if not np.isfinite(y).any():
        return None
    m = np.nanmax(y)
    if m <= 0:
        return None
    return float(np.asarray(s_grid)[np.nonzero(y >= m)[0][0]])


def corner_angle_ratio() -> float:
    """Largest rotation angle of one MC proposal divided by the step (the cube corner
    x, y, z = +-tan(step)/sqrt(12) of get_near_identity_point; the ratio does not depend on the
    step at the sizes used). The median proposal is 0.97 of the step, the corner 1.71."""
    import itertools

    from icenine.orientation_search import QuaternionGrid

    g = QuaternionGrid()
    step = math.radians(0.1)
    r = math.tan(step) / math.sqrt(12.0)
    best = 0.0
    for sx, sy, sz in itertools.product((-1.0, 1.0), repeat=3):
        q = g.get_near_identity_point(sx * r, sy * r, sz * r)
        best = max(best, math.degrees(2.0 * math.acos(min(1.0, abs(float(q[0]))))))
    return best / 0.1


def restart_jump_angles(
    box_rad: float, n: int, rng: np.random.Generator, grid_gen: Any
) -> np.ndarray:
    """Rotation angles (deg) of n restart jumps, drawn as MCOptimizer.optimize draws them: x, y, z
    ~ U(-box/2, box/2) passed straight to get_near_identity_point (the inherited C++ convention:
    no tan/sqrt(12) scaling, unlike the MC proposal)."""
    half = box_rad / 2.0
    ang = np.empty(n)
    for i in range(n):
        rx = rng.uniform(-half, half)
        ry = rng.uniform(-half, half)
        rz = rng.uniform(-half, half)
        dq = grid_gen.get_near_identity_point(rx, ry, rz)
        ang[i] = math.degrees(2.0 * math.acos(min(1.0, abs(float(dq[0])))))
    return ang
