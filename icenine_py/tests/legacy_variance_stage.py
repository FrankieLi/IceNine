"""Frozen copies of MCOptimizer.variance_minimizing_optimize BEFORE the restart fix (restart
offsets +-box/2 about the current global best instead of C++'s +-SubregionRadius about the initial
orientation; see MIGRATION_HISTORY, "C++ vs Python variance stage") and of MCOptimizer.optimize
BEFORE the RandomRestartZeroTemp port (per-step loop, restart rule n_since >= min_ergodic, restart
about the global best in +-box/2; "C++-faithful MC: reruns"). Only for tests that reproduce stored
experiment results recorded with the old behaviour; use by monkeypatching
MCOptimizer.variance_minimizing_optimize / MCOptimizer.optimize. Do not use for new work.
"""

import math
from typing import Dict, List, Optional

import numpy as np

from icenine.orientation_search import (
    MCOptimizer,
    SearchCandidate,
    _quat_misorientation_deg,
    _quat_multiply,
)
from icenine.sampling import matrix_to_quaternion, quaternion_to_matrix

_ = MCOptimizer, np


def legacy_variance_minimizing_optimize(
    self,
    initial_orientation: np.ndarray,
    search_box_side: float,
    max_mc_steps: int,
    successive_restarts: int,
    max_convergence_cost: float,
    convergence_variance: float,
    cost_fn_angular_resolution: float = math.radians(0.5),
) -> SearchCandidate:
    """
    Adaptive sampling MC with variance-based convergence.

    Uses adaptive subregion radius: shrinks on improvement, expands on failure.
    Each subregion runs zero-temp MC with Welford variance tracking.
    Budget extends dynamically if variance indicates rough landscape.

    Args:
        initial_orientation: Starting 3x3 rotation matrix
        search_box_side: Side length of search region (radians)
        max_mc_steps: Maximum total MC steps budget
        successive_restarts: Max random restarts
        max_convergence_cost: Cost threshold for convergence
        convergence_variance: Variance threshold for convergence
        cost_fn_angular_resolution: Angular resolution for step count calc
                                    (default 0.5 degrees, hardcoded in C++)

    Returns:
        Best SearchCandidate found

    C++ Reference: OrientationSearch.cpp:206-287 AdaptiveSamplingZeroTemp
    """
    # Initialize
    global_best_q = matrix_to_quaternion(initial_orientation)
    current_q = global_best_q.copy()
    global_best_info = self.cost_fn.evaluate(
        initial_orientation, self.voxel_vertices, self.phase_index
    )
    global_min_cost = global_best_info.cost

    # C++: SubregionRadius = tan(search_box_side) / sqrt(48)
    subregion_radius = math.tan(search_box_side) / math.sqrt(48.0)
    total_steps_taken = 0
    max_steps = max_mc_steps
    n_global_restarts = 0

    while total_steps_taken < max_steps:
        # Adaptive step count per subregion: ceil(radius/resolution)^2.7
        # C++ OrientationSearch.cpp:227-228
        if cost_fn_angular_resolution > 0:
            n_subregion_steps = int(math.ceil(subregion_radius / cost_fn_angular_resolution) ** 2.7)
        else:
            n_subregion_steps = 10
        n_subregion_steps = max(n_subregion_steps, 10)

        # Run zero-temp MC with variance tracking in subregion
        current_mat = quaternion_to_matrix(current_q)
        new_orient, new_cost, variance, new_info = self._zero_temp_with_variance(
            current_mat, subregion_radius, n_subregion_steps
        )

        # Extend budget if landscape is rough
        if variance > convergence_variance:
            max_steps += n_subregion_steps

        total_steps_taken += n_subregion_steps

        # Decision: did we improve?
        if new_cost >= global_min_cost:
            # No improvement — expand and restart
            subregion_radius = min(2.0 * subregion_radius, search_box_side)
            # Random restart from global best
            half_box = search_box_side / 2.0
            rx = self._rng.uniform(-half_box, half_box)
            ry = self._rng.uniform(-half_box, half_box)
            rz = self._rng.uniform(-half_box, half_box)
            restart_q = self._grid_gen.get_near_identity_point(rx, ry, rz)
            current_q = _quat_multiply(restart_q, global_best_q)
            n_global_restarts += 1
        else:
            # Improvement — update global best and shrink
            global_min_cost = new_cost
            global_best_q = matrix_to_quaternion(new_orient)
            global_best_info = new_info
            current_q = global_best_q.copy()
            subregion_radius *= 0.5

        # Convergence check: cost AND variance both below threshold
        if global_min_cost < max_convergence_cost and abs(variance) < convergence_variance:
            break

    return SearchCandidate(
        orientation=quaternion_to_matrix(global_best_q),
        cost=global_min_cost,
        overlap_info=global_best_info,
    )


def legacy_mc_optimize(
    self,
    initial_orientation: np.ndarray,
    angular_box_side: float,
    angular_step: float,
    max_mc_steps: int,
    max_restarts: int,
    max_convergence_cost: float = 0.0,
    trajectory: Optional[List[Dict]] = None,
) -> SearchCandidate:
    """
    Run zero-temperature MC optimization from initial orientation.

    Args:
        initial_orientation: Starting 3x3 rotation matrix
        angular_box_side: Side length of search box (radians)
        angular_step: Initial angular step size (radians)
        max_mc_steps: Maximum MC steps
        max_restarts: Maximum number of random restarts
        max_convergence_cost: Early stop if cost drops below this
        trajectory: Optional list; when provided, accepted-move records are
            appended as dicts with keys step, event_type, angular_step_deg,
            cur_step_rad. Non-breaking: ignored when None (default).

    Returns:
        Best SearchCandidate found

    C++ Reference: OrientationSearch.cpp:296-365 RandomRestartZeroTemp
    """
    # Convert initial orientation to quaternion
    best_q = matrix_to_quaternion(initial_orientation)
    optimal_q = best_q.copy()
    prev_best_q = best_q.copy()  # for trajectory angular-step computation

    # Evaluate initial cost
    best_info = self.cost_fn.evaluate(initial_orientation, self.voxel_vertices, self.phase_index)
    global_min_cost = best_info.cost
    current_cost = global_min_cost

    cur_step = angular_step
    n_restarts = 0
    n_steps_since_improve = 0

    # Ergodic step count: estimate how many steps to cover the search box
    min_ergodic = (
        max(1, int(2.0 * (angular_box_side / cur_step) ** 3)) if cur_step > 0 else max_mc_steps
    )

    for step in range(max_mc_steps):
        # Generate random perturbation
        # C++: radius = tan(step_size) / sqrt(12)
        radius = math.tan(cur_step) / math.sqrt(12.0) if cur_step > 0 else 0.01
        x = self._rng.uniform(-radius, radius)
        y = self._rng.uniform(-radius, radius)
        z = self._rng.uniform(-radius, radius)

        delta_q = self._grid_gen.get_near_identity_point(x, y, z)

        # Compose: trial = delta * optimal
        trial_q = _quat_multiply(delta_q, optimal_q)
        trial_mat = quaternion_to_matrix(trial_q)

        # Evaluate cost
        trial_info = self.cost_fn.evaluate(trial_mat, self.voxel_vertices, self.phase_index)

        if trial_info.cost < current_cost:
            # Accept improvement
            current_cost = trial_info.cost
            optimal_q = trial_q.copy()

            if current_cost < global_min_cost:
                global_min_cost = current_cost
                best_q = optimal_q.copy()
                best_info = trial_info
                n_steps_since_improve = 0

                # Halve step size on improvement
                cur_step *= 0.5
                min_ergodic = (
                    max(1, int(2.0 * (angular_box_side / cur_step) ** 3))
                    if cur_step > 0
                    else max_mc_steps
                )

                if trajectory is not None:
                    trajectory.append(
                        {
                            "step": step,
                            "event_type": "mc_accept",
                            "angular_step_deg": _quat_misorientation_deg(prev_best_q, best_q),
                            "cur_step_rad": cur_step,  # already halved
                        }
                    )
                    prev_best_q = best_q.copy()

                # Early convergence check
                if global_min_cost < max_convergence_cost:
                    break
        else:
            n_steps_since_improve += 1

        # Restart check
        if n_steps_since_improve >= min_ergodic:
            n_restarts += 1
            if n_restarts > max_restarts:
                break

            # Random restart within box
            half_box = angular_box_side / 2.0
            rx = self._rng.uniform(-half_box, half_box)
            ry = self._rng.uniform(-half_box, half_box)
            rz = self._rng.uniform(-half_box, half_box)
            restart_q = self._grid_gen.get_near_identity_point(rx, ry, rz)
            restart_q = _quat_multiply(restart_q, best_q)
            optimal_q = restart_q.copy()

            if trajectory is not None:
                trajectory.append(
                    {
                        "step": step,
                        "event_type": "mc_restart",
                        "angular_step_deg": _quat_misorientation_deg(best_q, optimal_q),
                        "cur_step_rad": angular_step,  # reset to original step size
                    }
                )

            # Re-evaluate at restart point
            restart_mat = quaternion_to_matrix(optimal_q)
            restart_info = self.cost_fn.evaluate(restart_mat, self.voxel_vertices, self.phase_index)
            current_cost = restart_info.cost

            # Reset step size
            cur_step = angular_step
            n_steps_since_improve = 0
            min_ergodic = (
                max(1, int(2.0 * (angular_box_side / cur_step) ** 3))
                if cur_step > 0
                else max_mc_steps
            )

    return SearchCandidate(
        orientation=quaternion_to_matrix(best_q),
        cost=global_min_cost,
        overlap_info=best_info,
    )
