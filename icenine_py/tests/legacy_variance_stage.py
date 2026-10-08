"""Frozen copy of MCOptimizer.variance_minimizing_optimize BEFORE the restart fix (restart offsets
+-box/2 about the current global best instead of C++'s +-SubregionRadius about the initial
orientation; see MIGRATION_HISTORY, "C++ vs Python variance stage"). Only for tests that reproduce
stored experiment results recorded with the old behaviour; use by monkeypatching
MCOptimizer.variance_minimizing_optimize. Do not use for new work.
"""

import math

import numpy as np

from icenine.orientation_search import MCOptimizer, SearchCandidate, _quat_multiply
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
