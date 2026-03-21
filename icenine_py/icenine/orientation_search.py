"""
Orientation search algorithms for reconstruction.

Implements discrete grid search and zero-temperature Monte Carlo optimization
for finding crystal orientations that best match experimental diffraction data.

C++ Reference:
    Src/DiscreteSearch.h/cpp   — Discrete grid search
    Src/ContinuousSearch.h     — MC optimization wrapper
    Src/OrientationSearch.cpp  — Zero-temp MC optimizer
    Src/Reconstructor.cpp      — Multi-level reconstruction loop
"""

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

from .cost_functions import OverlapInfo, VoxelCostFunction
from .sampling import (
    QuaternionGrid,
    generate_local_grid,
    matrix_to_quaternion,
    quaternion_to_matrix,
    _quat_multiply,
)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SearchCandidate:
    """
    A candidate orientation with associated cost.

    C++ Reference: SearchDetails.h SCandidate
    """

    orientation: np.ndarray  # 3x3 rotation matrix
    cost: float = 1.0
    overlap_info: Optional[OverlapInfo] = None

    def __lt__(self, other: "SearchCandidate") -> bool:
        return self.cost < other.cost


@dataclass
class SearchParameters:
    """
    Parameters controlling the orientation search.

    C++ Reference: SearchDetails.h SSearchParameter
    """

    local_grid_radius: float = math.radians(5.0)  # radians
    min_local_resolution: int = 0
    max_local_resolution: int = 3
    max_discrete_candidates: int = 100
    max_mc_steps: int = 3500
    mc_radius_scale_factor: float = 1.0
    successive_restarts: int = 2
    max_convergence_cost: float = 0.1
    max_deepening_hit_ratio: float = 0.7
    max_accepted_cost: float = 0.9

    @classmethod
    def from_config(cls, config) -> "SearchParameters":
        """Create from ConfigFile."""
        return cls(
            local_grid_radius=config.local_orientation_grid_radius,
            min_local_resolution=config.min_local_resolution,
            max_local_resolution=config.max_local_resolution,
            max_discrete_candidates=config.max_discrete_candidates,
            max_mc_steps=config.max_mc_steps,
            mc_radius_scale_factor=config.mc_radius_scale_factor,
            successive_restarts=config.successive_restarts,
            max_convergence_cost=config.max_convergence_cost,
            max_deepening_hit_ratio=config.max_deepening_hit_ratio,
            max_accepted_cost=config.max_accepted_cost,
        )


# ---------------------------------------------------------------------------
# Discrete search
# ---------------------------------------------------------------------------

def run_discrete_search(
    cost_fn: VoxelCostFunction,
    fz_orientations: np.ndarray,
    local_grid: np.ndarray,
    voxel_vertices,
    phase_index: int = 0,
) -> List[SearchCandidate]:
    """
    Evaluate all FZ orientation × local grid combinations.

    For each (FZ_orientation, local_perturbation):
        candidate = local_perturbation @ FZ_orientation
        cost = cost_fn.evaluate(candidate)
        if cost < 1.0: keep as candidate

    Args:
        cost_fn: VoxelCostFunction for evaluation
        fz_orientations: (N_fz, 3, 3) array of FZ orientation matrices
        local_grid: (N_local, 3, 3) array of local perturbation matrices
        voxel_vertices: Voxel triangle vertices (torch.Tensor, shape (3, 3))
        phase_index: Crystal phase index

    Returns:
        List of SearchCandidate with cost < 1.0, sorted by cost

    C++ Reference: DiscreteSearch.h LocalAngularInterpProcess
    """
    candidates = []

    for fz_idx in range(len(fz_orientations)):
        search_center = fz_orientations[fz_idx]

        for lg_idx in range(len(local_grid)):
            # Compose: candidate = local_grid @ search_center
            candidate_orientation = local_grid[lg_idx] @ search_center

            overlap_info = cost_fn.evaluate(
                orientation=candidate_orientation,
                voxel_vertices=voxel_vertices,
                phase_index=phase_index,
            )

            if overlap_info.peak_overlap > 0:
                candidates.append(SearchCandidate(
                    orientation=candidate_orientation,
                    cost=overlap_info.cost,
                    overlap_info=overlap_info,
                ))

    candidates.sort()
    return candidates


# ---------------------------------------------------------------------------
# Zero-temperature Monte Carlo optimizer
# ---------------------------------------------------------------------------

class MCOptimizer:
    """
    Zero-temperature Monte Carlo optimization for orientation refinement.

    Greedy descent with random restarts: accepts perturbations only if they
    reduce cost. Halves step size on improvement, restarts from random
    position if stuck.

    C++ Reference: OrientationSearch.cpp RandomRestartZeroTemp
    """

    def __init__(
        self,
        cost_fn: VoxelCostFunction,
        voxel_vertices,
        phase_index: int = 0,
        rng: Optional[np.random.Generator] = None,
    ):
        self.cost_fn = cost_fn
        self.voxel_vertices = voxel_vertices
        self.phase_index = phase_index
        self._grid_gen = QuaternionGrid()
        self._rng = rng or np.random.default_rng()

    def optimize(
        self,
        initial_orientation: np.ndarray,
        angular_box_side: float,
        angular_step: float,
        max_mc_steps: int,
        max_restarts: int,
        max_convergence_cost: float = 0.0,
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

        Returns:
            Best SearchCandidate found

        C++ Reference: OrientationSearch.cpp:296-365 RandomRestartZeroTemp
        """
        # Convert initial orientation to quaternion
        best_q = matrix_to_quaternion(initial_orientation)
        optimal_q = best_q.copy()

        # Evaluate initial cost
        best_info = self.cost_fn.evaluate(
            initial_orientation, self.voxel_vertices, self.phase_index
        )
        global_min_cost = best_info.cost
        current_cost = global_min_cost

        cur_step = angular_step
        n_restarts = 0
        n_steps_since_improve = 0

        # Ergodic step count: estimate how many steps to cover the search box
        min_ergodic = max(
            1, int(2.0 * (angular_box_side / cur_step) ** 3)
        ) if cur_step > 0 else max_mc_steps

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
            trial_info = self.cost_fn.evaluate(
                trial_mat, self.voxel_vertices, self.phase_index
            )

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
                    min_ergodic = max(
                        1, int(2.0 * (angular_box_side / cur_step) ** 3)
                    ) if cur_step > 0 else max_mc_steps

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

                # Re-evaluate at restart point
                restart_mat = quaternion_to_matrix(optimal_q)
                restart_info = self.cost_fn.evaluate(
                    restart_mat, self.voxel_vertices, self.phase_index
                )
                current_cost = restart_info.cost

                # Reset step size
                cur_step = angular_step
                n_steps_since_improve = 0
                min_ergodic = max(
                    1, int(2.0 * (angular_box_side / cur_step) ** 3)
                ) if cur_step > 0 else max_mc_steps

        return SearchCandidate(
            orientation=quaternion_to_matrix(best_q),
            cost=global_min_cost,
            overlap_info=best_info,
        )


# ---------------------------------------------------------------------------
# Convergence checks
# ---------------------------------------------------------------------------

def hit_ratio_converged(
    overlap_info: OverlapInfo,
    threshold: float,
) -> bool:
    """
    Check if hit ratio exceeds threshold.

    C++ Reference: ContinuousSearch.h HitRatioConvergenceFn
    """
    return overlap_info.hit_ratio >= threshold
