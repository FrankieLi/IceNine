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
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

from .cost_functions import OverlapInfo, VoxelCostFunction
from .sampling import (
    QuaternionGrid,
    generate_local_grid,
    get_misorientation,
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


def _spacing_filter(
    candidates: List[SearchCandidate],
    angular_radius: float,
    symmetry_quats: np.ndarray,
) -> List[SearchCandidate]:
    """
    Filter candidates by angular spacing: reject a candidate if a closer
    candidate with better cost already exists.

    This matches the C++ Acceptable() logic in DiscreteSearch.h:244-258.
    The partition scheme mirrors the C++ swap-and-advance pattern from
    GetSpacedCandidates (DiscreteSearch.h:364-398).

    Args:
        candidates: List of SearchCandidate (with cost and orientation set)
        angular_radius: Minimum angular separation (radians)
        symmetry_quats: (N_sym, 4) symmetry operator quaternions

    Returns:
        Filtered list of accepted candidates
    """
    if len(candidates) <= 1:
        return candidates

    # Convert orientations to quaternions for misorientation check
    quats = np.array([matrix_to_quaternion(c.orientation) for c in candidates])

    # Partition: accepted candidates go to front, rejected stay at back.
    # C++ iterates pCur from element 1, comparing against [pFirstGood, end).
    # pFirstGood advances when a candidate is NOT acceptable (swap to front).
    # At the end, [pFirstGood, pCur) contains accepted candidates.
    #
    # Rewritten as a simple accept/reject list for clarity.
    accepted = [candidates[0]]
    accepted_quats = [quats[0]]

    for i in range(1, len(candidates)):
        # Check if acceptable: no existing accepted candidate within angular_radius
        # that also has better (lower) cost
        acceptable = True
        for j in range(len(accepted)):
            mis = get_misorientation(accepted_quats[j], quats[i], symmetry_quats)
            if mis < angular_radius and accepted[j].cost < candidates[i].cost:
                acceptable = False
                break
        if acceptable:
            accepted.append(candidates[i])
            accepted_quats.append(quats[i])

    return accepted


def get_symmetry_quaternions(symmetry) -> np.ndarray:
    """
    Extract proper rotation quaternions from a CrystalSymmetry object.

    Filters to proper rotations only (det > 0), converts to quaternions.

    Args:
        symmetry: CrystalSymmetry object

    Returns:
        (N_sym, 4) array of quaternions [w, x, y, z]
    """
    matrices = symmetry.get_rotation_matrices()
    quats = []
    for m in matrices:
        m = np.asarray(m, dtype=np.float64)
        if np.linalg.det(m) > 0:
            quats.append(matrix_to_quaternion(m))
    return np.array(quats)


def run_discrete_search_spaced(
    global_cost_fn: VoxelCostFunction,
    local_cost_fn: VoxelCostFunction,
    fz_orientations: np.ndarray,
    local_grid: np.ndarray,
    voxel_vertices,
    angular_radius: float,
    symmetry_quats: np.ndarray,
    phase_index: int = 0,
) -> List[SearchCandidate]:
    """
    Discrete search with per-clique angular spacing filter.

    For each FZ orientation (clique):
    1. Evaluate all local grid perturbations with global cost fn (pixel_radius=3)
    2. Re-evaluate accepted candidates with local cost fn (pixel_radius=0)
    3. Apply spacing filter to remove angular-near duplicates

    This matches C++ GetSpacedCandidates (DiscreteSearch.h:316-410) + the
    re-evaluation in RunDiscreteSearch (DiscreteAdaptive.tmpl.cpp:50-103).

    Args:
        global_cost_fn: Cost function for initial screening (pixel_radius=3)
        local_cost_fn: Cost function for re-evaluation (pixel_radius=0)
        fz_orientations: (N_fz, 3, 3) FZ orientation matrices
        local_grid: (N_local, 3, 3) local perturbation matrices
        voxel_vertices: Voxel triangle vertices
        angular_radius: Spacing filter radius (radians) = search diameter
        symmetry_quats: (N_sym, 4) crystal symmetry quaternions
        phase_index: Crystal phase index

    Returns:
        List of SearchCandidate (sorted by cost), filtered by spacing
    """
    all_candidates = []

    for fz_idx in range(len(fz_orientations)):
        search_center = fz_orientations[fz_idx]

        # Phase 1: Screen with global cost fn (pixel_radius=3)
        clique_candidates = []
        for lg_idx in range(len(local_grid)):
            candidate_orientation = local_grid[lg_idx] @ search_center
            overlap_info = global_cost_fn.evaluate(
                orientation=candidate_orientation,
                voxel_vertices=voxel_vertices,
                phase_index=phase_index,
            )
            if overlap_info.peak_overlap > 0:
                clique_candidates.append(SearchCandidate(
                    orientation=candidate_orientation,
                    cost=0.0,  # will be set by local cost fn below
                ))

        if not clique_candidates:
            continue

        # Phase 2: Re-evaluate with local cost fn (pixel_radius=0)
        # C++: GetSpacedCandidates lines 352-358, cost = 1 - GetConfidence
        for cand in clique_candidates:
            info = local_cost_fn.evaluate(
                orientation=cand.orientation,
                voxel_vertices=voxel_vertices,
                phase_index=phase_index,
            )
            cand.cost = 1.0 - info.confidence if info.peak_on_detector > 0 else 1.0
            cand.overlap_info = info

        # Phase 3: Apply spacing filter within this clique
        # C++: GetSpacedCandidates lines 364-398
        filtered = _spacing_filter(clique_candidates, angular_radius, symmetry_quats)
        all_candidates.extend(filtered)

    all_candidates.sort()
    return all_candidates


# ---------------------------------------------------------------------------
# SO(3) distance helper (quaternion-based, no symmetry reduction)
# ---------------------------------------------------------------------------

def _quat_misorientation_deg(q1: np.ndarray, q2: np.ndarray) -> float:
    """Geodesic distance in degrees between two orientations (no symmetry).

    Uses the half-angle formula: dist = 2 * arccos(|q1 · q2|).
    """
    dot = float(np.abs(np.dot(q1, q2)))
    dot = min(1.0, dot)
    return float(np.degrees(2.0 * np.arccos(dot)))


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

                    if trajectory is not None:
                        trajectory.append({
                            "step": step,
                            "event_type": "mc_accept",
                            "angular_step_deg": _quat_misorientation_deg(prev_best_q, best_q),
                            "cur_step_rad": cur_step,  # already halved
                        })
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
                    trajectory.append({
                        "step": step,
                        "event_type": "mc_restart",
                        "angular_step_deg": _quat_misorientation_deg(best_q, optimal_q),
                        "cur_step_rad": angular_step,  # reset to original step size
                    })

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

    def _zero_temp_with_variance(
        self,
        initial_orientation: np.ndarray,
        subregion_radius: float,
        n_steps: int,
    ) -> Tuple[np.ndarray, float, float, "OverlapInfo"]:
        """
        Run zero-temp MC within a subregion, tracking cost variance via Welford.

        Args:
            initial_orientation: Starting 3x3 rotation matrix
            subregion_radius: Angular radius of search neighborhood (radians)
            n_steps: Number of MC steps to run

        Returns:
            (best_orientation, best_cost, variance, best_overlap_info)

        C++ Reference: OrientationSearch.cpp:142-198 ZeroTemperatureOptimizationWithVariance
        """
        best_q = matrix_to_quaternion(initial_orientation)
        optimal_q = best_q.copy()

        best_info = self.cost_fn.evaluate(
            initial_orientation, self.voxel_vertices, self.phase_index
        )
        best_cost = best_info.cost
        current_cost = best_cost

        # Welford online variance tracking
        n_sampled = 1
        mean = best_cost
        m2 = 0.0  # sum of squared deviations

        radius = math.tan(subregion_radius) / math.sqrt(12.0) if subregion_radius > 0 else 0.01

        for _ in range(n_steps):
            x = self._rng.uniform(-radius, radius)
            y = self._rng.uniform(-radius, radius)
            z = self._rng.uniform(-radius, radius)

            delta_q = self._grid_gen.get_near_identity_point(x, y, z)
            trial_q = _quat_multiply(delta_q, optimal_q)
            trial_mat = quaternion_to_matrix(trial_q)

            trial_info = self.cost_fn.evaluate(
                trial_mat, self.voxel_vertices, self.phase_index
            )
            trial_cost = trial_info.cost

            # Welford update (C++ OrientationSearch.cpp:185-187)
            n_sampled += 1
            delta = trial_cost - mean
            mean = mean + delta / n_sampled
            m2 += delta * (trial_cost - mean)

            if trial_cost < current_cost:
                current_cost = trial_cost
                optimal_q = trial_q.copy()
                if trial_cost < best_cost:
                    best_cost = trial_cost
                    best_q = optimal_q.copy()
                    best_info = trial_info

        variance = m2 / (n_sampled - 1) if n_sampled > 1 else -1.0
        return quaternion_to_matrix(best_q), best_cost, variance, best_info

    def variance_minimizing_optimize(
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
                n_subregion_steps = int(
                    math.ceil(subregion_radius / cost_fn_angular_resolution) ** 2.7
                )
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
            if (global_min_cost < max_convergence_cost and
                    abs(variance) < convergence_variance):
                break

        return SearchCandidate(
            orientation=quaternion_to_matrix(global_best_q),
            cost=global_min_cost,
            overlap_info=global_best_info,
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
