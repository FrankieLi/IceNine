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
from typing import Any, Dict, List, Optional, Tuple

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
# Optional geoopt dependency (for RiemannianAdamOptimizer)
# ---------------------------------------------------------------------------

_GEOOPT_AVAILABLE = False
try:
    import geoopt as _geoopt
    import torch as _torch

    _GEOOPT_AVAILABLE = True
except ImportError:
    pass


def _require_geoopt() -> None:
    if not _GEOOPT_AVAILABLE:
        raise RuntimeError(
            "geoopt is required for RiemannianAdamOptimizer. "
            "Install with: uv sync --extra riemannian"
        )


def _make_stiefel_param(R_init: np.ndarray) -> "_geoopt.ManifoldParameter":
    """Wrap a 3×3 rotation matrix as a geoopt Stiefel manifold parameter."""
    _require_geoopt()
    manifold = _geoopt.manifolds.Stiefel()
    return _geoopt.ManifoldParameter(_torch.from_numpy(R_init).float(), manifold=manifold)


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


LOCAL_OPTIMIZERS: Tuple[str, ...] = ("mc", "cma")


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
    # Hybrid optimizer fields (opt-in; default preserves MC-only behavior)
    use_hybrid_optimizer: bool = False
    adam_n_steps: int = 100
    adam_lr: float = 1e-4
    adam_scale: int = 2  # index into MultiScaleImageStack scales [1, 4, 8]
    # Local-refinement switch (opt-in; "mc" keeps the C++-parity MC + VarianceMinimizing path).
    # "cma" replaces the refinement MC (FindOptimal MC + VarianceMinimizing in
    # refine_from_candidates, the VarianceMinimizing call in local_optimization) with one local
    # CMA-ES run (CMAOptimizer) from the same start. Coarse search and quick MC are unchanged.
    local_optimizer: str = "mc"
    cma_sigma0_deg: float = 0.2  # initial CMA step (degrees of rotation vector)
    cma_max_evals: int = 1000  # total cost evaluations per CMA run, start included
    cma_popsize: Optional[int] = None  # None: the cma package default (4 + 3 ln 3 = 7)
    # BFS-only knobs (used by BFSReconstruction; never read by the no-start search).
    # cma_neighbor_max_evals: CMA budget of a BFS neighbour (and refit) local_optimization; seed
    #   voxels keep cma_max_evals. cma_retry_sigma0_deg: a rejected CMA neighbour is retried once
    #   from the same inherited start with this wider initial step (0 = no retry). Both act only
    #   when local_optimizer == "cma".
    # bfs_refit: opt-in port of the C++ LazyBFSClient::Refit pass over the voxels left REFIT.
    #   bfs_refit_conf: peak-overlap-ratio gate of that pass (None: the config's
    #   PartialResultAcceptanceConfidence).
    cma_neighbor_max_evals: int = 250
    cma_retry_sigma0_deg: float = 1.5
    bfs_refit: bool = False
    bfs_refit_conf: Optional[float] = None

    def __post_init__(self) -> None:
        if self.local_optimizer not in LOCAL_OPTIMIZERS:
            raise ValueError(
                f"local_optimizer must be one of {LOCAL_OPTIMIZERS}, got {self.local_optimizer!r}"
            )
        if not self.cma_sigma0_deg > 0:
            raise ValueError(f"cma_sigma0_deg must be > 0, got {self.cma_sigma0_deg}")
        if self.cma_max_evals < 2:
            raise ValueError(f"cma_max_evals must be >= 2, got {self.cma_max_evals}")
        if self.cma_popsize is not None and self.cma_popsize < 2:
            raise ValueError(f"cma_popsize must be >= 2 or None, got {self.cma_popsize}")
        if self.cma_neighbor_max_evals < 2:
            raise ValueError(
                f"cma_neighbor_max_evals must be >= 2, got {self.cma_neighbor_max_evals}"
            )
        if not self.cma_retry_sigma0_deg >= 0:
            raise ValueError(f"cma_retry_sigma0_deg must be >= 0, got {self.cma_retry_sigma0_deg}")
        if self.bfs_refit_conf is not None and not 0 <= self.bfs_refit_conf <= 1:
            raise ValueError(f"bfs_refit_conf must be in [0, 1] or None, got {self.bfs_refit_conf}")

    @classmethod
    def from_config(cls, config) -> "SearchParameters":
        """Create from ConfigFile."""
        refit_conf = getattr(config, "bfs_refit_conf", -1.0)  # < 0: unset
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
            # optional keys (LocalOptimizer, CMASigma0, CMAMaxEvals, CMAPopSize); the getattr
            # defaults are the "mc" path, so configs without the keys are unchanged
            local_optimizer=getattr(config, "local_optimizer", "mc"),
            cma_sigma0_deg=getattr(config, "cma_sigma0_deg", 0.2),
            cma_max_evals=getattr(config, "cma_max_evals", 1000),
            cma_popsize=(getattr(config, "cma_popsize", 0) or None),
            cma_neighbor_max_evals=getattr(config, "cma_neighbor_max_evals", 250),
            cma_retry_sigma0_deg=getattr(config, "cma_retry_sigma0_deg", 1.5),
            bfs_refit=bool(getattr(config, "bfs_refit", False)),
            bfs_refit_conf=(None if refit_conf < 0 else float(refit_conf)),
        )


@dataclass
class CMAResult(SearchCandidate):
    """SearchCandidate returned by CMAOptimizer, with the run's bookkeeping."""

    n_evals: int = 0  # cost evaluations used, the start included
    stop_reason: str = ""  # "max_evals", "converged_cost" or the cma package's stop key(s)


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
                candidates.append(
                    SearchCandidate(
                        orientation=candidate_orientation,
                        cost=overlap_info.cost,
                        overlap_info=overlap_info,
                    )
                )

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

    This matches the C++ Acceptable() logic in DiscreteSearch.h:244-258 and
    the swap-and-advance partition in GetSpacedCandidates (DiscreteSearch.h:
    364-398): each candidate is checked against the *entire remaining pool*
    [first_good, end) — not just previously-accepted candidates — so a
    lower-cost candidate later in the list can still reject an earlier one
    that hasn't been visited yet. Rejected candidates are swapped to the
    front; the surviving suffix [first_good, end) is the accepted set.

    Args:
        candidates: List of SearchCandidate (with cost and orientation set)
        angular_radius: Minimum angular separation (radians)
        symmetry_quats: (N_sym, 4) symmetry operator quaternions

    Returns:
        Filtered list of accepted candidates
    """
    if len(candidates) <= 1:
        return list(candidates)

    cands = list(candidates)
    quats = [matrix_to_quaternion(c.orientation) for c in cands]
    n = len(cands)

    first_good = 0
    cur = 1
    while cur < n:
        acceptable = True
        for j in range(first_good, n):
            if j == cur:
                continue
            mis = get_misorientation(quats[j], quats[cur], symmetry_quats)
            if mis < angular_radius and cands[j].cost < cands[cur].cost:
                acceptable = False
                break
        if not acceptable:
            cands[first_good], cands[cur] = cands[cur], cands[first_good]
            quats[first_good], quats[cur] = quats[cur], quats[first_good]
            first_good += 1
        cur += 1

    return cands[first_good:]


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
                clique_candidates.append(
                    SearchCandidate(
                        orientation=candidate_orientation,
                        cost=0.0,  # will be set by local cost fn below
                    )
                )

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

    @property
    def rng(self) -> np.random.Generator:
        """The generator driving this optimizer (shared with a CMAOptimizer seeded from it)."""
        return self._rng

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
                restart_info = self.cost_fn.evaluate(
                    restart_mat, self.voxel_vertices, self.phase_index
                )
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

            trial_info = self.cost_fn.evaluate(trial_mat, self.voxel_vertices, self.phase_index)
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
            if global_min_cost < max_convergence_cost and abs(variance) < convergence_variance:
                break

        return SearchCandidate(
            orientation=quaternion_to_matrix(global_best_q),
            cost=global_min_cost,
            overlap_info=global_best_info,
        )


# ---------------------------------------------------------------------------
# Local CMA-ES refinement
# ---------------------------------------------------------------------------


class CMAOptimizer:
    """
    Local CMA-ES refinement of an orientation (opt-in alternative to MC + VarianceMinimizing).

    The search variable is a rotation vector v (degrees) about the start R0: R = exp(v) R0, with
    v = 0 at the start and an isotropic initial step sigma0_deg. The start is evaluated first;
    the result is the lowest-cost orientation evaluated (first on ties), so it is never worse than
    the start. The evaluation budget is exact: no more than ``max_evals`` calls of the cost
    function (start included), and a generation is cut off, not completed, at the budget.

    Stopping: the cma package's flat-fitness, function-value and stagnation tests are switched
    off (the hard cost is quantised, so a small population is often flat or stalled); the run ends
    on the package's own x-tolerance (``tolx``, in degrees of v), on
    ``max_evals``, or when the cost falls below ``max_convergence_cost`` (0 disables it, the
    default; MC's default stop is not used because it would end the run well before the precision
    this optimizer reaches). These are the settings of the Phase B3 "cma_02" run.

    Random numbers: cma's seed is ``seed + 1`` (the package treats seed 0 as "random"). Pass
    ``seed`` to ``optimize`` for a fixed seed; otherwise one is drawn from ``rng`` (so a run is
    deterministic given the generator, as MCOptimizer is).

    Args:
        cost_fn: VoxelCostFunction (anything with ``evaluate(R, vertices, phase) -> info`` where
            ``info.cost`` is the scalar to minimise).
        voxel_vertices: Triangle vertices in sample frame, passed to cost_fn.
        phase_index: Crystal phase index.
        rng: Generator the per-run cma seeds are drawn from.
        sigma0_deg: Initial step (degrees). Default 0.2.
        max_evals: Evaluation budget per run, start included. Default 1000.
        popsize: Population size; None for the cma default (7 in 3 dimensions).
        tolx: cma x-tolerance stop. Default 1e-9.
        max_convergence_cost: Stop once the best cost is below this; 0 disables. Default 0.
    """

    def __init__(
        self,
        cost_fn: VoxelCostFunction,
        voxel_vertices: Any,
        phase_index: int = 0,
        rng: Optional[np.random.Generator] = None,
        sigma0_deg: float = 0.2,
        max_evals: int = 1000,
        popsize: Optional[int] = None,
        tolx: float = 1e-9,
        max_convergence_cost: float = 0.0,
    ):
        if not sigma0_deg > 0:
            raise ValueError(f"sigma0_deg must be > 0, got {sigma0_deg}")
        if max_evals < 1:
            raise ValueError(f"max_evals must be >= 1, got {max_evals}")
        self.cost_fn = cost_fn
        self.voxel_vertices = voxel_vertices
        self.phase_index = phase_index
        self._rng = rng if rng is not None else np.random.default_rng()
        self.sigma0_deg = float(sigma0_deg)
        self.max_evals = int(max_evals)
        self.popsize = None if popsize is None else int(popsize)
        self.tolx = float(tolx)
        self.max_convergence_cost = float(max_convergence_cost)

    def optimize(
        self,
        initial_orientation: np.ndarray,
        seed: Optional[int] = None,
    ) -> CMAResult:
        """
        Run CMA-ES from ``initial_orientation`` (3x3 rotation matrix).

        Args:
            initial_orientation: Starting rotation matrix.
            seed: Fixed seed (cma seed = seed + 1); None draws one from the generator.

        Returns:
            CMAResult: orientation and cost of the lowest-cost orientation evaluated, its
            overlap_info (from that evaluation: no extra call), n_evals and stop_reason.
        """
        import cma  # type: ignore[import-untyped]  # core dependency, imported lazily

        if seed is None:
            seed = int(self._rng.integers(0, 2**31 - 2))
        if seed < 0:
            raise ValueError(f"seed must be >= 0, got {seed}")
        R0 = np.asarray(initial_orientation, dtype=np.float64)

        n = 0
        best_info = self.cost_fn.evaluate(R0, self.voxel_vertices, self.phase_index)
        n += 1
        best_cost = float(best_info.cost)
        best_R = R0.copy()
        stop_reason = "max_evals"

        def converged() -> bool:
            return self.max_convergence_cost > 0 and best_cost < self.max_convergence_cost

        if converged():
            return CMAResult(best_R, best_cost, best_info, n_evals=n, stop_reason="converged_cost")
        if n >= self.max_evals:
            return CMAResult(best_R, best_cost, best_info, n_evals=n, stop_reason="max_evals")

        opts: Dict[str, Any] = dict(
            # cma would reseed the global legacy numpy RNG from an integer seed; NaN skips that and
            # the private RandomState below (seeded seed + 1, as the package would) is the stream
            seed=float("nan"),
            randn=np.random.RandomState(int(seed) + 1).randn,
            verbose=-9,
            tolfun=0.0,
            tolfunhist=0.0,
            tolflatfitness=10**9,
            tolstagnation=10**9,
            tolx=self.tolx,
        )
        if self.popsize is not None:
            opts["popsize"] = self.popsize
        es = cma.CMAEvolutionStrategy(np.zeros(3), self.sigma0_deg, opts)

        done = False
        while not done and not es.stop():
            X = es.ask()
            costs: List[float] = []
            for x in X:
                if n >= self.max_evals:
                    done = True
                    break
                R = Rotation.from_rotvec(np.radians(np.asarray(x))).as_matrix() @ R0
                info = self.cost_fn.evaluate(R, self.voxel_vertices, self.phase_index)
                n += 1
                c = float(info.cost)
                costs.append(c)
                if c < best_cost:
                    best_cost, best_R, best_info = c, R, info
                    if converged():
                        done, stop_reason = True, "converged_cost"
                        break
            if not done:
                es.tell(X, costs)
                if n >= self.max_evals:
                    done = True
        if stop_reason == "max_evals" and n < self.max_evals:
            stop_reason = "+".join(sorted(es.stop().keys())) or "cma_stop"
        return CMAResult(best_R, best_cost, best_info, n_evals=n, stop_reason=stop_reason)


# ---------------------------------------------------------------------------
# Hybrid Riemannian Adam + MC-restart optimizer
# ---------------------------------------------------------------------------


class RiemannianAdamOptimizer:
    """
    Hybrid optimizer: Riemannian Adam gradient descent with MC-style restarts.

    Uses geoopt's Stiefel manifold RiemannianAdam for fast gradient-based
    convergence near the basin, then applies random restarts (same strategy as
    MCOptimizer) when stuck. Hard VoxelCostFunction drives convergence decisions;
    DifferentiableCostFunction provides the gradient signal.

    Requires geoopt: uv sync --extra riemannian
    """

    _BETA1 = 0.9
    _BETA2 = 0.999

    def __init__(
        self,
        hard_cost_fn: VoxelCostFunction,
        diff_cost_fn,  # DifferentiableCostFunction
        voxel_vertices,
        phase_index: int = 0,
        rng: Optional[np.random.Generator] = None,
    ):
        _require_geoopt()
        self.hard_cost_fn = hard_cost_fn
        self.diff_cost_fn = diff_cost_fn
        self.voxel_vertices = voxel_vertices
        self.phase_index = phase_index
        self._grid_gen = QuaternionGrid()
        self._rng = rng or np.random.default_rng()

    def optimize(
        self,
        initial_orientation: np.ndarray,
        angular_box_side: float,
        n_steps: int = 100,
        lr: float = 1e-4,
        scale: int = 2,
        max_restarts: int = 2,
        max_convergence_cost: float = 0.0,
    ) -> "SearchCandidate":
        """
        Run hybrid Adam + MC-restart optimization from initial orientation.

        For each restart attempt:
        1. Initialize Stiefel manifold parameter from current orientation.
        2. Run n_steps of RiemannianAdam using differentiable cost (gradient signal).
        3. SVD re-orthogonalize the result (guards Stiefel float drift).
        4. Evaluate with hard cost function.
        5. If improved, update best; check convergence.
        6. If not converged and restarts remain, perturb within angular_box_side.

        Args:
            initial_orientation: Starting 3×3 rotation matrix
            angular_box_side: Side length of search box for restart perturbations (radians)
            n_steps: Number of Adam gradient steps per restart
            lr: Adam learning rate
            scale: Scale index into MultiScaleImageStack (0=full, 1=4×down, 2=8×down)
            max_restarts: Maximum number of MC-style restarts after Adam
            max_convergence_cost: Early stop if hard cost drops below this

        Returns:
            Best SearchCandidate found (same type as MCOptimizer.optimize)
        """
        best_orientation = initial_orientation.copy()
        best_info = self.hard_cost_fn.evaluate(
            best_orientation, self.voxel_vertices, self.phase_index
        )
        best_cost = best_info.cost
        current_orientation = initial_orientation.copy()

        for restart in range(max_restarts + 1):
            R = _make_stiefel_param(current_orientation)
            optimizer = _geoopt.optim.RiemannianAdam([R], lr=lr, betas=(self._BETA1, self._BETA2))

            for step in range(n_steps):
                optimizer.zero_grad()
                diff_info = self.diff_cost_fn.evaluate(
                    R,
                    self.voxel_vertices,
                    phase_index=self.phase_index,
                    scale=scale,
                )
                if diff_info.n_peaks == 0:
                    # No observable peaks at this orientation: diff_info.cost
                    # is a graph-connected zero (so .backward() doesn't raise)
                    # but its gradient is identically zero — further Adam
                    # steps here are wasted work, not a genuine convergence.
                    # Stop this restart's Adam loop now; the hard-cost restart
                    # logic below will perturb and try again.
                    break
                if diff_info.cost.requires_grad:
                    diff_info.cost.backward()
                    optimizer.step()

            # SVD re-orthogonalize (guards Stiefel float drift; cheap 3×3).
            # geoopt's Stiefel manifold is O(3), not SO(3), so U @ Vt can be a
            # reflection (det=-1); flip the last row of Vt to force det=+1.
            candidate_np = R.detach().numpy()
            U, _, Vt = np.linalg.svd(candidate_np)
            if np.linalg.det(U @ Vt) < 0:
                Vt[-1] *= -1
            candidate_np = U @ Vt

            hard_info = self.hard_cost_fn.evaluate(
                candidate_np, self.voxel_vertices, self.phase_index
            )
            if hard_info.cost < best_cost:
                best_cost = hard_info.cost
                best_orientation = candidate_np.copy()
                best_info = hard_info

            if best_cost < max_convergence_cost:
                break

            if restart < max_restarts:
                half_box = angular_box_side / 2.0
                rx = self._rng.uniform(-half_box, half_box)
                ry = self._rng.uniform(-half_box, half_box)
                rz = self._rng.uniform(-half_box, half_box)
                perturb_q = self._grid_gen.get_near_identity_point(rx, ry, rz)
                best_q = matrix_to_quaternion(best_orientation)
                current_orientation = quaternion_to_matrix(_quat_multiply(perturb_q, best_q))

        return SearchCandidate(
            orientation=best_orientation,
            cost=best_cost,
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
