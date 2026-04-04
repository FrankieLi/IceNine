"""
Reconstruction orchestrator — serial voxel-by-voxel orientation recovery.

Ports the trivial (non-BFS, non-parallel) reconstruction pipeline from C++:
  SerialReconstruction → BasicVoxelReconstructor → multi-level adaptive search

C++ Reference:
    Src/Reconstructor.h/cpp        — BasicVoxelReconstructor, ReconstructVoxel
    Src/SerialReconstruction.h      — Main loop
    Src/ReconstructionSetup.h       — Setup and data loading
"""

import math
import time
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np
import torch

from .config_file import ConfigFile
from .cost_functions import OverlapInfo, VoxelCostFunction
from .crystal_structure import CrystalStructure
from .detector import Detector
from .experiment_setup import XDMExperimentSetup
from .experimental_data import ExperimentalData
from .forward_simulation import ForwardSimulation
from .mic_file import MicFile, ReconstructionState
from .orientation_search import (
    MCOptimizer,
    RiemannianAdamOptimizer,
    SearchCandidate,
    SearchParameters,
    get_symmetry_quaternions,
    hit_ratio_converged,
    run_discrete_search,
    run_discrete_search_spaced,
)
from .sample import Sample
from .sampling import (
    generate_local_grid,
    generate_local_grid_multi_level,
    load_fundamental_zone_file,
    matrix_to_quaternion,
    quaternion_to_matrix,
)
from .simulation import Simulation
from .simulation_range import SimulationRange


# ---------------------------------------------------------------------------
# Convergence codes
# ---------------------------------------------------------------------------

class ConvergenceCode:
    NOT_CONVERGED = 0
    HIT_RATIO_CONVERGED = 1
    COST_CONVERGED = 2
    MAX_LEVEL_REACHED = 3


# ---------------------------------------------------------------------------
# ReconstructionSetup — holds all loaded data
# ---------------------------------------------------------------------------

@dataclass
class ReconstructionSetup:
    """
    Container for all data needed during reconstruction.

    C++ Reference: Src/ReconstructionSetup.h
    """

    config: ConfigFile
    exp_setup: XDMExperimentSetup
    exp_data: ExperimentalData
    fz_orientations: np.ndarray  # (N_fz, 3, 3) rotation matrices
    search_params: SearchParameters
    simulator: Simulation
    detector_list: List[Detector]
    range_map: SimulationRange
    sample: Sample
    structure_list: List[CrystalStructure]
    diff_cost_fn: Optional[object] = None  # DifferentiableCostFunction, set for hybrid optimizer


def setup_reconstruction(
    config: ConfigFile,
    exp_data: Optional[ExperimentalData] = None,
    fz_orientations: Optional[np.ndarray] = None,
) -> ReconstructionSetup:
    """
    Initialize all reconstruction components from config.

    Args:
        config: ConfigFile with all reconstruction parameters
        exp_data: Pre-loaded experimental data (if None, loads from files)
        fz_orientations: Pre-loaded FZ orientations (if None, loads from file)

    Returns:
        ReconstructionSetup with all components initialized
    """
    # Initialize experiment
    exp_setup = XDMExperimentSetup(config)
    exp_setup.initialize_experiment()

    detector_list = exp_setup.get_detector_list()
    range_map = exp_setup.get_range_to_index_map()

    # Load experimental data
    if exp_data is None:
        exp_data = ExperimentalData.from_ascii_files(config, exp_setup)

    # Load FZ orientations
    if fz_orientations is None:
        fz_file = config.fundamental_zone_filename
        if fz_file:
            fz_orientations = load_fundamental_zone_file(fz_file)
        else:
            raise ValueError("No FundamentalZoneFilename in config and no fz_orientations provided")

    # Initialize sample
    sample = Sample()
    exp_setup.initialize_sample(sample, detector_list[0])
    structure_list = sample.get_structure_list()

    # Initialize simulator
    simulator = Simulation(exp_setup)

    # Search parameters
    search_params = SearchParameters.from_config(config)

    return ReconstructionSetup(
        config=config,
        exp_setup=exp_setup,
        exp_data=exp_data,
        fz_orientations=fz_orientations,
        search_params=search_params,
        simulator=simulator,
        detector_list=detector_list,
        range_map=range_map,
        sample=sample,
        structure_list=structure_list,
    )


def build_diff_cost_fn(
    setup: ReconstructionSetup,
    downsample_factors: Optional[List[int]] = None,
    omega_window: int = 1,
):
    """
    Build a DifferentiableCostFunction from an existing ReconstructionSetup.

    This is the standard way to enable the hybrid RiemannianAdamOptimizer.
    Assign the result to setup.diff_cost_fn and set
    setup.search_params.use_hybrid_optimizer = True.

    Args:
        setup: Fully initialized ReconstructionSetup
        downsample_factors: Downsample levels for MultiScaleImageStack.
            Default [1, 4, 8] → scale indices 0, 1, 2 in SearchParameters.adam_scale.
        omega_window: Omega integration window passed to SparseImageStack

    Returns:
        DifferentiableCostFunction ready for use with RiemannianAdamOptimizer
    """
    from .differentiable_cost import DifferentiableCostFunction, MultiScaleImageStack

    if downsample_factors is None:
        downsample_factors = [1, 4, 8]

    ms = MultiScaleImageStack(
        setup.exp_data.to_sparse_image_stack(),
        downsample_factors,
        omega_window=omega_window,
    )
    return DifferentiableCostFunction(
        simulator=setup.simulator,
        detector_list=setup.detector_list,
        range_map=setup.range_map,
        image_stack=ms,
        sample=setup.sample,
        structure_list=setup.structure_list,
        eta_limit=setup.exp_setup.get_eta_limit(),
    )


# ---------------------------------------------------------------------------
# BasicVoxelReconstructor — single voxel reconstruction
# ---------------------------------------------------------------------------

class BasicVoxelReconstructor:
    """
    Reconstruct a single voxel's crystal orientation.

    Multi-level adaptive search:
    1. For each resolution level (coarse → fine):
       a. Discrete search: FZ orientations × local grid
       b. Quick MC optimization (20 steps, no restarts)
       c. Keep top candidates
       d. Full MC optimization with convergence check
       e. If converged, stop

    C++ Reference: Src/Reconstructor.cpp ReconstructVoxel
    """

    def __init__(self, setup: ReconstructionSetup):
        self.setup = setup
        self.params = setup.search_params

        # Pre-generate local grids at all resolution levels
        self._local_grids = {}
        for level in range(self.params.min_local_resolution,
                           self.params.max_local_resolution + 1):
            self._local_grids[level] = generate_local_grid(
                self.params.local_grid_radius, level
            )

    def reconstruct_voxel(
        self,
        voxel_vertices: torch.Tensor,
        phase_index: int = 0,
        rng: Optional[np.random.Generator] = None,
    ) -> SearchCandidate:
        """
        Reconstruct orientation for a single voxel.

        Args:
            voxel_vertices: Triangle vertices in sample frame, shape (3, 3)
            phase_index: Crystal phase index
            rng: Random number generator for MC optimization

        Returns:
            Best SearchCandidate found

        C++ Reference: Reconstructor.cpp:181-245 ReconstructVoxel
        """
        # Two cost functions matching C++ two-tier approach:
        # - Global search (discrete): pixel_radius=3, wider cost landscape
        #   C++ Reference: DiscreteAdaptive.tmpl.cpp:65 nPixelRadius=3
        # - Local search (MC): pixel_radius=0, exact triangle rasterization
        #   C++ Reference: Reconstructor.cpp LocalSearchCostFunctions
        eta_limit = self.setup.exp_setup.get_eta_limit()
        global_cost_fn = VoxelCostFunction(
            simulator=self.setup.simulator,
            detector_list=self.setup.detector_list,
            range_map=self.setup.range_map,
            exp_data=self.setup.exp_data,
            sample=self.setup.sample,
            structure_list=self.setup.structure_list,
            mode='hard',
            eta_limit=eta_limit,
            pixel_radius=3,
            max_q=5.0,
        )
        local_cost_fn = VoxelCostFunction(
            simulator=self.setup.simulator,
            detector_list=self.setup.detector_list,
            range_map=self.setup.range_map,
            exp_data=self.setup.exp_data,
            sample=self.setup.sample,
            structure_list=self.setup.structure_list,
            mode='hard',
            eta_limit=eta_limit,
            pixel_radius=0,
        )

        mc_optimizer = MCOptimizer(
            cost_fn=local_cost_fn,
            voxel_vertices=voxel_vertices,
            phase_index=phase_index,
            rng=rng,
        )

        best_candidate = SearchCandidate(orientation=np.eye(3), cost=1.0)
        converged = False

        # Compute angular step size for MC
        # C++: box_width = local_grid_radius / 2^max_local_resolution
        box_width = self.params.local_grid_radius / (
            2 ** self.params.max_local_resolution
        )
        mc_step = box_width * self.params.mc_radius_scale_factor

        for level in range(self.params.min_local_resolution,
                           self.params.max_local_resolution + 1):
            local_grid = self._local_grids[level]
            t_level = time.time()

            # Phase 1: Discrete search — search ALL FZ orientations × local_grid
            # at every level, matching C++ ReconstructVoxel (Reconstructor.cpp:212-215).
            # The C++ does NOT do adaptive narrowing — it always searches the full
            # FZ set with progressively finer local grids.
            n_evals = len(self.setup.fz_orientations) * len(local_grid)
            print(f"    Level {level}: discrete search "
                  f"({len(self.setup.fz_orientations)} FZ × {len(local_grid)} local "
                  f"= {n_evals} evals)", flush=True)
            candidates = run_discrete_search(
                cost_fn=global_cost_fn,
                fz_orientations=self.setup.fz_orientations,
                local_grid=local_grid,
                voxel_vertices=voxel_vertices,
                phase_index=phase_index,
            )
            t_discrete = time.time() - t_level

            if not candidates:
                print(f"    Level {level}: no candidates found ({t_discrete:.1f}s)",
                      flush=True)
                continue

            print(f"    Level {level}: {len(candidates)} candidates ({t_discrete:.1f}s), "
                  f"best cost={candidates[0].cost:.4f}", flush=True)

            # Phase 2: Quick MC optimization (20 steps, no restarts)
            t_mc = time.time()
            n_quick = min(len(candidates), self.params.max_discrete_candidates)
            quick_candidates = []
            for cand in candidates[:n_quick]:
                result = mc_optimizer.optimize(
                    initial_orientation=cand.orientation,
                    angular_box_side=box_width,
                    angular_step=mc_step,
                    max_mc_steps=20,
                    max_restarts=0,
                    max_convergence_cost=self.params.max_convergence_cost,
                )
                quick_candidates.append(result)

            # Phase 3: Sort and keep top N
            quick_candidates.sort()
            top_candidates = quick_candidates[:self.params.max_discrete_candidates]
            t_quick = time.time() - t_mc
            print(f"    Level {level}: quick MC on {n_quick} candidates ({t_quick:.1f}s), "
                  f"best cost={top_candidates[0].cost:.4f}", flush=True)

            # Phase 4: Full MC optimization
            t_full = time.time()
            for ci, cand in enumerate(top_candidates):
                result = mc_optimizer.optimize(
                    initial_orientation=cand.orientation,
                    angular_box_side=box_width,
                    angular_step=mc_step,
                    max_mc_steps=self.params.max_mc_steps,
                    max_restarts=self.params.successive_restarts,
                    max_convergence_cost=self.params.max_convergence_cost,
                )

                if result.cost < best_candidate.cost:
                    best_candidate = result

                # Check convergence
                if (result.overlap_info is not None and
                        hit_ratio_converged(result.overlap_info,
                                            self.params.max_deepening_hit_ratio)):
                    converged = True
                    break

            t_full_mc = time.time() - t_full
            t_total_level = time.time() - t_level
            print(f"    Level {level}: full MC on {ci + 1} candidates ({t_full_mc:.1f}s), "
                  f"best cost={best_candidate.cost:.4f}, "
                  f"level total={t_total_level:.1f}s"
                  f"{' CONVERGED' if converged else ''}", flush=True)

            if converged:
                break

        return best_candidate


# ---------------------------------------------------------------------------
# AdaptiveVoxelReconstructor — single voxel, adaptive narrowing
# ---------------------------------------------------------------------------

class AdaptiveVoxelReconstructor:
    """
    Reconstruct a single voxel's crystal orientation using adaptive refinement.

    Multi-level adaptive search with candidate narrowing:
    1. For each resolution level:
       a. Generate local grid at fixed resolution (levels 0-1) with shrinking diameter
       b. Discrete search: FZ candidates × local grid (global cost fn, pixel_radius=3)
       c. Re-evaluate candidates with local cost fn (pixel_radius=0)
       d. Quick MC (10 steps, 5 restarts) on all candidates
       e. Sort, keep top 1/4 as FZ candidates for next level
       f. Shrink diameter by ÷1.5, increment nQMax
    2. After all levels:
       a. FindOptimal: full MC on top candidates with convergence check
       b. VarianceMinimizing: refine best until variance < 0.02²
       c. Final overlap evaluation

    C++ Reference: Src/DiscreteAdaptive.tmpl.cpp:108-250 ReconstructVoxel
    """

    def __init__(self, setup: ReconstructionSetup):
        self.setup = setup
        self.params = setup.search_params
        # Eval count tracking (set during reconstruct_voxel)
        self._last_global_evals = 0
        self._last_local_evals = 0

    @property
    def last_eval_counts(self) -> tuple:
        """Return (global_evals, local_evals, total_evals) from most recent reconstruct_voxel call."""
        total = self._last_global_evals + self._last_local_evals
        return (self._last_global_evals, self._last_local_evals, total)

    def reconstruct_voxel(
        self,
        voxel_vertices: torch.Tensor,
        phase_index: int = 0,
        rng: Optional[np.random.Generator] = None,
    ) -> SearchCandidate:
        """
        Reconstruct orientation for a single voxel using adaptive refinement.

        Args:
            voxel_vertices: Triangle vertices in sample frame, shape (3, 3)
            phase_index: Crystal phase index
            rng: Random number generator for MC optimization

        Returns:
            Best SearchCandidate found

        C++ Reference: DiscreteAdaptive.tmpl.cpp:108-250
        """
        eta_limit = self.setup.exp_setup.get_eta_limit()

        # Local cost function (pixel_radius=0) for MC and re-evaluation
        # C++ Reference: DiscreteAdaptive.tmpl.cpp:120-121 LocalSearchCostFunctions
        self._local_cost_fn = local_cost_fn = VoxelCostFunction(
            simulator=self.setup.simulator,
            detector_list=self.setup.detector_list,
            range_map=self.setup.range_map,
            exp_data=self.setup.exp_data,
            sample=self.setup.sample,
            structure_list=self.setup.structure_list,
            mode='hard',
            eta_limit=eta_limit,
            pixel_radius=0,
        )

        # MC optimizer uses local cost function (always built; used in VarianceMinimizing
        # and as fallback when use_hybrid_optimizer is False)
        mc_optimizer = MCOptimizer(
            cost_fn=local_cost_fn,
            voxel_vertices=voxel_vertices,
            phase_index=phase_index,
            rng=rng,
        )

        # Hybrid optimizer for FindOptimal phase (opt-in via SearchParameters)
        use_hybrid = (
            self.params.use_hybrid_optimizer
            and self.setup.diff_cost_fn is not None
        )
        if use_hybrid:
            find_optimizer: RiemannianAdamOptimizer = RiemannianAdamOptimizer(
                hard_cost_fn=local_cost_fn,
                diff_cost_fn=self.setup.diff_cost_fn,
                voxel_vertices=voxel_vertices,
                phase_index=phase_index,
                rng=rng,
            )

        # Crystal symmetry for spacing filter
        # C++ DiscreteAdaptive.tmpl.cpp:88
        symmetry = self.setup.exp_setup.get_sample_symmetry()
        symmetry_quats = get_symmetry_quaternions(symmetry) if symmetry else None

        # Initialize search state
        # C++ DiscreteAdaptive.tmpl.cpp:139-141
        fz_orientations = self.setup.fz_orientations  # full FZ set initially
        diameter = self.params.local_grid_radius
        n_q_max = 5.0 + self.params.min_local_resolution

        candidates = []
        total_global_evals = 0

        for level in range(self.params.max_local_resolution + 1):
            t_level = time.time()

            # Generate local grid: always levels 0-1 with current diameter
            # C++ DiscreteAdaptive.tmpl.cpp:149
            local_grid = generate_local_grid_multi_level(diameter, 0, 1)

            # Global cost function with current nQMax
            # C++ DiscreteAdaptive.tmpl.cpp:65-66 nPixelRadius=3
            global_cost_fn = VoxelCostFunction(
                simulator=self.setup.simulator,
                detector_list=self.setup.detector_list,
                range_map=self.setup.range_map,
                exp_data=self.setup.exp_data,
                sample=self.setup.sample,
                structure_list=self.setup.structure_list,
                mode='hard',
                eta_limit=eta_limit,
                pixel_radius=3,
                max_q=n_q_max,
            )

            # Phase 1: Discrete search
            n_evals = len(fz_orientations) * len(local_grid)
            print(f"    Level {level}: discrete search "
                  f"({len(fz_orientations)} FZ × {len(local_grid)} local "
                  f"= {n_evals} evals, nQMax={n_q_max:.0f}, "
                  f"diameter={math.degrees(diameter):.2f}°)", flush=True)

            candidates = run_discrete_search_spaced(
                global_cost_fn=global_cost_fn,
                local_cost_fn=local_cost_fn,
                fz_orientations=fz_orientations,
                local_grid=local_grid,
                voxel_vertices=voxel_vertices,
                angular_radius=diameter,
                symmetry_quats=symmetry_quats,
                phase_index=phase_index,
            )
            t_discrete = time.time() - t_level

            # Increment nQMax for next level
            # C++ DiscreteAdaptive.tmpl.cpp:155
            n_q_max += 1

            if not candidates:
                print(f"    Level {level}: no candidates found ({t_discrete:.1f}s)",
                      flush=True)
                continue

            print(f"    Level {level}: {len(candidates)} candidates "
                  f"(discrete={t_discrete:.1f}s), "
                  f"best cost={candidates[0].cost:.4f}", flush=True)

            # Phase 3: Quick MC on all candidates (10 steps, 5 restarts)
            # C++ DiscreteAdaptive.tmpl.cpp:173-194
            t_mc = time.time()
            trial_radius = max(diameter / 3.0, math.radians(0.2))
            # C++ ContinuousSearch.h:70-73 CalculateSearchParameter
            # BoxWidth = localGridRadius / 2^localResolution
            # For trial search: localResolution = min_local_resolution
            trial_box_width = trial_radius / (
                2 ** self.params.min_local_resolution
            )
            trial_step = trial_box_width * self.params.mc_radius_scale_factor

            for cand in candidates:
                result = mc_optimizer.optimize(
                    initial_orientation=cand.orientation,
                    angular_box_side=trial_box_width,
                    angular_step=trial_step,
                    max_mc_steps=10,
                    max_restarts=5,
                    max_convergence_cost=self.params.max_convergence_cost,
                )
                cand.orientation = result.orientation
                cand.cost = result.cost
                cand.overlap_info = result.overlap_info
            t_quick = time.time() - t_mc

            # Phase 4: Shrink diameter, keep top 1/4
            # C++ DiscreteAdaptive.tmpl.cpp:196-205
            diameter /= 1.5
            candidates.sort()
            n_keep = max(1, len(candidates) // 4)
            fz_orientations = np.array([c.orientation for c in candidates[:n_keep]])

            # Accumulate global cost fn evals for this level
            total_global_evals += global_cost_fn.eval_count

            t_total = time.time() - t_level
            print(f"    Level {level}: quick MC ({t_quick:.1f}s), "
                  f"kept {n_keep}/{len(candidates)}, "
                  f"best cost={candidates[0].cost:.4f}, "
                  f"level total={t_total:.1f}s", flush=True)

        # After all levels: FindOptimal + VarianceMinimizing
        # C++ DiscreteAdaptive.tmpl.cpp:210-246
        if not candidates:
            return SearchCandidate(orientation=np.eye(3), cost=1.0)

        final_radius = max(diameter / 3.0, math.radians(0.2))
        final_box_width = final_radius / (2 ** self.params.min_local_resolution)
        final_step = final_box_width * self.params.mc_radius_scale_factor

        # FindOptimal: full MC (or hybrid Adam) on top candidates with convergence check
        # C++ ContinuousSearch.h:316-346
        n_final = min(
            len(candidates), self.params.max_discrete_candidates
        )
        optimizer_label = "hybrid Adam" if use_hybrid else "full MC"
        print(f"    FindOptimal: {optimizer_label} on {n_final} candidates", flush=True)
        t_final = time.time()

        best_candidate = SearchCandidate(orientation=np.eye(3), cost=1.0)
        converged = False
        for ci, cand in enumerate(candidates[:n_final]):
            if use_hybrid:
                result = find_optimizer.optimize(
                    initial_orientation=cand.orientation,
                    angular_box_side=final_box_width,
                    n_steps=self.params.adam_n_steps,
                    lr=self.params.adam_lr,
                    scale=self.params.adam_scale,
                    max_restarts=self.params.successive_restarts,
                    max_convergence_cost=self.params.max_convergence_cost,
                )
            else:
                result = mc_optimizer.optimize(
                    initial_orientation=cand.orientation,
                    angular_box_side=final_box_width,
                    angular_step=final_step,
                    max_mc_steps=self.params.max_mc_steps,
                    max_restarts=self.params.successive_restarts,
                    max_convergence_cost=self.params.max_convergence_cost,
                )
            if result.cost < best_candidate.cost:
                best_candidate = result
                # Convergence check: hit_ratio >= 1.0
                # C++ ContinuousSearch.h:276-286 HitRatioConvergenceFn with ratio=1.0
                if (result.overlap_info is not None and
                        hit_ratio_converged(result.overlap_info, 1.0)):
                    converged = True
                    break
        t_find = time.time() - t_final
        print(f"    FindOptimal: {ci + 1} evaluated ({t_find:.1f}s), "
              f"best cost={best_candidate.cost:.4f}"
              f"{' CONVERGED' if converged else ''}", flush=True)

        # VarianceMinimizing: refine best until variance < 0.02²
        # C++ DiscreteAdaptive.tmpl.cpp:232
        t_var = time.time()
        var_result = mc_optimizer.variance_minimizing_optimize(
            initial_orientation=best_candidate.orientation,
            search_box_side=final_box_width,
            max_mc_steps=self.params.max_mc_steps,
            successive_restarts=self.params.successive_restarts,
            max_convergence_cost=0.0,  # C++ sets this to 0 for final optimization
            convergence_variance=0.02 ** 2,
        )
        if var_result.cost < best_candidate.cost:
            best_candidate = var_result
        t_var_elapsed = time.time() - t_var
        print(f"    VarianceMin: ({t_var_elapsed:.1f}s), "
              f"cost={best_candidate.cost:.4f}", flush=True)

        # Final overlap evaluation
        # C++ DiscreteAdaptive.tmpl.cpp:234-242
        final_info = local_cost_fn.evaluate(
            orientation=best_candidate.orientation,
            voxel_vertices=voxel_vertices,
            phase_index=phase_index,
        )
        best_candidate.overlap_info = final_info
        best_candidate.cost = final_info.cost

        # Store eval counts for benchmarking
        self._last_global_evals = total_global_evals
        self._last_local_evals = local_cost_fn.eval_count

        return best_candidate

    def evaluate_overlap(
        self,
        orientation: np.ndarray,
        voxel_vertices: torch.Tensor,
        phase_index: int = 0,
    ) -> OverlapInfo:
        """
        Evaluate overlap info for a given orientation (for BFS acceptance checks).

        C++ Reference: DiscreteAdaptive.tmpl.cpp:255-275 EvaluateOverlapInfo
        """
        eta_limit = self.setup.exp_setup.get_eta_limit()
        local_cost_fn = VoxelCostFunction(
            simulator=self.setup.simulator,
            detector_list=self.setup.detector_list,
            range_map=self.setup.range_map,
            exp_data=self.setup.exp_data,
            sample=self.setup.sample,
            structure_list=self.setup.structure_list,
            mode='hard',
            eta_limit=eta_limit,
            pixel_radius=0,
        )
        return local_cost_fn.evaluate(
            orientation=orientation,
            voxel_vertices=voxel_vertices,
            phase_index=phase_index,
        )

    def local_optimization(
        self,
        voxel_vertices: torch.Tensor,
        phase_index: int,
        initial_orientation: np.ndarray,
        rng: Optional[np.random.Generator] = None,
    ) -> SearchCandidate:
        """
        MC-only optimization from a given starting orientation (for BFS neighbors).

        Runs variance-minimizing MC without any discrete search. This is the
        cheap path used when orientation is inherited from a fitted neighbor.

        Args:
            voxel_vertices: Triangle vertices in sample frame, shape (3, 3)
            phase_index: Crystal phase index
            initial_orientation: Starting 3x3 rotation matrix (from neighbor)
            rng: Random number generator

        Returns:
            Optimized SearchCandidate

        C++ Reference: DiscreteAdaptive.tmpl.cpp:280-317 LocalOptimization
        """
        eta_limit = self.setup.exp_setup.get_eta_limit()
        local_cost_fn = VoxelCostFunction(
            simulator=self.setup.simulator,
            detector_list=self.setup.detector_list,
            range_map=self.setup.range_map,
            exp_data=self.setup.exp_data,
            sample=self.setup.sample,
            structure_list=self.setup.structure_list,
            mode='hard',
            eta_limit=eta_limit,
            pixel_radius=0,
        )

        mc_optimizer = MCOptimizer(
            cost_fn=local_cost_fn,
            voxel_vertices=voxel_vertices,
            phase_index=phase_index,
            rng=rng,
        )

        # C++ ContinuousSearch.h:70-73: BoxWidth = localGridRadius / 2^localResolution
        box_width = self.params.local_grid_radius / (
            2 ** self.params.min_local_resolution
        )

        # Variance-minimizing MC
        # C++ DiscreteAdaptive.tmpl.cpp:305
        result = mc_optimizer.variance_minimizing_optimize(
            initial_orientation=initial_orientation,
            search_box_side=box_width,
            max_mc_steps=self.params.max_mc_steps,
            successive_restarts=self.params.successive_restarts,
            max_convergence_cost=self.params.max_convergence_cost,
            convergence_variance=0.02 ** 2,
        )

        # Final overlap evaluation
        # C++ DiscreteAdaptive.tmpl.cpp:307-314
        final_info = local_cost_fn.evaluate(
            orientation=result.orientation,
            voxel_vertices=voxel_vertices,
            phase_index=phase_index,
        )
        result.overlap_info = final_info
        result.cost = final_info.cost

        return result


# ---------------------------------------------------------------------------
# BFSReconstruction — breadth-first spatial propagation
# ---------------------------------------------------------------------------

class BFSReconstruction:
    """
    Breadth-first reconstruction with spatial orientation propagation.

    Full adaptive search on seed voxels, then MC-only local optimization
    on neighbors using inherited orientations. This is much faster than
    independent per-voxel search for spatially coherent microstructures.

    Algorithm:
    1. Shuffle all voxels (randomized seed order)
    2. Pick next unvisited voxel as seed
    3. Full AdaptiveVoxelReconstructor.reconstruct_voxel() on seed
    4. If seed quality >= threshold: mark FITTED, propagate to neighbors via BFS
    5. BFS loop: pop neighbor, run local_optimization (MC-only),
       accept if (hit_ratio / best_hit_ratio) > 0.9
    6. Repeat from step 2 until all voxels visited

    C++ Reference:
        BreadthFirstReconstructor.tmpl.cpp:112-188 Fit()
        ReconstructionStrategies.tmpl.cpp:334-356 InsertSeed()
    """

    def __init__(self, setup: ReconstructionSetup):
        self.setup = setup
        self.reconstructor = AdaptiveVoxelReconstructor(setup)

    def reconstruct_sample(
        self,
        output_mic: Optional[str] = None,
        max_voxels: Optional[int] = None,
        rng: Optional[np.random.Generator] = None,
    ) -> List[int]:
        """
        Reconstruct all voxels using BFS propagation.

        Args:
            output_mic: Path to save reconstructed .mic file (optional)
            max_voxels: Limit number of voxels to process (for testing)
            rng: Random number generator

        Returns:
            List of voxel indices in order they were processed
        """
        if rng is None:
            rng = np.random.default_rng()

        mic = self.setup.sample.get_mic()
        n_total = len(mic.voxels)
        n_process = min(n_total, max_voxels) if max_voxels else n_total

        # Initialize all voxels to NOT_VISITED
        # C++ ReconstructionStrategies.tmpl.cpp:57
        for v in mic.voxels:
            v.reconstruction_id = ReconstructionState.NOT_VISITED

        # Randomized seed order (C++ uses random_shuffle)
        voxel_order = list(range(n_process))
        rng.shuffle(voxel_order)

        print(f"BFS Reconstruction: {n_process} voxels", flush=True)
        start_time = time.time()
        all_processed = []
        n_seeds = 0
        n_fitted = 0
        n_refit = 0

        for seed_idx in voxel_order:
            if mic.voxels[seed_idx].reconstruction_id != ReconstructionState.NOT_VISITED:
                continue

            n_seeds += 1
            t_seed = time.time()
            print(f"\n  Seed #{n_seeds}: voxel {seed_idx}", flush=True)

            processed = self._fit_from_seed(mic, seed_idx, rng)
            all_processed.extend(processed)

            seed_fitted = sum(
                1 for i in processed
                if mic.voxels[i].reconstruction_id == ReconstructionState.FITTED
            )
            seed_refit = len(processed) - seed_fitted
            n_fitted += seed_fitted
            n_refit += seed_refit

            t_elapsed = time.time() - t_seed
            print(f"  Seed #{n_seeds} done: {len(processed)} voxels "
                  f"({seed_fitted} fitted, {seed_refit} refit) in {t_elapsed:.1f}s",
                  flush=True)

        total_time = time.time() - start_time
        print(f"\nBFS complete: {n_seeds} seeds, {n_fitted} fitted, "
              f"{n_refit} refit, {total_time:.1f}s total", flush=True)

        # Save output .mic file
        if output_mic is not None:
            mic.write(output_mic)
            print(f"Saved reconstructed mic to {output_mic}", flush=True)

        return all_processed

    def _fit_from_seed(
        self,
        mic: MicFile,
        seed_idx: int,
        rng: np.random.Generator,
    ) -> List[int]:
        """
        Full reconstruction on seed, then BFS propagation to neighbors.

        C++ Reference: BreadthFirstReconstructor.tmpl.cpp:112-188 Fit()
        """
        voxel = mic.voxels[seed_idx]
        vertices = _get_voxel_vertices(voxel)

        # Full adaptive reconstruction on seed
        t0 = time.time()
        result = self.reconstructor.reconstruct_voxel(
            voxel_vertices=vertices,
            phase_index=voxel.phase,
            rng=rng,
        )
        t_recon = time.time() - t0

        # Evaluate overlap
        # C++ BreadthFirstReconstructor.tmpl.cpp:136
        overlap_info = self.reconstructor.evaluate_overlap(
            result.orientation, vertices, voxel.phase
        )

        # Compute confidence and hit_ratio
        # C++ CostFunctions.cpp: GetConfidence = peak_overlap / peak_on_detector
        # C++ CostFunctions.cpp: GetHitRatio = pixel_overlap / pixel_on_detector
        confidence = (overlap_info.peak_overlap / overlap_info.peak_on_detector
                      if overlap_info.peak_on_detector > 0 else 0.0)
        hit_ratio = (overlap_info.pixel_overlap / overlap_info.pixel_on_detector
                     if overlap_info.pixel_on_detector > 0 else 0.0)

        # Update seed voxel
        voxel.orientation = result.orientation
        voxel.cost = result.cost
        voxel.confidence = confidence
        voxel.overlap_ratio = hit_ratio

        print(f"    Seed voxel {seed_idx}: cost={result.cost:.4f}, "
              f"hit_ratio={hit_ratio:.3f}, conf={confidence:.3f} ({t_recon:.1f}s)",
              flush=True)

        # Check acceptance threshold
        # C++ BreadthFirstReconstructor.tmpl.cpp:142-149
        min_accel = self.setup.config.min_acceleration_threshold
        if hit_ratio < min_accel:
            voxel.reconstruction_id = ReconstructionState.REFIT
            print(f"    Seed rejected (hit_ratio {hit_ratio:.3f} < "
                  f"threshold {min_accel:.3f})", flush=True)
            return [seed_idx]

        # Mark fitted, start BFS
        # C++ BreadthFirstReconstructor.tmpl.cpp:153-156
        voxel.reconstruction_id = ReconstructionState.FITTED
        solution = [seed_idx]

        bfs_queue: deque[int] = deque()
        self._insert_seed(mic, seed_idx, bfs_queue)

        best_conf = hit_ratio
        n_bfs = 0

        # BFS expansion loop
        # C++ BreadthFirstReconstructor.tmpl.cpp:157-184
        while bfs_queue:
            neighbor_idx = bfs_queue.popleft()

            # Skip already-fitted voxels (C++ Pop() skips FITTED)
            if mic.voxels[neighbor_idx].reconstruction_id == ReconstructionState.FITTED:
                continue

            neighbor = mic.voxels[neighbor_idx]
            n_vertices = _get_voxel_vertices(neighbor)
            n_bfs += 1

            # MC-only optimization from inherited orientation
            t_local = time.time()
            opt_result = self.reconstructor.local_optimization(
                voxel_vertices=n_vertices,
                phase_index=neighbor.phase,
                initial_orientation=neighbor.orientation,
                rng=rng,
            )
            t_local_elapsed = time.time() - t_local

            # Compute hit_ratio from result
            n_info = opt_result.overlap_info
            n_hit_ratio = (n_info.pixel_overlap / n_info.pixel_on_detector
                           if n_info and n_info.pixel_on_detector > 0 else 0.0)
            n_confidence = (n_info.peak_overlap / n_info.peak_on_detector
                            if n_info and n_info.peak_on_detector > 0 else 0.0)

            # Update neighbor voxel
            neighbor.orientation = opt_result.orientation
            neighbor.cost = opt_result.cost
            neighbor.confidence = n_confidence
            neighbor.overlap_ratio = n_hit_ratio

            # Track best quality
            # C++ BreadthFirstReconstructor.tmpl.cpp:162
            best_conf = max(n_hit_ratio, best_conf)

            # Acceptance check: 90% of best quality
            # C++ BreadthFirstReconstructor.tmpl.cpp:163
            if best_conf > 0 and (n_hit_ratio / best_conf) > 0.9:
                # Accept: mark fitted, propagate to neighbors
                neighbor.reconstruction_id = ReconstructionState.FITTED
                solution.append(neighbor_idx)
                self._insert_seed(mic, neighbor_idx, bfs_queue)
                print(f"    BFS #{n_bfs} voxel {neighbor_idx}: FITTED "
                      f"cost={opt_result.cost:.4f}, hit_ratio={n_hit_ratio:.3f} "
                      f"({t_local_elapsed:.1f}s)", flush=True)
            else:
                # Reject: mark for later re-fitting
                neighbor.reconstruction_id = ReconstructionState.REFIT
                solution.append(neighbor_idx)
                print(f"    BFS #{n_bfs} voxel {neighbor_idx}: REFIT "
                      f"cost={opt_result.cost:.4f}, hit_ratio={n_hit_ratio:.3f} "
                      f"(ratio={n_hit_ratio / best_conf:.3f} < 0.9) "
                      f"({t_local_elapsed:.1f}s)", flush=True)

        return solution

    def _insert_seed(
        self,
        mic: MicFile,
        voxel_idx: int,
        queue: deque,
    ) -> None:
        """
        Propagate orientation to unvisited neighbors and add to BFS queue.

        C++ Reference: ReconstructionStrategies.tmpl.cpp:334-356 InsertSeed()
        """
        voxel = mic.voxels[voxel_idx]
        # C++ uses GetNeighbors with the solution grid
        # Python uses KDTree-based neighbor lookup with 2x side_length radius
        radius = 2.0 * voxel.side_length
        neighbors = mic.get_neighbors(voxel_idx, radius=radius)

        for n_idx in neighbors:
            neighbor = mic.voxels[n_idx]
            if neighbor.reconstruction_id == ReconstructionState.NOT_VISITED:
                neighbor.reconstruction_id = ReconstructionState.VISITED
                neighbor.orientation = voxel.orientation.copy()  # PROPAGATION
                queue.append(n_idx)


# ---------------------------------------------------------------------------
# SerialReconstruction — main loop over all voxels
# ---------------------------------------------------------------------------

class SerialReconstruction:
    """
    Reconstruct all voxels in a sample sequentially.

    C++ Reference: Src/SerialReconstruction.h
    """

    def __init__(self, setup: ReconstructionSetup):
        self.setup = setup
        self.reconstructor = BasicVoxelReconstructor(setup)
        self.results: List[SearchCandidate] = []

    def reconstruct_sample(
        self,
        output_mic: Optional[str] = None,
        max_voxels: Optional[int] = None,
        rng: Optional[np.random.Generator] = None,
    ) -> List[SearchCandidate]:
        """
        Reconstruct all voxels in the sample.

        Args:
            output_mic: Path to save reconstructed .mic file (optional)
            max_voxels: Limit number of voxels to reconstruct (for testing)
            rng: Random number generator

        Returns:
            List of SearchCandidate results, one per voxel
        """
        mic = self.setup.sample.get_mic()
        voxel_list = mic.voxels

        if max_voxels is not None:
            voxel_list = voxel_list[:max_voxels]

        n_voxels = len(voxel_list)
        print(f"Reconstructing {n_voxels} voxels...", flush=True)

        self.results = []
        start_time = time.time()

        for idx, voxel in enumerate(voxel_list):
            t_voxel = time.time()
            elapsed = t_voxel - start_time
            print(f"  Voxel {idx + 1}/{n_voxels} (phase={voxel.phase}, "
                  f"elapsed={elapsed:.1f}s)", flush=True)

            # Get voxel vertices
            vertices = _get_voxel_vertices(voxel)

            # Reconstruct
            result = self.reconstructor.reconstruct_voxel(
                voxel_vertices=vertices,
                phase_index=voxel.phase,
                rng=rng,
            )

            self.results.append(result)
            voxel_time = time.time() - t_voxel
            oi = result.overlap_info
            hit = (oi.pixel_overlap / oi.pixel_on_detector
                   if oi and oi.pixel_on_detector > 0 else 0)
            print(f"  Voxel {idx + 1}/{n_voxels} done: "
                  f"cost={result.cost:.4f}, hit_ratio={hit:.3f}, "
                  f"time={voxel_time:.1f}s", flush=True)

            # Update voxel orientation with reconstructed result
            voxel.orientation = result.orientation

        elapsed = time.time() - start_time
        print(f"Reconstruction complete: {n_voxels} voxels in {elapsed:.1f}s",
              flush=True)

        # Save output .mic file
        if output_mic is not None:
            mic.save(output_mic)
            print(f"Saved reconstructed mic to {output_mic}")

        return self.results


def _get_voxel_vertices(voxel) -> torch.Tensor:
    """Extract triangle vertices from a voxel in sample frame."""
    x, y, z = voxel.position
    s = voxel.side_length
    sqrt3_half = 0.5 * math.sqrt(3.0)

    if voxel.points_up:
        vertices = torch.tensor([
            [x, y, z],
            [x + s, y, z],
            [x + s * 0.5, y + s * sqrt3_half, z],
        ], dtype=torch.float32)
    else:
        vertices = torch.tensor([
            [x, y, z],
            [x + s * 0.5, y - s * sqrt3_half, z],
            [x + s, y, z],
        ], dtype=torch.float32)

    return vertices
