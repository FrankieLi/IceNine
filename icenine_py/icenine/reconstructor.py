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
from typing import Any, Callable, Dict, List, Optional, Tuple, TYPE_CHECKING

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
    CMAOptimizer,
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

if TYPE_CHECKING:
    from .differentiable_cost import DifferentiableCostFunction


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
    diff_cost_fn: Optional["DifferentiableCostFunction"] = None  # set for hybrid optimizer


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
) -> "DifferentiableCostFunction":
    """
    Build a DifferentiableCostFunction from an existing ReconstructionSetup.

    This is the standard way to enable the hybrid RiemannianAdamOptimizer.
    Assign the result to setup.diff_cost_fn and set
    setup.search_params.use_hybrid_optimizer = True.

    Args:
        setup: Fully initialized ReconstructionSetup
        downsample_factors: Downsample levels for MultiScaleImageStack.
            Default [1, 4, 8] → scale indices 0, 1, 2 in SearchParameters.adam_scale.
        omega_window: Omega integration window passed to MultiScaleImageStack

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
        for level in range(self.params.min_local_resolution, self.params.max_local_resolution + 1):
            self._local_grids[level] = generate_local_grid(self.params.local_grid_radius, level)

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
            mode="hard",
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
            mode="hard",
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
        box_width = self.params.local_grid_radius / (2**self.params.max_local_resolution)
        mc_step = box_width * self.params.mc_radius_scale_factor

        for level in range(self.params.min_local_resolution, self.params.max_local_resolution + 1):
            local_grid = self._local_grids[level]
            t_level = time.time()

            # Phase 1: Discrete search — search ALL FZ orientations × local_grid
            # at every level, matching C++ ReconstructVoxel (Reconstructor.cpp:212-215).
            # The C++ does NOT do adaptive narrowing — it always searches the full
            # FZ set with progressively finer local grids.
            n_evals = len(self.setup.fz_orientations) * len(local_grid)
            print(
                f"    Level {level}: discrete search "
                f"({len(self.setup.fz_orientations)} FZ × {len(local_grid)} local "
                f"= {n_evals} evals)",
                flush=True,
            )
            candidates = run_discrete_search(
                cost_fn=global_cost_fn,
                fz_orientations=self.setup.fz_orientations,
                local_grid=local_grid,
                voxel_vertices=voxel_vertices,
                phase_index=phase_index,
            )
            t_discrete = time.time() - t_level

            if not candidates:
                print(f"    Level {level}: no candidates found ({t_discrete:.1f}s)", flush=True)
                continue

            print(
                f"    Level {level}: {len(candidates)} candidates ({t_discrete:.1f}s), "
                f"best cost={candidates[0].cost:.4f}",
                flush=True,
            )

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
            top_candidates = quick_candidates[: self.params.max_discrete_candidates]
            t_quick = time.time() - t_mc
            print(
                f"    Level {level}: quick MC on {n_quick} candidates ({t_quick:.1f}s), "
                f"best cost={top_candidates[0].cost:.4f}",
                flush=True,
            )

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
                if result.overlap_info is not None and hit_ratio_converged(
                    result.overlap_info, self.params.max_deepening_hit_ratio
                ):
                    converged = True
                    break

            t_full_mc = time.time() - t_full
            t_total_level = time.time() - t_level
            print(
                f"    Level {level}: full MC on {ci + 1} candidates ({t_full_mc:.1f}s), "
                f"best cost={best_candidate.cost:.4f}, "
                f"level total={t_total_level:.1f}s"
                f"{' CONVERGED' if converged else ''}",
                flush=True,
            )

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

    def __init__(self, setup: ReconstructionSetup, min_sin_eta: float = 0.0):
        """min_sin_eta: optional lower bound on |sin eta| of the diffracted beam, applied to every
        cost function reconstruct_voxel builds (default 0 = off, the C++ behaviour)."""
        self.setup = setup
        self.params = setup.search_params
        self.min_sin_eta = min_sin_eta
        # Eval count tracking (set during reconstruct_voxel)
        self._last_global_evals = 0
        self._last_local_evals = 0
        # evaluations of the most recent local_optimization call (BFS provenance)
        self.last_local_optimization_evals = 0
        # Diagnostics of the most recent call: per-level best candidate (orientation, cost) after
        # the quick MC, and FindOptimal's winner (candidate index, converged flag)
        self.last_level_best: List[Tuple[int, np.ndarray, float]] = []
        self.last_find_optimal: Dict[str, Any] = {}
        # Optional non-invasive recorder: called as recorder(event, data) with copies of the
        # search state ("discrete", "quick_mc", "find_candidate", "find_final", "variance",
        # "final"). It never touches the random stream, so results are bit-identical with or
        # without it. Default None = off.
        self.recorder: Optional[Callable[[str, Dict[str, Any]], None]] = None
        # Search knobs probed by scripts/findoptimal_robustness (defaults = the C++ behaviour)
        self.keep_fraction: float = 0.25  # fraction of candidates kept per level
        # also keep the top keep_fraction by the discrete-stage score (union with the post-MC rank)
        self.keep_union_discrete: bool = False
        self.n_q_start_offset: float = 0.0  # added to the initial n_q_max (5 + min resolution)
        self.global_pixel_radius: int = 3  # pixel radius of the coarse (global) cost function
        # Optional hook: extra_candidates(level, candidates) -> extra SearchCandidates (orientation
        # only; the quick MC fills their cost) added before the quick MC of that level (F1b).
        # Default off; see MIGRATION_HISTORY 'FindOptimal robustness'.
        self.extra_candidates: Optional[
            Callable[[int, List[SearchCandidate]], List[SearchCandidate]]
        ] = None
        # Optional hook: rank_key(level, candidates) -> array, lower = better, replaces the
        # post-quick-MC local cost as the sort key used for pruning (and the hand-off order)
        self.rank_key: Optional[Callable[[int, List[SearchCandidate]], np.ndarray]] = None

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
        self._check_local_optimizer()
        eta_limit = self.setup.exp_setup.get_eta_limit()

        # Local cost function (pixel_radius=0) for MC and re-evaluation
        # C++ Reference: DiscreteAdaptive.tmpl.cpp:120-121 LocalSearchCostFunctions
        self._local_cost_fn = local_cost_fn = self._make_local_cost_fn()
        mc_optimizer, find_optimizer = self._make_optimizers(
            local_cost_fn, voxel_vertices, phase_index, rng
        )

        # Crystal symmetry for spacing filter
        # C++ DiscreteAdaptive.tmpl.cpp:88
        symmetry = self.setup.exp_setup.get_sample_symmetry()
        # Empty (not None) when symmetry is NONE: reduce_to_fundamental_zone's
        # loop over symmetry_quats then no-ops, correctly leaving q unreduced.
        symmetry_quats = (
            get_symmetry_quaternions(symmetry) if symmetry is not None else np.zeros((0, 4))
        )

        # Initialize search state
        # C++ DiscreteAdaptive.tmpl.cpp:139-141
        fz_orientations = self.setup.fz_orientations  # full FZ set initially
        diameter = self.params.local_grid_radius
        n_q_max = 5.0 + self.params.min_local_resolution + self.n_q_start_offset

        candidates = []
        total_global_evals = 0
        self.last_level_best = []
        self.last_find_optimal = {}

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
                mode="hard",
                eta_limit=eta_limit,
                pixel_radius=self.global_pixel_radius,
                max_q=n_q_max,
                min_sin_eta=self.min_sin_eta,
            )

            # Phase 1: Discrete search
            n_evals = len(fz_orientations) * len(local_grid)
            print(
                f"    Level {level}: discrete search "
                f"({len(fz_orientations)} FZ × {len(local_grid)} local "
                f"= {n_evals} evals, nQMax={n_q_max:.0f}, "
                f"diameter={math.degrees(diameter):.2f}°)",
                flush=True,
            )

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
                print(f"    Level {level}: no candidates found ({t_discrete:.1f}s)", flush=True)
                continue

            # Note: candidates[0].cost here is 1 - confidence (peak-ratio metric
            # set in run_discrete_search_spaced), not 1 - quality (pixel/detector
            # ratio) used by the "best cost=" prints below once Phase 3 MC
            # overwrites .cost — the two are different scales, not comparable.
            print(
                f"    Level {level}: {len(candidates)} candidates "
                f"(discrete={t_discrete:.1f}s), "
                f"best cost={candidates[0].cost:.4f}",
                flush=True,
            )

            if self.extra_candidates is not None:
                # NOTE: these bypass the spacing filter of the discrete stage and enlarge the
                # list that n_keep is computed from below
                candidates = list(candidates) + list(self.extra_candidates(level, candidates))
            if self.recorder is not None:
                self.recorder(
                    "discrete",
                    dict(
                        level=level,
                        n_q_max=n_q_max - 1,
                        diameter=diameter,
                        R=np.stack([c.orientation for c in candidates]).copy(),
                        score=np.array([c.cost for c in candidates]),
                    ),
                )

            # Phase 3: Quick MC on all candidates (10 steps, 5 restarts)
            # C++ DiscreteAdaptive.tmpl.cpp:173-194
            t_mc = time.time()
            ids = {id(c): i for i, c in enumerate(candidates)}
            disc_scores = np.array([c.cost for c in candidates])
            trial_radius = max(diameter / 3.0, math.radians(0.2))
            # C++ ContinuousSearch.h:70-73 CalculateSearchParameter
            # BoxWidth = localGridRadius / 2^localResolution
            # For trial search: localResolution = min_local_resolution
            trial_box_width = trial_radius / (2**self.params.min_local_resolution)
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
            pre_sort = list(candidates)  # discrete-stage order (the quick MC mutated in place)
            if self.rank_key is not None:
                key = np.asarray(self.rank_key(level, candidates))
                assert len(key) == len(candidates), "rank_key must return one key per candidate"
                order = np.argsort(key, kind="stable")
                candidates = [candidates[i] for i in order]
            else:
                candidates.sort()
            self.last_level_best.append(
                (level, candidates[0].orientation.copy(), float(candidates[0].cost))
            )
            n_keep = max(1, int(len(candidates) * self.keep_fraction))
            kept_list = list(candidates[:n_keep])
            if self.keep_union_discrete:
                # top n_keep by the discrete-stage score as well; the candidate objects carry the
                # post-quick-MC orientation (it is the one kept, not the discrete-stage one)
                by_disc = np.argsort(disc_scores, kind="stable")[:n_keep]
                in_kept = {id(c) for c in kept_list}
                kept_list += [
                    pre_sort[int(i)] for i in by_disc if id(pre_sort[int(i)]) not in in_kept
                ]
            fz_orientations = np.array([c.orientation for c in kept_list])
            if self.recorder is not None:
                self.recorder(
                    "quick_mc",
                    dict(
                        level=level,
                        perm=np.array([ids[id(c)] for c in candidates]),  # sorted -> discrete idx
                        R=np.stack([c.orientation for c in candidates]).copy(),
                        cost=np.array([c.cost for c in candidates]),
                        n_keep=n_keep,
                    ),
                )

            # Accumulate global cost fn evals for this level
            total_global_evals += global_cost_fn.eval_count

            t_total = time.time() - t_level
            print(
                f"    Level {level}: quick MC ({t_quick:.1f}s), "
                f"kept {len(kept_list)}/{len(candidates)}, "
                f"best cost={candidates[0].cost:.4f}, "
                f"level total={t_total:.1f}s",
                flush=True,
            )

        # After all levels: FindOptimal + VarianceMinimizing + final overlap evaluation
        # C++ DiscreteAdaptive.tmpl.cpp:210-246
        if not candidates:
            return SearchCandidate(orientation=np.eye(3), cost=1.0)

        best_candidate = self.refine_from_candidates(
            candidates,
            voxel_vertices,
            phase_index,
            diameter,
            local_cost_fn=local_cost_fn,
            mc_optimizer=mc_optimizer,
            find_optimizer=find_optimizer,
        )

        # Store eval counts for benchmarking
        self._last_global_evals = total_global_evals
        self._last_local_evals = local_cost_fn.eval_count

        return best_candidate

    def _make_local_cost_fn(self) -> VoxelCostFunction:
        """Local cost function (pixel_radius=0) used for MC, re-evaluation and the final overlap."""
        return VoxelCostFunction(
            simulator=self.setup.simulator,
            detector_list=self.setup.detector_list,
            range_map=self.setup.range_map,
            exp_data=self.setup.exp_data,
            sample=self.setup.sample,
            structure_list=self.setup.structure_list,
            mode="hard",
            eta_limit=self.setup.exp_setup.get_eta_limit(),
            pixel_radius=0,
            min_sin_eta=self.min_sin_eta,
        )

    def _make_optimizers(
        self,
        local_cost_fn: VoxelCostFunction,
        voxel_vertices: torch.Tensor,
        phase_index: int,
        rng: Optional[np.random.Generator],
    ) -> Tuple[MCOptimizer, Optional[RiemannianAdamOptimizer]]:
        """The MC optimizer (always built: used by VarianceMinimizing and as the FindOptimal
        optimizer unless hybrid) and, when SearchParameters.use_hybrid_optimizer is set and
        setup.diff_cost_fn exists, the hybrid Adam optimizer for FindOptimal (else None)."""
        mc_optimizer = MCOptimizer(
            cost_fn=local_cost_fn,
            voxel_vertices=voxel_vertices,
            phase_index=phase_index,
            rng=rng,
        )
        find_optimizer: Optional[RiemannianAdamOptimizer] = None
        if self.params.use_hybrid_optimizer and self.setup.diff_cost_fn is not None:
            find_optimizer = RiemannianAdamOptimizer(
                hard_cost_fn=local_cost_fn,
                diff_cost_fn=self.setup.diff_cost_fn,
                voxel_vertices=voxel_vertices,
                phase_index=phase_index,
                rng=rng,
            )
        return mc_optimizer, find_optimizer

    def _check_local_optimizer(self) -> None:
        """local_optimizer='cma' and use_hybrid_optimizer both replace FindOptimal's MC; reject
        the combination up front (params only: before any search, whether or not a diff cost
        function exists)."""
        if self.params.local_optimizer == "cma" and self.params.use_hybrid_optimizer:
            raise ValueError(
                "local_optimizer='cma' and use_hybrid_optimizer are mutually exclusive"
            )

    def _make_cma_optimizer(
        self,
        local_cost_fn: VoxelCostFunction,
        voxel_vertices: torch.Tensor,
        phase_index: int,
        rng: Optional[np.random.Generator],
    ) -> CMAOptimizer:
        """The CMA-ES refiner for SearchParameters.local_optimizer == "cma" (per-run seeds are
        drawn from rng, so a run is deterministic given the generator)."""
        return CMAOptimizer(
            cost_fn=local_cost_fn,
            voxel_vertices=voxel_vertices,
            phase_index=phase_index,
            rng=rng,
            sigma0_deg=self.params.cma_sigma0_deg,
            max_evals=self.params.cma_max_evals,
            popsize=self.params.cma_popsize,
        )

    def refine_from_candidates(
        self,
        candidates: List[SearchCandidate],
        voxel_vertices: torch.Tensor,
        phase_index: int = 0,
        diameter: Optional[float] = None,
        rng: Optional[np.random.Generator] = None,
        local_cost_fn: Optional[VoxelCostFunction] = None,
        mc_optimizer: Optional[MCOptimizer] = None,
        find_optimizer: Optional[RiemannianAdamOptimizer] = None,
    ) -> SearchCandidate:
        """
        The final stage of reconstruct_voxel: FindOptimal, VarianceMinimizing, final overlap.

        Args:
            candidates: Candidates sorted best first (the coarse levels' hand-off in
                reconstruct_voxel; may be a single user-supplied orientation).
            voxel_vertices: Triangle vertices in sample frame, shape (3, 3)
            phase_index: Crystal phase index
            diameter: Search diameter (radians) the coarse levels ended with. The FindOptimal /
                VarianceMinimizing search box is
                max(diameter / 3, 0.2 deg) / 2^min_local_resolution, independent of how far a
                candidate is from the optimum. Default: the diameter
                reconstruct_voxel reaches after all its levels,
                local_grid_radius / 1.5^(max_local_resolution + 1).
            rng: Random generator, used only when the optimizers are built here.
            local_cost_fn, mc_optimizer, find_optimizer: the objects reconstruct_voxel built (so
                its eval counts and random stream continue); built here when omitted. Passing
                mc_optimizer without find_optimizer disables the hybrid (Adam) optimizer: full
                MC is used for FindOptimal.

        Returns:
            Best SearchCandidate, with final overlap info and cost; the identity orientation with
            cost 1.0 when candidates is empty (as reconstruct_voxel). When called directly,
            self.last_eval_counts is (0, local evals, local evals).

        C++ Reference: DiscreteAdaptive.tmpl.cpp:210-246
        """
        self._check_local_optimizer()
        if not candidates:
            return SearchCandidate(orientation=np.eye(3), cost=1.0)
        standalone = local_cost_fn is None
        if local_cost_fn is None:
            local_cost_fn = self._make_local_cost_fn()
        if mc_optimizer is None:
            mc_optimizer, find_optimizer = self._make_optimizers(
                local_cost_fn, voxel_vertices, phase_index, rng
            )
        use_hybrid = find_optimizer is not None
        use_cma = self.params.local_optimizer == "cma"
        cma_optimizer: Optional[CMAOptimizer] = None
        if use_cma:
            # seeds come from the generator the MC optimizer was built with (reconstruct_voxel's)
            cma_optimizer = self._make_cma_optimizer(
                local_cost_fn, voxel_vertices, phase_index, mc_optimizer.rng
            )
        if diameter is None:
            diameter = self.params.local_grid_radius / 1.5 ** (self.params.max_local_resolution + 1)

        final_radius = max(diameter / 3.0, math.radians(0.2))
        final_box_width = final_radius / (2**self.params.min_local_resolution)
        final_step = final_box_width * self.params.mc_radius_scale_factor

        # FindOptimal: full MC (or hybrid Adam) on top candidates with convergence check
        # C++ ContinuousSearch.h:316-346
        n_final = min(len(candidates), self.params.max_discrete_candidates)
        optimizer_label = "CMA-ES" if use_cma else ("hybrid Adam" if use_hybrid else "full MC")
        print(f"    FindOptimal: {optimizer_label} on {n_final} candidates", flush=True)
        t_final = time.time()

        best_candidate = SearchCandidate(orientation=np.eye(3), cost=1.0)
        converged = False
        best_ci = -1
        result: SearchCandidate
        for ci, cand in enumerate(candidates[:n_final]):
            if cma_optimizer is not None:
                result = cma_optimizer.optimize(cand.orientation)
            elif use_hybrid:
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
            if self.recorder is not None:
                self.recorder(
                    "find_candidate",
                    dict(
                        index=ci,
                        R_in=np.asarray(cand.orientation).copy(),
                        R_out=np.asarray(result.orientation).copy(),
                        cost=float(result.cost),
                    ),
                )
            if result.cost < best_candidate.cost:
                best_candidate = result
                best_ci = ci
                # Convergence check: hit_ratio >= 1.0
                # C++ ContinuousSearch.h:276-286 HitRatioConvergenceFn with ratio=1.0
                if result.overlap_info is not None and hit_ratio_converged(
                    result.overlap_info, 1.0
                ):
                    converged = True
                    break
        t_find = time.time() - t_final
        self.last_find_optimal = dict(
            winner_index=best_ci,
            n_evaluated=ci + 1,
            n_candidates=len(candidates),
            converged=converged,
        )
        print(
            f"    FindOptimal: {ci + 1} evaluated ({t_find:.1f}s), "
            f"best cost={best_candidate.cost:.4f}"
            f"{' CONVERGED' if converged else ''}",
            flush=True,
        )

        # VarianceMinimizing: refine best until variance < 0.02²
        # C++ DiscreteAdaptive.tmpl.cpp:232
        # (skipped with local_optimizer='cma': the CMA run per candidate replaces MC and this stage)
        if not use_cma:
            t_var = time.time()
            var_result = mc_optimizer.variance_minimizing_optimize(
                initial_orientation=best_candidate.orientation,
                search_box_side=final_box_width,
                max_mc_steps=self.params.max_mc_steps,
                successive_restarts=self.params.successive_restarts,
                max_convergence_cost=0.0,  # C++ sets this to 0 for final optimization
                convergence_variance=0.02**2,
            )
            if self.recorder is not None:
                self.recorder(
                    "variance",
                    dict(R=np.asarray(var_result.orientation).copy(), cost=float(var_result.cost)),
                )
            if var_result.cost < best_candidate.cost:
                best_candidate = var_result
            t_var_elapsed = time.time() - t_var
            print(
                f"    VarianceMin: ({t_var_elapsed:.1f}s), " f"cost={best_candidate.cost:.4f}",
                flush=True,
            )

        # Final overlap evaluation
        # C++ DiscreteAdaptive.tmpl.cpp:234-242
        final_info = local_cost_fn.evaluate(
            orientation=best_candidate.orientation,
            voxel_vertices=voxel_vertices,
            phase_index=phase_index,
        )
        best_candidate.overlap_info = final_info
        best_candidate.cost = final_info.cost
        if self.recorder is not None:
            self.recorder(
                "final",
                dict(R=np.asarray(best_candidate.orientation).copy(), cost=float(final_info.cost)),
            )

        if standalone:
            self._last_global_evals = 0
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
            mode="hard",
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
        cma_sigma0_deg: Optional[float] = None,
        cma_max_evals: Optional[int] = None,
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
            cma_sigma0_deg: local_optimizer == "cma" only: initial step override (default
                params.cma_sigma0_deg)
            cma_max_evals: local_optimizer == "cma" only: budget override (default
                params.cma_max_evals; BFS passes params.cma_neighbor_max_evals)

        Returns:
            Optimized SearchCandidate. The number of cost evaluations of the call (final
            overlap evaluation included) is left in ``last_local_optimization_evals``.

        C++ Reference: DiscreteAdaptive.tmpl.cpp:280-317 LocalOptimization
        """
        self._check_local_optimizer()
        eta_limit = self.setup.exp_setup.get_eta_limit()
        local_cost_fn = VoxelCostFunction(
            simulator=self.setup.simulator,
            detector_list=self.setup.detector_list,
            range_map=self.setup.range_map,
            exp_data=self.setup.exp_data,
            sample=self.setup.sample,
            structure_list=self.setup.structure_list,
            mode="hard",
            eta_limit=eta_limit,
            pixel_radius=0,
        )

        result: SearchCandidate
        if self.params.local_optimizer == "cma":
            # one local CMA-ES run from the inherited start replaces the variance-minimizing MC
            cma_opt = self._make_cma_optimizer(local_cost_fn, voxel_vertices, phase_index, rng)
            if cma_sigma0_deg is not None:
                if not cma_sigma0_deg > 0:
                    raise ValueError(f"cma_sigma0_deg must be > 0, got {cma_sigma0_deg}")
                cma_opt.sigma0_deg = float(cma_sigma0_deg)
            if cma_max_evals is not None:
                if cma_max_evals < 2:
                    raise ValueError(f"cma_max_evals must be >= 2, got {cma_max_evals}")
                cma_opt.max_evals = int(cma_max_evals)
            result = cma_opt.optimize(initial_orientation)
        else:
            mc_optimizer = MCOptimizer(
                cost_fn=local_cost_fn,
                voxel_vertices=voxel_vertices,
                phase_index=phase_index,
                rng=rng,
            )

            # C++ ContinuousSearch.h:70-73: BoxWidth = localGridRadius / 2^localResolution
            box_width = self.params.local_grid_radius / (2**self.params.min_local_resolution)

            # Variance-minimizing MC
            # C++ DiscreteAdaptive.tmpl.cpp:305
            result = mc_optimizer.variance_minimizing_optimize(
                initial_orientation=initial_orientation,
                search_box_side=box_width,
                max_mc_steps=self.params.max_mc_steps,
                successive_restarts=self.params.successive_restarts,
                max_convergence_cost=self.params.max_convergence_cost,
                convergence_variance=0.02**2,
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
        self.last_local_optimization_evals = int(local_cost_fn.eval_count)

        return result


# ---------------------------------------------------------------------------
# BFSReconstruction — breadth-first spatial propagation
# ---------------------------------------------------------------------------


@dataclass
class BFSVoxelRecord:
    """
    Per-voxel provenance of a BFS reconstruction (cheap: one small object per visited voxel).

    ``source`` is the final classification: "seed" (full search, accepted), "neighbor" (inherited
    start accepted at the first local fit), "neighbor_retry" (accepted only after the wider CMA
    retry), "refit" (fixed in the opt-in refit pass), "unresolved" (still REFIT at the end).
    ``n_evals`` and ``wall_s`` are cumulative over every attempt on the voxel (full search or
    local fit, retry, refit); the one overlap evaluation after a full search is not counted.
    """

    voxel_idx: int
    source: str = ""
    n_evals: int = 0
    wall_s: float = 0.0
    hit_ratio: float = 0.0
    seed_rejected: bool = False  # the full search ended below MinAccelerationThreshold
    local_rejected: bool = False  # a local fit failed the 0.9 acceptance test
    retried: bool = False  # the wider-sigma CMA retry ran
    retry_accepted: bool = False
    refit_tried: bool = False
    refit_mode: str = ""  # "local" (kept the local fit, C++ skip-discrete) or "full" (new search)


_Snapshot = Tuple[np.ndarray, float, float, float]  # orientation, cost, confidence, overlap_ratio


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
    7. Opt-in (SearchParameters.bfs_refit): one refit pass over the voxels left REFIT

    Opt-in extras (the defaults leave the C++-parity behaviour unchanged):
    - cma mode: neighbours use ``cma_neighbor_max_evals``; a rejected neighbour is retried once
      from the same inherited start with ``cma_retry_sigma0_deg`` (0 = off).
    - ``bfs_refit``: port of C++ LazyBFSClient::Refit (see ``_refit_pass``).

    After reconstruct_sample, ``records`` maps voxel index -> BFSVoxelRecord and ``stats`` holds
    the counters (seeds, retries, refit candidates/attempts/resolved, unresolved, evaluations
    and wall time split by seed / neighbour / refit).

    C++ Reference:
        BreadthFirstReconstructor.tmpl.cpp:112-188 Fit()
        BreadthFirstReconstructor.tmpl.cpp:86-103 Refit()
        ReconstructionStrategies.tmpl.cpp:334-356 InsertSeed()
    """

    def __init__(self, setup: ReconstructionSetup):
        self.setup = setup
        self.reconstructor = AdaptiveVoxelReconstructor(setup)
        self.records: Dict[int, BFSVoxelRecord] = {}
        self.stats: Dict[str, Any] = {}
        self._refit_candidates: List[int] = []
        self._phase = "neighbor"  # "neighbor" | "refit": which stats bucket local fits count in

    @staticmethod
    def _ratios(info: Optional[OverlapInfo]) -> Tuple[float, float]:
        """(hit_ratio, confidence) of an overlap evaluation.

        C++ CostFunctions.cpp: GetHitRatio = pixel_overlap / pixel_on_detector,
        GetConfidence = peak_overlap / peak_on_detector."""
        hit = (
            info.pixel_overlap / info.pixel_on_detector
            if info and info.pixel_on_detector > 0
            else 0.0
        )
        conf = (
            info.peak_overlap / info.peak_on_detector if info and info.peak_on_detector > 0 else 0.0
        )
        return hit, conf

    @staticmethod
    def _new_stats() -> Dict[str, Any]:
        """Zeroed counters (see class docstring)."""
        return dict(
            n_seeds=0,
            n_seed_rejected=0,
            n_neighbor_fits=0,
            n_neighbor_accepted_first=0,
            n_retry_attempted=0,
            n_retry_accepted=0,
            n_refit_candidates=0,
            n_refit_attempted=0,
            n_refit_local=0,
            n_refit_full=0,
            n_refit_resolved=0,
            n_unresolved=0,
            n_evals_seed=0,
            n_evals_neighbor=0,
            n_evals_refit=0,
            wall_seed_s=0.0,
            wall_neighbor_s=0.0,
            wall_refit_s=0.0,
        )

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

        self.records = {}
        self.stats = self._new_stats()
        self._refit_candidates = []
        self._phase = "neighbor"

        # Randomized seed order (C++ uses random_shuffle)
        voxel_order = list(range(n_process))
        rng.shuffle(voxel_order)

        print(f"BFS Reconstruction: {n_process} voxels", flush=True)
        start_time = time.time()
        all_processed = []
        n_seeds = 0

        for seed_idx in voxel_order:
            if mic.voxels[seed_idx].reconstruction_id != ReconstructionState.NOT_VISITED:
                continue

            n_seeds += 1
            t_seed = time.time()
            print(f"\n  Seed #{n_seeds}: voxel {seed_idx}", flush=True)

            processed = self._fit_from_seed(mic, seed_idx, rng)
            all_processed.extend(processed)

            seed_fitted = sum(
                1
                for i in processed
                if mic.voxels[i].reconstruction_id == ReconstructionState.FITTED
            )
            seed_refit = len(processed) - seed_fitted

            t_elapsed = time.time() - t_seed
            print(
                f"  Seed #{n_seeds} done: {len(processed)} voxels "
                f"({seed_fitted} fitted, {seed_refit} refit) in {t_elapsed:.1f}s",
                flush=True,
            )
        self.stats["n_seeds"] = n_seeds

        if self.setup.search_params.bfs_refit:
            self._refit_pass(mic, rng)

        for i in self._refit_candidates:
            if mic.voxels[i].reconstruction_id != ReconstructionState.FITTED:
                self.records[i].source = "unresolved"
        n_fitted = sum(
            1
            for i in all_processed
            if mic.voxels[i].reconstruction_id == ReconstructionState.FITTED
        )
        n_refit = len(all_processed) - n_fitted
        total_time = time.time() - start_time
        self.stats["n_unresolved"] = n_refit
        self.stats["n_voxels"] = len(all_processed)
        self.stats["wall_total_s"] = total_time
        print(
            f"\nBFS complete: {n_seeds} seeds, {n_fitted} fitted, "
            f"{n_refit} refit, {total_time:.1f}s total",
            flush=True,
        )

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
        rec = self.records[seed_idx] = BFSVoxelRecord(seed_idx, source="seed")

        # Full adaptive reconstruction on seed
        t0 = time.time()
        result = self.reconstructor.reconstruct_voxel(
            voxel_vertices=vertices,
            phase_index=voxel.phase,
            rng=rng,
        )
        t_recon = time.time() - t0
        n_ev = self.reconstructor.last_eval_counts[2]
        rec.n_evals += n_ev
        rec.wall_s += t_recon
        self.stats["n_evals_seed"] += n_ev
        self.stats["wall_seed_s"] += t_recon

        # Evaluate overlap
        # C++ BreadthFirstReconstructor.tmpl.cpp:136
        overlap_info = self.reconstructor.evaluate_overlap(
            result.orientation, vertices, voxel.phase
        )

        # Compute confidence and hit_ratio
        hit_ratio, confidence = self._ratios(overlap_info)

        # Update seed voxel
        voxel.orientation = result.orientation
        voxel.cost = result.cost
        voxel.confidence = confidence
        voxel.overlap_ratio = hit_ratio
        rec.hit_ratio = hit_ratio

        print(
            f"    Seed voxel {seed_idx}: cost={result.cost:.4f}, "
            f"hit_ratio={hit_ratio:.3f}, conf={confidence:.3f} ({t_recon:.1f}s)",
            flush=True,
        )

        # Check acceptance threshold
        # C++ BreadthFirstReconstructor.tmpl.cpp:142-149
        min_accel = self.setup.config.min_acceleration_threshold
        if hit_ratio < min_accel:
            voxel.reconstruction_id = ReconstructionState.REFIT
            rec.seed_rejected = True
            self.stats["n_seed_rejected"] += 1
            self._refit_candidates.append(seed_idx)
            print(
                f"    Seed rejected (hit_ratio {hit_ratio:.3f} < " f"threshold {min_accel:.3f})",
                flush=True,
            )
            return [seed_idx]

        # Mark fitted, start BFS
        # C++ BreadthFirstReconstructor.tmpl.cpp:153-156
        voxel.reconstruction_id = ReconstructionState.FITTED
        solution = [seed_idx]
        self._expand(mic, seed_idx, hit_ratio, rng, solution)
        return solution

    def _local_fit(
        self,
        voxel: Any,
        vertices: torch.Tensor,
        rng: np.random.Generator,
        rec: BFSVoxelRecord,
        sigma0_deg: Optional[float] = None,
    ) -> Tuple[SearchCandidate, float, float]:
        """One local_optimization from the voxel's current (inherited) orientation.

        In cma mode the budget is params.cma_neighbor_max_evals (seeds keep cma_max_evals);
        sigma0_deg overrides the start step (the retry). Adds evaluations and time to ``rec``
        and to the stats bucket of the current phase. Returns (result, hit_ratio, confidence)."""
        params = self.setup.search_params
        kwargs: Dict[str, Any] = {}
        if params.local_optimizer == "cma":
            kwargs["cma_max_evals"] = params.cma_neighbor_max_evals
            if sigma0_deg is not None:
                kwargs["cma_sigma0_deg"] = sigma0_deg
        t0 = time.time()
        opt_result = self.reconstructor.local_optimization(
            voxel_vertices=vertices,
            phase_index=voxel.phase,
            initial_orientation=voxel.orientation,
            rng=rng,
            **kwargs,
        )
        dt = time.time() - t0
        n_ev = self.reconstructor.last_local_optimization_evals
        rec.n_evals += n_ev
        rec.wall_s += dt
        self.stats[f"n_evals_{self._phase}"] += n_ev
        self.stats[f"wall_{self._phase}_s"] += dt
        hit, conf = self._ratios(opt_result.overlap_info)
        return opt_result, hit, conf

    def _expand(
        self,
        mic: MicFile,
        start_idx: int,
        best_conf: float,
        rng: np.random.Generator,
        solution: List[int],
        refit_pass: bool = False,
    ) -> None:
        """
        BFS expansion from a fitted voxel (the loop of C++ Fit() after the centre is accepted).

        Each popped neighbour gets a local fit from the inherited orientation and is accepted
        if hit_ratio / best_conf > 0.9. In cma mode a rejected neighbour is retried once with
        the wider cma_retry_sigma0_deg (same budget, same inherited start); the retry result is
        kept only if it passes the same test, otherwise the first fit stands. Rejected voxels are
        marked REFIT and become refit candidates. In the refit pass the expansion visits REFIT
        voxels (not NOT_VISITED ones); a still-rejected voxel keeps the better of its old and new
        fit by confidence (C++ Push keeps the higher fConfidence).
        """
        params = self.setup.search_params
        bfs_queue: deque[int] = deque()
        snapshots: Dict[int, _Snapshot] = {}
        self._insert_seed(mic, start_idx, bfs_queue, refit_pass, snapshots)
        retry_on = params.local_optimizer == "cma" and params.cma_retry_sigma0_deg > 0
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
            rec = self.records.setdefault(neighbor_idx, BFSVoxelRecord(neighbor_idx))
            if not refit_pass:
                self.stats["n_neighbor_fits"] += 1

            # local optimization from inherited orientation
            opt_result, n_hit_ratio, n_confidence = self._local_fit(neighbor, n_vertices, rng, rec)

            # Track best quality
            # C++ BreadthFirstReconstructor.tmpl.cpp:162
            best_conf = max(n_hit_ratio, best_conf)

            # Acceptance check: 90% of best quality
            # C++ BreadthFirstReconstructor.tmpl.cpp:163
            accepted = best_conf > 0 and (n_hit_ratio / best_conf) > 0.9
            first_accepted = accepted
            if not accepted:
                rec.local_rejected = True
                if retry_on:
                    # wider-start retry from the same inherited orientation (cma mode only)
                    self.stats["n_retry_attempted"] += 1
                    rec.retried = True
                    r_result, r_hit, r_conf = self._local_fit(
                        neighbor, n_vertices, rng, rec, sigma0_deg=params.cma_retry_sigma0_deg
                    )
                    r_best = max(r_hit, best_conf)
                    if r_best > 0 and (r_hit / r_best) > 0.9:
                        accepted, rec.retry_accepted = True, True
                        self.stats["n_retry_accepted"] += 1
                        opt_result, n_hit_ratio, n_confidence, best_conf = (
                            r_result,
                            r_hit,
                            r_conf,
                            r_best,
                        )

            # Update neighbor voxel
            old = snapshots.get(neighbor_idx)
            if accepted or old is None or n_confidence >= old[2]:
                neighbor.orientation = opt_result.orientation
                neighbor.cost = opt_result.cost
                neighbor.confidence = n_confidence
                neighbor.overlap_ratio = n_hit_ratio
            else:  # refit pass, rejected again, and the old fit was better: keep the old one
                neighbor.orientation, neighbor.cost, neighbor.confidence = old[0], old[1], old[2]
                neighbor.overlap_ratio = old[3]
            rec.hit_ratio = neighbor.overlap_ratio

            if accepted:
                # Accept: mark fitted, propagate to neighbors
                neighbor.reconstruction_id = ReconstructionState.FITTED
                if refit_pass:
                    rec.source = "refit"
                    self.stats["n_refit_resolved"] += 1
                else:
                    solution.append(neighbor_idx)
                    rec.source = "neighbor_retry" if rec.retry_accepted else "neighbor"
                    if first_accepted:
                        self.stats["n_neighbor_accepted_first"] += 1
                self._insert_seed(mic, neighbor_idx, bfs_queue, refit_pass, snapshots)
                print(
                    f"    BFS #{n_bfs} voxel {neighbor_idx}: FITTED "
                    f"cost={opt_result.cost:.4f}, hit_ratio={n_hit_ratio:.3f} "
                    f"({rec.wall_s:.1f}s)",
                    flush=True,
                )
            else:
                # Reject: mark for later re-fitting
                neighbor.reconstruction_id = ReconstructionState.REFIT
                rec.source = "unresolved"
                if not refit_pass:
                    solution.append(neighbor_idx)
                    self._refit_candidates.append(neighbor_idx)
                print(
                    f"    BFS #{n_bfs} voxel {neighbor_idx}: REFIT "
                    f"cost={opt_result.cost:.4f}, hit_ratio={n_hit_ratio:.3f} "
                    f"(ratio={n_hit_ratio / best_conf:.3f} < 0.9) "
                    f"({rec.wall_s:.1f}s)",
                    flush=True,
                )

    def _refit_pass(self, mic: MicFile, rng: np.random.Generator) -> None:
        """
        Opt-in port of the C++ LazyBFSClient::Refit (BreadthFirstReconstructor.tmpl.cpp:86-103).

        C++ semantics: the client receives a voxel with a stored orientation (a REFIT voxel of an
        earlier run) and (1) runs LocalOptimization from that orientation; (2) if the resulting
        confidence (peak overlap / peak on detector) >= PartialResultAcceptanceConfidence it calls
        Fit(v, skip_discrete_search=True): the locally optimized orientation is the centre (code
        PARTIAL); otherwise Fit(v, false): a full ReconstructVoxel search. (3) Fit then checks the
        centre's hit ratio against MinAccelerationThreshold (skipped only for a CONVERGED search),
        on success marks it FITTED and runs the BFS expansion (0.9 acceptance), else leaves REFIT.
        Each REFIT voxel is one work unit per run.

        Python port: after the main BFS, every voxel left REFIT (rejected seed or rejected
        neighbour) is visited once, in the order it became REFIT, unless an earlier refit's
        expansion already fixed it. Steps 1-3 are reproduced; the gate is params.bfs_refit_conf,
        else config.partial_result_acceptance_conf.

        Intentional differences from C++:
        - One pass inside the same run, not a second program run with a partial .mic as input
          (C++ RESTART_FIT resets every voxel of the partial result to NOT_VISITED; here only the
          REFIT voxels are revisited and FITTED ones are left alone).
        - The expansion from a refit centre visits REFIT neighbours only (all others are already
          visited), and a still-rejected voxel keeps the better of its old and new fit by
          confidence (the C++ Push rule) instead of the grid being reset per work unit.
        - A full-search centre is gated by the hit ratio even when the search converged (the
          Python seed path never kept the CONVERGED code).
        - The local fit uses the neighbour budget (cma_neighbor_max_evals) and no wider-sigma
          retry. Refit centres themselves are not retried.
        """
        params = self.setup.search_params
        gate = (
            params.bfs_refit_conf
            if params.bfs_refit_conf is not None
            else float(getattr(self.setup.config, "partial_result_acceptance_conf", 0.0))
        )
        min_accel = self.setup.config.min_acceleration_threshold
        candidates = [
            i
            for i in self._refit_candidates
            if mic.voxels[i].reconstruction_id != ReconstructionState.FITTED
        ]
        self.stats["n_refit_candidates"] = len(candidates)
        self._phase = "refit"
        print(f"\nRefit pass: {len(candidates)} REFIT voxels, gate {gate:.3f}", flush=True)

        for idx in candidates:
            voxel = mic.voxels[idx]
            if voxel.reconstruction_id == ReconstructionState.FITTED:
                continue  # fixed by an earlier refit's expansion
            rec = self.records[idx]
            rec.refit_tried = True
            self.stats["n_refit_attempted"] += 1
            vertices = _get_voxel_vertices(voxel)
            old: _Snapshot = (
                voxel.orientation.copy(),
                voxel.cost,
                voxel.confidence,
                voxel.overlap_ratio,
            )

            opt_result, hit, conf = self._local_fit(voxel, vertices, rng, rec)
            if conf >= gate:
                rec.refit_mode = "local"
                self.stats["n_refit_local"] += 1
                orientation, cost = opt_result.orientation, opt_result.cost
            else:
                rec.refit_mode = "full"
                self.stats["n_refit_full"] += 1
                t0 = time.time()
                full = self.reconstructor.reconstruct_voxel(
                    voxel_vertices=vertices, phase_index=voxel.phase, rng=rng
                )
                n_ev = self.reconstructor.last_eval_counts[2]
                dt = time.time() - t0
                rec.n_evals += n_ev
                rec.wall_s += dt
                self.stats["n_evals_refit"] += n_ev
                self.stats["wall_refit_s"] += dt
                info = self.reconstructor.evaluate_overlap(full.orientation, vertices, voxel.phase)
                hit, conf = self._ratios(info)
                orientation, cost = full.orientation, full.cost

            if conf >= old[2] or hit >= min_accel:
                voxel.orientation, voxel.cost = orientation, cost
                voxel.confidence, voxel.overlap_ratio = conf, hit
            else:  # not accepted and not better: keep the earlier fit
                voxel.orientation, voxel.cost, voxel.confidence, voxel.overlap_ratio = old
            rec.hit_ratio = voxel.overlap_ratio

            if hit < min_accel:
                voxel.reconstruction_id = ReconstructionState.REFIT
                continue
            voxel.reconstruction_id = ReconstructionState.FITTED
            rec.source = "refit"
            self.stats["n_refit_resolved"] += 1
            self._expand(mic, idx, hit, rng, [], refit_pass=True)

    def _insert_seed(
        self,
        mic: MicFile,
        voxel_idx: int,
        queue: deque,
        refit_pass: bool = False,
        snapshots: Optional[Dict[int, _Snapshot]] = None,
    ) -> None:
        """
        Propagate orientation to unvisited neighbors and add to BFS queue.

        With ``refit_pass`` the REFIT neighbours are enqueued instead (their previous fit is
        saved in ``snapshots`` first).

        C++ Reference: ReconstructionStrategies.tmpl.cpp:334-356 InsertSeed()
        """
        voxel = mic.voxels[voxel_idx]
        # C++ uses GetNeighbors with the solution grid
        # Python uses KDTree-based neighbor lookup with 2x side_length radius
        radius = 2.0 * voxel.side_length
        neighbors = mic.get_neighbors(voxel_idx, radius=radius)
        wanted = ReconstructionState.REFIT if refit_pass else ReconstructionState.NOT_VISITED

        for n_idx in neighbors:
            neighbor = mic.voxels[n_idx]
            if neighbor.reconstruction_id == wanted:
                if refit_pass and snapshots is not None:
                    snapshots[n_idx] = (
                        neighbor.orientation.copy(),
                        neighbor.cost,
                        neighbor.confidence,
                        neighbor.overlap_ratio,
                    )
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
            print(
                f"  Voxel {idx + 1}/{n_voxels} (phase={voxel.phase}, " f"elapsed={elapsed:.1f}s)",
                flush=True,
            )

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
            hit = oi.pixel_overlap / oi.pixel_on_detector if oi and oi.pixel_on_detector > 0 else 0
            print(
                f"  Voxel {idx + 1}/{n_voxels} done: "
                f"cost={result.cost:.4f}, hit_ratio={hit:.3f}, "
                f"time={voxel_time:.1f}s",
                flush=True,
            )

            # Update voxel orientation with reconstructed result
            voxel.orientation = result.orientation

        elapsed = time.time() - start_time
        print(f"Reconstruction complete: {n_voxels} voxels in {elapsed:.1f}s", flush=True)

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
        vertices = torch.tensor(
            [
                [x, y, z],
                [x + s, y, z],
                [x + s * 0.5, y + s * sqrt3_half, z],
            ],
            dtype=torch.float32,
        )
    else:
        vertices = torch.tensor(
            [
                [x, y, z],
                [x + s * 0.5, y - s * sqrt3_half, z],
                [x + s, y, z],
            ],
            dtype=torch.float32,
        )

    return vertices
