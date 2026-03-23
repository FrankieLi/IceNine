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
from .mic_file import MicFile
from .orientation_search import (
    MCOptimizer,
    SearchCandidate,
    SearchParameters,
    hit_ratio_converged,
    run_discrete_search,
)
from .sample import Sample
from .sampling import (
    generate_local_grid,
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
