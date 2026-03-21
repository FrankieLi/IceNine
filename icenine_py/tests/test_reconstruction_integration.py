"""
Integration tests for the full reconstruction pipeline.

Tests the end-to-end workflow:
  Forward Simulation → Experimental Data → Reconstruction → Orientation Verification

Uses Example2.ThreeVoxels as the test case. The ground truth orientations are
known from the input .mic file, and the forward simulation output serves as
synthetic experimental data.
"""

import math
from pathlib import Path

import numpy as np
import pytest
import torch

from icenine.config_file import ConfigFile
from icenine.cost_functions import OverlapInfo, VoxelCostFunction
from icenine.experiment_setup import XDMExperimentSetup
from icenine.experimental_data import ExperimentalData
from icenine.mic_file import MicFile
from icenine.orientation_search import MCOptimizer, SearchCandidate, run_discrete_search
from icenine.reconstructor import (
    BasicVoxelReconstructor,
    ReconstructionSetup,
    _get_voxel_vertices,
    setup_reconstruction,
)
from icenine.sample import Sample
from icenine.sampling import (
    generate_local_grid,
    load_fundamental_zone_file,
    matrix_to_quaternion,
    get_misorientation,
)
from icenine.simulation import Simulation
from icenine.symmetry import create_cubic_symmetry


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def project_root():
    return Path(__file__).parent.parent.parent


@pytest.fixture(autouse=True)
def chdir_to_project_root(project_root, monkeypatch):
    """Config files use relative paths from project root."""
    monkeypatch.chdir(project_root)


@pytest.fixture
def example_dir(project_root):
    return project_root / "Examples" / "Example2.ThreeVoxels"


@pytest.fixture
def sim_config(example_dir):
    """Load ThreeVoxels simulation config."""
    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
    if not config_path.exists():
        pytest.skip(f"Config not found: {config_path}")
    # chdir to example dir for relative path resolution
    import os
    os.chdir(example_dir)
    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = "3Grains.sim"
    return config


@pytest.fixture
def ground_truth_mic(example_dir):
    """Load ground truth orientations from input .mic file."""
    mic_path = example_dir / "SimInput" / "three_voxels.mic"
    if not mic_path.exists():
        pytest.skip(f"Ground truth .mic not found: {mic_path}")
    return MicFile.read(str(mic_path))


@pytest.fixture
def exp_data(example_dir):
    """Load forward simulation output as experimental data."""
    data_dir = example_dir / "ScatteringData_Python"
    if not data_dir.exists() or len(list(data_dir.glob("*.d*"))) < 360:
        pytest.skip(f"Forward sim output not found or incomplete: {data_dir}")
    return ExperimentalData.from_image_directory(
        directory=str(data_dir),
        basename="3Grains.sim",
        ext="d",
        serial_length=5,
        n_omega=180,
        n_detectors=2,
        num_rows=2048,
        num_cols=2048,
    )


@pytest.fixture
def cubic_symmetry_quats():
    """Get 24 proper cubic symmetry quaternions for misorientation checks."""
    sym = create_cubic_symmetry(4.0)
    matrices = sym.get_rotation_matrices()
    proper = [m for m in matrices if np.linalg.det(m) > 0]
    quats = np.array([matrix_to_quaternion(np.array(m)) for m in proper])
    return quats


# ============================================================================
# Test: Experimental data loading sanity
# ============================================================================

class TestExperimentalDataSanity:
    """Verify loaded experimental data has expected properties."""

    def test_data_dimensions(self, exp_data):
        """180 omega intervals, 2 detectors."""
        assert exp_data.n_omega_intervals == 180
        assert exp_data.n_detectors == 2

    def test_has_bright_pixels(self, exp_data):
        """Should have significant number of bright pixels from 3 voxels."""
        n_bright = exp_data.count_bright_pixels()
        assert n_bright > 1000, f"Expected >1000 bright pixels, got {n_bright}"


# ============================================================================
# Test: Cost function evaluation with known orientation
# ============================================================================

class TestCostFunctionWithGroundTruth:
    """Verify cost function gives good scores for correct orientations."""

    def test_correct_orientation_has_overlap(
        self, sim_config, ground_truth_mic, exp_data
    ):
        """
        A voxel evaluated at its ground truth orientation should produce
        non-zero pixel overlap with the experimental data.
        """
        # Initialize experiment
        exp_setup = XDMExperimentSetup(sim_config)
        exp_setup.initialize_experiment()

        detector_list = exp_setup.get_detector_list()
        range_map = exp_setup.get_range_to_index_map()

        # Initialize sample and simulator
        sample = Sample()
        exp_setup.initialize_sample(sample, detector_list[0])
        simulator = Simulation(exp_setup)
        structure_list = sample.get_structure_list()

        # Create cost function
        cost_fn = VoxelCostFunction(
            simulator=simulator,
            detector_list=detector_list,
            range_map=range_map,
            exp_data=exp_data,
            sample=sample,
            structure_list=structure_list,
            mode='hard',
        )

        # Evaluate at ground truth orientation for voxel 0
        voxel = ground_truth_mic.voxels[0]
        vertices = _get_voxel_vertices(voxel)

        overlap_info = cost_fn.evaluate(
            orientation=voxel.orientation,
            voxel_vertices=vertices,
            phase_index=voxel.phase,  # phase indexes directly into structure_list
        )

        # Ground truth orientation should produce some overlap
        assert overlap_info.pixel_on_detector > 0, (
            "Ground truth orientation should produce pixels on detector"
        )
        assert overlap_info.peak_on_detector > 0, (
            "Ground truth orientation should produce qualified peaks"
        )


# ============================================================================
# Test: Single-voxel reconstruction
# ============================================================================

class TestSingleVoxelReconstruction:
    """Test reconstructing a single voxel from synthetic data."""

    @pytest.fixture
    def recon_setup(self, sim_config, exp_data):
        """Set up reconstruction components (shared across tests)."""
        exp_setup = XDMExperimentSetup(sim_config)
        exp_setup.initialize_experiment()
        detector_list = exp_setup.get_detector_list()
        range_map = exp_setup.get_range_to_index_map()
        sample = Sample()
        exp_setup.initialize_sample(sample, detector_list[0])
        simulator = Simulation(exp_setup)
        structure_list = sample.get_structure_list()

        cost_fn = VoxelCostFunction(
            simulator=simulator,
            detector_list=detector_list,
            range_map=range_map,
            exp_data=exp_data,
            sample=sample,
            structure_list=structure_list,
            mode='hard',
        )
        return cost_fn

    def test_ground_truth_has_overlap(
        self, recon_setup, ground_truth_mic, cubic_symmetry_quats
    ):
        """
        The ground truth orientation should produce non-zero overlap and
        better cost than a random orientation.

        Note: For sparse data (3 voxels), absolute quality is low because
        most Bragg peaks land in empty detector regions. What matters is
        that ground truth has BETTER quality than wrong orientations.
        """
        cost_fn = recon_setup
        voxel = ground_truth_mic.voxels[0]
        vertices = _get_voxel_vertices(voxel)

        gt_info = cost_fn.evaluate(
            orientation=voxel.orientation,
            voxel_vertices=vertices,
            phase_index=voxel.phase,
        )

        # Ground truth must produce pixel overlap
        assert gt_info.pixel_overlap > 0, (
            "Ground truth should have non-zero pixel overlap"
        )
        assert gt_info.peak_overlap > 0, (
            "Ground truth should have non-zero peak overlap"
        )
        # Pixel hit ratio should be high (near 1.0)
        assert gt_info.hit_ratio > 0.5, (
            f"Ground truth should have high hit ratio, got {gt_info.hit_ratio:.3f}"
        )

        # Compare against a random orientation — should be worse
        rng = np.random.default_rng(42)
        from scipy.spatial.transform import Rotation
        random_orient = Rotation.random(random_state=42).as_matrix().astype(np.float32)
        rand_info = cost_fn.evaluate(
            orientation=random_orient,
            voxel_vertices=vertices,
            phase_index=voxel.phase,
        )

        assert gt_info.quality > rand_info.quality, (
            f"Ground truth quality ({gt_info.quality:.4f}) should exceed "
            f"random ({rand_info.quality:.4f})"
        )

    def test_mc_converges_from_perturbation(
        self, recon_setup, ground_truth_mic, cubic_symmetry_quats
    ):
        """
        Start from ground truth + small random perturbation (~1.5 deg),
        run MC optimization, verify it converges back close to ground truth.
        """
        cost_fn = recon_setup
        voxel = ground_truth_mic.voxels[0]
        vertices = _get_voxel_vertices(voxel)
        gt_orientation = voxel.orientation

        # Apply small perturbation (~1.5 degrees)
        from scipy.spatial.transform import Rotation
        perturbation = Rotation.from_rotvec(
            np.array([0.02, 0.01, -0.01])
        ).as_matrix()
        start_orientation = perturbation @ gt_orientation

        mc = MCOptimizer(
            cost_fn=cost_fn,
            voxel_vertices=vertices,
            phase_index=voxel.phase,
            rng=np.random.default_rng(42),
        )

        result = mc.optimize(
            initial_orientation=start_orientation,
            angular_box_side=math.radians(3.0),
            angular_step=math.radians(0.5),
            max_mc_steps=100,
            max_restarts=1,
            max_convergence_cost=0.001,
        )

        # Should find overlap
        assert result.cost < 1.0, (
            f"MC should find overlap, got cost={result.cost}"
        )

        # Should converge close to ground truth
        q_result = matrix_to_quaternion(result.orientation)
        q_truth = matrix_to_quaternion(gt_orientation)
        misori = get_misorientation(q_result, q_truth, cubic_symmetry_quats)
        misori_deg = math.degrees(misori)

        assert misori_deg < 10.0, (
            f"Expected misorientation < 10 deg, got {misori_deg:.2f} deg"
        )
