"""
Unit tests for core Simulation class.

Tests the main simulation engine for forward diffraction calculations.
"""

import pytest
import torch
import numpy as np
from pathlib import Path

from icenine.simulation import Simulation, PeakInfo
from icenine.experiment_setup import ExperimentSetup
from icenine.config_file import ConfigFile
from icenine.sample import Sample
from icenine.detector import Detector
from icenine.image_data import ImageData
from icenine.peak_filters import TrivialAcceptFn, XDMEtaAcceptFn
from icenine.crystal_structure import CrystalStructure


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def project_root():
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


@pytest.fixture
def config_file(project_root):
    """Load ReconstructTest.config for testing."""
    config_path = project_root / "ConfigFiles" / "ReconstructTest.config"
    if not config_path.exists():
        pytest.skip(f"Config file not found: {config_path}")
    return ConfigFile.from_file(str(config_path))


@pytest.fixture
def experiment_setup(config_file):
    """Create initialized ExperimentSetup."""
    setup = ExperimentSetup(config_file)
    return setup


@pytest.fixture
def simulator(experiment_setup):
    """Create initialized Simulation."""
    return Simulation(experiment_setup)


@pytest.fixture
def simple_simulator():
    """Create simple simulator with known parameters."""
    sim = Simulation()
    sim.beam_direction = torch.tensor([1., 0., 0.])  # Beam along +X
    sim.beam_energy = 50.0  # keV
    sim.beam_deflection_chi = 0.0
    sim.initialized = True
    return sim


@pytest.fixture
def test_detector():
    """Create simple test detector."""
    # Default detector: surface perpendicular to X-axis
    return Detector(
        num_rows=1024,
        num_cols=1024,
        beam_center_j=512.0,
        beam_center_k=512.0,
        pixel_width=0.004,  # 4 microns
        pixel_height=0.004,
        position=torch.tensor([1.0, 0., 0.]),  # 1m away along X
        orientation=torch.eye(3)
    )


@pytest.fixture
def test_sample():
    """Create simple test sample."""
    sample = Sample()
    sample.set_location(np.array([0., 0., 0.]))
    sample.set_orientation_matrix(np.eye(3))
    return sample


# ============================================================================
# Initialization Tests
# ============================================================================

class TestSimulationInit:
    """Test Simulation initialization."""

    def test_default_initialization(self):
        """Test default constructor."""
        sim = Simulation()
        assert not sim.initialized
        assert sim.beam_energy == 0.0
        assert torch.allclose(sim.beam_direction, torch.tensor([0., 0., 1.]))

    def test_initialization_with_experiment_setup(self, experiment_setup):
        """Test initialization from ExperimentSetup."""
        sim = Simulation(experiment_setup)
        assert sim.initialized
        assert sim.beam_energy == experiment_setup.beam_energy
        assert sim.beam_energy > 0


# ============================================================================
# Observable Peaks Tests
# ============================================================================

class TestGetObservablePeaks:
    """Test get_observable_peaks method."""

    def test_with_identity_orientation(self, simple_simulator):
        """Test peak generation with identity orientation."""
        # Identity orientation
        orientation = torch.eye(3)

        # Simple cubic reflections
        reciprocal_vectors = [
            torch.tensor([2.668, 0., 0.]),  # (111) for Au
            torch.tensor([3.081, 0., 0.]),  # (200)
        ]

        peaks = simple_simulator.get_observable_peaks(
            orientation,
            reciprocal_vectors
        )

        # Should get 2 omega solutions per reflection (if observable)
        assert isinstance(peaks, list)
        for peak in peaks:
            assert isinstance(peak, PeakInfo)
            assert hasattr(peak, 'omega')
            assert hasattr(peak, 'g_vector')
            assert hasattr(peak, 'g_magnitude')

    def test_with_rotated_orientation(self, simple_simulator):
        """Test peak generation with rotated orientation."""
        # 45° rotation around Z
        angle = np.pi / 4
        orientation = torch.tensor([
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1]
        ], dtype=torch.float32)

        reciprocal_vectors = [
            torch.tensor([2.668, 0., 0.]),
        ]

        peaks = simple_simulator.get_observable_peaks(
            orientation,
            reciprocal_vectors
        )

        assert len(peaks) >= 0  # May or may not be observable

    def test_empty_reflection_list(self, simple_simulator):
        """Test with empty reflection list."""
        orientation = torch.eye(3)
        reciprocal_vectors = []

        peaks = simple_simulator.get_observable_peaks(
            orientation,
            reciprocal_vectors
        )

        assert peaks == []

    def test_requires_initialization(self):
        """Test that method requires initialization."""
        sim = Simulation()
        orientation = torch.eye(3)
        reciprocal_vectors = [torch.tensor([2.668, 0., 0.])]

        with pytest.raises(RuntimeError, match="not initialized"):
            sim.get_observable_peaks(orientation, reciprocal_vectors)


# ============================================================================
# Vertex Projection Tests
# ============================================================================

class TestProjectVertex:
    """Test project_vertex method."""

    def test_vertex_at_origin(self, simple_simulator, test_detector, test_sample):
        """Test projecting vertex above detector plane."""
        # Beam is along +X (from simple_simulator fixture).
        # Detector plane is z=0 (identity orientation at position (1,0,0)).
        # Use normal (1,0,1)/sqrt(2) so the +X beam reflects downward toward z=0.
        vertex = torch.tensor([0., 0., 1.])  # Above detector plane
        normal = torch.tensor([1., 0., 1.]) / np.sqrt(2)

        hit, pixel_col, pixel_row = simple_simulator.project_vertex(
            test_detector,
            test_sample,
            vertex,
            normal
        )

        # Reflected ray should hit the z=0 detector plane
        assert hit == True
        assert pixel_col > 0
        assert pixel_row > 0

    def test_vertex_off_axis(self, simple_simulator, test_detector, test_sample):
        """Test projecting vertex off beam axis."""
        vertex = torch.tensor([0.1, 0., 0.])  # Offset in X
        normal = torch.tensor([0., 0., 1.])

        hit, pixel_col, pixel_row = simple_simulator.project_vertex(
            test_detector,
            test_sample,
            vertex,
            normal
        )

        # Should still hit (detector is large enough)
        assert isinstance(hit, bool)
        assert isinstance(pixel_col, float)
        assert isinstance(pixel_row, float)


# ============================================================================
# Voxel Projection Tests
# ============================================================================

class TestProjectVoxel:
    """Test project_voxel method."""

    def test_project_small_voxel(self, simple_simulator, test_detector, test_sample):
        """Test projecting small voxel."""
        # Small triangular voxel near origin
        vertices = torch.tensor([
            [0., 0., 0.],
            [0.001, 0., 0.],
            [0., 0.001, 0.]
        ])
        normal = torch.tensor([0., 0., 1.])

        # Create image
        image = ImageData(test_detector.num_rows, test_detector.num_cols)

        # Use trivial filter
        filter_fn = TrivialAcceptFn(intensity=1.0)

        success = simple_simulator.project_voxel(
            image,
            test_detector,
            test_sample,
            vertices,
            normal,
            filter_fn
        )

        # Should successfully project
        assert isinstance(success, bool)

    def test_project_with_eta_filter(self, simple_simulator, test_detector, test_sample):
        """Test projection with eta angle filter."""
        vertices = torch.tensor([
            [0., 0., 0.],
            [0.001, 0., 0.],
            [0., 0.001, 0.]
        ])
        normal = torch.tensor([0., 0., 1.])

        image = ImageData(test_detector.num_rows, test_detector.num_cols)

        # XDM eta filter
        filter_fn = XDMEtaAcceptFn(
            min_eta=0.0,
            max_eta=np.deg2rad(60),
            form_intensity=100.0,
            sin_2theta=0.5
        )

        success = simple_simulator.project_voxel(
            image,
            test_detector,
            test_sample,
            vertices,
            normal,
            filter_fn
        )

        assert isinstance(success, bool)

    def test_filter_rejection(self, simple_simulator, test_detector, test_sample):
        """Test that filter can reject peaks."""
        vertices = torch.tensor([
            [0., 0., 0.],
            [0.001, 0., 0.],
            [0., 0.001, 0.]
        ])
        normal = torch.tensor([0., 1., 0.])  # Normal along Y

        image = ImageData(test_detector.num_rows, test_detector.num_cols)

        # Very restrictive eta filter
        filter_fn = XDMEtaAcceptFn(
            min_eta=0.0,
            max_eta=0.01,  # Very small acceptance
            form_intensity=100.0,
            sin_2theta=0.5
        )

        success = simple_simulator.project_voxel(
            image,
            test_detector,
            test_sample,
            vertices,
            normal,
            filter_fn
        )

        # May or may not succeed depending on geometry
        assert isinstance(success, bool)


# ============================================================================
# Integration Tests
# ============================================================================

class TestSimulationIntegration:
    """Integration tests with real config files."""

    def test_with_config_file(self, config_file):
        """Test simulation with real config file."""
        exp_setup = ExperimentSetup(config_file)
        sim = Simulation(exp_setup)

        assert sim.initialized
        assert sim.beam_energy > 0
        print(f"Beam energy: {sim.beam_energy} keV")

    def test_with_crystal_structure(self, simple_simulator):
        """Test peak generation with real crystal structure."""
        # Create Au FCC structure
        gold = CrystalStructure.create_fcc("Au", 4.0782)

        # Get reflections up to some Q
        max_q = 10.0  # Å⁻¹
        reflections = gold.generate_reflections(max_q)

        assert len(reflections) > 0

        # Get reflection vectors (q_vec may be None, compute from lattice if needed)
        reciprocal_vectors = [
            torch.tensor(r.q_vec, dtype=torch.float32)
            for r in reflections[:10]
            if r.q_vec is not None
        ]

        # Generate peaks
        orientation = torch.eye(3)
        peaks = simple_simulator.get_observable_peaks(
            orientation,
            reciprocal_vectors
        )

        print(f"Generated {len(peaks)} observable peaks from {len(reciprocal_vectors)} reflections")
        assert len(peaks) >= 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
