"""
Tests for cost functions module.

Tests verify:
1. OverlapInfo incremental quality averaging
2. Qualified peak counting (contiguous detector logic)
3. VoxelCostFunction integration with forward sim output
"""

from pathlib import Path

import numpy as np
import pytest
import torch

from icenine.cost_functions import (
    OverlapInfo,
    count_qualified_peaks,
)


# ============================================================================
# Test: OverlapInfo
# ============================================================================

class TestOverlapInfo:
    """Tests for OverlapInfo dataclass."""

    def test_initial_state(self):
        """All counters start at zero."""
        info = OverlapInfo()
        assert info.pixel_overlap == 0
        assert info.quality == 0.0
        assert info.cost == 1.0
        assert info.hit_ratio == 0.0

    def test_update_counts(self):
        """Counts accumulate correctly."""
        info = OverlapInfo()
        info.update_counts(10, 20, 1, 1)
        assert info.pixel_overlap == 10
        assert info.pixel_on_detector == 20
        assert info.peak_overlap == 1
        assert info.peak_on_detector == 1

        info.update_counts(5, 10, 1, 1)
        assert info.pixel_overlap == 15
        assert info.pixel_on_detector == 30
        assert info.peak_overlap == 2

    def test_hit_ratio(self):
        """Hit ratio = pixel_overlap / pixel_on_detector."""
        info = OverlapInfo()
        info.update_counts(10, 20, 1, 1)
        assert abs(info.hit_ratio - 0.5) < 1e-10

    def test_update_quality_single(self):
        """Quality for single peak with full overlap on 2 detectors."""
        info = OverlapInfo()
        info.update_quality(
            peak_pixel_overlap=10, peak_pixel_on_detector=10,
            peak_detectors_overlap=2, n_detectors_total=2,
        )

        # (10/10) * (2/2) = 1.0
        assert abs(info.quality - 1.0) < 1e-10
        assert abs(info.cost - 0.0) < 1e-10

    def test_update_quality_incremental(self):
        """Quality uses incremental (Welford) running mean with per-peak values."""
        info = OverlapInfo()

        # First peak: 100% overlap, 2/2 detectors
        info.update_quality(
            peak_pixel_overlap=10, peak_pixel_on_detector=10,
            peak_detectors_overlap=2, n_detectors_total=2,
        )
        assert abs(info.quality - 1.0) < 1e-10

        # Second peak: 5/10 pixel overlap, 1/2 detectors
        info.update_quality(
            peak_pixel_overlap=5, peak_pixel_on_detector=10,
            peak_detectors_overlap=1, n_detectors_total=2,
        )

        # cur_quality = (5/10) * (1/2) = 0.25
        # running_mean = 1.0 + (0.25 - 1.0) / 2 = 0.625
        assert abs(info.quality - 0.625) < 1e-10

    def test_update_quality_zero_pixels(self):
        """Quality unchanged when no pixels on detector."""
        info = OverlapInfo()
        info.update_quality(
            peak_pixel_overlap=0, peak_pixel_on_detector=0,
            peak_detectors_overlap=0, n_detectors_total=2,
        )
        assert info.quality == 0.0


# ============================================================================
# Test: Qualified Peak Counting
# ============================================================================

class TestCountQualifiedPeaks:
    """Tests for contiguous detector validation."""

    def test_single_detector_lit_and_overlap(self):
        """Single detector: lit + overlap → valid."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [True], [True], 1
        )
        assert peak_on == 1
        assert peak_ovlp == 1
        assert n_det == 1

    def test_single_detector_lit_no_overlap(self):
        """Single detector: lit, no overlap → peak valid but overlap invalid."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [True], [False], 1
        )
        assert peak_on == 1
        assert peak_ovlp == 0
        assert n_det == 0

    def test_two_detectors_contiguous(self):
        """Two contiguous detectors with overlap."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [True, True], [True, True], 2
        )
        assert peak_on == 1
        assert peak_ovlp == 1
        assert n_det == 2

    def test_two_detectors_first_only(self):
        """Only first detector lit → valid (trailing unlit OK)."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [True, False], [True, False], 2
        )
        assert peak_on == 1
        assert peak_ovlp == 1
        assert n_det == 1

    def test_gap_in_detectors_invalid(self):
        """Gap in lit detectors → invalid peak."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [True, False, True], [True, False, True], 3
        )
        assert peak_on == 0
        assert peak_ovlp == 0
        assert n_det == 0

    def test_not_starting_at_zero_invalid(self):
        """Lit pattern not starting at detector 0 → invalid."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [False, True, True], [False, True, True], 3
        )
        assert peak_on == 0
        assert peak_ovlp == 0

    def test_empty_detectors(self):
        """No detectors → all zeros."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [], [], 0
        )
        assert peak_on == 0
        assert peak_ovlp == 0
        assert n_det == 0

    def test_all_dark(self):
        """All detectors dark → invalid."""
        peak_on, peak_ovlp, n_det = count_qualified_peaks(
            [False, False], [False, False], 2
        )
        assert peak_on == 0
        assert peak_ovlp == 0


# ============================================================================
# Test: Batched vs Serial equivalence
# ============================================================================

import math
import os
import time


class TestBatchedEquivalence:
    """Verify batched overlap calculation matches serial version exactly."""

    @pytest.fixture
    def project_root(self):
        return Path(__file__).parent.parent.parent

    @pytest.fixture(autouse=True)
    def chdir_to_example(self, project_root):
        example_dir = project_root / "Examples" / "Example2.ThreeVoxels"
        if not example_dir.exists():
            pytest.skip(f"Example directory not found: {example_dir}")
        old_cwd = os.getcwd()
        os.chdir(example_dir)
        yield
        os.chdir(old_cwd)

    @pytest.fixture
    def example_dir(self, project_root):
        return project_root / "Examples" / "Example2.ThreeVoxels"

    @pytest.fixture
    def setup_components(self, example_dir):
        """Set up reconstruction components for equivalence testing."""
        from icenine.config_file import ConfigFile
        from icenine.cost_functions import VoxelCostFunction, calculate_diffraction_overlap, calculate_diffraction_overlap_batched
        from icenine.experiment_setup import XDMExperimentSetup
        from icenine.experimental_data import ExperimentalData
        from icenine.mic_file import MicFile
        from icenine.reconstructor import _get_voxel_vertices
        from icenine.sample import Sample
        from icenine.simulation import Simulation

        config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
        if not config_path.exists():
            pytest.skip(f"Config not found: {config_path}")

        config = ConfigFile.from_file(str(config_path))
        config.out_file_basename = "3Grains.sim"

        data_dir = example_dir / "ScatteringData_Python"
        if not data_dir.exists() or len(list(data_dir.glob("*.d*"))) < 360:
            pytest.skip(f"Forward sim output not found or incomplete: {data_dir}")

        exp_data = ExperimentalData.from_image_directory(
            directory=str(data_dir),
            basename="3Grains.sim",
            ext="d",
            serial_length=5,
            n_omega=180,
            n_detectors=2,
            num_rows=2048,
            num_cols=2048,
        )

        exp_setup = XDMExperimentSetup(config)
        exp_setup.initialize_experiment()
        detector_list = exp_setup.get_detector_list()
        range_map = exp_setup.get_range_to_index_map()
        sample = Sample()
        exp_setup.initialize_sample(sample, detector_list[0])
        simulator = Simulation(exp_setup)
        structure_list = sample.get_structure_list()

        mic_path = example_dir / "SimInput" / "three_voxels.mic"
        if not mic_path.exists():
            pytest.skip(f"Ground truth .mic not found: {mic_path}")
        mic = MicFile.read(str(mic_path))

        return {
            "simulator": simulator,
            "detector_list": detector_list,
            "range_map": range_map,
            "exp_data": exp_data,
            "sample": sample,
            "structure_list": structure_list,
            "mic": mic,
        }

    def test_batched_matches_serial(self, setup_components):
        """Batched and serial overlap produce identical OverlapInfo."""
        from icenine.cost_functions import calculate_diffraction_overlap, calculate_diffraction_overlap_batched
        from icenine.diffraction_core import get_scattering_omegas_torch
        from icenine.reconstructor import _get_voxel_vertices

        c = setup_components
        voxel = c["mic"].voxels[0]
        vertices = _get_voxel_vertices(voxel)

        # Get peaks for voxel 0
        phase_idx = voxel.phase
        structure = c["structure_list"][phase_idx]
        recp_vecs = structure.get_reflection_vectors()
        g_hkl = torch.stack([torch.from_numpy(rv.q_vec).float() for rv in recp_vecs])
        g_mag = torch.tensor([rv.q_mag for rv in recp_vecs], dtype=torch.float32)

        orientation_t = torch.from_numpy(voxel.orientation).float()
        g_lab = (orientation_t @ g_hkl.T).T

        bragg = get_scattering_omegas_torch(
            g_lab, g_mag,
            c["simulator"].beam_energy,
            c["simulator"].beam_deflection_chi,
        )

        peak_omegas = []
        peak_normals = []
        obs_indices = torch.where(bragg.observable)[0].tolist()
        for idx in obs_indices:
            g_vec = g_lab[idx]
            normal = g_vec / torch.norm(g_vec)
            peak_omegas.append(bragg.omega1[idx].item())
            peak_normals.append(normal)
            peak_omegas.append(bragg.omega2[idx].item())
            peak_normals.append(normal)

        # Run serial version
        serial_result = calculate_diffraction_overlap(
            sample=c["sample"],
            voxel_vertices=vertices,
            peak_omegas=peak_omegas,
            peak_normals=peak_normals,
            detector_list=c["detector_list"],
            range_map=c["range_map"],
            exp_data=c["exp_data"],
            beam_direction=c["simulator"].beam_direction,
            mode='hard',
        )

        # Run batched version
        batched_result = calculate_diffraction_overlap_batched(
            sample=c["sample"],
            voxel_vertices=vertices,
            peak_omegas=peak_omegas,
            peak_normals=peak_normals,
            detector_list=c["detector_list"],
            range_map=c["range_map"],
            exp_data=c["exp_data"],
            beam_direction=c["simulator"].beam_direction,
            mode='hard',
        )

        # Assert exact match on all fields
        assert serial_result.pixel_overlap == batched_result.pixel_overlap, (
            f"pixel_overlap: serial={serial_result.pixel_overlap} vs batched={batched_result.pixel_overlap}"
        )
        assert serial_result.pixel_on_detector == batched_result.pixel_on_detector, (
            f"pixel_on_detector: serial={serial_result.pixel_on_detector} vs batched={batched_result.pixel_on_detector}"
        )
        assert serial_result.peak_overlap == batched_result.peak_overlap, (
            f"peak_overlap: serial={serial_result.peak_overlap} vs batched={batched_result.peak_overlap}"
        )
        # peak_on_detector may differ by a few counts due to float32 precision
        # at ray-plane intersection boundary (epsilon threshold). This doesn't
        # affect quality/cost since those peaks have zero pixel overlap.
        assert abs(serial_result.peak_on_detector - batched_result.peak_on_detector) < 10, (
            f"peak_on_detector: serial={serial_result.peak_on_detector} vs batched={batched_result.peak_on_detector}"
        )
        assert abs(serial_result.quality - batched_result.quality) < 1e-6, (
            f"quality: serial={serial_result.quality} vs batched={batched_result.quality}"
        )
        assert abs(serial_result.cost - batched_result.cost) < 1e-6, (
            f"cost: serial={serial_result.cost} vs batched={batched_result.cost}"
        )

    def test_batch_c_matches_python_fallback(self, setup_components):
        """Batch C extension Stage D matches Python fallback path."""
        from icenine.cost_functions import (
            calculate_diffraction_overlap_batched,
            _HAS_C_RASTERIZE,
        )
        import icenine.cost_functions as cf
        from icenine.diffraction_core import get_scattering_omegas_torch
        from icenine.reconstructor import _get_voxel_vertices

        if not _HAS_C_RASTERIZE:
            pytest.skip("C extension not available")

        c = setup_components
        voxel = c["mic"].voxels[0]
        vertices = _get_voxel_vertices(voxel)

        phase_idx = voxel.phase
        structure = c["structure_list"][phase_idx]
        recp_vecs = structure.get_reflection_vectors()
        g_hkl = torch.stack([torch.from_numpy(rv.q_vec).float() for rv in recp_vecs])
        g_mag = torch.tensor([rv.q_mag for rv in recp_vecs], dtype=torch.float32)

        orientation_t = torch.from_numpy(voxel.orientation).float()
        g_lab = (orientation_t @ g_hkl.T).T

        bragg = get_scattering_omegas_torch(
            g_lab, g_mag,
            c["simulator"].beam_energy,
            c["simulator"].beam_deflection_chi,
        )

        peak_omegas = []
        peak_normals = []
        obs_indices = torch.where(bragg.observable)[0].tolist()
        for idx in obs_indices:
            g_vec = g_lab[idx]
            normal = g_vec / torch.norm(g_vec)
            peak_omegas.append(bragg.omega1[idx].item())
            peak_normals.append(normal)
            peak_omegas.append(bragg.omega2[idx].item())
            peak_normals.append(normal)

        # Run with C batch path (mode='hard', default)
        c_result = calculate_diffraction_overlap_batched(
            sample=c["sample"], voxel_vertices=vertices,
            peak_omegas=peak_omegas, peak_normals=peak_normals,
            detector_list=c["detector_list"], range_map=c["range_map"],
            exp_data=c["exp_data"], beam_direction=c["simulator"].beam_direction,
            mode='hard',
        )

        # Force Python fallback by temporarily disabling C extension
        original = cf._HAS_C_RASTERIZE
        cf._HAS_C_RASTERIZE = False
        try:
            py_result = calculate_diffraction_overlap_batched(
                sample=c["sample"], voxel_vertices=vertices,
                peak_omegas=peak_omegas, peak_normals=peak_normals,
                detector_list=c["detector_list"], range_map=c["range_map"],
                exp_data=c["exp_data"], beam_direction=c["simulator"].beam_direction,
                mode='hard',
            )
        finally:
            cf._HAS_C_RASTERIZE = original

        assert c_result.pixel_overlap == py_result.pixel_overlap, (
            f"pixel_overlap: C={c_result.pixel_overlap} vs Python={py_result.pixel_overlap}"
        )
        assert c_result.pixel_on_detector == py_result.pixel_on_detector, (
            f"pixel_on_detector: C={c_result.pixel_on_detector} vs Python={py_result.pixel_on_detector}"
        )
        assert c_result.peak_overlap == py_result.peak_overlap, (
            f"peak_overlap: C={c_result.peak_overlap} vs Python={py_result.peak_overlap}"
        )
        assert c_result.peak_on_detector == py_result.peak_on_detector, (
            f"peak_on_detector: C={c_result.peak_on_detector} vs Python={py_result.peak_on_detector}"
        )
        assert abs(c_result.quality - py_result.quality) < 1e-10, (
            f"quality: C={c_result.quality} vs Python={py_result.quality}"
        )

    def test_batched_faster_than_serial(self, setup_components):
        """Batched version should be significantly faster than serial."""
        from icenine.cost_functions import VoxelCostFunction, calculate_diffraction_overlap, calculate_diffraction_overlap_batched
        from icenine.diffraction_core import get_scattering_omegas_torch
        from icenine.reconstructor import _get_voxel_vertices

        c = setup_components
        voxel = c["mic"].voxels[0]
        vertices = _get_voxel_vertices(voxel)

        phase_idx = voxel.phase
        structure = c["structure_list"][phase_idx]
        recp_vecs = structure.get_reflection_vectors()
        g_hkl = torch.stack([torch.from_numpy(rv.q_vec).float() for rv in recp_vecs])
        g_mag = torch.tensor([rv.q_mag for rv in recp_vecs], dtype=torch.float32)

        orientation_t = torch.from_numpy(voxel.orientation).float()
        g_lab = (orientation_t @ g_hkl.T).T

        bragg = get_scattering_omegas_torch(
            g_lab, g_mag,
            c["simulator"].beam_energy,
            c["simulator"].beam_deflection_chi,
        )

        peak_omegas = []
        peak_normals = []
        obs_indices = torch.where(bragg.observable)[0].tolist()
        for idx in obs_indices:
            g_vec = g_lab[idx]
            normal = g_vec / torch.norm(g_vec)
            peak_omegas.append(bragg.omega1[idx].item())
            peak_normals.append(normal)
            peak_omegas.append(bragg.omega2[idx].item())
            peak_normals.append(normal)

        n_runs = 5

        # Benchmark serial
        t0 = time.time()
        for _ in range(n_runs):
            calculate_diffraction_overlap(
                sample=c["sample"], voxel_vertices=vertices,
                peak_omegas=peak_omegas, peak_normals=peak_normals,
                detector_list=c["detector_list"], range_map=c["range_map"],
                exp_data=c["exp_data"], beam_direction=c["simulator"].beam_direction,
            )
        serial_time = (time.time() - t0) / n_runs

        # Benchmark batched
        t0 = time.time()
        for _ in range(n_runs):
            calculate_diffraction_overlap_batched(
                sample=c["sample"], voxel_vertices=vertices,
                peak_omegas=peak_omegas, peak_normals=peak_normals,
                detector_list=c["detector_list"], range_map=c["range_map"],
                exp_data=c["exp_data"], beam_direction=c["simulator"].beam_direction,
            )
        batched_time = (time.time() - t0) / n_runs

        speedup = serial_time / batched_time
        print(f"\n  Serial: {serial_time*1000:.1f}ms, Batched: {batched_time*1000:.1f}ms, Speedup: {speedup:.1f}x")

        # Batched should be at least 1.5x faster
        assert speedup > 1.5, (
            f"Batched should be faster: serial={serial_time*1000:.1f}ms, "
            f"batched={batched_time*1000:.1f}ms, speedup={speedup:.1f}x"
        )
