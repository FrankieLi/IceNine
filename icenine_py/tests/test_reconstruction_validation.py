"""
Comprehensive reconstruction validation: blind search on ThreeVoxels.

Runs the full reconstruction pipeline (discrete search → MC optimization)
on all 3 ThreeVoxels voxels starting from blind search (no hint of ground
truth), then verifies recovered orientations match ground truth within
numerical precision of the search grid.

Generates a persistent validation report at:
    icenine_py/tests/reconstruction_validation_report.txt
"""

import json
import math
import os
import time
from datetime import datetime
from pathlib import Path
from textwrap import dedent
from typing import List

import numpy as np
import pytest
import torch

from icenine.config_file import ConfigFile
from icenine.experimental_data import ExperimentalData
from icenine.mic_file import MicFile
from icenine.reconstructor import (
    BasicVoxelReconstructor,
    ReconstructionSetup,
    _get_voxel_vertices,
    setup_reconstruction,
)
from icenine.sampling import get_misorientation, matrix_to_quaternion
from icenine.symmetry import create_cubic_symmetry


# ============================================================================
# Helpers
# ============================================================================

def _format_matrix(m: np.ndarray) -> str:
    """Format a 3x3 matrix as a compact multi-line string."""
    rows = []
    for i in range(3):
        vals = " ".join(f"{m[i, j]:8.5f}" for j in range(3))
        rows.append(f"  [{vals}]")
    return "\n".join(rows)


def _rotation_axis_angle(m: np.ndarray):
    """Extract single-axis rotation representation (axis, angle) from rotation matrix."""
    from scipy.spatial.transform import Rotation
    r = Rotation.from_matrix(m)
    rotvec = r.as_rotvec()
    angle = np.linalg.norm(rotvec)
    if angle < 1e-12:
        axis = np.array([0.0, 0.0, 1.0])
    else:
        axis = rotvec / angle
    return axis, math.degrees(angle)


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture(scope="module")
def project_root():
    return Path(__file__).parent.parent.parent


@pytest.fixture(scope="module", autouse=True)
def chdir_to_example(project_root):
    """Config files use relative paths from example directory."""
    example_dir = project_root / "Examples" / "Example2.ThreeVoxels"
    if not example_dir.exists():
        pytest.skip(f"Example directory not found: {example_dir}")
    old_cwd = os.getcwd()
    os.chdir(example_dir)
    yield
    os.chdir(old_cwd)


@pytest.fixture(scope="module")
def example_dir(project_root):
    return project_root / "Examples" / "Example2.ThreeVoxels"


@pytest.fixture(scope="module")
def ground_truth_mic(example_dir):
    """Load ground truth orientations from input .mic file."""
    mic_path = example_dir / "SimInput" / "three_voxels.mic"
    if not mic_path.exists():
        pytest.skip(f"Ground truth .mic not found: {mic_path}")
    return MicFile.read(str(mic_path))


@pytest.fixture(scope="module")
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


@pytest.fixture(scope="module")
def sim_config(example_dir):
    """Load ThreeVoxels simulation config with moderate reconstruction params."""
    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"
    if not config_path.exists():
        pytest.skip(f"Config not found: {config_path}")
    config = ConfigFile.from_file(str(config_path))
    config.out_file_basename = "3Grains.sim"

    # Use MaxQ=8 for faster evaluation (~112 reciprocal vectors vs 868 at MaxQ=16)
    config.max_q = 8.0

    # Match C++ behavior: aggressive candidate pruning, moderate MC
    config.max_local_resolution = 3
    config.max_mc_steps = 200
    config.successive_restarts = 2
    config.max_discrete_candidates = 30
    return config


@pytest.fixture(scope="module")
def recon_setup(sim_config, exp_data):
    """Set up the full reconstruction pipeline (shared across all voxel tests)."""
    return setup_reconstruction(sim_config, exp_data=exp_data)


@pytest.fixture(scope="module")
def cubic_symmetry_quats():
    """24 proper cubic symmetry quaternions for misorientation calculation."""
    sym = create_cubic_symmetry(4.0)
    matrices = sym.get_rotation_matrices()
    proper = [m for m in matrices if np.linalg.det(m) > 0]
    return np.array([matrix_to_quaternion(np.array(m)) for m in proper])


@pytest.fixture(scope="module")
def report_path(project_root):
    return project_root / "icenine_py" / "tests" / "reconstruction_validation_report.txt"


# ============================================================================
# Validation test
# ============================================================================

from dataclasses import dataclass


@dataclass
class VoxelValidationResult:
    voxel_index: int
    ground_truth: np.ndarray
    reconstructed: np.ndarray
    misorientation_deg: float
    cost: float
    hit_ratio: float
    quality: float
    convergence_code: str
    elapsed_seconds: float
    gt_axis: np.ndarray
    gt_angle_deg: float
    recon_axis: np.ndarray
    recon_angle_deg: float


@pytest.mark.slow
class TestBlindReconstruction:
    """
    Run full blind reconstruction on all 3 ThreeVoxels voxels.
    Verify each recovers ground truth orientation within 2 degrees.
    Persist a detailed validation report.

    Run with: uv run pytest tests/test_reconstruction_validation.py -m slow -v
    """

    MISORIENTATION_THRESHOLD_DEG = 2.0

    def _reconstruct_voxel(self, reconstructor, voxel, rng, cubic_symmetry_quats):
        """Reconstruct one voxel and compute misorientation vs ground truth."""
        vertices = _get_voxel_vertices(voxel)
        gt_orientation = voxel.orientation.copy()

        t0 = time.time()
        result = reconstructor.reconstruct_voxel(
            voxel_vertices=vertices,
            phase_index=voxel.phase,
            rng=rng,
        )
        elapsed = time.time() - t0

        # Compute misorientation
        q_result = matrix_to_quaternion(result.orientation)
        q_truth = matrix_to_quaternion(gt_orientation)
        misori_rad = get_misorientation(q_result, q_truth, cubic_symmetry_quats)
        misori_deg = math.degrees(misori_rad)

        # Axis-angle representations
        gt_axis, gt_angle = _rotation_axis_angle(gt_orientation)
        recon_axis, recon_angle = _rotation_axis_angle(result.orientation)

        # Convergence code
        if result.overlap_info is not None and result.overlap_info.hit_ratio >= 0.8:
            conv_code = "HIT_RATIO_CONVERGED"
        elif result.cost < 0.0001:
            conv_code = "COST_CONVERGED"
        elif result.cost < 1.0:
            conv_code = "PARTIAL"
        else:
            conv_code = "NOT_CONVERGED"

        return VoxelValidationResult(
            voxel_index=-1,  # set by caller
            ground_truth=gt_orientation,
            reconstructed=result.orientation,
            misorientation_deg=misori_deg,
            cost=result.cost,
            hit_ratio=result.overlap_info.hit_ratio if result.overlap_info else 0.0,
            quality=result.overlap_info.quality if result.overlap_info else 0.0,
            convergence_code=conv_code,
            elapsed_seconds=elapsed,
            gt_axis=gt_axis,
            gt_angle_deg=gt_angle,
            recon_axis=recon_axis,
            recon_angle_deg=recon_angle,
        )

    def _generate_report(self, results: List[VoxelValidationResult], total_time: float) -> str:
        """Generate human-readable validation report."""
        lines = []
        lines.append("=" * 78)
        lines.append("IceNine Python Reconstruction Validation Report")
        lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"Test case: ThreeVoxels (3 voxels, copper, cubic symmetry)")
        lines.append(f"Data: Synthetic (forward simulation output as experimental data)")
        lines.append(f"Search: Blind discrete search over 4886 FZ orientations + MC refinement")
        lines.append(f"Config: MaxLocalResolution=3, MaxMCSteps=200, Restarts=2, Candidates=30")
        lines.append("=" * 78)
        lines.append("")

        all_pass = True
        for r in results:
            passed = r.misorientation_deg < self.MISORIENTATION_THRESHOLD_DEG
            if not passed:
                all_pass = False
            status = "PASS" if passed else "FAIL"

            lines.append(f"--- Voxel {r.voxel_index} [{status}] ---")
            lines.append(f"  Time: {r.elapsed_seconds:.1f}s")
            lines.append(f"  Convergence: {r.convergence_code}")
            lines.append(f"  Cost: {r.cost:.6f}  |  Quality: {r.quality:.6f}  |  Hit ratio: {r.hit_ratio:.3f}")
            lines.append(f"  Misorientation: {r.misorientation_deg:.4f} deg (threshold: {self.MISORIENTATION_THRESHOLD_DEG} deg)")
            lines.append("")
            lines.append(f"  Ground truth orientation (rotation matrix):")
            lines.append(_format_matrix(r.ground_truth))
            lines.append(f"  Ground truth axis-angle: axis=[{r.gt_axis[0]:.5f}, {r.gt_axis[1]:.5f}, {r.gt_axis[2]:.5f}], angle={r.gt_angle_deg:.4f} deg")
            lines.append("")
            lines.append(f"  Reconstructed orientation (rotation matrix):")
            lines.append(_format_matrix(r.reconstructed))
            lines.append(f"  Reconstructed axis-angle: axis=[{r.recon_axis[0]:.5f}, {r.recon_axis[1]:.5f}, {r.recon_axis[2]:.5f}], angle={r.recon_angle_deg:.4f} deg")
            lines.append("")

        lines.append("=" * 78)
        lines.append(f"SUMMARY")
        lines.append(f"  Total voxels: {len(results)}")
        lines.append(f"  Passed: {sum(1 for r in results if r.misorientation_deg < self.MISORIENTATION_THRESHOLD_DEG)}/{len(results)}")
        lines.append(f"  Total time: {total_time:.1f}s")
        lines.append(f"  Mean misorientation: {np.mean([r.misorientation_deg for r in results]):.4f} deg")
        lines.append(f"  Max misorientation: {np.max([r.misorientation_deg for r in results]):.4f} deg")
        lines.append(f"  Mean cost: {np.mean([r.cost for r in results]):.6f}")
        lines.append(f"  Overall: {'ALL PASS' if all_pass else 'SOME FAILED'}")
        lines.append("=" * 78)

        return "\n".join(lines)

    def test_blind_reconstruction_all_voxels(
        self, recon_setup, ground_truth_mic, cubic_symmetry_quats, report_path
    ):
        """
        Reconstruct all 3 voxels from blind search and verify orientation recovery.

        This is the main validation test. For each voxel:
        1. Start from blind search (no knowledge of ground truth)
        2. Run full multi-level discrete search + MC optimization
        3. Compare recovered orientation to ground truth via misorientation angle
        4. Assert misorientation < 2 degrees
        5. Generate detailed report
        """
        reconstructor = BasicVoxelReconstructor(recon_setup)
        rng = np.random.default_rng(42)

        results = []
        total_t0 = time.time()

        for idx, voxel in enumerate(ground_truth_mic.voxels):
            print(f"\n  Reconstructing voxel {idx}...")
            r = self._reconstruct_voxel(reconstructor, voxel, rng, cubic_symmetry_quats)
            r.voxel_index = idx
            results.append(r)
            print(f"  Voxel {idx}: misorientation={r.misorientation_deg:.4f} deg, "
                  f"cost={r.cost:.6f}, time={r.elapsed_seconds:.1f}s")

        total_time = time.time() - total_t0

        # Generate and save report
        report = self._generate_report(results, total_time)
        report_path.write_text(report)
        print(f"\n  Validation report saved to: {report_path}")
        print(report)

        # Assertions
        for r in results:
            assert r.misorientation_deg < self.MISORIENTATION_THRESHOLD_DEG, (
                f"Voxel {r.voxel_index}: misorientation {r.misorientation_deg:.4f} deg "
                f"exceeds threshold {self.MISORIENTATION_THRESHOLD_DEG} deg"
            )
            assert r.cost < 0.5, (
                f"Voxel {r.voxel_index}: cost {r.cost:.4f} too high (expected < 0.5)"
            )
            assert r.hit_ratio > 0.3, (
                f"Voxel {r.voxel_index}: hit_ratio {r.hit_ratio:.3f} too low (expected > 0.3)"
            )
