#!/usr/bin/env python3
"""
Single-peak diagnostic: trace one voxel + one reciprocal vector through
the full forward simulation pipeline, printing every intermediate value.

Use this to compare line-by-line against C++ output.

Usage:
    cd Examples/Example2.ThreeVoxels
    uv run python debug_single_peak.py
"""

import sys
from pathlib import Path
import numpy as np
import torch

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "icenine_py"))

from icenine.config_file import ConfigFile
from icenine.experiment_setup import XDMExperimentSetup
from icenine.forward_simulation import ForwardSimulation
from icenine.diffraction_core import (
    get_scattering_omegas_torch,
    get_reflected_ray_dir,
    build_reflected_ray,
    get_illuminated_pixel,
)
from icenine.constants import KEV_OVER_HBAR_C_IN_ANG


def main():
    import os
    example_dir = Path(__file__).parent
    os.chdir(example_dir)
    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"

    # Load config and experiment setup (same as ForwardSimulation.__init__ + simulate)
    config = ConfigFile.from_file(str(config_path))
    exp_setup = XDMExperimentSetup(config)
    exp_setup.initialize_experiment()

    detector_list = exp_setup.get_detector_list()

    # Initialize sample (same as ForwardSimulation.simulate_detector_images)
    from icenine.sample import Sample
    sample = Sample()
    exp_setup.initialize_sample(sample, detector_list[0])

    mic = sample.get_mic()
    structure_list = sample.get_structure_list()
    wavenumber = KEV_OVER_HBAR_C_IN_ANG * exp_setup.beam_energy
    beam_dir = torch.from_numpy(exp_setup.get_beam_direction()).float()

    print("=" * 70)
    print("Single-Peak Diagnostic")
    print("=" * 70)
    print(f"Beam energy: {exp_setup.beam_energy} keV")
    print(f"Beam direction: {beam_dir}")
    print(f"Wavenumber: {wavenumber:.6f} Å⁻¹")
    print(f"Num detectors: {len(detector_list)}")
    print(f"Num voxels: {len(mic.voxels)}")

    # Pick first voxel
    VOXEL_IDX = 0
    voxel = mic.voxels[VOXEL_IDX]

    print(f"\n--- Voxel {VOXEL_IDX} ---")
    print(f"  Position: {voxel.position}")
    print(f"  Side length: {voxel.side_length}")
    print(f"  Points up: {voxel.points_up}")
    print(f"  Generation: {voxel.generation}")
    print(f"  Phase: {voxel.phase}")
    print(f"  Orientation (matrix):\n{voxel.orientation}")

    # Get vertices
    fs = ForwardSimulation.__new__(ForwardSimulation)
    vertices = fs._get_voxel_vertices(voxel)
    print(f"\n  Vertices:")
    for i in range(3):
        print(f"    v{i}: {vertices[i].numpy()}")

    # Get reciprocal vectors for this voxel's phase
    crystal_structure = structure_list[voxel.phase]
    reciprocal_vectors = crystal_structure.get_reflection_vectors()
    print(f"\n  Num reciprocal vectors: {len(reciprocal_vectors)}")

    # Pick first reciprocal vector
    RECP_IDX = 0
    recp = reciprocal_vectors[RECP_IDX]
    voxel_orientation = torch.from_numpy(voxel.orientation).float()

    print(f"\n--- Reciprocal Vector {RECP_IDX} ---")
    print(f"  q_vec (crystal frame): {recp.q_vec}")
    print(f"  q_mag: {recp.q_mag}")
    print(f"  hkl: ({recp.h}, {recp.k}, {recp.l})")
    print(f"  intensity: {recp.intensity}")

    # Transform to lab frame
    g_hkl = torch.from_numpy(recp.q_vec).float()
    g_lab = voxel_orientation @ g_hkl
    print(f"  g_lab (after orientation): {g_lab.numpy()}")

    # Solve Bragg condition
    result = get_scattering_omegas_torch(
        g_lab.unsqueeze(0),
        torch.tensor([recp.q_mag]),
        exp_setup.beam_energy,
        exp_setup.get_beam_deflection_chi_laue(),
    )

    print(f"\n  Observable: {result.observable[0].item()}")
    if not result.observable[0]:
        print("  Peak not observable, trying next reciprocal vector...")
        # Try to find one that's observable
        for idx, rv in enumerate(reciprocal_vectors):
            g = voxel_orientation @ torch.from_numpy(rv.q_vec).float()
            r = get_scattering_omegas_torch(
                g.unsqueeze(0),
                torch.tensor([rv.q_mag]),
                exp_setup.beam_energy,
                exp_setup.get_beam_deflection_chi_laue(),
            )
            if r.observable[0]:
                RECP_IDX = idx
                recp = rv
                g_hkl = torch.from_numpy(recp.q_vec).float()
                g_lab = voxel_orientation @ g_hkl
                result = r
                print(f"  Found observable peak at index {idx}: hkl=({recp.h},{recp.k},{recp.l})")
                break
        else:
            print("  No observable peaks found!")
            return 1

    omega1 = result.omega1[0].item()
    omega2 = result.omega2[0].item()
    print(f"  Omega solutions: {np.degrees(omega1):.4f}°, {np.degrees(omega2):.4f}°")

    # Use first omega solution
    omega = omega1

    # Compute sin(2θ) and scattering direction
    sin_theta = recp.q_mag / (2.0 * wavenumber)
    sin_2theta = np.sin(2.0 * np.arcsin(sin_theta))
    scattering_dir = g_lab / torch.norm(g_lab)

    print(f"\n--- Projection at omega = {np.degrees(omega):.4f}° ---")
    print(f"  sin(theta): {sin_theta:.8f}")
    print(f"  sin(2theta): {sin_2theta:.8f}")
    print(f"  Scattering direction: {scattering_dir.numpy()}")

    # Save and rotate sample
    current_orientation = sample.get_orientation()
    print(f"  Sample orientation (before rotate): {current_orientation}")
    sample.rotate_z(omega)
    print(f"  Sample orientation (after rotate_z({np.degrees(omega):.4f}°))")

    # Get reflected ray direction
    reflected_dir = get_reflected_ray_dir(sample, scattering_dir, beam_dir)
    print(f"  Reflected ray direction: {reflected_dir.numpy()}")

    # Project each vertex onto each detector
    for det_idx, detector in enumerate(detector_list):
        print(f"\n  --- Detector {det_idx} ---")
        print(f"    Position: {detector._position.numpy()}")
        print(f"    Plane normal: {detector._detector_plane.normal.numpy()}")
        print(f"    Plane D: {detector._detector_plane.d.item():.6f}")
        print(f"    Lab basis J: {detector._lab_frame_basis_j.numpy()}")
        print(f"    Lab basis K: {detector._lab_frame_basis_k.numpy()}")
        print(f"    Coord origin (lab): {detector._lab_frame_coord_origin.numpy()}")

        all_hit = True
        pixels = []
        for vi in range(3):
            vertex = vertices[vi]
            ray = build_reflected_ray(sample, vertex, reflected_dir)
            print(f"\n    Vertex {vi}: {vertex.numpy()}")
            print(f"      Ray origin (lab): {ray.origin.numpy()}")
            print(f"      Ray direction: {ray.direction.numpy()}")

            hit, col, row = get_illuminated_pixel(detector, ray)
            print(f"      Hit: {hit.item()}, col={col.item():.2f}, row={row.item():.2f}")

            if not hit.item():
                all_hit = False
            else:
                pixels.append((col.item(), row.item()))

        if all_hit:
            print(f"\n    All vertices hit! Triangle pixels:")
            for i, (c, r) in enumerate(pixels):
                print(f"      v{i}: col={c:.2f}, row={r:.2f}")
        else:
            print(f"\n    Some vertices missed detector")

    # Restore sample orientation
    sample.set_orientation(*current_orientation)

    print("\n" + "=" * 70)
    print("Done. Compare these values with C++ output.")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
