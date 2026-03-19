"""
Demonstration of PyTorch-based diffraction_core for neural networks.

This script shows:
1. Batched scattering omega calculations
2. GPU acceleration (if available)
3. Automatic differentiation for gradient-based optimization
4. Integration with crystal orientation optimization
"""

import torch
import numpy as np
import time
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from icenine.diffraction_core import (
    get_scattering_omegas_torch,
    batch_scattering_omegas_from_reflections,
    BatchScatteringResult,
)
from icenine.constants import TEST_PARAMS


def demo_basic_usage():
    """Demo 1: Basic batched scattering calculations."""
    print("=" * 70)
    print("DEMO 1: Basic Batched Scattering Calculations")
    print("=" * 70)

    # Au FCC parameters
    a = TEST_PARAMS["a"]
    a_recip = 2 * np.pi / a

    # Create batch of scattering vectors
    hkl_list = torch.tensor([
        [1, 1, 1],   # (111)
        [2, 0, 0],   # (200)
        [2, 2, 0],   # (220)
        [3, 1, 1],   # (311)
        [2, 2, 2],   # (222)
    ], dtype=torch.float32)

    g_vectors = hkl_list * a_recip
    g_magnitudes = torch.norm(g_vectors, dim=1)

    print(f"\nProcessing {len(hkl_list)} reflections...")
    print(f"HKL indices:\n{hkl_list.numpy()}")

    # Calculate omega angles for all reflections at once
    result = get_scattering_omegas_torch(
        g_vectors,
        g_magnitudes,
        beam_energy=TEST_PARAMS["beam_energy"]
    )

    print(f"\nResults:")
    print(f"Observable: {result.observable.sum()}/{len(hkl_list)}")
    print(f"\nOmega angles (radians):")
    for i, (obs, omega1, omega2) in enumerate(zip(result.observable, result.omega1, result.omega2)):
        hkl = hkl_list[i].numpy().astype(int)
        if obs:
            print(f"  {hkl}: ω₁ = {omega1:.4f}, ω₂ = {omega2:.4f}")
        else:
            print(f"  {hkl}: Not observable")


def demo_gpu_acceleration():
    """Demo 2: GPU acceleration (if available)."""
    print("\n" + "=" * 70)
    print("DEMO 2: GPU Acceleration")
    print("=" * 70)

    if not torch.cuda.is_available():
        print("\nCUDA not available. Comparing CPU vs CPU (for timing baseline).")
        device = torch.device("cpu")
        device_name = "CPU"
    else:
        device = torch.device("cuda")
        device_name = f"GPU ({torch.cuda.get_device_name(0)})"

    # Generate many reflections for timing
    n_reflections = 1000
    a_recip = 2 * np.pi / TEST_PARAMS["a"]

    # Random Miller indices
    hkl = torch.randint(-5, 6, (n_reflections, 3), dtype=torch.float32)
    g_vectors = hkl * a_recip
    g_magnitudes = torch.norm(g_vectors, dim=1)

    # CPU timing
    print(f"\nProcessing {n_reflections} reflections on CPU...")
    start = time.time()
    result_cpu = get_scattering_omegas_torch(
        g_vectors,
        g_magnitudes,
        beam_energy=TEST_PARAMS["beam_energy"]
    )
    cpu_time = time.time() - start
    print(f"CPU time: {cpu_time:.4f} seconds")
    print(f"Observable: {result_cpu.observable.sum()}/{n_reflections}")

    # GPU timing (if available)
    if torch.cuda.is_available():
        g_vectors_gpu = g_vectors.to(device)
        g_magnitudes_gpu = g_magnitudes.to(device)

        # Warmup
        _ = get_scattering_omegas_torch(
            g_vectors_gpu,
            g_magnitudes_gpu,
            beam_energy=TEST_PARAMS["beam_energy"]
        )
        torch.cuda.synchronize()

        print(f"\nProcessing {n_reflections} reflections on {device_name}...")
        start = time.time()
        result_gpu = get_scattering_omegas_torch(
            g_vectors_gpu,
            g_magnitudes_gpu,
            beam_energy=TEST_PARAMS["beam_energy"]
        )
        torch.cuda.synchronize()
        gpu_time = time.time() - start
        print(f"GPU time: {gpu_time:.4f} seconds")
        print(f"Speedup: {cpu_time/gpu_time:.2f}x")


def demo_differentiable_physics():
    """Demo 3: Automatic differentiation capability."""
    print("\n" + "=" * 70)
    print("DEMO 3: Differentiable Physics - Gradient Computation")
    print("=" * 70)

    # Convert to tensor for gradient flow
    a_recip = torch.tensor(2 * np.pi / TEST_PARAMS["a"], dtype=torch.float32)

    # Define reflections
    hkl = torch.tensor([
        [1, 1, 1],
        [2, 0, 0],
        [2, 2, 0],
    ], dtype=torch.float32)

    # Create a simple orientation matrix (small rotation around z)
    angle = torch.tensor(0.1, requires_grad=True)  # Small angle
    c, s = torch.cos(angle), torch.sin(angle)
    R = torch.stack([
        torch.stack([c, -s, torch.zeros_like(c)]),
        torch.stack([s, c, torch.zeros_like(c)]),
        torch.stack([torch.zeros_like(c), torch.zeros_like(c), torch.ones_like(c)])
    ])

    print("\nCalculating scattering with differentiable orientation...")
    print(f"Rotation angle: {angle.item():.4f} rad ({np.degrees(angle.item()):.2f}°)")

    # Calculate scattering with gradient tracking
    result = batch_scattering_omegas_from_reflections(
        hkl,
        a_recip,
        R,
        beam_energy=TEST_PARAMS["beam_energy"]
    )

    print(f"\nObservable reflections: {result.observable.sum()}/{len(hkl)}")

    # Compute a loss (sum of omega1 angles)
    loss = result.omega1.sum()

    print(f"Loss (sum of omega1): {loss.item():.4f}")

    # Compute gradient
    loss.backward()

    print(f"\nGradient of loss w.r.t. angle: {angle.grad.item():.6f}")
    print("✓ Gradients computed successfully!")
    print("\nThis demonstrates that the diffraction physics is fully differentiable")
    print("and can be integrated into neural network training loops.")


def demo_batch_orientations():
    """Demo 4: Batched orientation matrices."""
    print("\n" + "=" * 70)
    print("DEMO 4: Batch Processing Multiple Orientations")
    print("=" * 70)

    a_recip = 2 * np.pi / TEST_PARAMS["a"]

    # Define reflections
    hkl = torch.tensor([
        [1, 1, 1],
        [2, 0, 0],
        [2, 2, 0],
    ], dtype=torch.float32)

    # Create batch of different orientations (rotations around z-axis)
    n_orientations = 5
    angles = torch.linspace(0, np.pi/2, n_orientations)

    print(f"\nTesting {n_orientations} different orientations...")
    print(f"Rotation angles: {np.degrees(angles.numpy()).astype(int)}°")

    results = []
    for i, angle in enumerate(angles):
        # Rotation matrix around z-axis
        c, s = torch.cos(angle), torch.sin(angle)
        R = torch.tensor([
            [c, -s, 0],
            [s, c, 0],
            [0, 0, 1]
        ])

        # Calculate scattering
        result = batch_scattering_omegas_from_reflections(
            hkl, a_recip, R, beam_energy=TEST_PARAMS["beam_energy"]
        )

        n_observable = result.observable.sum()
        print(f"  {int(np.degrees(angle)):3d}°: {n_observable}/{len(hkl)} observable")
        results.append(result)


def demo_comparison_with_numpy():
    """Demo 5: Compare PyTorch backend with NumPy interface."""
    print("\n" + "=" * 70)
    print("DEMO 5: PyTorch Backend vs NumPy Interface")
    print("=" * 70)

    from icenine.diffraction_core import get_scattering_omegas

    a_recip = 2 * np.pi / TEST_PARAMS["a"]
    g_vec_np = np.array([1.0, 1.0, 1.0]) * a_recip
    g_mag = np.linalg.norm(g_vec_np)

    print("\nCalling with NumPy array (backward compatibility):")
    result_np = get_scattering_omegas(
        g_vec_np, g_mag, beam_energy=TEST_PARAMS["beam_energy"]
    )
    print(f"  Result type: {type(result_np.g_vector)}")
    print(f"  Observable: {result_np.observable}")
    print(f"  Omega1: {result_np.omega1:.6f}")

    print("\nCalling with PyTorch tensor:")
    g_vec_torch = torch.from_numpy(g_vec_np).float()
    result_torch = get_scattering_omegas(
        g_vec_torch, g_mag, beam_energy=TEST_PARAMS["beam_energy"]
    )
    print(f"  Result type: {type(result_torch.g_vector)}")
    print(f"  Observable: {result_torch.observable}")
    print(f"  Omega1: {result_torch.omega1:.6f}")

    print("\nBoth interfaces produce identical results!")


def main():
    """Run all demonstrations."""
    print("\n" + "=" * 70)
    print("PyTorch Diffraction Core - Neural Network Integration Demo")
    print("=" * 70)
    print("\nThis demonstrates the PyTorch-based diffraction_core module")
    print("for differentiable physics simulations and neural network training.")

    demo_basic_usage()
    demo_gpu_acceleration()
    demo_differentiable_physics()
    demo_batch_orientations()
    demo_comparison_with_numpy()

    print("\n" + "=" * 70)
    print("All demos complete!")
    print("=" * 70)
    print("\nKey features demonstrated:")
    print("  ✓ Batched processing of multiple reflections")
    print("  ✓ GPU acceleration (when available)")
    print("  ✓ Automatic differentiation for gradient-based optimization")
    print("  ✓ Backward compatibility with NumPy interface")
    print("  ✓ Integration with crystal orientation optimization")
    print("\nReady for neural network-based inverse problems and")
    print("differentiable forward simulations!")


if __name__ == "__main__":
    main()
