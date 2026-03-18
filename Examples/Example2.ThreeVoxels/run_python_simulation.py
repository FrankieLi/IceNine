#!/usr/bin/env python3
"""
Run Python IceNine forward simulation for Example2.ThreeVoxels.

Compares against C++ IceNine output in ScatteringData/.
Python output goes to ScatteringData_Python/.

Usage:
    cd Examples/Example2.ThreeVoxels
    python run_python_simulation.py
"""

import sys
from pathlib import Path
import time

# Add icenine_py to Python path
REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "icenine_py"))

from icenine.config_file import ConfigFile
from icenine.forward_simulation import ForwardSimulation


def main():
    """Run Python forward simulation."""

    # Set working directory
    example_dir = Path(__file__).parent
    config_path = example_dir / "ConfigFiles" / "Example2.Simulation.config"

    print("=" * 70)
    print("Python IceNine Forward Simulation - Example2.ThreeVoxels")
    print("=" * 70)
    print(f"\nConfig file: {config_path}")
    print(f"Working directory: {example_dir}\n")

    # Load configuration
    print("Loading configuration...")
    config = ConfigFile.from_file(str(config_path))

    print(f"  Sample: {config.sample_filename}")
    print(f"  Structure: {config.structure_filename}")
    print(f"  Beam energy: {config.beam_energy} keV")
    print(f"  Max Q: {config.max_q} Å⁻¹")
    print(f"  Detectors: {config.num_detectors}")
    print(f"  Original output: {config.out_file_basename}")

    # Modify output directory to ScatteringData_Python/
    original_basename = config.out_file_basename
    config.out_file_basename = "ScatteringData_Python/3Grains.sim"

    print(f"  Modified output: {config.out_file_basename}")
    print()

    # Create simulator
    print("Initializing simulator...")
    simulator = ForwardSimulation(config)

    # Run simulation
    print("=" * 70)
    print("Starting forward simulation...")
    print("=" * 70)
    start_time = time.time()

    try:
        images = simulator.simulate_detector_images(
            output_dir=example_dir / "ScatteringData_Python"
        )

        elapsed = time.time() - start_time

        print("=" * 70)
        print("✓ Simulation completed successfully!")
        print("=" * 70)
        print(f"  Omega steps: {len(images)}")
        print(f"  Detectors: {len(images[0]) if images else 0}")
        print(f"  Total files: {len(images) * len(images[0]) if images else 0}")
        print(f"  Time elapsed: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
        print(f"  Output directory: ScatteringData_Python/")
        print()

        # Compare file counts
        cpp_dir = example_dir / "ScatteringData"
        python_dir = example_dir / "ScatteringData_Python"

        cpp_files = list(cpp_dir.glob("*.d*"))
        python_files = list(python_dir.glob("*.d*"))

        print("File count comparison:")
        print(f"  C++ output: {len(cpp_files)} files")
        print(f"  Python output: {len(python_files)} files")

        if len(cpp_files) == len(python_files):
            print("  ✓ File counts match!")
        else:
            print(f"  ⚠ Mismatch: {len(python_files) - len(cpp_files):+d} files")

        print()
        print("Next step: Run compare_outputs.py to validate numerical agreement")

    except Exception as e:
        print("=" * 70)
        print("✗ Simulation failed!")
        print("=" * 70)
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
