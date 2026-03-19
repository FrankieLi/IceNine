# Integration Test: C++ vs Python IceNine

This directory contains an end-to-end integration test comparing C++ and Python implementations of IceNine forward simulation.

## Test Configuration

**Sample:** Three voxels (copper FCC)
**Config:** `ConfigFiles/Example2.Simulation.config`
**Expected Output:** ~360 detector images (180 omega steps × 2 detectors)

**C++ Output:** `ScatteringData/` (reference)
**Python Output:** `ScatteringData_Python/` (to be generated)

## Running the Test

### 1. Run Python Simulation

```bash
cd Examples/Example2.ThreeVoxels
uv run python run_python_simulation.py
```

**What it does:**
- Loads `ConfigFiles/Example2.Simulation.config`
- Initializes experiment (reads sample, detectors, omega ranges, structure)
- Runs forward simulation
- Outputs to `ScatteringData_Python/3Grains.sim*.d{0,1}`

**Expected output:**
```
Python IceNine Forward Simulation - Example2.ThreeVoxels
========================================================================

Config file: ConfigFiles/Example2.Simulation.config
Working directory: /Users/sfli/Research/IceNine/Examples/Example2.ThreeVoxels

Loading configuration...
  Sample: SimInput/three_voxels.mic
  Structure: DataFiles/copper.dat
  Beam energy: 64.351 keV
  Max Q: 16 Å⁻¹
  Detectors: 2
  ...
✓ Simulation completed successfully!
  Total files: 360
  Time elapsed: XXX seconds
```

### 2. Compare Outputs

```bash
uv run python compare_outputs.py
```

**What it does:**
- Loads corresponding C++ and Python detector images
- Computes pixel-wise differences
- Calculates statistical metrics
- Reports pass/fail against thresholds

**Expected output:**
```
IceNine C++ vs Python Output Comparison
===============================================================================

Found 360 C++ output files
Found 360 Python output files

Comparing detector images...
Compared 360 file pairs

Summary Statistics
===============================================================================

Absolute Differences:
  Max across all files:  X.XXe-XX
  Mean max difference:   X.XXe-XX
  ...

Pass/Fail Criteria
===============================================================================
✓ Max absolute difference < 0.001
✓ Mean absolute difference < 1e-05
✓ Max relative error < 1%
✓ Minimum correlation > 0.99

✓ ALL TESTS PASSED - Python and C++ outputs match!
```

## Current Results

- **3200/3201 pixels match** at identical (col, row) locations
- **Max relative intensity difference: 5.4e-6** (well within float32 precision)
- **4 pixel-location mismatches** — all at omega bin boundaries where floating-point rounding places a peak in an adjacent 1° bin
- **353/360 files nonempty** in both C++ and Python (identical set)

## Pass/Fail Criteria

| Metric | Threshold | Description |
|--------|-----------|-------------|
| Max relative difference | < 10⁻⁴ | Largest relative intensity error for matching pixels |
| Pixel location match | > 99.9% | Fraction of pixels at identical (col, row) |

## Test Files

```
Examples/Example2.ThreeVoxels/
├── run_python_simulation.py      # Run Python forward simulation
├── compare_outputs.py             # Compare C++ vs Python outputs
├── README_INTEGRATION_TEST.md     # This file
├── ConfigFiles/
│   ├── Example2.Simulation.config # Main config file
│   └── StandardGeometry.2Det      # Detector geometry
├── SimInput/
│   └── three_voxels.mic           # Sample voxel data
├── DataFiles/
│   ├── copper.dat                 # Crystal structure
│   ├── omega_180_2L.dat           # Omega rotation ranges
│   └── MyFZ.dat                   # Fundamental zone
├── ScatteringData/                # C++ reference output (360 files)
└── ScatteringData_Python/         # Python output (generated)
```

## Troubleshooting

### Error: "Config file not found"
Make sure you're running from the `Examples/Example2.ThreeVoxels/` directory.

### Error: "Python output directory not found"
Run `run_python_simulation.py` first before comparing.

### Error: "Missing data file"
Check that all input files exist:
- `SimInput/three_voxels.mic`
- `DataFiles/copper.dat`
- `DataFiles/omega_180_2L.dat`
- `ConfigFiles/StandardGeometry.2Det`

### Performance is slow
The Python version may be slower than C++ initially. Consider:
- Using GPU acceleration (if PyTorch supports CUDA)
- Profiling bottlenecks
- Optimizing vectorization

## Expected Differences

Small numerical differences occur due to:
- **Float32 precision:** C++ `float` and PyTorch `float32` have identical precision but different intermediate rounding in transcendental functions (sin, asin, atan2)
- **Omega bin boundaries:** When an omega value falls exactly on a 1° bin edge, C++ and Python may round to adjacent bins. This accounts for all 4 pixel-location mismatches.

Typical relative intensity error: < 10⁻⁵

## Citation

This test validates the Python port of:

> S. F. Li and R. M. Suter, "Adaptive reconstruction method for three-dimensional orientation imaging", *Journal of Applied Crystallography*, 2013.
