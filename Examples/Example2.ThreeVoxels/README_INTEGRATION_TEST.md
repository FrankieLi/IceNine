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
python run_python_simulation.py
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
python compare_outputs.py
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

## Pass/Fail Criteria

| Metric | Threshold | Description |
|--------|-----------|-------------|
| Max absolute difference | < 10⁻³ | Largest pixel error across all images |
| Mean absolute difference | < 10⁻⁵ | Average pixel error |
| Max relative error | < 1% | Largest percentage error |
| Minimum correlation | > 0.99 | Pearson correlation coefficient |

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

Even with identical algorithms, small numerical differences may occur due to:
- **Floating point precision:** C++ uses `float`, Python uses `torch.float32`
- **Library differences:** Different implementations of transcendental functions
- **Rounding errors:** Accumulated over thousands of calculations

Typical acceptable differences: 10⁻⁶ to 10⁻⁵ (relative error)

## Citation

This test validates the Python port of:

> S. F. Li and R. M. Suter, "Adaptive reconstruction method for three-dimensional orientation imaging", *Journal of Applied Crystallography*, 2013.
