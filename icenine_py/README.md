# IceNine Python/PyTorch Implementation

Python port of IceNine crystallography and diffraction primitives with PyTorch acceleration.

## Overview

This is a test-driven migration of IceNine's core diffraction physics calculations from C++ to Python, focusing on:

- **Crystallography primitives**: Symmetry operations, reciprocal lattice generation
- **Diffraction physics**: Scattering vector calculations, Bragg condition, observable peaks
- **Validation**: Every function tested against C++ ground truth (tolerance < 1e-6)

## Project Structure

```
icenine_py/
├── icenine/               # Python package
│   ├── constants.py       # Physical constants (KEV_OVER_HBAR_C_IN_ANG, etc.)
│   ├── symmetry.py        # Crystal symmetry operations (pymatgen wrapper)
│   ├── crystal_structure.py  # Crystal structures and reciprocal lattice
│   └── diffraction_core.py   # Core diffraction calculations
├── tests/                 # Test suite
│   ├── test_symmetry.py   # Symmetry operation validation
│   ├── test_miller_indices.py  # Miller index generation validation
│   └── test_diffraction.py     # Diffraction calculation validation
├── cpp_harness/           # C++ ground truth generator
│   └── generate_test_data.cpp
└── cpp_outputs/           # C++ test data (JSON)
```

## Installation

```bash
# Install dependencies
pip install -e .

# Install with development tools
pip install -e ".[dev]"
```

## Dependencies

- **pymatgen**: Crystallography and symmetry operations
- **PyTorch**: Tensor operations and GPU acceleration
- **NumPy**: Array operations
- **SciPy**: Rotation and transformation utilities

## Testing

All tests compare Python outputs against C++ ground truth:

```bash
# Run all tests
pytest tests/

# Run specific test suite
pytest tests/test_symmetry.py -v

# Run with coverage
pytest --cov=icenine tests/
```

## Test Parameters

Using Gold (Au) FCC from `ConfigFiles/ReconstructTest.config`:

- Element: Au
- Lattice: FCC (a = 4.0782 Å)
- Space group: 225 (Fm-3m)
- Beam energy: 50.02099 keV
- Wavelength: 0.2478 Å
- Max Q: 8.0 Å⁻¹

## Validation Strategy

1. **C++ Test Harness**: Generate ground truth data from original IceNine
2. **Python Implementation**: Port functions using modern libraries
3. **Numerical Comparison**: Assert Python matches C++ within tolerance
4. **Validation Report**: Statistical analysis of differences

## Citation

Based on: S. F. Li and R. M. Suter, "Adaptive reconstruction method for three-dimensional orientation imaging", *Journal of Applied Crystallography*, 2013.
