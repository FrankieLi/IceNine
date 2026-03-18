# C++ Validation Status

**Date:** 2025-01-15
**Status:** ⚠️ **NOT YET VALIDATED AGAINST C++**

## Current Validation Status

### What We've Done ✅

1. **Unit tests**: 24/24 tests passing for internal consistency
2. **Algorithm verification**: Line-by-line comparison with C++ source code
3. **Formula matching**: All mathematical formulas match C++ implementation
4. **Euler angles**: Phase 1 validated euler_to_matrix_torch() against C++
5. **Roundtrip tests**: Lab ↔ detector ↔ pixel transformations verified

### What We Haven't Done ❌

**Direct numerical comparison**: Running identical inputs through both C++ and Python implementations and comparing outputs.

## Why This Matters

While we've verified the **algorithm** is correct by:
- Reading the C++ source code line-by-line
- Implementing the same formulas
- Verifying internal consistency

We haven't verified **numerical accuracy** by:
- Compiling the C++ code
- Running the same test inputs through both
- Comparing outputs digit-by-digit

## Validation Test Suite Created

I've created a validation test harness in `tests/`:

### Files

1. **`validate_cpp_detector.cpp`** - C++ test program
   - Creates detectors with various configurations
   - Tests coordinate transformations
   - Tests ray intersections
   - Outputs results in parseable format

2. **`validate_py_detector.py`** - Python test program
   - Mirrors the C++ tests exactly
   - Same inputs, same operations
   - Same output format

3. **`run_validation.sh`** - Automated comparison script
   - Compiles C++ test
   - Runs both C++ and Python tests
   - Compares outputs line-by-line
   - Reports differences

## How to Run Validation

```bash
cd /Users/sfli/Research/IceNine/icenine_py/tests

# Option 1: Automatic (recommended)
./run_validation.sh

# Option 2: Manual steps
# 1. Compile C++ test
g++ -std=c++11 \
    -I../../XDM++/libXDM \
    -I../../Src \
    validate_cpp_detector.cpp \
    ../../XDM++/libXDM/3dMath.cpp \
    ../../Src/Detector.cpp \
    -o validate_cpp_detector

# 2. Run C++ test
./validate_cpp_detector > cpp_detector_output.txt

# 3. Run Python test
python3 validate_py_detector.py > py_detector_output.txt

# 4. Compare
diff -u cpp_detector_output.txt py_detector_output.txt
```

## Expected Validation Tests

The validation suite tests:

### Test 1: Basic Detector at Origin
- Configuration:
  - Position: (100, 0, 0)
  - Orientation: Identity
  - Size: 2048 × 2048 pixels
  - Pixel size: 0.2 × 0.2 mm
  - Beam center: (1024, 1024) pixels

- Operations tested:
  - `lab_to_detector_coordinate()` at detector center
  - `lab_to_detector_coordinate()` at offset point
  - `lab_to_pixel()` transformations
  - `detector_to_lab_coordinate()` inverse
  - `pixel_to_lab_coordinate()` inverse
  - Basis vector retrieval

### Test 2: Rotated Detector
- Configuration:
  - Position: (150, 10, -5)
  - Orientation: Euler angles (10°, 5°, 0°)
  - Size: 1024 × 1024 pixels
  - Beam center: (512, 512) pixels

- Operations tested:
  - Coordinate transformations with rotation
  - Rotated basis vectors
  - Orientation matrix

### Test 3: Ray Intersection
- Configuration:
  - Detector at (100, 0, 0)
  - Ray from origin along +X axis

- Operations tested:
  - `intersect_ray()` for perpendicular ray
  - Hit point calculation
  - Hit pixel calculation

## Known Potential Differences

### 1. Floating Point Precision
- **C++**: Uses `Float` (typically `float` = 32-bit)
- **Python**: Uses `torch.float64` for validation (64-bit)
- **Impact**: Python should be MORE accurate

### 2. Pixel Rounding
- **C++**: Truncates to `int` in `ToRowPixel()`, `ToColPixel()`
- **Python**: Returns continuous `float` for differentiability
- **Impact**: ~0.5 pixel difference expected (documented in tests)

### 3. Matrix-to-Euler (Not Tested)
- **C++**: Has `GetEulerAngles()`
- **Python**: `matrix_to_euler_torch()` not implemented (not needed)
- **Impact**: None for current use cases

## Validation Results

**Status:** Not yet run

To validate:
```bash
cd /Users/sfli/Research/IceNine/icenine_py/tests
./run_validation.sh
```

Expected outcome:
- ✅ **Pass**: All outputs match within numerical precision
- ⚠️ **Minor differences**: Pixel rounding (expected, documented)
- ❌ **Fail**: Significant numerical differences (investigate)

## If Validation Fails

1. **Check compilation**: Ensure C++ compiles with correct flags
2. **Check precision**: Verify float vs double handling
3. **Check coordinate systems**: Verify J/K vs X/Y mapping
4. **Check Euler conventions**: Verify ZXZ vs ZYZ (C++ comment is misleading)
5. **Isolate failure**: Which test case fails? Which operation?

## Action Items

- [ ] Run validation test suite
- [ ] Document results
- [ ] Fix any discrepancies found
- [ ] Add validation to CI/CD pipeline (optional)

## References

- **C++ Source**: `Src/Detector.{h,cpp}`
- **Python Source**: `icenine_py/icenine/detector.py`
- **Phase 1 Validation**: `PHASE1_GEOMETRY_COMPLETE.md` (Euler angles validated)
- **Implementation Plan**: `DETECTOR_IMPLEMENTATION_PLAN.md`

---

**Note to maintainers**: This validation should be run at least once to verify correctness. It's especially important if:
- The C++ implementation is known to be production-tested
- Numerical accuracy is critical for reconstruction
- The Python code will be used in production

**Recommendation**: Run validation before using the Python implementation for scientific results.
