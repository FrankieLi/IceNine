# IceNine Python/PyTorch Migration History

Consolidated record of the C++ → Python/PyTorch migration. For current package usage, see [README.md](README.md).

## Migration Approach

- **Test-driven against C++ ground truth**: Every function validated against C++ output (tolerance < 1e-6). Ground truth generated via C++ test harnesses, stored as JSON in `cpp_outputs/`.
- **Dual-backend architecture**: `diffraction_core.py` supports both NumPy (backward compatible) and PyTorch (batched/GPU/autograd). Batched operations achieve 5-167x speedup over single-element processing.
- **Incremental porting**: Each module ported independently with its own test suite before integration.

## Completed Phases

### Phase 1: Constants & Symmetry
- Ported `PhysicalConstants.h` → `constants.py`
- Crystal symmetry via pymatgen wrapper (`symmetry.py`), validated against 24 cubic rotation matrices

### Phase 2: Crystal Structure & Diffraction Physics
- `crystal_structure.py`: Reciprocal lattice generation, Miller index enumeration (136 reflections → 9 unique after cubic symmetry reduction)
- `diffraction_core.py`: Scattering vectors, Bragg condition, observable peaks. Full PyTorch autograd support for gradient-based orientation optimization.

### Phase 3: MIC File I/O & Spatial Indexing
- `mic_file.py`: Read/write `.mic` files (Bunge Euler angles in degrees on disk, rotation matrices in radians in memory). Validated against all 37 repository `.mic` files (268K voxels).
- Spatial indexing via `scipy.spatial.cKDTree` (radius queries, k-nearest, boundary detection). 10-100x faster than pure Python KDTree with identical API.

### Phase 4: Geometry & Simulation Range
- `geometry.py`: Euler conversions, Plane/Ray classes, refactored from mic_file.py for reuse
- `simulation_range.py`: Omega range system for discontinuous data collection wedges

### Phase 5: Detector & ImageData
- `detector.py` (~700 lines): DetectorParameters dataclass, differentiable coordinate transforms, ray-detector intersection
- `image_data.py` (~1075 lines): Dense/sparse tensor modes, soft/hard rasterization, differentiable overlap computation

### Phase 6: Sample, ConfigFile, ExperimentSetup
- `sample.py`: Sample class with critical Euler convention handling (see gotchas below)
- `config_file.py`: Parser for 80+ parameters across 10 categories with automatic degree-to-radian conversion
- `experiment_setup.py`, `file_io.py`: Experiment configuration and file I/O utilities

### Phase 7: Simulation & Forward Simulation
- `simulation.py`: Core simulation engine
- `forward_simulation.py`: Forward model implementation
- `peak_filters.py`: Peak filtering utilities
- Integration test framework in `Examples/Example2.ThreeVoxels/`

## Critical Gotchas

### Detector Geometry (Root Causes of Forward Simulation Mismatch)

Two bugs in detector geometry caused zero pixel overlap between Python and C++ forward simulation output. Fixing them achieved ~86% pixel match.

**Bug 1: Detector plane normal computed incorrectly** (`detector.py:_calculate_image_plane`)
- **C++ method** (`Detector.cpp:211-244 CalculateImagePlane`): Defines the detector plane using 3 corner points `(1,0,0)`, `(0,1,0)`, `(0,0,0)` rotated by the orientation matrix and translated by position. The plane normal is `cross(p2-p1, p3-p1)`.
- **Python bug**: Was computing plane normal as `cross(lab_basis_j, lab_basis_k)`, which gives the wrong normal when J/K unit vectors are not the standard basis (e.g., for typical HEDM detectors where J=`(1,0,0)` and K=`(0,-1,0)`).
- **Fix**: Replicated the C++ 3-point plane construction exactly.

**Bug 2: Detector J/K unit vectors hardcoded incorrectly** (`detector.py:__init__`)
- **C++ method** (`DetectorFile.cpp:GetDetector`): Reads `JUnitVector` and `KUnitVector` from the detector file and passes them to `SetImageParameter`.
- **Python bug**: Hardcoded J=`(0,1,0)`, K=`(0,0,1)` instead of reading from file. The correct defaults matching C++ convention are J=`(1,0,0)`, K=`(0,-1,0)`.
- **Fix**: Added `j_unit_vector`/`k_unit_vector` parameters to `Detector.__init__()` with correct defaults, and pass parsed values from `file_io.py:DetectorInfo.get_detector()`.

**Related fix in `file_io.py`**: `LabFrameOrientation` Euler angles were being double-converted from degrees to radians — `euler_to_matrix()` internally converts degrees, but the code was passing already-converted radians.

### Fixed Bug: Sample orientation save/restore (degree/radian mismatch)
- `sample.get_orientation()` was returning Euler angles in **radians** but `sample.set_orientation()` expects **degrees**
- `rotate()`, `rotate_axis_angle()` were not updating `orientation_euler` after composition
- **Fix**: `orientation_euler` now stores degrees. `get_orientation()` returns degrees. `rotate()` and `rotate_axis_angle()` extract Euler angles after updating the matrix. `rotate_z()` does not update `orientation_euler` (hot-path, callers save/restore full matrix).

### Euler Angle Conventions
- **Voxel orientations**: Active ZXZ Bunge convention (standard crystallography)
- **Global sample orientation**: Passive Euler convention (`SetPassiveEulerMatrix` in C++), implemented as `passive_euler_matrix()` in Python. These are NOT the same — discovered through test failures.
- **Disk format**: `.mic` files store angles in **degrees**; conversion to radians happens at I/O boundaries.

### RNG Determinism
`GetRandomLocalGrid()` in C++ creates a local unseeded RNG per call, producing the same sequence each time. This behavior is preserved in the Python port — it may be intentional for reproducibility or an original oversight.

### Performance Requirements
Single-element physics operations are too slow for production. All diffraction calculations must use batched PyTorch operations for acceptable performance.

### Phase 8: Test Suite Hardening
Fixed all 12 test failures and 17 test errors to achieve a clean test suite (256 passed, 34 skipped, 0 failed, 0 errors).

**Bug fixes:**
- `symmetry.py`: `vectors_equivalent()` used `np.allclose` with default `rtol=1e-5`, causing the `tolerance` parameter to be ignored for large-magnitude vectors. Fixed by adding `rtol=0`.
- `experiment_setup.py`: Implemented `get_sample_symmetry()` (was returning `None`). Raises `NotImplementedError` for unsupported symmetry types (hexagonal).

**Test fixes:**
- `test_experiment_setup.py`: Added `chdir_to_project_root` fixture to match C++ CWD behavior for relative config paths. Tests that require data files (`omega_2L_100_cont.dat`) now skip gracefully.
- `test_simulation.py`: Fixed `get_reflections()` → `generate_reflections()` API mismatch, `r.g_vector` → `r.q_vec`.
- `test_experiment_setup.py`: Fixed `test_get_reciprocal_vector` to use actual beam direction from config (was hardcoded to `[0,0,1]`, config has `[1,0,0]`).
- `test_mic_file_validation.py`: Fixed histogram bins non-monotonic when `max(sizes)` is small.
- `test_diffraction.py`, `test_miller_indices.py`, `test_symmetry.py`: Added `pytest.skip` guards for missing C++ ground truth JSON files (previously caused `FileNotFoundError` errors).
- `test_miller_indices.py`, `test_symmetry.py`: Added `skipif` guards for optional `pytest-benchmark` dependency.

## Remaining Work

### Ported
- Forward simulation (serial + batched differentiable)
- Reconstruction: serial (BasicVoxelReconstructor), adaptive (AdaptiveVoxelReconstructor), BFS (BFSReconstruction)
- Cost functions (batched stages A-D with C extension)
- Search strategies: discrete grid search, zero-temperature MC, variance-minimizing MC

### Not Yet Ported
- Parallel processing (MPI-based server/client)
- Step size file reading (`_read_step_size_file` stub in `experiment_setup.py`)
- Crystal symmetry unique reflection list filtering

### Integration Testing — Pixel-Exact Match Achieved

C++ vs Python forward simulation comparison using three-voxel test case (360 detector images, 2 detectors, 180 omega steps, copper FCC).

**Final status:** 3200/3201 pixels match at identical locations with max relative intensity difference of 5.4e-6. Only 4 pixel-location mismatches remain, all caused by floating-point omega values landing on bin boundaries (rounding to adjacent 1° bin).

**Bug fix history (chronological):**
1. Detector plane normal computed incorrectly → fixed to match C++ 3-point construction
2. J/K unit vectors hardcoded wrong → fixed to read from detector file
3. `LabFrameOrientation` Euler angles double-converted → fixed to pass degrees to `euler_to_matrix()`
4. Voxel vertices used hardcoded 1mm right triangle → fixed to use actual equilateral triangle geometry from .mic file (side_length, points_up, position), matching C++ `MicIO.h:300-315`
5. **Rasterization mismatch (4.3x extra pixels)**: Python used soft/sigmoid rasterization. Replaced with C++-matching scanline rasterizer: float→int truncation matching `ToRowPixel`/`ToColPixel`, Sutherland-Hodgman polygon clipping with C++ boundary conventions, Bresenham edge table, scanline fill. Reduced from 4.3x → 1.36x pixel ratio.
6. **Multi-detector short-circuit**: Implemented C++ `GetProjectedVertices` behavior — if ANY vertex fails to intersect ANY detector plane, skip ALL detectors for that peak. (No effect on Example2 since all rays hit both detector planes.)
7. **C++ tokenizer bug (`Parser.cpp`)**: `Tokenize()` dropped the last line of files without trailing newlines. The mic file `three_voxels.mic` had no trailing `\n`, so the 3rd voxel was silently skipped — C++ only simulated 2 of 3 voxels. Fixed by treating end-of-buffer as implicit newline. This was the root cause of the remaining 843 extra pixels (Python correctly loaded all 3 voxels).

**Key implementation details (for future debugging):**
- `image_data.py:add_triangle_scanline()` — C++-matching rasterizer. Vertex coordinates are truncated (not rounded) to match `ToRowPixel`/`ToColPixel`. Clipping uses `x < x_max` (strictly less) for right/bottom and `x >= x_min` for left/top, matching C++ `SutherlandHodgman.h`.
- `simulation.py:project_voxel_multi_detector()` — Projects all 3 vertices onto all detectors before rasterizing. Short-circuits if any vertex misses any detector plane.
- `forward_simulation.py:_simulate_peaks()` — Triple loop: voxels × reflections × omega solutions. Uses `range_map.angle_to_wedge_index()` for omega bin mapping. Peak filter (`XDMEtaAcceptFn`) applies Lorentz-polarization correction: `I = I_form / (|sin(η)| × sin(2θ))`.
- Omega-to-bin mapping: `SimulationRange._build_index_list()` marks only the CENTER bin of each wedge, matching C++ `SimulationData.h:Set()`. For 180 individual 1° wedges, all 180 bins have valid indices.

### Performance Optimization — 6x Forward Simulation Speedup

Profiling the ManyGrains test (24,570 voxels, 112 reflections, 180 omega steps) revealed torch scalar operation overhead as the bottleneck. Each per-element torch tensor op has ~10-50µs Python overhead; with ~3,360 torch ops per voxel this dominated runtime at ~58ms/voxel (~22min total).

**Optimizations applied:**
1. **Batched Bragg solving**: Single `get_scattering_omegas_torch()` call per voxel with all reflections stacked as (N,3) tensor, instead of N individual calls
2. **Pre-computed detector geometry**: Detector plane normals, origins, unit vectors, and pixel parameters extracted as plain Python floats before the voxel loop
3. **Pure-Python inner loop**: Replaced torch scalar tensor ops with `math` module for rotation, ray-plane intersection, and lab-to-pixel conversion
4. **Functional rotation**: Compute Rz(omega) @ base_rotation as inline float math instead of mutating/restoring sample state per omega

**Result**: 9.3ms/voxel (down from 58.2ms), ~4.5min for full ManyGrains simulation.

**ManyGrains validation**: 99.97% pixel match rate (7,422,506 / 7,424,450), pixel count ratio 1.0000, total intensity ratio 1.000000, mean correlation 0.9996. Remaining ~0.03% mismatches are omega bin-boundary rounding (same as ThreeVoxels).

**Note on differentiability**: The inner loop now uses plain floats (not torch autograd). This is acceptable because forward simulation generates images for comparison and doesn't need gradients. See "Fully Batched Differentiable Pipeline" below for the autograd-preserving path.

**Diagnostic scripts in `Examples/Example2.ManyGrains/`:**
- `run_python_simulation.py` — runs full forward simulation (24,570 voxels)
- `compare_outputs.py` — statistical comparison with pass/fail criteria

**Diagnostic scripts in `Examples/Example2.ThreeVoxels/`:**
- `run_python_simulation.py` — runs full forward simulation
- `compare_outputs.py` / `compare_pixel_overlap.py` — pixel comparison
- `debug_single_peak.py` — traces one voxel + one reciprocal vector through full pipeline
- `debug_detector_geometry.py` — prints detector geometric properties
- See `Examples/Example2.ThreeVoxels/README_INTEGRATION_TEST.md`

### Fully Batched Differentiable Pipeline

Added `_simulate_peaks_batched()` as a dual path alongside the serial `_simulate_peaks()`. The batched pipeline processes ALL voxels simultaneously via large tensor operations, preserving torch autograd through stages 1-5 for future gradient-based orientation optimization.

**Architecture**: 7-stage pipeline with configurable `batch_size` chunking:
- **Stage 0**: Collect voxel orientations (V,3,3), vertices (V,3,3), detector geometry into contiguous tensors
- **Stage 1**: Batched Bragg solving — `bmm(orientations, g_hkl.T)` → `get_scattering_omegas_torch()` for all voxels×reflections at once
- **Stage 2**: Vectorized omega-to-wedge lookup via `SimulationRange.to_lookup_tensor()` + compact valid peaks with `torch.where()`
- **Stage 3**: Build `Rz(omega)` rotation matrices as (N,3,3) tensors, transform scattering directions, compute reflection vectors
- **Stage 4**: `batch_eta_filter()` — vectorized eta acceptance + Lorentz-polarization intensity correction
- **Stage 5**: Transform voxel vertices to lab frame, per-detector ray-plane intersection + pixel coordinate projection
- **Stage 6**: Sequential rasterization via existing `add_triangle_scanline()` (non-differentiable, accumulates into shared images)

**Differentiability boundary**: Stages 1-5 are fully differentiable via torch autograd. Stage 6 (scanline rasterization) is discrete/non-differentiable by design — gradients flow up to pixel coordinates but not through the integer rasterization step.

**API**: `simulate_detector_images(batched=True, batch_size=None)` selects the batched path. `batch_size=None` means all voxels in one chunk; set a value to limit memory for large samples.

**New helper functions**:
- `peak_filters.batch_eta_filter()` — vectorized eta filter + intensity computation
- `SimulationRange.to_lookup_tensor()` — omega-to-wedge index tensor for batched lookup

**Validation**:
- ThreeVoxels: pixel-exact match (3203/3203, max rel diff 2.1e-7)
- ManyGrains 100-voxel subset: 82,666/82,678 match, 12 bin-boundary mismatches (0.015%)
- ManyGrains full (24,570 voxels): 7,424,420 pixels in 3.3min (8.0ms/voxel)
- Autograd smoke tests: non-zero gradients verified through each differentiable stage

**Bug fixes during implementation**:
1. Omega bin lookup: `.long()` truncates negative floats like -0.001 to 0, falsely passing the `bin_idx >= 0` check. Fixed by checking `f_vals >= 0` on the float value before truncation.
2. Used float64 for omega-to-bin conversion to match serial path precision (reduces bin-boundary mismatches).

**Performance note**: For the current workload, rasterization (Stage 6) dominates at ~2.5M sequential `add_triangle_scanline` calls. The batched stages (1-5) add minimal overhead. Significant speedup would require a batched rasterizer or GPU parallelism.

**Regression test scripts**:
- `Examples/Example2.ThreeVoxels/test_batched_vs_serial.py` — pixel-exact comparison
- `Examples/Example2.ManyGrains/test_batched_vs_serial.py` — subset comparison + full benchmark

### Phase 8: Reconstruction Algorithm (Trivial/Serial)

Port of the C++ reconstruction pipeline: `SerialReconstruction` → `BasicVoxelReconstructor` → sequential voxel-by-voxel orientation search. No BFS, no MPI parallelism, no parameter optimization.

**New modules**:
- `sampling.py` (~570 lines): SO(3) uniform sampling via Sukharev grids (Yershova & LaValle, 2003). `CQuaternionGrid` port with barycentric-to-quaternion mapping on 4 hyperfaces of the upper hemi-hypersphere. Includes fundamental zone reduction (max |w| selection), misorientation computation, SLERP, structured/random local grid generation.
- `experimental_data.py` (~200 lines): Loads detector images from ASCII files or forward simulation output into `[omega_interval][detector]` array for pixel queries during reconstruction.
- `cost_functions.py` (~340 lines): `OverlapInfo` (Welford running mean quality metric), `VoxelCostFunction` (per-peak projection + overlap counting), `count_qualified_peaks` (contiguous detector validation). Supports both hard (discrete) and soft (differentiable) modes via existing `ImageData.get_triangle_overlap_property()`.
- `orientation_search.py` (~270 lines): `run_discrete_search` (FZ orientations × local grid perturbations), `MCOptimizer` (zero-temperature greedy descent with step halving and random restarts).
- `reconstructor.py` (~310 lines): `ReconstructionSetup` (loads all data from config), `BasicVoxelReconstructor` (multi-level adaptive: discrete → quick MC → filter → full MC → convergence), `SerialReconstruction` (voxel loop + .mic output).

**Architecture**:
```
SerialReconstruction.reconstruct_sample()
  └─ BasicVoxelReconstructor.reconstruct_voxel()
       └─ for level in [0..max_local_resolution]:
            1. run_discrete_search(fz_orientations × local_grid[level])
            2. MCOptimizer.optimize(candidates, 20 steps)  # quick
            3. Sort + keep top N
            4. MCOptimizer.optimize(top, full_params)       # full
            5. if hit_ratio_converged: break
```

**Key design decisions**:
- Direct port of Sukharev grid (no standard library equivalent for deterministic low-discrepancy SO(3) sampling)
- Uses existing `build_reflected_ray` + `get_illuminated_pixel` for projection (same code path as forward sim)
- FZ reduction uses proper rotation quaternions only (24 cubic ops filtered from 48 total)
- Cost function mode='hard' for C++ validation, mode='soft' for future gradient optimization

**Bugs found and fixed**:
- `reduce_to_fundamental_zone`: FZ reduction didn't enforce positive-w hemisphere after selecting the symmetry equivalent with max |w|. Added `if best[0] < 0: best = -best` (matches C++ `ToConvention`).
- `update_quality`: Used accumulated (running-sum) `self.pixel_overlap` / `self.pixel_on_detector` instead of per-peak values. C++ `UpdateQuality(oRHS, ...)` reads per-peak `oRHS.nPixelOverlap` / `oRHS.nPixelOnDetector`. Fixed by passing per-peak values as explicit arguments.
- Phase indexing: Phase 0 = empty space (no crystal), phase 1+ = fitted structures. `voxel.phase` maps directly into `structure_list` (no `- 1` offset).

**Integration tests** (`test_reconstruction_integration.py`, 5 tests):
- Forward sim output → experimental data → cost function → MC optimization → orientation convergence
- Tests: data dimensions, bright pixels, ground truth overlap (relative vs random), hit ratio, MC convergence from perturbed ground truth (< 10 deg misorientation)

**Test coverage**: 65 tests across 6 test files (all passing, 326 total suite):
- `test_sampling.py` (25 tests): Sukharev grid, SLERP, quaternion arithmetic, FZ reduction, misorientation
- `test_experimental_data.py` (8 tests): in-memory loading, ASCII file loading, pixel queries
- `test_cost_functions.py` (14 tests): OverlapInfo metrics, qualified peak counting
- `test_orientation_search.py` (8 tests): candidate sorting, parameters, convergence
- `test_reconstructor.py` (5 tests): voxel vertices, convergence codes
- `test_reconstruction_integration.py` (5 tests): end-to-end pipeline validation

### Performance: Cost Function Optimization

**Problem**: `get_triangle_overlap_property()` was the dominant bottleneck (73% of cost function wall time). For each peak projection, it allocated a full 2048×2048 temporary `ImageData`, rasterized a tiny triangle (~10-50 pixels), then called `torch.sum()` over all 4M pixels twice (overlap mask + lit mask). With ~1,700 peaks per cost evaluation and hundreds of evaluations per voxel, this caused massive memory churn and wasted computation.

**Profiling** (10 cost evaluations, ground truth orientation):
| Bottleneck | Time | % |
|---|---|---|
| `torch.sum()` on 2048×2048 masks | 16.9s | 43% |
| `get_triangle_overlap_property` overhead | 11.8s | 30% |
| `torch.zeros(2048, 2048)` temp alloc | 3.0s | 8% |

**Fix**: Rewrote `get_triangle_overlap_property` to compute barycentric coordinates only within the triangle's bounding box. No temporary image allocation — directly index into the existing experimental image at lit pixel positions.

**Result**: Integration tests 381s → 41s (**9.4x speedup**). Per-evaluation time 4.0s → 0.38s (**10.5x speedup**). `torch.sum()` and `torch.zeros()` completely eliminated from top bottlenecks.

**Remaining hotspots** (for future optimization): ray-plane intersection (45%), sample rotation (14%), detector coordinate conversion (13%) — all per-peak Python loops that could be batched in a future pass.

### Performance: Cost Function Optimization v2 (2026-03-22)

**Problem**: After the bounding-box optimization above, the Python cost function was still 209x slower than C++ (8,824 us vs 42 us). Per-operation profiling identified three bottlenecks:

| Bottleneck | Time | % of total | vs C++ |
|---|---|---|---|
| Eta filtering loop (Python for-loop, 2K iterations) | 4,137 us | 47% | 2,018x slower |
| Per-peak overlap (Sutherland-Hodgman + Bresenham + pixel lookup) | 377 us/peak | 40% | 2,690x slower |
| Stage A-C batched geometry | ~700 us | 8% | ~17x slower |

**Fix (3 priorities, implemented as task branches)**:

1. **Vectorize eta filtering** (`task/vectorize-eta-filter`): Replaced the Python for-loop (scalar `math.cos`/`math.sin`, per-iteration 3x3 tensor creation, `.item()` calls) with batched PyTorch operations. All 2K omega solutions processed in a single `torch.bmm` + `torch.atan2` pass. Files: `cost_functions.py` lines 718-760.

2. **C extension for triangle rasterization** (`task/c-extension-rasterizer`): Created `_rasterize.c` CPython C extension implementing `triangle_overlap()` and `pixel_radius_overlap()` in C. Ports Sutherland-Hodgman polygon clipping + Bresenham scanline fill with zero-copy numpy array access via `PyArray_GETPTR2`. Integrated into `image_data.py` (`get_triangle_overlap_property` hard mode) and `cost_functions.py` (pixel_radius mode) with Python fallback when C extension unavailable.

3. **Pixel-radius C fast path** (same task): Replaced Python nested dx/dy loop in Stage D with `_c_pixel_radius_overlap()` C call for the common `pixel_radius > 0` code path.

**Result**: `evaluate()` dropped from 8,824 us to 3,589 us (**2.46x speedup**). Per-peak overlap 377 us → 342 us. Batched overlap 37,134 us (serial) → 3,437 us (batched + C extension). All 328 tests pass, 34 skipped.

**Remaining gap**: Python 3,589 us vs C++ 42 us (~85x). The dominant remaining cost is the sequential Stage D loop.

#### Stage D: Sequential Per-Peak Overlap (the bottleneck)

`calculate_diffraction_overlap_batched()` is organized into four stages:

| Stage | What | Mode | Time |
|-------|------|------|------|
| A | Map omegas → wedge indices, filter invalid bins | Batched (numpy) | ~50 us |
| B | Batch rotation (Rz @ base_rot), reflection, vertex transform | Batched (torch.bmm) | ~300 us |
| C | Batch ray-detector intersection → pixel coordinates + hit masks | Batched (torch) | ~350 us |
| **D** | **Loop over M peaks × N detectors: look up experimental image, check pixel overlap, accumulate counters** | **Sequential (Python loop)** | **~2,900 us** |

Stage D (cost_functions.py ~line 523) iterates over each valid peak `p` in `range(M)` and each detector `det_idx` in `range(n_detectors)`. For each (peak, detector) pair it:

1. **Looks up the experimental image** via `exp_data.get_image(wedge_idx, det_idx)` — each peak maps to a different omega wedge, so a different experimental image
2. **Checks pixel overlap** using one of two modes:
   - `pixel_radius > 0`: Searches a ±radius square around the projected center pixel for any bright pixel (C extension fast path or Python fallback)
   - `pixel_radius == 0`: Full triangle rasterization via `get_triangle_overlap_property()` — Sutherland-Hodgman clip + Bresenham scanline + per-pixel overlap count (C extension fast path or Python fallback)
3. **Accumulates counters**: `detector_lit[det_idx]`, `spot_overlap[det_idx]`, `peak_pixel_overlap`, `peak_pixel_on_det`
4. **Calls `count_qualified_peaks()`** to determine if the peak qualifies (requires overlap on at least one detector)

**Why Stage D is inherently sequential**: Each peak projects onto a *different* experimental image (different omega wedge). The images cannot be stacked into a single tensor because they correspond to different physical measurements. The C++ handles this with a tight compiled loop (~0.2 us per peak); Python pays ~13 us per (peak × detector) iteration due to `.item()` calls, torch↔numpy conversions, and Python object dispatch.

**Future optimization paths**: (1) ~~Move entire Stage D loop to C extension~~ (done, see v3 below), (2) Cython/numba JIT, (3) ~~Pre-group peaks by wedge index~~ (done, see v3 below), (4) ~~Cache image conversions~~ (done, binary cache).

### Performance: Cost Function Optimization v3 — Batch C Extension (2026-03-22)

**Problem**: After v2, `evaluate()` was 3,589 us (~85x vs C++ 42 us). Stage D's sequential Python loop was the dominant bottleneck (~2,900 us of 3,589 us) due to M×N Python→C round trips, `.item()` calls, and repeated torch→numpy conversions.

**Fix (3 components)**:

1. **Binary image cache** (`ImageData._binary_cache`): Added `get_binary_numpy() -> np.ndarray` that returns a cached C-contiguous uint8 array where 1 = pixel > 0. Reconstruction only checks binary overlap, never intensity. The cache is lazily computed and invalidated by all mutator methods. Pre-populated at `ExperimentalData` load time via `prepare_for_reconstruction()`. Memory: 4x savings (uint8 vs float32).

2. **Batch C extension `stage_d_overlap()`** (`_rasterize.c`): Replaced the entire Stage D Python loop with a single C call. The function:
   - Receives all peak data as numpy arrays (wedge indices, detector hit masks, pixel coords/vertices)
   - Groups peaks by wedge index internally, processing all peaks against each binary image
   - Implements `count_qualified_peaks()` contiguity check and Welford running-mean quality update in C
   - Uses internal helpers `triangle_overlap_uint8()` and `pixel_radius_overlap_uint8()` that operate on binary images
   - Eliminates ~220 Python→C round trips per evaluation

3. **Differentiability documentation**: Added comment at Stage D explaining the intentional autograd break. Stages A-C preserve the PyTorch autograd graph; Stage D is inherently non-differentiable (binary pixel test, integer counting, contiguity validation). No code path calls `.backward()` — reconstruction uses discrete search + MC.

**Result**: `evaluate()` dropped from 3,589 us to 503 us (**7.1x speedup**). Batched overlap: 3,437 us → 364 us (**9.4x speedup**). Gap vs C++ reduced from 85x to ~12x.

| Component | v2 (us) | v3 (us) | Speedup |
|-----------|---------|---------|---------|
| Observable peaks (eta filter) | ~690 | ~690 | — |
| Stage D (overlap loop) | ~2,900 | ~170 | 17x |
| Data prep for C | 0 | ~50 | — |
| **evaluate() total** | **3,589** | **503** | **7.1x** |
| **Gap vs C++ (42 us)** | 85x | **12x** | — |

**Remaining gap**: The 12x gap is from Stages A-C (batched torch tensor ops vs C++ inline compiled code). Closing that would require moving A-C to C/CUDA — a separate effort.

**Files changed**: `image_data.py` (binary cache), `experimental_data.py` (eager cache population), `_rasterize.c` (batch `stage_d_overlap` + internal helpers), `cost_functions.py` (C batch path + Python fallback + autograd comment), `tests/test_cost_functions.py` (C-vs-Python equivalence test). All 329 tests pass.

### Validation: End-to-End Reconstruction (C++ vs Python)

**Date**: 2026-03-22

Ran identical reconstruction on ThreeVoxels example (MaxQ=8, 180 omega × 2 detectors, 4886 FZ orientations, 4 resolution levels 0-3, MaxMCSteps=200, SuccessiveRestarts=2) using same config (`ReconstructBenchmark.config`) and same experimental data (`ScatteringData/`).

**Timing:**

| Metric | C++ | Python |
|--------|-----|--------|
| Total wall time | 24.3s | 6637.5s |
| Data loading | ~1s | 1.3s |
| Reconstruction | ~23s | 6635.8s |
| Per-voxel average | ~8s | 2211.9s |
| Reconstruction slowdown | 1× | ~276× |

**Reconstruction quality:**

| Voxel | C++ Cost | Python Cost | Python Hit Ratio | Python Euler (recon) | Ground Truth Euler | Python Misori |
|-------|----------|-------------|------------------|----------------------|--------------------|---------------|
| 0 | 0.172 | 0.057 | 1.000 | (355.42, 5.19, 29.32) | (355.43, 5.19, 29.32) | 0.01° |
| 1 | 0.080 | 0.201 | 0.828 | (155.62, 45.18, 209.33) | (155.44, 45.18, 29.33) | ~0° (sym equiv) |
| 2 | 0.818 | 0.111 | 0.956 | (356.80, 3.70, 328.39) | (356.74, 3.70, 328.45) | 0.01° |

Python recovers all 3 orientations to within 0.01° of ground truth. Voxel 1's phi2 differs by 180° — this is a cubic symmetry equivalence. C++ finds a different local minimum for voxel 0 and fails on voxel 2 (cost=0.818).

**Timing breakdown (Python):** The dominant cost is level 3 discrete search: 4886 FZ × 512 local grid = 2,501,632 cost function evaluations per voxel. At ~400 us/eval, each level 3 takes ~1000s. MC optimization is negligible (<5s per voxel). All 3 voxels required all 4 levels before convergence.

**Progress statements**: Added `flush=True` progress reporting to `reconstructor.py` (per-level timing, candidate counts, convergence status) and `experimental_data.py` (image loading progress) for real-time visibility during long runs.

**Files changed**: `reconstructor.py` (progress statements), `experimental_data.py` (loading progress), `Examples/Example2.ThreeVoxels/run_python_reconstruction.py` (benchmark script).

### Cost Function C++/Python Equivalence Review

Systematic comparison of the C++ and Python cost function implementations. All code paths verified line-by-line.

**Equivalence table:**

| Aspect | C++ Source | Python Source | Match? |
|--------|-----------|---------------|--------|
| Scattering omegas (Bragg condition) | `Simulation.cpp:65-119` | `diffraction_core.py:59-194` | Yes |
| Observable peak generation | `Simulation.cpp:130-152` | `VoxelCostFunction.evaluate()` | Yes |
| Q-max filtering | `ExperimentSetup.cpp:SetDetectionLimit()` | `VoxelCostFunction.__init__` line 657 | Yes |
| Eta filtering | `Simulation.tmpl.cpp GetProjectedVertices` | `evaluate()` lines 743-750 | Yes (different location, same effect) |
| UpdateQuality (Welford mean) | `OverlapInfo.h:129-147` | `OverlapInfo.update_quality()` | Yes (fixed in Phase 8) |
| CountQualifiedPeaks | `OverlapInfo.tmpl.cpp:76-139` | `count_qualified_peaks()` | Yes |
| Quality = pixel_ratio × det_ratio | `OverlapInfo.h:136-138` | Lines 95-98 | Yes |
| PixelBasedPeakOverlapCounter | `OverlapInfo.h:336-382` | Lines 522-550 | Yes |
| ShapePixelOverlapCounter | `OverlapInfo.h:458-519` | Lines 551-569 | **Fixed** (detector_lit) |

**Bug found and fixed: `detector_lit` handling in triangle rasterization mode**

In `pixel_radius == 0` mode (full triangle rasterization), Python was setting `detector_lit[det_idx] = True` unconditionally whenever all 3 vertices intersected the detector plane (ray-plane hit). C++ (`ShapePixelOverlapCounter<SVoxel>`, OverlapInfo.h:498-517) only sets `bDetectorLit = true` when either:
- `nPixelOverlap > 0` (overlap found), OR
- At least one vertex is within detector bounds (`IsInBound`)

If all vertices project outside the detector image, C++ correctly leaves `bDetectorLit = false`. This affects `count_qualified_peaks` contiguous detector validation.

**Fix**: In both serial (`calculate_diffraction_overlap`, line 316) and batched (`calculate_diffraction_overlap_batched`, line 519) paths, moved `detector_lit` assignment after the overlap computation. Set it based on `n_lit_int > 0` (triangle has in-bounds pixels, equivalent to C++ `IsInBound`) or `n_overlap_int > 0`. The `pixel_radius > 0` branch was already correct (explicitly checks bounds via `found_in_bounds`).

**`GetConfidence` semantics (fixed)**: Python `OverlapInfo.confidence` now correctly returns `peak_overlap / peak_on_detector`, matching C++ `GetConfidence()`. Previously it returned `quality` (Welford mean).

### Q-max Behavior Clarification

Simulation data may use Q-max=16 Å⁻¹ while reconstruction uses Q-max=8 Å⁻¹. This does NOT degrade reconstruction quality. The cost function does not penalize under-observed peaks — peaks beyond max_q are filtered out at `VoxelCostFunction.__init__` time (line 657). Both C++ and Python handle this identically (C++ filters via `ExperimentSetup::SetDetectionLimit()`, Python via list comprehension on `rv.q_mag <= max_q`). Using fewer reciprocal vectors simply means fewer peaks contribute to the quality metric — faster evaluation but less discriminating.

### Cost Function Angular Sharpness

The cost function (pixel overlap quality) is extremely sharp in orientation space — quality drops rapidly within 0.2–0.5 degrees of the correct orientation. This is fundamental Bragg diffraction physics: diffraction peaks are narrow in angle, so even small orientation errors cause peaks to miss entirely. This sharpness motivates the multi-level adaptive search: coarse Sukharev grid sampling to find the correct basin, then fine MC refinement within it.

### Cost Function Validation — C++ vs Python (2026-03-22)

Apples-to-apples comparison using identical config (`ReconstructBenchmark.config`), identical data (`ScatteringData/`), and matched `eta_limit=86°`.

**Tool**: `icenine_py/benchmarks/validate_cost_function.py` (Python), `Src/CostFunctionBenchmark.cpp` (C++ with per-peak diagnostics).

**Results (voxel 2, ground truth orientation):**

| Metric | C++ | Python | Status |
|--------|-----|--------|--------|
| quality | 0.923077 | 0.916667 | 0.006 diff |
| pixel_overlap | 159 | 159 | Exact |
| pixel_on_detector | 159 | 160 | 1 pixel diff |
| peak_overlap | 52 | 52 | Exact |
| peak_on_detector | 52 | 52 | Exact |
| n_quality_points | 52 | 52 | Exact |

**Root cause of original 0.92 vs 0.67 gap (now resolved):**

1. **eta_limit not passed**: Python `VoxelCostFunction` defaulted to `π/2` (90°) instead of config's 86°. More peaks passed the eta filter, diluting quality. Fix: always pass `eta_limit=config.eta_limit`.
2. **Different data directory**: C++ reads `ScatteringData/` (no headers), Python was reading `ScatteringData_Python/` (with headers). Same pixel coordinates but different format.
3. **Different config file**: C++ benchmark used `ReconstructBenchmark.config` (MaxQ=8), Python used `Example2.Simulation.config` (MaxQ=16, but overridden by `max_q=8.0`).

**Remaining 0.006 quality difference**: A single peak at omega=-61.23° rasterizes 3 pixels in Python vs 2 in C++. One pixel sits on the triangle edge — borderline rounding in ray-plane intersection. This is within acceptable tolerance and does not indicate a formula bug.

### Adaptive BFS Reconstruction (2026-03-22 to 2026-03-23)

Ported the production C++ reconstruction pipeline to Python:

**Components ported:**

| C++ Source | Python Target | Description |
|-----------|---------------|-------------|
| `OrientationSearch.cpp:142-198` | `MCOptimizer._zero_temp_with_variance()` | Welford variance MC |
| `OrientationSearch.cpp:206-287` | `MCOptimizer.variance_minimizing_optimize()` | Adaptive sampling MC |
| `DiscreteAdaptive.tmpl.cpp:108-250` | `AdaptiveVoxelReconstructor.reconstruct_voxel()` | Multi-level adaptive search |
| `DiscreteAdaptive.tmpl.cpp:280-317` | `AdaptiveVoxelReconstructor.local_optimization()` | MC-only path for BFS neighbors |
| `BreadthFirstReconstructor.tmpl.cpp:112-188` | `BFSReconstruction._fit_from_seed()` | BFS spatial propagation |
| `ReconstructionStrategies.tmpl.cpp:334-356` | `BFSReconstruction._insert_seed()` | Orientation propagation to neighbors |
| `ReconstructionStrategies.h:257-263` | `ReconstructionState` | Voxel state machine |

**Key algorithm differences from BasicVoxelReconstructor:**

| Feature | Basic | Adaptive + BFS |
|---------|-------|----------------|
| FZ candidates between levels | Full FZ every level | Top 1/4 narrow between levels |
| Search diameter | Fixed | Shrinks ÷1.5 each level |
| nQMax | Fixed from config | Starts at 5 + min_level, increments |
| Quick MC | 20 steps, 0 restarts | 10 steps, 5 restarts |
| Final optimization | Full MC only | FindOptimal + VarianceMinimizing (σ²=0.02²) |
| Neighbor voxels | Independent full search | MC-only local_optimization from propagated orientation |
| Spatial propagation | None | BFS queue with 90% quality threshold |

**ThreeVoxels BFS benchmark:**

| Metric | Serial (Basic) | BFS (Adaptive) |
|--------|---------------|----------------|
| Seed voxels | 3 (all independent) | 2 (voxels 1, 2) |
| BFS neighbors | 0 | 1 (voxel 0, MC-only, 0.9s) |
| Total time | ~6638s | 252s |
| Speedup | 1× | ~26× |

Per-seed timing: 105s (voxel 2), 147s (voxel 1). The BFS benefit scales with sample size — for interior voxels in large samples (24K+), most get the cheap MC-only path (~1s each vs ~2200s for full search).

**Not yet ported**: C++ `Refit()` / `RESTART_FIT` — second-pass retry of REFIT voxels with LocalOptimization followed by full search if quality is still low.

### Bug Fix Batch (2026-03-23)

Five documented bugs fixed in one pass:

1. **Sample orientation degree/radian mismatch** (`sample.py`): `get_orientation()` returned radians but `set_orientation()` expects degrees — `set_orientation(*get_orientation())` silently produced wrong results. `rotate()` and `rotate_axis_angle()` didn't update `orientation_euler`. Fixed: `orientation_euler` now stores degrees, rotation methods extract Euler angles after matrix update, `rotate_z()` skips (hot-path, documented).

2. **OverlapInfo.confidence semantic mismatch** (`cost_functions.py`): `confidence` property returned `quality` (Welford mean) instead of `peak_overlap / peak_on_detector` (C++ `GetConfidence`). Fixed to match C++.

3. **Outdated "Remaining Work"** (`MIGRATION_HISTORY.md`): Listed reconstruction as "Not Yet Ported" despite being fully ported. Updated.

4. **_read_step_size_file print warning** (`experiment_setup.py`): Changed `print()` to `warnings.warn()` for proper Python warning semantics.

5. **Crystal symmetry TODO** (`experiment_setup.py`): Replaced bare TODO with explanation of why deferred (cubic symmetry doesn't need it).

### Adaptive Reconstruction Validation: C++ vs Python (Identical Algorithm) (2026-03-23)

The previous end-to-end reconstruction comparison (above) was **not an identical algorithm comparison**. It compared:
- **C++**: `DiscreteRefinement::ReconstructVoxel()` — the adaptive algorithm (FZ narrowing to top 1/4, shrinking diameter ÷1.5, nQMax incrementing, 10-step/5-restart quick MC, FindOptimal + VarianceMinimizing)
- **Python**: `BasicVoxelReconstructor` — a simpler algorithm (full FZ every level, fixed diameter, 20-step/0-restart MC)

This validation uses the **identical algorithm** on both sides: C++ `DiscreteRefinement` vs Python `AdaptiveVoxelReconstructor` (the actual port).

**Instrumentation:**
- **C++**: Added global evaluation counter `g_cost_eval_count` (defined in `CostFunctions.cpp`, incremented in `OverlapInfo.tmpl.cpp:CalculateDiffractionOverlap`). Counter reset and printed per-voxel in `SerialReconstruction.h` with `std::chrono` high-resolution timing. C++ `SerialReconstruction::ReconstructSample()` runs **both** `BasicVoxelReconstructor` and `DiscreteRefinement` on each voxel — eval counts are tracked separately.
- **Python**: Added `eval_count` to `VoxelCostFunction` (incremented per `evaluate()` call). Exposed via `AdaptiveVoxelReconstructor.last_eval_counts` property returning `(global_evals, local_evals, total_evals)`. Global cost function is recreated each level (different nQMax), so counts are accumulated across levels.

**Config**: `ReconstructBenchmark.config` (MaxQ=8, 180 omega × 2 detectors, 4886 FZ orientations, 4 levels 0-3, MaxMCSteps=200, SuccessiveRestarts=2, LocalGridRadius=5°, MCRadiusScaleFactor=1.0).

**Orientation Comparison:**

| Voxel | C++ Euler (adaptive) | Python Euler (adaptive) | Ground Truth Euler | C++ Misori | Py Misori |
|-------|---------------------|------------------------|--------------------|------------|-----------|
| 0 | (114.7, 87.5, 274.5) | (114.84, 87.47, 274.52) | (355.43, 5.19, 29.32) | ~90° (wrong basin) | ~90° (wrong basin) |
| 1 | (155.5, 45.2, 209.3) | (155.38, 45.18, 209.33) | (155.44, 45.18, 29.33) | ~0° (sym equiv) | ~0° (sym equiv) |
| 2 | (78.3, 48.6, 301.5) | (356.57, 3.70, 328.49) | (356.74, 3.70, 328.45) | FAILED (wrong) | 0.13° (correct) |

**Quality Metrics:**

| Voxel | C++ Cost | Py Cost | C++ Confidence | Py Confidence | C++ Hit Ratio | Py Hit Ratio |
|-------|----------|---------|----------------|---------------|---------------|-------------|
| 0 | 0.172 | 0.236 | 0.934 | 0.820 | 0.895 | 0.817 |
| 1 | 0.080 | 0.110 | 0.985 | 0.927 | 0.954 | 0.932 |
| 2 | 0.818 | 0.237 | 0.210 | 0.846 | 0.203 | 0.846 |

**Performance (Per-Voxel):**

| Voxel | C++ Time | Py Time | C++ Evals | Py Evals | C++ us/eval | Py us/eval |
|-------|----------|---------|-----------|----------|-------------|------------|
| 0 | 0.77s | 33.2s | 51,114 | 72,304 | 15.0 | 459.2 |
| 1 | 0.85s | 34.4s | 50,709 | 73,590 | 16.7 | 466.9 |
| 2 | 0.63s | 29.4s | 48,428 | 67,796 | 13.1 | 434.1 |

**Totals:**

| Metric | C++ | Python | Ratio |
|--------|-----|--------|-------|
| Total adaptive time | 2.25s | 97.0s | 43× |
| Total adaptive evals | 150,251 | 213,690 | 1.42× |
| Avg us/eval | 15.0 | 453.9 | 30× |

**Key findings:**

1. **Voxel 0**: Both C++ and Python find the same wrong local minimum (~90° from ground truth), confirming algorithm equivalence. This basin has good confidence (0.82-0.93) but is not the global optimum for this MaxQ=8 setting.

2. **Voxel 1**: Both find the same cubic symmetry equivalent orientation (~180° phi2 offset), confirming identical search behavior. C++ achieves slightly better quality metrics.

3. **Voxel 2**: Python **succeeds** (0.13° misori, cost=0.237) while C++ **fails** (cost=0.818). The difference is in candidate selection — Python evaluates more candidates (67K vs 48K evals) due to slight differences in floating-point intermediate values during discrete search, allowing it to discover the correct basin.

4. **Per-eval performance**: ~15 us (C++) vs ~450 us (Python) = **30× per-eval slowdown**. This is consistent with the Stage A-C batched torch operations overhead documented above (12× at the `evaluate()` level, plus Python-level overhead from the search/MC loops).

5. **Eval count difference**: Python performs ~42% more evaluations than C++. This is expected — floating-point differences in Bragg condition solving and eta filtering cause slightly different peaks to pass/fail, changing the number of candidates that enter MC refinement.

6. **Total wall time**: ~2.25s (C++) vs ~97s (Python) = **43× total slowdown** (lower than the previous 276× comparison because the adaptive algorithm does fewer evaluations than BasicVoxelReconstructor).

**Files changed**: `Src/CostFunctions.cpp` (global counter definition), `Src/OverlapInfo.tmpl.cpp` (counter increment), `Src/SerialReconstruction.h` (per-voxel timing + eval count reporting), `icenine_py/icenine/cost_functions.py` (eval_count tracking), `icenine_py/icenine/reconstructor.py` (expose eval counts), `Examples/Example2.ThreeVoxels/run_adaptive_reconstruction.py` (Python adaptive benchmark script).

### Spacing Filter Fix (2026-03-24)

**Problem**: Python did 42% more evaluations than C++ (213,690 vs 150,251). Root cause: C++ `GetSpacedCandidates` (DiscreteSearch.h:316-410) applies a per-clique angular spacing filter using `Acceptable()` (DiscreteSearch.h:244-258). Python was missing this filter entirely, causing candidate counts to balloon through adaptive levels (e.g., Level 3: 588 candidates in Python vs 4 in C++).

**The spacing filter**: For each FZ clique, after re-evaluating candidates with the local cost function (pixel_radius=0, cost = 1 - confidence), the filter walks the sorted candidate list and rejects any candidate if there exists an already-accepted candidate within `angular_radius` (= search diameter) that has a better (lower) cost. Misorientation is computed with crystal symmetry reduction.

**Implementation**: Three new functions in `orientation_search.py`:
- `_spacing_filter()` — core reject/accept logic matching C++ `Acceptable()`
- `get_symmetry_quaternions()` — extracts proper rotation quaternions from `CrystalSymmetry`
- `run_discrete_search_spaced()` — replaces `run_discrete_search` in adaptive path, combining global screening + local re-evaluation + per-clique spacing into one function

`reconstructor.py` `AdaptiveVoxelReconstructor.reconstruct_voxel()` updated to call `run_discrete_search_spaced()` instead of `run_discrete_search()` + separate re-evaluation phase.

**Results after fix:**

| Metric | C++ | Python (before) | Python (after) | After/C++ Ratio |
|--------|-----|-----------------|----------------|-----------------|
| Total time | 2.25s | 97.0s | 66.0s | 29× |
| Total evals | 150,251 | 213,690 | 157,726 | 1.05× |
| Avg us/eval | 15.0 | 453.9 | 418.7 | 28× |

Candidate counts per level now closely match C++. Voxel 2 now fails similarly to C++ (both cost ~0.8), confirming algorithmic alignment — the previous Python success was an artifact of excess candidates from the missing filter.

**Remaining bottleneck**: Per-eval gap is 28× (15 us C++ vs 419 us Python). This is pure computational overhead, not algorithmic difference.

### Differentiable Cost Function Infrastructure (2026-03-25)

New module `differentiable_cost.py` providing gradient-based orientation optimization infrastructure. Motivated by the goal of integrating the cost function into neural networks.

**Key components:**

| Class | Purpose |
|-------|---------|
| `ExperimentalImageStack` | Pre-stacks all (omega × detector) images into single contiguous tensor `(N, 1, H, W)` for batch `grid_sample`. Supports `binary=True` (0.0/1.0 matching existing pipeline) or `binary=False` (preserve intensities). Memory: ~5.6GB for 180×2×2048×2048 float32. |
| `MultiScaleImageStack` | Max-pool-downsampled image pyramid for coarse-to-fine optimization (later revised from an earlier Gaussian-blur/`F.conv2d` design — `max_pool2d` is morphological dilation for binary images and avoids the im2col memory blowup of chunked convolution). Widens angular basin from ~0.3° to ~2°+ at factor=8. |
| `DifferentiableCostFunction` | Replaces Stage D (sequential binary overlap counting) with differentiable bilinear sampling via `F.grid_sample`. Centroid point sampling (triangle centroid only, not full rasterization). |
| `DifferentiableOverlapInfo` | Dataclass with `quality`/`cost` tensors carrying `grad_fn` for backpropagation. |

**Architecture decisions:**

1. **Separate class, not refactored VoxelCostFunction**: DifferentiableCostFunction reimplements the A-C stages rather than sharing code, to avoid breaking the validated C++-matching pipeline.

2. **Centroid point sampling**: Each peak contributes one bilinear sample at `(v0+v1+v2)/3` instead of full triangle rasterization (~50 pixels). ~50x fewer samples, fully differentiable.

3. **Discrete routing**: Omega-to-wedge mapping (Stage A) and eta filtering use `.detach().numpy()` — no gradients through discrete routing. This is correct because the set of observable peaks doesn't change for small orientation perturbations.

4. **Multi-scale blur**: Broadens diffraction spots so gradient-based optimizers get non-zero signal even at ~2° offset (vs ~0.3° with unblurred images). Coarse-to-fine schedule: start blurry, progressively sharpen.

**Convenience method**: `ExperimentalData.to_image_stack(binary=True)` creates an `ExperimentalImageStack` directly.

**Tests (31 total, 29 passed, 2 skipped):**
- Phase 1 (20 tests): ExperimentalImageStack construction, binary/intensity modes, round-trip, flat indexing, batch gather, memory, sparse support. Gaussian kernel shape/normalization. MultiScaleImageStack blur preservation, peak spreading, shape matching.
- Phase 2 (11 tests): Gradient existence (`loss.backward()` → non-zero grads), finite-difference gradient check (autograd vs numerical Jacobian cosine similarity > 0.5), quality positive at ground truth (> 0.1), quality + cost = 1, correlation with hard cost function, wrong-orientation discrimination, all 3 ground truth voxels, invalid phase handling. Two blurred-scale tests skipped (require >16GB RAM).

**Files:** `differentiable_cost.py` (new, ~560 lines), `experimental_data.py` (added `to_image_stack()`), `tests/test_differentiable_cost.py` (new, 31 tests).

**Remaining work (Phases 3-4):**
- Phase 3: `GradientOrientationOptimizer` — axis-angle parameterization + coarse-to-fine Adam schedule
- Phase 4: Integration into `AdaptiveVoxelReconstructor` as optional gradient refinement after MC convergence

### MultiScaleImageStack Memory Optimization (2026-03-30)

Resolved ~9-11 GB memory usage during construction of multiple `MultiScaleImageStack` variants (needed for omega_window comparison benchmark). Three root causes and fixes:

**Root cause 1: `torch.maximum()` in `_omega_blend` allocates intermediate temp tensor each shift**
- Before: 3× tensor size overhead per shift (source, target, temp)
- Fix: `torch.maximum(a, b, out=a)` eliminates temp → 2× overhead
- `_prebuilt_downsampled` parameter

**Root cause 2: N independent `MultiScaleImageStack.__init__` calls each ran `_downsample_stack`**
- Each omega_window variant re-densified all 360 sparse frames independently
- Fix: `build_shared_base()` classmethod densifies once, `_prebuilt_downsampled` parameter shares result

**Root cause 3: `_downsample_stack` created new 16 MB dense tensor per frame**
- 360 frames × 16 MB = 5.6 GB allocator pool accumulation
- Fix: Pre-allocate `buf = torch.zeros(1, 1, H, W)`, reuse with `buf.zero_()`, wrapped in `torch.no_grad()` to prevent autograd graph accumulation

**Result**: Setup memory for 3 omega_window variants reduced from ~9-11 GB to ~2 GB for ManyGrains.

**New API** (added to `differentiable_cost.py`):
```python
# Build downsampled stacks once
shared_ds = MultiScaleImageStack.build_shared_base(image_stack, [1, 4, 8])
# Reuse across omega_window variants
for ow in [0, 1, 2]:
    ms = MultiScaleImageStack(image_stack, [1, 4, 8], omega_window=ow,
                              _prebuilt_downsampled=shared_ds)
```

**New benchmarks**:
- `benchmarks/bench_omega_window.py`: Compares quality vs misorientation for ω±0/1/2 at scale=2 (8× downsampled)
- `benchmarks/mem_profile_stack.py`: Progressive n_omega profiling (10→180) with `psutil` RSS measurement at each stage

### SparseImageStack (2026-03-30)

New class `SparseImageStack` in `differentiable_cost.py` stores only (row, col) pixel coordinates instead of dense float32 tensors. Loaded from `.d` files (ASCII, despite the extension) via `from_image_directory()`.

Memory comparison for ThreeVoxels 180 omegas × 2 detectors × 2048² images:
- Dense (`ExperimentalImageStack`): ~5.6 GB
- Sparse (`SparseImageStack`): ~12.5 KB (>400,000× smaller — binary diffraction images are extremely sparse)

The `_downsample_stack()` densifies frames on demand during construction; the downsampled result is stored dense since pooling fills in sparsity.

### Adam Gradient Optimization Benchmark (2026-03-30)

Benchmark `benchmarks/bench_gradient_optimization.py` tests Adam gradient descent on SO(3) via Lie algebra parameterization for orientation recovery. Sweeps: scale × omega_window × perturbation_deg.

**Parameterization**: `theta` (3-vector, axis-angle in radians) → `R = matrix_exp(skew(theta))` via `torch.matrix_exp`. Differentiable end-to-end through `DifferentiableCostFunction`.

**Result: Adam fails at all tested perturbation distances (1°, 2°, 5°)**

Root cause: `F.grid_sample` bilinear sampling of binary images gives **zero gradient in blob interiors** — only non-zero at the 1-pixel blob boundary. The gradient signal is structurally zero for the most important region (near ground truth orientation), making gradient descent unable to converge.

- Scale=1 (4× downsampled): Gradient non-zero up to ~0.3°; zero beyond
- Scale=2 (8× downsampled): Gradient up to ~0.5° with ω±1 blending; zero beyond
- Scale=0 excluded: Full-res densification takes ~277s per run; also zero basin past ~0.3°

**Output files**: `benchmarks/grad_opt_threevoxels.csv`, `benchmarks/grad_opt_manygrains.csv`, 6 PNG plots.

### CMA-ES Orientation Optimization Benchmark (2026-03-31)

Benchmark `benchmarks/bench_cmaes_optimization.py` tests CMA-ES (derivative-free Evolution Strategy) for orientation recovery, comparing hard cost vs diff cost (scale=2, ω±1).

**Algorithm**: CMA-ES on SO(3) using Lie algebra 3-vector parameterization. Initial step size `sigma0 = perturbation_rad` (full perturbation magnitude, not /3 — needed so initial population can spread back toward ground truth). All flat-landscape stop conditions disabled (`tolfun=0`, `tolfunhist=0`); relies on `tolx` for genuine convergence detection.

**ThreeVoxels results (3 voxels × 3 perturbations × 2 cost functions = 18 runs)**:

| | hard cost | diff_s2_ow1 |
|---|---|---|
| 1° | 0/3 recovered | 0/3 recovered |
| 2° | 0/3 recovered | 1/3 recovered (0.564°) |
| 5° | 0/3 recovered | 0/3 recovered |

**Hard cost**: 0/9 successful. Landscape completely flat outside ~0.5° basin — CMA-ES gets no signal, converges to random spurious orientations (40-165° misorientation).

**Diff cost**: 1/9 successful. Partial basin (~1-3°) in orientation space allows occasional recovery. However, false local optima from crystal symmetry (Cu has 24-fold cubic symmetry) trap CMA-ES at wrong orientations with quality 0.1-0.3 (vs ~0.4-0.7 at ground truth).

**Key finding**: Neither single-start Adam nor single-start CMA-ES can reliably recover orientations from perturbations ≥ 1°. The basin of attraction for both cost functions is too narrow relative to the search space. The original adaptive MC (random walk) approach is fundamentally more robust for this problem because it can make random moves without requiring gradient signal.

**Recommended directions**:
1. **Two-stage**: Use existing `AdaptiveMC` to get within ~0.5°, then apply gradient polishing (Adam or CMA-ES) — gradient IS reliable within the basin
2. **Multi-start CMA-ES**: Restart from N random initial orientations with a large sigma, take the best result — requires O(N) × 500 evaluations but doesn't need MC
3. **Distance field soft images**: Replace binary images with distance transform (distance to nearest bright pixel, stored float32) — gradient signal extends ~10-20px from blob, fixes zero-interior problem structurally

**ManyGrains results (20 voxels × 3 perturbations × 2 cost functions = 120 runs, 4.9 min)**:

| | hard cost | diff_s2_ow1 |
|---|---|---|
| 1° | 0/20 recovered | 1/20 recovered (vox 63: 0.063°) |
| 2° | 0/20 recovered | 1/20 recovered (vox 63: 0.056°) |
| 5° | 0/20 recovered | 0/20 recovered |

**Overall CMA-ES success rate**: Hard 0/120, Diff 3/120 (2.5%). The one voxel that consistently works (idx=63, φ1=322.68°, Φ=6.36°, φ2=43.72°, q=0.84) has an unusually favorable diff-cost landscape.

**Crystal symmetry false optima**: Diff cost frequently converges to orientations at 40-170° misorientation with quality 0.25-0.65 — these are genuine crystal symmetry equivalents (Cu has 24-fold cubic symmetry). Multi-start CMA-ES with symmetry-aware basin detection would be needed to handle this.

**Recommended directions** (confirmed by both benchmarks):
1. **Two-stage** (best near-term): Use existing `AdaptiveMC` to get within ~0.5°, then apply gradient polishing — gradient IS reliable within the basin
2. **Multi-start CMA-ES** with symmetry folding: Restart from each of the 24 cubic symmetry equivalents, take the best result
3. **Distance field soft images**: Replace binary images with distance transform — gradient extends ~10-20px from blobs, fixes zero-interior-gradient structurally

---

## Riemannian SO(3) Gradient Optimization Benchmark (2026-04-01)

**Branch**: `feature/differentiable-cost`
**Files**: `benchmarks/bench_riemannian_optimization.py` (NEW), `pyproject.toml` (added `riemannian` extra)

**Motivation**: The existing `euclidean_adam` optimizer from `bench_gradient_optimization.py` has two fundamental defects when optimizing over SO(3):
1. **Chart distortion**: Computing `d(cost)/d(theta)` (where `R = exp(skew(theta))`) mixes the Riemannian gradient with the Jacobian of `exp`, introducing distortion that grows as `R` drifts from `R_init`.
2. **Moment staleness**: Adam's first/second moment estimates accumulate in a fixed chart anchored at `R_init`. After many steps the moments reflect gradients at very different manifold points — they are never parallel-transported to the current `R`.

**Correct Riemannian approach**:
- Project Euclidean gradient `G = ∂cost/∂R` to tangent space T_R SO(3): `Ω = (R^T G − G^T R) / 2`
- Represent as 3-vector via `ω = [Ω₃₂, Ω₁₃, Ω₂₁]^T`
- Accumulate Adam moments in these so(3) coordinates (parallel transport is the identity under the left-trivialized flat connection — no moment rotation needed between steps)
- Retract via: `R_new = R · exp(-lr · skew(v_adam))` — stays exactly on SO(3)

**Three optimizers benchmarked**:

| Optimizer | Implementation | Description |
|---|---|---|
| `euclidean_adam` | Pure PyTorch | Baseline: Adam on θ ∈ ℝ³, R = matrix_exp(skew(θ)) |
| `riemannian_adam_manual` | Pure PyTorch | Manual Riemannian Adam: project grad, so(3) moments, matrix_exp retraction |
| `riemannian_adam_geoopt` | geoopt 0.5.1 | `geoopt.RiemannianAdam` on `Stiefel(3,3)` ≈ SO(3) |

Note: geoopt 0.5.x does not expose `SpecialOrthogonal`; `Stiefel(n=p=3)` is the orthogonal group O(3); geodesic retraction from a rotation matrix preserves det = +1 throughout.

**ThreeVoxels results (3 voxels × 3 perts × 2 scales × 3 ω-windows × 3 opts = 162 configs, n_steps=100, lr=0.01, 6.8 min)**:

Success defined as final misorientation < 0.5°:

| Optimizer | 1° pert | 2° pert | 5° pert |
|---|---|---|---|
| euclidean_adam | 7/18 | 2/18 | 1/18 |
| riemannian_adam_manual | 8/18 | 6/18 | 0/18 |
| riemannian_adam_geoopt | 8/18 | 2/18 | 3/18 |

Final misorientation (mean°):

| Optimizer | 1° pert | 2° pert | 5° pert |
|---|---|---|---|
| euclidean_adam | 1.80° | 3.25° | 5.09° |
| riemannian_adam_manual | 1.16° | 2.47° | 4.77° |
| riemannian_adam_geoopt | 1.81° | 2.19° | 5.79° |

**ManyGrains results (20 voxels × 3 perts × 2 scales × 3 ω-windows × 3 opts = 1080 configs, n_steps=100, lr=0.01, 8.1 min)**:

Success rate (< 0.5°):

| Optimizer | 1° pert | 2° pert | 5° pert |
|---|---|---|---|
| euclidean_adam | 24/120 (20%) | 25/120 (21%) | 4/120 (3%) |
| riemannian_adam_manual | 32/120 (27%) | 29/120 (24%) | 4/120 (3%) |
| riemannian_adam_geoopt | 37/120 (31%) | 29/120 (24%) | 5/120 (4%) |

Final misorientation (mean°):

| Optimizer | 1° pert | 2° pert | 5° pert |
|---|---|---|---|
| euclidean_adam | 1.91° | 2.40° | 5.99° |
| riemannian_adam_manual | 1.43° | 1.95° | 5.44° |
| riemannian_adam_geoopt | 1.31° | 1.96° | 5.66° |

**Key findings**:
- Both Riemannian optimizers outperform Euclidean Adam at 1° perturbation: +37% (geoopt) and +33% (manual) more successful recoveries on ManyGrains.
- At 2° perturbation the improvement is modest (+4%). At 5°, all methods fail equally — the basin of attraction problem dominates.
- The manual implementation confirms the theory: projecting gradients to T_R SO(3) at each step is more accurate than computing `d(cost)/d(theta)` in a fixed chart. The improvement is real but not dramatic — the gradient landscape is still too flat outside ~0.5° for reliable single-start optimization.
- `riemannian_adam_geoopt` and `riemannian_adam_manual` produce nearly identical results (mean 1.31° vs 1.43° at 1°), confirming the manual implementation is correct.

**geoopt API note**: `geoopt.manifolds.SpecialOrthogonal` does not exist in v0.5.1. Used `geoopt.manifolds.Stiefel()` instead; starting from SO(3) and using geodesic retraction keeps the matrix in SO(3) (det = +1 maintained throughout).

**Conclusions**: Riemannian structure helps at small perturbations but does not solve the flat-landscape / narrow-basin problem. The two-stage approach (AdaptiveMC → gradient polish) remains the recommended path.

---

## Riemannian SGD on SO(3) Benchmark (2026-04-01)

**Branch**: `feature/differentiable-cost`
**File**: `benchmarks/bench_sgd_optimization.py` (NEW)

**Motivation**: Compare SGD-family optimizers against Riemannian Adam on SO(3) to understand whether simpler mechanics (no second moment) are competitive, and to test whether Langevin noise (SGLD) provides probabilistic escape from the flat cost landscape.

**Five optimizers implemented** (all use geodesic retraction `R ← R·exp(−lr·Ω)`, same lr=0.01 as Adam benchmark):

| Optimizer | Implementation | Key difference from Adam |
|---|---|---|
| `riemannian_sgd_plain` | geoopt.RiemannianSGD, momentum=0 | No moments, constant lr |
| `riemannian_sgd_momentum` | geoopt.RiemannianSGD, β=0.9 | Heavy-ball momentum in so(3) |
| `riemannian_sgd_nesterov` | geoopt.RiemannianSGD, nesterov=True | Look-ahead momentum |
| `riemannian_sgd_cosine` | geoopt.RiemannianSGD + CosineAnnealingLR | lr_max=0.05 → lr_min=0.001 |
| `riemannian_sgld` | Pure PyTorch, Langevin noise | Isotropic noise on T_R SO(3), T annealing to 0 |
| `riemannian_adam_manual` | (baseline copy) | Adam, included for direct comparison |

**Results (ManyGrains, 20 voxels, n_steps=100, lr=0.01)**:

| Optimizer | 1° mean misori | 2° mean misori | 5° mean misori | 1° success (<0.5°) |
|---|---|---|---|---|
| riemannian_adam_manual | 1.43° | 1.95° | 5.44° | 32/120 (27%) |
| riemannian_sgld | 12.06° | 11.52° | 12.42° | 0/120 |
| riemannian_sgd_plain | 33.55° | 29.53° | 32.19° | 0/120 |
| riemannian_sgd_momentum | 124.36° | 127.25° | 127.09° | 0/120 |
| riemannian_sgd_nesterov | 124.15° | 121.65° | 127.89° | 0/120 |
| riemannian_sgd_cosine | 131.05° | 126.42° | 124.39° | 0/120 |

**Key finding: lr=0.01 is wildly incompatible with raw SGD on this problem.**

Adam's adaptive second-moment scaling effectively normalises the learning rate per-coordinate: `lr_eff ≈ lr / √m̂₂`. When gradients are small (which they are in the flat landscape most of the time), Adam inflates the effective lr to compensate, keeping steps bounded. Raw SGD applies lr directly to the gradient magnitude — when the cost landscape occasionally produces a large gradient near a blob boundary, the step `lr·Ω` overshoots and sends R far from its starting point. The momentum variants amplify this: a single large gradient gets accumulated into `m`, carried forward, and compounds over many steps → divergence to 100°+ misorientation.

**SGLD result**: Even with Langevin noise, lr=0.01 causes divergence (mean misorientation 12°). The noise helps slightly vs. plain SGD by providing random exploration, but the gradient component itself overshoots. A proper SGLD implementation would require a much smaller lr (≈10× smaller) to keep gradient steps bounded, with correspondingly larger noise scale to explore.

**Conclusion**: For this cost function, the lr chosen for Adam (0.01) is completely inappropriate for SGD without adaptive scaling. A fair comparison would require lr-tuning per optimizer (e.g. lr≈0.001 for SGD). The key takeaway is that **Adam's adaptive scaling is not just a convenience — it is essential for this problem** because:
1. Gradients are nearly zero almost everywhere (flat landscape), so the effective lr from Adam's `1/√m̂₂` normalization compensates correctly
2. Near blob boundaries where gradients are large, Adam's `√m̂₂` denominator attenuates steps, preventing overshoot
3. Raw SGD has no such self-regulation and diverges on this problem at standard learning rates

**Future direction**: Tune lr for SGD variants (try 0.0001–0.001) and re-run to get a fair comparison. The cosine schedule variant is most promising as it starts with a larger lr for exploration and decays to a small lr for refinement.

---

## Riemannian SGD Fair Comparison (lr=0.001 for SGD, lr=0.01 for Adam) (2026-04-01)

**Branch**: `feature/differentiable-cost`
**New files**: `sgd_opt_{example}_sgdlr1e3.csv`, `sgd_opt_*_sgdlr1e3.png`
**Code change**: Added `--sgd-lr` / `--adam-lr` argparse flags + `optimizer_lrs` dict threading through `sweep_voxel_multi` and `run_example` for per-optimizer lr override.

Re-ran all 6 optimizers with SGD variants at lr=0.001 (10× smaller than Adam's lr=0.01).

**ManyGrains results (20 voxels, n_steps=100, adam lr=0.01, sgd lr=0.001)**:

| Optimizer | 1° mean | 2° mean | 5° mean | 1° success (<0.5°) |
|---|---|---|---|---|
| riemannian_adam_manual (lr=0.01) | 1.43° | 1.95° | 5.44° | 32/120 (27%) |
| riemannian_sgld (lr=0.001) | 1.66° | 2.22° | 5.49° | 18/120 (15%) |
| riemannian_sgd_plain (lr=0.001) | 2.45° | 2.53° | 4.96° | 1/120 (<1%) |
| riemannian_sgd_momentum (lr=0.001) | 26.11° | 21.67° | 18.80° | 0/120 |
| riemannian_sgd_nesterov (lr=0.001) | 24.36° | 19.64° | 16.18° | 0/120 |
| riemannian_sgd_cosine (lr=0.001) | 26.11° | 21.67° | 18.80° | 0/120 |

**Key findings**:

1. **SGLD is now competitive**: At lr=0.001, SGLD (1.66° mean) approaches Adam (1.43° mean) and achieves 15% success rate. The Langevin noise provides some useful exploration that pure gradient descent misses. However, Adam still wins by a significant margin — its 27% success rate is nearly 2× SGLD's 15%.

2. **Momentum/Nesterov/cosine still diverge at lr=0.001**: These optimizers still reach 20–26° mean misorientation. The momentum accumulation (β=0.9 means ~10 steps of gradient are summed) amplifies any large gradient encountered, and the effective accumulated step is still too large even at lr=0.001. A fair lr for momentum would be ~0.0001.

3. **Plain SGD is stable but slow**: At lr=0.001, plain SGD reaches mean 2.45° (vs Adam 1.43°) with near-zero success rate. It moves too slowly to reach the basin in 100 steps — Adam's adaptive scaling effectively gives it ~10× larger effective lr in flat regions.

4. **cosine SGD = plain SGD here**: The cosine schedule starts at lr=0.001 and decays to lr_min=0.001/50=~0.00002 — effectively too slow throughout. A fair cosine run would need lr_max ≈ 0.005-0.01 for SGD.

**Overall conclusion**: Adam's adaptive second-moment normalization is doing critical work on this problem:
- In the flat landscape (most of the time): `lr_eff = lr/√m̂₂` inflates the effective lr when gradients are consistently small, moving faster than plain SGD at the same nominal lr
- Near blob boundaries: `√m̂₂` attenuates large gradient spikes, preventing divergence
- No single fixed lr for SGD captures both behaviors simultaneously

**Best single-start optimizer ranking**: riemannian_adam_geoopt ≈ riemannian_adam_manual > riemannian_sgld (lr=0.001) > riemannian_sgd_plain (lr=0.001) >> momentum variants

---

## Comprehensive HP Sweep — Gradient Methods vs. MC (2026-04-01)

**Branch**: `feature/differentiable-cost`
**Files added**:
- `benchmarks/bench_hp_sweep.py` (NEW) — full HP sweep benchmark
- `icenine/orientation_search.py` (MODIFIED) — `MCOptimizer.optimize()` gains optional `trajectory` parameter

**Motivation**: Prior benchmarks fixed hyperparameters and showed large lr sensitivity. This sweep exhaustively studies all hyperparameters for each optimizer family and compares against `MCOptimizer` to establish a definitive head-to-head comparison. Also captures optimization trajectories (angular step sizes per optimizer step) to understand search dynamics.

### Optimizers and HP grids

| Optimizer | HP Grid | Total configs |
|---|---|---|
| `riemannian_adam_geoopt` | lr ∈ {1e-4,...,0.1} × n_steps ∈ {100,200,500} × beta1 ∈ {0.9,0.95} | 42 |
| `riemannian_adam_manual` | Same as geoopt | 42 |
| `riemannian_sgd_plain` | lr ∈ {1e-4,...,0.01} × n_steps ∈ {100,200,500} | 15 |
| `riemannian_sgd_momentum` | lr ∈ {1e-5,...,1e-3} × momentum ∈ {0.5,0.9,0.99} (n_steps=200) | 15 |
| `riemannian_sgld` | lr ∈ {1e-4,...,0.01} × T_init ∈ {0.001,0.01,0.1} × n_steps ∈ {100,200,500} | 45 |
| `mc_optimizer` | max_mc_steps ∈ {100,500,1000,3500} × restarts ∈ {0,2,5} × step_frac ∈ {0.25,0.5,1.0} | 36 |

Total: 195 HP configs per voxel × perturbation. All gradient runs use `scale=2, omega_window=1`.

### Scale

- ThreeVoxels: all 3 qualifying voxels × 3 perturbations × 195 HP configs = ~1,755 runs
- ManyGrains: 100 voxels × 3 perturbations × 195 HP configs = ~58,500 runs (overnight)

### Output CSV schema

Two CSVs per example:
- `hp_sweep_{example}.csv` — one row per run with all HP, timing, memory, starting orientation, ground truth, final misorientation
- `hp_sweep_trajectory_{example}.csv` — subsampled step records (linked via `run_id` foreign key)

Gradient trajectory: every `TRAJ_SUBSAMPLE=10` steps, records angular step size and misorientation from ground truth.
MC trajectory: every accepted global improvement and every restart event — records angular step size and current step size (shows adaptive zoom behavior).

### Key modifications to `orientation_search.py`

Added `trajectory: Optional[List[Dict]] = None` parameter to `MCOptimizer.optimize()`:
- On each global improvement: appends `{step, event_type="mc_accept", angular_step_deg, cur_step_rad}` (cur_step already halved after improvement)
- On each restart: appends `{step, event_type="mc_restart", angular_step_deg, cur_step_rad=angular_step}`
- Non-breaking: default `trajectory=None` preserves original behavior
- Added `_quat_misorientation_deg(q1, q2)` helper for quaternion geodesic distance

### Benchmark usage

```bash
cd icenine_py

# Smoke test
uv run python benchmarks/bench_hp_sweep.py --smoke-test --example threevoxels

# Full ThreeVoxels (gradient optimizers, ~30 min)
uv run python benchmarks/bench_hp_sweep.py --example threevoxels --optimizer gradient

# MC only (fewer configs but MC is slow per run)
uv run python benchmarks/bench_hp_sweep.py --example threevoxels --optimizer mc

# Full ManyGrains overnight run
nohup uv run python benchmarks/bench_hp_sweep.py --example manygrains > hp_sweep_manygrains.log 2>&1 &

# Check progress
wc -l benchmarks/hp_sweep_manygrains.csv
```

### Results — Completed 2026-04-03

Both benchmarks completed successfully after fixing a `requires_grad` crash (see below).

**ManyGrains (58,500 runs, ~13.4 CPU-hours total):**

| Optimizer | Best HP | 1° success | 2° success | 5° success | Time/run |
|-----------|---------|-----------|-----------|-----------|----------|
| riemannian_adam_geoopt | lr=1e-4, n=100, β₁=0.9 | 96% | 41% | 0% | 0.44s |
| riemannian_adam_manual | lr=1e-4, n=200, β₁=0.9 | 94% | 49% | 0% | 0.58s |
| riemannian_sgd_plain | lr=1e-4, n=100 | 88% | 52% | 0% | 0.30s |
| riemannian_sgd_momentum | lr=1e-5, n=200, m=0.5 | 94% | 51% | 0% | 0.59s |
| riemannian_sgld | lr=1e-4, n=100, T=0.01 | 92% | 50% | 0% | 0.29s |
| mc_optimizer | n=3500, restarts=2, step=0.5 | 92% | 40% | 6% | 3.19s |

**ThreeVoxels (1,755 runs):** Same qualitative pattern; all methods succeed at 1° (3/3), most succeed at 2° (2/3), none at 5°.

**Key findings:**
1. Riemannian Adam (lr=1e-4, n=100) achieves 96% success at 1° vs. 92% for MC — and is 7× faster (0.44s vs. 3.19s)
2. Hard LR cliff: lr ≥ 0.05 → 0% success for all Adam variants. Optimal range: lr ∈ [1e-4, 1e-3]
3. More steps don't help at optimal lr: n=100 → 96%, n=200 → 94%, n=500 → 88%
4. Only MC achieves non-zero 5° success (6%) — gradient methods are trapped in flat landscape beyond ~2°
5. SGD variants surprisingly competitive at 2° perturbation (52%) vs Adam (41–49%)

**Bug fixed during run:**
In all 5 gradient `run_one_*` functions, replaced `if step > 0:` guard on `info.cost.backward()` with `if step > 0 and info.cost.requires_grad:`. When lr is very large (e.g. 0.1), R diverges → cost evaluates to a constant with `requires_grad=False` → `.backward()` throws `RuntimeError`. The guard skips the update step gracefully and continues recording.

**Output files:**
- `benchmarks/hp_sweep_manygrains.csv` — 58,500 runs (15 MB)
- `benchmarks/hp_sweep_threevoxels.csv` — 1,755 runs
- `benchmarks/hp_sweep_trajectory_threevoxels.csv` — step-by-step trajectories (3 MB)
- `benchmarks/hp_sweep_trajectory_manygrains.csv` — trajectories (105 MB, excluded from git via .gitignore)
- `benchmarks/hp_sweep_*.png` — LR sensitivity, n_steps sensitivity, optimizer comparison, trajectory plots

## Hybrid Riemannian Adam + MC-Restart Optimizer (2026-04-03)

**Branch**: `feature/differentiable-cost`
**Motivation**: HP sweep showed Riemannian Adam (lr=1e-4, n=100) achieves 96% success at 1° perturbation in 0.44s vs. 92% for MC in 3.19s. Replace the `FindOptimal` phase of `AdaptiveVoxelReconstructor` with a gradient-based optimizer that still falls back to random restarts when stuck.

### Architecture

**`RiemannianAdamOptimizer`** (new class in `orientation_search.py`):
- Uses `DifferentiableCostFunction` for gradient signal (geoopt Stiefel manifold parameter)
- Uses `VoxelCostFunction` for convergence decisions (hard binary overlap)
- MC-style restarts: random perturbation of best quaternion when stuck (same as `MCOptimizer`)
- SVD re-orthogonalization after each Adam loop (guards Stiefel float drift)
- Same `requires_grad` guard on `.backward()` as bench_hp_sweep.py bug fix

**Eval counts per restart** (n_steps=100, max_restarts=2):
- Hybrid: ~4 hard evals + 303 differentiable evals
- MC equivalent: ~3500 hard evals

### Files modified

| File | Change |
|------|--------|
| `icenine/orientation_search.py` | Added `_GEOOPT_AVAILABLE` guard, `_require_geoopt()`, `_make_stiefel_param()` helpers; added `use_hybrid_optimizer`, `adam_n_steps`, `adam_lr`, `adam_scale` to `SearchParameters`; added `RiemannianAdamOptimizer` class after `MCOptimizer` |
| `icenine/reconstructor.py` | Added `diff_cost_fn: Optional[object] = None` to `ReconstructionSetup`; added `RiemannianAdamOptimizer` to imports; added `build_diff_cost_fn()` module-level helper; updated `AdaptiveVoxelReconstructor.reconstruct_voxel()` FindOptimal block with hybrid/MC dispatch |
| `tests/test_orientation_search.py` | Added 4 tests: `test_requires_geoopt_error`, `test_returns_search_candidate`, `test_zero_restarts_single_hard_eval_after_adam`, `test_search_params_hybrid_defaults` |
| `benchmarks/bench_hybrid_optimizer.py` | New benchmark script: head-to-head Hybrid Adam vs MC, 5 perturbations × 20 voxels, with time breakdown (t_adam_sec + t_hard_eval_sec), 3 output plots |

### Enabling hybrid optimizer in reconstruction

```python
from icenine.reconstructor import setup_reconstruction, build_diff_cost_fn

setup = setup_reconstruction(config)
setup.diff_cost_fn = build_diff_cost_fn(setup)
setup.search_params.use_hybrid_optimizer = True
setup.search_params.adam_n_steps = 100
setup.search_params.adam_lr = 1e-4

reconstructor = AdaptiveVoxelReconstructor(setup)
```

### Benchmark usage

```bash
cd icenine_py
uv sync --extra riemannian

# Smoke test (3 voxels, 2 perturbations)
uv run python benchmarks/bench_hybrid_optimizer.py --smoke-test --example threevoxels

# Full ThreeVoxels (20 voxels, 5 perturbations)
uv run python benchmarks/bench_hybrid_optimizer.py --example threevoxels
```

### Outputs

- `benchmarks/bench_hybrid_{example}.csv` — one row per (voxel, perturbation, optimizer)
- `benchmarks/bench_hybrid_success_rate_{example}.png` — success rate vs perturbation
- `benchmarks/bench_hybrid_wall_time_{example}.png` — wall time with stacked Adam/hard-eval breakdown
- `benchmarks/bench_hybrid_scatter_{example}.png` — per-run MC vs hybrid misorientation scatter

## Toy Orientation NN — Theory Phase (2026-09-28)

**Branch**: `feature/nn-orientation-toy`
**Scope of the merged PR**: documentation only. The v0 prototype code stays local
and uncommitted for a follow-up feature branch (see "Prototype status" below).

**Motivation**: explore replacing the iterative orientation search
(`MCOptimizer`, `RiemannianAdamOptimizer`) with a neural network that predicts a
voxel's orientation offset from its detector data in one forward pass, starting
with local refinement around a known nominal orientation on simulated
Example2.ThreeVoxels data. A v0 smoke test exposed that the problem needed to be
formulated properly first; this phase did that.

### Documents added (`icenine_py/docs/`, render with pandoc, no `-N`)

| File | Content |
|------|---------|
| `omega_peak_width_derivation.md` | Rocking width of a peak in a rotation scan, Δω ≈ α/\|sin η\| (the Lorentz factor's 1/\|sin η\| part), and worst-case sensitivity of ω* to orientation, exactly 1/\|sin η\| rad/rad. Verified against `get_scattering_omegas_torch` on 5,100 (reflection, branch) pairs. |
| `nn_inverse_problem_formulation.md` | Why the frame/pixel-integrated, thresholded forward model has no inverse; what a supervised network learns instead (posterior mean with squared loss, covariance with Gaussian NLL); geometry of frame and pixel boundaries in orientation space; closed-form spot-motion Jacobian Γ_p (ring + parallax terms); angular-resolution estimates; references. Numerical checks throughout. |

### Key results

- **The target is a posterior, not an inverse.** The thresholded measurement map
  is piecewise constant in the orientation offset δ, so a region of orientations
  gives identical data. Squared loss converges to E[δ|D]; Gaussian NLL adds the
  covariance, which must be a full 3×3 matrix (the consistent set is
  anisotropic). The training perturbation distribution is the prior. The MMSE is
  a floor for every method, including MC and Riemannian Adam.
- **Frame boundaries** are level sets of ω*_p(δ); to first order parallel planes
  Δω_f·|sin η_p| apart. Holds within 0.5° for the 96% of Example2 peaks with
  |sin η| ≥ 0.3; strongly curved for near-axis peaks, which carry the most ω
  information. Rotation about the rotation axis shifts every ω* by exactly the
  rotation angle (all orders).
- **Spot motion** = sliding along the instantaneous Debye–Scherrer ring (cone
  with apex at the voxel's current lab position) + parallax (the voxel moves with
  the stage when ω* shifts). The ring accumulated over a scan is not a single
  conic. Closed-form Jacobian matches simulator finite differences to 1.5×10⁻⁴
  (median) for all of voxel 0's recorded peaks (364 at r_⊥ = 12 µm) at r_⊥ = 12–500 µm.
- **Resolution** (independent-quantisation model, idealised):
  σ_z ≈ 1/√(12P(1/Δω_f² + r_⊥²/(2a²))), σ_⊥ ≈ (a/d)/(κ√(6P)). About the
  rotation axis: frame-limited near the axis, parallax-limited beyond
  r_⊥ ≈ √2·a/Δω_f (≈120 µm for 1.48 µm pixels and 1° frames). Example2 voxel 0:
  σ_z ≈ 0.016°, σ_⊥ ≈ 0.0012° for the 364 recorded peaks (real HEDM is ~0.1°). Anisotropy ≈ 12 near the
  axis, ≈ 3–4 at 0.5 mm.

### Gotchas discovered

- **Peak count**: Example2's config has `MaxQ 16`. 790 (reflection, branch) pairs
  hit voxel 0's detector *plane*, but only 364 land on the 2048×2048 pixel grid
  and are ever recorded (corrected 2026-09-29; the first version of this section
  said 790). The printed "Max Q calculated: 24.79" is the geometric limit before
  the config cap. Reconstruction typically uses Q_max = 8.
- **Structure list index**: Example2's voxels have `phase = 1`;
  `structure_list[0]` is empty.
- **v0 evaluation was invalid**: the HP-sweep baseline is indexed by starting
  perturbation, v0 reported final-error thresholds, and there was no
  predict-nominal baseline.
- **v0 cannot see rotation about the rotation axis** for near-axis voxels: its
  windows drop the frame index, and pixels see that rotation only through
  parallax.
- **Only one detector used**: the prototype's ROI definition assigns each peak
  to the first detector it hits (all 364 of voxel 0's peaks are on detector 0);
  the experiment (and the C++ model) use both.
- **Simulator shortcuts**: noise-free intensities vary smoothly with δ and leak
  information real, thresholded data don't carry.
- **Sample translation** is zero and the base sample rotation is the identity
  in Example2, so solver frame = lab frame.

### Prototype status (local, uncommitted)

`icenine/orientation_nn.py` (ROI definition, single-frame windowed renderer,
perturbation sampler, dataset, loss/metric), `icenine/toy_orientation_model.py`
(flatten → FC 512/256/128 → quaternion), `scripts/generate_toy_orientation_dataset.py`,
`scripts/train_toy_orientation_nn.py`, `tests/test_orientation_nn.py` (14 tests,
passing), `.gitignore` entry `scripts/*.pt`, and the numerical-check scripts
behind `nn_inverse_problem_formulation.md` (`scripts/checks/`). The script behind
the derivation note's §7.1 check was lost with a temporary directory; its method
is described in that section and is easy to reproduce.

### Remaining plan

**Stage 0 — fix the evaluation.** Output a rotation offset δ plus a Cholesky
covariance, trained with Gaussian NLL (β-NLL if unstable). Report the ẑ and
perpendicular error components separately against the closed-form estimates.
Add a predict-nominal baseline and an exact Bayes baseline (sample the prior,
keep samples whose exact frame/pixel assignments match; gives the MMSE floor).
Report by starting-perturbation magnitude like `bench_hp_sweep.py`. Threshold or
noise the inputs. Re-run the v0 smoke test under this evaluation.

**Stage 1 — forward model and inputs** (in `orientation_nn.py`; leave the
C++-validated `ForwardSimulation` and cost functions unchanged). Record each peak
on every detector it hits. Use the reconstruction Q_max. Add per-peak metadata
(η, θ, nominal frame, ring tangent). Model the rocking width α/|sin η| as a box
split across frames (α = 0 must reproduce the single-frame renderer). Size each
peak's window from its spot-motion Jacobian and the frame drift β_max/|sin η|.
Test voxels at several distances from the axis. Keep peak sets and storage small
(Q_max, stratified subsets, sparse storage).

**Stage 2 — baselines without learning.** Predict-nominal; exact Bayes; Gauss–
Newton on (frame index, spot position) with the analytic Jacobians; HP-sweep MC
and Riemannian Adam at matching starting perturbations.

**Stage 3 — pixel network.** Shared per-peak encoder over frames × window plus
per-peak context (|sin η|, θ, detector, ring direction, position), masked pooling,
MLP head → δ + Cholesky factor. Compare per axis with Stage 2 and the floor.

**Stage 4 — realism.** Noise, intensity variation, detector point-spread and
partial pixel coverage, α from real data, multiple voxels and orientations, then
real data.

**Open decisions.**
- D1: α for simulation (sweep {0, 0.01°, 0.03°, 0.1°} until measured; config
  `BeamEnergyWidth 0.5` has unconfirmed units).
- D2: training prior = typical residual error after the coarse Sukharev search.
- D3: near-axis peaks vs window size (drop, per-peak ω-window, or small β_max).
- D4: a real dataset with known orientations to measure frame spread vs
  1/|sin η| and estimate α.
- D5: use the reconstruction Q_max (8 Å⁻¹) for the toy (recommended).

## Toy Orientation NN — Stage 0: Fix the Evaluation (2026-09-29)

**Branch**: `feature/nn-orientation-stage0`
**Goal**: replace the invalid v0 evaluation with one that can tell whether a network
learned anything: offset + covariance outputs, per-axis errors by perturbation
magnitude, and predict-nominal and exact-Bayes reference rows. Re-run v0 under it.

### What was added

| File | Content |
|------|---------|
| `icenine/orientation_eval.py` | `BatchedObserver`: float64, batched re-implementation of `_simulate_peaks`' per-peak maths. Returns for B candidate offsets: which ROI peaks are recorded, their frame, and their spot vertices. Agrees with the simulator on presence and frame index for every peak checked and on centroids to 5e-4 px (the test asserts 1e-3), ~1 ms/candidate. `ExactBayes`: posterior of the offset given the thresholded data, by importance sampling restricted to the set reproducing the observed frames and lit-pixel sets (`lit_pixel_set` reproduces the rasteriser's truncate, clip, round and fill exactly). `error_summary`: RMS about the stage axis (z) and perpendicular (x, y), median misorientation, success below 0.5° (the `bench_hp_sweep.py` criterion) and 0.1°. Rotation-vector helpers and priors. |
| `icenine/orientation_nn.py` | `cholesky_from_raw`, `gaussian_nll_loss` (with optional β-NLL), `spot_overlaps_grid`. `OrientationDataset` now returns offsets in degrees; windows are stored thresholded (uint8). |
| `icenine/toy_orientation_model.py` | `ToyOffsetNet`: same trunk as v0, head outputs offset (3) + Cholesky factor (6). |
| `scripts/generate_toy_orientation_dataset.py` | Train set from a ball prior (radius 2.5°); test set at fixed magnitudes 0.25/0.5/1/2° in random directions; thresholded windows. |
| `scripts/exact_bayes_baseline.py` | Exact posterior for every test case (`--frames-only` for the frames-only floor). |
| `scripts/train_toy_orientation_nn.py` | `--head offset` (NLL) or `--head quat` (v0), best-validation checkpoint, per-magnitude table vs predict-nominal and exact Bayes, saved predictions. |
| `tests/test_orientation_eval.py` | 22 tests: rotation conventions, observer vs simulator, exact ẑ shift of every ω*, the pixel-grid rule, `lit_pixel_set` vs the real rasteriser on 400 random triangles straddling the grid edge, NLL vs `MultivariateNormal`, Bayes posterior tightness, per-axis metrics. |
| `benchmarks/toy_orientation_stage0/` | Bayes posteriors, predictions and result tables behind the numbers below. |

### Correction found during Stage 0: most "ROI peaks" were never recorded

`define_roi_set`, the renderer and the first version of the observer counted a peak
whenever its ray hit the detector *plane*. The simulator's rasteriser clips spots
to the 2048×2048 pixel grid, so a spot outside it produces no data. For voxel 0,
**426 of the 790 ROI peaks (54%) are entirely off the grid** (spot rows reach 3,094)
and only **364** are recorded. Consequences of the first version, all corrected:

- v0's windows for those 426 peaks contained spots a real detector would never
  record, so its inputs were partly fictitious.
- The exact-Bayes floors used their frames and exact vertices as data and were
  too tight by 2–3×.
- Every "790 peaks" number in the theory docs (resolution checks, the 5,100 vs
  2,352 count in the ω-width note; the correct total over the three voxels is 1,083)
  used the wrong set.

The fix is one rule shared by the renderer and the observer (`spot_overlaps_grid`:
the bounding box of the rasteriser-truncated vertices meets the grid) plus the
rasteriser's own clip in the pixel-set comparison. It was caught by directly checking
the statement that no spots touch the detector edge, which failed. After
the fix the exact-Bayes sampler is calibrated to within sampling error on all three
axes (mean squared error of the posterior mean over mean posterior variance
0.99, 1.02, 1.02; it was 0.94, 0.98, 0.78 before).

### Results (voxel 0, 364 recorded peaks, 1,350 training samples, 30 epochs, one seed, 30 test cases per magnitude)

RMS error in degrees, z = about the stage axis, ⊥ = perpendicular; success = misorientation < 0.5°.

| \|δ\| | predict-nominal z / ⊥ (success) | v0 quaternion head z / ⊥ (success) | offset head z / ⊥ (success) | exact Bayes z / ⊥ |
|---|---|---|---|---|
| 0.25° | 0.149 / 0.142 (100%) | 0.386 / 0.129 (87%) | 0.518 / 0.107 (60%) | 1.4e-3 / 6e-5 |
| 0.50° | 0.308 / 0.279 (53%) | 0.454 / 0.123 (77%) | 0.576 / 0.100 (50%) | 2.0e-3 / 7e-5 |
| 1.00° | 0.636 / 0.546 (0%) | 0.410 / 0.145 (73%) | 0.501 / 0.111 (63%) | 1.9e-3 / 6e-5 |
| 2.00° | 1.094 / 1.184 (0%) | 0.562 / 0.221 (30%) | 0.668 / 0.266 (30%) | 1.7e-3 / 8e-5 |

- **Neither network beats the trivial baseline about z until |δ| ≥ 1°** (RMS z
  error is worse than predicting nominal at 0.25° and 0.5°), and at 0.25° both have
  a lower success rate than predicting nominal (100%). Perpendicular to the axis
  both are better than nominal at every magnitude (1.1–1.3× at 0.25°, 2–5× at
  0.5–2°). Both are far above the noise-free floor: about 2–5×10² times about z and
  1.5–3×10³ times perpendicular.
- **Head comparison is not conclusive**: the quaternion head is ahead in this run,
  but with one seed and n = 30 per bin the difference is within noise. The offset
  head's covariance is roughly calibrated (mean squared Mahalanobis distance
  2.0–4.1 against 3.0; predicted σ 0.19–0.50°).
- **Overfitting**: about 191M parameters (364×32×32×512 in the first layer) on 1,350 samples. The offset head's
  validation NLL was best at epoch 10 and noisy afterwards while training NLL kept
  falling; the best-validation checkpoint is what is reported.
- Not comparable to the HP-sweep success rates (Riemannian Adam 96% at 1°, MC 92%):
  different voxels and data, and the network is trained for this one voxel and
  nominal orientation.

### Exact noise-free floors (details in `docs/nn_inverse_problem_formulation.md` §3.6.7)

Median posterior σ (degrees), frames + lit pixels: ⊥ 5.5e-5, z 1.2e-3; frames only:
⊥ 3e-3, z 2.1e-3. These are 22× (⊥) and 13× (z) below the independent-quantisation
estimates of §3.6.5 (√P = 19), close to the 1/P scaling of the set-membership
picture, and independent of the offset size over 0.25–2°. The frames-only z cell
width √12·σ_z = 7.2e-3° is within 1.3× of the Vernier prediction 2Δω_f/(P+1).

### Findings and gotchas

- **The 32×32 windows are too small for the prior.** On average 34% of the 364
  windows are empty in the training set, and the empty fraction correlates 0.98
  with the perpendicular offset (0.36 per degree) but not with the z offset. Spots
  drift ~16 px/deg along their ring against a 16 px half-width. Most of what v0
  learns perpendicular to the axis is "which spots have left their window".
- **v0's information about z comes from the scan limits, not from frames.** Peaks
  near the +90° end of the ω range drop out in a way that tracks the z offset
  (correlation −0.83, about 1.2 lost peaks per degree; only 3 peaks lie within 3°
  of the −90° end); interior peaks are lost only through the perpendicular offset
  (correlation with z 0.01). That gives well under 1° resolution from ~11 peaks and
  would not exist for a 360° scan.
- **Noise-free floors are not a usable target.** 1.2e-3° (4 arcseconds) is far
  below anything a trained network, or a real experiment, will approach. They are a
  check on the physics and the sampler; meaningful comparisons are predict-nominal
  now and the Gauss–Newton baseline and the HP-sweep optimizers next.
- **`torch.linalg.norm` on length-3 vectors was ~10× slower than an explicit sum of
  squares** for large batches; the observer uses the latter (3.7× overall).
- Running the scripts changes the working directory (`setup_example` chdirs to the
  example folder), so every path argument must be resolved before the call.
- **Code review (Opus): APPROVE** with no critical issues; its warnings were
  addressed (hard-coded paths removed from `scripts/checks/`, a failed Bayes case now
  returns NaN instead of the truth, Bayes/test alignment asserted, claims softened,
  new files Black-formatted). Not done: mypy strict annotations for the new modules
  (19 errors, the rest of the repo is not clean either), an end-to-end test against
  `ForwardSimulation` output, and restoring the working directory after
  `setup_example`.
- The old v0 dataset/checkpoint files (2.6 GB, gitignored) were deleted: the
  dataset format changed.

### Approximations in the exact Bayes baseline

Only the ROI peaks (recorded at the nominal orientation) count as data; peaks that
would newly appear at the true offset are ignored (slightly wider posterior).
Overlaps between different peaks' lit pixels are ignored. A spot that overlaps the
grid's bounding box but clips to nothing is treated as recorded with an empty pixel
set (0.07% of present spots over 300 prior draws; none at nominal). Peaks are
required to hit only their home detector, whereas the serial simulator drops a peak
on every detector if any vertex misses any detector; for Example2 nothing changes,
but Stage 1 (both detectors) must apply the all-detectors rule. The simulator
works in float32 and the observer in float64 (a vertex within ~1e-4 px of an integer
can truncate differently). The network's windows are not exactly the data Bayes
conditions on: their origin is fractional (`nominal - 16`) and they clip to
themselves rather than to the detector grid, so the Bayes floor is a generous
reference, not the floor for the network's exact inputs.

### Plan status

Stage 0 is done apart from β-NLL training (implemented and unit-tested, not run).
Stage 1 (forward model with frames and both detectors, windows sized from the
spot-motion Jacobian, Q_max) is next. The window-size finding above, and the
open decisions D1–D5 in the previous section, apply.

## Toy Orientation NN — Stage 1: Network Inputs (in progress, started 2026-09-29)

**Branch**: `feature/nn-orientation-stage1`
**Goal**: fix what the network is given as input, which Stage 0 identified as the
main limitation (no frame index, spots drifting out of 32×32 windows, one
detector).

### Decisions (agreed 2026-09-29)

- **D1, rocking width α**: start at α = 0 (today's single-frame physics); add
  α ∈ {0.01°, 0.03°, 0.1°} once the frame input works.
- **D2, training prior**: measure the typical residual error after the coarse
  search by running the existing reconstructor on Example2; keep the 2.5° ball
  until then.
- **D3, near-axis peaks**: ±4 frames around the nominal frame, and drop peaks with
  |sin η| < 0.3 for now; revisit with per-peak window lengths later.
- **D4, real data for α**: open; not needed until α matters.
- **D5, Q_max**: use 8 Å⁻¹, the reconstruction value.

### D2 measured: what the coarse search hands to FindOptimal

`scripts/checks/measure_coarse_residual.py` runs `AdaptiveVoxelReconstructor` with
`ReconstructQ8.config` (Q_max = 8; 5° grid radius shrinking over 4 levels) on the three
Example2 voxels against `ScatteringData_Python`, 5 seeds each, and records the
orientations passed to FindOptimal (`benchmarks/toy_orientation_stage1/coarse_residual.npz`).

- **Voxel 2**: the coarse search fails on every seed (a known failure shared with C++);
  the best hand-off is ~60° off.
- **Voxels 0 and 1**: a good candidate reaches FindOptimal in 9 of 10 runs, with error
  0.07–1.09° (median about 0.35°), almost entirely rotation about the stage axis
  (0.2–0.7°; about 0.04° perpendicular; one case 0.98° about x). This is the anisotropy
  predicted in `docs/nn_inverse_problem_formulation.md` §3.6.
- The top-ranked candidate is often a wrong ~54° solution (4 of 10 runs) even when a
  good one is passed on, and FindOptimal then sometimes keeps the wrong one (voxel 1,
  seeds 0 and 1: final errors 54° and 42° from good hand-offs). A network should be
  evaluated per candidate, as FindOptimal is.
- **Decision**: training prior = uniform ball of radius 1°, which covers the good
  hand-offs; test magnitudes 0.1, 0.25, 0.5, 1°. With this prior, 32×32 windows and
  ±4 frames keep every spot that is present (max drift 16 px, 3 frames).

### Results: items 1–3 (2026-09-29)

Setup: voxel 0, Q_max = 8, both detectors, near-axis spots dropped (|sin η| < 0.3):
113 spots. Frame-coded exact windows, 32×32, ±4 frames. Prior: 1° ball. 10,000
samples (9,000 train, 1,000 validation), 30 epochs, one seed; 30 test cases per
magnitude. RMS error in degrees (z = about the stage axis, ⊥ = perpendicular);
success = misorientation < 0.1°.

| \|δ\| | predict-nominal z / ⊥ | offset head z / ⊥ (success) | offset head, no frame channel z / ⊥ (success) | quaternion head z / ⊥ (success) | exact Bayes z / ⊥ |
|---|---|---|---|---|---|
| 0.10° | 0.055 / 0.059 | 0.032 / 0.007 (100%) | 0.164 / 0.007 (40%) | 0.034 / 0.023 (100%) | 0.009 / 0.0004 |
| 0.25° | 0.148 / 0.143 | 0.027 / 0.008 (100%) | 0.192 / 0.009 (50%) | 0.047 / 0.026 (97%) | 0.008 / 0.0003 |
| 0.50° | 0.249 / 0.306 | 0.037 / 0.016 (100%) | 0.141 / 0.009 (37%) | 0.049 / 0.030 (90%) | 0.007 / 0.0003 |
| 1.00° | 0.584 / 0.574 | 0.058 / 0.019 (90%) | 0.190 / 0.026 (27%) | 0.083 / 0.038 (50%) | 0.007 / 0.0004 |

- **The new inputs work.** The offset head's median misorientation is 0.02–0.05°
  (Stage 0: 0.35–0.62°), 100% of cases are within 0.5° at every magnitude, and it beats
  predict-nominal about z at every magnitude. It is within about 4–8× of the noise-free
  floor about z and 20–50× perpendicular.
- **The frame channel is what fixes z.** Hiding it (same data, same model) leaves the
  perpendicular error unchanged but raises the z error to 0.14–0.19°, worse than
  predict-nominal at 0.1° and 0.25°. Pixels constrain perpendicular rotation, frames
  constrain rotation about the axis, as the docs predict.
- **Offset head vs quaternion head**: the offset head is better perpendicular (2–3×) and
  its covariance is roughly calibrated (mean squared Mahalanobis distance 1.4–4.1
  against 3.0).
- **Confounds vs Stage 0**: several things changed at once (prior 1° vs 2.5°,
  Q_max 8, both detectors, exact windows, 10,000 vs 1,500 samples). Only the frame
  ablation isolates a single factor. One seed.
- **Not yet comparable to FindOptimal**: the network is trained for one voxel and one
  nominal orientation on noise-free data. On the same Example2 data FindOptimal's final
  errors were 0.04–0.2° when it succeeded.

### Planned work, in order

1. Peak definition: record each peak on every detector whose pixel grid its spot
   overlaps, applying the simulator's rule that a peak is dropped on all detectors
   if any spot vertex misses any detector plane; configurable Q_max.
2. Frame index in the input: windows of a few frames × pixels around each peak's
   nominal frame.
3. Windows sized from each peak's spot-motion Jacobian and the perturbation range.
4. Frame-spread renderer (α/|sin η| profile); α = 0 must reproduce the current
   output exactly.
5. Re-run the Stage 0 evaluation on the new inputs and compare.

## Toy Orientation NN — Stage 2: Baselines Without Learning (2026-09-29)

**Branch**: `feature/nn-orientation-stage1`

### What was added

- `icenine/orientation_baselines.py`, `scripts/gauss_newton_baseline.py`:
  `CentroidGaussNewton` fits the orientation offset to each spot's recorded frame
  (central rotation angle, variance Δω²/12) and lit-pixel centroid (variance 1/12 px²)
  by damped Gauss–Newton on the batched observer. It uses exactly the information in
  the frame-coded windows, and the ROI peak identities are known.
- `scripts/optimizer_baselines.py`: MC (`MCOptimizer`, 3500 steps, 2 restarts, step 0.5,
  as in the sweep's headline config) and geoopt Riemannian Adam (lr 1e-4, 100 steps,
  scale 2, ω window 1) at the test set's magnitudes and directions, using the existing
  Python-simulated Example2 images (generated at the voxel's ground-truth orientation)
  with the search started from a displaced orientation, as in `bench_hp_sweep.py`.
  Reflections limited to |g| ≤ 8 Å⁻¹ like the networks.

### Results (voxel 0, 120 test cases, 30 per magnitude; RMS error in degrees, success = misorientation < 0.1°)

| \|δ\| | Gauss–Newton z / ⊥ (success) | Riemannian Adam z / ⊥ (success) | MC z / ⊥ (success) | exact Bayes z / ⊥ |
|---|---|---|---|---|
| 0.10° | 0.045 / 0.003 (100%) | 0.071 / 0.018 (87%) | 0.058 / 0.022 (87%) | 0.009 / 0.0004 |
| 0.25° | 0.033 / 0.004 (100%) | 0.143 / 0.018 (30%) | 0.150 / 0.072 (37%) | 0.008 / 0.0003 |
| 0.50° | 0.044 / 0.003 (97%) | 0.250 / 0.018 (20%) | 0.280 / 0.248 (3%) | 0.007 / 0.0003 |
| 1.00° | 0.045 / 0.003 (100%) | 0.584 / 0.027 (13%) | 0.424 / 0.496 (0%) | 0.007 / 0.0004 |

Cost: Gauss–Newton 7 ms per case (χ²/dof 0.91, so the quantisation noise model fits),
Adam 0.23 s, MC 1.8 s.

- **Gauss–Newton is a strong baseline.** With no learning it reaches ⊥ 0.003° and z
  0.03–0.045° at every magnitude. Its z error is 1.5× its own predicted σ (0.027°):
  frame quantisation errors are not independent across spots, so the independent-error
  model of `docs/nn_inverse_problem_formulation.md` §3.6.4 is optimistic about z.
- **Adam and MC leave the stage-axis component (largely) uncorrected on these voxels.**
  Adam's z error is about the starting z offset (0.58° at 1°, the RMS of a 1° offset's
  z component), while its perpendicular error is 0.02°. Reproduced with
  `bench_hp_sweep.py`'s own perturbation generator and runner (30 random axes, voxel 0):
  Adam success (< 0.5°) 60% at 1° for lr 1e-4 and 60–67% for lr 1e-3 to 1e-2.

### Reconciliation with the earlier, more extensive study

The 96% (Adam) and 92% (MC) figures for 1° in the HP-sweep section above are from the
**ManyGrains** study (58,500 runs, 100 voxels), not from Example2.ThreeVoxels, and the
two samples differ in exactly the way `docs/nn_inverse_problem_formulation.md` §3.6
predicts matters:

| | ThreeVoxels | ManyGrains sweep (100 selected voxels) |
|---|---|---|
| distance from the rotation axis r⊥ | 12 µm (all three voxels) | median 370 µm, 95% beyond 120 µm (min 75 µm) |
| voxel side | 0.75–1.5 µm (≈ 1 pixel) | 9.4 µm (≈ 6 pixels) |

The docs predict that pixels see rotation about the stage axis only through parallax,
which becomes the main source of that information beyond r⊥ ≈ √2·a/Δω_f ≈ 120 µm. So a
pixel-overlap cost should determine the stage-axis component for the ManyGrains voxels
and barely at all for Example2's, which is what the two studies show. Checks against the
records:

- The earlier ThreeVoxels benchmark (lr 0.01) recorded Riemannian Adam at 1° as 8/18
  successes (44%), consistent with 43% here.
- `hp_sweep_threevoxels.csv` (1,755 runs of the extensive sweep, same three voxels), for
  the ManyGrains headline configs, at 1°: Adam final misorientations 0.57 / 0.63 / 0.22°
  and MC 0.85 / 0.79 / 0.38° for voxels 0 / 1 / 2 (success 1 of 3 each). At 2° and 5°
  the headline Adam config succeeds on none of the three.
- The "all methods succeed at 1° (3/3)" statement for ThreeVoxels above is a best-of-sweep
  statement (the best of 42 or 36 HP configs per voxel). For fixed configs, 9 of 42 Adam
  configs and 0 of 36 MC configs succeed on at least two of the three voxels at 1°.

**Not established**: the per-run ManyGrains CSV (`hp_sweep_manygrains.csv`) was never
committed and is not on this machine (only the log and the trajectory CSV are), and
ManyGrains' `ScatteringData_Python` is absent, so the dependence on r⊥ (versus voxel size,
overlapping neighbours or the different reflection set) has not been measured. The two
samples differ in more than one way. The conclusion above applies to Example2's near-axis
voxels; it is not evidence against the ManyGrains results.

### Caveats on comparing Gauss–Newton with Adam and MC

- Gauss–Newton is given the peak identities (which predicted spot is which, from the
  nominal orientation) and the frame index, and starts inside the basin. Adam and MC must
  associate predicted spots with the image through pixel overlap, so they solve a harder
  problem; the comparison shows what the information in the frame-coded windows supports,
  not that the optimizers are worse algorithms.
- The optimizers were run with the data at the ground truth and the start displaced; the
  networks and Gauss–Newton see data displaced from a known nominal. These are the same
  local problem to first order, not identical.

### Open

- Test the near-axis explanation directly: forward-simulate a few ManyGrains voxels (or
  move Example2's voxel to r⊥ ≈ 400 µm), and rerun Adam and MC at 1° by r⊥; recover or
  regenerate `hp_sweep_manygrains.csv`.

## Toy Orientation NN — Stage 3: Iterating on the Set Network (2026-09-29)

Branch `feature/nn-orientation-stage1`. Plan: improve `PeakSetNet` (Step 1), repeat on a
voxel far from the rotation axis (Step 2), train on ~30 voxels (Step 3). Training runs on the
Apple GPU (`--device mps`; CPU vs MPS one-epoch losses with `--seed 0` agree: -0.87531 vs
-0.87529). Evaluation stays float64 on the CPU.

### Step 1: four rounds on voxel 0 (r⊥ = 12 µm), same data and test set as Stage 2

Median misorientation angle (deg) at |δ| = 0.1 / 0.25 / 0.5 / 1.0 (30 cases each), from the
saved `benchmarks/toy_orientation_stage2/res_*.json` (one clean run each). "set" is the Stage 2
run (lr 1e-3, 30 epochs, no schedule).

| run | median angle | z RMS | ⊥ RMS | net σ_z | mean Mahalanobis² |
|---|---|---|---|---|---|
| set (Stage 2) | 0.045/0.146/0.190/0.581 | 0.054/0.152/0.242/0.585 | 0.015/0.016/0.018/0.019 | 0.51/0.50/0.47/0.42 | 0.8/0.9/1.4/3.7 |
| r1 optimisation | 0.051/0.141/0.192/0.538 | 0.054/0.149/0.248/0.585 | 0.017/0.019/0.017/0.024 | 0.46/0.46/0.45/0.44 | 2.1/2.5/2.6/5.3 |
| r2 r1 with `--no-frame` | 0.052/0.142/0.191/0.547 | 0.054/0.148/0.247/0.585 | 0.019/0.016/0.019/0.026 | 0.45/0.45/0.45/0.45 | 2.1/2.0/2.4/4.6 |
| r3 + measurement features | 0.044/0.140/0.194/0.529 | 0.052/0.149/0.241/0.576 | 0.008/0.008/0.008/0.010 | 0.58/0.56/0.49/0.36 | 2.5/2.3/2.2/4.3 |
| r4 r3 with mean + sum pooling | 0.065/0.165/0.283/0.451 | 0.041/0.033/0.045/0.057 | 0.042/0.117/0.207/0.356 | 0.04/0.04/0.04/0.04 | 1.9/2.1/3.1/4.0 |
| fc (Stage 1) | 0.023/0.024/0.036/0.049 | 0.032/0.027/0.037/0.058 | 0.007/0.008/0.016/0.019 | | |
| Gauss–Newton | 0.040/0.024/0.035/0.036 | 0.045/0.032/0.044/0.045 | 0.003/0.003/0.003/0.003 | | |
| exact Bayes | 0.007/0.004/0.003/0.004 | 0.009/0.008/0.007/0.007 | 0.000/0.000/0.000/0.000 | | |

Rounds (all lr/clip/cosine/60 epochs from r1 on: `--lr 3e-4 --clip 1.0 --cosine`):

1. **Optimisation.** Gradient clipping, cosine decay and a lower lr remove the erratic
   validation loss (monotone to epoch 60) but change nothing in accuracy: z is still
   predict-nominal.
2. **Frame ablation (`--no-frame`).** Identical to r1 in every column, so the conv-encoded set
   net never used the frame channel; its z output is the prior.
3. **Explicit per-peak measurement features** (`measurement_features` in
   `toy_orientation_model.py`: present flag, lit count, mean frame offset, lit-pixel centroid
   relative to the window centre; they equal `extract_measurements()`, tested). This fixes the
   perpendicular axes (⊥ RMS 0.017-0.024 → 0.008-0.010, better than fc and the ⊥ of Stage 1's
   fc at 0.5/1 deg) but z is still the prior (σ_z ≈ 0.5, it knows it does not know z).
4. **Sum pooling** (mean and sum instead of mean and max). This fixes z (RMS 0.04-0.06, on par
   with fc and Gauss–Newton) but loses the ⊥ accuracy, particularly y (σ_y 0.16-0.46).

**Findings.** The stage axis is a global quantity: δ_z ≈ -(mean over peaks of the frame
residual), so it needs a pooled sum of per-peak frame measurements, and max pooling (which
discards it) plus a mean of ReLU features was not found by optimisation. Conversely max
pooling is what gives the ⊥ components. Neither r3 nor r4 alone meets the target (median ≤ fc
and ≤ 1.5x Gauss–Newton at every bin): r3 is at the prior in z; r4 is at 0.28-0.45 deg in
⊥ at large δ. The rounds were capped at four, so the obvious combination (mean + max + sum
pooling, `--pool all`, added to the code after round 4 and covered by the padding test) was
not run on voxel 0; it is tried on the Step 2 data.

### Step 2: a voxel far from the rotation axis

Generator changes: `--example {threevoxels,manygrains}` (stored in the dataset meta as `example`,
with `r_perp_um` and `side_um`); `gauss_newton_baseline.py` and `exact_bayes_baseline.py` read the
example from the dataset meta instead of hard-coding ThreeVoxels.

**Voxel**: `Examples/Example2.ManyGrains/SimInput/rand_500grains_1mm_inFZ.mic`, voxel **77**
(first candidate near 400 µm): r⊥ = **398.7 µm**, side 9.38 µm (triangle), **118 ROI peaks**
(Q_max 8, both detectors, |sin η| ≥ 0.3). Dataset `scripts/toy_orientation_stage3_far_*.pt`
(10 000 train from the 1° ball, 4 x 30 test at 0.1/0.25/0.5/1.0°; 96.8-98.1 % of spots inside
their windows). Exact Bayes (7 min) and Gauss–Newton (1 s) ran on the test set; MC/Adam were not
run (ManyGrains has no detector images). The Bayes row is summarised from
`far_test_bayes.npz` with `scripts/summarize_bayes_npz.py` (the trainings started before it
finished, so their tables have no Bayes row). All networks: lr 3e-4, clip 1, cosine, batch 64,
`--device mps`; fc 30 epochs, set 60. Results in `benchmarks/toy_orientation_stage3/far_res_*.json`.

Median angle (deg) at |δ| = 0.1 / 0.25 / 0.5 / 1.0, z and ⊥ RMS, mean Mahalanobis²:

| method | median angle | z RMS | ⊥ RMS | net σ (x,y,z) at 0.25 | Mahalanobis² |
|---|---|---|---|---|---|
| predict-nominal | 0.100/0.250/0.500/1.000 | 0.055/0.148/0.249/0.584 | 0.059/0.143/0.306/0.574 | | |
| fc | 0.009/0.008/0.009/0.020 | 0.009/0.007/0.006/0.016 | 0.004/0.004/0.006/0.017 | 0.004, 0.005, 0.008 | 3.7/2.7/4.1/11.1 |
| set, mean+max pool | 0.040/0.113/0.157/0.450 | 0.068/0.143/0.212/0.489 | 0.013/0.024/0.039/0.096 | 0.098, 0.054, 0.460 | 2.1/2.1/2.3/6.3 |
| set, mean+sum pool | 0.078/0.156/0.294/0.431 | 0.021/0.022/0.029/0.044 | 0.069/0.128/0.212/0.340 | 0.039, 0.237, 0.028 | 1.6/2.7/3.3/4.1 |
| **set, mean+max+sum pool** | 0.015/0.013/0.015/0.018 | 0.015/0.014/0.017/0.019 | 0.010/0.010/0.009/0.016 | 0.009, 0.008, 0.015 | 3.3/3.7/3.4/5.6 |
| Gauss–Newton | 0.012/0.011/0.009/0.008 | 0.013/0.010/0.010/0.008 | 0.003/0.003/0.003/0.003 | | |
| exact Bayes | 0.001/0.001/0.001/0.001 | 0.001/0.001/0.001/0.001 | 0.000/0.000/0.000/0.000 | | |

Findings.

- Parallax helps every method that can use it, as expected from the near-axis analysis: on
  voxel 0 (r⊥ = 12 µm) Gauss–Newton had z RMS 0.045° and exact Bayes 0.007-0.009°; at
  r⊥ = 399 µm they are 0.008-0.013° and 0.001°. fc median error is 0.009-0.020° (voxel 0:
  0.023-0.049°).
- **The pooling combination that fixes both axes is mean + max + sum** (the combination not tried
  in Step 1, `--pool all`): 0.013-0.018° at every magnitude. The single-pool nets reproduce the
  Step 1 split: mean+max is at the prior in z even here (σ_z 0.46), mean+sum loses ⊥ (σ_y 0.24).
- The target "median ≤ fc and ≤ 1.5x Gauss–Newton at every bin" is not met by the set net at
  0.1-0.5° (fc 0.008-0.009; GN 0.009-0.012; set 0.013-0.015); at 1° it beats fc (0.018 vs 0.020).
  Calibration: Mahalanobis² 3.3-3.7 for δ ≤ 0.5° (calibrated is 3), 5.6 at 1°; fc 2.7-4.1, but
  11.1 at 1°.
- Exact Bayes is ~10x below all methods here: the information in the thresholded data is far
  from exhausted (on voxel 0 the gap was 4-6x).

### Step 3: 30 voxels, one network (held-out voxels test generalization)

**Generator** (`--n-voxels 30 --voxel-seed 0 --r-max-um 500`): 30 target radii evenly spaced over
r⊥ = 0-500 µm, one random ManyGrains voxel near each (distinct orientations, ≥ 40 peaks; ROI
sets of 105-130 peaks, Q_max 8, both detectors, |sin η| ≥ 0.3, ±4 frames, 1° prior).
Every 5th voxel in r⊥ order (offset 2: r⊥ = 32, 126, 215, 301, 372, 457 µm) is **held out**: it
has no training data. Train: 24 voxels x 500 prior offsets = 12 000 samples; test: every voxel
x 4 magnitudes (0.1/0.25/0.5/1.0°) x 10 offsets = 1 200 cases (960 in-distribution offsets on
training voxels, 240 on held-out voxels). Storage: windows `(N, M_max = 130, 32, 32)` uint8
(zero = absent), per-voxel `context (30, 130, 16)` gathered by `voxel_id`, `R_nom (30, 3, 3)`;
train file 1.6 GB, test 0.16 GB (`scripts/toy_orientation_stage3_multi_*.pt`, gitignored);
generation 3 min. `PeakSetNet` needed no change (accepts a `(B, M, D)` context); tests added for
the context gather and for invariance to zero-padded peaks. No exact Bayes here (too slow for 1 200
cases and 30 voxels; GN is the reference). `fc` is not applicable (fixed peak count).

**Network**: `PeakSetNet(pool="all")` (measurement features, mean+max+sum pooling), lr 3e-4, clip 1,
cosine, batch 64, `--device mps`. 60 epochs takes 7 min; a 150-epoch run is also reported
(best validation at epoch 136). Results from
`benchmarks/toy_orientation_stage3/multi_res_set_all{,_150ep}.json`.

Median angle (deg) at |δ| = 0.1 / 0.25 / 0.5 / 1.0:

| test set | method | median angle | z RMS | ⊥ RMS | Mahalanobis² |
|---|---|---|---|---|---|
| in-distribution (24 voxels, 240 per bin) | predict-nominal | 0.100/0.250/0.500/1.000 | 0.057/0.147/0.282/0.601 | 0.058/0.143/0.292/0.565 | |
| | set, 60 epochs | 0.030/0.030/0.032/0.038 | 0.039/0.041/0.046/0.049 | 0.009/0.009/0.010/0.014 | 4.0/3.5/3.2/4.1 |
| | set, 150 epochs | 0.029/0.026/0.029/0.032 | 0.038/0.039/0.041/0.042 | 0.007/0.006/0.007/0.011 | 4.1/3.6/3.5/5.7 |
| | Gauss–Newton | 0.012/0.015/0.012/0.014 | 0.020/0.024/0.022/0.022 | 0.003/0.003/0.003/0.004 | |
| held-out voxels (6, 60 per bin) | set, 60 epochs | 0.036/0.043/0.046/0.055 | 0.030/0.041/0.047/0.055 | 0.019/0.020/0.020/0.028 | 19.5/21.6/17.4/17.8 |
| | set, 150 epochs | 0.034/0.037/0.039/0.053 | 0.031/0.043/0.051/0.054 | 0.018/0.016/0.018/0.024 | 31.1/24.5/26.3/34.3 |
| | Gauss–Newton | 0.014/0.012/0.012/0.010 | 0.017/0.020/0.020/0.016 | 0.004/0.003/0.004/0.004 | |

**Error vs r⊥** (median over voxels of each voxel's median angle over all magnitudes, 60-epoch run;
6/7/5/6 training voxels and 2/1/2/1 held-out voxels per bin):

| r⊥ bin (µm) | set, train voxels | set, held-out voxels | Gauss–Newton (train voxels) |
|---|---|---|---|
| 0-130 | 0.029 | 0.040 | 0.022 |
| 130-260 | 0.031 | 0.043 | 0.013 |
| 260-390 | 0.033 | 0.045 | 0.011 |
| 390-510 | 0.036 | 0.056 | 0.008 |

Across all 30 voxels the correlation of the voxel median error with r⊥ is +0.37 for the network
and -0.63 for Gauss–Newton.

**Findings.**

- One network trained on 24 voxels works on all of them: 0.026-0.038° median in-distribution,
  and 0.034-0.055° on the six voxels it never saw (about 1.2-1.6x the in-distribution error, with
  the largest degradation at 1° and at large r⊥). The per-voxel median error over all magnitudes
  (offsets up to 1°) is 0.024-0.056° for every one of the 30 voxels.
- Relative to per-voxel Gauss–Newton it is 2.5-4x worse, and the gap widens with r⊥ because
  Gauss–Newton exploits the parallax (its error falls from 0.022° to 0.008° with r⊥) and the
  network does not: its z RMS stays at 0.04° and its voxel error rises slightly with r⊥. The
  single-voxel specialist for voxel 77 (Step 2) reached 0.013-0.018° there, so the multi-voxel net
  gives up ~2x at large r⊥ in exchange for generality.
- Mean Mahalanobis² is 3.2-4.1 in distribution (nearly calibrated) but 17-34 on held-out voxels:
  the predicted covariance is overconfident by ~2-3x in σ on unseen voxels, and 150 epochs makes
  it worse (24-34). Do not trust the predicted σ for new voxels without recalibration.
- 150 epochs (best validation at epoch 136) improves in-distribution by ~10 % and held-out barely;
  the network is limited by the 24 training voxels / 12 000 samples and by its parameter
  budget, not by optimisation time.

**Caveats.** The validation split for early stopping is a random 10 % of the training samples, so
it is in-distribution and cannot detect the held-out degradation. Six held-out voxels give
per-radius-bin conclusions of only 1-2 voxels each. The held-out voxels come from the same
sample and distribution as the training voxels; nothing here says how the network behaves on a
different detector geometry or crystal structure. Test offsets are 10 per voxel and magnitude,
so per-voxel medians rest on 40 cases.

**Open**: (i) per-peak weights or attention pooling and per-peak Jacobian-times-residual features
to let the net reproduce Gauss–Newton's parallax use at large r⊥; (ii) more training voxels and
peaks per voxel, or fine-tuning on a new voxel; (iii) recalibrating the covariance on held-out
voxels; (iv) a one-voxel simulated image to run Adam/MC at r⊥ ≈ 400 µm (still missing for ManyGrains).
