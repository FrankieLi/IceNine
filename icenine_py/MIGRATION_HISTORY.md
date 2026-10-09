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

## Toy Orientation NN — Stage 1: Network Inputs (done 2026-09-29; items 3–4 deferred)

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
- The MC baseline's search box is set from the true perturbation magnitude (1.5 x |delta|), as
  in the HP-sweep protocol, so MC is told |delta|; Gauss-Newton is not.
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
-0.87529). The network runs on the training device for inference too (float32); predictions are
moved to the CPU and the error statistics are computed there in float64.

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
`far_test_bayes.npz` with `scripts/summarize_bayes_npz.py --bayes ... --test ... --out ...` (the trainings started before it
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

### Stage 3 — NLL diagnostic (2026-09-29/30)

Question: why does `PeakSetNet` (mean+max pooling, measurement features; round r3) fix the
perpendicular axes on voxel 0 but leave the stage axis z at the prior? Hypothesis: a Gaussian-NLL
pathology. The NLL gradient on the mean is Σ⁻¹(μ − y), so once σ_⊥ ≈ 0.01-0.02° and σ_z ≈ 0.5°
the z gradient is ~1000x weaker than the ⊥ ones and z never leaves "don't know" (Seitzer et al.
2022). Test: change only the loss. Data, architecture, optimiser as r3 (voxel 0, 9000 train,
`--arch set --pool meanmax --lr 3e-4 --clip 1 --cosine --epochs 60`, batch 32, `--device mps`,
seed 0), from `benchmarks/toy_orientation_stage3/nll/res_*.json` (one clean run each).

New code: `--loss {nll,decoupled,mse,mse-then-cov,mse-then-nll}` and `--mse-scale` in
`scripts/train_toy_orientation_nn.py`; `mse_deg_loss` and `decoupled_nll_loss` in
`icenine/orientation_nn.py` (MSE on the mean, equal weight per axis, in units of (0.1°)², plus the
NLL of the covariance at stopgrad(mean), so the mean gets exactly the MSE gradient);
`FrameProbeNet` and `--arch probe` in `icenine/toy_orientation_model.py`. The per-epoch log now
prints validation RMS (z, ⊥) and mean predicted σ (z, ⊥). The existing `--beta-nll` weights each
sample by one scalar, so it cannot rebalance z against ⊥ within a sample; it was run anyway.

Median angle (deg) / z RMS / ⊥ RMS / predicted σ_z / mean Mahalanobis², each at |δ| = 0.1 / 0.25 /
0.5 / 1.0 (30 cases each):

| run | median angle | z RMS | ⊥ RMS | net σ_z | Mahalanobis² |
|---|---|---|---|---|---|
| r3 (NLL, meanmax) | 0.044/0.140/0.194/0.529 | 0.052/0.149/0.241/0.576 | 0.008/0.008/0.008/0.010 | 0.58/0.57/0.49/0.36 | 2.5/2.3/2.2/4.3 |
| E1 β-NLL 0.5 | 0.099/0.227/0.383/0.807 | 0.074/0.147/0.255/0.585 | 0.046/0.120/0.213/0.356 | 0.49/0.48/0.46/0.43 | 0.1/0.4/1.1/4.0 |
| E1 β-NLL 1.0 | 0.095/0.226/0.387/0.792 | 0.073/0.146/0.249/0.586 | 0.049/0.118/0.207/0.350 | 0.49/0.48/0.46/0.43 | 0.5/0.8/1.7/5.6 |
| E2 decoupled, meanmax | 0.028/0.023/0.028/0.034 | 0.026/0.024/0.027/0.040 | 0.012/0.010/0.012/0.014 | 0.029/0.028/0.029/0.029 | 3.6/2.5/2.7/4.4 |
| E3 MSE 30 epochs, then decoupled | 0.025/0.028/0.032/0.033 | 0.019/0.027/0.030/0.041 | 0.013/0.015/0.014/0.015 | 0.029/0.029/0.030/0.031 | 2.7/4.0/3.4/4.8 |
| E5 decoupled, `--pool all` | 0.020/0.026/0.026/0.030 | 0.016/0.024/0.027/0.040 | 0.011/0.011/0.009/0.011 | 0.027/0.028/0.028/0.028 | 3.2/3.4/2.4/4.3 |
| E4 probe (frame + dω*/dδ only, MSE) | 0.085/0.081/0.076/0.087 | 0.023/0.035/0.038/0.046 | 0.057/0.063/0.053/0.060 | (fixed) | |
| fc (Stage 1) | 0.023/0.024/0.036/0.049 | 0.032/0.027/0.037/0.058 | 0.007/0.008/0.016/0.019 | | |
| Gauss–Newton | 0.040/0.024/0.035/0.036 | 0.045/0.032/0.044/0.045 | 0.003/0.003/0.003/0.003 | | |

Findings.

- **The hypothesis is confirmed.** Changing only the loss takes z from the prior (RMS 0.05-0.58°,
  σ_z ≈ 0.5°) to 0.024-0.040° with the same architecture, data and pooling. E2 and E3 reach
  median 0.023-0.034° at every magnitude: at or below fc at 0.5 and 1° (0.036, 0.049) and at
  0.25-1° on par with Gauss–Newton (0.024-0.036). E3's MSE-only phase already had validation z
  RMS 0.044° at epoch 20, before any covariance was fitted (`log_msecov.txt`), so the
  information was reachable with mean pooling and max pooling in place; only the NLL gradient
  scaling was withholding it.
- **Per-sample β-NLL does not help, as expected.** Its final-epoch validation z RMS is 0.43 and
  0.42° (`log_beta05.txt`, `log_beta10.txt`), the prior. The table rows for E1 are worse than
  r3 only because the best-validation checkpoint is epoch 1-2: the β weight `det(cov)^(β/3)`
  changes as σ shrinks, so the validation loss is not comparable across epochs and checkpoint
  selection picks the start. At epoch 60 their ⊥ RMS is 0.008° (as r3), z unchanged.
- **The frame information survives mean pooling.** The probe (a per-peak MLP of the mean frame
  offset and dω*/dδ, mean-pooled, no pixels) recovers z to 0.023-0.046° RMS; it has no pixel
  information, so ⊥ stays at 0.05-0.06° as designed.
- **σ_z does not balloon; it collapses correctly.** In the NLL runs σ_z stays at 0.4-0.6°
  throughout (a calibrated "don't know", Mahalanobis² 2-4) while σ_⊥ falls to 0.007-0.010°.
  With the decoupled loss σ_z falls to ≈ 0.03° along with the error (Mahalanobis² 2.4-4.4).
- **Costs.** ⊥ RMS is 0.009-0.014° against r3's 0.008° (a mild trade for balancing the axes) and
  the 1° z RMS is 0.040°. Exact Bayes remains 5-10x lower (0.004-0.007°).
- `--pool all` with the decoupled loss (E5) is the best voxel-0 run (median 0.020-0.030°), a
  modest gain over meanmax (0.023-0.034°): with the loss fixed, sum pooling is no longer needed.

Far voxel (ManyGrains 77, r⊥ = 399 µm), decoupled loss, batch 64, 60 epochs, one run each:

| run | median angle | z RMS | ⊥ RMS | net σ_z | Mahalanobis² |
|---|---|---|---|---|---|
| NLL, `--pool all` (Step 2) | 0.015/0.013/0.015/0.018 | 0.015/0.014/0.017/0.019 | 0.010/0.010/0.009/0.016 | 0.015 (at 0.25) | 3.3/3.7/3.4/5.6 |
| decoupled, meanmax | 0.019/0.021/0.023/0.024 | 0.014/0.015/0.021/0.023 | 0.011/0.014/0.013/0.035 | 0.016/0.016/0.017/0.019 | 2.5/3.4/3.5/8.4 |
| decoupled, `--pool all` | 0.014/0.016/0.020/0.024 | 0.016/0.014/0.019/0.021 | 0.008/0.010/0.010/0.017 | 0.012/0.013/0.015/0.018 | 4.5/3.7/3.6/5.5 |
| fc | 0.009/0.008/0.009/0.020 | 0.009/0.007/0.006/0.016 | 0.004/0.004/0.006/0.017 | | |
| Gauss–Newton | 0.012/0.011/0.009/0.008 | 0.013/0.010/0.010/0.008 | 0.003/0.003/0.003/0.003 | | |

On the far voxel the plain-NLL set net with `--pool all` was already good, and the decoupled loss
does not improve it (0.014-0.024 vs 0.013-0.018°, slightly worse at 0.5-1°), but it does let
plain `meanmax` pooling (which was at the prior in z there too: σ_z 0.46, z RMS 0.07-0.49°)
reach 0.019-0.024°. So the loss is what was missing for meanmax; sum pooling was compensating
for it. It is not a further gain over sum pooling.

Caveats: one seed per run; test sets of 30 cases per magnitude; no run of a full network with
MSE only (E3's first phase and the probe cover it); the mixing scale (0.1°) was not tuned; the
β-NLL rows use the epoch 1-2 checkpoint (see above). Not tried: decoupled loss on the 30-voxel
multi-voxel data (where the network does not exploit parallax; the loss is the first thing to
change there).

### Stage 3 — multi-voxel retrain with the decoupled loss (2026-09-30)

Question: does the decoupled loss (MSE on the mean + NLL of the covariance at stopgrad(mean),
`--mse-scale 0.1`) fix the 30-voxel network's gap to Gauss–Newton (GN), parallax use and held-out
overconfidence? Same data as Step 3. New: `--val-voxels N` in `scripts/train_toy_orientation_nn.py`
(`split_by_voxel` in `icenine/orientation_nn.py`, tested): N of the 24 training voxels, one per
r⊥ stratum drawn with `--seed`, are removed from training and used only for early stopping /
best-checkpoint selection (the test file's six held-out voxels are never used for selection).
With `--seed 0` and N = 4 the validation voxels are **6, 11, 19, 24 (r⊥ = 99, 181, 326, 410 µm)**;
training is then 20 voxels / 10 000 samples, validation 2 000 samples. Default (`--val-voxels 0`)
is still the random 10 % sample split.

All runs: `--arch set --head offset --device mps --lr 3e-4 --clip 1 --cosine --batch-size 64
--epochs 60` (the settings of the previous multi runs), seed 0, one run each. Results from
`benchmarks/toy_orientation_stage3/multi_{log,res,pred}_<run>.*`. R3 (150 epochs) was **not** run:
60 epochs is not under-trained (see below). R5 was added to separate the loss from the val split.
Median angle (deg) at |δ| = 0.1/0.25/0.5/1.0; "in-dist" = the 24 training-set voxels (960 test
cases; for the voxel-val runs 4 of these 24 voxels were unseen in training), "held-out" = the 6
test-file voxels (240 cases).

| run (loss, pool, val split; best epoch) | group | median angle | z RMS | ⊥ RMS | Mahalanobis² |
|---|---|---|---|---|---|
| previous: NLL, all, sample-val, 60 ep | in-dist | 0.030/0.030/0.032/0.038 | 0.039/0.041/0.046/0.049 | 0.009/0.009/0.010/0.014 | 4.0/3.5/3.2/4.1 |
| | held-out | 0.036/0.043/0.046/0.055 | 0.030/0.041/0.047/0.055 | 0.019/0.020/0.020/0.028 | 19.5/21.6/17.4/17.8 |
| previous: NLL, all, sample-val, 150 ep | in-dist | 0.029/0.026/0.029/0.032 | 0.038/0.039/0.041/0.042 | 0.007/0.006/0.007/0.011 | 4.1/3.6/3.5/5.7 |
| | held-out | 0.034/0.037/0.039/0.053 | 0.031/0.043/0.051/0.054 | 0.018/0.016/0.018/0.024 | 31.1/24.5/26.3/34.3 |
| R1: decoupled, all, voxel-val (ep 12) | in-dist | 0.040/0.045/0.053/0.071 | 0.041/0.047/0.053/0.060 | 0.018/0.020/0.026/0.041 | 1.6/1.9/2.2/3.3 |
| | held-out | 0.040/0.047/0.057/0.068 | 0.033/0.040/0.055/0.062 | 0.022/0.026/0.032/0.038 | 1.8/2.8/3.7/3.5 |
| R2: decoupled, meanmax, voxel-val (ep 20) | in-dist | 0.049/0.048/0.055/0.070 | 0.057/0.054/0.056/0.063 | 0.019/0.019/0.023/0.036 | 2.7/2.4/2.5/3.9 |
| | held-out | 0.051/0.060/0.060/0.060 | 0.053/0.058/0.061/0.066 | 0.022/0.026/0.026/0.029 | 2.9/3.4/3.6/3.5 |
| R4: NLL, all, voxel-val (ep 18) | in-dist | 0.051/0.045/0.050/0.061 | 0.054/0.053/0.056/0.066 | 0.019/0.018/0.019/0.026 | 2.5/2.2/2.0/2.6 |
| | held-out | 0.045/0.048/0.054/0.070 | 0.047/0.050/0.055/0.065 | 0.024/0.021/0.028/0.038 | 3.2/3.0/4.7/5.6 |
| R5: decoupled, all, sample-val (ep 58) | in-dist | 0.033/0.031/0.035/0.039 | 0.037/0.038/0.041/0.041 | 0.013/0.012/0.014/0.019 | 3.1/3.0/3.3/3.8 |
| | held-out | 0.036/0.043/0.050/0.057 | 0.033/0.041/0.053/0.052 | 0.017/0.020/0.027/0.032 | 4.8/7.2/11.1/10.4 |
| Gauss–Newton | in-dist | 0.012/0.015/0.012/0.014 | 0.020/0.024/0.022/0.022 | 0.003/0.003/0.003/0.004 | |
| | held-out | 0.014/0.012/0.012/0.010 | 0.017/0.020/0.020/0.016 | 0.004/0.003/0.004/0.004 | |

Per-voxel median error (all magnitudes) vs r⊥, correlation over all 30 voxels (net) and median
over voxels in four r⊥ bins (0-130 / 130-260 / 260-390 / 390-510 µm, the 20-24 training voxels):

| run | corr(err, r⊥) | bin medians | median of per-voxel medians: 20 train / 4 val / 6 held-out voxels |
|---|---|---|---|
| previous NLL, all, 60 ep | +0.37 | 0.028/0.031/0.033/0.036 | 0.030 / 0.035 / 0.044 |
| R1 decoupled, all, voxel-val | +0.27 | 0.043/0.050/0.050/0.049 | 0.048 / 0.067 / 0.053 |
| R2 decoupled, meanmax, voxel-val | +0.42 | 0.047/0.059/0.058/0.058 | 0.055 / 0.062 / 0.056 |
| R4 NLL, all, voxel-val | +0.38 | 0.043/0.049/0.053/0.057 | 0.052 / 0.064 / 0.052 |
| R5 decoupled, all, sample-val | +0.21 | 0.033/0.032/0.035/0.035 | 0.034 / 0.035 / 0.044 |
| Gauss–Newton | -0.63 | 0.022/0.013/0.011/0.008 | 0.012 / 0.019 / 0.014 |

Answers.

- **(a) Gap to GN: not closed.** The best decoupled run (R5, same protocol as before) is 0.031-0.039°
  in-distribution against 0.026-0.038° for plain NLL: no change. GN is 0.012-0.015° (2.5-3x
  better). On voxel 0 the loss took z from the prior to GN level; on 30 voxels z was already
  learned with plain NLL and `--pool all` (z RMS 0.04°), so there was no NLL pathology left to fix.
- **(b) Parallax: not used.** No net run shows error falling with r⊥ (corr +0.21 to +0.42, bin
  medians flat or rising; GN -0.63, 0.022 → 0.008°). z RMS is 0.037-0.041° (R5), 2x GN's 0.020°; no
  run has z RMS below 0.02° at large r⊥. The loss was not the limiting factor; the pooled
  representation is (see Open items in Step 3).
- **(c) Held-out overconfidence: calibration improves, accuracy does not.** With the voxel-wise val
  split the held-out Mahalanobis² is 1.8-3.7 (R1), 2.9-3.6 (R2) and 3.2-5.6 (R4) against 17-34
  before, i.e. calibrated. The loss is not what did it: R4 (plain NLL) calibrates as well as R1.
  R5 (decoupled loss, sample-val) is between: 4.8-11.1. So the gain comes from the checkpoint
  selection. But the voxel-val checkpoints are early (epoch 12/20/18 of 60; the val-voxel loss is
  noisy, z RMS on the 4 val voxels bounces between 0.07 and 0.086 over epochs 5-60 while the
  train loss keeps falling) and those nets are *worse* everywhere, by 1.3-2x
  (in-dist 0.040-0.071 vs 0.033-0.039 for R5; held-out 0.040-0.070 vs 0.036-0.057). The
  overconfidence was partly the sharpening of σ during late training (σ_z 0.06 at epoch 12 vs
  0.04 at epoch 58 on the val voxels; R1 log), which is not accompanied by better error on unseen
  voxels. What a covariance head can learn from 20-24 voxels is limited: it cannot tell a new voxel
  from a training voxel, so it fits the in-distribution error level; the held-out (and val-voxel)
  error is ~1.2-1.6x larger and the Mahalanobis² rises to 5-11 (R5) once σ has shrunk to the
  training-voxel level. The voxel-val runs avoid this by stopping at a wider σ, at the price of
  accuracy.
- **Pooling.** meanmax (R2) is 10-20 % worse than `all` (R1) under the decoupled loss (z RMS 0.053-0.066 vs
  0.033-0.062); the far-voxel result (loss makes sum pooling unnecessary) does not carry over to
  the multi-voxel net.
- **Recommendation.** For accuracy use R5's protocol (sample-val, late checkpoint); for honest σ on
  new voxels either recalibrate on a held-out set of voxels or use voxel-val, but the noisy
  4-voxel validation loss picks a poor checkpoint. A smoother selection criterion (validation
  z/⊥ RMS instead of the training loss, more validation voxels, or averaging the last epochs) is
  the obvious next step.

Caveats: one seed per run; only 6 held-out voxels (1-2 per r⊥ bin); the "in-dist" test group of
the voxel-val runs contains the 4 val voxels, which those nets never trained on (the third column
of the second table separates them: the 4 val voxels are the worst-off group in R1/R4, 0.064-0.067);
the val voxels influence checkpoint selection, so they are not a clean test either. Best-epoch
selection on a noisy loss confounds "voxel-wise split" with "early checkpoint"; a last-epoch
evaluation of R1/R4 was not run. R3 (150 epochs) was not run.

## Toy Orientation NN — Architecture (parallax) (2026-09-30)

**Architecture branch summary (Steps 1-4 and 3b).** Putting the physics into the network (a learned Gauss–Newton layer, `GNLayerNet`) removes the
3x gap of the pooled set net on clean multi-voxel data: it is at GN's accuracy on seen voxels, 1.1-1.7x GN on unseen voxels (see the
re-check under Step 2), with GN's falling error vs r⊥ and Mahalanobis² 2-4. Detector pairing (Step 3) neither helps nor hurts on clean data.
With realistic data (deliberately simulated nuisances added to exact synthetic windows, not faulty data; Step 4: neighbour/twin spots, spurious blobs, missing spots, edge jitter) plain GN degrades ~20x (median 0.22-0.33° vs
0.012-0.014°) and a Huber-robust GN cuts that by about a third, but a GNLayerNet *trained on realistic windows* reaches 0.06-0.08° (median) on the same data, 2-3x better
than robust GN, while a net trained on clean data only is not better than robust GN and barely better than plain GN. On pixel-level noise alone (no distractor spots) robust GN is as good as or
better than the net. Pairing gives at most a few per cent on realistic data (Step 3b fixed an inert pairing MLP in `GNLayerNet`; the fixed pairing is still within seed spread of the unpaired net). Remaining errors on distractor data are heavy-tailed (perp RMS 0.05° vs 0.004°
clean). Details and caveats below.

Branch `feature/nn-orientation-arch`. Goal (plan `warm-plotting-platypus`): an architecture whose
multi-voxel error falls with r⊥ like Gauss–Newton (GN), i.e. that uses parallax, and that survives
realistic nuisance signal. Data and protocol as in Stage 3 step 3 (30 ManyGrains voxels, 24 train /
6 held out, 0.1/0.25/0.5/1.0°); new runs use `--val-voxels 4` (voxels 6, 11, 19, 24 [Correction 2026-10-01: seed 0 only; seed 1 uses 3, 11, 20, 29]), `--ema 0.998
--checkpoint ema` (EMA weights at the last epoch instead of best-epoch selection), 60 epochs, lr 3e-4,
clip 1, cosine, batch 64, mps, decoupled loss, seeds 0 and 1. Numbers are means of the two seeds,
from `benchmarks/toy_orientation_arch/*res*.json` (tables: `*_summary.txt`,
`scripts/summarize_results.py`).

### Step 1: diagnostics without training

New code: `extract_measurements(..., detectors=[..])` (GN on one detector), `CentroidGaussNewton.information`
(J^T W J), `.solve_linear` (one undamped GN step from nominal), `orientation_eval.pair_index` /
`nominal_offsets`, `scripts/arch_diagnostics.py`, `scripts/dataset_problems.py`,
`scripts/make_dataset_aux.py` (sidecar with per-peak exact nominal offsets and pair index, so the
existing gitignored datasets are reused without re-rendering), `--aux/--subpixel` in the trainer.
Outputs: `diag_{stage1,far,multi}.json`, `diag_multi_summary.txt`.

**Pairing.** ROI sets have 105-130 entries (mean 119: 63 on detector 0, 56 on detector 1). Entries
sharing (reflection, ω-branch) on the two detectors are the same ray: 50-62 pairs per voxel, i.e.
94 % of entries are paired (voxel 0: 53 pairs of 113; voxel 77: 53 of 118). On average 111 of 119
entries per test sample have their partner also recorded.

**GN per detector, 30 multi-voxel test voxels** (median over voxels in r⊥ bins; median angle / z RMS / ⊥ RMS, degrees):

| r⊥ bin (µm) | both detectors | detector 0 only | detector 1 only |
|---|---|---|---|
| 0-130 | 0.019 / 0.032 / 0.0029 | 0.022 / 0.032 / 0.0052 | 0.022 / 0.032 / 0.0032 |
| 130-260 | 0.013 / 0.017 / 0.0031 | 0.015 / 0.019 / 0.0052 | 0.018 / 0.024 / 0.0034 |
| 260-390 | 0.011 / 0.012 / 0.0033 | 0.014 / 0.015 / 0.0048 | 0.011 / 0.013 / 0.0037 |
| 390-510 | 0.0076 / 0.0064 / 0.0044 | 0.011 / 0.0079 / 0.0054 | 0.0099 / 0.010 / 0.0044 |
| corr(error, r⊥) | -0.63 | -0.58 | -0.76 |
| corr(z RMS, r⊥) | -0.84 | -0.79 | -0.88 |

- The stage axis z is the error. GN's z RMS falls 0.032 → 0.006° from r⊥ = 0 to 500 µm (corr -0.84)
  while ⊥ RMS is flat (0.003-0.004°). This is the parallax: the Fisher sigma_z from J^T W J is 0.023° → 0.007°
  (corr -0.96 with r⊥) and the weakest eigenvector of J^T W J is ≥ 95 % z at every radius.
- Each detector alone reproduces the parallax (z RMS corr -0.79/-0.88). Both detectors together
  improve ⊥ by ~1.5x (0.003 vs 0.005 on detector 0) and z by 0-20 %, not by a factor: the two
  detectors are mostly redundant for z on these data (ΔL ≈ 2 mm is small next to the rotation-axis lever arm).
- **Conditioning:** cond(J^T W J) = 48 (r⊥ < 130 µm) → 5 (r⊥ > 390 µm) in sigma-normalised units:
  well posed from pixels and frames alone at every radius, poorer near the axis because the pixel
  columns of J for z scale with r⊥.
- **Linearisation is not a limit:** one undamped GN step from nominal (`solve_linear`) reaches the converged GN
  error (median angle 0.0194/0.0138/0.0120/0.0076 vs 0.0191/0.0127/0.0110/0.0076 by bin, corr -0.63). So a network layer
  that solves the linear normal equations at the nominal point can in principle reach GN.

**Sub-pixel bias.** `measurement_features` measures relative to the window centre; the exact nominal
centroid differs by the fractional position (0-1 px, a fixed per-peak number) and the nominal crossing
is off the frame centre by -0.5..0.5 frames. A test confirms that feature minus `nominal_offsets` is
exactly (lit centroid - exact nominal centroid). Retraining the multi-voxel set net with the corrected
measurement (`--subpixel`) versus the same protocol without it (2 seeds each; median angle / ⊥ RMS at
0.1/0.25/0.5/1.0°):

| run | group | median angle | z RMS | ⊥ RMS | Mahalanobis² |
|---|---|---|---|---|---|
| set, window-centre measurement | in-dist | 0.034/0.035/0.040/0.044 | 0.039/0.042/0.047/0.046 | 0.015/0.015/0.017/0.025 | 3.6/4.1/4.5/5.6 |
| | held-out | 0.038/0.042/0.053/0.058 | 0.036/0.042/0.054/0.054 | 0.019/0.022/0.026/0.032 | 5.5/7.4/10.2/11.0 |
| set, `--subpixel` | in-dist | 0.033/0.035/0.033/0.043 | 0.039/0.039/0.038/0.041 | 0.011/0.013/0.014/0.023 | 3.1/3.6/3.4/5.1 |
| | held-out | 0.029/0.034/0.039/0.053 | 0.040/0.037/0.038/0.049 | 0.012/0.015/0.019/0.029 | 4.0/4.8/6.0/9.3 |
| Gauss–Newton | in-dist / held-out | 0.012-0.015 / 0.010-0.014 | 0.020-0.024 / 0.016-0.020 | 0.003-0.004 | |

Correlation of per-voxel error with r⊥: +0.29 for both nets (GN -0.63). **Result:** the bias was a real
but small contributor: ⊥ RMS improves by 25-35 % and held-out median by 10-25 %, calibration improves
(held-out Mahalanobis² 4-9 vs 5.5-11); z RMS (0.04°) and the r⊥ trend do not change, so it does not
explain the 3x gap to GN. The baseline here (EMA weights, voxel-wise validation) reproduces the earlier
best set net (held-out 0.036/0.043/0.050/0.057 -> 0.038/0.042/0.053/0.058), so the comparison is like for like.
Caveats: 2 seeds; per-seed spread not tabulated (only the mean is quoted); 6 held-out voxels.

### Step 2: learned Gauss–Newton layer (`GNLayerNet`)

**What was built** (`icenine/toy_orientation_model.py`, `--arch gn`). A shared per-peak encoder (window
convolutions + context + measurement features, as `PeakSetNet`) emits, per present peak, a
reliability weight for each of its three measured rows (centroid column, centroid row, frame; a
softplus, initialised to 1) and a correction to the measurement (initialised to 0). The measurement is
y = (lit centroid - exact nominal centroid, frame - exact nominal crossing) in units of the
quantisation sigma (`--subpixel` features, `nominal_offsets`); the Jacobian J comes from the context
(Γ = context[:, 0:6] × 20 px/deg, dω*/dδ = context[:, 6:9] converted to frames/deg). The pooled normal
equations A = Σ JᵀWJ + 1e-3·I, δ = A⁻¹ Σ JᵀW(y + dy) are solved in closed form (3×3 adjugate inverse and
Cholesky in plain tensor ops: differentiable and runs on MPS, which has no `linalg.solve`). The
returned covariance is D A⁻¹ D with a learned diagonal D (calibration). `--gn-iters K` unrolls K
IRLS-style rounds in which the head also sees asinh of the residual y - Jδ and the current δ. The
initial state (W=1, dy=0) **is** one undamped GN step from nominal (`CentroidGaussNewton.solve_linear`).
Trained with the decoupled loss (same protocol as Step 1 runs: 60 epochs, EMA weights, `--val-voxels 4`
for multi-voxel; single-voxel runs use the default random 10 % sample validation).
Tests: layer == weighted least squares; W=1/dy=0 reproduces `solve_linear` on real windows
(|Δδ| < 2e-3°, covariance diag within 3 %); permutation and padding invariance (with and without pairing);
finite gradients.

**Single-voxel sanity** (30 test cases per magnitude; median angle in degrees at 0.1/0.25/0.5/1.0°; one run):

| voxel | method | median angle | z RMS | ⊥ RMS | Mahalanobis² |
|---|---|---|---|---|---|
| 0 (r⊥ 12 µm) | set net (Stage 3, decoupled, all) | 0.020/0.026/0.026/0.030 | 0.016/0.024/0.027/0.040 | 0.011/0.011/0.009/0.011 | 3.2/3.4/2.4/4.3 |
| | **GNLayerNet K=1** | 0.014/0.018/0.018/0.025 | 0.019/0.030/0.024/0.037 | 0.0028/0.0035/0.0033/0.0039 | 2.1/3.8/2.7/3.8 |
| | Gauss–Newton | 0.040/0.024/0.035/0.036 | 0.045/0.032/0.044/0.045 | 0.0027/0.0035/0.0031/0.0033 | |
| 77 (r⊥ 399 µm) | set net (Stage 3, decoupled, all) | 0.014/0.016/0.020/0.024 | 0.016/0.014/0.019/0.021 | 0.008/0.010/0.010/0.017 | 4.5/3.7/3.6/5.5 |
| | **GNLayerNet K=1** | 0.0098/0.0078/0.0086/0.0095 | 0.010/0.012/0.009/0.012 | 0.0025/0.0028/0.0028/0.0036 | 2.9/3.2/2.4/2.8 |
| | Gauss–Newton | 0.012/0.011/0.009/0.008 | 0.013/0.010/0.010/0.008 | 0.0031/0.0032/0.0034/0.0035 | |

(`single_v0_gn.*`, `single_v77_gn.*`.) The net reaches GN's ⊥ accuracy immediately (0.003°), and is
at or below GN in z (on voxel 0 better than GN, because it learns corrections to the quantised
measurements; GN's voxel-0 median at 0.1° is dominated by few cases). Exact Bayes for voxel 0 is
0.004-0.007°, so ~3x headroom remains.

**Multi-voxel, 30 voxels** (2 seeds, means; in-dist = 24 training-set voxels of which 4 are the
validation voxels the net never trained on; held-out = 6 test-file voxels). From
`multi_res_gn_k{1,3}_s{0,1}.json`, `step2_summary.txt`:

| run | group | median angle | z RMS | ⊥ RMS | Mahalanobis² |
|---|---|---|---|---|---|
| set net (Step 1 baseline, window-centre) | in-dist | 0.034/0.035/0.040/0.044 | 0.039/0.042/0.047/0.046 | 0.015/0.015/0.017/0.025 | 3.6/4.1/4.5/5.6 |
| | held-out | 0.038/0.042/0.053/0.058 | 0.036/0.042/0.054/0.054 | 0.019/0.022/0.026/0.032 | 5.5/7.4/10.2/11.0 |
| set net, `--subpixel` | held-out | 0.029/0.034/0.039/0.053 | 0.040/0.037/0.038/0.049 | 0.012/0.015/0.019/0.029 | 4.0/4.8/6.0/9.3 |
| **GNLayerNet K=1** | in-dist | 0.013/0.015/0.012/0.013 | 0.020/0.023/0.021/0.022 | 0.006/0.005/0.004/0.004 | 3.1/3.2/3.1/3.2 |
| | held-out | 0.017/0.016/0.015/0.013 | 0.023/0.022/0.022/0.019 | 0.007/0.006/0.005/0.004 | 3.7/4.0/3.9/3.5 |
| **GNLayerNet K=3** | in-dist | 0.013/0.015/0.012/0.015 | 0.021/0.023/0.021/0.023 | 0.003/0.003/0.003/0.004 | 3.3/3.4/3.1/3.2 |
| | held-out | 0.018/0.016/0.017/0.016 | 0.024/0.024/0.023/0.022 | 0.003/0.003/0.004/0.005 | 4.2/4.0/3.8/3.9 |
| Gauss–Newton | in-dist | 0.012/0.015/0.012/0.014 | 0.020/0.024/0.022/0.022 | 0.003/0.003/0.003/0.004 | |
| | held-out | 0.014/0.012/0.012/0.010 | 0.017/0.020/0.020/0.016 | 0.004/0.003/0.004/0.004 | |

Per-voxel median error vs r⊥ (30 voxels, mean of 2 seeds of the correlation): K=1 **-0.72**, K=3 **-0.64**,
GN -0.63, set nets +0.29. Median per-voxel error: K=1 0.0136 in-dist / 0.0188 held-out; K=3 0.0131 / 0.0169; GN
0.0128 / 0.0143; set net 0.037 / 0.048.

**Finding.** The plan's success criteria are met on clean data: held-out median is 1.0-1.5x GN **in the 2-seed mean (this overstates it per seed, see the re-check below)** (0.013-0.018 vs 0.010-0.014°; the plan allowed 1.5x), the error
falls with r⊥ like GN's, the z RMS (0.02°) equals GN's, and Mahalanobis² is 3.1-4.2 in-dist *and* held-out (the set net's held-out was
5.5-11). Putting the physics in the architecture removes the 3x gap; the val-voxel z RMS (unseen voxels, epoch-60 EMA; from the training logs) drops from 0.074° (set net; 0.045° with --subpixel) to 0.020-0.023° (seed 0). The K=3 variant trades nothing on clean data
and cleans up ⊥ (0.003 vs 0.005-0.007°); K=1's ⊥ is noisier because the freed weight/correction head can trade ⊥ against z.
On clean data iterations add nothing (the linear model is already exact enough); their purpose is outlier rejection (Step 4).
Caveats: 2 seeds, the mean of the four |δ| bins is quoted with no spread; 6 held-out voxels; the held-out gap to GN (1.0-1.5x) is small but
consistently positive at 0.1° (0.017-0.018 vs 0.014). The linear-at-nominal Jacobian ignores curvature: at 1° the layer matches GN (0.013-0.016 vs 0.010-0.014).
Not verified: behaviour on a different geometry/structure.

**Re-check of the "Gauss–Newton parity" claim (Step 4 session, from the saved per-seed jsons).** Ratio of the net's median angle to
GN's, mean of the four |δ| bins, per seed (`multi_res_gn_k1_s*`, `multi_res_gn_k3_s*`; the Step 3 lr 1e-4 and paired runs for comparison):

| run | in-dist s0 / s1 | held-out s0 / s1 | corr(voxel error, r⊥) s0 / s1 |
|---|---|---|---|
| K=1, lr 3e-4 | 0.93 / 1.12 | 1.15 / 1.43 | -0.64 / -0.79 |
| K=3, lr 3e-4 | 1.17 / 0.91 | **1.74** / 1.10 | -0.43 / -0.85 |
| K=3, lr 1e-4 | 0.92 / 1.19 | 1.24 / 1.38 | -0.75 / -0.47 |
| K=3 paired, lr 1e-4 | 0.93 / 1.04 | 1.13 / 1.40 | -0.70 / -0.50 |

(GN: corr -0.63.) So: on the 24 training-set voxels the net is at parity (0.9-1.2x GN); on the 6 unseen voxels it is 1.1-1.7x GN
and the "within 1.5x" criterion fails for one of eight runs (K=3 s0, 1.74x; seed spread alone moves K=3 from 1.10x to 1.74x);
the sign of the r⊥ trend holds in every run (-0.43 to -0.85) but its strength varies by seed. **The accurate statement is
"close to GN on seen voxels, 1.1-1.7x GN on unseen voxels, with GN's parallax trend"**, not exact parity. The 3x gap of the set
net is removed in every run (set net held-out ~3.5x). The commit subject of e6f4072 ("reaches Gauss–Newton parity on 30 voxels")
is therefore an overstatement for held-out voxels; the numbers above are the record.

### Step 3: detector pairing (`--pairing`)

**Implementation (a deviation from the plan's "reflection token").** Instead of a separate token holding both
windows, every entry keeps its own window encoding and is given its partner's: `pair_index` (same
reflection and ω-branch, other detector; `orientation_eval.pair_index`, stored in the aux sidecar, padded
per voxel) gathers the partner's encoding; `f ← f + MLP([f, f_partner, has_partner])`. In the IRLS iterations the weight head
additionally sees asinh of the partner's residual (y - Jδ) and the flag, which is the L1-L2 consistency check:
a spot whose partner disagrees with the current δ is down-weighted. Both detector rows enter the normal equations through their own J rows, as before.
The pair MLP's output layer is zero-initialised (the network starts as the unpaired one); a first attempt without that was unstable
(`multi_res_gnpairv1_k3_s0.*`: per-voxel median error 0.029 vs 0.013 unpaired; training loss oscillating), a
learning rate of 3e-4 with zero-init gave 0.016/0.019 and lr 1e-4 fixed it. To compare like with like, the unpaired K=3 net was re-run at lr 1e-4.
Tests: `pair_index` links exactly the same-reflection, same-branch, other-detector entry (symmetric, equal nominal ω), incl. unpaired/3-detector
cases; GNLayerNet permutation/padding invariance holds with pairing (partner indices permuted consistently).

[Correction 2026-10-01: the pair-encoding MLP described above never trains. `GNLayerNet.forward` applies
`f = f + relu(pair2(relu(pair1(...))))` with `pair2` zero-initialised (`toy_orientation_model.py` lines 398-402, 448), so the
pre-activation of the ReLU is exactly 0, ReLU'(0) = 0, and `pair1`/`pair2` receive exactly zero gradient forever (verified by backward;
`head`/`cov` parameters do get gradient). In every paired run (`gnpair*`) pairing therefore acted only through the head's
asinh(partner residual) input and the has-partner flag; 12,480 of the 115,393 parameters are inert. The statement that zero-initialisation
"fixed" the instability is better read as zero-initialisation switching the mixing path off; only the unstable `gnpairv1` had live mixing.
The Step 3 and Step 4 conclusions on pairing ("no clear benefit") therefore test only the partner-residual input, not encoder-level mixing
of the two detectors. The unit tests missed it: the invariance test re-initialises the zero parameters and the gradient test checks only
finiteness. See `docs/orientation_nn_design.md` Section 3.4. No code change was made in this correction.]

**Clean data, 30 voxels, K=3, 2 seeds** (`step3_summary.txt`; median angle / z RMS / ⊥ RMS at 0.1/0.25/0.5/1.0°):

| run | group | median angle | z RMS | ⊥ RMS | Mahalanobis² |
|---|---|---|---|---|---|
| paired, lr 1e-4 | in-dist | 0.011/0.014/0.012/0.014 | 0.018/0.024/0.023/0.023 | 0.003/0.004/0.004/0.005 | 2.8/3.2/3.1/3.1 |
| | held-out | 0.013/0.015/0.016/0.016 | 0.020/0.022/0.025/0.021 | 0.003/0.004/0.004/0.005 | 4.4/3.9/4.2/4.5 |
| unpaired, lr 1e-4 | in-dist | 0.013/0.016/0.013/0.013 | 0.021/0.025/0.022/0.023 | 0.003/0.003/0.003/0.004 | 2.9/3.2/2.8/3.2 |
| | held-out | 0.016/0.016/0.015/0.015 | 0.027/0.026/0.023/0.021 | 0.003/0.003/0.004/0.004 | 3.7/3.5/3.5/3.6 |
| unpaired, lr 3e-4 (Step 2) | held-out | 0.018/0.016/0.017/0.016 | 0.024/0.024/0.023/0.022 | 0.003/0.003/0.004/0.005 | 4.2/4.0/3.8/3.9 |
| Gauss–Newton | held-out | 0.014/0.012/0.012/0.010 | 0.017/0.020/0.020/0.016 | 0.004/0.003/0.004/0.004 | |

corr(voxel error, r⊥): paired -0.60, unpaired -0.61 (lr 1e-4) / -0.64 (lr 3e-4); GN -0.63. **Finding:** on clean data
pairing neither helps nor hurts beyond seed noise (held-out median 0.013-0.016 vs 0.015-0.016; in-dist 0.011-0.014 vs 0.013-0.016; the
paired net's ⊥ at 1° is slightly worse, 0.005 vs 0.004). This is expected: on clean data the two detectors' rows in the normal equations already
carry the pair information (Step 1: detectors are largely redundant for z), and there are no outliers to reject. The value of pairing, if any,
is in the realistic-data test (Step 4 below). Caveat: 2 seeds; differences of ~0.002° are within seed spread (not tabulated).

### Step 4: distractors and noise

**Terminology: realistic data.** "Realistic" windows are exact synthetic windows with deliberately simulated complexity added to approach real detector data: genuine spots of neighbouring voxels and a Σ3 twin overlapping the target (`neighbours`, more complex rather than noisier) and detector noise (`noise`: missing spots, threshold jitter, hot pixels, spurious blobs); `all` = both; "clean" = neither. Orientation labels are always exact. "Realism-trained" = trained on windows made realistic on the fly (`--realistic-train all`). Earlier versions called this "corruption" (old flags/names remain as aliases; filenames with `corr` = realism-trained); it never meant faulty data. Full definition: `docs/orientation_nn_design.md` Section 2.7.

**Realism layer (verified against `orientation_eval.RealismConfig` / `make_realistic_windows` and the generator).** Applied to the frame-coded windows,
independently per entry:
- *neighbours* (distractor layer, rendered at dataset-generation time by `render_distractor_windows`): the spots of the 2 nearest ManyGrains mic voxels (≈9 µm away) and a Σ3
  twin of the nearest one (`--neighbors 2 --twin`, defaults `--neighbor-p 0.5 --neighbor-sigma-deg 0.3`): each source follows the target's perturbed orientation plus its own random misorientation
  (0.3° total), and is present in a given sample with probability 0.5. Its lit pixels (with their frame codes) are drawn into every target window they overlap; the
  target's own pixels always win (`combine_windows`), so the target's spots are never altered (tested).
- *noise* (`p_miss`, `p_flip`, `p_hot`, `p_blob` = 0.1, 0.05, 0.05, 0.1): a whole spot is not recorded (window zeroed, 10 %); each lit pixel is dropped and each 4-neighbour of a lit pixel is lit
  with probability 0.05 (threshold jitter at the spot edge); one isolated hot pixel (5 %); one spurious 2-4 px blob at a random frame (10 %).
- Variants: `clean`, `neighbours` (layer only), `noise` (pixel-level only), `all` (both). Test windows are made realistic deterministically (`make_realistic_dataset`, fixed seed), so GN and the nets see identical inputs; training windows are made realistic on the fly (`--realistic-train all`).
Not implemented (not in the plan's list either): a per-sample random severity of the nuisances, partial-spot occlusion, intensity/threshold jitter beyond the edge flips, noise correlated across the two detectors.

**Data.** Same 30 voxels, selection and protocol as Step 2/3 but a new dataset (`scripts/toy_orientation_arch_dis_{train,test}.pt`, gitignored, 3.2 GB / 0.3 GB, seed 42; 12000 train / 1200 test
samples; `dis_windows` holds the distractor layer). The `--aux` sidecar of the Step 2/3 data is reused (same voxels, ROI sets and contexts; checked equal. The clean windows and offsets of the test file are identical to `toy_orientation_stage3_multi_test.pt`, the generator seed being the same).

**Compared.** (a) plain GN; (b) Huber-robust GN (`CentroidGaussNewton(huber_c=c)`: IRLS-weighted Levenberg–Marquardt on the normalised residuals); (c) GNLayerNet K=3 trained with `--realistic-train all`;
(d) the same with `--pairing`; (e) GNLayerNet K=3 trained on the *clean* windows of the same file (to separate "has an IRLS head" from "was trained on realistic windows"). Nets: lr 1e-4, 60 epochs, EMA 0.998, `--val-voxels 4` (validation windows made realistic too), decoupled loss, seeds 0 and 1:
```
uv run python scripts/train_toy_orientation_nn.py --head offset --loss decoupled --device mps --clip 1 --cosine --batch-size 64 --epochs 60 --ema 0.998 --checkpoint ema \\
  --lr 1e-4 --arch gn --gn-iters 3 [--pairing] [--realistic-train all] --train scripts/toy_orientation_arch_dis_train.pt --test scripts/toy_orientation_arch_dis_test.pt \\
  --aux scripts/toy_orientation_stage3_multi_aux.pt --val-voxels 4 --eval-variants clean,neighbours,noise,all --seed S --results-json ... --save-predictions ...
```
(`dis_res_gn_k3_corr_s{0,1}`, `dis_res_gnpair_k3_corr_s{0,1}`, `dis_res_gn_k3_cleantrain_s{0,1}`.) GN: `scripts/gauss_newton_baseline.py --realistic V [--huber c]`. All numbers below are produced by
`scripts/summarize_arch_step4.py` from the saved predictions (`step4_summary.{txt,json}`: every cell, per-bin medians, per-seed values), not from the trainer's own tables.
The clean-variant predictions `dis_res_gn_k3_cleantrain_s{0,1}.npz` are byte-identical to `multi_res_gn_k3_lr1e-4_s{0,1}.npz` (the same deterministic run, evaluated again on the realistic test sets); the `_all`/`_neighbours`/`_noise` predictions exist only under the `cleantrain` name, so both sets of files are kept.
**Padding caveat.** The realism layer also hits padded entries of the multi-voxel arrays in all reported runs (`--mask-padding` was not used): they have J = 0, so they do not move the estimate delta, but hot pixels/blobs can make them look "present", so they enter the pooled covariance features and the count n; the Gauss-Newton baseline slices `[:n_pk]`, so the comparison is slightly asymmetric against the net. Use `--mask-padding` for future runs.

**Huber threshold: chosen on a validation split, not on the test set.** 200 samples from four *training-set* voxels (6, 11, 19, 24; `step4_huber_val.txt`), windows made realistic with `all`:
median angle plain GN 0.095°, Huber c = 0.5 / **1** / 2 / 3 / 5: 0.0336 / **0.0325** / 0.0328 / 0.0369 / 0.0426° (on clean validation windows: 0.0179 / 0.0158 / 0.0125 / 0.0134 / 0.0135 vs GN 0.0135°).
c = 1 minimises the realistic-validation median and c = 1 and 2 are within 1 %; c = 1 is used for all variants, with a clean-data cost (below). (An earlier session had picked c per test variant by looking at test results; those
numbers are not used. Logs for c = 0.5-5 on test variants are in the directory but are not part of any table.) Checked afterwards: with `--max-iter 200` (22 % of c = 1 fits hit the default cap of 30 on `all`) every fit converges and the median
error is unchanged (0.156-0.174° vs 0.14-0.22° in-dist/held-out above; `dis_*_huber1it200_*`), so the iteration cap is not what limits the robust fit.

**Results** (test set: 30 voxels, 24 in-dist / 6 held-out, 40 samples each = 10 at each |δ| of 0.1/0.25/0.5/1.0°; 2 seeds for nets, means shown, per-seed in brackets).
Median angle is pooled over the four |δ| bins (it does not depend on |δ| in any row; per-bin values are in `step4_summary.txt`); in-dist = 24 training-set voxels (4 of them
the network's validation voxels), held-out = 6 voxels never seen. RMS and "< 0.1°" (fraction of cases) are held-out; Mah² in-dist / held-out; corr = per-voxel median error vs r⊥ (30 voxels; GN on clean -0.63).

| test set | method | median in-dist | median held-out [seeds] | < 0.1° | z / ⊥ RMS | Mah² | corr(err, r⊥) [seeds] |
|---|---|---|---|---|---|---|---|
| clean | plain GN | 0.0131 | 0.0119 | 1.00 | 0.018 / 0.004 | | -0.63 |
| | Huber GN (c=1) | 0.0125 | 0.0144 | 1.00 | 0.023 / 0.004 | | -0.73 |
| | GNLayerNet, clean-trained | 0.0138 | 0.0152 [0.0143, 0.0161] | 1.00 | 0.025 / 0.004 | 3.0 / 3.6 | -0.75, -0.47 |
| | GNLayerNet, realism-trained | 0.0185 | 0.0193 [0.0176, 0.0210] | 1.00 | 0.028 / 0.004 | 2.1 / 2.2 | -0.72, -0.71 |
| | paired, realism-trained | 0.0185 | 0.0186 [0.0185, 0.0186] | 1.00 | 0.027 / 0.005 | 2.3 / 2.3 | -0.73, -0.79 |
| neighbours | plain GN | 0.2253 | 0.3428 | 0.20 | 0.641 / 0.100 | | -0.11 |
| | Huber GN (c=1) | 0.1375 | 0.2223 | 0.33 | 0.285 / 0.070 | | +0.13 |
| | GNLayerNet, clean-trained | 0.1934 | 0.2706 [0.2916, 0.2496] | 0.23 | 0.513 / 0.091 | 1556 / 1951 | -0.14, +0.04 |
| | GNLayerNet, realism-trained | 0.0576 | 0.0719 [0.0723, 0.0715] | 0.66 | 0.068 / 0.051 | 3.0 / 2.9 | -0.18, -0.09 |
| | paired, realism-trained | 0.0594 | 0.0710 [0.0702, 0.0718] | 0.65 | 0.068 / 0.052 | 3.0 / 2.7 | -0.10, -0.11 |
| noise | plain GN | 0.0545 | 0.0573 | 0.85 | 0.066 / 0.028 | | -0.80 |
| | Huber GN (c=1) | **0.0177** | **0.0215** | 1.00 | 0.030 / 0.005 | | -0.72 |
| | GNLayerNet, clean-trained | 0.0704 | 0.0741 [0.0746, 0.0736] | 0.71 | 0.076 / 0.040 | 223 / 272 | -0.49, -0.27 |
| | GNLayerNet, realism-trained | 0.0241 | 0.0242 [0.0227, 0.0258] | 0.99 | 0.036 / 0.006 | 3.0 / 3.0 | -0.61, -0.71 |
| | paired, realism-trained | 0.0233 | 0.0231 [0.0228, 0.0233] | 1.00 | 0.035 / 0.006 | 3.0 / 3.0 | -0.69, -0.76 |
| all | plain GN | 0.2228 | 0.3311 | 0.20 | 0.617 / 0.100 | | -0.14 |
| | Huber GN (c=1) | 0.1405 | 0.2239 | 0.32 | 0.288 / 0.071 | | +0.14 |
| | GNLayerNet, clean-trained | 0.2033 | 0.2757 [0.3016, 0.2498] | 0.16 | 0.493 / 0.094 | 1462 / 1796 | -0.17, +0.04 |
| | GNLayerNet, realism-trained | **0.0615** | **0.0821** [0.0827, 0.0815] | 0.64 | 0.071 / 0.054 | 2.9 / 2.6 | -0.17, -0.12 |
| | paired, realism-trained | **0.0616** | **0.0762** [0.0734, 0.0789] | 0.62 | 0.071 / 0.055 | 2.8 / 2.5 | -0.12, -0.12 |

Per-bin held-out medians on `all` (|δ| = 0.1/0.25/0.5/1.0°): plain GN 0.307/0.353/0.316/0.333; Huber 0.185/0.249/0.224/0.233; realism-trained net 0.064/0.084/0.073/0.086; paired 0.064/0.083/0.076/0.086; clean-trained net 0.264/0.288/0.256/0.278. In-dist seed values for the realism-trained nets agree to ≤ 0.004° (`all`: 0.0614/0.0615 unpaired, 0.0625/0.0608 paired).

**Findings.**
1. *Distractor spots break GN and the robust variant helps only partly.* Plain GN goes from 0.012-0.013° to 0.22-0.34° (≈20x; only 20-34 % of cases within 0.1°); Huber GN cuts that by about a third (0.14° in-dist, 0.22° held-out) and loses the parallax
   trend (corr +0.13/+0.14). Pixel-level noise alone costs plain GN 4x (0.055°) and Huber GN recovers nearly all of it (0.018-0.022°).
2. *A GNLayerNet trained on realistic windows is 2-3x better than robust GN on distractor data* (`all`: 0.0615° in-dist / 0.082° held-out vs 0.1405 / 0.224; 64 % vs 32 % of held-out cases below 0.1°) and its
   covariance stays calibrated (Mahalanobis² 2.5-3.0, in-dist and held-out). The comparison is fair in the sense that the robust baseline's only free parameter was set on separate validation voxels, but note the net sees the nuisance family it
   is tested on (see caveats). Where the *only* nuisance is pixel noise the robust GN is as good or better (0.0215 vs 0.0242° held-out, 2 seeds, every seed worse than Huber): the plan's criterion "matches or beats robust GN" is **met for neighbours/all, not for noise**.
3. *Realism-aware training is what matters, not the IRLS head.* The same architecture trained on clean windows is no better than plain GN on `neighbours`/`all` (0.20-0.28° vs 0.22-0.33°; worse than Huber GN) and its Mahalanobis² explodes (1500-2000: confidently wrong), and on `noise` it is worse than plain GN (0.074 vs 0.057°). 
   Realism-aware training costs accuracy on clean data (held-out 0.0176-0.0210° vs 0.0143-0.0161° clean-trained vs 0.012° GN: 1.2-1.4x worse than the clean-trained net). Robust GN also costs ≈ 20 % on clean held-out (0.0144 vs 0.0119°) at c = 1.
4. *Generalisation to unseen voxels degrades under distractors.* Clean: held-out ≈ in-dist (0.0193 vs 0.0185°). `all`: 0.082 vs 0.0615° (1.3x); `neighbours` 0.072 vs 0.058° (1.25x). The gap is not specific to the net: plain GN's is 1.5x (0.331 vs 0.223°) and Huber GN's 1.6x (0.224 vs 0.141°), and on clean data GN's held-out error is *lower* than in-dist, so under distractors the 6 held-out voxels are simply harder, and the net's gap (1.3x) is no worse than GN's. It is still a gap that the clean data did not show.
5. *Pairing: no clear benefit.* `neighbours`, `noise`, in-dist `all`, clean: differences ≤ 0.001° (inside the seed spread). Only held-out `all` shows paired < unpaired (0.0734/0.0789 vs 0.0827/0.0815, ≈ 7 %, both paired seeds below both unpaired seeds), with 2 seeds and 6 voxels that is suggestive, not established. The partner-residual consistency check does not buy the large gain one might expect; one hypothesis (untested) is that a distractor spot from a neighbour voxel 9 µm away
   lands consistently on both detectors, so the pair check cannot reject it.
6. *Parallax is lost under distractors for every method.* Error-vs-r⊥ correlation drops from -0.6…-0.8 (clean, noise) to -0.1…-0.18 (neighbours, all) for plain GN and for the nets alike, and is positive for Huber GN. The remaining error is heavy-tailed: ⊥ RMS 0.05° vs 0.004° on clean (median 0.06-0.08°, but 34-38 % of held-out cases above 0.1°).
   The net's gain over Huber GN is largest in z (RMS 0.071 vs 0.288°, 4x) and smaller in ⊥ (0.054 vs 0.071°, 1.3x): the ⊥ tail, not the z error, is what the net leaves.

**Success criteria of the plan (realistic data).** "Beats plain GN": yes on all three realistic sets (≥ 2.4x). "Matches or beats robust GN": yes on neighbours/all (2-3x), no on noise-only (Huber 0.9x of the net's error). "Mahalanobis² 2-5 held-out": yes (2.5-3.0) for the realism-trained nets, no for clean-trained ones.

**Caveats / not verified.**
- 2 seeds, 6 held-out voxels, 10 samples per voxel and |δ|; the seed spread of the held-out median is up to 0.003° (clean), 0.0055° (realism-trained nets on `all`: paired 0.0734 vs 0.0789) and 0.04° (clean-trained net on `all`); no confidence intervals.
- Train and test nuisances come from the same family and parameters (same p's, same neighbour set per voxel, same twin construction); the test distractors are new random draws, not new *kinds* of nuisance. The net's advantage over robust GN may shrink out of family (different rates, other grains, intensity effects); not tested. A real-data check is the real test.
- The neighbour layer uses at most 3 sources per voxel from the same mic; a voxel with a dense neighbourhood or genuinely overlapping twin spots may be harder.
- Robust GN is one fixed construction (Huber on normalised residuals, one start at nominal, LM, c = 1); a different robust loss, multi-start, or explicit spot-to-source assignment could narrow the gap. Threshold selection used only 4 voxels × 50 samples.
- Net lr/K were taken from Step 3 (lr 1e-4, K = 3), not re-tuned for realistic data.
- Not analysed: what the weight head learned (does it down-weight distractor entries?); the dependence on the number of distractor pixels in a window; the tail failures.
- The seed-0 realistic unpaired net was trained earlier (before the interruption); a rerun with the same command reproduced its first 19 epochs exactly, so it is the same configuration; its result was kept.

**Deferred.** Out-of-family nuisance sweeps (rates, intensity); a weight-head analysis; per-sample severity as training augmentation; more seeds for the pairing question; Step 1-3 style diagnostics on realistic data; real data.

### Step 3b: pairing fix (2026-10-01)

**Bug.** See the 2026-10-01 correction under Step 3: with `pair2` zero-initialised and `f = f + relu(pair2(relu(pair1(...))))`, the pair MLP got exactly zero gradient forever.

**Fix** (`GNLayerNet.forward`): `f = f + hm * pair2(relu(pair1(cat[f, f_partner, hm])))`. No ReLU after the zero-initialised `pair2`, so step 0 is still the unpaired network and the update is trainable; masked by the has-partner flag `hm`, so unpaired entries are unchanged. Gradient reaches `pair2` as soon as `head2`/`cov2` (also zero-initialised) have left zero, and `pair1` once `pair2` has. No gate or LayerNorm was added: the simplest variant trained stably at lr 1e-4 in all four runs (no NaN, no divergence; the loss curves look like the unpaired ones), so the gate variant was not tried. Parameter count unchanged (115,393). The old behaviour is not kept as an option.

**Tests** (`tests/test_orientation_baselines.py::TestGNLayer`): `gradients_flow_and_are_finite[gn, gn_paired, set]` now requires a non-zero finite gradient for every parameter of `GNLayerNet` (unpaired, paired, T=3) and `PeakSetNet` after three small steps off the zero-initialised layers; `pair_mlp_trains_from_the_zero_init`; `paired_net_equals_unpaired_net_at_initialisation` (extra head inputs zeroed, shared weights copied; `pair1` irrelevant while `pair2` = 0; path live once `pair2` != 0). The gradient tests and `pair_mlp_trains...` fail on the pre-fix code (checked by stashing the model change). The permutation/padding tests still pass with pairing.

**Reruns** (same data, settings and seeds as the original runs: `--head offset --arch gn --gn-iters 3 --pairing --loss decoupled --device mps --lr 1e-4 --clip 1 --cosine --batch-size 64 --epochs 60 --ema 0.998 --checkpoint ema --val-voxels 4`; clean: `scripts/toy_orientation_stage3_multi_{train,test}.pt`, `--aux scripts/toy_orientation_stage3_multi_aux.pt`, `--extra gn=benchmarks/toy_orientation_stage3/multi_pred_gauss_newton.npz`; realistic: `scripts/toy_orientation_arch_dis_{train,test}.pt`, same aux, `--realistic-train all --eval-variants clean,neighbours,noise,all`; seeds 0, 1). Outputs `benchmarks/toy_orientation_arch/{multi_res_gnpairfix_k3_lr1e-4,dis_res_gnpairfix_k3_corr}_s{0,1}.*`; tables `pairfix_clean_summary.txt` (`summarize_results.py`) and `step4_pairfix_summary.{txt,json}` (`summarize_arch_step4.py --runs ... "paired-fixed=dis_res_gnpairfix_k3_corr"`). The old runs are untouched. (The first attempt at the realistic runs was killed by a background time limit at epoch 35 and rerun from scratch; the reported runs are the complete ones.)

**Clean multi-voxel, T=3, lr 1e-4** (30 voxels; medians at |δ| = 0.1/0.25/0.5/1.0°; in-dist includes the 4 validation voxels; per-seed values `s0 | s1`):

| run | group | median angle | z RMS | ⊥ RMS | Mahalanobis² |
|---|---|---|---|---|---|
| paired, fixed | in-dist s0 | 0.014/0.016/0.013/0.014 | 0.021/0.026/0.023/0.024 | 0.003/0.003/0.003/0.004 | 2.7/3.1/2.9/3.3 |
| | in-dist s1 | 0.013/0.015/0.014/0.014 | 0.019/0.024/0.023/0.022 | 0.003/0.003/0.004/0.004 | 3.0/3.4/3.1/3.0 |
| | held-out s0 | 0.014/0.015/0.015/0.013 | 0.019/0.021/0.020/0.021 | 0.003/0.004/0.004/0.004 | 3.1/3.2/3.4/3.6 |
| | held-out s1 | 0.016/0.017/0.016/0.014 | 0.022/0.023/0.022/0.019 | 0.003/0.003/0.004/0.005 | 3.6/3.7/3.3/3.3 |
| paired, inert (old) | held-out s0 | 0.012/0.014/0.014/0.013 | 0.017/0.020/0.022/0.018 | 0.003/0.003/0.004/0.004 | 3.6/3.3/3.9/4.2 |
| | held-out s1 | 0.014/0.016/0.017/0.019 | 0.024/0.025/0.028/0.023 | 0.004/0.004/0.004/0.007 | 5.1/4.4/4.4/4.9 |
| unpaired | held-out s0 | 0.015/0.016/0.013/0.015 | 0.022/0.021/0.021/0.021 | 0.003/0.003/0.003/0.004 | 3.4/3.1/3.4/4.0 |
| | held-out s1 | 0.017/0.016/0.017/0.015 | 0.033/0.031/0.025/0.021 | 0.004/0.004/0.004/0.004 | 4.1/3.9/3.6/3.2 |
| Gauss-Newton | held-out | 0.014/0.012/0.012/0.010 | 0.017/0.020/0.020/0.016 | 0.004/0.003/0.004/0.004 | |

Two-seed means (in-dist | held-out): fixed 0.013/0.015/0.013/0.014 | 0.015/0.016/0.015/0.013; inert 0.011/0.014/0.012/0.014 | 0.013/0.015/0.016/0.016; unpaired 0.013/0.016/0.013/0.013 | 0.016/0.016/0.015/0.015. Held-out Mahalanobis² (means over bins): fixed 3.3-3.5, inert 3.9-4.5, unpaired 3.5-3.7. Per-voxel median error, in-dist / held-out: fixed 0.0145 / 0.0166, inert 0.0143 / 0.0145, unpaired 0.0138 / 0.0163, GN 0.0127 / 0.0143. corr(voxel median error, r⊥) (s0, s1): fixed -0.49, -0.77 (mean -0.63); inert -0.70, -0.50; unpaired -0.75, -0.47; GN -0.63.

**Realism-trained nets** (T=3, lr 1e-4; held-out median angle pooled over the four bins, `[seed 0, seed 1]`; in-dist and the per-bin values are in `step4_pairfix_summary.txt`):

| test set | unpaired | paired, inert (old) | paired, fixed | Huber GN | plain GN |
|---|---|---|---|---|---|
| clean | 0.0193 [0.0176, 0.0210] | 0.0186 [0.0185, 0.0186] | 0.0199 [0.0214, 0.0183] | 0.0144 | 0.0119 |
| neighbours | 0.0719 [0.0723, 0.0715] | 0.0710 [0.0702, 0.0718] | 0.0706 [0.0739, 0.0673] | 0.2223 | 0.3428 |
| noise | 0.0242 [0.0227, 0.0258] | 0.0231 [0.0228, 0.0233] | 0.0257 [0.0262, 0.0253] | 0.0215 | 0.0573 |
| all | 0.0821 [0.0827, 0.0815] | 0.0762 [0.0734, 0.0789] | 0.0775 [0.0780, 0.0770] | 0.2239 | 0.3311 |

In-dist medians (fixed vs unpaired vs inert): clean 0.0188 / 0.0185 / 0.0185, neighbours 0.0585 / 0.0576 / 0.0594, noise 0.0250 / 0.0241 / 0.0233, all 0.0629 / 0.0615 / 0.0616. Held-out Mahalanobis² (fixed): clean 2.5, neighbours 2.9, noise 3.3, all 2.8. Held-out per-bin medians on `all`: fixed 0.068/0.081/0.082/0.079, inert 0.064/0.083/0.076/0.086, unpaired 0.064/0.084/0.073/0.086. corr(voxel median error, r⊥) (s0, s1): clean -0.58, -0.74; neighbours -0.12, -0.16; noise -0.57, -0.62; all -0.08, -0.10 (unpaired -0.17/-0.12 on `all`, inert -0.12/-0.12).

**Findings.**
(a) *Clean data: encoder-level pairing neither helps nor hurts.* Held-out/in-dist medians, z/⊥ RMS and Mahalanobis² of the fixed net lie inside the spread of the unpaired and inert-paired runs (differences ≤ 0.003°, comparable to the seed-to-seed spread of one configuration, e.g. unpaired held-out 0.0144 vs 0.0181 per-voxel). Held-out ⊥ RMS is 0.003-0.005° in every run, so the ray-direction difference between the two detectors does not improve ⊥ beyond what the two detectors' rows in the normal equations already provide; z RMS (0.019-0.024° held-out) is not better than GN's 0.016-0.020°. The fixed net's held-out Mahalanobis² (3.3-3.5) is slightly lower than the inert (3.9-4.5) and unpaired (3.5-3.7) ones, a small effect in two seeds. The gap to GN (held-out 1.0-1.5x) remains.
(b) *Realistic data: no separation.* Fixed pairing is within seed spread of the unpaired net on `neighbours` (0.0706 vs 0.0719), `noise` (slightly worse, 0.0257 vs 0.0242, 0.0231 inert) and `all` (0.0775 vs 0.0821). The earlier "suggestive" gain of the inert-paired net on held-out `all` (0.0762, 7 % below unpaired) is reproduced in size by the fixed net (0.0775, 6 % below) with both seeds (0.0780, 0.0770) below both unpaired seeds (0.0827, 0.0815) but also within 0.004° of the inert ones; this is as before suggestive, not established (2 seeds, 6 held-out voxels), and it is not larger with live mixing, so the partner-residual input and the encoder mixing give the same, small amount. The untested hypothesis stands: a neighbour spot lands consistently on both detectors and cannot be rejected by a pair-consistency check.
(c) *Error vs r⊥:* unchanged. Clean -0.63 (fixed) vs -0.60 (inert), -0.61 (unpaired), GN -0.63, with seed spread (-0.49 to -0.77) larger than the differences; under distractors the trend is still lost (-0.08 to -0.16 on `neighbours`/`all`).

**Caveats.** 2 seeds, 6 held-out voxels, one lr (1e-4) and T = 3 (lr 3e-4 and the gate variant were not run); differences of ~0.003° are inside the seed spread. The weights of the fixed runs are not saved, so that the pair MLP actually moved away from zero is inferred from the tests (gradient flows) and from the training curves differing from the inert runs, not measured directly. `gnpairv1`'s instability was not reproduced with the zero-initialised linear update at lr 1e-4. Conclusion for Step 3/4: the "no clear benefit of pairing" finding now holds for live encoder-level mixing as well, on this simulated geometry.

### Review fixes (2026-10-01)

Guards, tests and hygiene after the branch review; no behaviour change, the saved benchmark numbers are unaffected.
- `BatchedObserver` asserts `range_width > 0` and a contiguous, increasing `range_index` over valid omega bins (the Jacobian uses `|range_width|`, `nominal_offsets` the signed width and wedge index = bin order).
- Silent no-ops now assert: `--eval-variants neighbours/all` and `gauss_newton_baseline.py --realistic neighbours/all` without `dis_windows`; `--subpixel` without `--aux`.
- `summarize_arch_step4.py` keys the cached `dis_test_meta.npz` to the test file (path + size) and rebuilds on mismatch.
- `make_realistic_windows` / `make_realistic_dataset` take an optional `valid` mask (default None = apply it to padded entries too, as in all reported runs; the mask was not used in them).
- `scripts/make_dis_val.py` recreates `toy_orientation_arch_dis_val.pt` (checked sample-for-sample against the existing file).
- `ema_update` / `select_checkpoint_state` factored out of the trainer (`orientation_nn.py`); `inv3` clamps |det| to 1e-30; small script cleanups; `.gitignore` for `*_done.flag`.
- New tests (`TestReviewFixes`): EMA and checkpoint selection, distractor frame filter/coding, default ridge vs the undamped GN step, padded-entry realism and the mask, `inv3`.

### Review fixes (pre-merge)

Second review pass on PR #29; no saved benchmark number changed and nothing was rerun.
- Trainer: `--save-predictions` naming now goes through `variant_path` (stem + suffix, default `.npz`; previously a path without `.npz` made every realistic variant overwrite the clean predictions). New `--mask-padding` (default off, so reported runs stay reproducible) passes the valid-entry mask to the train, validation and test realism layer (`make_realistic_windows`). `--subpixel` with `--arch probe` asserts (it was silently ignored); the `--aux` table is loaded once; the results json is always `{variant: rows}` (`summarize_results.py` reads this and the older single-variant layout, and derives the |delta| bin keys from the file). "saved ..." messages print the path as given, not the absolute one.
- `GNLayerNet` warns when called without `aux["nom_off"]` (the window-centre fallback reintroduces the sub-pixel bias). `gauss_newton_baseline.py` checks the renderer before applying the realism layer to windows.
- `CentroidGaussNewton.solve` now uses `_linearize` (bit-identical results; `_linearize` returns the kept-spot mask too); `information()` is documented as J^T J without Huber weights, unlike the covariance returned by `solve`.
- Type annotations on the remaining untyped functions of the branch, `main()` guards in the summarize scripts, long lines wrapped, the Huber-choice predictions (`dis_pred_huber*.npz`, `val/val_*.npz`) committed so the threshold choice is reproducible from the repo. Summaries of the committed results (`summarize_results.py`, `summarize_diagnostics.py`, `summarize_arch_step4.py`) are identical before and after.
- Tests: `tests/test_toy_scripts.py` (prediction naming, padding mask, both json layouts, the GN warning).

- **Terminology rename (2026-10-03): "corruption" → "realistic".** "Corrupted data" suggests something is wrong; here it means deliberately simulated, more realistic windows (genuine spots of neighbouring voxels and a Σ3 twin, plus detector noise), with exact labels. Code and docs renamed; the old names remain as aliases; nothing was re-run and no committed result or log was edited (filenames with `corr`, e.g. `dis_res_gn_k3_corr_s*`, `dis_res_gnpairfix_k3_corr_*`, mean realism-trained; the run label `corr-trained=` in the `summarize_arch_step4.py` commands is kept so the committed `step4_pairfix_summary.*` reproduce byte for byte). Summaries of the committed data are identical before and after.

| Old | New |
|---|---|
| `CorruptionConfig` | `RealismConfig` |
| `corrupt_windows` | `make_realistic_windows` |
| `corrupt_dataset` | `make_realistic_dataset` |
| trainer `--corrupt-train` | `--realistic-train` (dest `realistic_train`) |
| `gauss_newton_baseline.py --corrupt` | `--realistic` |
| corruption-aware / corrupted-trained / corr-trained | realism-aware / realism-trained |
| `corr` in result filenames | realism-trained (files not renamed) |

Variant names (`clean`, `neighbours`, `noise`, `all`) and the `RealismConfig` fields (`p_miss`, `p_flip`, `p_hot`, `p_blob`) are unchanged. The old Python names are module-level aliases in `orientation_eval.py`, and the old flags are extra option strings of the same arguments (tested in `TestReviewFixes`).


## Toy Orientation NN — Perturbation sweep (2026-10-04)

**Question.** How does the network's error grow as the nominal orientation gets further from the truth? Randomly sampled voxels of the 500-grain sample are each tested on their own: the nominal is the voxel's true `.mic` orientation rotated by a random rotation of angle r, the net sees windows rendered exactly as in the dataset generator, and its error is the angle between its answer and the true orientation. No other baselines (the truth is known).

**Plan / what was built.**
1. `train_toy_orientation_nn.py --save-model PATH` (new): saves the final weights plus `model_kwargs` (constructor arguments), `no_frame`, `frame_half_width`, `window_size`, `realistic_train`, `mask_padding`, `seed`, and the trainer's argparse namespace. Weights were never saved before. Files go to `scripts/toy_orientation_sweep_model_{clean,realistic}_s{0,1}.pt` (gitignored by `scripts/*.pt`).
2. Retrained the four headline nets (unpaired `GNLayerNet`, T=3, lr 1e-4, 60 epochs, EMA 0.998, `--val-voxels 4`, seeds 0 and 1) with the recorded commands (design doc Section 9): clean-trained = `multi_res_gn_k3_lr1e-4` on `toy_orientation_stage3_multi_{train,test}.pt`; realism-trained = `dis_res_gn_k3_corr` (`--realistic-train all`, `toy_orientation_arch_dis_{train,test}.pt`); both with `--aux scripts/toy_orientation_stage3_multi_aux.pt`. Outputs went to a scratch directory, not over the committed results.
3. `scripts/perturbation_sweep.py` (new, tests in `tests/test_perturbation_sweep.py`).

**Reproduction check.** All four retrains reproduce the committed predictions bit-for-bit (`pred_deg` and `chol`, max abs diff 0.0, for the clean test set and, for the realism-trained nets, all four variants; `benchmarks/toy_orientation_sweep/retrain_reproduction.txt`). One caveat: two trainings run concurrently on the one MPS GPU diverged from the committed logs in the 6th digit by epoch 4 (GPU non-determinism under sharing), so they were killed and the four runs redone one at a time; run alone, they are exact.

**Protocol.**
- Voxels: 50 drawn uniformly at random (`--voxel-seed 0`, a seeded permutation of the eligible voxels, taking the first usable ones) from the 24,570 voxels, with r_perp <= 500 um, excluding the 30 voxels of `stage3_multi`/`arch_dis` (read from the dataset files; the same 30 in both), and requiring >= 40 ROI peaks at the true orientation plus a valid window spec. 20,601 eligible, none rejected. They span r_perp 57-499 um (terciles 57-246 / 249-371 / 381-499 um; 17/16/17 voxels), 46 distinct grains, 35 within 1.25 voxel pitch of another grain, 4 in a grain that also contains one of the 30 dataset voxels. Per-voxel index, position, r_perp, grain id, boundary flag and cross-grain-neighbour count are in the raw npz.
- Perturbation: for each (voxel, r) 20 directions with uniformly random axes; delta (length exactly r) and nominal R_nom = exp(-[delta]x) R_true, so exp([delta]x) R_nom = R_true, the dataset convention (`offsets_to_matrices`). The net predicts delta; error = angle of exp([delta_hat]x) R_nom R_true^T. Radii 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5 deg (10 radii x 20 directions x 50 voxels = 10,000 cases per window variant).
- Windows per case: its own nominal (so its own ROI set, `define_roi_set` + `|sin eta| >= 0.3`, observer, window spec, context, sub-pixel offsets), Q_max 8, both detectors, 32x32, K=4, alpha=0, rendered at the true orientation. `clean` = windows only. `all` = generator's distractor layer (2 nearest mic neighbours within 30 um + a Sigma3 twin of the nearest; each source present with probability 0.5 at an independent 0.3 deg rms offset) plus `make_realistic_dataset("all")` with a fixed seed per (voxel, radius), and the valid (padding) mask (the committed runs did not use it).
- Distractor sources: the sources sit at their own mic orientations (the target's truth is its mic orientation, so a same-grain neighbour is the target's own grain). In the generator the sources follow the target's delta from its nominal; here the target's "delta" is absorbed in the nominal, so each source is the neighbour's mic orientation plus the 0.3 deg rms Gaussian offset only.
- Failure rule: a case fails at pass 1 if the perturbed nominal has no ROI peak, fewer than 40 ROI entries, an invalid window spec, or fewer than `--min-present 20` windows with a lit pixel (counts in `all` include distractor/noise pixels). Failures are counted per radius/variant and excluded from the medians (not silently dropped: see the failure rows and `n_ok`). Padding to 130 entries.
- One-shot = pass 1. Iterated = 3 passes; each pass rebuilds ROI set, windows and context at the previous estimate, re-renders the SAME true orientation with the same distractor draws and the same realism seed (the noise realisation is the same draw but entry indices belong to the new ROI set, so it is not the same pixel noise on the same spot). A case that cannot continue (ROI/spec/too few present) keeps its last estimate; stops are counted.
- Metrics per (variant, model, r, pass): median and RMS angle, RMS z and perp (rotation-vector components of the error, sample frame), fraction < 0.1 deg, fraction improved (error < r), mean Mahalanobis^2, n_ok, and the same split by r_perp tercile.

**Commands.**
```
cd icenine_py
# retrain (one at a time; add --save-model), e.g. realism-trained seed 0:
uv run python scripts/train_toy_orientation_nn.py --head offset --arch gn --gn-iters 3 --loss decoupled --device mps \
  --lr 1e-4 --clip 1 --cosine --batch-size 64 --epochs 60 --ema 0.998 --checkpoint ema --val-voxels 4 --realistic-train all \
  --eval-variants clean,neighbours,noise,all --seed 0 --train scripts/toy_orientation_arch_dis_train.pt \
  --test scripts/toy_orientation_arch_dis_test.pt --aux scripts/toy_orientation_stage3_multi_aux.pt \
  --save-model scripts/toy_orientation_sweep_model_realistic_s0.pt --save-predictions <scratch>/real_s0.npz
# sweep + summary + plot (10 workers, CPU, 1340 s wall time)
M=scripts/toy_orientation_sweep_model
uv run python scripts/perturbation_sweep.py run --workers 10 --out-dir benchmarks/toy_orientation_sweep \
  --models clean_s0=${M}_clean_s0.pt clean_s1=${M}_clean_s1.pt realistic_s0=${M}_realistic_s0.pt realistic_s1=${M}_realistic_s1.pt
uv run python scripts/perturbation_sweep.py summarize --out-dir benchmarks/toy_orientation_sweep   # re-make txt/json/png
```
Outputs in `benchmarks/toy_orientation_sweep/`: `perturbation_sweep_raw.npz` (per-case arrays `[voxel, radius, direction, variant, model, pass]`), `perturbation_sweep_summary.{txt,json}` (all metrics, per seed, per tercile, failure and stop counts), `perturbation_sweep.png`, `retrain_reproduction.txt`, `run.log`.

**Results** (median angular error in degrees, mean of the two seeds; 1000 cases per radius, 940 at r = 5 deg in clean windows; one-shot -> iterated x3; per-seed values and everything else in the summary files):

| r (deg) | clean-trained, clean windows | realism-trained, clean windows | clean-trained, `all` windows | realism-trained, `all` windows |
|---|---|---|---|---|
| 0.05 | .0125 -> .0124 | .0122 -> .0121 | .274 -> .278 | .062 -> .065 |
| 0.1 | .0125 -> .0124 | .0131 -> .0121 | .259 -> .263 | .059 -> .063 |
| 0.25 | .0123 -> .0124 | .0137 -> .0121 | .243 -> .251 | .064 -> .064 |
| 0.5 | .0116 -> .0124 | .0156 -> .0121 | .253 -> .264 | .063 -> .065 |
| 0.75 | .0122 -> .0124 | .0177 -> .0121 | .245 -> .261 | .067 -> .065 |
| 1 | .0130 -> .0124 | .0208 -> .0121 | .241 -> .266 | .071 -> .066 |
| 1.5 | .0194 -> .0125 | .0281 -> .0122 | .277 -> .272 | .088 -> .061 |
| 2 | .0303 -> .0125 | .0382 -> .0121 | .413 -> .280 | .171 -> .072 |
| 3 | .0606 -> .0125 | .0662 -> .0121 | 1.11 -> .291 | 1.01 -> .079 |
| 5 | .168 -> .0125 | .273 -> .0120 | 3.85 -> .800 | 4.24 -> 1.97 |

Failures at pass 1 (of 1000 cases): clean windows 0 at every radius except r = 5 deg, 60 (all "fewer than 20 present windows"); `all` windows 0 except 6 at r = 5 deg. Early stops in later passes: none in clean windows; in `all` windows 6 / 1 / 1 / 5 case-stops for clean_s0 / clean_s1 / realistic_s0 / realistic_s1 (all at r = 5 deg except one spec failure at 0.5 deg).

Reading (what the data support, with 50 voxels and 2 seeds):
- Inside the training ball (r <= 1 deg) the one-shot errors on clean windows are 0.012-0.021 deg and the realism-trained error rises slowly with r (0.012 -> 0.021 deg), while the clean-trained one is flat (0.012-0.013). Beyond 1 deg the one-shot error grows roughly like r^1.5-2 for both (about 0.03 at 2 deg, 0.06 at 3 deg, 0.17-0.27 at 5 deg), i.e. it stays well below r (error < r in 100% of clean-window cases) but loses sub-0.1 deg accuracy (fraction < 0.1 deg: 0.80 at 3 deg, 0.10-0.20 at 5 deg).
- Re-centring fixes that: after 2-3 passes the clean-window error is at its floor (~0.012 deg, independent of r, all cases < 0.1 deg for every r including 5 deg) because the second pass starts within the training ball. The floor is the same as the r = 0.05 deg error. 5 deg is extrapolation (trained in a 1 deg ball), and on clean windows that extrapolated first step is still good enough to bring the estimate inside the ball.
- With realistic (`all`) windows the error floor is set by the distractors, not by r: clean-trained ~0.25 deg (mean Mahalanobis^2 ~1,000-2,000, i.e. confidently wrong), realism-trained ~0.06-0.07 deg (Mahalanobis^2 ~3). These floors agree with the committed held-out numbers (`.2033`/`.0615`, 30-voxel test). The realism-trained net stays at its floor up to r = 1.5 deg one-shot (0.088), rises to 0.17 at 2 deg and fails at >= 3 deg one-shot (1.0 deg at 3; 4.2 at 5, about the same as no correction); iterating recovers it to 0.079 at 3 deg but not at 5 deg (median 1.97 after 3 passes, 24% of cases < 0.1 deg; wide spread between seeds: 1.30 / 2.63). For the clean-trained net the ~0.25 floor is not improved by iterating.
- Iterating does not help inside 1 deg on realistic windows (differences <= 0.003 deg at r <= 0.75; clean-trained: slightly worse, 0.241 -> 0.266 at 1 deg).
- r_perp: on clean windows the error falls with r_perp (floor 0.023 / 0.012 / 0.010 deg for the three terciles, the parallax effect again); on realistic windows no clear trend. Voxels near a grain boundary (35 of 50) are not worse than the others (r = 1: 0.065 vs 0.09 deg for the realism-trained net); with 15 interior voxels this is not a finding.

**Caveats.** (1) Trained within a 1 deg ball, so r > 1 deg is extrapolation. (2) Single-voxel windows, not a full-sample render: neighbours/twin come from the generator's source model only (2 neighbours + twin, same-family nuisances as in training), not from a rendering of everything that actually overlaps. (3) 50 voxels x 20 directions per radius, 2 seeds per model type; standard errors/intervals were not computed, voxel-to-voxel variation is not separated from direction-to-direction. (4) The failure rule counts windows with any lit pixel, so noise and distractors can mask a low peak count in `all` windows (6 vs 60 failures at 5 deg). (5) Iteration reuses distractor draws and the realism seed but pixel noise is attached to entry indices of a rebuilt ROI set. (6) Only the median is plotted; RMS, z/perp components and the fractions are in the summary. (7) The padding mask was applied here but not in the committed runs, so the `all` numbers are not pixel-identical to the committed test protocol, though the floors agree. (8) The perturbed nominal's ROI set and window centres differ from the dataset's (they are defined at the perturbed nominal, as a real hand-off would), a slight protocol difference from the datasets, where the ROI set is defined at the truth.

### Comparison with existing optimizers (2026-10-05)

**Question.** On exactly the sweep's cases, how do the project's non-learned optimizers compare with the network? Cases are reproduced, not re-drawn: `perturbation_sweep.case_draws` (factored out of `sweep_voxel`) gives the same seeded directions and distractor draws, the voxels come from the sweep's raw npz, and every task asserts that its per-case ROI counts and pass-1 failure codes (both variants, hence the realism windows' present counts) equal the sweep's. Error = angle of R_est R_true^T, the sweep's metric (no crystal-symmetry reduction; irrelevant for r <= 5 deg). 50 voxels x 10 radii x 20 directions = 1000 cases per radius and variant, 500 tasks, all run for every method (no reduction was needed: the projection from the 1-voxel pilot was 2 h; the full run took 1 h 49 min on 10 workers). Summaries use the cases where the net ran at pass 1 (n = 1000 per radius, 940 at r = 5 clean, 994 at r = 5 realistic), the same set for every column.

**Methods** (`scripts/optimizer_sweep.py`, existing code and settings):
- **MC**: `MCOptimizer` on the binary pixel-overlap cost via `bench_hp_sweep.run_one_mc`: 3500 steps, 2 restarts, step fraction 0.5, search box 1.5 r. **It is told r through the box (an advantage the net does not have).** Note that `MCOptimizer` restarts after `2 (box/step)^3` steps without improvement and stops after the second restart, so when it sees no signal it ends after ~50 steps; the seed in `optimizer_baselines.py` has no effect (unseeded `default_rng`), here the generator is made deterministic per case.
- **Adam**: `run_one_riemannian_adam_geoopt` on the differentiable cost, scale 2 (8x max-pool), omega window 1, lr 1e-4, 100 steps. Not told r.
- **GN / Huber GN (c = 1)**: `CentroidGaussNewton` on the net's windows, one-shot at the perturbed nominal and iterated x3 with the sweep's re-centring (same pass loop, same re-rendering and gating). Runtime counts measurement extraction + solve, not window rendering.
- All: |g| <= 8, both detectors, the net's eligibility. The cost functions gained an optional `min_sin_eta` (default 0 = unchanged) so `|sin eta| >= 0.3` and the config's EtaLimit (86 deg, not the cost default 90) apply to MC and Adam as to the net's ROI set.

**What each method sees.**
- Net / GN: 32x32 windows around the ROI spots of the case's perturbed nominal, +-4 frames, frame-coded; distractors and noise exist only inside windows.
- MC / Adam (new, no full-sample image exists): per-case detector images made from pixel sets, the lit pixels the simulator's rasteriser would give (`lit_pixel_set`). `clean` = every spot of the target voxel at its TRUE orientation (|g| <= 8, all frames, both detectors; no |sin eta| cut in the image, the cost applies it). `all` = that plus (i) every active distractor source (2 neighbours + the Sigma3 twin, the sweep's draws) over the WHOLE detector, plus (ii) the sweep's realism edits transplanted pixel by pixel: the difference between the net's windows before and after `make_realistic_dataset` (pixels dropped or grown, hot pixels, blobs; a window zeroed as a missing spot removes all its pixels, target and distractor) is applied at the window's detector position and frame. The edits are made on the pass-1 windows, so they follow that nominal's window grid. Checks (`tests/test_optimizer_sweep.py`, probe on one voxel): the target-only image gives hard-cost quality exactly 1.0 at the truth (54 peaks); every pixel of the net's realistic windows is in the realistic image and every removed pixel was present before the edit; the coarse image stack equals `MultiScaleImageStack` (scale, omega window 1) exactly. So the optimizers see the whole detector (all of the voxel's spots, more than the net's ROI windows) but the same random nuisances inside the windows. Differences: pixel noise outside the net's windows does not exist (there is none to inject: the noise is defined per window), noise pixels that fall off the detector or outside the frame range are dropped, spots whose frame has left the +-4 window (large r) are still in the image, and the image holds pixel sets of the simulator's own spots with no intensity/blur.
- Cost at the true orientation on the exact (clean) images (stored per case, `q_true`): the hard cost's quality is 1.0 for the first checked voxel, but its median over the 50 voxels is 0.93 (per-voxel medians 0.90-0.97, 1.0 for only that one voxel). The images are exact there (checked on voxels 3108 and 2910: pixel overlap = pixels on the detector and every peak overlaps), so I attribute the ceiling to the cost's own per-peak detector factor (quality_i = pixel ratio x detectors-overlapping / n_det in `OverlapInfo.update_quality`, which is below 1 for spots recorded on one detector), not to the images; I did not trace it further. On `all` images the median is 0.76; Adam's soft scale-2 cost is 0.66 on both.

**Results**: median angular error (deg), mean of the two net seeds for the net columns (the `clean-trained` net is in the summary files). `x3` = iterated.

Clean data:

| r (deg) | MC | Adam | GN | GN x3 | Huber | Huber x3 | net (realism-trained) one-shot | net x3 |
|---|---|---|---|---|---|---|---|---|
| 0.05 | .0274 | .134 | .0130 | .0127 | .0123 | .0126 | .0122 | .0121 |
| 0.1 | .0573 | .133 | .0128 | .0127 | .0123 | .0126 | .0131 | .0121 |
| 0.25 | .136 | .119 | .0122 | .0127 | .0120 | .0126 | .0137 | .0121 |
| 0.5 | .288 | .145 | .0127 | .0127 | .0118 | .0126 | .0156 | .0121 |
| 0.75 | .424 | .235 | .0130 | .0127 | .0123 | .0126 | .0177 | .0121 |
| 1 | .608 | .344 | .0135 | .0127 | .0122 | .0126 | .0208 | .0121 |
| 1.5 | 1.11 | .575 | .0217 | .0127 | .0149 | .0126 | .0281 | .0122 |
| 2 | 1.64 | .875 | .0328 | .0127 | .0202 | .0126 | .0382 | .0121 |
| 3 | 2.26 | 2.17 | .0641 | .0127 | .0406 | .0126 | .0662 | .0121 |
| 5 | 4.50 | 4.80 | .127 | .0122 | .0868 | .0123 | .273 | .0120 |

Realistic (`all`) data:

| r (deg) | MC | Adam | GN | GN x3 | Huber | Huber x3 | net (realism-trained) one-shot | net x3 |
|---|---|---|---|---|---|---|---|---|
| 0.05 | .0282 | .269 | .274 | .284 | .196 | .196 | .0622 | .0652 |
| 0.1 | .0566 | .246 | .249 | .263 | .171 | .172 | .0589 | .0626 |
| 0.25 | .137 | .201 | .249 | .258 | .173 | .177 | .0637 | .0644 |
| 0.5 | .287 | .187 | .261 | .271 | .175 | .177 | .0632 | .0654 |
| 0.75 | .422 | .276 | .256 | .266 | .184 | .190 | .0665 | .0646 |
| 1 | .576 | .373 | .258 | .280 | .179 | .186 | .0707 | .0664 |
| 1.5 | 1.05 | .613 | .306 | .289 | .182 | .188 | .0880 | .0613 |
| 2 | 1.65 | .911 | .464 | .290 | .187 | .206 | .171 | .0720 |
| 3 | 2.45 | 2.20 | 1.11 | .294 | .278 | .218 | 1.01 | .0789 |
| 5 | 4.38 | 4.83 | 3.80 | .876 | 2.64 | .232 | 4.24 | 1.97 |

RMS, fraction < 0.1 deg, fraction improved, runtimes, the clean-trained net, and every paired comparison are in `benchmarks/toy_orientation_sweep/optimizer_sweep_summary.{txt,json}`; plot `perturbation_sweep_vs_optimizers.png`; per-case raw arrays `[voxel, radius, direction, variant, method, pass]` in `optimizer_sweep_raw.npz` (3 MB, committed), log `optimizer_run.log`.

Paired comparisons (identical cases, realism-trained net iterated x3 has the smaller error; seed mean, fraction of 1000 cases):

| r (deg) | clean: vs MC | vs Adam | vs GN x3 | vs Huber x3 | realistic: vs MC | vs Adam | vs GN x3 | vs Huber x3 |
|---|---|---|---|---|---|---|---|---|
| 0.05 | .76 | .99 | .48 | .49 | .26 | .95 | .95 | .82 |
| 0.5 | .99 | 1.00 | .48 | .49 | .91 | .90 | .96 | .80 |
| 1 | .99 | 1.00 | .48 | .48 | .96 | .96 | .94 | .80 |
| 2 | .99 | 1.00 | .49 | .49 | .97 | .99 | .93 | .79 |
| 3 | 1.00 | 1.00 | .49 | .49 | .94 | .94 | .85 | .72 |
| 5 | .99 | 1.00 | .48 | .49 | .70 | .81 | .43 | .26 |

Median runtime per case (s, data preparation excluded): MC 2.4, Adam 0.30, GN x3 0.02-0.03, Huber x3 0.04-0.13; the net was not timed here (one batched forward pass per case, three for x3).

Reading (50 voxels x 20 directions, 2 net seeds; no confidence intervals):
- Clean data: GN and Huber GN reach the quantisation floor (0.012-0.013 deg) from every start up to 5 deg once iterated (one-shot degrades like the net's: 0.13 deg at 5 deg), and the net iterated x3 sits on the same floor; the paired win rate against GN x3 is 0.48 (a coin flip at the floor). So on clean data the net offers no accuracy gain over the existing centroid GN, only (potentially) speed/robustness elsewhere. MC and Adam, on the sharp pixel-overlap/coarse costs, do not converge to the floor: MC's median error is about 0.5-0.6 r (it moves toward the truth in 77-89% of cases but stops short; < 0.1 deg only for r <= 0.1), Adam's is 0.13 deg at r = 0.05 and 0.34 at 1 deg.
- Realistic data: plain GN is dominated by the distractors (0.26-0.29 deg floor, like the clean-trained net; paired win rate of the realism-trained net 0.93-0.96 for r <= 3); Huber GN halves it (0.17-0.22) but stays about 3x above the realism-trained net (0.065); the net wins 80% of paired cases up to r = 2. At r = 5 deg the iterated Huber GN (0.23) clearly beats the iterated net (1.97): the net wins only 26% of paired cases there. MC (told r) is the only method that beats the net at the smallest radii (r = 0.05: 0.028 vs 0.065, net wins 26% of cases; r = 0.1: equal, 45%), because its error scales with r while the net's floor is constant; from r = 0.25 the net is better (76% -> 96%).
- Adam (HP-sweep-best setting) is worse than the net everywhere (net wins >= 80% of cases at every radius, both variants). Its lr 1e-4 x 100 steps limits the travel to about 0.6 deg, so it cannot correct r >~ 1 deg (errors 2-5 deg at r = 3-5), and for r <= 0.1 it moves away from the truth (error 0.13-0.27 deg > r: the scale-2 cost's optimum is offset from the truth, its value at the truth being 0.66).

**Caveats.** (1) MC is told r. (2) The optimizers' images are an approximation of "the same data": whole-detector distractors, target spots at |g| <= 8 with no cut, window-level realism edits transplanted, but the cost sees all eligible peaks at the candidate orientation, not the net's fixed ROI set; the costs at the truth are below 1 even on exact images (see above), so their ceiling is not 1 and, on `all` images, the optimum need not sit exactly at the truth. (3) The Adam cost has a systematic coordinate offset at scale 2 (`differentiable_cost.py` comments), which sets its small-r error. (4) The settings are the HP-sweep protocol's, not tuned per radius: Adam with a larger lr or more steps, or MC with more steps, would do better at large r; this compares the existing protocol, not the best achievable. (5) GN runtimes exclude window rendering; the net's were not measured. (6) Same single-voxel caveats as the sweep; sources use the generator's model.

**Caveat on identical noise (added after review).** The optimizers' realistic images carry only the pass-1 realism edits. The net's passes 2-3 re-render the windows, and their pixel noise attaches to the entry indices of the rebuilt ROI set, so "same nuisances" holds strictly for the one-shot comparison only; the iterated-net comparison on realistic data is therefore not on identical noise. The one-shot paired numbers in the summaries are strictly like-for-like. The committed raw npz files and logs contain absolute local paths (not secrets).


### Comparison with multi-level reconstruction (FindOptimal) (2026-10-05)

**Question.** How does the project's multi-level adaptive reconstruction (`AdaptiveVoxelReconstructor`: coarse levels -> FindOptimal -> VarianceMinimizing -> final overlap) do on the sweep's single-voxel cases? Two experiments, `scripts/findoptimal_sweep.py` (summary and plot: `findoptimal_sweep_summary.py`), on the same 50 voxels and the same per-voxel clean and realistic (`all`) detector images the optimizer sweep builds (its pixel-set helpers; nothing simulates the whole sample; each task asserts the case is aligned with the network sweep).

**Refactor (behaviour unchanged).** The part of `reconstruct_voxel` after the coarse levels is now `AdaptiveVoxelReconstructor.refine_from_candidates(candidates, voxel_vertices, phase_index, diameter, ...)`, which `reconstruct_voxel` itself calls with the objects it built (same cost function, same MC optimizer and random stream, same eval counts); helpers `_make_local_cost_fn` / `_make_optimizers` hold the construction code that moved. `AdaptiveVoxelReconstructor(setup, min_sin_eta=0.0)` passes the new cost option to every cost function it builds (default off). Diagnostics `last_level_best` and `last_find_optimal` are recorded. Proof: `tests/test_findoptimal_refactor.py::test_reconstruct_voxel_unchanged_by_refactor` compares a deterministic `reconstruct_voxel` run (ThreeVoxels voxel 0, 3-orientation FZ set, 1 level, 150 MC steps, seed 7) with values recorded from the pre-refactor code (orientation error vector to 1e-9 and cost); `test_refine_from_candidates_standalone_converges` runs the new method on its own. The pre-refactor code was run twice and gave identical numbers (deterministic).

**Settings (both experiments).** `SearchParameters.from_config` of `Examples/Example2.ThreeVoxels/ConfigFiles/ReconstructQ8.config` (there is no ManyGrains reconstruction config; the geometry/detector files of the two examples are identical): grid radius 5 deg, `MinLocalResolution 0` / `MaxLocalResolution 3` (4 levels, diameter 5 -> 0.99 deg, divided by 1.5 per level), `MaxMCSteps 200`, `MCRadiusScaleFactor 0.4`, `SuccessiveRestarts 2`, `MaxConvergenceCost 1e-4`, `MaxDiscreteCandidates 30`, hard (binary pixel) cost, pixel radius 3 for the coarse screening and 0 afterwards, the level Q_max 5, 6, 7, 8 (`nQMax = 5 + level`), the local cost function with the simulator's reflection list at Q_max 8 (config MaxQ overridden to 8 as in the sweep), the ManyGrains physics, `MyFZ.dat`, `|sin eta| >= 0.3` through the new `min_sin_eta` (the reconstructor passes it through, so it is applied exactly as for the net, MC and Adam), EtaLimit 86 deg. No hybrid/Adam optimizer (FindOptimal is full MC).

**Experiment A (global, no starting guess, so no r).** `reconstruct_voxel` from scratch, once per voxel and variant (50 x 2 = 100 runs), rng seeded per voxel. The image of a voxel is the one the optimizer sweep gave its case at the smallest radius, direction 0 (r = 0.05 deg: the realism windows sit at the truth). It sees the whole detector (every spot of the target with |g| <= 8 [+ distractors and the realism edits in `all`]), not the net's ROI windows.

**Experiment B (local, depends on r).** `refine_from_candidates` with the case's perturbed nominal as the single candidate: FindOptimal (full MC, convergence check) + VarianceMinimizing + final evaluation, as `reconstruct_voxel` runs them. **Search size: reconstruct_voxel's FindOptimal does not use r.** Its box is `max(d/3, 0.2 deg) / 2^MinLocalResolution` with d the diameter after the levels (5 deg / 1.5^4 = 0.99 deg), i.e. 0.329 deg, MC step 0.4 x 0.329 = 0.132 deg, 200 steps, 2 restarts, with a convergence stop at `hit_ratio >= 1` (never reached here: 0/100 in A, <= 2.5% in B); the same box for r = 0.05 and r = 5 deg. B was run on all 20 directions x 10 radii x 50 voxels x 2 variants = 20,000 cases (the timing pilot, 2 voxels x radii 0.05 and 5 x 20 directions + A on 2 voxels, projected 1.05 h on 10 workers, far below the 4 h limit; the actual run took 5 min (A) + 32 min (B)); images exactly as the optimizer sweep builds them per case (seed per case); error summaries on the cases where the net could run at pass 1 (n = 1000 per radius and variant, 940 clean and 994 realistic at r = 5 deg), the same set as the other columns. A start that FindOptimal cannot improve is returned as is only if some candidate scores cost < 1; when MC finds no overlap at all (cost >= 1) `reconstruct_voxel`'s code returns the identity matrix (its initial best candidate), and that is what is scored (see the fallback table below).

**Error metric.** Misorientation angle to the true `.mic` orientation reduced by cubic symmetry: min over the 24 proper cubic operators S of angle(R_est S R_true^T) (`findoptimal_sweep.misorientation_deg`, operators from `icenine.symmetry.create_cubic_symmetry`; tested against `sampling.get_misorientation` and on a symmetry-equivalent rotation in `tests/test_findoptimal_sweep.py`). The sweep's other methods report the plain angle of R_est R_true^T; at r <= 5 deg their errors are far below 45 deg, and no cubic operator other than the identity is closer than 90 deg, so reduction cannot change them. Checked on 20,000 random cases per method (R_est rebuilt from the stored error vector and the voxel's mic orientation): max |reduced - plain| = 0 and no case lowered, for all four network models, MC, Adam, GN and Huber GN (largest plain error anywhere: 25.6 deg (MC), network 14.2 deg, Adam/GN/Huber 6.0-6.6 deg). For FindOptimal B likewise max |reduced - plain| = 0 at every radius; in A they differ for 16/50 (clean) and 13/50 (realistic) voxels (an equivalent orientation was returned).

**Results A (50 voxels per variant).**

| | clean | realistic |
|---|---|---|
| median / RMS / max error (deg) | 0.069 / 35.1 / 60.0 | 0.052 / 31.1 / 60.0 |
| fraction < 0.1 deg | 0.52 | 0.60 |
| **wrong (> 1 deg)** | **21/50** | **17/50** |
| median / RMS / max over the voxels it got right | 0.028 / 0.050 / 0.150 (n = 29) | 0.025 / 0.047 / 0.114 (n = 33) |
| runtime per voxel (median) | 26.8 s | 29.2 s |
| cost evals (median): global, pixel radius 3 / local | 44,474 / 4,656 | 44,748 / 7,790 |
| FindOptimal winner = best coarse candidate | 46/50 | 46/50 |

For comparison the net iterated x3 has a median 0.0121 (clean) and 0.062-0.079 (realistic, r <= 3 deg) from a start that is already within r; MC 0.027 at r = 0.05; GN x3 0.0127 / 0.26-0.29. A needs no starting guess, but its right answers (0.025-0.03 deg median) are 2x worse than the net's clean floor, and about 40% of the voxels are wrong: the median over all voxels is not a meaningful accuracy number here, which is why the plot also shows the median over the correct ones (dashed).

All 38 wrong solutions (clean 21, realistic 17) are coincidence-site-lattice (CSL) neighbours of the truth or close to one. Method (`findoptimal_sweep_summary.classify_csl`): crystal-frame misorientation M = R_true^T R_final S over the 24 cubic operators, smallest angle, axis folded into the cubic fundamental family; the lowest-Sigma CSL (odd Sigma <= 29, ideal angle/axis from tan(theta/2) = sqrt(N)/m) whose deviation (smallest angle between M and the ideal rotation over the cubic group on both sides) is within the Brandon limit 15 deg / sqrt(Sigma) is assigned. Counts (clean + realistic): **Sigma3 25** (60 deg about <111>, axes within 0.3 deg of <111>, deviation from the ideal 60 deg <= 0.41 deg), **Sigma7 5** (38.2 deg about <111>), **Sigma5 4** (36.9 deg about <100>), **Sigma11 1** (voxel 2489 clean, 50.45 deg about <110>), **Sigma17 1** (voxel 970 realistic: 29.6 deg, axis 3.8 deg from <100>, Sigma17a = 28.07 deg; deviation 2.45 deg against a Brandon limit of 3.64 deg: borderline), and **2 matching no CSL with Sigma <= 29** (voxel 16182 realistic: 37.9 deg, axis 25.5 deg from <111>; voxel 1053 realistic: 48.5 deg, axis 17.8 deg from <110>). No Sigma13. So the dominant failure is the Sigma3 twin (as in `docs/todo_findoptimal_wrong_candidates.md`), and the other relations are the low-Sigma CSLs, which share a subset of reflections with the truth: consistent with the known wrong-candidate issue. **In every wrong case the cost at the true orientation is lower than at the returned answer** (21/21 and 17/17: 0.0-0.09 vs 0.64-0.92 clean, 0.18-0.35 vs 0.42-0.86 realistic), so the failure is in the search (coarse levels keep the CSL orientation, FindOptimal polishes it), not in a cost that prefers the wrong orientation. In 3 of the 38 (voxels 16182 and 121 clean, 24474 realistic) an earlier level's best candidate was within 10 deg of the truth (4.8-8.2 deg) and the last level's was not, i.e. the coarse ranking dropped the right basin; in one more (voxel 2489 clean) the last level's best candidate was 8.7 deg off but FindOptimal returned another candidate (the Sigma11 one); in the other 34 no level's best candidate was within 10 deg of the truth. FindOptimal's winner is the best coarse candidate (index 0) in 46/50 runs per variant (clean: 17/21 wrong and 29/29 right; realistic: 17/17 wrong and 29/33 right): it rarely overrides the coarse ranking.

**Results B** (`findoptimal_sweep_summary.txt` has the RMS, fractions improved, per-radius details, one-shot net and clean-trained columns). Median error (deg); net = realism-trained, mean of 2 seeds, iterated x3; MC is told r; n = 1000 per cell (940 at r = 5 clean, 994 realistic).

| r (deg) | clean: FindOptimal | MC | GN x3 | net x3 | realistic: FindOptimal | MC | GN x3 | net x3 |
|---|---|---|---|---|---|---|---|---|
| 0.05 | .0197 | .0274 | .0127 | .0121 | .0189 | .0282 | .284 | .0652 |
| 0.1 | .0257 | .0573 | .0127 | .0121 | .0248 | .0566 | .263 | .0626 |
| 0.25 | .0426 | .136 | .0127 | .0121 | .0370 | .137 | .258 | .0644 |
| 0.5 | .0557 | .288 | .0127 | .0121 | .0552 | .287 | .271 | .0654 |
| 0.75 | .0748 | .424 | .0127 | .0121 | .0739 | .422 | .266 | .0646 |
| 1 | .0892 | .608 | .0127 | .0121 | .0905 | .576 | .280 | .0664 |
| 1.5 | .234 | 1.11 | .0127 | .0122 | .223 | 1.05 | .289 | .0613 |
| 2 | 1.18 | 1.64 | .0127 | .0121 | .869 | 1.65 | .290 | .0720 |
| 3 | 3.07 | 2.26 | .0127 | .0121 | 2.98 | 2.45 | .294 | .0789 |
| 5 | 23.9 | 4.50 | .0122 | .0120 | 6.02 | 4.38 | .876 | 1.97 |

Fraction of cases < 0.1 deg (FindOptimal; clean / realistic): .985/.990 (0.05), .958/.961 (0.1), .861/.862 (0.25), .707/.720 (0.5), .584/.570 (0.75), .522/.533 (1), .339/.339 (1.5), .164/.172 (2), .008/.012 (3), 0/0 (5). Median runtime of `refine_from_candidates` per case (s, clean; realistic similar): 1.64 (0.05), 1.56, 1.23, 0.73, 0.53, 0.49 (1), 0.37, 0.33, 0.30, 0.25 (5), against 2.4 for the MC protocol of the optimizer sweep; median cost evaluations 2,300 down to 320. With the fixed 0.329 deg box the number of evaluations falls as r grows (presumably because the search stops sooner when it finds no improvement; not examined).

Paired: fraction of identical cases where the realism-trained net (iterated x3) has the smaller error than FindOptimal (seed mean; one-shot in the summary):

| r (deg) | 0.05 | 0.1 | 0.25 | 0.5 | 0.75 | 1 | 1.5 | 2 | 3 | 5 |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | .63 | .72 | .81 | .86 | .88 | .89 | .93 | .97 | 1.00 | 1.00 |
| realistic | .24 | .29 | .38 | .49 | .58 | .62 | .75 | .86 | .95 | .91 |

Fallback when MC finds no overlap (cost >= 1, identity returned): the fraction of such cases is 0 for r <= 1 deg, then 0.002 / 0.033 / 0.193 / 0.532 (clean) and 0 / 0.015 / 0.129 / 0.427 (realistic) at r = 1.5 / 2 / 3 / 5 deg. If those results are replaced by the unchanged input (error = r) the medians are 3.00 / 5.00 (clean, r = 3 / 5) and 2.98 / 5.00 (realistic) and the RMS 6.9 / 11.1 and 6.5 / 11.7; the tables above are the reconstructor's actual output.

Reading (50 voxels x 20 directions per radius, one realisation of the images; no confidence intervals):
- Used as a local refiner (B), FindOptimal's error grows slowly with r (clean: 0.020 at r = 0.05, 0.043 at 0.25, 0.089 at 1 deg): its 0.329 deg box cannot undo more than about 1-1.5 deg; at r >= 2 deg median errors are 0.9-1.2, then 3 (about the start) and, at r = 5, the search often finds nothing (43-53% identity). Within r <= 1 deg 96-100% of cases improve (error < r, clean, r >= 0.1; 86% at 0.05), but it does not reach the net's/GN's floor (0.012 clean): fraction < 0.1 deg falls from 0.99 (r = 0.05) to 0.52 (r = 1).
- On clean data the net iterated x3 is better than FindOptimal in 63% (r = 0.05) to 100% of cases (r >= 3), GN x3 sits at the same floor as the net. On realistic data the picture is mixed at small r: FindOptimal's median error (0.019-0.055 deg for r <= 0.5) is below the realism-trained net's 0.063-0.065, and the net wins only 24-49% of paired cases up to r = 0.5, 58-62% at 0.75-1 deg, 75-86% at 1.5-2 and >= 90% at >= 3 deg. The same pattern as MC's, presumably for the same reason: its error scales with the start's distance while the net has a distractor-set floor. It is also much better than GN x3 / Huber x3 on realistic data below r ~ 1 deg (0.02-0.09 vs 0.17-0.28); I did not test why (it uses the whole-detector pixel overlap rather than windows).
- A (global) and B are different things: B's start is within r of the truth; A's has to find the basin. A's right answers (n = 29 clean, 33 realistic) have a median of 0.025-0.028 deg, about B's error at r ~ 0.1 deg (0.025-0.026), but 34-42% of the voxels are lost to CSL neighbours that the coarse search prefers.

**Caveats.** (1) Images are the optimizer sweep's approximation of the data (whole-detector pixel sets of the simulator's spots, window-level realism edits transplanted); the reconstructor sees all eligible peaks at the candidate orientation, the net only its ROI windows. (2) A uses one image per voxel (case r = 0.05, direction 0) and one seed per voxel; 50 voxels, so 21/50 vs 17/50 wrong is not a difference between variants (a different realisation could change several voxels; 11 voxels are wrong in both variants: 7336, 16182, 20206, 3142, 13884, 2489, 9484, 970, 9421, 24474, 17771). The Q_max 8 / `min_sin_eta` choices follow the net; the reconstruction's own default (no filter, Q_max from the config) was not run. (3) The search settings are ReconstructQ8's (200 MC steps, not the HP-sweep 3500), and FindOptimal here has no knowledge of r; a start-aware box would do better at r > 1 deg (not tested; MC with a 1.5 r box in the optimizer sweep does not either). (4) The identity return when MC finds nothing is the reconstructor's behaviour, kept; see the fallback numbers. (5) Runtimes are CPU wall time per case on 10 single-thread workers, not benchmarked in isolation; pilot per-case mean 0.97 s, full-run medians in the table. (6) The Brandon-criterion assignment is by lowest Sigma within tolerance and uses Sigma <= 29: at the borderline (voxel 970) or for the 2 unmatched cases a different Sigma list could change the label.

**Commands.**
```
cd icenine_py
uv run python scripts/findoptimal_sweep.py pilot --workers 10                  # timing pilot, log in findoptimal_pilot.log
nohup sh -c "uv run python scripts/findoptimal_sweep.py run-a --workers 10 && uv run python scripts/findoptimal_sweep.py run-b --workers 10 --n-dirs 20" > benchmarks/toy_orientation_sweep/findoptimal_run.log 2>&1 &
uv run python scripts/findoptimal_sweep.py summarize   # findoptimal_sweep_summary.{txt,json}, perturbation_sweep_vs_findoptimal.png
```
Outputs in `benchmarks/toy_orientation_sweep/`: `findoptimal_a_raw.npz`, `findoptimal_b_raw.npz` (1.6 MB), `findoptimal_sweep_summary.{txt,json}`, `perturbation_sweep_vs_findoptimal.png`, `findoptimal_pilot.log`, `findoptimal_run.log`. Per-task results are cached in `scripts/findoptimal_sweep_cache/` (gitignored).

**Caveat on identical noise (added after review).** The optimizers' realistic images carry only the pass-1 realism edits. The net's passes 2-3 re-render the windows, and their pixel noise attaches to the entry indices of the rebuilt ROI set, so "same nuisances" holds strictly for the one-shot comparison only; the iterated-net comparison on realistic data is therefore not on identical noise. The one-shot paired numbers in the summaries are strictly like-for-like. The committed raw npz files and logs contain absolute local paths (not secrets).


## FindOptimal robustness — wrong candidates and CSL traps (2026-10-06)

**Question.** Full multi-level reconstruction from scratch (`AdaptiveVoxelReconstructor.reconstruct_voxel`) is wrong (> 1 deg) in 21/50 clean and 17/50 realistic single-voxel runs on the 500-grain sample; the wrong answers are mostly CSL relatives of the truth (25 Sigma3, 5 Sigma7, 4 Sigma5, 1 Sigma11, ...). In every wrong case the local cost at the truth is lower than at the returned answer, so the search fails, not the cost. (1) Where is the truth lost, and what cheap fixes work? (2a) Would a low-cost classifier make FindOptimal more robust? (2b) What training data would that need?

**Plan.**

Data: per-voxel images exactly as `scripts/findoptimal_sweep.py` experiment A (no full-sample simulation). The existing 50 voxels plus 150 more random ManyGrains voxels (seeded; excluding the 30 NN-dataset voxels; recorded) = 200 voxels. Variants: clean and realistic (`all`). 3 reconstruct seeds per (voxel, variant) for E0. Pilot timing first; budget about 6 h wall on 10 CPU workers (1 torch thread each, nohup, resumable caches in a gitignored directory); cut realistic seeds or voxel count if over budget, and say so. Error = cubic-symmetry-reduced misorientation; wrong = > 1 deg; CSL label by the Brandon criterion (Sigma <= 29).

- **E0 diagnosis.** A non-invasive recorder hook in `reconstruct_voxel` (default behaviour bit-identical, checked against the golden test) records per level every candidate's orientation, discrete-stage score, post-quick-MC local cost and kept/pruned flag, plus FindOptimal's per-candidate results. Per run: min error over candidates per level, rank of the best truth-basin candidate (1 and 3 deg basins). Wrong runs are classified S1 (truth basin never among level-0 candidates), S2 (present but pruned at level k; record whether ranking by the local cost would have kept it), S3 (present at hand-off but FindOptimal returns a trap). Also wrong rate per seed, voxels wrong in all 3 seeds vs some, dependence on r_perp, boundary proximity, and the number of reflections the truth shares with its CSL relatives.
- **E1 cheap fixes**, each on the same runs (wrong rate with 95% CI, median error of right answers, extra cost evaluations, runtime): F1 CSL-relative check after the search (post hoc from E0 final answers); F1b CSL expansion at each level; F2 keep more (top 1/2, more discrete candidates); F3 coarse-cost settings (n_q_max start 8, pixel_radius 1); F4 best of N seeds by local cost. F1b-F3 on a subset (all wrong runs plus equal right runs) if needed.
- **E2 classifier feasibility (2a).** Candidate dataset (harvested E0 candidates labelled basin / near / CSL class / other, plus synthetic CSL relatives and basin orientations with residuals drawn from the E0 coarse-level error distribution). Features from one cost evaluation per candidate (hit fractions by Q shell, {hkl} family, detector; peak counts; pixel_radius 0 and 3; CSL-aware shared-vs-nonshared reflection hit fractions). Logistic regression and gradient-boosted trees, split by voxel; cross-variant tests; ROC-AUC vs the local-cost baseline; precision at top; end-to-end use as a pruning reranker and as the final chooser vs F1; learning curve at 25/50/100/150 voxels.
- **E3 training-data requirements (2b)** as a written section backed by the E2 learning curve.

Deliverables: `scripts/findoptimal_robustness/`, tests, `benchmarks/findoptimal_robustness/`, four plots, `docs/findoptimal_robustness_report.md`, status update of `docs/todo_findoptimal_wrong_candidates.md`.

### Completion summary (2026-10-06)

Report: `docs/findoptimal_robustness_report.md`; scripts: `scripts/findoptimal_robustness/` (README there); results and the four plots: `benchmarks/findoptimal_robustness/`; tests: `tests/test_findoptimal_robustness.py`, `tests/test_findoptimal_refactor.py` (recorder hook bit-identical; suite 537 passed, 34 skipped). `AdaptiveVoxelReconstructor` gained optional attributes `recorder`, `keep_fraction`, `keep_union_discrete`, `n_q_start_offset`, `global_pixel_radius`, `extra_candidates`, `rank_key`; defaults reproduce the old behaviour exactly. scikit-learn added to the dev extra. Compute: E0 about 57 min on 10 workers, everything about 4 h.

**Diagnosis (200 voxels x 3 seeds, clean and realistic; seed 0 reproduces experiment A exactly).** Wrong (> 1 deg): 204/600 = 34.0% clean (CI 30.3-37.9), 143/600 = 23.8% realistic (20.6-27.4); the cost at the truth is below the cost of the returned answer in 347/347 wrong runs. The truth basin (within 3 deg) is present at level 0 in 98-99% of runs; it is lost by the **pruning rank** (S2: present, then not among the best quarter by the post-quick-MC local cost) in 88% (clean) and 83% (realistic) of the wrong runs; S1 (never a level-0 candidate) 6% / 3.5%; S3 (reaches FindOptimal, trap returned) 4% / 13%. Mechanism: a basin candidate after the quick MC is still 1-2 deg off, the cost is sharp to 0.2-0.5 deg, so its cost is 0.94-0.99 (cost at the exact truth 0.06 / 0.27) and it loses a ranking among near-1 costs to a wrong candidate (Sigma3 62-66% of the wrong answers, then Sigma5, Sigma7). FindOptimal receives only 2-4 candidates (the set shrinks 4x per level), so `max_discrete_candidates` never binds. No C++/Python difference in candidate selection found (C++ `floor(n/4)` can be 0). Failure is strongly voxel-dependent (28 voxels wrong in all 3 seeds, 7.9 expected if independent), no consistent dependence on r_perp or grain-boundary proximity (the runs are not independent; clean not-near-boundary 41.1%, n=141, against near 31.8%, n=459; realistic first r_perp quartile 0.31 against 0.21-0.22), weakly related to the fraction of reflections shared with Sigma3/Sigma7 relatives (Spearman about +0.1).

**Fixes (wrong rate clean / realistic).** Baseline seed 0 (n=200): 33.5% / 26.0%. F1 CSL-relative check after the search (all Sigma<=29 relatives, quick MC, then FindOptimal on the best 3 per Sigma subset: `f1_run.py` refines the union of the top 3 of the Sigma<=11 subset and the top 3 of the Sigma<=29 subset, up to 6 orientations, and `summarize_fixes.f1_answer` uses the top 3 of the chosen subset): 3.5% / 3.3% (n=600 each), 0 of the 853 right runs broken, +11% evaluations, +10% time. F1b relatives added at every level: 0.0% / 1.5% (n=200), +26-28% evaluations (+10.9 s on 24.5 / 26.5 s, i.e. +41-44% time). F2a keep 1/2: about 11% / 5.5% implied (+15-23%); F2b union with the discrete-score rank: 11% / 9.6%. F3a coarse levels at full Q_max: 15% / 12% (+30-45%, breaks some right runs); F3b coarse pixel radius 1: 56% / 40% (worse). F4 best of 2 / 3 seeds: 20.5% / 14.0% clean, 11.0% / 6.5% realistic (2-3x cost). F2/F3 were run on a wrong+right subset (seed 0) with implied overall rates; F3c not run.

**Classifier (2a).** One forward pass per candidate (62 features, 6 ms), logistic regression and gradient-boosted trees, grain-disjoint folds, 360k candidates: ROC-AUC basin vs trap 0.997 (trees) against 0.963 for the local cost; pruning recall 0.985 against 0.900; the cost features and the CSL-aware features add almost nothing, the tolerant (radius 1, 3) hit fractions carry the gain. End to end, as the pruning rank: 4.0% / 4.5% wrong (seed 0), +4-6% evaluations; with F1 2.0% / 3.5%; not better than F1 (1.5% / 4.0%) or F1b (0.0% / 1.5%). As the final chooser it adds nothing (the cost already separates basin from trap; worse on realistic data). Learning curve: 25 voxels give most of the pruning gain, 100-150 the full AUC. **Training data (2b):** truth labels from simulation, search-generated hard negatives plus synthetic CSL relatives, 4.9% positives, grain-disjoint evaluation, domain-shift risks listed in the report.

**Recommendation.** Add F1 (Sigma<=11 gets most of it); consider F1b; keep more candidates per level; do not shrink the coarse pixel tolerance; the classifier is not needed for CSL failures. **Deviations:** F2/F3 on a subset; F3c not run; first classifier analysis used voxel-disjoint folds, redone grain-disjoint after finding 28 grains spanning folds (headline AUC and pruning recall within 0.004; larger differences elsewhere: ablation AUC B without tolerant features 0.853 against 0.868, final precision train-clean/test-realistic 0.831 against 0.847, end-to-end C-a+F1 realistic 3.5% against 2.0%, C-b final realistic 4.2% against 5.0%; both kept).

**Decision (2026-10-06, project owner).** The CSL-relative fixes F1 and F1b stay OFF by default so the Python port keeps matching the C++ reconstruction's behaviour; they remain opt-in and may be re-enabled later. Evidence (`benchmarks/findoptimal_robustness/fixes_summary.txt`): baseline wrong rate 34% clean / 24% realistic; F1 (Sigma<=11) 3.8% / 3.5% for about +5% evaluations, broke 0 right runs; F1b 0% / 1.5% (seed 0, n=200) for +26-28% evaluations (+41-44% time). How to enable: F1b is the attribute `AdaptiveVoxelReconstructor.extra_candidates` (hook `extra_candidates(level, candidates)`; set it to `make_expander(top_k=3, max_sigma=11)` from `scripts/findoptimal_robustness/fixes_run.py`, which is what `apply_knobs` does for `csl_expand=3`). F1 is not an attribute of the reconstructor: it is the post-search step in `scripts/findoptimal_robustness/f1_run.py` (`task()` generates the relatives with `csl.csl_relatives` and runs the quick MC on each; `refine_one` is only the refine step, `refine_from_candidates`; the answer is the lowest final local cost); applying it in production means calling that logic after `reconstruct_voxel`. What would justify turning them on: real data showing CSL traps (Sigma3/5/7 twins returned instead of the truth), or a C++ port of the same fix so both implementations stay equal.

**Known limitations.**
- The E2 classifier is trained with a 3 degree basin label (`fit` in `scripts/findoptimal_robustness/e2_models.py`) but is mostly evaluated at 1 degree.
- Trained classifier models live only in the gitignored cache (`scripts/findoptimal_robustness/cache/models/`); they must be retrained to be used.
- All results use per-voxel images with at most 3 distractor sources, not full-sample renders.
- F2/F3 are seed-0 subset estimates (wrong+right subset with implied overall rates).

## Hybrid NN finisher, Q_max-8 cost proxy, run-time profiling

Status: **complete** (started and finished 2026-10-06; completion summary at the end of this section). Branch `feature/nn-hybrid-proxy-profiling` with sub-task branches `feature/nn-hybrid-proxy-profiling-hybrid`, `-proxy` and `-profiling` (git cannot hold `feature/x/task` while `feature/x` exists, so the `/task` form of the gitflow rules was replaced by `-task`). Full plan: (plan file, deleted at completion). No changes to `icenine/reconstructor.py` are planned; evaluation counts and stage times come from wrappers in the scripts.

**Task 1: hybrid network -> seeded FindOptimal** (`scripts/nn_hybrid/`, `benchmarks/nn_hybrid/`).
Goal: does a network estimate, passed to `refine_from_candidates` (FindOptimal), beat the network alone and FindOptimal alone? Reuses the 50 voxels x 10 radii x 20 directions x 2 variants cases of the perturbation / FindOptimal sweeps, B's image and B's MC seed so results pair with `findoptimal_b_raw`. Pipelines: H0 FindOptimal alone; N1/N3 net x1 / x3 alone; H1 net x1 -> FO; H3 net x3 -> FO; H3c net x3 -> FO with a covariance-sized box b = clip(3 sigma_max, 0.329, 2.0 deg) (H3c: run, dropped, D3); H3m net x3 -> MC with box 3 sigma_max (H3m: run, dropped, D4); HG Huber GN x3 -> net x3 -> FO (r = 1.5, 2, 3, 5); post hoc rules (H3 falling back to HG when pass-1 sigma_max > tau; min cost of H0 and H3). Primary net `realistic_s0` on both variants; seed replicate `realistic_s1` (realistic, H3 only).
Success (realistic): H3 median <= 0.035 deg and fraction < 0.1 deg >= 0.9 for r <= 3; wrong <= 1% for r <= 3; H3 beats N3 in >= 70% paired for r <= 3; H3 beats H0 in >= 70% for r >= 0.75; at r <= 0.5 H3 not worse than H0 (median within 0.005 deg, win >= 45%); HG median <= 0.1 deg at r = 5.
Decision points: H1 within 0.005 deg of H3 at r <= 2 -> recommend x1; H3 stuck at the net's ~0.06 deg floor -> compare cost at truth vs at the result, try FO from the H0 result (min-cost rule); H3c not better than H3 where b > 0.329 -> drop the covariance box; H3m within 0.01 deg of H3 at under half the time -> H3m is the recommended finisher. Risk 8: if clean-data criteria fail, run H3 with `clean_s0` on clean data.
Budget: ~3.5-4 h on 10 workers.

**Task 2: Q_max-8 cost proxy as the pruning key** (`scripts/coarse_proxy/`, `benchmarks/coarse_proxy/`). Train low-Q (Q4/Q5) feature models on the E2 candidate set, labelled by "basin cost" (Q8 cost the candidate would reach if refined), as the `rank_key` for pruning; grain-disjoint folds; offline AUC / pruning recall, then end to end on 200 voxels x 2 variants vs baseline, F1, F1b, E2 rerank. Decision points D1 (Q5-only recall < 0.95 -> F-lowQ + cost8 only), D2 (recall at keep 1/8 >= 0.98 -> also run keep 1/8; row (ii) keep 1/8: not run, D2), D3 (regression vs classifier). Success: rerank wrong <= E2 rerank (4.0% / 4.5%) at less proxy time; or keep-1/8 wrong <= 5% at evaluations <= baseline; or proxy + F1 wrong <= F1 at fewer evaluations than F1b. Budget ~3 h.

**Task 3: does the NN reduce total run time?** (`scripts/nn_hybrid/stage_timer.py`, profiling scripts). Use cases U0 (no start), U1 (start within 0.1 deg), U2 (start 0.5-3 deg); stage timers and evaluation counters by wrapping; single-worker interleaved timing with `OMP/MKL_NUM_THREADS=1`; cProfile and torch.profiler. "Helps" = >= 20% lower median wall time at matched accuracy (U0), beats H0 on time and accuracy (U1), >= 1.5x faster than H0 + fallback or the only method with wrong <= 1% (U2). Budget ~3.5 h single worker + 0.5 h profiling.

Total budget ~10-11 h wall. Caveats carried in every table: per-voxel images (at most 3 distractor sources), not full-sample renders; noise is only strictly paired for pass 1 + finisher (H1 is the strictly paired hybrid).

### Task 1 results: hybrid network -> seeded FindOptimal (2026-10-06)

Code (reviewed and revised after the first write-up; numbers below are those of the regenerated summaries): `scripts/nn_hybrid/` (`run.py`, `summary.py`, `stage_timer.py`, `run_all.sh`), tests `tests/test_nn_hybrid.py`, results `benchmarks/nn_hybrid/` (`summary.txt` / `summary.json` are the full tables; `summary_s1_replicate.txt`, `summary_clean_s0_fallback.txt`, `pilot_summary.txt`, `wall_times.json`, raw `*_raw.npz`). Nothing in `icenine/` changed. Cases: 50 voxels x 10 radii x 20 directions x 2 variants = 20,000 per pipeline; metrics over the sweep's pass-1 success mask (9,000 cases per variant at r <= 3, and at r = 5 940 clean / 994 realistic: 9,940 clean + 9,994 realistic = 19,934 in all; the 60 clean and 6 realistic r = 5 cases whose pass 1 failed are fall-backs, FindOptimal from the start, counted in none of the tables; they have median 41.9 / 47.8 deg and 100% wrong, H0 and H3 identical). Error = cubic-reduced misorientation; wrong = > 1 deg; "realistic" is the sweep variant `all`. Net: `realistic_s0` on both variants.

**Checks.** The net passes of the copied pass loop reproduce `perturbation_sweep_raw` err_angle to 2.4e-07 deg (19,934 cases); H0 (FindOptimal alone, our images and seeds) reproduces `findoptimal_b_raw` R_final bit-for-bit in 20,000/20,000 cases, with equal evaluation counts. `tests/test_findoptimal_refactor.py` still passes. Gotcha found: `findoptimal_b_raw` (and the new raws) store the voxel axis sorted by voxel index, the perturbation / optimizer / A raws in sweep order; everything is re-indexed.


**Realistic data, median error (deg, cubic-reduced)**

| r (deg) | H0 | N3 | H1 | H3 | H3c | H3m | HG |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.0189 | 0.0677 | 0.0207 | 0.0205 | 0.0214 | 0.0307 | - |
| 0.1 | 0.0248 | 0.0623 | 0.0205 | 0.0197 | 0.0208 | 0.0313 | - |
| 0.25 | 0.037 | 0.0652 | 0.0203 | 0.02 | 0.0219 | 0.0301 | - |
| 0.5 | 0.0552 | 0.0667 | 0.0208 | 0.0212 | 0.0217 | 0.0304 | - |
| 0.75 | 0.0739 | 0.0658 | 0.0212 | 0.02 | 0.0209 | 0.0301 | - |
| 1 | 0.0905 | 0.0666 | 0.0207 | 0.0203 | 0.0207 | 0.0313 | - |
| 1.5 | 0.223 | 0.062 | 0.0236 | 0.0207 | 0.0213 | 0.0294 | 0.0212 |
| 2 | 0.869 | 0.073 | 0.034 | 0.0218 | 0.022 | 0.0332 | 0.0215 |
| 3 | 2.98 | 0.0803 | 0.0945 | 0.0214 | 0.0225 | 0.0407 | 0.0199 |
| 5 | 6.02 | 1.3 | 4.3 | 0.175 | 0.132 | 1.25 | 0.0199 |

**Realistic data, fraction < 0.1 deg**

| r (deg) | H0 | N3 | H1 | H3 | H3c | H3m | HG |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.990 | 0.640 | 0.973 | 0.978 | 0.957 | 0.816 | - |
| 0.1 | 0.961 | 0.664 | 0.974 | 0.980 | 0.960 | 0.836 | - |
| 0.25 | 0.862 | 0.644 | 0.974 | 0.976 | 0.961 | 0.842 | - |
| 0.5 | 0.720 | 0.654 | 0.979 | 0.978 | 0.963 | 0.839 | - |
| 0.75 | 0.570 | 0.647 | 0.970 | 0.975 | 0.966 | 0.832 | - |
| 1 | 0.533 | 0.650 | 0.979 | 0.973 | 0.960 | 0.830 | - |
| 1.5 | 0.339 | 0.664 | 0.953 | 0.971 | 0.950 | 0.836 | 0.978 |
| 2 | 0.172 | 0.613 | 0.844 | 0.962 | 0.942 | 0.801 | 0.971 |
| 3 | 0.012 | 0.568 | 0.508 | 0.916 | 0.908 | 0.740 | 0.973 |
| 5 | 0.000 | 0.238 | 0.016 | 0.472 | 0.480 | 0.308 | 0.958 |

**Realistic data, wrong rate (> 1 deg, %)**

| r (deg) | H0 | N3 | H1 | H3 | H3c | H3m | HG |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 0.1 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 0.25 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 0.5 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 0.75 | 0.10 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 1 | 0.90 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 1.5 | 15.90 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 2 | 46.80 | 0.70 | 1.10 | 0.10 | 0.30 | 0.70 | 0.00 |
| 3 | 94.60 | 6.70 | 19.40 | 4.10 | 3.90 | 6.60 | 0.00 |
| 5 | 100.00 | 52.92 | 95.17 | 43.56 | 40.34 | 52.21 | 2.31 |

**Clean data, median error (deg, cubic-reduced)**

| r (deg) | H0 | N3 | H1 | H3 | H3c | H3m | HG |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.0197 | 0.0111 | 0.0105 | 0.0105 | 0.0105 | 0.00847 | - |
| 0.1 | 0.0257 | 0.0111 | 0.0114 | 0.0104 | 0.0104 | 0.00895 | - |
| 0.25 | 0.0426 | 0.0111 | 0.0123 | 0.0105 | 0.0105 | 0.00871 | - |
| 0.5 | 0.0557 | 0.0111 | 0.0131 | 0.0105 | 0.0105 | 0.009 | - |
| 0.75 | 0.0748 | 0.011 | 0.0151 | 0.0105 | 0.0105 | 0.00865 | - |
| 1 | 0.0892 | 0.0111 | 0.0165 | 0.0105 | 0.0105 | 0.0089 | - |
| 1.5 | 0.234 | 0.011 | 0.0178 | 0.0105 | 0.0105 | 0.00887 | 0.0104 |
| 2 | 1.18 | 0.011 | 0.0181 | 0.0105 | 0.0105 | 0.00835 | 0.0105 |
| 3 | 3.07 | 0.0111 | 0.0206 | 0.0107 | 0.0107 | 0.00833 | 0.0105 |
| 5 | 23.9 | 0.0111 | 0.0365 | 0.0107 | 0.0107 | 0.00915 | 0.0105 |

**Clean data, fraction < 0.1 deg**

| r (deg) | H0 | N3 | H1 | H3 | H3c | H3m | HG |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.985 | 1.000 | 1.000 | 1.000 | 1.000 | 0.998 | - |
| 0.1 | 0.958 | 1.000 | 1.000 | 0.999 | 0.999 | 1.000 | - |
| 0.25 | 0.861 | 1.000 | 1.000 | 0.999 | 0.999 | 0.998 | - |
| 0.5 | 0.707 | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 | - |
| 0.75 | 0.584 | 1.000 | 0.998 | 1.000 | 1.000 | 1.000 | - |
| 1 | 0.522 | 1.000 | 0.997 | 1.000 | 1.000 | 0.998 | - |
| 1.5 | 0.339 | 1.000 | 0.997 | 0.999 | 0.999 | 0.997 | 0.999 |
| 2 | 0.164 | 1.000 | 0.990 | 0.999 | 0.999 | 0.998 | 1.000 |
| 3 | 0.008 | 1.000 | 0.965 | 1.000 | 1.000 | 1.000 | 1.000 |
| 5 | 0.000 | 1.000 | 0.871 | 0.999 | 0.999 | 0.999 | 0.999 |

**Clean data, wrong rate (> 1 deg, %)**

| r (deg) | H0 | N3 | H1 | H3 | H3c | H3m | HG |
|---|---|---|---|---|---|---|---|
| 0.05 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 0.1 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 0.25 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 0.5 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 0.75 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 1 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | - |
| 1.5 | 18.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 2 | 52.70 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 3 | 96.10 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 5 | 100.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |


Paired win rates (realistic, ties half; full tables incl. MC and Huber in `summary.txt`): H3 vs N3 0.79-0.81 for every r <= 3 (0.68 at r = 5); H3 vs H0 0.49, 0.55, 0.68, 0.75, 0.81, 0.83, 0.87, 0.93, 0.98, 0.88 for r = 0.05 ... 5; H3 vs H1 about 0.50 for r <= 1.5, 0.64 (r = 2), 0.79 (r = 3), 0.85 (r = 5); HG vs H3 0.50 / 0.51 / 0.53 / 0.78 at r = 1.5 / 2 / 3 / 5. Mean finisher cost evaluations per case (realistic): H0 3060 at r = 0.05, 2260 at 0.25, 1300 at 0.5, 696 at 1.5, 373 at 5; H1 about 2900 (r <= 1), 2360 (r = 2), 1290 (r = 3); H3 2890-3070 (r <= 3), 1730 (r = 5); H3c about 2420; H3m about 2240; HG 2940-3620.

**Time (10 workers, contended, mean seconds per realistic case; all radii).** Net stage: prepare 0.13 + render 0.06 + forward pass 0.015, about 0.2 s per case (about 0.15 s without the harness rendering). Finisher: H0 0.93 s (1358 evaluations), H1 1.67 s (2438), H3 1.93 s (2829), H3c 1.60 s (2344), H3m 1.58 s (2313), HG 2.18 s (3222) plus 0.25 s Huber GN. About 94% of the finisher time is cost evaluations (about 0.65 ms each: H3 1.82 s of evaluations in a 1.93 s finisher). Time inclusive of its evaluations: FindOptimal 0.11 s and VarianceMinimizing 1.82 s for H3. VarianceMinimizing time depends on how close the start is to the truth (H0: 1.98 s at r = 0.05 falling to 0.15 s at r = 5; H3: about 1.9 s at every r <= 3); `converged` is False in every case of every pipeline, H0 included, and `refine_from_candidates` always runs VarianceMinimizing, so no convergence explanation is offered. On the radii common to all pipelines (r = 1.5, 2, 3, 5; `summary.txt`) the finisher takes H0 0.35 s (511 evaluations), H1 1.16 s (1700), H3 1.84 s (2682), H3c 1.53 s (2227), H3m 1.64 s (2398), HG 2.18 s (3222), so HG costs +0.34 s over H3 there (plus 0.25 s of Huber GN, a figure that includes the harness's own prepare / render of GN passes 2-3). The stored H3m times include harness work (pixel grouping and the start / truth quality evaluations of `run_mc_adam`, and its stage appears under the old name `find_optimal`); `run.py` now times only the MC (`o["seconds"]`, stage `mc_optimize`): on one voxel the same bit-identical H3m result took 0.65 s against 0.79 s stored (about 18% less), which biases D4's time ratio against H3m, but D4 is decided on accuracy. Pipeline wall times (10 workers, `wall_times.json`): H3 66 min, H0 36, H1 61, H3c 11 (copies H3 where the box is not larger), H3m 43, HG 32 (4 radii), s1 replicate 37 (realistic only), clean_s0 fall-back about 30; total about 5 h, longer than the planned 3.5-4 h because the hybrid finishers start near the truth and run about twice as long as H0 on average.

**Success criteria (realistic, `realistic_s0`).**
- C1 H3 median <= 0.035 deg and fraction < 0.1 deg >= 0.9 for r <= 3: **MET** (medians 0.0197-0.0218 deg, fraction 0.916-0.980; worst r = 3 with 0.916).
- C2 H3 wrong <= 1% for r <= 3: **NOT MET** only at r = 3 (4.10%, interval 3.0-5.5); r = 2 is 0.10%, all r <= 1.5 are 0.00%. (N3 has 6.7% at r = 3, H1 19.4%.)
- C3 H3 beats N3 in >= 70% paired for r <= 3: **MET** (0.79-0.81).
- C4 H3 beats H0 in >= 70% for r >= 0.75: **MET** (0.81, 0.83, 0.87, 0.93, 0.98, 0.88).
- C5 at r <= 0.5 H3 not worse than H0 (median within 0.005 deg, win >= 45%): **MET** (median H3 - H0 = +0.002, -0.005, -0.017, -0.034 deg at r = 0.05, 0.1, 0.25, 0.5; win 0.49, 0.55, 0.68, 0.75). At r = 0.05 H0 is marginally better (0.0189 against 0.0205 deg).
- C6 HG median <= 0.1 deg at r = 5: **MET** (0.0199 deg; 2.31% wrong, interval 1.5-3.4; H3 at r = 5: median 0.175 deg, 43.6% wrong).
Seed replicate (`realistic_s1`, H3 only): C3-C5 met; C1 not met (median 0.0237 deg but fraction < 0.1 deg 0.887 at r = 3 against 0.9); C2 not met at r = 3 (7.4%) and r = 2 (1.1%), 0.2% at r = 1.5; medians 0.0207-0.0237 deg for r <= 3; wins vs N3 0.79-0.81, vs H0 0.79-0.94 for r >= 0.75. So the r = 3 corner is borderline for C1 and fails C2 in both seeds, and r = 2 is 0.1% / 1.1% wrong.
Clean data (H3): C1, C2, C4, C5, C6 met (median 0.0104-0.0107 deg, wrong 0.00% at every radius, H3 beats H0 0.91-1.00 for r >= 0.75). C3 reads "not met" (0.52) only because H3 and N3 agree within the 0.002 deg tie band in 93% of cases: the net alone already reaches 0.011 deg on clean data and FindOptimal adds nothing. Because of that formal miss the risk-8 fall-back was run: H3 with `clean_s0` on clean data gives the same picture (median 0.0102-0.0103 deg, wrong 0.00% everywhere, wins vs H0 0.91-1.00, vs N3 0.51-0.52 = ties; `summary_clean_s0_fallback.txt`), so the clean-data conclusion does not depend on the network.

**Decision points.**
- D1 (H1 within 0.005 deg of H3 at r <= 2 -> recommend x1): median H1 - H3 = +0.0003, +0.0008, +0.0002, -0.0004, +0.0012, +0.0004, +0.0029, +0.0122 deg for r = 0.05 ... 2, so within 0.005 up to r = 1.5 but not at r = 2; H1 is also worse in the tail (r = 2: 1.10% wrong against 0.10%; r = 3: 19.4% against 4.1%; fraction < 0.1 deg 0.844 against 0.962 at r = 2). **Decision: keep x3** (x1 is equivalent for r <= 1.5 and saves two net passes, about 0.2 s of preparation and rendering, but the finisher time is unchanged).
- D2 (H3 stuck at a floor): the net's floor is 0.065 deg; H3 reaches 0.020 deg, but 25-30% of its cases stay above 0.035 deg. The finisher's result has a higher cost than the truth in about 96% of ALL cases (H3: 95.3-96.6% per radius for r <= 3, 95.5% at r = 0.05, 95.8% at r = 3; H0: 93.8% at r = 0.05 rising to 100% at r = 3), not only in the tail (97.5-100%), so the finisher generally stops before the cost minimum; the cause (step size, variance stopping rule, or MC stochasticity) is not diagnosed (the earlier "flat plateau" reading was unsupported). The min-cost rule (b) of H0 and H3 (needs both finishers, about 1.5x the H3 finisher cost: H0 adds on average 0.93 s and 1358 evaluations, 3060 at r = 0.05) lowers the median to 0.0137-0.0214 deg (0.0137 against 0.0205 at r = 0.05, 0.0181 against 0.0200 at r = 0.75, 0.0214 unchanged at r = 3), H0 being chosen in 59% of cases at r = 0.05, 12% at r = 1.5, 2% at r = 3, and leaves the wrong rate unchanged (r = 3: 4.00% against 4.10%). The plan's extra FindOptimal from the H0 result was not run (the min-cost rule is the cheaper version of it, since the H0 run is already available in the study).
- D3 (H3c vs H3 where b > 0.329 deg): b exceeds the default box in 160-192 of 1000 cases per radius at r <= 3 (540 of 994 at r = 5). Pooled over those 2093 cases H3c is worse: median 0.0494 against 0.0339 deg, H3c wins 43%. Per-radius medians 0.034-0.041 against 0.022-0.030 deg; wrong rates equal or slightly lower only at r = 3 (10.9% against 12.0% in that subset) and r = 5 (61.5% against 67.4%), higher at r = 2 (1.64% against 0.55%). **Decision: drop the covariance box.**
- D4 (H3m within 0.01 deg of H3 at under half the time): realistic: median H3m - H3 = +0.009 to +0.019 deg and time 0.75-0.80 of H3 (not under half), fraction < 0.1 deg 0.82-0.84 against 0.97-0.98 for r <= 1.5 (0.801 against 0.962 at r = 2, 0.740 against 0.916 at r = 3): **no, H3m is not recommended** (the stored H3m times include harness work, see Time; removing it would give about 0.62x, still not under half). Clean data: H3m is 0.0014-0.0023 deg better at 0.33-0.39x the time, because on clean data FindOptimal and the variance loop add nothing; that is a property of the clean data, not a recommendation for the realistic case.
- Post hoc (a), falling back to HG where the pass-1 predicted sigma_max exceeds tau: sigma_max does not flag the failures. tau = 0.1 deg selects HG in 22% / 46% / 82% / 99% of cases at r = 1.5 / 2 / 3 / 5 and gives 0.00% / 0.00% / 1.00% / 2.92% wrong; tau >= 0.2 selects HG in 1-5% of cases at r <= 3 and leaves r = 3 at 4.1% and (tau = 0.2) r = 5 at 37.6% wrong. HG alone is 0.00% at r <= 3 and 2.31% at r = 5, so a sigma_max rule is no better than using HG whenever a large start (r >= 2) is possible.

**Findings.** (1) Net x3 -> FindOptimal (H3) has a median of about 0.020 deg at every r <= 3, 3x lower than the net alone (0.062-0.080) and lower than FindOptimal alone for r >= 0.25 (0.037-2.98 deg); H0 is marginally better at r = 0.05 (0.0189 against 0.0205 deg) and H3 slightly better at r = 0.1 (0.0197 against 0.0248, win 0.55). (2) Its remaining tail is a net failure (r = 3: 4.1% wrong, r = 5: 43.6%), which FindOptimal's 0.33 deg box cannot repair. (3) Huber Gauss-Newton x3, then net x3, then FindOptimal (HG) removes the tail: median 0.020-0.021 deg and 0.00% wrong for r = 1.5-3, 2.31% wrong at r = 5, for +0.25 s of GN and, on the common radii 1.5-5, +0.34 s finisher time over H3 (2.18 s against 1.84 s); it ties H3 at r = 1.5-2 (win 0.50 / 0.51), is better at r = 3 (0.53) and r = 5 (0.78). (4) With a start known to be within 0.1 deg, H0 and H3 cost the same number of evaluations (about 3000), H0 is marginally more accurate at r = 0.05 and H3 slightly more accurate at r = 0.1; H3's evaluation count stays near 2900 while H0's falls with r (1300 at 0.5, 696 at 1.5), so for starts this close FindOptimal alone is the cheaper choice for the same accuracy and the hybrid pays only where H0's accuracy degrades (win rates 0.68-0.75 at r = 0.25-0.5, Task 3 quantifies the time on one worker).

**Caveats.** Per-voxel images (at most 3 distractor sources), not full-sample renders. Noise pairing: the finisher always sees case B's image; the net's pass 1 sees the same noise, but passes 2-3 (and the HG net stage) re-render the windows, so only H1 is strictly paired; H3 / HG results are paired in voxel, direction and truth but not in the noise of the net's later passes. Times are 10-worker contended times; single-worker times belong to Task 3. The second net seed (`realistic_s1`) was run for H3 only. H3m uses the optimizer sweep's MC (hard pixel cost, 3500 steps x 2 restarts), not the reconstructor's cost function. The sub-task branch is `feature/nn-hybrid-proxy-profiling-hybrid`, because git cannot hold `feature/nn-hybrid-proxy-profiling` and `feature/nn-hybrid-proxy-profiling/hybrid` together.

### Task 2 results: Q_max-8 cost proxy as the pruning key (2026-10-06)

Code: `scripts/coarse_proxy/` (`lowq_dataset.py`, `labels.py`, `models.py`, `endtoend.py`, `timing.py`, `summary.py`, shared `base.py`), the optional `q_max` of `scripts/findoptimal_robustness/features.py` (default None = the E2 features, bit-identical; tested against the stored E2 table), an optional `--cache-root` of `f1_run.py`, tests `tests/test_coarse_proxy.py`. Results: `benchmarks/coarse_proxy/` (`summary.txt` / `summary.json` hold every number below; `offline_metrics.*`, `label_validation.json`, `qlevels.json`, `timing.json`). Caches (gitignored): `scripts/coarse_proxy/cache/`. Nothing in `icenine/` changed. Branch `feature/nn-hybrid-proxy-profiling-proxy` (the plan's `.../proxy` name is blocked by the existing parent ref).

**Q levels (verified, `qlevels.json`).** `FeatureExtractor.q_levels` = 3.015, 3.481, 4.923, 5.773, 6.029, 6.962, 7.587, 7.784 with 8, 6, 12, 24, 8, 6, 24, 24 reflections (112 in all). |q| <= 3: no reflection (Q3 dropped); <= 4: {111}, {200} (14 reflections); <= 5: adds {220} (26). Exactly the plan's expectation. The cut-off is `max_q` in Angstrom^-1, as in `VoxelCostFunction`.

**Label choice (`labels.py`).** Regression target `y_bcost` = the Q8 local cost a candidate would reach if refined: within 3 deg of the truth -> cost at the truth after `refine_from_candidates([truth])`; within 3 deg of an exact Sigma <= 29 relative -> cost after a short refinement of that relative (quick MC of the last level, then one MC of up to 200 steps, box 0.329 deg; 270 relatives per case, about 20 s per case on one worker); otherwise censored = the candidate's own cost (harvested candidates; synthetic ones without a label). Never the candidate's own cost for a truth-basin or relative-basin candidate. Counts over the 360k E2 candidates: truth basin 17,688; relative basin 199,721; censored 134,367; unlabelled 8,220. Classification labels: y3 (< 3 deg), y1 (< 1 deg), level-matched y_lm (y3 at levels 0-2 and synthetic, y1 at level 3 and FindOptimal results). The label values themselves are nearly bimodal: the truth-basin cost is 0.068 (clean; 10-90% over the 200 cases 0.043-0.095) / 0.263 (realistic; 0.207-0.308), the short-refined relatives have a median cost of 0.963 / 0.969 (median over cases of the per-case median) and the censored candidates are near 1, so the regression is mostly a basin detector. The censored branch was validated on 2,000 random censored harvested candidates (5 per case) refined with the same short procedure: **Spearman 0.861 pooled**, which passes the 0.8 threshold; but ranking happens within a case and level, and there the picture is borderline: mean within case 0.790 (< 0.8), by level 0.810 / 0.872 / 0.783 for L0 / L1 / L2 (L2 below 0.8). The regression label was kept on the pooled criterion; D3 shows the classifiers tie it offline (0.977 vs 0.977 pooled recall), so the fallback to the classifiers would not have changed the offline recall. The refinement lowers the cost only slightly (median 0.986 -> 0.979, lower for 75.8%), i.e. the censored labels are close to their refined values in rank but are almost all near 1.

**Offline protocol.** 4 grain-disjoint folds (`e2_models.fold_of`; tested), HistGradientBoosting (classifiers: the E2 model; regressors: same hyper-parameters), trained on clean + realistic together and on every candidate incl. the synthetic ones (built from the truth, so training-only information), tested on held-out grains. Sets A (contested harvested), C (basin vs harvested near-CSL) and the pruning recall D / final precision E are harvested candidates; set B is synthetic and reported as such. D = per (voxel, variant, seed, level) group, is the best basin candidate (< 3 deg; < 1 deg at level 3) among the max(1, int(n x frac)) best by the score. E2 full 62-feature GBT reproduces the earlier pooled L0-2 recall (0.985). "hand" = mean of overall hit0, hit_any, hit1, hit3 (untrained). "+c8" adds the free post-quick-MC Q8 local cost, "cache" = low-Q columns of the cached full E2 table (needs the full pass, so not deployable), "lowq" = a new pass over the |q| <= Q reflections plus costs at max_q = Q. rho = Spearman of the score with -y_bcost (harvested labelled candidates); for classifiers it is only a sanity number.

**Offline, test on clean data (held-out grains).**

| model | AUC A | AUC B | AUC C | recall 1/4: L0 | L1 | L2 | L3@1deg | L0-2 | recall 1/8: L0 | L1 | L2 | L3@1deg | L0-2 | E | rho |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cost8 | 0.961 | 0.851 | 0.966 | 0.885 | 0.840 | 0.905 | 0.966 | 0.876 | 0.735 | 0.781 | 0.814 | 0.966 | 0.773 | 0.986 | 0.558 |
| cost_q5 | 0.913 | 0.840 | 0.931 | 0.727 | 0.766 | 0.871 | 0.981 | 0.781 | 0.574 | 0.687 | 0.820 | 0.981 | 0.682 | 0.982 | 0.390 |
| hand | 0.990 | 0.900 | 0.990 | 0.941 | 0.929 | 0.959 | 1.000 | 0.942 | 0.875 | 0.906 | 0.902 | 1.000 | 0.893 | 0.995 | 0.587 |
| cache5/clf3 | 0.990 | 0.917 | 0.989 | 0.959 | 0.937 | 0.959 | 0.996 | 0.952 | 0.906 | 0.920 | 0.946 | 0.996 | 0.922 | 1.000 | 0.143 |
| cache5/reg | 0.997 | 0.921 | 0.996 | 0.931 | 0.937 | 0.966 | 0.992 | 0.943 | 0.903 | 0.917 | 0.949 | 0.992 | 0.921 | 1.000 | 0.409 |
| lowq4/reg | 0.988 | 0.919 | 0.990 | 0.926 | 0.915 | 0.959 | 1.000 | 0.932 | 0.860 | 0.906 | 0.956 | 1.000 | 0.903 | 1.000 | 0.357 |
| lowq5/reg | 0.999 | 0.918 | 0.999 | 0.941 | 0.943 | 0.959 | 0.996 | 0.947 | 0.903 | 0.917 | 0.956 | 0.996 | 0.923 | 1.000 | 0.459 |
| lowq4+c8/reg | 0.998 | 0.912 | 0.998 | 0.964 | 0.952 | 0.959 | 1.000 | 0.959 | 0.944 | 0.926 | 0.953 | 1.000 | 0.940 | 1.000 | 0.562 |
| lowq5+c8/reg | 0.999 | 0.910 | 0.999 | 0.977 | 0.957 | 0.963 | 1.000 | 0.966 | 0.957 | 0.929 | 0.956 | 1.000 | 0.947 | 1.000 | 0.610 |
| lowq5+c8/clf3 | 0.994 | 0.916 | 0.994 | 0.980 | 0.954 | 0.959 | 1.000 | 0.965 | 0.954 | 0.937 | 0.956 | 1.000 | 0.949 | 1.000 | 0.304 |
| lowq5+c8/clf_lm | 0.995 | 0.916 | 0.995 | 0.982 | 0.954 | 0.963 | 1.000 | 0.967 | 0.949 | 0.943 | 0.956 | 1.000 | 0.949 | 1.000 | 0.304 |
| e2full/clf3 | 1.000 | 0.938 | 1.000 | 0.995 | 0.963 | 0.969 | 1.000 | 0.977 | 0.990 | 0.957 | 0.966 | 1.000 | 0.972 | 1.000 | 0.277 |
| e2full/reg | 1.000 | 0.931 | 1.000 | 0.992 | 0.977 | 0.969 | 0.996 | 0.981 | 0.987 | 0.963 | 0.966 | 0.996 | 0.973 | 1.000 | 0.642 |

Groups: n(L0..L3) = [392, 351, 295, 264], n(L0-2) = 1038, E groups 219, set sizes A 34747, C 42742

**Offline, test on realistic data (held-out grains).**

| model | AUC A | AUC B | AUC C | recall 1/4: L0 | L1 | L2 | L3@1deg | L0-2 | recall 1/8: L0 | L1 | L2 | L3@1deg | L0-2 | E | rho |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cost8 | 0.963 | 0.850 | 0.970 | 0.907 | 0.905 | 0.964 | 0.934 | 0.923 | 0.756 | 0.847 | 0.879 | 0.921 | 0.824 | 0.949 | 0.706 |
| cost_q5 | 0.923 | 0.846 | 0.940 | 0.759 | 0.815 | 0.918 | 0.950 | 0.826 | 0.590 | 0.738 | 0.885 | 0.943 | 0.729 | 0.949 | 0.551 |
| hand | 0.984 | 0.894 | 0.988 | 0.947 | 0.962 | 0.982 | 0.934 | 0.963 | 0.889 | 0.935 | 0.952 | 0.921 | 0.923 | 0.949 | 0.702 |
| cache5/clf3 | 0.987 | 0.916 | 0.991 | 0.952 | 0.981 | 0.988 | 0.946 | 0.973 | 0.912 | 0.967 | 0.970 | 0.937 | 0.948 | 0.959 | 0.115 |
| cache5/reg | 0.988 | 0.917 | 0.993 | 0.960 | 0.978 | 0.991 | 0.896 | 0.975 | 0.932 | 0.973 | 0.967 | 0.871 | 0.956 | 0.878 | 0.426 |
| lowq4/reg | 0.992 | 0.912 | 0.994 | 0.937 | 0.962 | 0.991 | 0.972 | 0.962 | 0.889 | 0.959 | 0.985 | 0.959 | 0.942 | 0.980 | 0.408 |
| lowq5/reg | 0.995 | 0.915 | 0.996 | 0.960 | 0.973 | 0.997 | 0.978 | 0.975 | 0.942 | 0.967 | 0.988 | 0.975 | 0.964 | 0.983 | 0.531 |
| lowq4+c8/reg | 0.995 | 0.906 | 0.996 | 0.977 | 0.970 | 0.988 | 0.972 | 0.978 | 0.960 | 0.965 | 0.979 | 0.968 | 0.967 | 0.966 | 0.659 |
| lowq5+c8/reg | 0.995 | 0.909 | 0.996 | 0.990 | 0.978 | 0.994 | 0.972 | 0.987 | 0.972 | 0.967 | 0.985 | 0.965 | 0.974 | 0.976 | 0.674 |
| lowq5+c8/clf3 | 0.994 | 0.914 | 0.995 | 0.982 | 0.975 | 0.997 | 0.981 | 0.984 | 0.967 | 0.975 | 0.985 | 0.975 | 0.975 | 0.990 | 0.359 |
| lowq5+c8/clf_lm | 0.993 | 0.913 | 0.994 | 0.987 | 0.978 | 0.991 | 0.981 | 0.985 | 0.967 | 0.975 | 0.982 | 0.975 | 0.974 | 0.973 | 0.360 |
| e2full/clf3 | 0.994 | 0.936 | 0.996 | 0.997 | 0.986 | 0.994 | 0.972 | 0.993 | 0.992 | 0.984 | 0.985 | 0.972 | 0.987 | 0.976 | 0.312 |
| e2full/reg | 0.994 | 0.931 | 0.996 | 0.997 | 0.989 | 0.994 | 0.972 | 0.994 | 0.990 | 0.984 | 0.982 | 0.965 | 0.985 | 0.969 | 0.688 |

Groups: n(L0..L3) = [398, 367, 331, 317], n(L0-2) = 1096, E groups 295, set sizes A 36641, C 54081

(pooled over both test variants: lowq5+c8/reg recall L0-2 0.977 at 1/4 and 0.961 at 1/8, E2 full/clf3 0.985 and 0.980, cost8 0.900 and 0.799; `offline_metrics.txt`.)

**Decision points (pooled held-out).**
- D1: F-cache Q5 GBT pruning recall (keep 1/4, L0-2) is 0.963 (classifier; 0.960 regressor) >= 0.95, so the low-Q signal is **not** mostly lost (E2 full 0.985). The best deployable proxy, F-lowQ Q5 + cost8 (regression), is 0.977, a gap of 0.008 to E2, which is not more than 0.01, so the "adds nothing beyond E2" branch was **not** taken; F-lowQ + cost8 was carried forward as the deployable proxy.
- D2: recall at keep 1/8 for the best deployable model is 0.963 (lowq5+c8/clf3; 0.961 for the regressor) < 0.98, so D2 does **not** hold and the keep_fraction = 1/8 row (ii) was **not run**. (Untested: what keep 1/8 would do end to end.)
- D3: regression 0.977 vs best classifier 0.977 for lowq5+c8 (0.969 vs 0.969 for lowq4+c8), within 0.005, so regression is carried (the end-to-end row uses `lowq5+c8` / `reg`).

**End to end** (200 voxels x 2 variants, seed 0, same images and rng as E0, model of the fold without the voxel's grain; per-voxel images; wall time is the 10-worker, contended `reconstruct_voxel` time and differs from the single-worker accounting). "evals" = global + local cost evaluations of the reconstructor; "incl. proxy" adds the scored candidates times the measured single-worker cost of one scored candidate in Q8-local-evaluation equivalents (proxy 2.95, E2 4.22). Row (i) = proxy rerank, (iii) = proxy rerank + F1 (`f1_run.py --source p_i --cache-root ...`).

| row | clean wrong | fixed / broken | med err right | evals global / local | incl. proxy | wall |
|---|---|---|---|---|---|---|
| E0 baseline | 67/200 = 0.335 [0.273, 0.403] | - | 0.0301 | 44478 / 4711 | 49190 | 24.5 s |
| F1 | 3/200 = 0.015 [0.005, 0.043] | 64 / 0 | 0.0272 | 44478 / 10044 | 54523 | 27.2 s |
| F1b | 0/200 = 0.000 [0.000, 0.019] | 67 / 0 | 0.0200 | 45492 / 17519 | 63011 | 35.4 s |
| E2 rerank | 8/200 = 0.040 [0.020, 0.077] | 59 / 0 | 0.0285 | 44505 / 5986 | 51528 | 26.6 s |
| E2 rerank + F1 | 4/200 = 0.020 [0.008, 0.050] | 63 / 0 | 0.0285 | 44505 / 10798 | 56339 | 29.1 s |
| (i) proxy rerank | 13/200 = 0.065 [0.038, 0.108] | 55 / 1 | 0.0281 | 44494 / 5739 | 50941 | 25.5 s |
| (iii) proxy + F1 | 0/200 = 0.000 [0.000, 0.019] | 67 / 0 | 0.0282 | 44494 / 10646 | 55848 | 27.9 s |

| row | realistic wrong | fixed / broken | med err right | evals global / local | incl. proxy | wall |
|---|---|---|---|---|---|---|
| E0 baseline | 52/200 = 0.260 [0.204, 0.325] | - | 0.0278 | 44778 / 7628 | 52406 | 26.5 s |
| F1 | 8/200 = 0.040 [0.020, 0.077] | 44 / 0 | 0.0257 | 44778 / 13007 | 57785 | 29.2 s |
| F1b | 3/200 = 0.015 [0.005, 0.043] | 50 / 1 | 0.0176 | 45789 / 20143 | 65931 | 37.3 s |
| E2 rerank | 9/200 = 0.045 [0.024, 0.083] | 46 / 3 | 0.0228 | 44824 / 9625 | 56106 | 29.8 s |
| E2 rerank + F1 | 7/200 = 0.035 [0.017, 0.070] | 47 / 2 | 0.0228 | 44824 / 14743 | 61223 | 32.6 s |
| (i) proxy rerank | 15/200 = 0.075 [0.046, 0.120] | 40 / 3 | 0.0285 | 44805 / 9111 | 55046 | 28.2 s |
| (iii) proxy + F1 | 4/200 = 0.020 [0.008, 0.050] | 50 / 2 | 0.0296 | 44805 / 14278 | 60214 | 30.9 s |

Seeds 1-2 for row (i) (paired with E0 seeds 0-2, 600 runs per variant): clean 42/600 = 0.070 [0.052, 0.093] (fixed 167, broken 5), realistic 33/600 = 0.055 [0.039, 0.076] (fixed 117, broken 7); evaluations incl. proxy 50975 / 55156 against baseline 49196 / 52409. Per-run change of row (i) vs the baseline (seed 0): global +16 / +27, local +1027 / +1483, scoring +708 / +1130 equivalents, total +3.6% / +5.0%, wall +0.9 s / +1.7 s (E2 rerank: +4.8% / +7.1%, +2.1 s / +3.3 s). The scored candidates per run are 240 / 384. Where row (i)'s extra local evaluations go (measured, `eval_split.py`: 20 voxels x 2 variants, seed 0, stage timer patches, the patched runs reproduce the unpatched `R_final` bit for bit; mean per run, baseline -> proxy):

| stage (local cost evaluations unless noted) | clean | realistic |
|---|---|---|
| discrete search, levels 1-3 (local) | 180 -> 212 (+32) | 324 -> 372 (+48) |
| discrete search, pixel-radius-3 evaluations (all levels) | 44476 -> 44490 (+14) | 44818 -> 44842 (+24) |
| quick MC, levels 1-3 | 647 -> 744 (+97) | 1127 -> 1308 (+181) |
| FindOptimal | 475 -> 956 (+481) | 940 -> 2020 (+1080) |
| VarianceMinimizing | 1239 -> 1750 (+511) | 2192 -> 1995 (-197) |
| total local, reconstructor | 4607 -> 5728 (+1121) | 7983 -> 9094 (+1111) |
| candidates handed to FindOptimal (evaluated) | 2.6 -> 4.8 | 5.0 -> 10.4 |
| candidates at levels 1 / 2 / 3 | 44.9 / 11.3 / 2.6 -> 47.6 / 15.2 / 4.8 | 77.1 / 20.4 / 5.0 -> 81.5 / 27.1 / 10.4 |
| proxy's own low-Q cost evaluations (not in the reconstructor count) | 239 local + 239 radius-3 | 401 + 401 |

So in these 40 runs most of the extra local evaluations are in FindOptimal (+481 clean, +1080 realistic; more candidates reach it) and, on clean data, VarianceMinimizing (+511; on realistic data it uses 197 fewer); quick MC (+97 / +181) and the discrete stage of the later levels (+32 / +48) add the rest, and level 0 is unchanged. (Stage split as measured; no cause is claimed beyond it.)

**Cost accounting (one worker, threads = 1, median of 1000 candidates, `timing.json`).** One Q8 local evaluation 0.520 ms (the unit); one Q5 radius-3 global evaluation 0.387 ms (0.75); F-lowQ Q4 pass 1.42 ms (2.73; geometry only 0.682 ms = 1.31); F-lowQ Q5 pass 1.52 ms (2.93; geometry only 0.742 ms = 1.43); full E2 pass 2.19 ms (4.21); GBT prediction 0.0064 ms per candidate (batch 200; 0.012). So the proxy is 1.4x cheaper than the E2 pass as timed per candidate, but that E2 pass recomputes the Q8 local cost, which inside `rank_key` is free (`c.cost`): a deployable E2 costs about 4.21 - 1.00 = 3.2 equivalents, so the real advantage is about 1.1x, not the 5-10x one might expect from 26 of 112 reflections: most of the time is per-pass Python overhead (peak table, hit lookup, two cost evaluations), which the plan flagged as risk 6 (batching the pass over a level's candidates was not done). Per run, each row with its own scored count (proxy 240 / 384, E2 rerank 246 / 393 candidates): proxy 0.37 s / 0.59 s, E2 as timed 0.54 s / 0.86 s, deployable E2 (without the free Q8 cost) 0.41 s / 0.66 s (`summary.txt` [7]).

**Success criteria.**
- Label validation Spearman >= 0.8: met on the pooled value (0.861); the within-case mean is 0.790 and L2 0.783, so borderline below 0.8 where the ranking actually happens. The classifiers tie the regression offline (D3), so the fallback would not change the offline recall.
- Informational D1 sub-check, not a plan success criterion: best deployable proxy within 0.01 recall of E2. Borderline: gap 0.008 against `e2full/clf3` pooled, but 0.011 pooled against the best E2 variant (`e2full/clf_lm` 0.988), 0.017 clean and 0.006 realistic. The deployable model was chosen as the best of 6 held-out results, which is mildly optimistic.
- Rerank wrong <= E2 rerank (4.0% / 4.5%) at less proxy time than E2: **not met** on the wrong rate (proxy 6.5% / 7.5%; the Wilson intervals overlap; paired exact McNemar on the same runs: discordant 6 vs 1, p = 0.125 clean, 12 vs 6, p = 0.238 realistic, i.e. not significant), met on proxy time (0.37 vs 0.54 s clean, 0.59 vs 0.86 s realistic). Across 3 seeds the proxy is 7.0% / 5.5% wrong.
- (ii) keep 1/8 wrong <= 5% at evaluations <= baseline: not evaluated (D2 did not hold). Row (i) itself fails this test's thresholds (wrong 6.5% / 7.5%; evaluations incl. proxy +3.6% / +5.0%).
- (iii) [cubic-specific: F1 uses the Sigma relatives of cubic symmetry] proxy + F1 wrong <= F1 seed 0 (1.5% / 4.0%) at fewer evaluations than F1b: **met** (0/200 and 4/200 wrong against 3 and 8; evaluations incl. proxy 55,848 / 60,214 against F1b 63,011 / 65,931, i.e. -11.4% / -8.7%; against F1 alone 54,523 / 57,785 it is +2.4% / +4.2% more). The wrong counts are small (a few runs; intervals overlap with F1's and E2 + F1's 4 / 7), so on seed 0 this is a "not worse" result, not a demonstrated improvement over F1. Paired exact McNemar on seed 0, proxy + F1 against F1: discordant 0 vs 3, p = 0.25 (clean) and 3 vs 7, p = 0.344 (realistic); against E2 + F1: 0 vs 4, p = 0.125 and 3 vs 6, p = 0.508. Over three seeds (row (iii) F1 on the seeds 1-2 runs of row (i), `f1_run.py --source-seeds 1 2`; 600 runs per variant): proxy + F1 wrong 0/600 = 0.000 [0.000, 0.006] (clean) and 9/600 = 0.015 [0.008, 0.028] (realistic) against F1 21/600 = 0.035 [0.023, 0.053] and 20/600 = 0.033 [0.022, 0.051]; paired discordant 0 vs 21 (p = 9.5e-07) and 4 vs 15 (p = 0.019). Evaluations incl. proxy 55,894 / 60,266 against F1 54,485 / 57,813 (+2.6% / +4.2%); F1b exists for seed 0 only (63,011 / 65,931), so the 3-seed comparison is against F1 and the F1b comparison stays seed 0. So over three seeds proxy + F1 has a lower wrong rate than F1 at about +3-4% evaluations; the seed-0 result alone was not significant.

**Domain shift check (risk 3).** Harvesting the proxy run's own candidates (all levels, seed 0) and recomputing the recall (L0-2, keep 1/4) with the same fold models: 0.986 (clean) and 0.990 (realistic), against 0.960 / 0.984 for the baseline candidates of the same voxels offline. The recall did not drop, so by the plan's rule (> 0.01) **no retraining** was done. Caveat: the proxy-run recall is conditional on a basin candidate surviving to that level, which can only raise it; it does not capture runs that lost the basin earlier. The end-to-end wrong rate is higher than the E2 rerank although the offline gap is 0.008; compounding the per-level recalls (L0 / L1 / L2 pooled 0.984 / 0.968 / 0.979 for the proxy against 0.996 / 0.975 / 0.982 for E2 full) gives about 0.933 against 0.954, consistent with the observed 6.5% / 7.5% against 4.0% / 4.5%; this arithmetic assumes the levels are independent and ignores level 3 and FindOptimal, so it is only a consistency check and was not tested as the cause.

**Framing (general vs cubic-specific).** The general, symmetry-agnostic result is the proxy rerank: it lowers the wrong rate from 33.5% / 26% to 6.5% / 7.5% at +5% time (the E2 rerank to 4.0% / 4.5%). Proxy + F1 is a cubic Sigma-relative fix (CSL traps matter mainly for cubic), so it is not the headline.

**Findings.** (1) The raw cost is a poor pruning key (recall 0.900 pooled at keep 1/4, 0.799 at 1/8), the raw low-Q cost worse (0.80 / 0.71 for Q5). (2) The untrained hand key (overall hit rates, 0.953) already beats the cost, and the trained low-Q models beat both: F-lowQ Q5 alone 0.962 (0.944 at 1/8), + the free Q8 cost 0.977 (0.961), E2 full 0.985 (0.980). (3) The free post-quick-MC Q8 cost adds about 0.015 recall to the low-Q features (0.962 -> 0.977), more at Q4 (0.947 -> 0.969). (4) Regression on y_bcost and classification are equal in recall to within 0.005; the regressor is better at Spearman against y_bcost (0.65 vs 0.34) but that number is dominated by the near-constant label of the many non-basin candidates. (5) End to end the proxy rerank fixes 55 / 40 of the baseline's wrong runs (67 / 52) and breaks 1 / 3, with 13 / 15 wrong left; combined with F1 it reaches F1's accuracy or better at about F1's evaluations (3 seeds: 0.0% / 1.5% wrong against F1's 3.5% / 3.3%).

**Caveats.** Per-voxel images, not full-sample renders (at most 3 distractor sources). Label leakage: the synthetic candidates and all y_bcost labels use the truth; they are in training only, the metrics of sets A, D, E and the end-to-end runs are harvested (non-synthetic) and set B (synthetic) is reported separately. Domain shift: the model is trained on candidates of baseline searches and used inside a changed search (checked above, no drop by the recall measure). Wall times at 10 workers are contended and are not a speed claim; the evaluation-equivalent accounting uses single-worker medians. The F1 rows' wall time is the evaluation-scaled approximation of `summarize_fixes`. Untested: keep_fraction 1/8 end to end (D2), batching the feature pass, adding the level as a feature, other structures / geometries, full-sample images.

**Deviations from the plan.** (a) Row (ii) not run (D2 did not hold). (b) Q3 dropped (no reflections). (c) The sub-task branch is `feature/nn-hybrid-proxy-profiling-proxy` (the slash form is blocked by an existing ref). (d) `f1_run.py` got `--cache-root` so that `--source` can point to the proxy cache, and `--source-seeds` for seeds 1-2 of a source (defaults unchanged). (e) The relative short refinement is quick MC + one MC of up to 200 steps with `successive_restarts` = 2 restarts (the FindOptimal MC settings without VarianceMinimizing); the truth-basin label uses the full `refine_from_candidates` as the plan says (the short procedure applied to the truth gave exactly the same cost in all 400 cases; the cost is a ratio of integer pixel counts, so this says the counts agree, not that the orientations do). (f) The label-validation draw is 5 candidates per case (stratified by case) rather than a uniform draw over the 360k candidates, restricted to the censored harvested ones.

### Task 3 results: does the network (or the proxy) reduce the total run time? (2026-10-06)

Code: `scripts/profiling/` (`prof_common.py` isolation record, interleaving, stage instrumentation; `prof_u0.py`; `prof_seeded.py`; `prof_profile.py`; `prof_overhead.py`, `prof_summary.py`; `run_all.sh`), tests `tests/test_profiling.py`, results `benchmarks/profiling/` (`summary.txt` / `summary.json` are the full tables; `u0_runs.json.gz`, `seeded_runs.json.gz` the per-run records; `cprofile_*.txt`, `torch_profiler_net.txt`, `mps_vs_cpu.json`). Branch `feature/nn-hybrid-proxy-profiling-profiling` (merged into the feature branch, 8531239). Nothing in `icenine/` changed; all wrappers are installed from the scripts by patching module attributes where they are looked up (`icenine.reconstructor.run_discrete_search_spaced`, `MCOptimizer.optimize` split quick MC / FindOptimal / other MC by `max_mc_steps`, `variance_minimizing_optimize`, `VoxelCostFunction.evaluate` calls and time by role global / local / mc, `perturbation_sweep.prepare_nominal`, `render_windows`, `render_distractor_windows`, `make_realistic_dataset`, `decode_windows`, GN `extract_measurements` and `solve`). `tests/test_profiling.py` checks that one `reconstruct_voxel` run is bit-identical with and without the wrappers (orientation, cost, evaluation counts; the counted `evaluate` calls equal the reconstructor's own counts) and that every patch is removed again.

**Isolation (recorded in `summary.json`, section isolation; raw records in `scripts/profiling/cache`, not committed).** Single worker (spawn pool of 1), `OMP_NUM_THREADS = MKL_NUM_THREADS = 1`, `torch.set_num_threads(1)`, mains power ("AC Power", battery 100% charged), Apple M2 Max (8 P + 4 E cores), no thermal warning level recorded. A background load average of about 1.4 existed before the jobs started (editor, browser, system daemons, probably another agent session). 1-minute load at the start of each kept task (after the re-timing below): U0 median 1.37, maximum 1.92; U1/U2 median 1.41, maximum 1.94 (the timing job itself contributes about 1 to the load). No kept task has a non-ignorable process above 20% CPU in its start record (the busiest are listed in `summary.json`). Before the re-timing, 32 of the 118 tasks (5 of 48 U0, 27 of 70 U1/U2) had a competing process above 20% CPU or a 1-minute load above 2.0 at their start (browser renderers up to 109.6%, `mobileassetd` 99.6%, `mediaanalysisd` 55.7% and 99.4%, `Code Helper (Renderer)` 58.4%, `git` 58.6%, `cloudd` 57.8%); they were re-timed, see the re-timing note below. Discard rule: a task whose 1-minute load at its start was >= about 3.3 was discarded and re-run. Five U1/U2 tasks (loads 3.3 to 10.0: browser renderer, `siriactionsd`) were re-run after the first pass; their original timings were uniformly slower than the re-run: pooled over the five tasks (40 cases per pipeline, 24 for HG) the means were H0 0.547 vs 0.499 s, H1 1.761 vs 1.605 s, H3 2.001 vs 1.815 s, MC told r 1.996 vs 1.829 s, HG 2.183 vs 1.984 s, i.e. every pipeline 9.1-10.2% slower, so the discard favours no pipeline. After review two more tasks that had started next to a `python3` process at 99.3% (v7336 r0 rep1, load 1.84) and `mediaanalysisd` at 139.4% (v16269 r2 rep0, load 2.63) were re-run on a quiet machine (no python / uv process, load 1.35 / 1.40 / 1.84): every pipeline of the two tasks was 9-13% faster than in the original runs (ratios old / new 1.09 to 1.13, again uniform over H0, H1, H3, MC told r), the errors were identical, and the re-run values replaced the originals in all tables here (the medians of the tables changed in the third digit at most, e.g. the r = 0.25 row; the verdict outcomes did not change). Pipelines are interleaved per case in a rotated order (every pipeline runs first equally often); repeats: 3 timings on 4 of the 40 U0 cases (all clean variant; per-case (max-min)/median of the three timings: median 0.712%, 90th percentile 1.78%) and on 10 of the 50 U1/U2 (voxel, radius) tasks, directions 0-1 (40 of the 400 cases per pipeline); per-case medians over the repeats are used in the tables (means of the per-case medians in the U2 verdicts). An untimed warm-up of every pipeline precedes the single-worker runs. Wrapper overhead (`scripts/profiling/prof_overhead.py`, `overhead.json`): 1.71 microseconds per wrapped `evaluate` call on a dummy method, i.e. 47,394 calls x 1.71 us = 0.45% of a 18.1 s baseline run (voxel 7336, clean); a paired run, plain 18.14 / 18.10 s vs wrapped 18.34 / 18.36 / 18.35 s, gives x1.0129 (+1.3%).

**Re-timing (2026-10-06).** *What was flagged.* The owner asked for the remaining affected tasks to be re-timed. A task was flagged if its stored start record (`isolation` in the task's raw file) showed a process other than this job's own python / uv ancestors and `kernel_task` above 20% CPU, or a 1-minute load above 2.0. 5 of 48 U0 tasks (v13884 realistic, v16905 clean, v17451 clean, v18541 realistic, v24550 clean: 10.4%) and 27 of 70 U1/U2 tasks (38.6%; v16182 r1-r4, r8, r9 and r2 rep1, v16269 r0, r1, r3-r6 (r6 also rep1), r8, v2910 r3 (rep0, rep1), r5, r6, r8 (rep0, rep1), r9, v3108 r7, r9 rep2, v7336 r2, r4, r9) were flagged. The rule set before the re-run was to re-time a whole use case if more than 60% of its tasks were flagged; neither reached that, so only the flagged tasks were re-run (estimated about 10 min for U0 and 25 min for U1/U2 from the stored per-task times plus the warm-up; measured 589 s and 1,311 s for the first pass). *How.* The same code, cases, seeds, rotation index, repeats, one worker, `OMP_NUM_THREADS = MKL_NUM_THREADS = torch threads = 1` and warm-up; the only code change is a default-off `--only <task file stems>` option of `prof_u0.py` and `prof_seeded.py`. Before each run the machine was checked (1-minute load < 2.0, no other `uv run` / python job, mains power, no other process above 20% CPU in two `ps` samples 1 s apart) and the run was held until it passed (a browser, `spotlight`, `XprotectService` and photo-analysis daemons delayed the starts, by about 40 minutes for the first U1/U2 pass); the per-task isolation record was stored again. The busy-process check of the first U0 pass was broken (a script error made it pass without evaluating), and in that pass 3 U0 tasks (loads 4.16, 2.73 and 2.56, the last with a browser renderer at 86.5%) and, in the first U1/U2 pass (whose gate worked; background daemons became busy after it passed), 4 U1/U2 tasks (loads 2.47 and 2.62; `Code Helper (Renderer)` 66.7% with `WindowServer` 52.6% and `airportd` 43.7%; `corespotlightd` 50.8%; `dasd` 24.3%) were flagged again by the same rule. They were re-run in a second pass (3 U0 and 4 U1/U2 tasks, task-start loads 1.18 to 1.37) and the 2 tasks flagged again there (U0 v17451 with `networkserviceproxy` 56.2%, U1/U2 v2910 r8 rep1 with `sysmond` 26.8%) in a third (task-start loads 1.22 and 1.51). The quiet gates passed at loads 1.38 (pass 1 U1/U2), 1.13 and 1.29 (pass 2), 0.95 and 1.78 (pass 3). 9 U0 and 32 U1/U2 task runs were made in all; the timings of the last run of each flagged task replaced the originals in the stored runs (old ones kept in the gitignored `cache/*/retimed_old/`). *Reproduction.* U0: R_final of every run is bit-identical to the old run and to the stored earlier-task run (`=stored`, 5 tasks x 6 pipelines). U1/U2: R_final, error, converged and fall-back flags are identical to the old run for all 196 cases of H0, H1, H3 and MC told r and 84 of HG (so the H3 / HG wrong flags agree in 100% of the re-timed cases; the rates against the full Task 1 runs are the 98.5% and 99.4% of the Checks paragraph, unchanged). *Old / new.* Median over cases of old time / new time: U0 baseline 0.967, F1 0.984, F1b 0.972, E2 0.974, proxy 0.968, proxy + F1 0.975 (per-task median old/new 0.954 to 0.999, i.e. the re-run U0 tasks took 0.1% to 4.8% longer than the originals); U1/U2 H0 1.011, H1 1.008, H3 1.009, MC told r 1.013, HG 1.003 (per task 0.944 to 1.054; per case 0.83 to 1.12). The change is small and about the same for all pipelines (within 1.7% for U0 and 1% for U1/U2), so the original timings of the flagged tasks were not biased toward any pipeline; the direction differs between U0 (new slower) and U1/U2 (new slightly faster) and the cause was not tested. *Effect.* The U0 medians and means moved by at most 0.5%; the U1/U2 band medians by at most 1.8% and the per-radius medians by up to 3.9% (H0 at r = 0.75, 0.015 s); the 10-90% ranges moved more; the contention factors did not change (U0 1.25, U1/U2 1.21). The largest changes are U1 clean H1 vs H0 +5.8% to +7.2% and H3 vs H0 +9.7% to +11.1%, the H0 + fallback comparator rows (+0.1 s, through the baseline reconstruction 21.8 s to 21.9 s realistic), and the U0 F1 and proxy medians by 0.1 s. **No verdict changed** (U0, U1, U2 pooled, U2 per radius including r = 5, and the cubic F1b-fallback comparator rows are the same true / false as before; the numbers in this section follow the regenerated `summary.txt` / `summary.json`).

**Checks.** Every U0 run reproduces its stored earlier run bit for bit: baseline (E0), F1b, E2 rerank, proxy row (i) and the two F1 post hoc steps (28/28 clean and 20/20 realistic runs per pipeline, repeats included; `=stored` column in `summary.txt`), so the accuracy of the U0 timing subset (20 voxels x 2 variants, the first 20 of the 200 voxels) is exactly that of the stored runs on those voxels. U1/U2: H0 reproduces the Task 1 error on all 400 cases (identical), H1 on all 400; H3 and HG are exactly equal on 50.5% and 50.6% of their cases (passes 2-3 are re-rendered on a batch of one with another noise stream, so only pass 1 is identical) with the same wrong flag in 98.5% and 99.4%; MC told r (one-shot optimizer sweep, stored unreduced) does not reproduce its stored result: exactly equal in 9.5% of the cases (median |difference| 0.000441 deg) and the wrong flag differs in 8.2% of the cases (it is timed, but its accuracy comes from the stored sweep, not from the timing run). The accuracy columns of the tables below ("full run") come from the full Task 1 runs (all 50 voxels x 20 directions) and the Task 2 end-to-end summary, not from the small timing subset.

**U0 (no start), single worker, wall seconds per voxel, median [10-90%], 20 voxels per variant, seed 0.** F1 and proxy + F1 are the parent run plus the post hoc step (same run).

| pipeline | clean s | realistic s | evaluations clean (global / local / all incl. proxy and F1) | wrong clean / realistic on the 20 voxels | wrong in the 200-voxel run (Task 2) clean / realistic |
|---|---|---|---|---|---|
| baseline | 19.8 [19-20.9] | 21.7 [19.9-24] | 4.45e4 / 4.61e3 / 4.91e4 | 9/20, 7/20 | 33.5%, 26% |
| F1 (Sigma <= 29) | 23 [22.5-23.6] | 24.5 [23.2-27.2] | 4.45e4 / 1.04e4 / 5.49e4 | 0/20, 2/20 | 1.5%, 4% |
| F1b | 27.5 [27-28] | 28.9 [27.9-31.2] | 4.55e4 / 1.78e4 / 6.33e4 | 0/20, 0/20 | 0%, 1.5% |
| E2 rerank | 21.2 [20.2-21.9] | 23.1 [21.5-27.4] | 4.45e4 / 5.96e3 / 5.1e4 | 1/20, 2/20 | 4%, 4.5% |
| proxy rerank (lowq5+c8 regression) | 20.8 [20-21.8] | 22.9 [20.6-25.7] | 4.45e4 / 5.73e3 / 5.07e4 | 1/20, 3/20 | 6.5%, 7.5% |
| proxy + F1 | 23.5 [22.9-24.4] | 25.8 [23.8-29.3] | 4.45e4 / 1.07e4 / 5.57e4 | 0/20, 0/20 | 0%, 2% (3 seeds 0%, 1.5%) |

(The "all incl. proxy" column counts the cost-function calls of the proxy's low-Q features as well; the reconstructor's own counts are the global and local columns. The 10-worker wall times of Task 2 were 24.5 s baseline and 27.9 s proxy + F1 clean: single-worker times are 16-19% lower, see the contention factor.)

**U0 stage breakdown** (mean exclusive seconds per run, clean; realistic in `summary.txt`). Baseline 19.8 s: cost evaluations 17.1 s global (86.2%) + 2.39 s local (12.1%); everything else together (discrete-search loop, quick-MC bookkeeping, FindOptimal and VarianceMinimizing code outside `evaluate`) 0.34 s (1.7%). Realistic baseline: 17.3 s global (79.0%) + 4.15 s local (19.0%). F1b: 17.5 s global + 9.21 s local (33.5%). Proxy: set_image + features + prediction 5.95e-06 + 0.183 + 0.00371 s per run (clean), 0.308 s realistic; E2 0.277 s / 0.478 s. The F1 post hoc step costs 3.19 s clean / 3.05 s realistic (94% of it in local evaluations).

**U0 verdicts** (plan rule: helps vs a comparator if the median wall time is >= 20% lower at matched accuracy, or the wrong rate is lower with non-overlapping Wilson intervals at <= 10% extra time; matched = wrong rate <= comparator's upper bound and median error of right answers within +0.005 deg; seed-0 n = 200 accuracy, 3-seed n = 600 for proxy + F1 against baseline and F1; judged per variant, "helps" only if it holds in both). **Headline: no U0 pipeline reduces the run time.** Every alternative to the baseline costs time (+5% to +39%); the only "helps" verdicts are accuracy gains at <= 10% extra time. Framing: F1, F1b and proxy + F1 are cubic-specific fixes (they use Sigma relatives); the symmetry-agnostic U0 result is the rerank against the baseline: the E2 and proxy rerank lower the wrong rate from 33.5% / 26% to 4% / 4.5% and 6.5% / 7.5% at +5% to +7% time, and nothing is faster. F1b is not used as the accuracy target of this headline.
- E2 rerank vs baseline: **helps (accuracy, +6.6% to +7.0% time; not faster)**: wrong 4% [2.04, 7.69] / 4.5% [2.39, 8.33] vs 33.5% [27.3, 40.3] / 26% [20.4, 32.5], time +7.0% / +6.6%. Proxy rerank vs baseline: **helps (accuracy, +5.2% to +5.5% time; not faster)**: wrong 6.5% [3.84, 10.8] / 7.5% [4.6, 12] at +5.2% / +5.5%.
- Cubic-specific rows: F1 vs baseline: not (+16.0% / +13.0% time, above the 10% allowance). F1b vs baseline: not (+38.9% / +33.3%). F1b vs F1: not (+19.7% / +17.9% time for 1.5% vs 0% and 4% vs 1.5% wrong with overlapping intervals). Proxy + F1 vs F1 (3 seeds): clean 0% [0, 0.636] vs 3.5% [2.3, 5.29] at +2.4% time (helps), realistic 1.5% [0.791, 2.83] vs 3.33% [2.17, 5.09] at +5.3% (intervals overlap: not); overall **not** (it needs both variants). Proxy + F1 vs baseline: not (+18.8% / +19.0% time). Proxy + F1 vs F1b: not: it is 14.4% / 10.7% faster (23.5 vs 27.5 s clean, 25.8 vs 28.9 s realistic) but the median error of right answers is 0.0282 / 0.0296 vs 0.0200 / 0.0176 deg, outside the +0.005 deg tolerance (wrong 0% vs 0% clean, 2% vs 1.5% realistic).
- No pipeline is >= 20% faster than any comparator at matched accuracy. Per-variant details and all other pairs: `summary.txt`.


**U1 / U2 (a start exists), single worker, production-equivalent wall seconds per case** (400 cases per pipeline: first 5 sweep voxels x 10 radii x 4 directions x 2 variants; HG at r = 1.5, 2, 3, 5 only; both variants pooled; median [10-90%]):

| r (deg) | H0 | H1 | H3 | MC told r | HG |
|---|---|---|---|---|---|
| 0.05 | 1.3 [0.993-2.99] | 1.45 [0.991-2.55] | 1.66 [1.14-2.16] | 1.99 [1.97-2.02] | - |
| 0.1 | 1.42 [0.725-2.43] | 1.61 [0.987-2.86] | 1.5 [1.06-2.98] | 1.98 [1.96-2.01] | - |
| 0.25 | 1.14 [0.65-1.87] | 1.32 [0.861-2.96] | 1.44 [0.904-2.22] | 1.96 [1.84-2] | - |
| 0.5 | 0.41 [0.251-1.17] | 1.36 [0.994-2.63] | 1.44 [0.881-2.17] | 1.97 [1.93-2.03] | - |
| 0.75 | 0.374 [0.259-1.12] | 1.57 [0.855-2.21] | 1.65 [0.954-3] | 1.98 [1.93-2.04] | - |
| 1 | 0.327 [0.257-1.17] | 1.4 [0.902-2.84] | 1.63 [0.975-2.4] | 1.99 [1.86-2.02] | - |
| 1.5 | 0.303 [0.241-1.11] | 1.38 [0.822-2.16] | 1.52 [0.985-2.64] | 1.98 [1.86-2.03] | 1.81 [1.16-2.75] |
| 2 | 0.248 [0.237-0.334] | 1.11 [0.496-2.25] | 1.43 [0.927-2.53] | 1.98 [1.9-2.01] | 1.65 [1.24-4.08] |
| 3 | 0.242 [0.239-0.268] | 1.07 [0.3-2.13] | 1.55 [0.919-2.39] | 1.97 [0.631-2.01] | 1.76 [1.36-2.65] |
| 5 | 0.189 [0.175-0.248] | 0.307 [0.221-1.28] | 1.17 [0.231-2.17] | 1.99 [0.0326-2.05] | 1.83 [0.229-3.73] |

Use-case bands (median [10-90%] production time s; mean evaluations; wrong in the full Task 1 run with Wilson 95% interval and median error in deg; subset wrong in `summary.txt`):

| band, variant | pipeline | time s | evals | wrong (full run) | median err (full run) |
|---|---|---|---|---|---|
| U1 (r 0.05, 0.1) clean | H0 | 1.32 [0.981-1.72] | 2.31e3 | 0% [0, 0.192] | 0.0227 |
| | H1 | 1.41 [1.07-1.9] | 2.52e3 | 0% | 0.011 |
| | H3 | 1.46 [1.12-1.99] | 2.42e3 | 0% | 0.0105 |
| | MC told r | 1.99 [1.96-2.02] | 3.33e3 | 0% | 0.0381 |
| U1 realistic | H0 | 1.67 [0.734-3.44] | 3.27e3 | 0% | 0.0219 |
| | H1 | 1.61 [0.95-3.66] | 3.39e3 | 0% | 0.0205 |
| | H3 | 1.77 [1.05-3.37] | 3.19e3 | 0% | 0.0201 |
| | MC told r | 1.98 [1.96-2.02] | 3.5e3 | 0% | 0.0381 |
| U2 (r 0.5-3) clean | H0 | 0.285 [0.24-0.766] | 722 | 28% [26.8, 29.1] | 0.19 |
| | H1 | 1.36 [0.899-1.95] | 2.37e3 | 0% [0, 0.064] | 0.0164 |
| | H3 | 1.53 [0.986-1.97] | 2.41e3 | 0% | 0.0105 |
| | HG (r 1.5-3) | 1.69 [1.13-2.09] | 2.34e3 | 0% [0, 0.128] | 0.0105 |
| | MC told r | 1.97 [1.9-2.02] | 3.37e3 | 37.6% | 0.659 |
| U2 realistic | H0 | 0.268 [0.24-1.2] | 898 | 26.4% [25.3, 27.5] | 0.183 |
| | H1 | 1.34 [0.34-2.93] | 2.66e3 | 3.42% [2.99, 3.91] | 0.0265 |
| | H3 | 1.61 [0.882-3.07] | 3.11e3 | 0.7% [0.518, 0.945] | 0.0208 |
| | HG (r 1.5-3) | 1.96 [1.33-4.08] | 3.43e3 | 0% [0, 0.128] | 0.0208 |
| | MC told r | 1.98 [1.87-2.03] | 3.35e3 | 37.5% | 0.651 |
| U2 stress (r 5) clean | H0 | 0.177 [0.175-0.244] | 339 | 100% | 23.9 |
| | H1 | 0.741 [0.22-1.34] | 1.35e3 | 0% [0, 0.407] | 0.0365 |
| | H3 | 1.32 [0.216-2.11] | 1.94e3 | 0% | 0.0107 |
| | HG | 1.68 [0.224-2.24] | 2.01e3 | 0% | 0.0105 |
| U2 stress (r 5) realistic | H0 | 0.238 [0.175-0.255] | 381 | 100% | 6.02 |
| | H1 | 0.295 [0.224-0.318] | 403 | 95.2% | 4.3 |
| | H3 | 0.415 [0.341-2.53] | 1.64e3 | 43.6% [40.5, 46.7] | 0.175 |
| | HG | 2.11 [1.45-4.52] | 3.8e3 | 2.31% [1.55, 3.45] | 0.0199 |

(The full-run accuracy is over the pass-1 successes of the sweep; at r = 5 the clean variant has 60 and the realistic variant 6 of 1000 fall-back cases in which the network has no estimate: they are wrong for every pipeline and are not in these wrong rates. In the timing subset they are 5 of 20 clean cases at r = 5, which is why its clean wrong counts for H1 / H3 / HG are 5/20.)

**Network time per case** (median over all cases; measured = including the harness rendering of every pass, production-equivalent = prepare_nominal (ROI set, observer, window spec), `decode_windows` and the forward pass only): H1 0.0568 vs 0.0498 s, H3 0.168 vs 0.148 s, HG (GN included) 0.378 vs 0.35 s. Stage breakdown, mean exclusive seconds per case: H3: prepare 0.132, decode 0.0009, forward 0.0137, harness render 0.0220 (render 0.0005 + physics 0.0093 + distractors 0.0100 + realism 0.0023), finisher cost evaluations 1.42 s (evaluate_local) + 0.0776 VarianceMinimizing + 0.0042 FindOptimal bookkeeping; H1: prepare 0.0445, forward 0.0048, evaluate_local 1.33; H0: evaluate_local 0.706; HG: prepare 0.258, GN solve 0.0627 + extract 0.0027, forward 0.0134, evaluate_local 1.55; MC told r: evaluate 1.80 s. The pass-1 windows of every net pipeline are taken from the sweep batch (stand for the experimental windows) and are not charged.

**U1 verdicts.** H1 vs H0: not (clean: +7.2% time, median error 0.011 vs 0.0227 deg; realistic: -3.4% time, median error 0.0205 vs 0.0219 deg, wrong 0% vs 0%). H3 vs H0: not (+11.1% / +6.3% time). Under the literal rule (any lower median error counts as "beats on accuracy", without the 0.002 deg tie used here) U1 realistic H1 would be "helps" (-3.4% time, 0.0205 vs 0.0219 deg) but that -3.4% is on 40 cases with a 10-90% range of 0.95-3.66 s, i.e. within noise, and the clean variant fails on time (+7.2%); so the overall verdict stays no and H0 remains the recommended finisher for a start within 0.1 deg. MC told r (1.99 s, median error 0.0381) is slower and less accurate than H0.

**U2 verdicts.** **The networks help only where H0 fails (r >= 1.5 deg); for a start within 1 deg H0 is faster at an equal wrong rate (but with a 2.6-4.5x higher median error), so U2 is a conditional result: the pooled r = 0.5-3 verdict depends on the uniform weighting of the radii 0.5, 0.75, 1, 1.5, 2, 3 in the sweep design.** Per radius (means of per-case times, ideal-trigger H0 + fallback = H0 mean time + H0's wrong rate x the baseline reconstruction's mean time (19.8 s clean, 21.9 s realistic); the fallback answer is a full reconstruction, itself wrong 33.5% / 26%):
- r = 0.5, 0.75, 1: H0 + fallback takes 0.44-0.81 s (wrong 0% to 0.34%) against 1.44-2.43 s for H1 and H3, so the networks are 1.96-3.4x slower (speed ratio 0.30-0.51x clean, 0.34-0.47x realistic; H3 0.30x / 0.47x at r = 0.5, 0.47x / 0.36x at r = 1) at an equal wrong rate (H0's median error is 2.6-4.5x higher: 0.055-0.090 deg against H3's about 0.020 deg, and its fraction < 0.1 deg is 0.53-0.72 against 0.97-0.98): no help by the plan's rule. H0's median time at these radii is 0.33-0.41 s.
- r = 1.5: H0 + fallback 4 / 4.05 s (wrong 6.03% / 4.13%, H0 itself 18% / 15.9%): H1 1.33 / 1.57 s (3x / 2.58x faster), H3 1.48 / 1.8 s (2.71x / 2.25x), HG 1.64 / 2.13 s (2.43x / 1.9x), all with 0% wrong: helps. r = 2: H0 + fallback 10.7 / 10.5 s: H1 7.94x / 8.49x, H3 7.4x / 5.27x, HG 6.56x / 3.96x: helps. r = 3: H0 + fallback 19.3 / 21 s; H3 12.7x / 11.9x faster with 0% / 4.1% wrong and H1 14.4x / 22.4x with 0% / 19.4% wrong (better than the fallback's 32.2% / 24.6% but not accurate), HG 11.3x / 9.35x with 0% / 0%: by the plan's rule all three pass (the rule compares with the best non-NN option), but only HG is accurate on realistic data.
- Pooled r = 0.5-3 (H1, H3), comparators H0 (mean 0.413 / 0.519 s, wrong 28% / 26.4%), the baseline reconstruction and H0 + fallback (5.95 / 6.3 s, wrong 9.37% / 6.86%): H3 mean 1.51 / 1.92 s, wrong 0% / 0.7% [0.518, 0.945], 3.94x / 3.28x faster than H0 + fallback, 13.1x / 11.4x faster than a full reconstruction: **helps** (symmetry-agnostic comparators). H1: mean 1.39 / 1.57 s, wrong 0% / 3.42%, 4.27x / 4.01x, 14.2x / 13.9x: helps. HG (timed at r = 1.5-3, comparators recomputed on those radii: H0 mean 0.321 / 0.358 s, wrong 55.6% / 52.4%, H0 + fallback 11.3 / 11.8 s, wrong 18.6% / 13.6%): HG mean 1.66 / 2.34 s, wrong 0% / 0%, 6.83x / 5.06x faster than H0 + fallback, 11.9x / 9.35x faster than full reconstruction: helps.
- Cubic-specific comparator: a fallback to F1b instead of the baseline reconstruction (H0 + F1b fallback, ideal trigger; F1b mean 27.5 / 29.5 s, wrong 0% / 1.5% from Task 2) takes 8.11 / 8.32 s (pooled r = 0.5-3) with wrong 0% / 0.396%, so H3 is 5.37x / 4.33x faster. The "only method at wrong <= 1%" route then fails (the fallback is at <= 1%), and the >= 1.5x route still holds, but the accuracy-match test also fails for H1 (3.42% wrong) and for realistic H3 (0.7% wrong vs an upper bound of the 0.396% fallback): with this comparator the pooled r = 0.5-3 verdict becomes **not** for H1 and H3 (clean passes, realistic fails); HG on r = 1.5-3 still helps (0% / 0% wrong; 9.43x / 6.77x faster). The fallback is an ideal oracle in both cases; a real trigger would favour the networks more.
- Stress r = 5: only **HG helps** (wrong 0% / 2.31%; mean 1.41 / 2.57 s; 14.1x / 8.6x faster than H0 + fallback; 14x / 8.51x faster than full reconstruction). H1 and H3 do not (realistic wrong 95.2% and 43.6%); their clean results use the full-run wrong rates over pass-1 successes and exclude the 6% clean fall-back cases (wrong for every pipeline), so with those counted none of the pipelines is at <= 1% wrong in the clean variant.
- Speed-up vs full reconstruction (mean per-voxel time of the baseline over the mean per-case time): H0 is roughly 25-100x faster but wrong at r >= 1.5; H1 and H3 reach 11-23x and HG 8.3-12x (r = 1.5-3, both variants) and the hybrids use 2.3e3-3.4e3 cost evaluations per case vs 4.9e4-5.3e4 for the full reconstruction (about 15-20x fewer). Means are used for the U2 verdicts and these ratios (S1 of the review: H0 + fallback is an expectation, so the pipelines are compared as means); the tables report medians.


**Profilers.** `cprofile_u0_*.txt`, `cprofile_seeded_*.txt` (top 30 by cumulative time, 5 runs per pipeline; function names without directories). Share of the total self time (tottime) under cProfile, baseline: icenine/ Python code 60.5%, C builtins and methods 28.9%, scripts 5.31%, third-party Python 4.9%; the seeded pipelines are within 58-62% / 27-28% / 6.9-8.0% / 3.6-4.5% (H3: 60.8% / 27.7% / 7.54% / 3.73%). cProfile adds a per-call cost, so the Python share is overstated; the stage timer (no per-line tracing) puts 98.3% of a baseline run (clean) in `VoxelCostFunction.evaluate`. The "physics" and the Python overhead inside `evaluate` were not separated. Outside `evaluate` the search machinery (discrete-search loop, MC and FindOptimal bookkeeping) is 1.7% of a baseline run and the `rank_key` proxy 0.9-2.0% (0.18-0.48 s per run). Top cumulative entries of every pipeline are `evaluate` and its callees. `torch.profiler` (CPU, net x3 stage, 6 cases, realistic and clean, radius 2 deg; `torch_profiler_net.txt`), record_function ranges, total over 6 cases: net 1510 ms, prepare 1170 ms (18 calls, 65.0 ms per call), render (harness) 241 ms (distractors 148, physics 75, realism 15), forward 97.2 ms (18 calls, 5.40 ms per call), decode 7.87 ms (0.437 ms per call).

**Forward pass, MPS vs CPU** (realistic_s0, torch 2.8.0, median ms per call [10-90%], 200 repeats): batch 1: CPU 1 thread 4.4 [4.38-4.42], CPU 8 threads 4.65 [4.58-4.78], MPS resident 3.8 [3.66-3.93], MPS with host-device transfer 4.46 [4.43-4.55]; batch 20: CPU 1 thread 77.6 [77.3-77.9], CPU 8 threads 75.1 [74.5-75.5], MPS resident 3.4 [3.24-3.68], MPS with transfer 8.3 [7.85-8.42]. Outputs agree (max |difference| 7.75e-07 / 9.54e-07 deg in the mean, 4.47e-08 / 2.53e-07 in the Cholesky factor). The forward pass is 0.0048-0.0138 s of a 1.4-1.7 s hybrid case, so the batch-20 MPS speed-up over one CPU thread (22.8x resident, 9.3x with host-device transfer; batch 1: 1.16x resident, 0.99x with transfer) would not change the run time of a single case.

**Contention (10 workers vs 1 worker, same cases; median of per-run time ratios).** U0: 1.25 (20 voxels x 6 pipelines; per pipeline 1.24-1.25, 10-90% 1.18-1.29); U1/U2: 1.21 (880 matched pipeline runs; per pipeline 1.20-1.22, 10-90% 1-1.32). A 10-worker run therefore delivers about 10 / 1.25 = 8.0x (U0) and 10 / 1.21 = 8.3x (U1/U2) the single-worker throughput of the timed portion. The 10-worker wall: 288 s for the 20 U0 tasks (20 x 6 pipelines), 185 s for the 50 U1/U2 tasks (2 directions).

**Caveats.** Per-voxel images (at most 3 distractor sources), not full-sample renders; detector-window rendering is harness work and excluded from "production" times (the pass-1 windows are not charged at all); the image attach (`attach_images`) is excluded for every pipeline; the F1 post hoc step is the `f1_run.py` procedure exactly (quick MC of all Sigma <= 29 relatives, FindOptimal on the union of the top 3 of Sigma <= 11 and of Sigma <= 29), so it refines up to six candidates although the Sigma <= 29 answer only needs the top 3 (the same code gave the accuracy it is paired with); H3 and HG passes 2-3 re-render on a batch of one with a different noise stream than Task 1's batch of 20; the proxy and E2 models are those of Task 2; the "H0 + fallback" trigger is ideal; the timing subsets are small (20 voxels per variant for U0, 5 voxels for U1/U2) so the wrong rates of the timing subset are only a consistency check; Apple P/E cores, thermals and background daemons are not controlled (load recorded above); one machine; the five load-spike re-runs and the 32 re-timed tasks ran at different times than the rest.

**Deviations from the plan.** (a) "Proxy keep 1/8" is dropped from U0 (not run in Task 2, D2 failed). (b) U2's "full reconstruction ignoring the start" is not timed separately: its time and accuracy are the U0 baseline's (the start is ignored, so the cost does not depend on r). (c) H3m is not in the table; "MC told r" is `optimizer_sweep.run_mc_adam` with the true r from the nominal. (d) H1 is kept as the strictly paired hybrid. (e) The repeats are on 4 clean U0 cases only (the first of every ten). (f) F1 and proxy + F1 are not run as independent pipelines but as a post hoc step on the parent run in the same process right after it (their time is the parent's plus the step's); so F1's and baseline's times are not independent samples. (g) Five U1/U2 tasks were re-run after a load spike (see isolation), and 32 tasks (5 U0, 27 U1/U2) with a competing process or a load above 2.0 at their start were re-timed (see the re-timing note). (h) cProfile is run on 5 runs per pipeline at radii 0.1, 0.5, 1.5, 2, 3 (HG at 1.5, 2, 3, 5, 1.5), voxel 0, direction 0, variants alternating.


### Completion summary (2026-10-06)

All three tasks are merged into `feature/nn-hybrid-proxy-profiling`. Nothing in `icenine/` changed; the new code is in `scripts/nn_hybrid/`, `scripts/coarse_proxy/` and `scripts/profiling/` (tests `test_nn_hybrid.py`, `test_coarse_proxy.py`, `test_profiling.py`, including bit-identity checks against the stored runs). Numbers below are from `benchmarks/nn_hybrid/summary*.txt|json`, `benchmarks/coarse_proxy/summary.txt|json` (and `eval_split.json`) and `benchmarks/profiling/summary.txt|json` (and `overhead.json`). All results use per-voxel images (at most 3 distractor sources), not full-sample renders. "Realistic" data = neighbour/twin overlap plus detector noise.

**Decisions**

| Item | Decision | Reason |
|---|---|---|
| Net passes (D1, Task 1) | keep net x3 | x1 is equivalent for r <= 1.5 but worse at r = 2 (1.10% vs 0.10% wrong) and r = 3 |
| H3c, covariance-sized FindOptimal box (D3) | dropped | worse than H3 where the box exceeds 0.329 deg (median 0.0494 vs 0.0339 deg) |
| H3m, net -> MC finisher (D4) | dropped | median +0.009 to +0.019 deg worse, time 0.75-0.80 of H3, not under half |
| sigma_max-triggered HG fallback | no gain | no better than using plain HG whenever a large start is possible |
| Regression vs classifier label (Task 2 D3) | regression carried | recall 0.977 vs 0.977 |
| keep 1/8 end to end (Task 2 D2) | not run | offline recall 0.963 < 0.98 |
| Q3 features | dropped | no reflection with \|q\| <= 3 |
| Proxy retraining | none | no domain-shift drop in recall (0.986 / 0.990) |
| Sub-branch naming | `-task` instead of `/task` | git cannot hold `feature/x/task` with `feature/x` |
| `icenine/` | unchanged | all instrumentation is by patching from the scripts |

**Headline results**

| Task | Result |
|---|---|
| 1. Hybrid net -> FindOptimal | H3 median about 0.020 deg for r <= 3 (net alone 0.062-0.080 deg). Wrong rate (realistic): H3 4.1% at r = 3 and 43.6% at r = 5; HG 0% for r <= 3 and 2.3% at r = 5. Criterion C2 fails at r = 3 in both seeds. On clean data FindOptimal adds nothing after the net (net alone 0.011 deg). About 96% of finisher results have a higher cost than the truth (cause not diagnosed). |
| 2. Q_max-8 cost proxy | Offline pruning recall (keep 1/4, L0-2, pooled): cost 0.900, hand-made key 0.953, low-Q5 + cost8 0.977, E2 0.985. End to end, proxy rerank 6.5% / 7.5% wrong against E2 rerank 4.0% / 4.5% (not significant, McNemar p = 0.125 / 0.238; baseline 33.5% / 26%); proxy scoring is about 1.1x cheaper than a deployable E2. Proxy + F1 over 3 seeds: 0% / 1.5% wrong against F1 3.5% / 3.3% (cubic-specific, Sigma relatives). |
| 3. Run-time profiling | U0: no pipeline is faster; alternatives cost +5% to +39%, so the reranks are accuracy-only wins (+5.2% to +5.5% proxy, +6.6% to +7.0% E2). U1: use H0. U2: the hybrids help only at r >= 1.5 (H3 3.3-3.9x faster than H0 + ideal fallback pooled over r = 0.5-3, about 12-14x faster than a full reconstruction); at r = 5 only HG helps. 98.3% of a baseline run is in `VoxelCostFunction.evaluate`. MPS does not matter per case. Contention factor (10 workers vs 1) 1.21-1.25. |

**Recommendations**
- No start (U0): use a symmetry-agnostic rerank, the proxy or the E2 rerank (accuracy gain at +5% to +7% time, no speed-up). F1 / proxy + F1 are cubic-specific add-ons.
- Start within 0.1 deg: H0 (FindOptimal alone).
- Start that may be 1.5 deg or more off: HG (or H3 when r <= 3 and a few percent wrong at r = 3 is acceptable); for r = 5 only HG.
- Speed work belongs in `VoxelCostFunction.evaluate` (98.3% of the time), not in the network (0.15-0.38 s of a 1.4-2.1 s hybrid case).

**Open follow-ups**
- Batch the low-Q feature pass (most of its time is per-pass Python overhead).
- keep 1/8 end to end (not run, D2).
- Why the finisher stops above the truth's cost (about 96% of results).
- Full-sample renders (denser realism).
- A real fallback trigger for H0 (the study used an ideal oracle).
- Candidate de-duplication (`docs/todo_coarse_cost_proxy_and_dedup.md`, Idea 2).
- The three ideas in `docs/todo_future_ideas_nn_active_fourier.md` (no-start NN, active imaging, Fourier / resolution theory).
- The level as a proxy feature.
- Adopt the shared stats / doc_tables / preflight helpers from `feature/dev-tooling` once it merges.

## Developer tooling (2026-10-06)

**Plan (short form).** Turn rules and checks that were repeated in LLM prompts into code that runs without an LLM, and shrink the LLM prompts to planning, interpretation, logic review and writing. Deliverables: (A) a Bash PreToolUse hook script (bare python/pytest/pip, blanket `git add`/`git commit -a`); (B) a committed but not activated git pre-commit hook with Python-implemented checks; (C) task-branch scripts and skills (`feature/<parent>-<task>` naming, because git cannot hold `feature/x/task` beside `feature/x`); (D) `scripts/common/stats.py`; (E) doc-table generation and sync; (F) a number audit; (G) a timing preflight; (H) job status and checkpoint; (I) a TODO scaffold; (J) an implementer agent and a claims-audit section for the reviewer; (K) README, this section and CLAUDE.md. Built in a separate worktree (`IceNine-tooling`) while a single-worker timing job ran in the main tree; targeted tests only under `nice` while it ran.

### Completion summary (2026-10-06)

All of A-K implemented on `feature/dev-tooling`; tests in `tests/test_dev_tooling.py` (99 passed after the review fixes; full suite 639 passed, 34 skipped with the main venv). Scripts are listed in the README section "Developer tooling". The git hooks are not activated; the main session runs `scripts/dev/install_hooks.sh` after the merge (and after the other feature's agent has committed, since `core.hooksPath` is shared by all worktrees). Not changed: `scripts/nn_hybrid`, `scripts/coarse_proxy`, `scripts/profiling`, `icenine/`; adopting `stats`, `doc_tables` and `preflight` there is a later step.

**Number audit on the FindOptimal robustness section above.** `audit_numbers.py` extracts 128 numbers and matches 125 against `benchmarks/findoptimal_robustness/*`; the three unmatched are real non-source numbers (the test count 537, the derived 853 right runs, and the derived 347/347). The calibration line shows the limit (reported per precision bucket: 0 dp 70%, 1 dp 91%, 2 dp 86%, 3+ dp 89%): shifting every number by 3 in its last printed digit still matches about 80% overall, because the source files contain thousands of values and low-precision decimals (0.2, 0.94) match almost anything. The audit catches wrong-magnitude and invented numbers, not small transcription errors; counts and high-precision values are the strongest checks.

## Follow-ups (2026-10-06)

Plan (short form). Branch `feature/followups` off develop; one sub-task branch per task
(`scripts/dev/start_task.sh <name>`), merged back after review. Results go in a subsection per task.
FindOptimal F1/F1b are cubic-specific.

- **T1 `env`**: reproducible environment for the golden tests. A fresh `uv pip install` picked
  Python 3.12 / numpy 2.4 / torch 2.10 and 4 golden bit-identity tests in
  `tests/test_findoptimal_refactor.py` fail at 1e-9. Find where outputs diverge (float drift vs
  branch flip); make `uv sync` from `uv.lock` the documented setup (`.python-version` 3.9); verify in
  a temporary clone; record numpy/torch/python versions next to the golden data and fail clearly on
  mismatch (no loosened tolerance, no silent skip).
- **T2 `stats-adoption`**: use `scripts/common/stats.py` (and `preflight.py` where the record format
  stays readable) in the study scripts; regenerated summaries must be byte-identical.
- **T3 `lowq-batch`**: batch the F-lowQ Q5 feature pass over candidates and peaks; identical
  features (<= 1e-12); time per candidate at batch sizes 1, 50, 200 on a single worker; report the
  new proxy cost in evaluation equivalents.
- **T4 `keep-eighth`**: proxy rerank with keep_fraction 1/8 (and 1/6) end to end, 200 voxels x 2
  variants, seed 0, paired with baseline and keep-1/4 proxy; then a 20-voxel single-worker timing.
  "Useful" means wrong <= keep-1/4 proxy Wilson upper bound at fewer total evaluations than baseline.
- **T5 `finisher-diagnosis`**: why the finisher stops above the truth's cost (path barrier, longer or
  finer continuation, pixel granularity, stopping rule). Observations only; no reconstructor change.

Not in this round: full-sample renders and the ideas in `docs/todo_future_ideas_nn_active_fourier.md`.

### T1 `env`: locked environment for the golden tests (2026-10-06)

**Where the golden outputs diverge.** One golden case (ThreeVoxels voxel 0, `reconstruct_voxel`, seed 7)
was run with the recorder attached in the main environment (Python 3.9.6, numpy 2.0.2, torch 2.8.0,
scipy 1.13.1) and in a fresh `uv pip install -e ".[dev]"` environment (Python 3.12.12, numpy 2.5.3,
torch 2.14.1, scipy 1.18.1; the resolver picks newer versions than the 3.12/2.4/2.10 first seen).
Every stage event has the same sequence, the same costs (maximum cost difference 0.0, final cost
0.8180327868852459 in both) and the same discrete-stage orientation (difference 0). Rotation matrices first
differ in the quick-MC stage (maximum entry difference 6.1e-9) and carry through to the final
orientation (3.8e-7 deg in the third rotation-vector component, against the 1e-9 test tolerance). So this
case shows small floating-point drift (about 6e-9 in matrix entries), not a branch flip in the MC
search. The source of the drift was not isolated; `scipy.spatial.transform.Rotation` conversions on 1000
random vectors are bit-identical across the two environments, which rules out that part only. The other
3 tests in the file pass in both environments; the 4 failing tests are the 4 that compare against
GOLDEN. This is one case, so other cases could still differ differently.

**Change.** `icenine_py/.python-version` (3.9); `uv sync --extra dev` is the documented setup in
`icenine_py/README.md`. `GOLDEN_ENV` next to the golden data in `tests/test_findoptimal_refactor.py`
records python 3.9 (major.minor), numpy 2.0.2, torch 2.8.0, scipy 1.13.1; the 4 golden tests take a
`golden_env` fixture that fails (not skips) with "golden values recorded with ...; running ... Run
`uv sync --extra dev`" on a mismatch. Tolerances unchanged.

**Verification.** In a temporary clone outside the repo, `env -u VIRTUAL_ENV uv sync --extra dev`
produced Python 3.9.6, numpy 2.0.2, torch 2.8.0, scipy 1.13.1 and left `uv.lock` unchanged; the 7 tests
in `test_findoptimal_refactor.py` pass. In the drifted 3.12 environment the 4 golden tests fail with
the version message. Main tree full suite: 663 passed, 34 skipped, 1 failed
(`test_dev_tooling.py::test_githook_and_installer_exist_and_are_not_activated`, which asserts
`core.hooksPath` is unset; it is now set to `.githooks` in this checkout, so the failure is the hooks
having been activated, not this change; fixed in d27575a). A later module-name clash between two
`test_nn_hybrid` summary imports in the full suite was fixed in 402fbd0 (summary loaded by path).

**Not done here.** The CLAUDE.md setup lines (`uv pip install -e ".[dev]"` to `uv sync --extra dev`) are
blocked by the pre-commit check on CLAUDE.md changes (override not allowed); the main session should make that
edit.

### T2 `stats-adoption`: shared statistics helpers in the study scripts (2026-10-06)

`scripts/findoptimal_robustness/common.py` and `scripts/profiling/prof_common.py` now compute their Wilson
interval with `scripts/common/stats.wilson` (the 3-tuple `(rate, lo, hi)` is rebuilt at the wrapper, so every
call site is unchanged). `scripts/nn_hybrid/summary.py` uses `stats.win_rate` (tie band `TIE_DEG = 0.002`);
`scripts/coarse_proxy/summary.py` uses `stats.paired_discordant` and `stats.mcnemar_exact`.
`prof_summary.py` delegates its Wilson to `prof_common` / the nn_hybrid summary; its only change is the gzip fix below.

Clamp fix accepted (main session decision): `nn_hybrid/summary.py` also uses `stats.wilson`, which clamps to
[0, 1]. The only change in the regenerated summaries is `wrong_hi` 1.0000000000000002 to 1.0 (one ulp, k = n):
4 values, one in each of `nn_hybrid/summary.json`, `summary_s1_replicate.json`, `summary_clean_s0_fallback.json`
and `profiling/summary.json` (the last reuses the nn_hybrid code); the `.txt` files are unchanged. These four
`.json` files were committed with the new value; all other summaries stay byte-identical.

Deterministic gzip: `prof_summary.write_gz` now writes with no mtime and no file name in the header, and
`compact()` iterates its key set in sorted order (it iterated a `set`, so the key order of `seeded_runs.json.gz` /
`u0_runs.json.gz` varied with the hash seed). The decompressed JSON is equal to the committed one as data
(key order differs once); the two `.gz` were re-committed once. Two further regenerations gave identical files
and a clean `git status`.

Kept local, on purpose:
- `nn_hybrid/summary.py::reorder`: it fills voxels absent from the source with NaN/0 (pilot on a subset);
  `stats.reorder` raises on a missing id.
- `prof_common.isolation_record`: not moved to `preflight.py`. The stored records (`isolation_*.json`, the
  `isolation` entry in every `v*.json`) have keys `time`, `battery`, `thermal`, `processes_over_5pct_cpu`,
  `cores`, `note`, which `prof_summary.isolation_lines` reads and `preflight.preflight()` does not produce, so the
  format would not stay readable.

Regenerated and confirmed byte-identical to the committed files (except the four `wrong_hi` values above):
- `benchmarks/nn_hybrid/summary.txt|json`, `summary_s1_replicate.txt|json` (`--model realistic_s1 --variant all
  --pipes H3`), `summary_clean_s0_fallback.txt|json` (`--model clean_s0 --variant clean --pipes H3`)
- `benchmarks/coarse_proxy/summary.txt|json`
- `benchmarks/profiling/summary.txt|json`
- `benchmarks/findoptimal_robustness/`: `fixes_summary.txt|json` (`summarize_fixes.py`), `e2_endtoend.txt|json`
  (`e2_summary.py`), `e0_summary.txt|json` and `e0_runs.npz` (`analyze_e0.py`)

Not regenerated: `nn_hybrid/pilot_summary.txt` (pilot subset; its command line was not recorded, not tried),
`coarse_proxy/offline_metrics.*`, `e0_dependence.*`, `e2_ablation.*`, `e2_results*.json`, `f1_remaining.txt`,
`s2_mechanism.txt` (do not use the replaced helpers, or need a model fit / training run; unit tests cover the
helpers). 

### T3 `lowq-batch`: batched F-lowQ Q5 feature pass (2026-10-06)

**Change.** `FeatureExtractor.features_batch(Rs, vertices, phase, with_cost=True)` in
`scripts/findoptimal_robustness/features.py` scores all candidates of one level in one pass: the geometry
(Bragg angles, reflection, ray-detector intersection, lit-pixel lookups) is vectorised over candidates and peaks,
the aggregates are integer counts per candidate (`np.bincount`, exact in float64), and the two cost columns reuse
the same geometry (the per-candidate path recomputes it twice) and call the existing C stage-D routine once per
candidate and radius, with the rows of each candidate in the cost function's order (the quality is a running
mean). Cost-function evaluation counters are advanced as the per-candidate path advances them. Chunks of 256
candidates; without the C stage-D extension it falls back to the per-candidate loop. `features` is unchanged and
stays the default and the reference; `scripts/coarse_proxy/endtoend.py --batched` opts in to the batched pass
(`proxy_features(..., batched=True)`). Nothing in `icenine/` changed.

**Equality.** Two details were needed for bit-identity: a batched `torch.matmul` of the rotations with the
reciprocal vectors differs from the 2-D matmul in the last bit (52346 of 271152 elements in one test), so
`g_lab` is built with one 2-D matmul per candidate, and it must keep the strided `(R @ g.T).T` layout (a
contiguous copy changes the last bit of the omegas). Intermediate states on the 5500 ad hoc candidates (the
Q_max-5 columns were identical throughout; only the unrestricted extractor with 112 reflections differed): batched
matmul, 5 rows differed (up to 0.037 in a hit rate, 0.003 in a cost); per-candidate matmul but contiguous layout,
2 rows (cost columns). With both, `features_batch` equals `features` exactly (`np.array_equal`,
all columns including both costs) on the E2 candidates of 6 random cases (5500 candidates, Q_max 5 and the
unrestricted extractor; ad hoc check) and on the 5000 candidates of the timing run (25 cases, Q_max 5:
5000 / 5000 rows identical; `timing_batch.json`). `tests/test_coarse_proxy.py::test_features_batch_equals_per_candidate`
covers Q_max 5 and None, batch sizes 1, 7 and the whole set, the evaluation counters, `with_cost=False` and the empty
batch on one voxel's E2 candidates. Not tested: other detector geometries or crystal phases than the study's.

**Timing (single-worker, all thread counts 1, `require_quiet()` preflight saved as
`benchmarks/coarse_proxy/timing_batch_preflight.json`).** Q_max-5 pass with both low-Q costs; 25 random cases x 200
candidates, 3 repetitions, methods interleaved in each repetition; per candidate = total / 200; median over cases.
One Q_max-8 local evaluation (the unit) was 0.529 ms in the same run. "proxy" adds the GBT prediction
(0.0064 ms per candidate, batch 200, from `timing.json`).

<!-- table:t3_lowq_batch_timing -->
| path | batch | ms | speedup | equiv | proxy |
|---|---|---|---|---|---|
| per-candidate (reference) | 1 | 1.555 | 1.0 | 2.94 | 2.95 |
| batched | 1 | 1.074 | 1.4 | 2.03 | 2.04 |
| batched | 50 | 0.151 | 10.3 | 0.28 | 0.30 |
| batched | 200 | 0.133 | 11.6 | 0.25 | 0.26 |
<!-- /table:t3_lowq_batch_timing -->

(Batch size 1 is the batched code called once per candidate: its 1.07 ms is the fixed per-call overhead. Geometry
only, without the two costs: 0.734 ms per candidate per-candidate, 0.124 ms batched at 200.)

**New proxy cost.** At batch 200 the F-lowQ Q5 pass is 0.25 evaluation equivalents (2.94 per candidate in this run, 2.93 in
Task 2), the proxy with its GBT prediction 0.26 (2.95), a factor 11.6; at batch 50, 0.28 and 0.30. A deployable E2 pass was estimated
at about 3.2 equivalents (Task 2 accounting, not re-timed here), so the proxy is now about 12x cheaper than that,
not 1.1x. Scoring the 240 / 384 candidates of one proxy run would take about 32 / 51 ms at the batch-200 rate
(computed from the per-candidate time, not measured end to end), against 0.37 / 0.59 s in Task 2. The levels of a
real run have fewer candidates than 200 per call; the realized gain end to end depends on the batch sizes of
`rank_key` (timed in T4: proxy about 0.05 s per run, single worker). The GBT time added at every batch size is the batch-200
per-candidate cost, which understates the cost at batch 1; levels 2-3 of a run have 1.3-7 candidates per call. The T4 wall timing
bounds this effect (proxy about 0.05 s per run), so the verdict holds.

Tests: `tests/test_coarse_proxy.py` 11 passed (2 new); full suite 666 passed, 34 skipped, 1 deselected, 0 failed.
Audit of this section: the numbers not found in the source JSONs are unit conversions (the JSON is in seconds, the
text in ms), the ad hoc equality check counts above (not saved), and figures quoted from Task 2 (3.2, 0.59).

### T4 `keep-eighth`: proxy rerank at keep 1/8 and 1/6, end to end (2026-10-06)

**Setup.** `scripts/coarse_proxy/endtoend.py run --batched` (set `lowq5+c8`, regression target, the proxy rerank is
symmetry-agnostic; no F1/F1b rows were run), 200 voxels x 2 variants, seed 0, 10 workers, keep_fraction 1/8 (tag
`p_k8`) and 1/6 (`p_k6`), paired with the stored baseline (E0) and the stored keep-1/4 proxy row (`p_i`).
`keep_eighth.py` (`summary`, `timing`, `timing-summary`) builds the tables; `endtoend.py` now also stores the
batch size of every `rank_key` call (`call_levels`, `call_sizes`). Nothing in `icenine/` changed.

**The batched pass changes nothing but the time.** On 6 runs (3 voxels x 2 variants, keep 1/8) the default path and
`--batched` gave identical `R_final`, `cost_final` and evaluation counts. A full batched re-run of keep 1/4 (`p_k4b`)
reproduced all 400 stored `p_i` runs exactly (`R_final`, global and local evaluations, scored candidates). The
keep-1/4 row below is the stored one; its call sizes come from `p_k4b`.

**Candidates per `rank_key` call.** Four calls per run, one per level (mean [min-max] over the 400 runs; `keep_eighth.txt`):
level 0 220.9 [95-562] candidates in all rows; keep 1/4: level 1 63.3 [26-165], level 2 20.7 [7-61], level 3 7.1 [1-24];
keep 1/6: 42.4 [19-111], 9.4 [3-29], 1.8 [1-9]; keep 1/8: 32.0 [15-90], 5.1 [1-18], 1.3 [1-4]. Only the first call is a large batch; the later
levels are small, and the batched cost per candidate is much higher there (T3: 2.04 eq at batch 1, 0.30 at 50, 0.26 at 200).
The proxy's own cost below is each call's batch size times the T3 cost interpolated (log batch size) at that
size, summed per run; it is an estimate from the T3 single-worker timing, not measured in these runs.

**End to end** (paired on the same 200 voxels per variant; `peq` = proxy evaluation equivalents per run at the batched
cost; `tot` = global + local evaluations of the reconstructor + `peq`; the reconstructor's own counts are in `glob` / `loc`).

<!-- table:t4_keep_eighth -->
| variant | row | wrong | rate | ci | med | glob | loc | scored | peq | tot |
|---|---|---|---|---|---|---|---|---|---|---|
| clean | E0 baseline | 67 | 0.335 | [0.273, 0.403] | 0.0301 | 44478 | 4711 | 0 | n/a | 49190 |
| clean | proxy keep 1/4 | 13 | 0.065 | [0.038, 0.108] | 0.0281 | 44494 | 5739 | 240 | 81 | 50313 |
| clean | proxy keep 1/6 | 18 | 0.090 | [0.058, 0.138] | 0.0359 | 44282 | 4473 | 212 | 72 | 48828 |
| clean | proxy keep 1/8 | 22 | 0.110 | [0.074, 0.161] | 0.0296 | 44196 | 4325 | 201 | 68 | 48590 |
| realistic | E0 baseline | 52 | 0.260 | [0.204, 0.325] | 0.0278 | 44778 | 7628 | 0 | n/a | 52406 |
| realistic | proxy keep 1/4 | 15 | 0.075 | [0.046, 0.120] | 0.0285 | 44805 | 9111 | 384 | 118 | 54034 |
| realistic | proxy keep 1/6 | 23 | 0.115 | [0.078, 0.167] | 0.0378 | 44465 | 6676 | 337 | 103 | 51244 |
| realistic | proxy keep 1/8 | 44 | 0.220 | [0.168, 0.282] | 0.0358 | 44323 | 6010 | 317 | 98 | 50430 |
<!-- /table:t4_keep_eighth -->

Exact McNemar (discordant counts are only-first-wrong / only-second-wrong):

| comparison | clean | realistic |
|---|---|---|
| keep 1/6 vs baseline | 5 / 54, p = 1.9e-11 | 10 / 39, p = 3.9e-05 |
| keep 1/8 vs baseline | 6 / 51, p = 5.7e-10 | 28 / 36, p = 0.38 |
| keep 1/6 vs keep 1/4 | 6 / 1, p = 0.125 | 12 / 4, p = 0.077 |
| keep 1/8 vs keep 1/4 | 10 / 1, p = 0.012 | 37 / 8, p = 1.5e-05 |

Pooled over both variants (400 runs): keep 1/4 wrong 28 (0.070, Wilson 0.049-0.099), keep 1/6 41 (0.102,
0.076-0.136), keep 1/8 66 (0.165, 0.132-0.205), baseline 119. Keep 1/6 vs keep 1/4: 18 / 5, p = 0.011; keep 1/8 vs keep 1/4: 47 / 9,
p = 2.6e-07. Both rows are better than the baseline (1/6: 15 / 93, p = 6.4e-15; 1/8: 34 / 87, p = 1.6e-06).

**"Useful"** (wrong <= the keep-1/4 proxy's Wilson upper bound, at fewer total evaluations than the baseline):
- Per variant: keep 1/6 is useful in both (clean 0.090 <= 0.108, 48828 < 49190; realistic 0.115 <= 0.120,
  51244 < 52406). Keep 1/8 is not useful in either: clean 0.110 vs the bound 0.108 (22/200, a miss by 0.002; total 48590 < 49190),
  realistic 0.220 vs 0.120.
- Pooled over the 400 runs: neither is useful. Keep 1/6 is close (0.102 vs the pooled bound 0.099; total 50036 vs
  baseline 50798, -1.5%); keep 1/8 is not (0.165; total 49510, -2.5%). The plan does not say per variant or pooled; both are shown.
- **Decision (main session).** The pooled verdict is the one that counts: it uses all 400 runs, and the per-variant pass
  for keep 1/6 is within 0.005-0.018 of the bound on one seed. Neither keep 1/6 nor keep 1/8 is adopted. Against keep 1/4,
  keep 1/6 is wrong in 18 more paired cases and right in 5 more (p = 0.011), and keep 1/8 is wrong in 47 more and right in
  9 more (p = 2.6e-7). Measured on a single worker, keep 1/6 saves only 2.5% of wall time and keep 1/8 only 3.3%.
  Tightening the keep fraction is therefore not a route to a run-time saving in the no-start case. Time is spent in the
  global discrete search, which these settings barely change (global evaluations, pooled, -0.6% at keep 1/6 and -0.8% at keep 1/8).
  The keep-1/4 proxy rerank stays the recommendation, for accuracy.
- The keep-1/4 row itself is not useful by this definition: it uses more evaluations than the baseline (+3.1% realistic at the
  batched proxy cost, +2.3% clean), consistent with the earlier T3 result.
- The reconstructor's own evaluation count falls with the keep fraction (global + local, clean 49190 baseline -> 50233 / 48756 / 48521;
  realistic 52406 -> 53916 / 51140 / 50333), but the number of wrong answers rises, and the median error of the right answers is larger for 1/6
  (0.036 clean, 0.038 realistic) and for realistic 1/8 (0.036) than for keep 1/4 (0.028) and the baseline (0.028-0.030) (deg; not tested further). Realistic keep 1/8 broke 28 baseline-right cases (vs 3 for keep 1/4).

**Single-worker timing** (first 20 voxels x 2 variants = 40 cases, seed 0, batched pass, `require_quiet()` preflight saved as
`benchmarks/coarse_proxy/keep_eighth_preflight.json`, all thread counts 1; arms rotated per case so none is always first;
the keep-1/4, 1/6, 1/8 and baseline arms were all run because 1/6 is useful per variant and 1/8 is close in the clean variant).
Every proxy run reproduced the stored 10-worker run exactly (40/40 for each arm); the baseline arm differs from the proxies in all 40.

<!-- table:t4_keep_eighth_timing -->
| arm | wall | speed | ev | dev | proxy | msev | same |
|---|---|---|---|---|---|---|---|
| base | 21.81 | 1.000 | 50942 | +0.0 | 0.00 | 0.428 | 40 |
| p_i | 22.53 | 0.968 | 52077 | +2.2 | 0.05 | 0.433 | 0 |
| p_k6 | 21.25 | 1.026 | 49844 | -2.2 | 0.05 | 0.426 | 0 |
| p_k8 | 21.09 | 1.034 | 49522 | -2.8 | 0.05 | 0.426 | 0 |
<!-- /table:t4_keep_eighth_timing -->

`wall` = mean `reconstruct_voxel` seconds per run (including the proxy), `speed` = baseline total time / arm total time, `ev` =
mean global + local evaluations per run, `dev` = change against the baseline in percent, `msev` = wall per evaluation (ms).
The saving in time follows the saving in evaluations: 0.426 ms per evaluation for keep 1/6 and 1/8 against 0.428 for the baseline (proxy 0.05 s per
run in all arms, about 0.2% of the run). Keep 1/6 saves 2.5% of the time and 2.2% of the evaluations, keep 1/8 3.3% and 2.8%; keep 1/4 costs 3.3% more time.
This is one 40-case subset with no interval on the time difference, and these 40 cases are a subset of the 400 (their evaluation
counts differ from the 200-voxel means above).

**Reading.** At keep 1/6 the proxy rerank gives a wrong rate that stays inside the keep-1/4 interval (per variant) for a 2-3% saving of time;
the pooled wrong rate is slightly above the keep-1/4 upper bound, and 1/8 is clearly worse. The time saving is small
(about 3%), consistent with global evaluations barely changing (-0.6% / -0.8%); the stage split of the wall time was not measured. Not tested: other keep fractions, other seeds
(the stored keep-1/4 row has seeds 1-2; these new rows do not), whether the extra wrong answers at 1/8 are candidates lost at the pruning step
(the harvest step was not run on the new rows).

### T5 `finisher-diagnosis`: why the finisher stops above the truth's cost (2026-10-07)

**Setup.** `scripts/finisher_diagnosis/diagnose.py run --workers 10` (wall 2852 s, 54 tasks; contended 10-worker run, no
timing claims are made), `diag_summary.py` builds the tables and `summary.json` (`benchmarks/finisher_diagnosis/`). Cases are the
realistic ("all") single-voxel cases of the Task 1 pipelines H3 (net x3 start; 200 cases) and H0 (perturbed nominal start;
100 cases), radii 0.05-3 deg (9 radii), stratified by radius (H3: 22-23 per radius over 4 random voxels; H0: 11-12 over 2
voxels), seeded selection. The finisher is re-run with the Task 1 seed through the unmodified `refine_from_candidates`, with
`LoggedMC` (a `MCOptimizer` subclass in the script whose `optimize` / `variance_minimizing_optimize` are copies that also
record the stopping rule; the same random stream). All 300 re-runs reproduce the Task 1 `R_final` exactly (max abs difference
0.0). Nothing in `icenine/` changed. The config's `max_mc_steps` is 200 (not 3500), `successive_restarts` 2; the default
box is 0.329 deg. "Gap" = cost of the result minus cost at the truth; "above" = gap > 1e-9. Intervals are Wilson 95%.

**Result vs truth.** The result's cost is above the truth's in 191/200 H3 and 100/100 H0 cases (H3 96%, CI 92-98%; H0 CI 96-100%).
In the 9 H3 cases at or below (5 strictly below) there is no gap.

<!-- table:t5_overview -->
| set | n | result above truth | gap q25/50/75 | angle res-truth q25/50/75 (deg) |
|---|---|---|---|---|
| H3 | 200 | 191/200 (96%, 92-98) | 0.010/0.028/0.050 | 0.013/0.023/0.043 |
| H0 | 100 | 100/100 (100%, 96-100) | 0.033/0.103/0.497 | 0.032/0.089/0.409 |
<!-- /table:t5_overview -->

<!-- table:t5_by_radius -->
| set | radius (deg) | n | result above truth | gap q50 | angle res-truth q50 (deg) | mc_long reaches truth cost | vm_long reaches truth cost |
|---|---|---|---|---|---|---|---|
| H3 | 0.05 | 23 | 21/23 | 0.036 | 0.028 | 5/23 | 2/23 |
| H3 | 0.1 | 23 | 23/23 | 0.023 | 0.025 | 1/23 | 1/23 |
| H3 | 0.25 | 22 | 22/22 | 0.035 | 0.032 | 1/22 | 0/22 |
| H3 | 0.5 | 22 | 20/22 | 0.025 | 0.020 | 3/22 | 2/22 |
| H3 | 0.75 | 22 | 19/22 | 0.026 | 0.017 | 5/22 | 3/22 |
| H3 | 1 | 22 | 21/22 | 0.026 | 0.018 | 3/22 | 1/22 |
| H3 | 1.5 | 22 | 22/22 | 0.029 | 0.025 | 3/22 | 1/22 |
| H3 | 2 | 22 | 21/22 | 0.025 | 0.022 | 2/22 | 1/22 |
| H3 | 3 | 22 | 22/22 | 0.031 | 0.027 | 1/22 | 0/22 |
| H0 | 0.05 | 12 | 12/12 | 0.023 | 0.016 | 1/12 | 0/12 |
| H0 | 0.1 | 11 | 11/11 | 0.034 | 0.014 | 0/11 | 0/11 |
| H0 | 0.25 | 11 | 11/11 | 0.043 | 0.038 | 1/11 | 0/11 |
| H0 | 0.5 | 11 | 11/11 | 0.058 | 0.057 | 0/11 | 0/11 |
| H0 | 0.75 | 11 | 11/11 | 0.126 | 0.172 | 0/11 | 0/11 |
| H0 | 1 | 11 | 11/11 | 0.192 | 0.212 | 4/11 | 1/11 |
| H0 | 1.5 | 11 | 11/11 | 0.326 | 0.232 | 0/11 | 0/11 |
| H0 | 2 | 11 | 11/11 | 0.709 | 1.750 | 0/11 | 0/11 |
| H0 | 3 | 11 | 11/11 | 0.810 | 3.353 | 0/11 | 0/11 |
<!-- /table:t5_by_radius -->

H3: the gap (median 0.028) and the angle from result to truth (median 0.023 deg) are about the same at all radii, i.e. the
finisher ends about 0.02 deg from the truth whatever the net's start error was. H0: gap and angle grow with the radius
(at 2-3 deg the result is a median 1.75 / 3.35 deg from the truth); that is a different failure (not reaching the basin) and
the continuations below do not help there.

**(a) Cost along the straight path result -> truth** (41 uniform points plus 7 near the result at t = 0.001-0.1; 'barrier' =
highest interior cost minus the result's cost).

<!-- table:t5_geodesic -->
| set | path cost rises above result | barrier q25/50/75 | a point below result | a point <= truth cost | monotone |
|---|---|---|---|---|---|
| H3 | 73/200 (36%, 30-43) | 0.000/0.000/0.001 | 191/200 (96%, 92-98) | 117/200 (58%, 52-65) | 25/200 (12%, 9-18) |
| H0 | 38/100 (38%, 29-48) | 0.000/0.000/0.000 | 100/100 (100%, 96-100) | 36/100 (36%, 27-46) | 15/100 (15%, 9-23) |
<!-- /table:t5_geodesic -->

The path cost rises above the result's in 36% (H3) and 38% (H0) of the cases, but the rise exceeds the gap itself in only
11/200 (H3) and 1/100 (H0); the median barrier is 0. In 191/200 (H3) and 100/100 (H0) cases some interior point of the path
has a lower cost than the result, and in 117/200 (H3) / 36/100 (H0) some interior point has cost at or below the truth's. The
path is monotone non-increasing in 25/200 (H3) and 15/100 (H0). So no barrier of the size of the gap lies on the straight
line; it is a rough cost, not a smooth ramp.

**(c) Granularity** (12 random axes per angle, at the truth and at the result; 'same' = cost equal within 1e-9, 'lower' = below the
cost at the centre).

<!-- table:t5_granularity -->
| set | angle (deg) | truth: same cost | truth: lower | truth: median |d| | result: same cost | result: lower | result: median |d| |
|---|---|---|---|---|---|---|---|
| H3 | 0.001 | 0.18 | 0.07 | 0.0018 | 0.17 | 0.28 | 0.0015 |
| H3 | 0.003 | 0.01 | 0.01 | 0.0070 | 0.02 | 0.27 | 0.0049 |
| H3 | 0.01 | 0.00 | 0.00 | 0.0249 | 0.00 | 0.18 | 0.0171 |
| H3 | 0.02 | 0.00 | 0.00 | 0.0529 | 0.00 | 0.09 | 0.0396 |
| H3 | 0.05 | 0.00 | 0.00 | 0.1422 | 0.00 | 0.02 | 0.1193 |
| H3 | 0.1 | 0.00 | 0.00 | 0.2903 | 0.00 | 0.00 | 0.2561 |
| H3 | 0.2 | 0.00 | 0.00 | 0.4807 | 0.00 | 0.00 | 0.4433 |
| H0 | 0.001 | 0.19 | 0.06 | 0.0018 | 0.33 | 0.16 | 0.0008 |
| H0 | 0.003 | 0.01 | 0.02 | 0.0072 | 0.15 | 0.19 | 0.0027 |
| H0 | 0.01 | 0.00 | 0.00 | 0.0263 | 0.09 | 0.18 | 0.0109 |
| H0 | 0.02 | 0.00 | 0.00 | 0.0548 | 0.06 | 0.14 | 0.0222 |
| H0 | 0.05 | 0.00 | 0.00 | 0.1469 | 0.03 | 0.09 | 0.0738 |
| H0 | 0.1 | 0.00 | 0.00 | 0.2921 | 0.02 | 0.06 | 0.1655 |
| H0 | 0.2 | 0.00 | 0.00 | 0.4808 | 0.02 | 0.03 | 0.3157 |
<!-- /table:t5_granularity -->

At the truth no rotation of 0.01 deg or more lowered the cost, and at 0.001 / 0.003 deg 7% / 1% of the draws did (H3). The
median change of the cost for a 0.01 deg rotation at the truth (0.025 H3, 0.026 H0) is close to the median H3 gap (0.028): a
gap of this size corresponds to about 0.01 deg of rotation at the truth. The cost steps are fine against the gap: the smallest nonzero difference
between distinct cost values seen among each case's 84 truth-neighbour evaluations has a median of 2.4e-5 (H3) and 2.4e-5 (H0)
(an upper bound on the quantum, from a small sample). At the result, rotations of 0.001-0.02 deg lowered the cost in 28% to 9%
(H3) of the draws (H0: 14-19% at 0.001-0.02 deg), so the result is not a local minimum on that scale in 175/200 (H3) and 80/100 (H0) cases (sampled: 60 neighbours at 0.001-0.05 deg, `res_local_min_0.05deg`). The cost at the truth changes by
a median 0.0018 for a 0.001 deg rotation, i.e. already rough at that scale. Median counts at the truth / result (pixel overlap,
pixels on detector, peaks overlapping, peaks on detector, peaks counted) are in `summary.json` (`counts_*`).

**(d) Stopping rules in the default finisher run.**

<!-- table:t5_stopping -->
| set | FO stop codes (budget/restarts/cost) | FO steps q50 | FO accepts q50 | FO last accept q50 | FO final step q50 (deg) | VM steps q50 | VM lowers cost | VM final radius q50 (deg) |
|---|---|---|---|---|---|---|---|---|
| H3 | 122/78/0 | 200 | 8 | 24 | 0.0008 | 2245 | 107/200 (54%, 47-60) | 0.329 |
| H0 | 97/3/0 | 200 | 6 | 28 | 0.0021 | 375 | 80/100 (80%, 71-87) | 0.329 |
<!-- /table:t5_stopping -->

FindOptimal MC (200 steps, step halved at every accepted improvement): the stop code counts are step budget / restarts exhausted /
cost converged. In H3 the restarts-exhausted stop fired in 78/200 runs, and these are exactly the 78 runs with no accepted step
(0 accepts in 78/200 H3 and 3/100 H0 runs); otherwise the run ended on the 200-step budget with a median of 8 (H3) / 6 (H0)
accepts, the last accept at step 24 (H3) / 28 (H0) (medians), and a final step size of 0.0008 (H3) / 0.0021 (H0) deg (median). The
hit-ratio "converged" flag never fired (0/200, 0/100) and the cost-convergence stop (`max_convergence_cost`) never fired. VarianceMinimizing
has `max_convergence_cost` 0, so it can only stop on its step budget, which is extended by one subregion (10 steps) for every
subregion whose cost variance is above 0.02^2 = 0.0004; it ended after a median 2245 (H3) / 375 (H0) steps (budget 200 plus the
extensions), with the subregion radius back at the box (0.329 deg) in most runs, and it lowered the cost of the FindOptimal
result in 107/200 (54%, 47-60%) H3 and 80/100 (80%, 71-87%) H0 runs. Which candidate the reconstructor returned as final was not logged,
so "VM supplied the final result" is not claimed.

**(b) Continuations from the result.** Components called directly from the script with changed parameters: `mc_long` (FindOptimal MC,
same box and step, 25x the steps = 5000, 10 restarts), `mc_smallstep` (step / 10, 10000 steps, 5 restarts), `mc_smallbox` (box and step / 4,
10000 steps, 5 restarts), `vm_long` (VarianceMinimizing, same box, budget 5000 steps), `vm_smallbox` (box / 4, budget 5000), `rerun_default`
(the whole default finisher again from the result, another seed) and `from_truth` (the default finisher started at the truth). The
VarianceMinimizing runs were stopped at 50000 steps if the variance rule had not ended them (a safety cap that is not in the reconstructor;
`vm_capped` in `summary.json`: reached in 139/200 `vm_long` and 198/200 `vm_smallbox` H3 runs, 65/100 and 83/100 H0 runs, so those two rows
are mostly capped runs, not rule-terminated ones; the default re-runs were never capped). 'Gap closed' is the median over cases with
gap > 0 of (result cost - end cost) / gap; 'toward truth' is the decrease of the angle to the truth.

<!-- table:t5_continuations -->
| set | continuation | improves | reaches truth cost | gap closed (median) | toward truth q25/50/75 (deg) | moved q50 (deg) |
|---|---|---|---|---|---|---|
| H3 | mc_long | 53/200 (26%, 21-33) | 24/200 (12%, 8-17) | 0.00 | -0.000/0.000/0.001 | 0.000 |
| H3 | mc_smallstep | 194/200 (97%, 94-99) | 52/200 (26%, 20-32) | 0.61 | 0.002/0.006/0.011 | 0.011 |
| H3 | mc_smallbox | 137/200 (68%, 62-75) | 29/200 (14%, 10-20) | 0.69 | 0.000/0.006/0.023 | 0.016 |
| H3 | vm_long | 113/200 (56%, 50-63) | 11/200 (6%, 3-10) | 0.24 | -0.000/0.000/0.028 | 0.015 |
| H3 | vm_smallbox | 186/200 (93%, 89-96) | 42/200 (21%, 16-27) | 0.93 | 0.005/0.016/0.035 | 0.022 |
| H3 | rerun_default | 67/200 (34%, 27-40) | 13/200 (6%, 4-11) | 0.00 | -0.000/0.000/0.010 | 0.000 |
| H3 | from_truth | 0/200 (0%, 0-2) | 200/200 (100%, 98-100) | nan | -0.000/-0.000/-0.000 | 0.000 |
| H0 | mc_long | 65/100 (65%, 55-74) | 6/100 (6%, 3-12) | 0.05 | -0.000/0.005/0.092 | 0.058 |
| H0 | mc_smallstep | 88/100 (88%, 80-93) | 12/100 (12%, 7-20) | 0.12 | 0.001/0.007/0.011 | 0.013 |
| H0 | mc_smallbox | 76/100 (76%, 67-83) | 10/100 (10%, 6-17) | 0.15 | 0.000/0.018/0.034 | 0.030 |
| H0 | vm_long | 78/100 (78%, 69-85) | 1/100 (1%, 0-5) | 0.76 | 0.000/0.045/0.227 | 0.073 |
| H0 | vm_smallbox | 90/100 (90%, 83-94) | 14/100 (14%, 9-22) | 0.96 | 0.008/0.043/0.205 | 0.060 |
| H0 | rerun_default | 60/100 (60%, 50-69) | 3/100 (3%, 1-8) | 0.16 | -0.000/0.016/0.136 | 0.044 |
| H0 | from_truth | 0/100 (0%, 0-4) | 100/100 (100%, 96-100) | nan | -0.000/-0.000/-0.000 | 0.000 |
<!-- /table:t5_continuations -->

The default finisher started at the truth never lowered the cost (0/200, 0/100). Reaching exactly the truth's cost is rare for
any single continuation (H3 at most 52/200 = 26%, CI 20-32%, for `mc_smallstep`; H0 at most 14/100), but the cost drops a lot: `vm_smallbox` closes
a median 93% (H3) / 96% (H0) of the gap (but it sits at the 50,000-step cap in 198/200 and 83/100 runs, about 19x / 81x the default evaluations) and `mc_smallstep` improves the cost in 97% (H3) of the runs, and
ends a median 0.006 deg (`vm_smallbox`, H3) from the truth against 0.023 deg at the result. The best of the five continuations per case
closes a median 99% of the gap and reaches the truth's cost in 96/200 (48%, 41-55%) (H3) and 35/100 (35%, 26-45%) (H0) cases (`best_of_five` in `summary.json`). The same-step `mc_long` improves the cost
in only 53/200 H3 runs with a median movement of 0.000 deg, although its budget is 5000 steps (25x the default 200). Logged
(`mc_long_log` in `summary.json`): it stopped on restarts exhausted in 154/200 H3 runs (H0: 44/100), with a median of 0 accepts
(H0: 4.5) and 21 restarts (H0: 5.5); it ran a median of only 821 of its 5000 steps in H3 (H0: 5000), i.e. most H3 runs used about 4x the
default budget, and its final step was still the initial 0.132 deg (H3 median; H0 0.006 deg). Why it finds no first accept was not tested.
`mc_smallstep` (1/10 of the step) improves 194/200 H3 runs. In H0 the
continuations move the orientation toward the truth by a median of at most 0.045 deg and mostly do not reach the truth's cost at 0.75 deg and above
(by-radius table: `vm_long` reaches it in 1/11 at 1 deg and in none elsewhere above 0.5 deg; `mc_long` in 4/11 at 1 deg and in none elsewhere above 0.5 deg).

**Reading (observations only).** In 96% of the H3 cases the finisher's result is about 0.02 deg from the truth with a cost 0.03 above
it; the cost along the straight line to the truth has no barrier as big as that gap, a 0.01 deg rotation at the truth changes the cost by
about as much as the gap, the finisher's FindOptimal step ends 0.001-0.002 deg and the result is not a local minimum on the 0.001-0.05 deg
scale in 175/200 (H3) and 80/100 (H0) cases (sampled, 60 neighbours), and continuing with a smaller step or a smaller VarianceMinimizing box lowers the cost to within a small remainder of the truth's in many
cases. This does not show why the default parameters stop where they do. Not tested: whether the minimum of the realistic-data cost is at
the truth (the truth's cost is itself noisy: the cost at the truth is above the cost of some nearby orientation in 30% of H3 cases at 0.05 deg,
see `truth_best_within_0.05deg` in `summary.json`), other box and step ratios, the effect of the 200-step MC budget on the result, and symmetry
equivalents of the truth (the work is symmetry-agnostic). Angles are to the stored truth; the largest is 1.9 deg (H3) and 32.2 deg (H0). Below 45 deg the unreduced cubic angle equals the
reduced one (the smallest cubic symmetry rotation is 90 deg) and the cost is symmetry-invariant, so no reported number changes; for lower symmetries the bound is smaller (e.g. 30 deg for 6-fold).

**Effect estimate for a possible fix (diagnostics only, not implemented).** Appending a VarianceMinimizing pass with a quarter-size box after
the default finisher (`vm_smallbox`) reduces the H3 gap by a median 93% and moves the result a median 0.016 deg toward the truth, reaching the truth's cost in 42/200 (21%, 16-27%) cases
(14/100 for H0). It reached the 50,000-step cap in 198/200 H3 runs (83/100 H0; median 50,000 steps) against a median default finisher of
2,629 (H3) / 616 (H0) evaluations, i.e. about 19x (H3) / 81x (H0) the default evaluation count, in essentially every run
(`vm_smallbox_cost` in `summary.json`); a cheaper variant (a fixed 5000-step budget) has not been estimated here. An MC pass with step / 10 (`mc_smallstep`, 10000 evaluations) closes a median 61% of the H3 gap and
reaches the truth's cost in 26%. Neither was tested for its effect on the reconstruction's final success rate, only on the cost and the angle to the truth in these
cases.

### Completion summary (Follow-ups, 2026-10-07)

**Status: complete.** T1-T5 merged into `feature/followups`; the branch review's doc and test fixes are applied (no new experiments).

**Decisions.**
- The Wilson clamp is accepted (4 ulp values).
- Keep 1/6 and keep 1/8 are not adopted, on the pooled verdict (41/400 and 66/400 wrong against the bound 0.099; McNemar against keep 1/4: p = 0.011 and 2.6e-7).
- The keep-1/4 proxy rerank stays, for accuracy.
- The golden tests fail rather than skip in other environments (they skip only when the data is absent).
- `nn_hybrid` keeps its local reorder.
- `isolation_record` stays out of preflight.

**Headline per task.**
- T1: the drift is at float level (6e-9), not a branch flip.
- T2: regeneration is byte-identical and the gzip is deterministic.
- T3: exact equality of the features; 11.6x at batch 200 (0.26 evaluation equivalents).
- T4: about 2.5-3.3% wall time saved, at an accuracy cost.
- T5: the result is about 0.023 deg from the truth with gap 0.028; no barrier on the straight path is as large as the gap; the default MC step is too coarse at the end; `vm_smallbox` closes 93% of the gap at about 19x the default evaluations (at the 50,000-step cap in 198/200 runs).

**Recommendations.**
- Test a cheaper fixed-budget small-box VM or step/10 finisher pass for its effect on reconstruction success.
- Use `--batched` in production proxy runs.

**Open items.**
- The CLAUDE.md `uv sync` edit (owner).
- Seeds 1-2 for the keep rows.
- Harvest on keep 1/8.
- Whether the realistic-cost minimum sits at the truth (the truth is not the local best in 59/200).
- The H0 wrong-basin failures at 2-3 deg.
- The stage split of wall time.
- `pilot_summary.txt` was not regenerated.
- Full-sample renders and the three future ideas (owner direction).

Final full suite on `feature/followups`: 673 passed, 34 skipped, 0 failed (707 collected).

## Bug fix: MCOptimizer C++ parity (2026-10-09)

Branch `feature/mc-cpp-parity-fix` (off develop at 9863007): the two places where `MCOptimizer` did not follow C++, brought over from `feature/finisher-mc-study` (commits 2c579c9 and 19933d5) without any study code. Library change only in `icenine/orientation_search.py`; nothing else in `icenine/` changed.

| # | Function | C++ (`OrientationSearch.cpp`) | Python before |
|---|---|---|---|
| 1 | `variance_minimizing_optimize` restart | offsets uniform in +-SubregionRadius (the radius just doubled, capped at the box), applied to the INITIAL orientation (lines 248-255) | offsets uniform in +-box/2 about the current global best |
| 2 | `optimize` (quick MC and the FindOptimal MC) | `RandomRestartZeroTemp` (lines 297-363, with `ZeroTemperatureOptimization` 100-135): fixed blocks of `nMinErgodicSteps = int(2 (box/step)^3)` computed once; the step halves only after a block ends strictly below the global best; a failed block restarts at `delta * initial`, `delta` from U(-r, r)^3 with `r = tan(box)/sqrt(48)`, step reset; stop after more than `SuccessiveRestarts` failures in a row | one step at a time, halving at every improvement of the global best with the threshold recomputed, restart about the global best with +-box/2 |

Implementation: `MCOptimizer.optimize` is the block loop, `_mc_block` its inner block, `last_run` the bookkeeping of the last call (stop reason, blocks, accepts, restarts); the trajectory records one event per block (`mc_accept`, `mc_restart`). `_on_block_end` is an empty per-block hook (kept so the study scripts' subclasses work unchanged on the feature branch). The CMA optimizer, the `rng` property and the extra `SearchParameters` keys of the feature branch are not on this branch.

**Verification** (done on the feature branch, full record there): debug prints in a scratch build of the C++ `RandomRestartZeroTemp` and `AdaptiveSamplingZeroTemp` traced block lengths, step halving and reset, and restart radii with zero block-rule violations (all 93 Python FindOptimal calls checked as well). On the 6 interior seed voxels of the full-sample seed-cost diagnosis (one seed each, n = 6, voxel 15901 timed contended) the variance stage takes, after fix 1 and measured before the item-2 port, a mean of 1,609 evaluations in Python against 1,295 steps in C++ (C++ counts steps, which are about 10% fewer than evaluations), down from 707,732 evaluations before fix 1.

**Changed meanings.** The trajectory has one record per block (`step` is the cumulative step count at the block end; a restart record's `angular_step_deg` is measured from the previous best), where it had one record per accepted step; the stored `hp_sweep_trajectory_*.csv` and the description of the old records earlier in this file use the old meaning. In `last_run` and the finisher-diagnosis `LOG_KEYS`: `since_improve` is the number of successive failed blocks (was trial steps since the last improvement), `min_ergodic` is the fixed initial block length (was the last recomputed threshold), `steps_run` counts trial steps only. The stored `benchmarks/finisher_diagnosis/summary.json` uses the old meanings.

**Stale results on develop.** The MC results recorded on develop with the old optimizer (HP sweep, finisher diagnosis, FindOptimal robustness) predate this fix and are not re-run here.

**Effect on existing behaviour.** This changes the default classic reconstruction path (the quick MC and FindOptimal MC draw different random numbers and run a different block structure), so results are not bit-identical to before and no earlier random-draw pin is valid. Tests: new `tests/test_variance_stage_parity.py` and `tests/test_mc_cpp_parity.py` (scripted semantics tests that fail on the old code); `GOLDEN` in `tests/test_findoptimal_refactor.py` (4 tests) re-recorded (old values in the comment); `test_logged_mc_matches_mcoptimizer` (3 cases) passes again after `scripts/finisher_diagnosis/diagnose.py` `LoggedMC` was changed to log the library `optimize` and to use the new variance restart; `test_findoptimal_alone_reproduces_stored_experiment_b` (skipped when the sweep raws are absent) patches in the frozen pre-fix methods of `tests/legacy_variance_stage.py`, since the stored experiment predates both fixes. No `cpp_outputs`-based check changed. For the reruns of earlier MC results with the corrected optimizer, see `feature/finisher-mc-study`, MIGRATION_HISTORY subsections "C++ vs Python variance stage" and "C++-faithful MC: reruns".
