# IceNine Python/PyTorch Implementation

Python port of IceNine for synchrotron X-ray diffraction — forward simulation and reconstruction of crystal grain orientations. Forward simulation validated pixel-exact against the C++ implementation.

## Installation

```bash
cd icenine_py
uv pip install -e ".[dev]"
```

Dependencies: numpy, torch, pymatgen, scipy (see `pyproject.toml`).

## Quick Start: Forward Simulation

```bash
cd Examples/Example2.ThreeVoxels
uv run python run_python_simulation.py
```

This simulates X-ray diffraction from a 3-voxel copper sample (64.35 keV beam, 180 omega steps × 2 detectors) and writes 360 detector images to `ScatteringData_Python/`.

### Programmatic Usage

```python
from icenine.config_file import ConfigFile
from icenine.forward_simulation import ForwardSimulation

config = ConfigFile.from_file("path/to/experiment.config")
simulator = ForwardSimulation(config)

# Serial path (reference implementation, non-differentiable)
images = simulator.simulate_detector_images(output_dir="output/")

# Batched path (differentiable through stages 1-5, torch autograd preserved)
images = simulator.simulate_detector_images(output_dir="output/", batched=True)

# Batched with memory chunking (for large samples)
images = simulator.simulate_detector_images(
    output_dir="output/", batched=True, batch_size=5000
)
# images[omega_index][detector_index] is an ImageData object
```

### What ForwardSimulation Does

For each voxel in the sample:
1. Loads crystal structure and generates reciprocal lattice vectors (filtered by MaxQ and MinAmplitudeFraction)
2. For each reciprocal vector, solves the Bragg condition for omega angles
3. Maps each omega to an experimental wedge via `SimulationRange`
4. Rotates the sample to the omega angle
5. Projects the voxel triangle onto all detectors (Sutherland-Hodgman clipping + scanline rasterization)
6. Accumulates intensity with Lorentz-polarization correction: `I = I_form / (|sin(η)| × sin(2θ))`

## Package Modules

| Module | Description |
|--------|-------------|
| `forward_simulation.py` | Main simulation loop — generates detector images from a sample |
| `simulation.py` | Core engine — Bragg condition solving, vertex projection, voxel rasterization |
| `experiment_setup.py` | Reads config, initializes detectors/sample/omega ranges |
| `config_file.py` | Parser for `.config` files (80+ parameters, auto degree→radian conversion) |
| `detector.py` | Detector geometry, coordinate transforms, ray-plane intersection |
| `image_data.py` | Detector image storage, scanline triangle rasterizer |
| `sample.py` | Sample with orientation, translation, crystal structures |
| `mic_file.py` | Read/write `.mic` voxel grid files (Bunge Euler angles) |
| `crystal_structure.py` | Unit cell, reciprocal lattice, reflection vector generation |
| `diffraction_core.py` | Scattering omega calculation, reflected rays (PyTorch batched) |
| `peak_filters.py` | Eta-angle acceptance filter with Lorentz-polarization correction (`batch_eta_filter()` for vectorized ops) |
| `simulation_range.py` | Omega range system for discontinuous data collection wedges |
| `geometry.py` | Euler conversions, Plane/Ray classes |
| `symmetry.py` | Crystal symmetry operations (pymatgen wrapper) |
| `constants.py` | Physical constants |
| `file_io.py` | Detector file, crystal structure file, and omega file I/O |
| `reconstructor.py` | Reconstruction orchestrators — serial, adaptive, and BFS spatial propagation |
| `orientation_search.py` | Discrete grid search + zero-temperature MC optimization |
| `cost_functions.py` | Overlap computation between simulated projections and experimental data (batched Stages A-C + sequential Stage D) |
| `_rasterize.c` | CPython C extension for fast triangle rasterization and pixel overlap (Sutherland-Hodgman + Bresenham) |
| `experimental_data.py` | Load experimental detector images for reconstruction |
| `sampling.py` | SO(3) uniform sampling via Sukharev grids (Yershova & LaValle) |

## Quick Start: Reconstruction

Reconstruction recovers crystal orientations from experimental (or synthetic) detector images. The workflow is: load a config → load experimental data → search orientation space → output a `.mic` file with fitted orientations.

### Using Synthetic Data (Forward Sim → Reconstruct)

The simplest way to test reconstruction is with synthetic data from a forward simulation:

```python
import os
from icenine.config_file import ConfigFile
from icenine.experimental_data import ExperimentalData
from icenine.reconstructor import setup_reconstruction, SerialReconstruction

# Work from the example directory (config uses relative paths)
os.chdir("Examples/Example2.ThreeVoxels")

config = ConfigFile.from_file("ConfigFiles/Example2.Simulation.config")

# Point to forward simulation output as "experimental" data
config.out_file_basename = "3Grains.sim"
exp_data = ExperimentalData.from_image_directory(
    directory="ScatteringData_Python",
    basename="3Grains.sim",
    ext="d",
    serial_length=5,
    n_omega=180,
    n_detectors=2,
    num_rows=2048,
    num_cols=2048,
)

# Set up reconstruction (loads FZ orientations, detectors, sample, etc.)
setup = setup_reconstruction(config, exp_data=exp_data)

# Reconstruct all voxels and save result
recon = SerialReconstruction(setup)
results = recon.reconstruct_sample(
    output_mic="reconstructed.mic",
    max_voxels=3,  # limit for testing; remove for full sample
)

# Inspect results
for i, r in enumerate(results):
    print(f"Voxel {i}: cost={r.cost:.4f}, hit_ratio={r.hit_ratio:.3f}, "
          f"convergence={r.convergence_code}")
```

### Using Real Experimental Data

```python
config = ConfigFile.from_file("path/to/experiment.config")

# ExperimentalData loads from the InfileBasename/InfileExtension in config
setup = setup_reconstruction(config)

recon = SerialReconstruction(setup)
results = recon.reconstruct_sample(output_mic="output.mic")
```

### What Reconstruction Does

For each voxel in the sample grid:
1. **Discrete search**: Evaluates all FZ orientations × local grid perturbations
2. **Quick MC**: Runs short Monte Carlo optimization (20 steps) on top candidates
3. **Filter**: Keeps the best N candidates by cost
4. **Full MC**: Runs full MC optimization with restarts and convergence checking
5. **Adaptive deepening**: If not converged, refines the local grid and repeats

The search is multi-level adaptive — it starts with a coarse orientation grid and progressively refines around promising candidates until the cost function converges.

### BFS Reconstruction (Recommended for Large Samples)

For spatially coherent microstructures, BFS reconstruction is much faster than independent per-voxel search. It does a full adaptive search on a seed voxel, then propagates the orientation to neighbors via BFS, using cheap MC-only optimization:

```python
from icenine.reconstructor import setup_reconstruction, BFSReconstruction

setup = setup_reconstruction(config, exp_data=exp_data)
recon = BFSReconstruction(setup)
processed = recon.reconstruct_sample(output_mic="bfs_result.mic")
```

BFS algorithm:
1. **Seed selection**: Pick next unvisited voxel (randomized order)
2. **Full search**: `AdaptiveVoxelReconstructor.reconstruct_voxel()` on seed (~100-150s)
3. **Propagate**: Copy seed orientation to all unvisited neighbors
4. **BFS loop**: For each neighbor, run `local_optimization()` (MC-only, ~1s)
5. **Accept/reject**: If neighbor quality > 90% of best, mark FITTED and propagate; else mark REFIT
6. **Repeat**: Until all voxels visited

For the ThreeVoxels test case: 252s (BFS) vs 6638s (serial) — 26× faster. The speedup is even larger for big samples where most voxels are interior neighbors.

### Key Config Parameters for Reconstruction

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| `FundamentalZoneFilename` | SO(3) sampling grid file | `DataFiles/MyFZ.dat` |
| `LocalOrientationGridRadius` | Local search radius (degrees) | 5 |
| `MinLocalResolution` / `MaxLocalResolution` | Adaptive deepening levels | 0 / 3-5 |
| `MaxMCSteps` | Monte Carlo steps per candidate | 300-3500 |
| `SuccessiveRestarts` | MC random restarts | 2-3 |
| `MaxConvergenceCost` | Cost threshold to stop MC early | 0.0001-0.01 |
| `MaxAcceptedCost` | Cost threshold to accept result | 0.9 |
| `MaxDiscreteCandidates` | Top candidates from discrete search | 50-100 |

### Single-Voxel Reconstruction

For debugging or testing, you can reconstruct a single voxel directly:

```python
import torch
from icenine.reconstructor import setup_reconstruction, BasicVoxelReconstructor, _get_voxel_vertices

setup = setup_reconstruction(config, exp_data=exp_data)
reconstructor = BasicVoxelReconstructor(setup)

# Get voxel from sample
mic = setup.sample.get_mic()
voxel = mic.voxels[0]
vertices = _get_voxel_vertices(voxel)

result = reconstructor.reconstruct_voxel(
    voxel_vertices=vertices,
    phase_index=voxel.phase,
)

print(f"Cost: {result.cost:.4f}")
print(f"Hit ratio: {result.hit_ratio:.3f}")
print(f"Orientation:\n{result.orientation}")
```

## Config File Format

Forward simulation requires a `.config` file specifying:

```
BeamEnergy           64.351           # keV
BeamDirection        0  0  1          # unit vector
MaxQ                 16               # Å⁻¹, max scattering vector magnitude
EtaLimit             86               # degrees, max eta for peak acceptance
MinAmplitudeFraction 0.25             # filter reflections below this fraction of max intensity
SampleFilename       SimInput/three_voxels.mic
StructureFilename    DataFiles/copper.dat
DetectorFilename     ConfigFiles/StandardGeometry.2Det
OmegaFilename        DataFiles/omega_180_2L.dat
OutFileBasename      3Grains.sim
OutFileExtension     d
OutFileSerialLength  5
```

See `Examples/Example2.ThreeVoxels/ConfigFiles/Example2.Simulation.config` for a complete example.

## Testing

```bash
cd icenine_py
uv run pytest tests/ -v                                 # all tests (326 passed, 34 skipped)
uv run pytest tests/test_simulation.py                   # specific module
uv run pytest tests/test_reconstruction_integration.py   # reconstruction end-to-end (~40s)
uv run pytest --cov=icenine tests/                       # with coverage
```

### Forward Simulation Integration Test (C++ vs Python)

```bash
cd Examples/Example2.ThreeVoxels
uv run python run_python_simulation.py   # generate Python output
uv run python compare_outputs.py         # compare against C++ reference
```

Expected: 3200/3201 pixels match at identical locations, max relative intensity difference < 1e-5.

### Reconstruction Integration Tests

`test_reconstruction_integration.py` runs 5 end-to-end tests using the ThreeVoxels example:
- Loads forward simulation output as synthetic experimental data
- Evaluates cost function at ground truth orientation (verifies overlap)
- Compares ground truth quality against random orientations
- Runs MC optimization from a perturbed starting point
- Verifies convergence to within 10° of ground truth

## Validation Status

### Forward Simulation

Validated pixel-exact against C++ on two test cases:

**ThreeVoxels** (3 voxels, 2 detectors, 180 omega steps):
- 3200/3201 pixels match, max relative intensity difference 5.4e-6
- 4 pixel-location mismatches at omega bin boundaries (floating-point rounding)

**ManyGrains** (24,570 voxels, 2 detectors, 180 omega steps):
- 99.97% pixel match rate (7,422,506 / 7,424,450)
- Pixel count ratio 1.0000, intensity ratio 1.000000

**Serial vs Batched** paths produce identical output (0.015% bin-boundary mismatches due to float32 precision in Bragg solver).

### Reconstruction

#### End-to-End Comparison (C++ vs Python)

Identical reconstruction on ThreeVoxels (MaxQ=8, 180 omega × 2 detectors, 4886 FZ orientations, 4 resolution levels):

| Metric | C++ | Python |
|--------|-----|--------|
| Total wall time | 24.3s | 6637.5s |
| Data loading | ~1s | 1.3s |
| Reconstruction | ~23s | 6635.8s |
| Per-voxel average | ~8s | 2211.9s |
| Reconstruction slowdown | 1× | ~276× |

**Per-voxel results (BasicVoxelReconstructor):**

| Voxel | C++ Cost | Python Cost | Python Euler (reconstructed) | Ground Truth Euler | Python Misori |
|-------|----------|-------------|------------------------------|--------------------|---------------|
| 0 | 0.172 | 0.057 | (355.42, 5.19, 29.32) | (355.43, 5.19, 29.32) | 0.01° |
| 1 | 0.080 | 0.201 | (155.62, 45.18, 209.33) | (155.44, 45.18, 29.33) | ~0° (sym equiv) |
| 2 | 0.818 | 0.111 | (356.80, 3.70, 328.39) | (356.74, 3.70, 328.45) | 0.01° |

Note: This comparison used Python `BasicVoxelReconstructor` vs C++ `DiscreteRefinement` — **different algorithms**. See the identical-algorithm comparison below.

#### Identical Algorithm Comparison (AdaptiveVoxelReconstructor)

Same config and data, using the **identical algorithm**: C++ `DiscreteRefinement` vs Python `AdaptiveVoxelReconstructor`. Both sides instrumented with exact evaluation counters.

| Metric | C++ | Python | Ratio |
|--------|-----|--------|-------|
| Total time | 2.25s | 97.0s | 43× |
| Total evals | 150,251 | 213,690 | 1.42× |
| Avg us/eval | 15.0 | 453.9 | 30× |

Per-voxel: Both find the same orientations for voxels 0 and 1 (same local minima). For voxel 2, Python succeeds (0.13° misori) while C++ fails (cost=0.818) due to candidate count differences from floating-point divergence. See [MIGRATION_HISTORY.md](MIGRATION_HISTORY.md) for full per-voxel tables.

#### Unit-Level Validation

- Ground truth orientations produce high overlap (hit ratio > 0.5, quality > random)
- MC optimizer converges from 1.5° perturbation to within 10° of ground truth
- Integration tests run in ~40s (optimized from ~380s via bounding-box overlap computation)
- Cost function evaluate(): 503 us (~12x vs C++ 42 us, optimized via batch C extension + binary image cache)

### Cost Function Pipeline (`calculate_diffraction_overlap_batched`)

The batched cost function is organized into four stages:

| Stage | What | Mode |
|-------|------|------|
| A | Map peak omegas → wedge indices, filter invalid | Numpy vectorized |
| B | Batch rotation/reflection, vertex transform to lab frame | PyTorch batched (bmm) |
| C | Batch ray-detector intersection → pixel coordinates | PyTorch batched |
| D | Per-peak overlap: rasterize + count against experimental images | Batch C extension (`stage_d_overlap`) |

Stage D processes all M peaks × N detectors in a single C call via `_rasterize.c:stage_d_overlap()`. It uses pre-cached uint8 binary images (`ImageData.get_binary_numpy()`) and implements triangle overlap, pixel-radius search, contiguity validation, and Welford quality aggregation entirely in C. A Python fallback path is available when the C extension is not compiled.

## Citation

S. F. Li and R. M. Suter, "Adaptive reconstruction method for three-dimensional orientation imaging", *Journal of Applied Crystallography*, 2013.
