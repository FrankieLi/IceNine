# IceNine Python/PyTorch Implementation

Python port of IceNine for synchrotron X-ray diffraction — forward simulation and reconstruction of crystal grain orientations. Forward simulation validated pixel-exact against the C++ implementation.

## Installation

```bash
cd icenine_py
uv sync --extra dev
```

`uv sync` builds `.venv` from the tracked `uv.lock` with the Python in `.python-version` (3.9), so
everyone gets the same numpy/torch/scipy versions (2.0.2 / 2.8.0 / 1.13.1). Do not use
`uv pip install`: it ignores the lock, resolves newer versions, and the golden bit-identity tests in
`tests/test_findoptimal_refactor.py` then fail (they check the recorded versions and say so).

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
| `differentiable_cost.py` | Gradient-capable cost infrastructure: `SparseImageStack`, `MultiScaleImageStack`, `DifferentiableCostFunction` |
| `orientation_nn.py` | Toy orientation-NN support: ROI peak set, windowed renderer (`define_roi_set`, `render_local_windows`, `spot_overlaps_grid`), perturbation samplers, dataset, quaternion loss and offset Gaussian NLL (`cholesky_from_raw`, `gaussian_nll_loss`), `mse_deg_loss` and `decoupled_nll_loss` (MSE on the mean, NLL on the covariance at stopgrad(mean)), `split_by_voxel` (validation voxels, one per r_perp stratum) |
| `orientation_eval.py` | Stage 0 evaluation: batched float64 `BatchedObserver` (presence, frame, spot vertices for many candidate offsets), `ExactBayes` posterior, `lit_pixel_set`, per-axis `error_summary`, rotation-vector helpers; frame-coded window helpers `WindowSpec`, `render_windows`, `decode_windows`, plus `combine_windows` / `make_realistic_windows` for distractors and the realism layer |
| `toy_orientation_model.py` | `ToyOrientationNet` (v0 quaternion head), `ToyOffsetNet` (offset + Cholesky covariance head) and `PeakSetNet` (shared per-peak encoder + explicit measurement features + mean/max/sum pooling; any number of peaks, order-invariant, per-voxel context), `GNLayerNet` (learned Gauss-Newton layer, current architecture) and `FrameProbeNet` (frame-only diagnostic net) |
| `orientation_baselines.py` | Stage 2 non-learning baseline: `extract_measurements` (frame index and lit-pixel centroid per spot; optional per-detector mask) and `CentroidGaussNewton` (plain or Huber-robust `huber_c`, `information` = J^T W J, `solve_linear` = one undamped step from nominal) |

Architecture (parallax) additions to the toy orientation NN (see MIGRATION_HISTORY.md, "Toy Orientation NN — Architecture (parallax)"):

- `toy_orientation_model.GNLayerNet` (`--arch gn`): a learned Gauss-Newton layer. A shared per-peak encoder emits reliability weights and measurement corrections; the measurement (centroid/frame minus the *exact* nominal prediction) and the Jacobian from the per-peak context are combined by pooled normal equations (closed-form 3x3 solve, MPS-safe), with the covariance `D A^-1 D`. `--gn-iters T` unrolls IRLS rounds, `--pairing` mixes each entry's encoding with its other-detector partner's through a zero-initialised linear residual MLP (trainable since the 2026-10-01 fix; before it, the ReLU after the zero-initialised layer made that MLP inert) and gives the head the partner entry's residual and a has-partner flag (see `docs/orientation_nn_design.md` Section 3.4). `PeakSetNet` takes the same exact-nominal measurement with `--subpixel`.
- `orientation_eval.nominal_offsets` / `pair_index` build the per-peak exact nominal offsets and the detector-pair index; `scripts/make_dataset_aux.py` writes them as a sidecar (`--aux`) for an existing dataset. `scripts/arch_diagnostics.py` runs the per-detector / conditioning diagnostics. `scripts/make_dis_val.py` rebuilds the 200-sample validation set used to choose the Huber threshold (`toy_orientation_arch_dis_val.pt`).
- Realism: `orientation_eval.render_distractor_windows` (neighbour-voxel and Sigma3-twin spots as a separate layer), `combine_windows` (target pixels always win), `RealismConfig`, `make_realistic_windows`, `make_realistic_dataset` (missing spots, edge-pixel jitter, hot pixels, spurious blobs). "Realistic" throughout means deliberately simulated complexity added to exact synthetic windows: genuine neighbour-voxel and twin spots overlapping the target plus detector noise, labels always exact (`docs/orientation_nn_design.md` Section 2.7; earlier called "corruption", old names kept as aliases). The generator's `--neighbors N --twin` writes the layer as `dis_windows`; the trainer takes `--realistic-train`, `--eval-variants`, `--extra-variant`; `scripts/gauss_newton_baseline.py` takes `--realistic`, `--huber`.
- Trainer: `--mask-padding` (keep zero-padded peak entries all-zero on realistic data; off in the reported runs), `--ema` / `--checkpoint {best,ema,last}` (EMA weights instead of noisy best-epoch selection); the results json gains a `summary` with corr(per-voxel error, r_perp). `scripts/summarize_results.py` averages seeds into tables. `scripts/summarize_arch_step4.py` scores plain/Huber Gauss-Newton and the nets on clean and realistic test sets from saved predictions (same code, same samples), with per-seed medians and the error-vs-r_perp correlation.

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

## Differentiable Cost Function

`differentiable_cost.py` provides gradient-capable infrastructure for orientation optimization experiments. It wraps the same physics as `VoxelCostFunction` but uses `F.grid_sample` (bilinear) instead of hard binary rasterization, making it differentiable via PyTorch autograd.

### Classes

**`SparseImageStack`** — memory-efficient image storage. Stores only (row, col) pixel coordinates instead of dense float32 tensors. For typical binary diffraction data: ~12.5 KB (ThreeVoxels) vs 5.6 GB dense — over 400,000× smaller.

```python
from icenine.differentiable_cost import SparseImageStack

image_stack = SparseImageStack.from_image_directory(
    directory="ScatteringData_Python",
    basename="3Grains.sim", ext="d", serial_length=5,
    n_omega=180, n_detectors=2, num_rows=2048, num_cols=2048, binary=True,
)
print(f"{image_stack.memory_bytes / 1024:.1f} KB")  # 12.5 KB
```

**`MultiScaleImageStack`** — multi-resolution image pyramid with optional omega blending. Downsamples by max-pooling at factors [4×, 8×, ...]; optionally blends each frame with its ±`omega_window` neighbors to widen the angular basin.

To construct multiple omega_window variants efficiently (avoids re-densifying all frames per variant):

```python
from icenine.differentiable_cost import MultiScaleImageStack

# Densify downsampled stacks once
shared_ds = MultiScaleImageStack.build_shared_base(image_stack, [1, 4, 8])

# Reuse across omega_window variants (~2 GB total vs ~9-11 GB without sharing)
ms_ow0 = MultiScaleImageStack(image_stack, [1, 4, 8], omega_window=0, _prebuilt_downsampled=shared_ds)
ms_ow1 = MultiScaleImageStack(image_stack, [1, 4, 8], omega_window=1, _prebuilt_downsampled=shared_ds)
```

**`DifferentiableCostFunction`** — drop-in replacement for `VoxelCostFunction` with autograd support. Uses centroid point sampling (one `grid_sample` per peak) instead of full rasterization — ~50× fewer samples, fully differentiable.

```python
from icenine.differentiable_cost import DifferentiableCostFunction
import torch

diff_fn = DifferentiableCostFunction(
    simulator=simulator, detector_list=detector_list, range_map=range_map,
    image_stack=ms_ow1, sample=sample, structure_list=structure_list,
)

# Evaluate at scale=2 (8× downsampled)
R = torch.eye(3, requires_grad=True)
info = diff_fn.evaluate(R, vertices, phase_index=0, scale=2)
print(f"quality={info.quality:.4f}")  # torch.Tensor with grad_fn
info.cost.backward()                  # backprop through cost
print(R.grad)
```

### Orientation Optimization Results

#### Experiment Setup

All gradient/optimizer benchmarks use a controlled orientation-recovery protocol applied to two datasets:

**Datasets:**
- **ThreeVoxels** (`Examples/Example2.ThreeVoxels/`): 3-voxel copper polycrystal. All 3 voxels used (hard quality > 0.1). Ground truth Bunge Euler angles: (355.4°, 5.2°, 29.3°), (155.4°, 45.2°, 29.3°), (356.7°, 3.7°, 328.4°).
- **ManyGrains** (`Examples/Example2.ManyGrains/`): 500-voxel copper polycrystal. 20 voxels selected by scanning up to 500 candidates and drawing a random subset with hard-cost quality > 0.1 (RNG seed=42).

**Physics setup:**
- Crystal: copper (FCC, a=3.61 Å, space group 225, 24-fold cubic symmetry)
- Beam energy: 64.351 keV (monochromatic), direction (0,0,1)
- Max Q: 16 Å⁻¹ (simulation); 8 Å⁻¹ effective during reconstruction (filters weak peaks)
- Eta limit: 86° (azimuthal acceptance)
- Omega range: 0°–179° in 1° steps = 180 frames
- Detectors: 2 flat-panel detectors at ~3.36 cm and ~5.39 cm from sample, 2048×2048 pixels, 14.8 μm pixel pitch
- Each voxel produces ~15–30 observable peaks across both detectors × 180 frames

**Perturbation construction:**
Each voxel has one starting orientation per perturbation size, constructed as:
```
R_start = R_perturb(axis, angle) @ R_ground_truth
```
where `axis` is a uniformly random unit vector drawn from `np.random.default_rng(seed=42)` (one axis per perturbation size, same axis used for all (scale, omega_window, optimizer) combinations for that voxel), and `angle` is the perturbation size in {1°, 2°, 5°}. The perturbation is a left-action rotation that shifts the starting orientation away from ground truth by exactly that geodesic distance on SO(3). The same random axis sequence is reused across all benchmark scripts for direct comparison.

**Cost function variants swept:**
- `scale`: 0 = full-resolution 2048² images (sharp basin), 1 = 4× downsampled 512² (wider basin), 2 = 8× downsampled 256² (widest basin). All gradient benchmarks sweep scales {1, 2}.
- `omega_window`: integer ω±k — each detector frame is morphologically dilated by taking the max with its ±k neighbors in omega before evaluation. Values {0, 1, 2} swept. Widens the angular coverage from 1° to (2k+1)° per frame, broadening the cost function basin.

**Convergence metric:** Geodesic misorientation between final R and ground truth: `arccos((tr(R_gt^T R_final) − 1) / 2)`. Success threshold: < 0.5°.

**Per-benchmark result counts:** Each cell in results tables below counts independent optimizer runs over all (voxel, scale, omega_window) combinations: 3 voxels × 2 scales × 3 ω-windows = 18 configs for ThreeVoxels; 20 voxels × 2 scales × 3 ω-windows = 120 configs for ManyGrains.

---

**Adam gradient descent** (`benchmarks/bench_gradient_optimization.py`): Euclidean Adam on θ ∈ ℝ³ (R = exp(skew(θ))), lr=0.01, n_steps=100. **Fails at all perturbation distances.** Root cause: `grid_sample` bilinear sampling of binary images gives zero gradient in blob interiors — only the 1-pixel blob boundary carries gradient signal. The gradient basin is ~0.3–0.5°, so Adam cannot converge from 1° or larger starting offsets.

**CMA-ES** (`benchmarks/bench_cmaes_optimization.py`, maxiter=500, σ₀=perturbation_rad):

| | Hard cost | Diff cost (scale=2, ω±1) |
|---|---|---|
| ThreeVoxels, 1° | 0/3 | 0/3 |
| ThreeVoxels, 2° | 0/3 | 1/3 (→0.56°) |
| ThreeVoxels, 5° | 0/3 | 0/3 |
| ManyGrains, 1° (20v) | 0/20 | 1/20 (→0.06°) |
| ManyGrains, 2° (20v) | 0/20 | 1/20 (→0.06°) |
| ManyGrains, 5° (20v) | 0/20 | 0/20 |

**Hard cost** fails because the landscape is completely flat outside the ~0.5° basin — CMA-ES receives no signal and drifts to random orientations (40–165° final misorientation). **Diff cost** occasionally succeeds when a run happens to sample the narrow basin, but is mostly trapped by crystal symmetry false optima (Cu has 24-fold cubic symmetry; symmetry-equivalent orientations achieve quality 0.25–0.65 at 40–170° misorientation).

**Riemannian Adam** (`benchmarks/bench_riemannian_optimization.py`, n_steps=100, lr=0.01): Euclidean Adam on θ ∈ ℝ³ has two defects: (1) chart distortion — `d(cost)/d(θ)` mixes the Riemannian gradient with the Jacobian of exp, growing as R drifts from R_start; (2) moment staleness — Adam moments accumulate in a fixed chart anchored at R_start, never parallel-transported. The Riemannian variants project gradients to T_R SO(3) at every step and retract via `R ← R·exp(-lr·Ω_adam)`, staying exactly on SO(3):

| Optimizer | ManyGrains 1° | ManyGrains 2° | ManyGrains 5° |
|---|---|---|---|
| euclidean_adam (baseline, θ ∈ ℝ³) | 24/120 (20%) | 25/120 (21%) | 4/120 (3%) |
| riemannian_adam_manual (pure PyTorch) | 32/120 (27%) | 29/120 (24%) | 4/120 (3%) |
| riemannian_adam_geoopt (geoopt Stiefel) | 37/120 (31%) | 29/120 (24%) | 5/120 (4%) |

Riemannian structure gives +37% more successes at 1° perturbation. Manual and geoopt variants agree closely, confirming correctness. At 5° all methods fail equally — flat landscape dominates.

**Riemannian SGD** (`benchmarks/bench_sgd_optimization.py`): Same setup as Riemannian Adam. At equal lr=0.01, all SGD variants diverge (momentum/nesterov/cosine reach 100–130° misorientation due to gradient spike accumulation). At the fair lr=0.001 (10× smaller, Adam still at 0.01):

| Optimizer | ManyGrains 1° success | 1° mean misori |
|---|---|---|
| riemannian_adam_manual (lr=0.01) | 32/120 (27%) | 1.43° |
| riemannian_sgld (lr=0.001, Langevin noise, T annealing→0) | 18/120 (15%) | 1.66° |
| riemannian_sgd_plain (lr=0.001) | 1/120 (<1%) | 2.45° |
| riemannian_sgd_momentum/nesterov/cosine (lr=0.001) | 0/120 | 20–26° |

SGLD is the best SGD variant due to Langevin noise providing probabilistic exploration. No fixed SGD lr achieves what Adam's `lr_eff ≈ lr/√m̂₂` does: automatic acceleration in flat regions and automatic attenuation of boundary spikes.

**Recommended approaches** (in order of simplicity):
1. **Two-stage MC + gradient polish**: Use existing `AdaptiveMC` to reach within ~0.5°, then apply Riemannian Adam — gradient signal IS reliable inside the basin
2. **Multi-start CMA-ES with symmetry folding**: Restart from all 24 cubic symmetry equivalents, take the best result
3. **Distance field soft images**: Replace binary images with distance transform (distance to nearest bright pixel, float32) — extends gradient signal ~10–20px beyond each blob, eliminating the zero-interior-gradient problem structurally

---

### Comprehensive HP Sweep — Gradient Methods vs. Monte Carlo

`benchmarks/bench_hp_sweep.py` performs a systematic hyperparameter sweep over all gradient optimizer families and the existing `MCOptimizer`, comparing them head-to-head on both the ThreeVoxels (3 voxels) and ManyGrains (100 voxels) datasets.

**Scope:** 195 HP configurations × 3 perturbation sizes × 100 voxels (ManyGrains) = 58,500 independent optimizer runs.
**Also recorded:** subsampled optimization trajectories (angular step size + misorientation from ground truth every 10 gradient steps; every accepted MC move).

#### Optimizer Families and HP Grids

| Optimizer | HPs swept | Configs |
|-----------|-----------|---------|
| `riemannian_adam_geoopt` | lr ∈ {1e-4..0.1} × n_steps ∈ {100,200,500} × β₁ ∈ {0.9,0.95} | 42 |
| `riemannian_adam_manual` | same grid | 42 |
| `riemannian_sgd_plain` | lr × n_steps | 15 |
| `riemannian_sgd_momentum` | lr × momentum (n_steps=200) | 15 |
| `riemannian_sgld` | lr × T_init × n_steps | 45 |
| `mc_optimizer` | max_steps × restarts × angular_step_frac | 36 |

#### Results — Best HP Config per Optimizer

Success threshold: final misorientation < 1°. Results shown for two datasets.

**ManyGrains (100 voxels per cell):**

| Optimizer | Best HP | 1° success | 2° success | 5° success | Time/run |
|-----------|---------|-----------|-----------|-----------|----------|
| riemannian_adam_geoopt | lr=1e-4, n=100, β₁=0.9 | **96%** | 41% | 0% | 0.44s |
| riemannian_adam_manual | lr=1e-4, n=200, β₁=0.9 | 94% | **49%** | 0% | 0.58s |
| riemannian_sgd_plain | lr=1e-4, n=100 | 88% | **52%** | 0% | 0.30s |
| riemannian_sgd_momentum | lr=1e-5, n=200, m=0.5 | 94% | 51% | 0% | 0.59s |
| riemannian_sgld | lr=1e-4, n=100, T=0.01 | 92% | 50% | 0% | 0.29s |
| mc_optimizer | n=3500, restarts=2, step=0.5 | 92% | 40% | **6%** | 3.19s |

**ThreeVoxels (3 voxels — qualitative; cell values are integer counts 0/1/2/3):**

| Optimizer | 1° success | 2° success | 5° success | Time/run |
|-----------|-----------|-----------|-----------|----------|
| riemannian_adam_geoopt | 3/3 | 2/3 | 0/3 | 1.7s |
| riemannian_adam_manual | 3/3 | 2/3 | 0/3 | 1.1s |
| riemannian_sgd_plain | 3/3 | 2/3 | 0/3 | 1.1s |
| riemannian_sgd_momentum | 3/3 | 2/3 | 0/3 | 2.2s |
| riemannian_sgld | 3/3 | 2/3 | 0/3 | 1.1s |
| mc_optimizer | 3/3 | 2/3 | 0/3 | 1.2s |

ThreeVoxels results are consistent with ManyGrains: same qualitative pattern, same optimal lr=1e-4, all methods fail at 5°. (3-voxel counts are insufficient for significance; use ManyGrains for quantitative comparison.)

**Key findings:**

1. **Gradient methods beat MC at small perturbations.** At 1° perturbation, best gradient optimizer (Riemannian Adam geoopt, lr=1e-4) achieves 96% success vs. 92% for MC — and is **7× faster** (0.44s vs. 3.19s). This is because at 1° the starting point is already near the basin, and gradient descent finds it efficiently.

2. **Hard LR cliff.** All Adam variants fail completely at lr ≥ 0.05 (0% success at 1° perturbation). The optimal range is lr ∈ [1e-4, 1e-3]. SGD and SGLD have similar cliffs at lr ≥ 5e-3 and lr ≥ 1e-2 respectively. Staying well below the cliff is the most impactful single HP choice.

3. **More steps don't help.** For Adam at the optimal lr=1e-4: n=100 → 96%, n=200 → 94%, n=500 → 88%. Diminishing returns set in quickly; extra steps can even hurt when the optimizer overshoots.

4. **Only MC succeeds at 5°.** At 5° perturbation, all gradient methods completely fail (0% success, stuck in flat landscape ≫0.5° from basin). MC achieves 6% — also poor, but it's the only method with any 5° successes, because random walk can occasionally land near the basin.

5. **SGD variants are surprisingly competitive.** Plain Riemannian SGD (lr=1e-4, n=100, 0.30s) achieves 88% success and **52% at 2°** (highest among all). The slower, simpler algorithm can do better at 2° because SGD's lack of momentum means it doesn't overshoot narrow basins at that range.

6. **Geoopt vs. manual Adam agree closely.** The geoopt Stiefel manifold retraction and the manual `R ← R·exp(−lr·Ω_adam)` retraction give nearly identical results (96% vs. 94% at 1°), confirming implementation correctness.

#### Trajectory Data

Each run also logs the optimization trajectory. Format: `(step, event_type, angular_step_deg, misori_from_gt_deg, quality)`. Gradient methods record every 10th step; MC records every accepted global improvement and every restart. Trajectory CSV: `benchmarks/hp_sweep_trajectory_{example}.csv` (linked to main CSV by `run_id`).

#### Output Plots

Eight PNG files are generated (four types × two datasets — ThreeVoxels and ManyGrains):

**`hp_sweep_lr_sensitivity_{example}.png`**
Grid of box plots — one subplot per gradient optimizer (MC excluded; it has no LR parameter). X-axis: learning rate (log scale). Y-axis: distribution of final misorientation (°) across all runs at that LR, aggregated over all voxels, perturbations, and n_steps values. Each box shows the median (center line), interquartile range IQR = Q75−Q25 (box edges), 1.5×IQR whiskers, and individual outliers as dots. Reveals the hard LR cliff: distributions shift from narrow and low (converged) to wide and high (diverged) at a specific learning rate threshold. Also shows bimodality — when two modes exist (some runs converging, others failing) at the same LR.

**`hp_sweep_nsteps_sensitivity_{example}.png`**
Grid of box plots — one subplot per optimizer (including MC). X-axis: number of optimization steps (n_steps for gradient methods; max_mc_steps for MC). Y-axis: distribution of final misorientation (°). Box statistics same as above (median, IQR, 1.5×IQR whiskers, outliers). Shows whether more steps improve results: at the optimal LR, gradient methods exhibit flat or worsening distributions beyond n=100, while MC shows the expected steady improvement.

**`hp_sweep_optimizer_comparison_{example}.png`**
Three-panel box plot — one panel per perturbation size (1°, 2°, 5°). X-axis: optimizer family. Y-axis: distribution of final misorientation (°) for each optimizer's best HP configuration (the HP with lowest mean misorientation for that optimizer). Box statistics same as above (median, IQR = Q75−Q25, 1.5×IQR whiskers, outlier dots). Shows the full distribution of outcomes at each optimizer's ceiling performance — distinguishing whether low mean is driven by a tight, reliably converging distribution or by a bimodal mix of successes and failures.

**`hp_sweep_trajectory_{example}.png`**
Step-size trajectory plot. X-axis: event index — the sequential count of recorded optimization events (one event per 10 gradient steps; one event per accepted MC move or restart). Y-axis: angular step size (°) — the geodesic distance on SO(3) between consecutive recorded states. Solid line = median angular step size across all runs for that optimizer at each event index; shaded band = interquartile range (IQR = Q25 to Q75, i.e. the middle 50% of the run distribution). Gradient methods show smooth monotonic decay as the optimizer converges; MC shows irregular bursts — large random steps on accepted improvements followed by smaller steps as the local optimum is refined.

#### Running the Benchmark

```bash
cd icenine_py
uv sync --extra riemannian
uv run python benchmarks/bench_hp_sweep.py --example threevoxels
uv run python benchmarks/bench_hp_sweep.py --example manygrains
uv run python benchmarks/bench_hp_sweep.py --smoke-test --example threevoxels  # quick test
uv run python benchmarks/bench_hp_sweep.py --plots-only --example manygrains   # regenerate plots only
```

### Hybrid Riemannian Adam + MC-Restart Optimizer

The HP sweep showed Riemannian Adam (lr=1e-4, n=100) beats MC at small perturbations while being much faster. `RiemannianAdamOptimizer` (in `icenine/orientation_search.py`) productionizes this as a drop-in replacement for `MCOptimizer` in the `FindOptimal` phase of `AdaptiveVoxelReconstructor`: it uses `DifferentiableCostFunction` (geoopt Stiefel manifold parameter) for gradient signal, `VoxelCostFunction` (hard binary overlap) for convergence decisions, and falls back to MC-style random-restart perturbation when stuck. SVD re-orthogonalization runs after each Adam loop to guard against Stiefel float drift.

**Enabling it:**

```python
from icenine.reconstructor import setup_reconstruction, build_diff_cost_fn

setup = setup_reconstruction(config)
setup.diff_cost_fn = build_diff_cost_fn(setup)
setup.search_params.use_hybrid_optimizer = True
setup.search_params.adam_n_steps = 100
setup.search_params.adam_lr = 1e-4

reconstructor = AdaptiveVoxelReconstructor(setup)
```

**Head-to-head results** (`benchmarks/bench_hybrid_optimizer.py`, 20 voxels × 5 perturbations, ManyGrains):

| Perturbation | Hybrid success | MC success | Hybrid speedup |
|---|---|---|---|
| 0.5° | 100% | 90% | 2.5× |
| 1.0° | 60% | 45% | 1.8× |
| 2–3° | (MC wins, consistent with the HP sweep — gradient methods lack signal this far from the basin) | | |

```bash
cd icenine_py
uv sync --extra riemannian
uv run python benchmarks/bench_hybrid_optimizer.py --smoke-test --example threevoxels  # 3 voxels, 2 perturbations
uv run python benchmarks/bench_hybrid_optimizer.py --example threevoxels               # full run
```

Outputs: `benchmarks/bench_hybrid_{example}.csv` (per voxel/perturbation/optimizer), plus success-rate, wall-time (Adam vs. hard-eval breakdown), and MC-vs-hybrid scatter plots.

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
uv run pytest tests/ -v                                 # all tests (~328 passed, ~34 skipped)
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

## Toy Orientation NN Scripts (`scripts/`)

Prototype pipeline for learning local orientation refinement (status and plan in
[MIGRATION_HISTORY.md](MIGRATION_HISTORY.md)). Requires the Example2.ThreeVoxels files.

```bash
cd icenine_py
uv run python scripts/generate_toy_orientation_dataset.py --n-train 1500 --test-per-bin 30   # train/test .pt (gitignored)
uv run python scripts/exact_bayes_baseline.py --test scripts/toy_orientation_stage0_test.pt --out benchmarks/toy_orientation_stage0/test_bayes.npz
uv run python scripts/train_toy_orientation_nn.py --head offset --train scripts/toy_orientation_stage0_train.pt \
    --test scripts/toy_orientation_stage0_test.pt --bayes benchmarks/toy_orientation_stage0/test_bayes.npz
```

Stages 1-3 (frame-coded windows, baselines, set network, other voxels):

```bash
G="--renderer observer --max-q 8 --detectors all --min-sin-eta 0.3 --frame-half-width 4 --prior-radius 1 --test-magnitudes 0.1 0.25 0.5 1.0"
# one voxel of ManyGrains (r_perp 399 um), or 30 voxels spanning r_perp 0-500 um (24 train / 6 held out)
uv run python scripts/generate_toy_orientation_dataset.py --example manygrains --voxel-index 77 $G --n-train 10000 --tag stage3_far
uv run python scripts/generate_toy_orientation_dataset.py --example manygrains --n-voxels 30 --per-voxel-train 500 --per-voxel-test 10 $G --tag stage3_multi
uv run python scripts/gauss_newton_baseline.py --test scripts/toy_orientation_stage3_multi_test.pt --out benchmarks/toy_orientation_stage3/multi_pred_gauss_newton.npz
uv run python scripts/train_toy_orientation_nn.py --arch set --pool all --head offset --device mps --lr 3e-4 --clip 1 --cosine --epochs 60 \
    --train scripts/toy_orientation_stage3_multi_train.pt --test scripts/toy_orientation_stage3_multi_test.pt \
    --extra "gauss-newton=benchmarks/toy_orientation_stage3/multi_pred_gauss_newton.npz"
```

Multi-voxel datasets store padded windows `(N, M_max, W, W)`, a per-voxel `context` table
`(V, M_max, 14 + n_detectors)` and a `voxel_id` per sample; the training script gathers the context per batch and
reports in-distribution and held-out-voxel tables plus a per-voxel median error vs r_perp.
`--example`, `--voxel-index` and `--n-voxels` select the sample; the baseline scripts read the
example from the dataset. Datasets (up to 3.2 GB) are gitignored; small `npz`/`json` results are in
`benchmarks/toy_orientation_stage{1,2,3}/`.

For multi-voxel data `--val-voxels N` holds out N training voxels (one per r_perp stratum, chosen with
`--seed`) for early stopping instead of a random 10 % of samples; `--loss decoupled` trains the mean
with MSE and the covariance with NLL at stopgrad(mean).

`--arch probe` trains the frame-probe network (`FrameProbeNet`). `scripts/optimizer_baselines.py` runs the MC and Riemannian-Adam baselines at the test perturbation sizes, `scripts/gauss_newton_baseline.py` the centroid Gauss-Newton baseline (reports a convergence status per case), and `scripts/summarize_bayes_npz.py --bayes B.npz --test T.pt --out R.json` tabulates an exact-Bayes run.

Commands for the `GNLayerNet` (`--arch gn`) and Step 4 (distractor/realism layer) experiments, including the data, baseline and summary scripts, are in `docs/orientation_nn_design.md` Section 9.

`train_toy_orientation_nn.py --save-model PATH` saves the final weights, the constructor arguments and the input settings (loaded by `load_model` in `scripts/perturbation_sweep.py`). `scripts/perturbation_sweep.py` is the single-voxel perturbation sweep: 50 random voxels (not among the dataset's 30), each started at a random rotation of angle r (0.05-5 deg) from its true orientation, one-shot and re-centred, on clean and realistic windows; `run` renders and evaluates (CPU worker pool), `summarize` rebuilds the tables, JSON and plot from the raw npz. Results are in `benchmarks/toy_orientation_sweep/` (see `MIGRATION_HISTORY.md`, "Perturbation sweep"); tests in `tests/test_perturbation_sweep.py`.
`scripts/optimizer_sweep.py` runs the existing optimizers (MC, Riemannian Adam, plain/Huber Gauss-Newton) on the same cases and puts their error next to the network's (`optimizer_sweep_summary.py` makes the tables and `perturbation_sweep_vs_optimizers.png`); per-case detector images are built from pixel sets (target spots at the truth, plus the sweep's distractors and realism edits for the realistic variant); tests in `tests/test_optimizer_sweep.py`. `VoxelCostFunction` / `DifferentiableCostFunction` gained an optional `min_sin_eta` filter (default off).
`scripts/findoptimal_sweep.py` runs the multi-level `AdaptiveVoxelReconstructor` (ReconstructQ8 search settings, Q_max 8, `|sin eta| >= 0.3`) on the same voxels and images: A = `reconstruct_voxel` from scratch (no starting guess), B = only its final stage (FindOptimal + VarianceMinimizing) started at each case's perturbed orientation; errors are cubic-symmetry reduced (`findoptimal_sweep_summary.py` makes `findoptimal_sweep_summary.{txt,json}` and `perturbation_sweep_vs_findoptimal.png`). `AdaptiveVoxelReconstructor` gained `refine_from_candidates` (the final stage factored out of `reconstruct_voxel`, behaviour unchanged) and an optional `min_sin_eta`; tests in `tests/test_findoptimal_refactor.py` and `tests/test_findoptimal_sweep.py`.

`scripts/findoptimal_robustness/` (own `README.md`) investigates why the multi-level reconstruction returns CSL relatives of the truth: instrumented runs of `reconstruct_voxel` (recorder hook) on 200 voxels x 3 seeds x clean/realistic, a CSL-relative check after the search (F1), search changes (F1b, F2, F3), and a one-pass candidate classifier (features, voxel/grain-disjoint evaluation, end-to-end reranking). `AdaptiveVoxelReconstructor` gained the optional attributes `recorder`, `keep_fraction`, `keep_union_discrete`, `n_q_start_offset`, `global_pixel_radius`, `extra_candidates` and `rank_key` (all default to the old behaviour). Results: `benchmarks/findoptimal_robustness/`; report: `docs/findoptimal_robustness_report.md`; tests in `tests/test_findoptimal_robustness.py`. scikit-learn (dev extra) is needed for the classifier scripts.

`scripts/nn_hybrid/` runs the hybrid network -> FindOptimal study on the same sweep cases (own header docstrings): `run.py` (pipelines H0 FindOptimal alone, H1/H3 net x1/x3 -> `refine_from_candidates`, H3c covariance-sized box, H3m net -> MC, HG Huber GN -> net -> FindOptimal; `pilot`, `run`; `run_all.sh` runs the full set one pipeline at a time), `summary.py` (metrics per radius x variant x pipeline, success criteria, decision points), `stage_timer.py` (context-manager stage timer plus per-instance `VoxelCostFunction.evaluate` counters, installed by patching module attributes; shared with the run-time profiling). Caches in `scripts/nn_hybrid/cache/` (gitignored), results in `benchmarks/nn_hybrid/`; tests in `tests/test_nn_hybrid.py`. Nothing in `icenine/` changed. Note that `findoptimal_b_raw.npz` and the hybrid raws store the voxel axis sorted by voxel index, the sweep raws in sweep order.

`scripts/coarse_proxy/` tests a low-Q_max cost proxy as the pruning key of `reconstruct_voxel` (the `rank_key` hook) on the FindOptimal-robustness cases (own header docstrings): `lowq_dataset.py` (F-lowQ features for every E2 candidate at Q = 4 and 5, via the new optional `q_max` of `findoptimal_robustness/features.py`; `qlevels` records the |q| families), `labels.py` (the basin-cost regression target `y_bcost`, the y3 / y1 / level-matched classification labels, their validation), `models.py` (grain-disjoint 4-fold GBT models, AUC / pruning-recall / final-precision / Spearman tables, decision points D1-D3, optional retraining with harvested candidates), `endtoend.py` (full `reconstruct_voxel` runs with the proxy key, and `harvest` of the proxy run's own candidates), `timing.py` (single-worker cost accounting), `timing_batch.py` (batched against per-candidate low-Q pass, single worker, preflight gated), `summary.py` (all tables, criteria and decisions; prints "not evaluated" for anything not run). `findoptimal_robustness/f1_run.py` gained the optional `--cache-root` and `--source-seeds`; `eval_split.py` measures where the proxy run's extra cost evaluations go (stage-timer counters, bit-identity asserted). Caches in `scripts/coarse_proxy/cache/` (gitignored), results in `benchmarks/coarse_proxy/`; tests in `tests/test_coarse_proxy.py`. `FeatureExtractor.features_batch` (findoptimal_robustness/features.py) scores many candidates in one vectorised pass with features exactly equal to `features` (per-candidate path stays the default; `endtoend.py --batched` opts in; `endtoend.py` also stores the batch size of every `rank_key` call). `keep_eighth.py` runs the T4 summary (keep 1/8 and 1/6 against the baseline and the keep-1/4 row, proxy cost at the batched per-call cost, the "useful" verdict) and the single-worker interleaved timing of a voxel subset (preflight gated). Nothing in `icenine/` changed.

`scripts/profiling/` measures where the run time goes and whether the network or the proxy reduces it (own header docstrings; `prof_common.py` holds the isolation record, the interleaving order and the stage instrumentation `Instrument`, which installs the `stage_timer.py` wrappers by patching module attributes where they are looked up; `prof_u0.py` times the full-reconstruction pipelines (baseline, F1, F1b, E2 rerank, proxy rerank, proxy + F1) single-worker and interleaved on the E0 images, asserting bit-identity with the stored runs; `prof_seeded.py` times H0 / H1 / H3 / MC-told-r / HG per case (U1 / U2), with the network stage on one case at a time and "measured" vs "production-equivalent" (without harness rendering) time; `prof_profile.py` runs cProfile (5 runs per pipeline), torch.profiler (`record_function` around the net sub-stages) and the forward pass on MPS vs CPU; `prof_overhead.py` measures the wrapper overhead; `prof_summary.py` writes the tables, the contention factor and the "helps" verdicts; `run_all.sh` chains the timing runs). Caches in `scripts/profiling/cache/` (gitignored), results in `benchmarks/profiling/`; tests in `tests/test_profiling.py`. Nothing in `icenine/` changed.

`scripts/finisher_diagnosis/` diagnoses why `refine_from_candidates` ends above the truth's cost on realistic sweep cases (own header docstrings): `diagnose.py` (`pilot`, `run`; re-runs the finisher with the Task 1 seed through the unmodified reconstructor with `LoggedMC`, a `MCOptimizer` subclass that records the stopping rule and iteration, then measures the cost along the geodesic result -> truth, the cost granularity at the truth and result, and continuations with longer runs / smaller steps / smaller boxes), `diag_summary.py` (tables and `summary.json`). Caches in `scripts/finisher_diagnosis/cache/` (gitignored), results in `benchmarks/finisher_diagnosis/`; tests in `tests/test_finisher_diagnosis.py`. Nothing in `icenine/` changed.

`scripts/cost_sensitivity/` (finisher/MC study, Phase A; own header docstrings): `landscape.py` (`pilot`, `run`; the cost landscape around the truth on the T5 cases, clean and realistic variants of the same voxel/radius/direction: 400 random directions x 10 radii from 0.0005 to 0.1 degrees plus the finisher's step), `resolution.py` (centroid-quantisation scale of a voxel's orientation from the geometry (a scale, not a bound); `fisher_quantisation`, `crb`), `summary.py` (plateau statistics, sampled minimum, cost-vs-angle correlation, gap-closure statistics -> `benchmarks/cost_sensitivity/summary.json` and `tables.md`). Cache in `scripts/cost_sensitivity/cache/` (gitignored); tests in `tests/test_cost_sensitivity.py`. Nothing in `icenine/` changed.

`scripts/checks/` holds the numerical checks behind the derivations in `docs/`. Results
of the Stage 0 run are in `benchmarks/toy_orientation_stage0/`.

## Design Notes (`docs/`)

Derivations that underpin ongoing work, written in Markdown with LaTeX math
(render with `pandoc <file>.md -o <file>.pdf`, without `-N`):

- [`docs/omega_peak_width_derivation.md`](docs/omega_peak_width_derivation.md) —
  angular width of a diffraction peak in a rotation scan (Δω ≈ α/|sin η|) and
  the sensitivity of the Bragg crossing to orientation (1/|sin η|).
- [`docs/nn_inverse_problem_formulation.md`](docs/nn_inverse_problem_formulation.md) —
  what a neural network for orientation refinement should approximate given
  frame/pixel-integrated data (the posterior, not an inverse), the spot-motion
  Jacobian, and angular-resolution estimates in terms of pixel size, frame width,
  number of frames and distance from the rotation axis.
- [`docs/orientation_nn_design.md`](docs/orientation_nn_design.md) — design of the
  current toy orientation network: data pipeline, `GNLayerNet`, training, evaluation, tests and results.

- [`docs/batched_lowq_pass.md`](docs/batched_lowq_pass.md) — how the batched low-Q feature pass
  (`FeatureExtractor.features_batch`, `--batched`) works, its two bit-identity pitfalls, the equality
  test and the timing at batch 1/50/200.

- [`docs/findings_orientation_search_2026-10.md`](docs/findings_orientation_search_2026-10.md) —
  consolidated findings, recommendations, decisions and open questions of the orientation-search work
  of 2026-10-01 to 2026-10-07 (network, FindOptimal robustness, hybrids, cost proxy, profiling,
  follow-ups).

See [`docs/orientation_nn_design.md`](docs/orientation_nn_design.md) for the status and
current results of this work, and [MIGRATION_HISTORY.md](MIGRATION_HISTORY.md) ("Toy
Orientation NN — Theory Phase" onward) for the step-by-step record.

## Developer tooling

Rules and checks that used to be repeated in LLM prompts live in code. Paths are from the repository root; run Python through `uv run --project icenine_py python ...`. The git pre-commit hook is committed but **not activated**; see `scripts/dev/install_hooks.sh`.

| Script | Purpose and usage |
|---|---|
| `.claude/hooks/check-bash-command.sh` | Claude Code PreToolUse hook: denies bare `python*`/`pytest`/`pip*` (also behind `nohup`, `nice`, `env`, `time`, `sudo`, subshells), `git add -A/--all/./://*/-u`, `git commit -a/-am/--all/--no-verify`, and `ALLOW_*` assignments (handles `&&`, `;`, pipes, newlines, backslash continuations; ignores quoted text and heredoc bodies). Wired in `.claude/settings.json`. |
| `.githooks/pre-commit` | Runs `scripts/dev/precommit_check.py` on staged changes: forbidden paths (`*.pt`, `*.pkl`, `*.joblib`, `cache/`, big `.npz`, `CLAUDE.md`/`.claude/**`), files over 10 MB, absolute home paths, Black on new `.py` files, lines over 100 characters, "realistic" terminology. Overrides: `ALLOW_CLAUDE_CONFIG=1`, `ALLOW_LARGE=1`, `ALLOW_ABS_PATHS=1`, `noqa: realistic`. |
| `scripts/dev/install_hooks.sh` | `git config core.hooksPath .githooks`; prints how to undo it. |
| `scripts/dev/start_task.sh <task> [--push]` | Create `feature/<parent>-<task>` off the current feature branch and record its parent. |
| `scripts/dev/finish_task.sh [--no-tests] [--yes] [--push] [-m msg]` | Run the Python suite and/or C++ build if relevant files changed, merge into the parent with `--no-ff`, delete the task branch. |
| `icenine_py/scripts/common/stats.py` | `wilson`, `mcnemar_exact`, `paired_discordant`, `win_rate`, `reorder` (id-checked re-indexing), `misorientation_deg_cubic`. Usage: put `scripts/common` on `sys.path`, `import stats`. |
| `icenine_py/scripts/common/doc_tables.py` | `markdown_table(rows, columns, formats)` and `write_tables(path, {name: md})` (marker-delimited blocks). |
| `scripts/dev/sync_doc_tables.py --doc D.md --tables T.md [--check]` | Copy generated blocks into the matching `<!-- table:NAME -->` markers; `--check` exits 1 on a difference. |
| `scripts/dev/audit_numbers.py --doc D.md --section "Heading" --sources FILES... [--strict]` | List numbers in a doc section that no source value matches (precision-aware; x100, /100 and sibling k/n ratios only for percentages and fractions). `-v` shows the nearest source value and file:key; `--per-file` requires each paragraph or table row to match within one source file. Advisory; prints the chance-match rate per precision bucket. A match is not verification: low-precision decimals match almost anything. |
| `icenine_py/scripts/common/preflight.py`, `scripts/dev/timing_preflight.py [--require] [--json out]` | Record load, power source, thread settings and busy processes before a timing run; `require_quiet` refuses a busy machine. |
| `scripts/dev/job_status.sh`, `scripts/dev/checkpoint.sh [note]` | Read-only snapshot of running jobs, recent logs, `.done` markers and git state; the checkpoint saves it to `.claude/checkpoints/`. |
| `scripts/dev/new_todo.sh <slug> "<title>" ["<subtitle>"]` | Scaffold `icenine_py/docs/todo_<slug>.md`. |
| `.claude/agents/implementer.md`, `.claude/agents/code-reviewer.md` | Slim LLM agents: an implementer (Sonnet) for approved plans, and the reviewer with a "Claims audit" section that runs the scripts above first. |

Tests: `tests/test_dev_tooling.py`.

## Citation

S. F. Li and R. M. Suter, "Adaptive reconstruction method for three-dimensional orientation imaging", *Journal of Applied Crystallography*, 2013.
