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

`uv sync` is exact: it uninstalls every extra you do not list. Use `uv sync --extra dev` for the library and the tests (the geoopt tests skip), and `uv sync --extra dev --extra riemannian --extra benchmarks` for the study scripts and the hybrid optimizer.

Dependencies: numpy, torch, pymatgen, scipy, cma (see `pyproject.toml`).

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
| `diffraction_core.py` | Scattering omega calculation, reflected rays (PyTorch batched); `BraggConstants` and `scattering_omegas_fixed_beam` (S1: the Bragg solve with the call-independent tensors hoisted, equal to `get_scattering_omegas_torch` on observable rows) |
| `peak_filters.py` | Eta-angle acceptance filter with Lorentz-polarization correction (`batch_eta_filter()` for vectorized ops) |
| `simulation_range.py` | Omega range system for discontinuous data collection wedges |
| `geometry.py` | Euler conversions, Plane/Ray classes |
| `symmetry.py` | Crystal symmetry operations (pymatgen wrapper) |
| `constants.py` | Physical constants |
| `file_io.py` | Detector file, crystal structure file, and omega file I/O |
| `reconstructor.py` | Reconstruction orchestrators — serial, adaptive, and BFS spatial propagation |
| `orientation_search.py` | Discrete grid search + zero-temperature MC optimization (`MCOptimizer.optimize` is the C++ `RandomRestartZeroTemp`: fixed-length blocks, restart about the initial orientation, stop on successive failed blocks; the VarianceMinimizing restart follows C++ too); opt-in `CMAOptimizer` (local CMA-ES on a rotation vector about the start, `SearchParameters.local_optimizer = "cma"`) and the hybrid `RiemannianAdamOptimizer` |
| `cost_functions.py` | Overlap computation between simulated projections and experimental data (batched Stages A-C + sequential Stage D); `frame_of_omegas` (stage-A lookup) and `FrameTable` (S1: the binary frames for the C stage, one table per cost function, filled on first use) |
| `_rasterize.c` | CPython C extension for fast triangle rasterization and pixel overlap (Sutherland-Hodgman + Bresenham) |
| `experimental_data.py` | Load experimental detector images for reconstruction |
| `sampling.py` | SO(3) uniform sampling via Sukharev grids (Yershova & LaValle) |
| `batched_screen.py` | Batched level screen (default on since 2026-10-13) (`ScreenEvaluator`): many (FZ point, local-grid) orientations per call, exact (`BatchedScreen`) or with C++'s first-order omega interpolation from the FZ centre (`CppOmegaInterp`) |
| `batched_cost.py` | Batched local cost (`BatchedCostEvaluator`, S3): many orientations per call; the CPU compact variant equals `VoxelCostFunction.evaluate` in every `OverlapInfo` field, the fixed-shape variant runs on `mps` or `cuda` (not bit-equal) for batches of at least `GPU_MIN_BATCH`; default rule, `CppOverlapRule` and (S1 Step 4) `CppEdgeRules`; shares the cost function's `FrameTable` |
| `batched_quick_mc.py` | Lockstep quick MC (`run_batched_quick_mc`, config key `BatchedQuickMC`, opt-in): all candidates of a level step together through one batched cost call, with draw option B (a fixed slice of 33 uniforms per candidate), so results differ from the sequential quick MC's |
| `candidate_dedup.py` | Candidate overlap and de-duplication (Task C, config key `CandidateDedup`, opt-in): `neighbour_lists` / `greedy_nms` (symmetry-reduced radius search of orientation matrices with a `cKDTree`: representative with the largest |w|, boundary copies, both quaternion signs; same groups as `scripts/findoptimal_robustness/csl.py::dedup_orientations`, tested), `component_labels`, `merged_clique_labels`, `DedupSpec` (the value of the key) and the level-loop stages `keep_before_quick_mc` (clique filter, dedup before the quick MC) and `distinct_order` (kept quarter counted over distinct candidates). `SearchCandidate.clique` is the index of the FZ point that produced a candidate (default -1) |
| `spot_pairing.py` | Geometric spot pairing across two detectors (Task D, no search over orientations): `line_origin` (the line through a detector-0 and a detector-1 spot centroid meets the layer z = 0; origin and its linearised tolerance), `SampleHull`, `PairingParams` (1 pixel centroid tolerance on each detector, optional arm A terms on detector 1: registration, voxel size x (L1/L0 - 1), arc motion r d_omega / 2; optional ring-angle and shape gates), `pair_frame` (batched over all spot pairs of a frame, greedy assignment), `attribute_peaks` (`cKDTree` radius attribution of paired peaks to voxel centroids). Pure geometry; images are read by the scripts. `ForwardSimulation.peak_hook` (default `None`) is the opt-in provenance hook used by `scripts/phase_d/render_provenance.py` |
| `smooth_cost.py` | Smooth ranking costs (Task E, config key `SmoothCostRank`, opt-in): `SmoothSpec` (the value of the key), `SmoothCost` (batched, one image set; lit-pixel `cKDTree`s built through S1's `FrameTable`, Bragg solve through `scattering_omegas_fixed_beam`) with the control arm S (binary pixel test at the three projected vertices and the centroid) and the arms E1 (distance-to-lit-pixel kernel), E2 (omega blending: largest frame credit times a kernel of the distance to that frame's omega interval), E3 (soft peak weights at the detector edge, eta limit and omega range ends; credit over the detectors a spot reaches), E4 (observed-ray residual: spots paired across the detectors with `spot_pairing`, attached to the voxel by Task D's footprint rule; robust angle-plus-omega residual, autograd gradient in `e4_cost_torch`); `evaluate_many`, `evaluate_specs` (several specs, shared queries), `level_widths` (the `auto` widths of a level), `SmoothRanker` / `make_ranker` (the reconstructor's ranking hook; a new `SmoothCost` is built in every `reconstruct_voxel` call, 16.6 s on the full sample: no cache per image set yet, Task G) |
| `outward_search.py` | Searches that start from an orientation a few degrees from the answer (Task F; library only, not wired into the reconstructor, no config key): `RayGaussNewton` (F4: Gauss-Newton on Task E's observed-ray residual, the E4 matching of each predicted peak to the nearest observed ray attached to the voxel; exact Jacobian by autograd through `SmoothCost.predict`, one batched backward pass; robust loss `r^2 / (r^2 + tau^2)` by reweighting, Levenberg-Marquardt damping and a trust region; the strict cost is evaluated at the start and the end only and decides `accepted`), `ShellSearch` (F5: deterministic equal-area shells on SO(3) about the start with spacing `spacing_deg`; per shell a batched pixel-radius screen, then E4 on the best `n_smooth`, then the strict cost on the best `n_strict`; stops at the first strict cost under `accept_cost`), `shell_rotations` / `shell_size` / `shell_directions` (the shell geometry), `screen_fn_from` / `smooth_fn_from` (callables from `ScreenEvaluator` and `SmoothCost`), `angle_deg`, `rotate_left`. Both count strict, screen and smooth evaluations separately. |
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
| `LocalOptimizer` | Optional, Python only: `mc` (default; C++ parity) or `cma` (local CMA-ES for the refinement stages) | `mc` |
| `CMASigma0` | Optional: CMA initial step (degrees of rotation vector); used only with `LocalOptimizer cma` | 0.2 |
| `CMAMaxEvals` | Optional: cost evaluations per CMA run, start included | 1000 |
| `CMAPopSize` | Optional: CMA population (0 = the `cma` package default, 7) | 0 |
| `CMANeighborMaxEvals` | Optional, `cma` only: CMA budget of a BFS neighbour fit; seed voxels keep `CMAMaxEvals`. **In `cma` mode this changes the Phase C BFS (neighbours 1000 -> 250 evaluations)**; set 1000 to keep the Phase C budget | 250 |
| `CMARetrySigma0` | Optional, `cma` only: a rejected BFS neighbour is retried once with this wider start step (degrees); 0 = no retry. **On by default in `cma` mode**, so the Phase C `cma` BFS is not reproduced unless this is 0 | 1.5 |
| `BFSRevisitRefit` | Optional, Python only: 1 = an expansion that reaches a voxel left REFIT fits it again from its own orientation (the multi-client behaviour of the C++ run) | 0 |
| `BFSRevisitMax` | Optional: local fits per voxel while it is REFIT (revisit and restart pass together), >= 1; with `BFSReseed` on the counters of a cluster are reset before each reseed, so the cap applies per reseed | 3 |
| `BFSRestartPass` | Optional, Python only: 1 = after the BFS, FITTED voxels bordering REFIT voxels push their orientation onto them (C++ restart semantics) | 0 |
| `BFSReseed` | Optional, Python only: 1 = when no NOT_VISITED voxel is left, each connected cluster of REFIT voxels gets one full no-start search on one voxel (nearest the cluster centroid, ties lowest index); accepted -> expansion with revisit; a rejected full search (ordinary seed or reseed) makes that voxel failed (never searched again, unresolved) but does not close its cluster, which is reseeded on another voxel. Off keeps the C++ seed rule and the output unchanged | 0 |
| `CppOverlapRule` | Optional, Python only: 1 = the literal C++ `CountQualifiedPeaks` overlap rule (a lit, non-overlapping later detector invalidates a peak's overlap) in the cost function. (The voxel centroid, as C++ `VoxelToCenter`, is projected by the pixel-radius-3 coarse screen for every rule, with or without this key since 2026-10-12; vertex 0 is only the explicit opt-out `SearchParameters.cpp_center_screen=False` / `VoxelCostFunction(center_screen=False)`, not a config key.) Off keeps the original Python rule. On makes the level-0 candidate counts and cost values match C++ (0.995-1.003x and no -0.013 offset) but lowers E0 accuracy (see MIGRATION_HISTORY "Level-0 candidates"); needs `icenine/_rasterize.c` rebuilt (`HAS_CPP_RULE`; from `icenine_py/`: `uv run --no-project --with setuptools --with numpy python setup.py build_ext --inplace`, since `uv sync --reinstall-package icenine` does not rebuild it), otherwise the slow Python path is used with a one-time warning | 0 |
| `CppEdgeRules` | Optional, Python only, opt-in: 1 = three C++ conventions at the edges in the cost function: an omega below the first frame has no frame (the default lookup truncates toward zero and compares it with frame 0, about 0.3% of the peaks of a 180 x 1 degree scan), a detector coordinate below 0 is pixel -1, and the detector-lit rule is exactly C++'s (`OverlapInfo.h:501-514`): a detector is lit iff the pixel overlap is > 0 or one of the three vertices is in bounds. There is no "pixels produced" clause, so a triangle that covers image pixels with all three vertices out of bounds stays unlit, and a triangle that produced no pixel (last image row or column) with a vertex in bounds is lit (the default lights on any produced pixel and is unchanged). Exact for the cost and the peak overlap; the pixel counts (`pixel_on_detector`) are the rasterizer's and need not match C++ (they match on 1,985 of 2,230 orientations of the complete cost). With the key the Python peak overlap and peak-on-detector counts equal C++ on 2,229 (screen) and 2,230 (complete cost) of the 2,230 orientations of the Task A fixed-orientation test, the same as with the earlier superset rule (a rerun gave identical rows: no peak of the test has a triangle with pixels and no vertex in bounds; MIGRATION_HISTORY "Task A" and `benchmarks/eval_reduction/A/step0/`; the findings document is written in Task G); not studied together with `RayOverlapRule`; needs the `_rasterize` extension rebuilt after 2026-10-13 (`HAS_VERTEX_LIT` = 2; an older build warns and falls back to the slow path). Independent of `CppOverlapRule`. | 0 |
| `CppOmegaInterp` | Optional, Python only, opt-in: 1 = the level screen takes the omegas of the local-grid orientations from the FZ centre as C++ `LocalAngularInterpProcess` does (`icenine/batched_screen.py`): the observable reflections and their two omegas are solved once per FZ point, each local-grid rotation `L` rotates the reflection direction exactly and shifts the omega to first order, `omega + (Rz(omega)[0].G' - x_lab)/(Rz(omega)[1].G')`; a shift above 1 rad drops the peak; the observable set is frozen at the centre. With `CppOverlapRule 1` and `CppEdgeRules 1` the per-FZ-point level-0 pass sets equal the C++ trace on 4,885 of 4,886 points (4,466 with exact omegas). The re-cost, quick MC and later stages are unchanged. Falls back, with a `RuntimeWarning`, to the per-orientation screen at a level that uses the ray rule. Independent of the other keys. | 0 |
| `BatchedScreen` | Optional, Python only, default 1 (2026-10-13; 0 = the per-orientation screen): 1 = the level screen scores 32 FZ points x the whole local grid per call (`SearchParameters.screen_chunk`) instead of one orientation per call. Same arithmetic as the unbatched screen, unless `CppOmegaInterp` is 1; no difference detected on the level-0 screen of voxel 15901 (identical `peak_overlap` on all 43,974 pairs), while a float32 last-bit difference for another batch shape remains possible (`icenine/batched_screen.py` docstring). Needs the `_rasterize` extension. Peak RSS of one no-start search of voxel 15901 (default configuration, one run each): 880.5 MiB with the key, 872.3 MiB without (+8.2 MiB), same evaluations and result. Where the cost function is not supported (ray rule at a level, no `_rasterize` extension) the reconstructor falls back to the per-orientation screen, with a `RuntimeWarning` only if `BatchedScreen` was written in the config or `CppOmegaInterp` is 1, silently under the implicit default. See MIGRATION_HISTORY "Task A". | 1 |
| `BatchedQuickMC` | Optional, Python only, opt-in, default 0: 1 = the quick MC (10 trial steps with 5 restarts, on every candidate of a level) runs all candidates of the level in lockstep, one batched local-cost call per step (`icenine/batched_cost.py`, `icenine/batched_quick_mc.py`), instead of one candidate at a time. **Not bit-equal to 0 by design** (draw option B: candidate i owns the fixed slice of 33 uniforms at stream positions 33 i .. 33 i + 32, whatever the other candidates' outcomes; the sequential draws depend on how many earlier blocks failed), so compare wrong counts, not results. Counted evaluations are the same (12 per candidate). Measured (contended machine): quick MC 4.5-4.7 times cheaper per evaluation on the CPU, 1.7-1.9 times shorter no-start search; E0 clean wrong 107/600 against 111/600 with 0 (no detectable difference). Used for a plain `MCOptimizer` with one block (all shipped configs); under the ray rule only the levels that use it (levels 0-2 by default, each with its own ray cost function) run the sequential quick MC and the later levels, which use the strict local cost, run in lockstep (mixed mode); a build without the `_rasterize` extension (or one older than `CppOverlapRule`) and study subclasses fall back to the sequential path at every level; each fallback gives a `RuntimeWarning`. `CppEdgeRules` runs in lockstep too (S1 Step 4; its results differ from a run without the key, as the edge rule changes the cost) and falls back only with an extension older than `HAS_VERTEX_LIT` 2. Needs the `_rasterize` extension. See MIGRATION_HISTORY "S3: batched quick MC". | 0 |
| `SearchDevice` | Optional, Python only, default `cpu` (`cpu`, `mps`, `cuda`): the device of the batched quick-MC cost, used only with `BatchedQuickMC 1`. On `cpu` the batched cost equals `evaluate`; on a GPU (batches of at least 64 orientations) it does not (the last bits of the transcendental functions), so GPU results are compared by wrong counts: on 50 E0 runs `mps` gave the same wrong flags as `cpu`. An unavailable device raises. `cuda` is untested. | cpu |
| `CandidateDedup` | Optional, Python only, opt-in, default `off`: removes duplicate candidates in the coarse levels (`icenine/candidate_dedup.py`, MIGRATION_HISTORY "Task C"). The value is a comma-separated list without spaces: `clique=cpp:K` (before the quick MC keep the best K candidates, by discrete score, of each FZ point), `clique=merged:K[@L]` (the same per connected component of FZ points whose candidates lie within L x diameter of each other, default L = 1; at L = 1 the component of level 0 is almost every FZ point), `pre=T` (before the quick MC, a greedy scan in score order drops a candidate that lies within T x diameter, symmetry-reduced, of a better-scored candidate that was kept), `post=T` (after the quick MC, before the ranking cut, count the kept quarter over candidates that are distinct within T x diameter by the same greedy scan in rank order; passed-over ones are ordered after the distinct ones, also in the FindOptimal hand-off), `quarter=orig|filtered` (the size of the kept quarter is counted from the discrete stage's list, default, or from the list left by the last filter). `diameter` is the level's local-grid radius. Applies at every level; both quick-MC paths (sequential and `BatchedQuickMC`) see the filtered list. The recorder's "discrete" record is unchanged; with the key on the "quick_mc" record gains `dd_n_discrete`, `dd_n_mc`, `dd_n_distinct`. Key off: the default path, bit-identical. | off |
| `SmoothCostRank` | Optional, Python only, opt-in, default `off`: ranks the post-quick-MC candidates of the coarse levels for the quarter cut by a smooth cost (`icenine/smooth_cost.py`, MIGRATION_HISTORY "Task E") instead of the strict cost; the strict cost stays on every candidate and is used everywhere else (quick MC, level 3, FindOptimal, the final cost). The value is the arms joined by `+` (`S`, `E1`, `E2`, `E3`, `E4`) and optional comma-separated `key=value` items, no spaces: `sigma` (E1 kernel in pixels), `omega` (E2 reach in frames past a frame boundary), `tau` (E4 scale in degrees), each a number or `auto` (default: from the level's quick-MC trial radius r = max(diameter / 3, 0.2 deg): sigma = r x detector-0 pixels per radian / 4, omega = max(1, r / frame width), tau = r / 4), `levels` (levels 0 .. levels-1 are ranked, default 3), `ramp`, `ramp_eta`, `ramp_omega` (E3 ramps: 10 px, 1 deg, 0.5 frame). E4 needs lab = Rz(omega) sample, the beam along +x and two detectors. Example: `SmoothCostRank E1+E2+E3`. At a smooth-ranked level the reconstructor's `last_level_best` entry and its "best cost=" print take the best candidate by the strict cost. Key off: the default path, bit-identical (nothing is built; tested against a search recorded before Task E, `tests/data/smooth_keyoff_ref.npz`). | off |
| `RayOverlapRule` | Optional, Python only, opt-in: `off` (default), `R<n>` (isotropic tolerance of n pixels on detector 1) or `A` (anisotropic window per peak). The ray-consistency overlap rule (`icenine/ray_overlap.py`): once a peak's predicted spot is bright on detector 0, the observed ray from the voxel through the bright detector-0 pixels is cast to detector 1 and decides peaks whose predicted detector-1 spot is not bright (cases a-f in the module docstring). **Two costs differ by design:** the rule is used only in the coarse search levels 0-2 (pixel-radius screen, re-cost, quick MC, ranking cut); FindOptimal, VarianceMinimizing, the final cost and the confidence use the strict predicted-spot rule (the C++ rule when `CppOverlapRule` is 1). The coarse levels run the C++ rule and the centroid screen as well. Setting `RayOverlapRule` without `CppOverlapRule 1` gives a `UserWarning` (coarse levels use the ray rule, while level 3, FindOptimal and the final cost use the lenient rule; not studied). The rule covers levels 0-2 only; if `max_local_resolution` < 3 every coarse level uses it. Needs `_rasterize.c` rebuilt (`HAS_RAY_RULE`). Arm A's registration term is `DetectorSpacingDeviation`. Case f (predicted spot off detector 1, observed ray on dark pixels) is always forgiven; the earlier `RayOverlapF1` key (reject) was removed 2026-10-12 because it wrongly rejected about 3% of real peaks and showed no benefit (may revisit; MIGRATION_HISTORY 'Ray-consistency overlap rule'). The observed-ray check runs in C (`_rasterize.ray_flags`; rebuild the extension). | off |

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
uv sync --extra dev --extra riemannian
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
uv sync --extra dev --extra riemannian
uv run python benchmarks/bench_hybrid_optimizer.py --smoke-test --example threevoxels  # 3 voxels, 2 perturbations
uv run python benchmarks/bench_hybrid_optimizer.py --example threevoxels               # full run
```

Outputs: `benchmarks/bench_hybrid_{example}.csv` (per voxel/perturbation/optimizer), plus success-rate, wall-time (Adam vs. hard-eval breakdown), and MC-vs-hybrid scatter plots.

### Local CMA-ES refinement (opt-in)

`CMAOptimizer` (in `icenine/orientation_search.py`) refines an orientation with the `cma` package (a core dependency). The search variable is a rotation vector v (degrees) about the start, `R = exp(v) R0`, with initial step `sigma0_deg` (default 0.2); the result is the lowest-cost orientation evaluated (never worse than the start), as a `CMAResult` (a `SearchCandidate` with `n_evals` and `stop_reason`). `max_evals` (default 1000, start included) is exact. cma's flat-fitness, function-value and stagnation stops are off (the hard cost is quantised); only its x-tolerance, `max_evals` and an optional `max_convergence_cost` (default off) end a run. The cma seed is `seed + 1`; without `seed` it is drawn from the generator, so a run is deterministic given the generator. Settings are those of the Phase B3 `cma_02` method (MIGRATION_HISTORY, "Phase B3 results" and "Phase C results").

**Switch.** `SearchParameters.local_optimizer` is `"mc"` (default, bit-identical to the C++-parity path) or `"cma"`. It is set in code, or in the config with the optional keys `LocalOptimizer cma`, `CMASigma0 <deg>`, `CMAMaxEvals <n>`, `CMAPopSize <n>` (absent keys mean `mc`; existing configs are unchanged; the C++ program rejects a config that contains these keys ("syntax not recognized", Src/ConfigFile.cpp), so use a Python-only config or set the fields in code). The fields are `local_optimizer`, `cma_sigma0_deg`, `cma_max_evals`, `cma_popsize`.

```python
setup = setup_reconstruction(config)           # a config with "LocalOptimizer cma" ...
setup.search_params.local_optimizer = "cma"    # ... or set it in code
setup.search_params.cma_sigma0_deg = 0.2
setup.search_params.cma_max_evals = 1000
BFSReconstruction(setup).reconstruct_sample(...)   # or AdaptiveVoxelReconstructor(setup)
```

Call sites that switch with `"cma"` (all of them are refinements of one start):
- `AdaptiveVoxelReconstructor.refine_from_candidates` (the end of `reconstruct_voxel`, so BFS seed voxels too): for each candidate, FindOptimal's MC **and** the final VarianceMinimizing pass are replaced by one CMA run from that candidate; the candidate loop, the "hit ratio 1.0" early exit, the choice of the lowest cost and the final overlap evaluation are unchanged.
- `AdaptiveVoxelReconstructor.local_optimization` (BFS neighbours, `_fit_from_seed`): the variance-minimizing MC from the inherited orientation is replaced by one CMA run from it; the acceptance test (hit ratio over 0.9 of the best) is unchanged.

Not switched: the coarse discrete search and the quick MC on its candidates (10 steps, 5 restarts; they rank and prune candidates for the next level, they are not a refinement of the final answer), and `BasicVoxelReconstructor` (the serial C++-parity reconstructor, not used by BFS). `local_optimizer="cma"` with `use_hybrid_optimizer` raises `ValueError`. A CMA run uses its whole budget unless x-tolerance stops it (it did not in the Phase C runs), so the cost of a run is `CMAMaxEvals` evaluations per candidate; CMA from sigma0 0.2 deg stays local (B3 far starts, `cma_02`: wrong in 1/11 at 1.5 deg, 4/11 at 2 deg, 10/11 at 3 deg), so `CMASigma0` may need to be larger for neighbours that cross a grain boundary.

### BFS revisit of REFIT voxels, CMA neighbour budget and retry, provenance (opt-in)

Python only (the C++ program rejects these config keys). With `local_optimizer = "mc"` and the options below off, the BFS is bit-identical to before (golden test on a scripted 196-voxel, 4-grain grid recorded from the parent commit). **In `cma` mode the defaults are not the Phase C behaviour**: `cma_neighbor_max_evals` is 250 instead of 1000 and `cma_retry_sigma0_deg` is 1.5 (retry on); set 1000 and 0 to get Phase C's neighbours back.

How C++ treats REFIT voxels, and what this class does:
- In a normal C++ LazyBFS run the server never re-seeds a REFIT voxel (`Pop` skips everything that is not NOT_VISITED). A REFIT voxel is fitted again only because every client keeps its own grid: when another client's expansion (often the voxel's true grain) reaches it while it is still NOT_VISITED there, it is fitted from that client's orientation, and the server's `Push` keeps the higher-confidence result (even over a FITTED one). `LazyBFSClient::Refit` is the restart path only (`RESTART_FIT`, a partial result; every stored voxel is a centre; gated by `PartialResultAcceptanceConfidence`, a key that turns restart on in a C++ config). The Python BFS has one client, so before these options a REFIT voxel was never revisited; that matched C++ with a single client.
- **`bfs_revisit_refit` / `BFSRevisitRefit 1`**: an expansion also enqueues REFIT neighbours it has not tried yet (a per-expansion set stands in for the client grid), at most `bfs_revisit_max` (default 3) local fits per voxel over the whole run (with `bfs_reseed` on, the counters of a cluster are reset before each reseed of it and stay reset if that reseed is rejected, so the cap applies per reseed; `BFSVoxelRecord.n_revisits_before_reset` keeps the discarded counts). The voxel inherits the expanding voxel's orientation and gets the usual local fit (with the cma retry) under the 0.9 relative test; accepted voxels become FITTED and expand; a rejected voxel keeps the better of its old and new fit by confidence (Push rule). Only local fits are made. FITTED voxels are not revisited (C++ would refit them and keep the higher confidence).
- **`bfs_restart_pass` / `BFSRestartPass 1`**: the restart semantics on this run's own result. FITTED voxels that border a REFIT voxel are expansion centres (index order), with revisit enabled and the same cap. There is no centre re-optimization and no full-search fallback. A full search for REFIT voxels is the separate option `bfs_reseed` (next item).
- **`bfs_reseed` / `BFSReseed 1`**: after the seed loop, the connected clusters of REFIT voxels (radius 2 x side length, as the expansion; failed voxels left out) are reseeded one at a time: one full no-start search on the voxel nearest the cluster centroid (ties: lowest index). Accepted by the seed test -> expand with revisit (forced on; the cluster's revisit counters are reset so the cap does not block the new orientation); then the clusters are recomputed and the loop repeats. Finality is per voxel: a voxel whose full search was rejected (ordinary seed or reseed) is failed, never searched again, and ends `unresolved` unless a later revisit fits it; it does not close its cluster, the rest is reseeded on other voxels. No retry, no lowered floor. `BFSReconstruction.reseed_log` has per reseed: cluster size, voxel, accepted, seconds, evaluations; `stats` has `n_reseeds`, `n_reseed_accepted`, `n_reseed_rejected`, `n_reseed_final_voxels` (voxels failed by a rejected reseed), `n_evals_reseed`, `wall_reseed_s` (the full searches are not counted in `n_evals_seed`). Off by default; with it off the output is the same as before (`tests/data/bfs_grid_golden.npz`).
- **`cma_neighbor_max_evals` / `CMANeighborMaxEvals`** (default 250, `cma` only): the budget of BFS neighbour and revisit fits; seeds keep `cma_max_evals`.
- **`cma_retry_sigma0_deg` / `CMARetrySigma0`** (default 1.5, 0 = off, `cma` only): a neighbour that fails the 0.9 test is fitted once more from the same inherited start with the wider step and the same budget; kept only if it passes the test, else the first fit stands.
- **Provenance and stats**: after `reconstruct_sample`, `bfs.records[idx]` is a `BFSVoxelRecord` (`source` in `seed | neighbor | neighbor_retry | revisit | restart | unresolved`, cumulative `n_evals` and `wall_s`, `hit_ratio`, `seed_rejected`, `local_rejected`, `retried`, `retry_accepted`, `n_revisits`, `n_revisits_before_reset`, `reseeded`; `source == "seed"` includes accepted reseeds) and `bfs.stats` the counters (seeds, `n_neighbor_fits`, `n_neighbor_accepted_first`, `n_retry_attempted/accepted_{neighbor,revisit,restart}`, `n_{revisit,restart}_{attempted,accepted,rejected}`, `n_revisit_capped`, `n_restart_candidates`, `n_unresolved`, evaluations and wall seconds split `seed / neighbor / revisit / restart`). "Accepted" is the relative hit-ratio test, not correctness against a truth; score errors against the truth by voxel index. Evaluations and times are cumulative over all attempts on a voxel; the overlap evaluation after a full search is not counted. Recording is one small object per voxel; on the scripted grid the whole 196-voxel run, fake fits included, takes about 0.01 s, so the bookkeeping is small next to any real fit (not measured separately).

```python
setup.search_params.local_optimizer = "cma"
setup.search_params.cma_neighbor_max_evals = 250
setup.search_params.cma_retry_sigma0_deg = 1.5
setup.search_params.bfs_revisit_refit = True    # or "BFSRevisitRefit 1" in a Python-only config
bfs = BFSReconstruction(setup)
bfs.reconstruct_sample(rng=np.random.default_rng(0))
unresolved = [i for i, r in bfs.records.items() if r.source == "unresolved"]
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

`scripts/mc_mechanism/` (finisher/MC study, Phase B1; own header docstrings): `mc_trace.py` (`pilot`, `run`; the T5 cases, clean and realistic, through the unmodified `refine_from_candidates` with `TracedMC`, a `LoggedMC` that records per step the step, `min_ergodic`, steps since the last improvement, accepts, restarts and the distance to the truth; plus the improvement-probability curve P(a proposal of step s lowers the cost) at the finisher result, the truth and the truth rotated by 0.005-0.05 deg), `mc_end_curve.py` (the same curve at the MC stage's own output), `summary.py` (mechanism checks, lock-in bound, curves -> `benchmarks/mc_mechanism/summary.json` and `tables.md`), `mc_helpers.py` (`mc_proposals`, which draws trial rotations exactly like `MCOptimizer.optimize`, and small pure helpers). Caches in `scripts/mc_mechanism/cache/` (gitignored); tests in `tests/test_mc_mechanism.py` and `tests/test_sweep_audit.py`. Nothing in `icenine/` changed.

`scripts/sweep_audit/` (finisher/MC study, Phase B2; no optimizer is run): `audit.py` re-analyses the April 2026 HP sweep (the per-run lines of `benchmarks/hp_sweep_manygrains.log`, linked to voxel positions in the `.mic` file and to `hp_sweep_trajectory_manygrains.csv`; `hp_sweep_threevoxels.csv`) and `bench_hybrid_*.csv`: reproduction of the recorded headline (final error < 1.0 deg), the precision of the "successes", success by distance from the rotation axis, each method at a matched budget, MC stopping, the hybrid benchmark, and a measured-versus-deployed table -> `benchmarks/sweep_audit/summary.json` and `tables.md`; `audit_helpers.py` (log-line parser, file-block splitter, quantiles). Results are in MIGRATION_HISTORY "Phase B1/B2 results".

`scripts/finisher_bench/` (finisher/MC study, Phase B3; own header docstrings): `optimizers.py` (`CountingCost`, which counts cost evaluations exactly, tracks the best orientation, snapshots it at checkpoints and stops a method at its budget; MC with local restarts, a (1+1)-ES with a success-rate step rule, Nelder-Mead and CMA-ES on the rotation vector, VarianceMinimizing in a small box, plain `MCOptimizer`), `bench.py` (`pilot`, `run`; the 300 T5 cases and 1000 perturbation-sweep cases, clean and realistic, every method on the same `VoxelCostFunction`, plus centroid Huber Gauss-Newton on the net's windows, the hybrid Riemannian Adam and the unmodified default finisher), `finisher_summary.py` (win rates, sign / voxel-clustered / McNemar tests, Wilson intervals -> `benchmarks/finisher_bench/summary.json`, `tables.md`), `finisher_timing.py` (single-worker, preflight-gated timing on 40 cases -> `timing.json`, `preflight.json`). Caches in `scripts/finisher_bench/cache/` (gitignored); tests in `tests/test_finisher_bench.py`. Nothing in `icenine/` changed.

`scripts/cma_finisher/` (finisher/MC study, Phase C; own header docstrings): `validate_lib.py` (the library `CMAOptimizer` against B3's stored `cma_02` results on 50 T5 H3 realistic cases, plus `refine_from_candidates` with `local_optimizer="cma"`), `e0_cma.py` (`pilot`, `run`, `timing`, `summary`; no-start `reconstruct_voxel` on the 200-voxel E0 set with `"mc"` and `"cma"` back to back, the `"mc"` run checked bit for bit against the stored E0 seed-0 run; `timing` is single-worker and preflight-gated -> `benchmarks/cma_finisher/`). Caches in `scripts/cma_finisher/cache/` (gitignored); tests in `tests/test_cma_optimizer.py`. Phase C changed `icenine/` (`orientation_search.py`, `reconstructor.py`, `config_file.py`; opt-in, default unchanged) and moved `cma` from the benchmarks extra to the core dependencies.

`scripts/phase_d/` (Phase D data of the finisher/MC study; own header docstrings): `make_neworient.py` (new orientations for the 500-grain sample -> `Examples/Example2.ManyGrains/SimInput/rand_500grains_1mm_neworient_s0.*`), `render_full.py` (full-sample forward simulation, clean and realistic detector-noise images; noise `--noise nodrop` by default, `--noise toy` for the old 10% drop; 10 workers -> `Examples/Example2.ManyGrains/ScatteringData_PhaseD/`, gitignored), `sanity_checks.py`, `frame_split_check.py`, `memmap_check.py`, `bfs_timing_probe.py`, `make_configs.py` (-> `configs/ReconstructPhaseD_{mc,cma,cma_noretry}_{clean,realistic,realistic_q16,realistic_nodrop}.config`, all reading the grid-only `..._grid.mic`), `summary.py`. The opt-in low-memory loader `ExperimentalData.from_binary_memmap` / `write_binary_stack` (`icenine/experimental_data.py`, uint8 `.npy` stack, hard cost only) is documented there. Record in MIGRATION_HISTORY, "Phase D data"; tests in `tests/test_phase_d.py`.

Realistic noise without the random spot drop is the default (Q8; owner D-20261014-I6, 2026-10-14; `docs/todo_q8_noise_nodrop.md`, MIGRATION_HISTORY "Q8"): the old noise erased each spot at random with probability 0.1, whatever its intensity, which is unphysical, so a caller that names no variant now gets no drop (edge jitter, hot pixels and blobs kept). Every named variant keeps its meaning, so the existing images, caches and run records stay valid and are reproduced by name. Phase D: `noise.NoiseParams()` and `add_detector_noise` without parameters are `nodrop` (`p_miss` 0, `noise.DEFAULT_NOISE`); `NoiseParams.named("toy")` is the old 10% drop. `render_full.py` without `--noise` renders `nodrop` into `<tag>/realistic_nodrop/` and writes `results/render_<tag>_nodrop.json`; `--noise toy` writes `<tag>/realistic/` and the tracked `results/render_<tag>.json`. It refuses to write into a `clean/` or `realistic*/` directory that already holds files unless `--overwrite` is given (for an existing tag use `--variants realistic`). The `ReconstructPhaseD_{mc,cma,cma_noretry}_realistic.config` arms still read the old 10%-drop images (`full/realistic`, `full_q16/realistic`); `make_configs.py` also writes the `*_realistic_nodrop.config` twins (not run: D-20261014-I7). E0 per-voxel images: `RealismConfig()` has `p_miss` 0; `RealismConfig.named("noise" | "all")` keep the drop (`orientation_eval.DROP_P_MISS` 0.1) and `"noise_nodrop" | "all_nodrop"` are the same without it. The E0 harness (`findoptimal_robustness/common.py`) renders and runs `clean` and `realistic_nodrop` by default (`C.DEFAULT_VARIANTS`: `e0_run.py` and `e0_cpprule.py run`); `C.VARIANTS` stays `clean, all` as the stored variants that the downstream readers (`analyze_e0.py`, the E2/F1 scripts, `e0_cpprule.py summary`) read by name (`all` = the old realistic images), and `C.ALL_VARIANTS` has all three. `e0_run.py images --variants all realistic_nodrop` builds the cached case images and `e0_cpprule.py run --variant realistic_nodrop --meta-dir <dir>` runs an arm on them. Other study scripts keep their own `VARIANTS = ["clean", "all"]` (the stored results they read; not changed). Q8 scripts: `q8_sanity.py` (three variants at the truth on the Phase D images: cost and rejected real peaks per overlap rule), `q8_e0_truth.py` (the same on the E0 images), `run_q8_step4.sh` (the six E0 arms on three image variants), `q8_e0_summary.py` (wrong counts, evaluations and paired splits), `q8_step7.py` (are the noisy images easier: level-0 candidates and the cost width at the degree scale, clean against no drop, per group of the Step 4 pairing) and `run_q8_offv0.sh` and `q8_lenient_check.py` (the lenient rule with the vertex-0 screen, `e0_cpprule.py --arm offv0`, to explain why the stored lenient counts differ); results in `benchmarks/q8_nodrop/` and `benchmarks/q8_nodrop/step7/`. Tests: `tests/test_phase_d.py` (the defaults are no-drop; `nodrop` frame equals `toy` outside the missed spots; `render_full.py` refuses to overwrite), `tests/test_orientation_eval.py::TestRealismNoDrop`, `tests/test_findoptimal_robustness.py`, `tests/test_q8_step7.py` (the Step 7 helpers).

Phase D pilot (`scripts/phase_d/`): `make_region.py` (the 2,000-voxel disc region -> `benchmarks/phase_d_pilot/region*`), `pilot_bfs.py` (one BFS arm on the region, `BFSRevisitRefit 1`, one process, one thread, memmap images, rng seed 0), `score_bfs.py` (symmetry-reduced misorientation scoring: Wilson and grain-bootstrap intervals, grain-clustered sign-flip test, boundary/interior and near/far axis splits, grains found/lost/fragmented, errors by provenance source), `summarize_pilot.py` (-> `benchmarks/phase_d_pilot/summary.json`, `tables.md`); tests in `tests/test_phase_d_scoring.py`. Results and go/no-go: MIGRATION_HISTORY, "Phase D pilot results".

Ray-consistency overlap rule (`scripts/phase_d/ray_overlap_study.py`, `ray_overlap_summary.py`, `icenine/ray_overlap.py`, tests `tests/test_ray_overlap.py`): opt-in, default off; used only for the coarse levels 0-2 while FindOptimal, VarianceMinimizing and the final cost stay strict (two costs by design). `study.py pilot|run --step 1|2`: step 1 peak-level filtering power and false rejections at the truth, near-truth and random orientations; step 2 the coarse stage levels 0-2 without FindOptimal, per arm. 2026-10-12 additions: `ray_overlap_summary30.py` (30-voxel coarse rerun, env `RAY_COARSE_DIR=coarse30`, survival at 3 / 2 / 1 deg and the pre-MC rank of the closest candidate), `ray_overlap_e2e_summary.py` (step 3: `e0_cpprule.py run --arm rayR2|rayA --msin 0` and `seed_diag.py run --opt ray_rule=R2|A`). Results and the run record: `benchmarks/phase_d_seed_diag/ray_overlap/` (`RUN_LOG.md`, `summary.json`, `tables.md`); MIGRATION_HISTORY, "Ray-consistency overlap rule".

Centroid start point everywhere (2026-10-12, `scripts/phase_d/centroid_step3_summary.py`; env `RAY_VARIANTS=clean`, `RAY_PEAK_DIR`, `RAY_COARSE_DIR`, `RAY_VERTEX0=1` in `ray_overlap_study.py`; `e0_cpprule.py --tag centroid --variant clean`; `seed_diag.py --tag base_centroid`): clean-image re-run of the peak level, the 30-voxel coarse stage, E0 and six full-sample voxels with the centroid screen for every rule; run script and record in `benchmarks/phase_d_seed_diag/centroid/` (`RUN_LOG.md`, `summary.json`, `tables.md`); MIGRATION_HISTORY, "Centroid start point everywhere".

Evaluation reduction, Task A (`scripts/eval_reduction/`, tests `tests/test_cpp_edge_rules.py`, `tests/test_batched_screen.py`, `tests/test_eval_reduction_tools.py`): `parity_inputs.py`, `parity_compare.py`, `parity_summary.py`, `parity_common.py`, `pixelless_check.py` (fixed-orientation C++ vs Python parity per peak, against a scratch C++ driver whose patch is `benchmarks/eval_reduction/A/step0/cpp_parity_patch.py.txt`); `eval_ledger.py` (per-level, per-stage evaluation ledger of a no-start search, `--opt NAME=VALUE`, `summary`); `cpp_ledger_parse.py`, `ledger_attribution.py`, `hit_ratio_check.py` (Python ledger against the C++ per-stage counts, scratch patch `step2/cpp_ledger_patch.py.txt`); `screen_check.py`, `interp_vs_cpp.py`, `screen_timing.py` (preflight-gated), `step3_summary.py`, `e0_interp_summary.py` (`e0_cpprule.py` arms `edge`, `interp`, `lenint`; both `interp` and `lenint` include `CppEdgeRules`), `screen_memory.py` (peak RSS of a no-start search with `BatchedScreen` off and on), `collect_tables.py`. Switches: `CppEdgeRules` and `CppOmegaInterp` opt-in and off, `BatchedScreen` on by default since 2026-10-13 (config table above, `icenine/batched_screen.py`). Results and the run record: `benchmarks/eval_reduction/A/` (`RUN_LOG.md`); MIGRATION_HISTORY, "Task A". The `_rasterize` extension must be rebuilt after pulling (`HAS_VERTEX_LIT`; a stale build with `CppEdgeRules` warns and falls back to the slow path). Task A-exact adds `slow_path_check.py` (Python slow path against the C path with `CppEdgeRules` on and off; `--stub` checks the detector-lit rule on dark images; `--out` saves the run, `benchmarks/eval_reduction/A/a_exact/`).

Evaluation reduction, Task B (`scripts/eval_reduction/`, no library change): `cost_edges.py` (per-peak trace of one cost evaluation, made by wrapping three module functions of `icenine.cost_functions` for one call and calling the C overlap stage once per peak; event classification and attribution of every change of the cost), `cost_edges_run.py` (rays from the truth, 20 directions, eight cost settings x rules; `pilot`, `run --kind full|e0`), `cost_edges_slices.py` (2-D slices through the truth), `cost_edges_summary.py`, `cost_edges_tables.py`, `cost_edges_figures.py`, `cost_edges_diagrams.py`, `peak_rejects.py` (the 3 real peaks the C++ rule rejects at the truth), `interp_levels.py` and `interp_levels_summary.py` (pass-set parity of the Python screen with `CppOmegaInterp` against C++ at levels 0 to 3; scratch C++ patch `benchmarks/eval_reduction/B/levels/cpp_levels_patch.py.txt`). Results: `benchmarks/eval_reduction/B/` (`RUN_LOG.md`, `tables.md`, `figures/`), `docs/findings_cost_edges_2026-10.md`; MIGRATION_HISTORY, "Task B".

Evaluation reduction, Task C (candidate overlap and de-duplication; library: `icenine/candidate_dedup.py`, config key `CandidateDedup`, tests `tests/test_candidate_dedup.py`, `tests/test_candidate_dedup_levels.py`): `scripts/eval_reduction/candidate_probe.py` (a probe around one no-start search that collects, per level, the candidates, their cliques, the evaluations of each quick-MC call and the post-MC orientations from the existing recorder hooks), `candidate_overlap.py` (step 1: `run-full`, `run-e0`, `analyze`; no library change), `candidate_overlap_figure.py` (one level's candidates in orientation space), `candidate_arms.py` (step 3 on the E0 set: `run --arm NAME --spec SPEC`, `summary`; the summary also writes the level-0 trace table `c_e0_level0` and, when the stored runs under `scripts/phase_d/cache/` exist, the baseline-against-stored-runs table `c_e0_reference`), `candidate_summary.py` (ledger and E0 arms side by side, the arm chosen for Task G) and `candidate_arms_figure.py`; `eval_ledger.py --opt candidate_dedup=SPEC` gives the per-stage ledger of an arm. Runners, arm list and results: `benchmarks/eval_reduction/C/` (`arms.txt`, `run_step1.sh`, `run_e0_arms.sh`, `run_ledger_arms.sh`, `run_full30.sh`, `RUN_LOG.md`, `step1/`, `step3/`, `figures/`; the raw per-run arrays are in the gitignored `cache/`). Findings: MIGRATION_HISTORY, "Task C".

Evaluation reduction, Task D (spot pairing across detectors and peak-to-voxel attribution; library: `icenine/spot_pairing.py`, tests `tests/test_spot_pairing.py`): `scripts/phase_d/render_provenance.py` (sibling of `render_full.py`: the same clean render through a recording image; writes per-triangle (voxel, reflection, omega, component) provenance to `<tag>/provenance/provenance.npz`; `--check` compares every frame with `clean/`), `scripts/eval_reduction/pairing_common.py` (geometry, provenance and truth loaders; refuses set-ups where lab = Rz(omega) sample does not hold), `pairing_run.py` (10 presets of the geometric pairing plus 2 binary-centroid runs, about 2 s each), `pairing_summary.py` (step 1 component table and step 2 precision, recall, origin error, rejections), `pairing_attribution.py attach|separate` (step 3: radius and footprint attribution, and truth against wrong orientations with the pixel cost), `pairing_separation_summary.py`, `pairing_figures.py` (`pairing_geometry.png`, `origin_error_map.png`, `attribution.png`, `sharpness.png`). Results: `benchmarks/eval_reduction/D/` (`RUN_LOG.md`, `step1/`, `step2/`, `step3/`, `figures/`, `doc_tables_all.md`; the arrays in `cache/` are gitignored). Findings: `docs/findings_spot_pairing_2026-10.md` and MIGRATION_HISTORY, "Task D".

Evaluation reduction, Task E (smoothing the cost; library: `icenine/smooth_cost.py`, config key `SmoothCostRank`, tests `tests/test_smooth_cost.py`): `scripts/eval_reduction/smooth_landscape.py` (measurement 1: the arms along Task B's rays, `pilot`, `run --kind full|e0 [--set fixed|L0]`; the strict cost is read from Task B's cache), `smooth_landscape_summary.py` (steps per degree, share of the variation in steps, step types from Task B's events, minima, basin radius, the ray figures), `smooth_rank_probe.py` (measurements 2 and 3 through Task C's probe: `run --kind e0|full --rule old|cpp --spec off|SPEC --tag T`; with `--spec off` every arm re-ranks the recorded post-MC lists of levels 0-2), `smooth_rank_summary.py rank|e2e|timing` (basin survival at the quarter cut; paired end-to-end comparison with the broken-voxel trace; the timing tables), `smooth_timing.py` (measurement 4: ms per evaluation, preflight-gated), `smooth_voxel_timing.py` (whole-voxel wall time, key off against on, preflight-gated or labelled contended), `smooth_e4_truth.py` (E4's cost at the truth split into peaks with and without an attached ray) and `smooth_keyoff_fixture.py` (records `tests/data/smooth_keyoff_ref.npz` with the library of a given commit). Results: `benchmarks/eval_reduction/E/` (`RUN_LOG.md`, `landscape/`, `rank/`, `e2e/`, `timing/`, `e4_truth/`, `figures/`; arrays in `cache/` are gitignored). Findings: `docs/findings_cost_smoothing_2026-10.md` and MIGRATION_HISTORY, "Task E".

Evaluation reduction, Task F (searching outward from a start; library: `icenine/outward_search.py`, tests `tests/test_outward_search.py`): `scripts/eval_reduction/outward_bench.py` (`pilot`, `run --kind e0|full`, `timing`: every method from starts at 0.1 to 5 deg in 20 directions; counts strict, screen and smooth evaluations, wall time, the error to the truth against the evaluations so far), `outward_summary.py e0|full` (capture tables with Wilson intervals, the choice rule, the figures, the policy "F4 first, F5 then F4 when F4 is rejected"), `outward_pilot.py schedule|spacing` (the settings of F4 and F5, chosen on four voxels the benchmark does not use). Results: `benchmarks/eval_reduction/F/` (`RUN_LOG.md`, `summary_*.json`, `tables_*.md`, `figures/`, `pilot/`, `timing/`; arrays in `cache/` are gitignored). Findings: `docs/findings_outward_search_2026-10.md` and MIGRATION_HISTORY, "Task F".

Speed parity (`scripts/speed_parity/`, Example2.ThreeVoxels, single thread; the first four scripts changed nothing in `icenine/`, S3 below added two modules): `prof_eval.py` (us per `evaluate` call and aten ops per call for the screen, re-cost and final arms; `--lp "<arm>"` adds a line profile, run with `uv run --with line_profiler`), `stage_split.py` (one call split into the Bragg solve, the torch geometry, the two Python loops and the stage-D C call), `search_split.py` (a no-start `reconstruct_voxel` split into time inside and outside `evaluate`), `translation_check.py` (rigid, linear spot translation under a voxel move; derivation in `docs/spot_translation_derivation.md`). S3 (batched quick MC): `s3_step0.py` (go/no-go measurement), `s3_timing.py` and `run_s3_timing.sh` (whole-search and quick-MC timings, key off against on, CPU and MPS; MPS crossover; preflight-gated), `s3_summary.py` (timing and E0 A/B summary, doc tables), `s3_draw_check.py` (the `uniform` mapping and vectorised trial-rotation checks); results `benchmarks/s3_steps1_5/`, tests `tests/test_batched_cost.py` and `tests/test_batched_quick_mc.py`, A/B arms via `scripts/findoptimal_robustness/e0_run.py` with a tag of their own (`--tag _s3off`, `--batched-quick-mc --tag _s3on`; options `--no-center-screen` and `--no-batched-screen` for the baseline check; a non-default run without a tag is refused, because `cache/e0` is the published Diagnosis cache). S1 (per-call loops out of `evaluate`): `s1_parity.py` (`dump`: every `OverlapInfo` field of 15 cost modes on four voxels to an `.npz`; `compare`: exact `==` of two dumps, including the orientation sets and indices, a shape or dtype mismatch counted as a difference; `compare-e0`: two E0 caches; `batched`: `BatchedCostEvaluator` against a dump; `--fixture` writes the small fixture `tests/data/s1_evaluate_ref.npz`, made on the base commit and read by `tests/test_s1_evaluate.py`, which compares it exactly only on the platform and torch version recorded in the fixture's meta and skips elsewhere), `run_s1_timing.sh` (base commit in a git worktree and S1 interleaved: `prof_eval.py`, `stage_split.py`, `search_split.py`; stops without timing when the preflight is not quiet unless `S1_ALLOW_CONTENDED=1`), `s1_summary.py` (doc tables); results `benchmarks/s1_loops/`. The library side is two helpers in `cost_functions.py` (`frame_of_omegas`, the stage-A lookup, and `FrameTable`, one table of binary frames per cost function, `VoxelCostFunction.frame_table()`; the batched cost shares both) and `BraggConstants` / `scattering_omegas_fixed_beam` in `diffraction_core.py` (the call-independent parts of the Bragg solve, `get_scattering_omegas_torch` unchanged); `evaluate` is bit-identical to before. Findings: MIGRATION_HISTORY, "Python speed parity" (the first four scripts), "S3: batched quick MC with a GPU option" (S3) and "S1: per-call loops out of evaluate" (S1). Open work for the whole port: `docs/ROADMAP.md`.

Seed-cost diagnosis and the variance-stage parity fix (`scripts/phase_d/seed_diag*.py`, `seed_diag_variance_parity*.py`, `tests/test_seed_diag.py`, `tests/test_variance_stage_parity.py`): `MCOptimizer.variance_minimizing_optimize` restarts as C++ does (offsets +-SubregionRadius about the initial orientation); see MIGRATION_HISTORY, "C++ vs Python variance stage".

`MCOptimizer.optimize` (the quick MC of the coarse levels and the FindOptimal MC) is a port of C++ `RandomRestartZeroTemp`: blocks of `int(2 (box/step)^3)` trial steps at a fixed step, each from the current state (`_mc_block`, C++ `ZeroTemperatureOptimization`); a block that does not end strictly below the global best restarts at `delta * initial` (`delta` from U(-r, r)^3, `r = tan(box)/sqrt(48)`) with the step reset, a better one becomes the global best and current state and halves the step; stop on a global cost below `max_convergence_cost`, then on more than `max_restarts` failed blocks in a row. The step budget counts trial steps (each block also evaluates its start); `last_run` holds the stop reason, blocks, accepts and restarts of the last call, and the trajectory records one event per block. Scripted semantics test `tests/test_mc_cpp_parity.py`; C++ trace and evaluation comparison `scripts/phase_d/mc_parity.py` (`benchmarks/phase_d_seed_diag/mc_parity/`); reruns of the earlier MC results in `scripts/mc_mechanism/summary_faithful.py`, `scripts/finisher_bench/mc_faithful_summary.py`, `scripts/cma_finisher/e0_cma.py --tag`, `validate_lib.py --tag`; see MIGRATION_HISTORY, "C++-faithful MC: reruns".

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

- [`docs/findings_finisher_mc_study_2026-10.md`](docs/findings_finisher_mc_study_2026-10.md) —
  findings of the finisher and MC study and the first full-sample BFS pilot (2026-10-07 to
  2026-10-09): cost sensitivity, the two C++ port fixes and what they changed, CMA-ES, the
  full-sample data, the BFS seeding gap, and which earlier claims are superseded.

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
| `scripts/dev/finish_task.sh [--no-tests] [--yes] [--push] [-m msg] [--next task]` | Run the Python suite and/or C++ build if relevant files changed, merge into the parent with `--no-ff` (inside the parent's worktree if it lives elsewhere), delete the task branch, detach and remove any other worktree holding it, optionally push and start the next task. A dirty `CLAUDE.md` is ignored by the cleanliness check and never staged. |
| `icenine_py/scripts/common/stats.py` | `wilson`, `mcnemar_exact`, `paired_discordant`, `win_rate`, `reorder` (id-checked re-indexing), `misorientation_deg_cubic`. Usage: put `scripts/common` on `sys.path`, `import stats`. |
| `icenine_py/scripts/common/doc_tables.py` | `markdown_table(rows, columns, formats)` and `write_tables(path, {name: md})` (marker-delimited blocks). |
| `icenine_py/scripts/common/gzio.py` | Transparent reads of result files stored gzip-compressed: `load_json`, `read_text`, `open_text`, `exists`, `resolve` take the logical name (`x.json`) and open `x.json` or `x.json.gz`, the plain file winning if both exist; `glob_paths(pattern)` matches both; `write_json_like` writes back keeping the stored compression (reproducible gzip bytes). The large Phase D result JSON (`benchmarks/phase_d_seed_diag/runs/*_full_v*.json` of the tags `base_centroid`, `lenientnow`, `cpprule`, `rayR2`, `rayA`, `varfix`, and `level0/py_v15901.json`) is stored as `.json.gz`; every reader in `scripts/phase_d/` goes through this module. To compress another result file: `gzip -9 -n FILE` and make its readers use `gzio`. |
| `scripts/dev/sync_doc_tables.py --doc D.md --tables T.md [--check]` | Copy generated blocks into the matching `<!-- table:NAME -->` markers; `--check` exits 1 on a difference. |
| `scripts/dev/audit_numbers.py --doc D.md --section "Heading" --sources FILES... [--strict]` | List numbers in a doc section that no source value matches (precision-aware; x100, /100 and sibling k/n ratios only for percentages and fractions). `-v` shows the nearest source value and file:key; `--per-file` requires each paragraph or table row to match within one source file. Advisory; prints the chance-match rate per precision bucket. A match is not verification: low-precision decimals match almost anything. |
| `icenine_py/scripts/common/preflight.py`, `scripts/dev/timing_preflight.py [--require] [--json out]` | Record load, power source, thread settings and busy processes before a timing run; `require_quiet` refuses a busy machine. |
| `scripts/dev/job_status.sh` | Read-only snapshot of running jobs, recent logs, `.done` markers and git state. |
| `scripts/dev/verify_handback.sh <parent> <report> [--branch B] [--fixes REVIEW] [--no-tests] [--agent I\|N]` | Mechanical gate on an implementer hand-back before `/review-task` (not a review; at most about 30 lines). Prints `ok`, `FLAG` (the main session decides) or `FAIL` (exit 1): commits since the parent (none = FAIL) and the hashes the report names against the branch; uncommitted changes in the branch's worktree outside `CLAUDE.md` and `.claude/**`; changed files by class (shared Track E/S files, other `icenine_py/icenine/**`, `Src/**`, `pyproject.toml`/`uv.lock`, `tests/data/*golden*`: FLAG; committed `DECISIONS.md`, `CLAUDE.md`, `.claude/**`: FAIL; `ROADMAP.md` hunks outside the caller's own rows and mailbox: FLAG); the eight report fields of `implementer.md` (missing = FAIL) and report length; unlabelled timing lines in Results (heuristic); with `--fixes` every number of the review's "Required fixes" must appear as a fix number in the report; and the `tests/...py[::name]` paths (or bare `test_*.py[::name]` names, taken as `icenine_py/tests/<name>` and flagged when absent) named on the report's Tests field, run in the branch's worktree (`VERIFY_HANDBACK_PYTEST_CMD` overrides `uv run pytest -q`). Own ROADMAP territory comes from `git config --local icenine.agent` (`I`: Track S, row H12, "From Ilium Agent"; `N`: Track E, "From Newt Agent") or `--agent`; with neither it refuses to run (exit 2). The report and `--fixes` paths may be relative to the current directory. |
| `scripts/dev/read_on_fetch.sh [--agent I\|N] [--other BRANCH] [--since REF] [--mark] [--save]` | ROADMAP protocol rule 2 as a read-only report: fetch `origin`, list the other agent's commits and changed files since the last-read point, diff its mailbox ("From Newt Agent" for `I`, "From Ilium Agent" for `N`), list its `D-*-N*` / `D-*-I*` decisions that `shmem query "Owner decision <ID>" --source user` does not find, with ready `shmem add` lines, and print `git merge origin/<other>` (never run). The agent letter comes from `git config --local icenine.agent` or `--agent`; with neither it refuses. The other branch defaults to `feature/eval-reduction` (`I`) or `feature/python-speed-parity` (`N`). Last-read point: `--since`, else `git config --local icenine.<branch>.readOnFetch` (written by `--mark` only), else the merge base with HEAD. `--save` runs the `shmem add` lines (`READ_ON_FETCH_SHMEM_CMD` overrides `shmem`); an unreachable shmem exits 3 and is reported, not read as "not saved". |
| `scripts/dev/checkpoint.sh [--memory-file F] [note]` | Write `.claude/checkpoints/<stamp>.md` (gitignored, under 80 lines): branch, parent, HEAD, commits ahead, status, worktrees, running jobs, plans in `~/.claude/plans/`, latest `.claude/reports/`, plus three sections to fill by hand. With `--memory-file` it keeps one pointer line in that memory file; it warns about running jobs. |
| `scripts/dev/new_todo.sh <slug> "<title>" ["<subtitle>"]` | Scaffold `icenine_py/docs/todo_<slug>.md`. |
| `scripts/dev/log_decision.sh --row R --title T --decision D --where W [--supersedes ID] [--date YYYY-MM-DD] [--agent I\|N] [--recorder NAME] [--trailer LINE] [--no-push] [--no-fetch] [--dry-run]` | Record one owner decision in `docs/DECISIONS.md` on the feature branch above the current task branch, without switching your checkout: it commits in the worktree that holds that branch (refused if its `DECISIONS.md` is dirty), or in a temporary one. Target: follow `branch.<b>.parent` from the current branch to the first branch with no recorded parent, which must be `feature/*` (a `feature/*` branch with no parent is itself the target); a chain ending at `develop`, `master`, `main` or another non-`feature/*` branch, and a detached HEAD, are refused. ID `D-YYYYMMDD-<agent><n>`, n counting from 1 per day and agent, read from the target branch tip. Agent letter from `--agent` or `git config --local icenine.agent` (no default); trailer from `--trailer` or `LOG_DECISION_TRAILER` (required unless `--dry-run`). Refuses when `origin` has commits the target lacks, when `origin` is unreachable (`--no-fetch` skips that check with a warning), on a duplicate title that day, a `--supersedes` that is not an existing heading ID, or a held lock (the lock names the holder's pid; `rmdir` it if that process is gone). Then `shmem add` (`LOG_DECISION_SHMEM_CMD` overrides) and push. Exit codes: 0 recorded, 1 refused or failed with nothing recorded, 2 bad usage, 3 partial failure (the commit was made, but shmem or the push failed: the command to rerun is printed). |
| `scripts/dev/roadmap_row.py --row ID (--status TEXT \| --state WORDS \| --show) [--check] [--doc F] [--agent I\|N]` | Show or edit the last (status) cell of one ROADMAP track-table row, found by its first cell in the `# Track X:` tables (exactly one match, else exit 2). `--status` replaces the cell, `--state` replaces only its first `**bold**` span; rows with and without a closing `\|` keep their form; text with a newline or an unescaped `\|` is refused. Prints the old (`-`) and new (`+`) row; `--check` writes nothing and exits 1 on a change. Editing is limited to the agent's own rows (constant `OWN_ROWS`: Ilium Agent `I` = Track S plus H12; Newt Agent `N` = Track E; any other row exits 2); agent from `--agent` or `git config --local icenine.agent`, no default; `--show` needs neither. Never commits. |
| `.claude/agents/implementer.md`, `.claude/agents/code-reviewer.md` | Slim LLM agents: an implementer (Sonnet) for approved plans, and the reviewer with a "Claims audit" section that runs the scripts above first. Both read `.claude/agents/study-conventions.md` first (standing rules, so briefs can say "follow study-conventions"), write the full report to `.claude/reports/<YYYYMMDD-HHMM>-<task>.md` (gitignored) and return at most about 250 words. |
| `.claude/skills/{handoff,resume,review-task,merge-task}` | `/handoff` runs `checkpoint.sh`, fills Done / Waiting on the owner / Next step and says it is safe to `/clear`; `/resume` reads only the newest checkpoint, its plan and `git status`/`log -5`; `/review-task` launches the opus reviewer, whose required fixes go to its report file; `/merge-task [next]` wraps `finish_task.sh --yes --push [--next]`. |
| `scripts/dev/make_brief.sh --plan P --task T [--section H] [--fixes REVIEW] [--trailer L] [--decisions F]` | Print the implementer brief for an approved plan: "Follow study-conventions", plan and section, the approval (the last DECISIONS.md entry that names the plan file on a path boundary and says "approved", not "not approved" or "unapproved"; refuses when none, or when a later entry names the plan and withdraws it or says `Supersedes:` that ID), the branch commands for the implementer's own worktree (`git checkout -b` and the parent record for a new task; `git checkout` of the freed task branch for fixes), scope (the plan's `**Library changes:**` paragraph, else no change under `icenine_py/icenine/`), the plan's `**Stop rule:**` paragraphs (else "stop and report"), the process limit from `git config --local icenine.maxHeavyProcs` (refuses when unset; a `**Process limit exception:**` paragraph replaces the standard line), the fixed "Refusals" rule (a refused tool or hook action is reported, never redone another way) and "Interruptions" rule (commit as you go), the commit trailer (`--trailer` or `MAKE_BRIEF_TRAILER`; refuses when neither is set) and the report path. A marker paragraph starts at column 0 or after a `- ` bullet, runs to the next blank line, is ignored inside code fences and is quoted whole. `--section` must match exactly one heading (an exact heading wins) and only names the section: the approval and the markers are read from the whole plan. With `--fixes` the brief applies a review report's required fixes; the report may be in another worktree and is printed as an absolute path. |
| `.claude/skills/{decide,implement-plan,apply-fixes}` | `/decide` records an owner decision with `log_decision.sh` (dry run, then real run; `roadmap_row.py` if a status changes); `/implement-plan <plan> <task> [section]` runs `make_brief.sh`, launches the implementer in the background and then `verify_handback.sh`; `/apply-fixes <review report> <plan> <task>` does the same with `--fixes` for a fresh implementer. Each skill only calls its script and shows the output. |
| `icenine_py/scripts/plot_mic.py --ref R.mic [--ref-index I.npy] [--z Z] --mic A.mic ... --labels ... --out F.png [--stats S.json]` | Orientation maps (cubic fundamental-zone Rodrigues vector as RGB, fixed scale) and misorientation maps against the reference (blue 0° to red 62.8°, the cubic maximum) for triangular `.mic` grids. Pairs voxels by index when the order matches, otherwise by nearest centroid on the same grid; unmatched reference voxels are grey. Python 3 port of `MicAnalyzer/Python/xdm_toolkit/xdm_plot.py`. |
| `icenine_py/scripts/mic_maps_report.py SPEC.json --out-dir D` | Build an Artifact-ready page (`index.html`, `maps/*.png`, `stats.json`) of several `plot_mic` figures from a JSON spec (format in its docstring). |
| `.claude/skills/mic-maps` | `/mic-maps`: choose the reference and runs, write the spec, build, check the wrong counts against scored numbers, and publish the page as an Artifact. For layers of 500 voxels or more. |

Tests: `tests/test_dev_tooling.py`, `tests/test_dev_tooling_decisions.py`, `tests/test_dev_tooling_handback.py`, `tests/test_dev_tooling_brief.py`, `tests/test_plot_mic.py`, `tests/test_gzio.py`.

## Citation

S. F. Li and R. M. Suter, "Adaptive reconstruction method for three-dimensional orientation imaging", *Journal of Applied Crystallography*, 2013.
