# C++ template config against the Python E0 parameters (Step 4a)

Python E0: search parameters from `Examples/Example2.ThreeVoxels/ConfigFiles/ReconstructQ8.config`
(`SearchParameters.from_config`, `ReconstructionSetup(config=...)` in `findoptimal_sweep.init_worker`), physics from
`Examples/Example2.ManyGrains/ConfigFiles/Example2.Simulation.config` with `MaxQ` overridden to 8
(`setup_example(max_q=8)`), `min_sin_eta` 0.3, images from the stored pixel keys.
C++: `cpp_template.config` (the Phase D C++ config of the earlier parity runs).

Method: both files normalised (comments dropped, `key=value`, sorted) and diffed with `diff`; the Simulation
config was diffed against the template the same way.

| Key group | Python E0 | C++ template | Match |
|---|---|---|---|
| `LocalOrientationGridRadius` (grid diameter 5 deg), `MinLocalResolution` 0, `MaxLocalResolution` 3 (4 levels, nQMax 5..8) | Q8 config | same | yes |
| `MaxQ` 8 | `setup_example(max_q=8)` | `MaxQ 8` | yes |
| `MaxDiscreteCandidates` 30, `MaxAcceptedCost` 0.9, `MaxConvergenceCost` 0.0001, `MaxMCSteps` 200, `MCRadiusScaleFactor` 0.4, `SuccessiveRestarts` 2 | Q8 config | same | yes |
| `EtaLimit` 86, `BeamEnergy` 64.351, `BeamDirection`, `BeamHeight`, `BeamEnergyWidth`, `BeamDeflectionChiLaue`, `DetectorFilename` (`StandardGeometry.2Det`), `NumDetectors` 2, `DetectorSpacing`, `RotationRangeFilename` (`omega_180_2L.dat`), `FundamentalZoneFilename` (`MyFZ.dat`), `StructureFilename` (`copper.dat`), `SampleSymmetry` Cubic, `MinAmplitudeFraction` | Simulation / Q8 config | same | yes (no diff line; the FZ, omega-range and copper files of the two examples are identical) |
| `MaxInitSideLength`, `MinSideLength` | 0.004 in the Q8 config, but the voxel geometry used by Python comes from the ManyGrains mic (side 0.009375) | 0.009375 | yes in effect: the search does not read the Q8 value, the voxel vertices do |
| Pixel radius of the coarse screen 3, local cost radius 0 | `global_pixel_radius` 3, `_make_local_cost_fn` | `nPixelRadius = 3` (`DiscreteAdaptive.tmpl.cpp:65`), local cost without pixel radius | yes |
| Keep 1/4 per level | `int(n * 0.25)`, at least 1 | `floor(n / 4)`, can be 0 | differs only when n < 4 (not met at levels 0 to 2 in these cases) |
| Peak filter `|sin eta| >= 0.3` | `min_sin_eta` 0.3 in the cost function | none (no key, no code) | **differs; C++ cannot be set**: Python is re-run with `min_sin_eta` 0 (Step 4c) |
| Images | pixel keys (binary) | `.d` frames written from the same keys (intensity constant 1000) | checked in 4b |
| Random number generator | numpy `default_rng(10000 + position)` | `std::mt19937`, fixed default seed 5489 | differs by design; C++ now also run with seeds 1 and 2 (patch 4b) |

Diff output (Q8 config against the template, excluding the placeholders): only `InfileBasename`, `OutStructureBasename`,
`OutfileBasename`, `SampleFilename`, `LazyBFS` (Python-only key absent in the template), `RunReconstruction`
(C++ run switch) and the two side-length keys.

Keys C++ reads but the single-voxel reconstruction does not use (BFS and parameter-optimisation keys:
`Optimization*`, `DetectionLimit*`, `ConsistencyError`, `BraggFilterTolerance`, `MaxDeepeningHitRatio`,
`MinAccelerationThreshold`, `SampleRadius`, `SampleCenter`, `ParamMC*`, `Cooling*`, `Thermalize*`,
`SearchVolumeReductionFactor`): not checked line by line for use; none of them is a search parameter of the
`ReconstructVoxel` path that was read for this task (`DiscreteAdaptive.tmpl.cpp`). Python-only keys with no C++
equivalent: `LocalOptimizer`, `CMA*`, `BFS*`, `CppOverlapRule`, `min_sin_eta`.

Result: no mismatch in the search or physics keys; the two differences are `min_sin_eta` (Python only) and the random
stream.
