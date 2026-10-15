---
title: "TODO: translation screen, correctness and loop ordering (ROADMAP S4)"
subtitle: "Is a voxel's spot set the translated spot set of another voxel, and in which order should the level-0 grid search run over orientations, voxels and detector frames?"
date: "2026-10-14"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Picked by the owner 2026-10-14 (D-20261014-I19); sized down and parallelised (D-20261014-I20); awaiting approval.** Track S,
Ilium Agent. Task branch `feature/python-speed-parity-s4-translation-screen`
(`scripts/dev/start_task.sh`), off `feature/python-speed-parity` after reading Newt Agent's mailbox
(ROADMAP "Between agents"). Implementer: Opus 5.5 at high effort (owner, D-20261014-I19). Agent
review when done (`/review-task`). **No library change**: scripts, a smoke test and docs only.

The first plan (2026-10-13) had a Step 0 go/no-go: a per-voxel table-driven screen saves at most
about 1% of a no-start search, because in the BFS only the 350-504 seed voxels run one. The owner's
pick replaces that gate with two questions, so the 1% ceiling no longer decides the task:

1. **Correctness.** When a voxel is translated, are its diffraction peaks and lit pixels exactly the
   original voxel's pixels moved by the per-peak transform of the translation?
2. **Speed.** Today each voxel runs its own level-0 grid search and then its own optimisation. In the
   new arrangement the grid search runs one orientation (or block) at a time over every voxel
   position, or one detector frame at a time. Which loop ordering gives the best throughput, and what
   does that mean for the optimisation that follows?

# Facts the plan rests on

- **The transform.** For a fixed orientation and peak $p$, the continuous pixel coordinates of any
  voxel point are $\mathbf{u}_p(\mathbf{x}) = A_p\mathbf{x} + \mathbf{c}_p$, with
  $A_p = \Pi Q_p R_z(\omega_p) B$ (`docs/spot_translation_derivation.md`). $\omega_p$, the frame, the
  eta filter and observability do not involve $\mathbf{x}$. So a translation $\Delta\mathbf{x}$ moves
  each spot rigidly by $A_p\Delta\mathbf{x}$, with a different shift per peak and per detector. The
  earlier check reached about $10^{-3}$ px, which is float32 rounding near 1000 px.
- **Pixels.** Lit pixels come from flooring and rasterising. Write the shift as
  $\mathbf{m}_p + \boldsymbol{\varphi}_p$, with $\mathbf{m}_p$ integer and
  $\boldsymbol{\varphi}_p \in [0,1)^2$. Then $\mathcal{R}(T + \Delta\mathbf{u}_p) =
  \mathcal{R}(T + \boldsymbol{\varphi}_p) + \mathbf{m}_p$ exactly. A pure pixel-set translation
  ($\mathcal{R}(T) + \text{round}(\Delta\mathbf{u}_p)$) is exact only when the shift is an integer. For
  a fractional shift, expect differences on the spot boundary. A spot can also leave the detector.
  Voxels come in up and down triangles, so shapes compare only between voxels of the same type. The
  centroid screen (default since 2026-10-12) has no shape.
- **The screen** (`icenine/batched_screen.py`, `_rasterize.c` `pixel_radius_overlap_uint8`). For each
  peak and detector it takes the centroid pixel $\lfloor\mathbf{u}\rfloor$. The detector is "lit" if
  that pixel is on the image, and "bright" if any lit pixel lies in the $(2r+1)^2$ box around it
  ($r = 3$). The C++ qualified-peak rule then turns the per-detector flags into `peak_overlap`, and an
  orientation passes when `peak_overlap > 0` (`run_discrete_search_spaced`). "Bright" equals one
  lookup in the frame image dilated by a $(2r+1)^2$ box, clipped at the image border. This has to be
  checked exactly (Step 1e).
- **What is voxel-independent.** Only level 0 is: 4,886 FZ points times 9 local-grid points gives
  43,974 orientations. The screen uses Q-max $5 + {}$`min_local_resolution`, pixel radius 3. The
  re-cost of the passers (local cost function, triangle mode, Q-max 8) also uses these orientations,
  so its $(f_p, A_p, \mathbf{c}_p)$ are voxel-independent too, and only the set of passers differs per
  voxel. The quick MC, FindOptimal and levels 1 and above propose new orientations per voxel.
- **Test bed.** Phase D clean full sample, here:
  `Examples/Example2.ManyGrains/SimInput/rand_500grains_1mm_neworient_s0_grid.mic`, 24,570 voxels in
  one layer. Images: `ScatteringData_PhaseD/full/clean`, 180 frames of 1 degree times 2 detectors.
  Config: `scripts/phase_d/configs/ReconstructPhaseD_cma_clean.config`. Studies stay on clean images
  (D-20261014-I5). Example2.ThreeVoxels serves for fast checks and the smoke test.
- **Rough sizes (to be measured).** At about 30 peaks per orientation and 2 detectors, a level-0
  screen of every voxel is about $24{,}570 \times 43{,}974 \times 30 \times 2 \approx 6.5\times10^{10}$
  peak tests. At 0.5-2 ns per test that is 30-130 s on one thread. Today it would be about 7.5
  CPU-hours: about 1.1 s per voxel with `BatchedScreen`, a ThreeVoxels figure, to be re-measured on
  Phase D. A full pass/count matrix of 43,974 by 24,570 is about $1.1\times10^9$ entries (1.1 GB as
  int8), so it has to be tiled or written out sparse. Dilated uint8 bitmaps of all 360 images need
  about 1.5 GB at 2048 x 2048, bit-packed about 190 MB.

# Candidate orderings

Every ordering computes the same quantity, `peak_overlap[o, v]` (or only the pass bit), for
orientations $o$ and voxels $v$.

| | Outer loop | Inner loop | Needs a stored table | Locality and parallelism |
|---|---|---|---|---|
| O1 | voxel | all orientations, all peaks | yes, $(f, A, \mathbf{c})$ for all 43,974 orientations | today's order, table-driven (the first plan). Each voxel walks every frame: poor image locality. One voxel per worker |
| O2 | orientation (block of $b$) | all voxels | no: peaks computed per block and discarded | the owner's "one orientation at a time over every voxel position". Each peak touches one frame region V times. Parallel over voxels (vector, GPU) and over blocks (threads) |
| O3 | frame and detector | all $(o, p)$ landing in that frame, times voxels | yes, sorted by frame | the owner's "detector positions" read as frames. Each image is loaded once. The accumulator over $(o, v)$ is the cost, so it is tiled over voxels |
| O4 | lit (dilated) pixel of a frame | $(o, p)$ in that frame: back-project the pixel to the voxels it covers and vote | yes, sorted by frame | "detector positions" read as Hough voting: a scatter, with atomics on a GPU. Wins only if lit pixels in a layer's footprint are much rarer than voxels |
| O5 | orientation block times voxel tile (tiled O2/O3) | the other | O2 form: no | the usual answer for gather-heavy kernels: block and tile sized to L1/L2. Sweep both |

Each ordering gets a cost model before it is built: operations, image bytes touched, working set and
the parallel axis. The measurements in Step 3 confirm or correct the model.

# Steps

**Step 0: baselines and the reference set.** Before any timing, run `scripts/dev/timing_preflight.py`
and record the load. Label every timing "quiet" or "contended" (D-20261014-I16). Build the level-0
global cost function and orientation set exactly as `reconstruct_voxel` does for the Phase D clean
config: same flags, `n_q_max`, pixel radius and centroid screen, read from the setup, not copied by
hand. Then measure:

- `BatchedScreen` per voxel, run as the reconstructor runs it (chunk 32), on 10 seeded random Phase D
  voxels: wall time split into torch geometry and the per-orientation stage-D loop, peaks per
  orientation after the eta filter and frame lookup, and **screen passers per voxel**. The passers set
  the re-cost load in Step 4.
- Image size, lit-pixel density per frame, raw and dilated by $r$, and how many frames occur.
- The re-cost time per voxel (`evaluate` of each passer with the local cost function), and the time of
  the rest of a no-start search, on 3 voxels. Reuse the S1/S3 numbers where they apply and say which.
- **Reference:** `peak_overlap` of all 43,974 orientations from `BatchedScreen`, for the 3 ThreeVoxels
  voxels and 20 seeded random Phase D voxels (both triangle types), saved under
  `benchmarks/s4_translation/`. Every later step is compared with it. It is an equality reference,
  not a timing, so it runs in up to 6 worker processes in the background.

**Step 1: correctness of translation (owner question 1).** Script
`scripts/speed_parity/s4/translation_correctness.py`. Orientations: 100 random over the FZ, plus 20
within 1 degree of the truth for 3 voxels. Position pairs: same triangle type, displacements from one
voxel step up to the sample diameter (1 mm).

- **1a, continuous.** Per peak and detector, compare the projected vertices and centroid computed
  directly (a float64 re-implementation of stages B and C, and the existing float32 path) with
  $\mathbf{u}_p(\mathbf{x}) + A_p\Delta\mathbf{x}$. Expected: at most about $10^{-9}$ px in float64
  and about $10^{-3}$ px in float32. Also check that the peak sets (observable reflection, omega
  solution, frame, eta filter) are identical, which is expected exactly. Report the largest deviation
  and where it occurs.
- **1b, lit pixels from the forward model.** Use the simulator's rasteriser, the one that renders the
  Phase D images and is validated against C++ (find it from `scripts/phase_d/render_full.py` and
  `orientation_eval.lit_pixel_set`). Rasterise one voxel at $\mathbf{x}$ and at $\mathbf{x}'$, then per
  peak and detector compare $S(\mathbf{x}')$ with:
  (i) $S(\mathbf{x}) + \mathbf{m}_p$ for a displacement built so that one chosen peak's shift is an
  integer (solve the in-layer $2\times2$ system). Expected exactly equal for that peak.
  (ii) $\mathcal{R}(T + \boldsymbol{\varphi}_p) + \mathbf{m}_p$, re-rasterising only the sub-pixel
  phase, for every peak and general displacements. Expected exactly equal.
  (iii) $S(\mathbf{x}) + \text{round}(A_p\Delta\mathbf{x})$, a pure pixel-set translation. Expected to
  differ on the boundary. Report the symmetric-difference counts per spot, as a distribution relative
  to the spot size, and confirm that every differing pixel touches the spot boundary.
  List spots that leave or enter the detector separately; they are expected, not failures.
- **1c, the cost-level effect of (iii).** On the Phase D images, for the 20 reference voxels and their
  re-cost passers, compare the final cost (triangle mode, `pixel_overlap`, `peak_overlap`, quality) of
  exact per-voxel rasterisation with the pure pixel-set translation. This decides whether "rasterise
  each spot shape once per (orientation, peak, triangle type), translate per voxel" is usable for the
  re-cost or only as a pre-filter.
- **1d, optional C++ cross-check.** If a Release build of `IceNine` takes minutes, simulate a
  single-voxel `.mic` at two positions (mode `s`) and compare spot pixel sets as in 1b(ii). Otherwise
  skip it and say so.
- **1e, the dilated-bitmap screen test.** For every pixel of 4 frames (both detectors, including the
  image border), the dilated-bitmap lookup and the centre bounds test must equal
  `pixel_radius_overlap_uint8`'s `(in_bounds, bright)` exactly.

**Stop rule:** if 1a in float64, 1b(i), 1b(ii) or 1e fails outside the listed exceptions, stop and
report. The premise of S4 is then wrong.

**Step 2: the table-driven screen and its equality.** Script `scripts/speed_parity/s4/screen_kernel.py`
(names may change):

- **Peaks.** Per orientation block, run the torch operations of `ScreenEvaluator._counts` up to the
  frame lookup unchanged, in float32, so the peak sets match. Batch shape can still flip a float32 last
  bit (MIGRATION_HISTORY "Task A"); count any frame that changes. Then compute $A_p$ and
  $\mathbf{c}_p$ in float64 from those directions, together with the parallel flag and the sign of the
  ray parameter (constant per peak, since the whole sample lies on the beam side of the detector).
- **Kernel.** $\lfloor A_p\mathbf{x} + \mathbf{c}_p\rfloor$ with `_pixel_coordinate`'s rule (the C++
  negative-coordinate rule, and the `CppEdgeRules` form if the config turns it on), the hit flag, the
  lit and bright flags from the dilated bitmap, the qualified-peak rule with the config's
  `cpp_overlap_rule`, then `peak_overlap[o, v]`.
- **Equality** with the Step 0 reference: per orientation, `peak_overlap` and pass. List every
  differing orientation with its voxel, the coordinate that flipped, and that coordinate's float64
  distance to the nearest pixel boundary. Do the same for a float32 table to measure its flip rate.

**Step 3: the loop-ordering study (owner question 2).** Script `scripts/speed_parity/s4/orderings.py`.
Build O1, O2, O3 and O5 on Step 2's kernel semantics. Check each against Step 2's output on the
reference voxels before timing it.

- **Sample size.** Time every ordering on one contiguous patch of 2,000 voxels, not the whole layer.
  Check that time grows linearly with the voxel count (500, 1,000 and 2,000 voxels), then project to
  24,570. Only the best ordering gets one measured full-layer run.

- **Implementation.** Vectorised numpy and torch first: numpy fixes the loop order to its vector axis,
  which is fine for O2. Loop order matters for O1, O3 and O5, so where the vectorised form hides it,
  write a scratch C kernel in `scripts/speed_parity/s4/`, compiled with the system `cc` and loaded
  through `ctypes`. No change to `icenine/_rasterize.c` and no new dependency.
- **O4** gets its cost model only, unless that model predicts at least 2x over the best of O2, O3 and
  O5. In that case build it too.
- **Sweeps.** For O2, a block size $b$ of one FZ point (9 orientations), 32 points, and larger. For O3
  and O5, voxel tile sizes. Bitmaps as uint8 and bit-packed, keeping only the frames that occur. Pass
  only (stop at the first qualified peak) against full `peak_overlap` counts, keeping the counts for
  the equality checks.
- **Threads and device.** CPU at 1, 2, 4 and up to the core count. Threads go over orientation blocks
  for O2, and over frames or voxel tiles for O3. Bitmaps are shared and read-only; at most 6 processes
  if processes are used on this 16 GB machine. Torch MPS for O2: MPS has no float64, so use
  $\mathbf{u} = \mathbf{u}_{\text{ref}} + A_p(\mathbf{x} - \mathbf{x}_{\text{ref}})$ with
  $\mathbf{u}_{\text{ref}}$ computed in float64 on the CPU, and count the flips.
- **Metrics per ordering.** Peak tests per second; (orientation, voxel) pairs per second; projected
  full-sample level-0 wall time for 24,570 voxels times 43,974 orientations; peak memory; table build time (O1, O3); the cost model's prediction next
  to each measurement. Also the passer output: count and memory of the sparse (voxel, orientation)
  list for the full sample.

**Running the steps in parallel (D-20261014-I20).** Correctness work and timing work are kept apart:

- Run concurrently, in the background: the Step 0 reference (up to 6 worker processes), all of Step 1,
  and the Step 2 equality runs. Write Step 2 and Step 3 code while they run.
- Run alone, one at a time, with nothing else running on the machine: every timing (the Step 0
  timings, Step 3 and the Step 4 re-cost timing). A timing taken next to another job is contended and
  the comparison between orderings loses its meaning.
- One implementer, not two. A second implementer would need the same Step 0 reference, and its jobs
  would contend with the timings; background jobs in one session give the same overlap.
- Memory: at most 6 heavy processes at once on this 16 GB machine; share the dilated bitmaps
  read-only (threads, or a memory map across processes).

**Step 4: what follows the screen (the "optimisation follows" half).** One measurement and an
analysis; no reconstruction run in this task.

- **Measure.** Re-cost the passers in orientation-major order: for each orientation, all voxels that
  passed. Take the vertices from $(A_p, \mathbf{c}_p)$ at Q-max 8 and the vertex offsets
  $A_p\boldsymbol{\delta}_i$, and rasterise exactly per voxel, through the existing stage-D C call or
  in a batch. Compare with today's per-voxel `evaluate` of each passer, for time and for equality of
  cost and confidence on the reference voxels.
- **Project the whole-run time** from the measured screen and re-cost and the existing per-voxel
  quick MC and level numbers (S1 and S3) for:
  (a) today's BFS: seeds run a no-start search, the others start from their neighbours;
  (b) one global screen and re-cost for every voxel, then per-voxel refinement from each candidate
  list;
  (c) a hybrid: BFS as today, but every voxel's level-0 candidates come from one global pass. Seeds,
  and voxels whose inherited orientation fails, start from that list instead of a fresh no-start
  search.
  Give each stage's share after S4. State plainly if the refinement after the screen dominates, which
  is likely, since it is per voxel and gains nothing from translation.
- **Recommend:** the ordering (screen and re-cost), CPU or MPS, and a sketch of a library path (a new
  module behind an opt-in key, no change to existing paths, the files it would touch). This is a
  proposal for the owner, not built here.

**Step 5: docs and test.** A MIGRATION_HISTORY section "S4: translation correctness and loop
ordering". Results and a summary in `benchmarks/s4_translation/`. The ROADMAP S4 row and this plan's
status. A smoke test `tests/test_s4_translation.py` on ThreeVoxels, under 15 s: the float64 affine law,
the dilated bitmap equal to `pixel_radius_overlap`, and the table screen equal to `BatchedScreen` on a
few FZ points. The README developer-tooling list gets the new script directory.

# Deliverables

- **Question 1, answered with numbers:** the largest float64 and float32 deviation of the affine law;
  integer-shift and phase-re-rasterised pixel-set equality (expected exact); the error of a pure
  pixel-set translation and its effect on the final cost; dilated-bitmap equality; the listed
  exceptions (detector edge, float boundary).
- **Question 2, answered with numbers:** a table of orderings with measured throughput, full-sample
  level-0 time, memory and scaling with threads (and MPS). The best ordering, and why, from the cost
  model and the measurements. The projected whole-run effect for (a), (b) and (c).
- **Equality:** the new screen's pass sets equal `BatchedScreen`'s on the reference voxels, apart from
  listed float-boundary flips, each with its distance to the pixel boundary.

# Risks

- **The refinement dominates.** Once the screen costs seconds for the whole sample, the re-cost of the
  passers and the per-voxel optimisation set the run time. Step 4 measures the re-cost and projects
  the rest rather than assuming a gain.
- **Float flips** at pixel boundaries: the float64 table against today's float32 path, and float32 on
  MPS. These are listed, not hidden.
- **Memory:** the $(o, v)$ accumulator (tile it), dilated bitmaps (bit-pack them), the O1 and O3 table
  (about 0.2 GB in float64; measure it).
- **Contended machine:** the preflight has not been quiet lately (load 2.6-2.8). Timings use
  interleaved repeats and carry the label. Ratios between orderings measured in the same run are more
  robust than absolute times.
- **One layer:** Phase D is one layer. Layer-to-layer effects are out of scope.

# Scope

The level-0 screen and its re-cost; the default level-0 rules of the Phase D clean config; centroid
screen. Out of scope: `CppOmegaInterp`, though the interpolated omegas are also orientation-only and
could be added later; the ray rule (`BatchedScreen` does not support it either); noisy images; any
change to `icenine/` or to a default.

# Relation to other work

- ROADMAP Track S: S5 (a BFS neighbour's first evaluation as a translation) uses the same
  coefficients and would follow naturally from Step 4(c). S2 (fused C evaluation) and S3 (batched
  quick MC) are the levers for the refinement that Step 4 projects. S6 (full-sample timing) is
  partly covered, because the Phase D images are now on this machine.
- `docs/spot_translation_derivation.md` (the affine law), MIGRATION_HISTORY "Task A" (`BatchedScreen`
  equality and timing), "S1" and "S3" (per-evaluation and quick MC costs).
- Track E: unaffected. No shared library file is touched; the mailbox gets a heads-up when the task
  merges.
