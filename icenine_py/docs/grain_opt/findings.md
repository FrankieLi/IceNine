# Full-grain optimization: findings as we go

Branch `feature/grain-opt` (from `feature/eval-reduction` d11b1a63). Plan: the "Full-grain
optimization plan" doc (claude.ai artifact B6TkLDJRwrMPsC9JkDtgLf). Module
`icenine/grain_fit.py`, tests `tests/test_grain_fit.py`, scripts `scripts/grain_opt/`, results
`benchmarks/grain_opt/`. Sample `grad` (497 grains, linear gradients 0 to 1 deg), Q-max 8.

## Data

The `grad` image stack was not on disk (the worktrees STATE.md lists are gone). It was rendered
again with the recorded command (`render_full.py --tag grad --variants clean --mic
gradient_truth.mic`, 4 workers, 58 s) into this worktree's git-ignored
`Examples/Example2.ManyGrains/ScatteringData_PhaseD/grad/`: 9,444,080 lit pixels and 104,912
spots, equal to the stored record (`benchmarks/autoresearch/ar3/gradient/grad_stack.json`).
`scripts/grain_opt/build_store.py` wrote the binary stack `stacks/grad_clean.npy` (9,444,080
non-zero, as recorded) and a sparse clean-intensity store `grad/intensity_clean.npz`.

## Shmem lookup (step 2)

Memory holds no other statistical method for matching intensity distributions. The main objective
stays the per-reflection Sinkhorn distance of `docs/todo_intensity_and_distribution_losses.md`
(section 2a), with geometric coverage as the prediction.

## Work item (1): grouping, regions, renderer

**Renderer.** `SpotModel` is `SmoothCost.predict` vectorised over voxels (each voxel its own
orientation and triangle; test: equal rows, frames and vertex pixels to 1e-9 on ThreeVoxels).
Each projected triangle is cut into n_sub^2 equal sub-triangles; their centroids carry the
projected area. At the ground-truth orientations of 20 `grad` grains (4 per level, sizes at quantiles 0.1
to 0.95), the share of coverage mass that falls on a lit pixel of the stock render, in the
predicted frame (pixel = floor of the coordinate), is a median 0.996 (min 0.952, 40 grain x
detector values). Mass that lands lit only in an adjacent frame is at most 0.0046. Test on
ThreeVoxels: above 0.95.

**Grouping.** Flood fill with neighbour misorientation below 2 deg: the ground-truth map gives exactly the
497 true grains (24,570 of 24,570 voxels in their true grain). The painted map gives 497 groups
holding 24,400 voxels in their group's majority true grain; 115 true grains are split (wrong
painted voxels form their own small groups).

**Region.** Per reflection, omega branch and detector: the start spot's bounding box padded by 6
pixels and one frame. Three filters were needed (found on 5 grains, `probe_objective.py`):

1. Observed connected components (8-connected in a frame; frames are not joined, which merged
   other grains' spots in the next frame) that come within 3 pixels and one frame of the start
   footprint. Other grains' spots in the window are dropped.
2. A ratio filter: observed lit pixels / pixels touched by the start prediction in [0.8, 1.6].
   Long thin streaks overlapping another grain's spot sit near 2. On the 20 grains this keeps
   44 to 99 of 105 to 221 spots (median about 85).
3. Unbalanced Sinkhorn (reach 3 pixels) for the remaining contamination. At the ground truth the mean
   divergence fell from 14 to 31 px^2 (balanced, before filters) to 0.4 to 0.7 px^2.

With these, on 5 probe grains the Sinkhorn divergence of a whole-grain rotation away from the
ground truth grows smoothly: about +10% at 0.05 deg, +50% at 0.1 deg, 3x at 0.2 deg, 10x at 0.5 deg
(least squares: similar shape). One evaluation takes 0.1 to 0.4 s (one thread).

## GO1: objective bias (ideas I1, I2, I17, I21)

Question: where does each objective's minimum sit against the ground truth? On the two level-0
grains 381 (21 voxels, 90 spots) and 302 (40 voxels, 93 spots), each grain has one ground-truth
orientation. CMA over the grain rotation only (3 numbers), started at the ground truth, 120
evaluations, sigma 0.03 deg (`probe_bias.py`; rows in `benchmarks/grain_opt/probe/`). Rule
(written before the runs): an objective qualifies if its optimum is below 0.01 deg from the
ground truth on both grains.

Reading the table: `offset` is the optimum's distance to the ground truth (deg); `drop` is how
much lower the objective is there than at the ground truth; `strict` is the mean strict cost at
the optimum (ground truth: 0.0806 on 381, 0.0522 on 302).

| objective | offset 381 | offset 302 | drop 381 | drop 302 | strict 381 | strict 302 |
|---|---|---|---|---|---|---|
| Sinkhorn, reach 3 px, blur 1 | 0.0330 | 0.0407 | 0.71% | 0.21% | 0.1048 | 0.0700 |
| Sinkhorn, balanced, blur 1 | 0.0947 | 0.0127 | 4.37% | 0.10% | 0.1651 | 0.0792 |
| Sinkhorn, reach 3, blur 0.5 | 0.0308 | 0.0409 | 0.55% | 0.18% | 0.1041 | 0.0699 |
| + pixel-box dilation, blur 1 | 0.0337 | 0.0408 | 0.77% | 0.22% | 0.1049 | 0.0701 |
| + pixel-box dilation, blur 0.5 | 0.0337 | 0.0407 | 0.67% | 0.20% | 0.1044 | 0.0700 |
| Sinkhorn, intensity-weighted | 0.0074 | 0.0263 | 0.09% | 0.03% | 0.1015 | 0.0692 |
| least squares, per-spot scale | 0.0117 | 0.0423 | 0.64% | 3.79% | 0.1025 | 0.0731 |
| strict cost (grain mean) | 0.0060 | 0.0082 | (none) | (none) | 0.0815 | 0.0528 |

- Only the strict cost qualifies. Every binary-mass Sinkhorn variant ends 0.031 to 0.041 deg from
  the ground truth with an objective under 1% lower: the bottom is flat at that scale, not
  shifted. Blur and the pixel-box dilation change nothing; balanced transport is worse.
- Intensity weighting sharpens the bottom on grain 381 (0.0074 deg) but not on 302 (0.026 deg).
  Least squares is in between.
- Moving 0.03 to 0.04 deg raises the strict cost by 30 to 35%; the strict cost resolves what the
  smooth objectives do not. The likely reason: the observed spot is a set of whole lit pixels and
  a voxel's spot moves well under a pixel for 0.03 deg; a transport distance on pixel centres
  averages that away, while the overlap count changes as edge pixels flip. Not tested (I4).
- Consequence for GO2: Sinkhorn can only be a basin finder; the finish uses the strict cost.

## Smoke test before GO1 (superseded)

A two-grain stage-A test on Sinkhorn (60 CMA evaluations) lowered the objective while moving
voxels away from the ground truth (median 0.016 to 0.050 deg on grain 381, 0.032 to 0.081 deg on
306); the strict-cost gate rejected both. GO1 above explains it.

## GO2: optimizer Pareto on G20 (ideas I5, I6, I9, I10, I28, I29)

Question: how many objective evaluations per voxel buy how much accuracy, and which chain of
stages reaches the strict-cost optimum cheapest? Benchmark G20 (20 `grad` grains, painted start,
start-map flood-fill groups, 1,082 voxels; 853 of them have an AR3 C record). Pipelines of
`run_pipeline.py`: A = grain-field CMA on an objective, V = per-voxel CMA on the strict cost
(sigma 0.02 deg). One objective evaluation scores every voxel of the grain once (grain fits) or
one voxel (per-voxel fits); chained stages add up. Tables: `benchmarks/grain_opt/go2/`.

Reading the table: on the 853 voxels with a C record; `eval/voxel` objective evaluations per
voxel; `ms/voxel` wall time of the fit stages per voxel (contended: Newt's runs share the
machine); `> 0.1 deg` pooled count with Wilson 95% interval; `median` and `p95` in degrees.

| arm | eval/voxel | ms/voxel | > 0.1 deg | median | p95 |
|---|---|---|---|---|---|
| painted start | 0 | 0 | 84 (8.0 to 12.0%) | 0.0349 | 0.135 |
| AR3 per-voxel CMA, C (sigma 0.05) | 61 | 6.5 | 18 (1.3 to 3.3%) | 0.0181 | 0.055 |
| per-voxel CMA, sigma 0.02, 20 | 22 | 3.4 | 45 (4.0 to 7.0%) | 0.0235 | 0.102 |
| per-voxel CMA, sigma 0.02, 61 | 62 | 6.6 | 34 (2.9 to 5.5%) | 0.0165 | 0.083 |
| grain field, strict, 50 | 63 | 1.8 | 20 (1.5 to 3.6%) | 0.0148 | 0.043 |
| grain field, strict, 100 | 113 | 3.1 | 3 (0.1 to 1.0%) | 0.0146 | 0.035 |
| grain field, strict, 200 | 213 | 5.2 | 3 (0.1 to 1.0%) | 0.0111 | 0.033 |
| grain field, strict, 400 | 413 | 11.1 | 3 (0.1 to 1.0%) | 0.0062 | 0.022 |
| quadratic field, strict, 200 | 207 | 5.7 | 20 (1.5 to 3.6%) | 0.0170 | 0.049 |
| field 100 + voxel 10 | 128 | 5.6 | 3 | 0.0151 | 0.034 |
| field 100 + voxel 30 | 149 | 7.2 | 3 | 0.0130 | 0.032 |
| field 100 + voxel 61 | 175 | 9.4 | 3 | 0.0114 | 0.030 |
| Sinkhorn 150, field 100, voxel 20 | 298 | 2,287 | 3 | 0.0145 | 0.033 |

On the 10-grain subset (494 C voxels): strict field 100 + voxel 20: 1 above 0.1 deg, median
0.0153 (5.6 ms/voxel); least squares 150 first: median 0.0136 (112 ms/voxel); intensity-weighted
Sinkhorn 150 first: median 0.0155 (1,391 ms/voxel).

- The grain field on the strict cost is the frontier: at 100 evaluations per voxel it beats C on
  voxels above 0.1 deg (3 against 18, intervals apart) and on p95 (0.035 against 0.055), at half
  of C's wall time. The 3 left are the 3 wrong voxels (above 1 deg) of the paint; no local fit
  reaches them. More evaluations lower the median (0.0146 at 100, 0.0062 at 400).
- A per-voxel polish after the field lowers the median a little (0.0151 to 0.0114 for 10 to 61
  more per voxel); 300 more grain evaluations do better (0.0062).
- The quadratic field (5 more coefficients per component) is worse at 200 than the linear field at
  100: CMA spends its budget on directions `grad` does not need (its gradients are linear).
- Per-voxel CMA alone is worse than the field at every budget here: each voxel's strict cost is
  noisy (few pixels per voxel), and the field pools all of a grain's spots.
- Smooth objectives as basin finders add nothing from this start: the Sinkhorn stage was kept on
  3 of 20 grains (gate), the end result equals strict alone, at about 700 times the wall time;
  least squares was kept on 7 of 10 and lowers the median slightly (0.0136 against 0.0153) at 20
  times the wall time. The painted start is already inside the strict basin; GO4 tests starts
  that are not.

## GO4: robustness across kinds of orientation variation (ideas I13, I6)

Question: does the GO2 frontier hold when the voxel-to-voxel variation is random (`rand`) or
spatially correlated (`corr3`, length 3 voxel sides), 0.5 to 5 deg RMS? Start: each voxel at its
grain's orientation before deformation (a perfect paint that knows nothing of the variation);
true grain ids; 20 grains per sample (5 per level), 1,043 voxels. Tables `benchmarks/grain_opt/go4/`.

Reading the table: share of voxels more than 0.1 deg from the ground truth (pooled over the 5
grains of the level) and the median error (deg), at the start and after each arm. Arms: field =
`Astrict/100`; + fine = field then per-voxel CMA 61 at step 0.02 deg; + coarse = the same at
step 0.3 deg; coarse only = per-voxel CMA 61 at 0.3 deg; LS first = `Als/150` then field + coarse.

| sample, level | start | field | + fine | + coarse | coarse only | LS first |
|---|---|---|---|---|---|---|
| corr3, 0.5 | 99%, 0.414 | 94%, 0.341 | 69%, 0.224 | 32%, 0.075 | 36%, 0.080 | 32%, 0.071 |
| corr3, 1 | 100%, 0.878 | 98%, 0.652 | 91%, 0.554 | 59%, 0.154 | 71%, 0.183 | 53%, 0.108 |
| corr3, 2 | 100%, 1.773 | 100%, 1.686 | 99%, 1.706 | 93%, 1.054 | 95%, 0.959 | 93%, 0.929 |
| corr3, 5 | 100%, 4.347 | 100%, 4.409 | 100%, 4.374 | 100%, 4.261 | 100%, 4.207 | 99%, 4.191 |
| rand, 0.5 | 97%, 0.443 | 97%, 0.460 | 90%, 0.452 | 71%, 0.314 | 72%, 0.363 | 68%, 0.334 |
| rand, 1 | 100%, 0.869 | 99%, 0.867 | 95%, 0.863 | 85%, 0.854 | 85%, 0.853 | 86%, 0.815 |
| rand, 2 | 100%, 1.818 | 99%, 1.860 | 98%, 1.877 | 95%, 1.862 | 95%, 1.812 | 97%, 1.907 |
| rand, 5 | 100%, 4.478 | 100%, 4.483 | 100%, 4.475 | 99%, 4.535 | 100%, 4.476 | 98%, 4.684 |

- From a start that ignores the variation, every arm recovers only voxels within about 0.5 deg of
  their start, the width of the strict cost's basin. On `corr3` at 0.5 deg the share above 0.1
  deg falls from 99% to 32%; at 1 deg to 53% (LS first); at 2 and 5 deg almost nothing is
  recovered. On `rand` nothing above 0.5 deg is recovered.
- The linear field alone barely moves either sample (it cannot represent the variation), so the
  GO2 winner does not carry over: robustness needs per-voxel freedom with a wider reach.
- A coarse per-voxel step (0.3 deg) beats the fine one (0.02 deg) everywhere here; least squares
  as a basin finder helps a little on `corr3` (median 0.154 to 0.108 at 1 deg).
- The start is pessimistic: a real paint seeds many grain entries in deformed grains (AR3 step 3:
  24 per level-1 grain), so voxels start nearer their own orientation. GO4b tests strategies that
  do not need that: continuation from fitted neighbours and a wider per-voxel search (I33, I34),
  and a real paint start (I35).

## GO4b: reaching far starts (ideas I33, I34)

Question: which per-voxel strategy recovers deformed grains whose voxels start 0.5 to 5 deg from
their ground truth? Same samples, start, membership and benchmarks as GO4. Stages: `L/n` per-voxel
CMA step ladder 1.0 / 0.3 / 0.05 deg with n/3 evaluations per leg; `N/30/0.3` neighbour
continuation (breadth-first from the best 10%, candidates = own and fitted neighbours'
orientations, then CMA 30 at 0.3 deg; two sweeps). Tables `benchmarks/grain_opt/go4b/`.

Reading the table: share of voxels more than 0.1 deg from the ground truth and median error
(deg) per level; `eval/voxel` objective evaluations per voxel; `ms/voxel` wall time (contended).

| arm | eval/voxel | ms/voxel | 0.5 deg | 1 deg | 2 deg | 5 deg |
|---|---|---|---|---|---|---|
| corr3, GO4 best (LS + field + coarse) | 340 | 57 | 32%, 0.071 | 53%, 0.108 | 93%, 0.929 | 99%, 4.191 |
| corr3, L/61 | 70 | 9 | 49%, 0.097 | 71%, 0.219 | 86%, 0.598 | 97%, 4.136 |
| corr3, L/183 | 196 | 18 | 9%, 0.021 | 11%, 0.018 | 47%, 0.080 | 91%, 3.498 |
| corr3, N/30 | 87 | 8 | 36%, 0.077 | 53%, 0.108 | 76%, 0.433 | 98%, 4.132 |
| corr3, L/61 then N/30 | 157 | 16 | 24%, 0.057 | 44%, 0.090 | 63%, 0.166 | 93%, 3.418 |
| corr3, LS + field + L/61 + N/30 | 433 | 70 | 23%, 0.057 | 39%, 0.080 | 65%, 0.219 | 94%, 3.348 |
| rand, GO4 best (field + coarse) | 177 | 9 | 71%, 0.314 | 85%, 0.854 | 95%, 1.862 | 99%, 4.535 |
| rand, L/61 | 70 | 10 | 75%, 0.370 | 88%, 0.830 | 91%, 1.903 | 100%, 4.526 |
| rand, L/183 | 196 | 18 | 51%, 0.123 | 66%, 0.786 | 81%, 1.989 | 94%, 4.424 |
| rand, N/30 | 87 | 8 | 72%, 0.302 | 82%, 0.882 | 92%, 2.046 | 98%, 4.844 |
| rand, L/61 then N/30 | 157 | 15 | 64%, 0.264 | 76%, 0.824 | 80%, 1.973 | 94%, 4.803 |
| rand, LS + field + L/61 + N/30 | 433 | 57 | 67%, 0.268 | 73%, 0.770 | 82%, 1.883 | 93%, 4.851 |

- The step ladder with enough budget per leg is the step forward: L/183 (61 per leg) takes corr3
  to 9%, 11% and 47% above 0.1 deg at 0.5, 1 and 2 deg (medians 0.021, 0.018, 0.080), from 32%,
  53% and 93% for GO4's best, at less than half its evaluations. With 20 per leg (L/61) it is
  worse than GO4's best: the budget per leg matters more than the arm.
- Neighbour continuation helps when neighbours are close (corr3) and does little on rand.
- rand stays hard (51% at 0.5 deg with L/183). It is a search problem, not an objective problem:
  on 150 voxels per level of both samples, the strict cost at the ground truth is lower than at
  the fit's end for practically every voxel left above 0.1 deg (median strict cost at the ground
  truth 0.065 to 0.096; at the end 0.10 to 0.71). So the strict cost after a fit flags the voxels
  that need a wider search.
- Next (I36, I37): a gated wider search for voxels whose strict cost stays high, and a larger
  ladder (a 3 deg first leg, more evaluations) for the 2 and 5 deg levels.

## GO4c: larger ladders (idea I37; I36 blocked)

Question: does more reach and budget per voxel recover the 2 and 5 deg levels? Same setup as
GO4b. Ladders `L/244/3` and `L/366/3` use four legs (3, 1, 0.3, 0.05 deg; 61 and 91 evaluations
per leg). Tables `benchmarks/grain_opt/go4c/`. Reading the table: share of voxels more than 0.1
deg from the ground truth, median error (deg) and share wrong (above 1 deg) per level.

| arm | eval/voxel | ms/voxel | 0.5 deg | 1 deg | 2 deg | 5 deg |
|---|---|---|---|---|---|---|
| corr3, L/183 (GO4b) | 196 | 18 | 9%, 0.021, 0% | 11%, 0.018, 0% | 47%, 0.080, 19% | 91%, 3.50, 89% |
| corr3, L/244/3 | 261 | 26 | 9%, 0.020, 0% | 15%, 0.018, 0% | 40%, 0.038, 11% | 77%, 2.33, 74% |
| corr3, L/366/3 | 401 | 35 | 3%, 0.012, 0% | 9%, 0.009, 0% | 31%, 0.011, 10% | 72%, 2.15, 69% |
| corr3, L/183 then N/30 | 283 | 26 | 6%, 0.020, 0% | 10%, 0.018, 0% | 40%, 0.058, 10% | 83%, 2.21, 75% |
| rand, L/183 (GO4b) | 196 | 18 | 51%, 0.123, 2% | 66%, 0.786, 40% | 81%, 1.99, 74% | 94%, 4.42, 92% |
| rand, L/244/3 | 261 | 26 | 50%, 0.104, 2% | 64%, 0.781, 39% | 74%, 1.90, 70% | 88%, 5.17, 87% |
| rand, L/366/3 | 401 | 36 | 48%, 0.047, 3% | 53%, 0.507, 36% | 77%, 1.89, 75% | 83%, 5.35, 83% |
| rand, L/183 then N/30 | 283 | 26 | 47%, 0.082, 2% | 61%, 0.757, 38% | 71%, 1.82, 64% | 89%, 4.81, 86% |

- More budget per leg keeps paying on corr3: L/366/3 leaves 3%, 9% and 31% of voxels above 0.1
  deg at 0.5, 1 and 2 deg (medians near 0.01 deg); at 5 deg 72% remain (69% wrong).
- rand gains less: the median at 0.5 deg falls to 0.047 deg but about half the voxels stay above
  0.1 deg. Its per-voxel landscape is rougher (neighbours' spots overlap and differ), and the
  ladder's first legs still start from the grain centre.
- The 5 deg level (both samples) needs a global search per voxel; I36 (gated seed path) waits for
  a cached ray index.

## GO5: strain in the forward model (ideas I14, I15, I16, I30)

Built by an implementer agent (branch `feature/grain-opt-strain`, merged after tests and an
independent review). `SpotModel` takes an optional per-voxel crystal-frame strain (6 numbers):
g = O M g0 with M = (I + eps)^-1, so the Bragg angle, omega and sin(2 theta) come from the strained
|g|. Tests pin the convention: a pure dilation d gives |g| = |g0|/(1 + d) and a smaller 2-theta; a
uniaxial e11 equals replacing every reciprocal vector by M g0 and differs from applying M in the
sample frame (added after review); a shear matches finite differences. The stock renderer has the
same option (`apply_deformation`, `render_full.py --strain`), default off and bit-identical;
rendered strained spots match `SpotModel` within 0.05 px. The `.mic` file holds M in columns 13 to
18, optionally written with 9 digits (M is held as float32 on read: about 6e-8 near 1). After
review, `GrainProblem` re-expresses the strain when it swaps a member orientation for a symmetry
equivalent (eps' = S^T eps S); before, strain fitting on such members would have been wrong.

Sample `strain`: the `grad` ground truth plus one constant strain tensor per grain (tensor RMS
sqrt(tr(eps^2)/3) of 1e-4, 5e-4 or 1e-3; 154, 156 and 187 grains). Rendered clean (84 s).

Reading line: for 5 grains at their ground-truth orientation and strain, the median strain step
(one component moved at a time) that raises a cost by 10%, and the spot shift per 1e-3 strain.

| component | coverage off lit pixels | intensity-weighted Sinkhorn | pixels per 1e-3 (det0 / det1) |
|---|---|---|---|
| e11 | 9e-5 | 2.2e-3 | 0.20 / 0.30 |
| e22 | 1.0e-4 | 1.5e-3 | 0.20 / 0.30 |
| e33 | 1.1e-4 | 8.7e-4 | 0.24 / 0.36 |
| e23 | 5.6e-5 | 8.8e-4 | 0.34 / 0.51 |
| e13 | 7.6e-5 | 9.0e-4 | 0.35 / 0.52 |
| e12 | 5.9e-5 | 1.2e-3 | 0.30 / 0.44 |

- In this near-field geometry a 1 px spot shift needs about 2e-3 to 5e-3 strain (0.2 to 0.5 px per
  1e-3): pixel/distance (4.4e-4, 2.7e-4) divided by a lever arm of order 2 tan(theta). Far field
  (HEXRD) gains by a much larger distance; here strain below about 1e-3 moves spots by well under
  a pixel.
- The Sinkhorn cost is flat to strain (a 10% rise needs about 1e-3 per component), as it was flat
  to sub-0.04 deg rotations (GO1). The coverage cost reacts near 1e-4, but those steps sit at the
  2e-5 grid floor of the scan: read them as at most 1e-4, not as a resolution.
- Leaving the strain out raises the Sinkhorn cost by 0.3% (a 1e-4 grain), 1.5% and 7.4% (5e-4)
  and 14% and 44% (1e-3): strain of 1e-3 is visible to the objective; 1e-4 hardly.
- Next (G1c): fit strain per grain with the strict cost (or the coverage cost) on `strain`, with
  orientation free, and measure recovered strain against the ground truth per level; per-grain
  fitting pools many spots, as the grain field did for orientation (GO2).

## GO5b: per-grain strain fit (idea I38)

Question: can a per-grain fit recover the elastic strain of sample `strain` from clean images?
Benchmark S12: 12 grains (4 per strain level, 25 to 70 voxels). Cost: the GO5 coverage cost (1
minus the share of predicted coverage mass on lit pixels). Parameters: grain field (9 numbers)
and one crystal-frame strain (6 numbers), CMA, strain from 0. Script `fit_strain.py`; rows in
`benchmarks/grain_opt/go5b/*.json`. Timing contended (Newt running), about 30 ms per evaluation;
the three arms took 2 to 4 minutes each.

Reading line: median over 4 grains of the strain error (tensor RMS of fitted minus ground-truth
strain), as a multiple of the level in brackets; "short" counts grains whose final cost is above
the cost at the ground truth (the optimizer stopped before the ground truth's cost).

| arm | 1e-4 | 5e-4 | 1e-3 | short | orientation median |
|---|---|---|---|---|---|
| `E@GT` strain only, orientation at ground truth, 300 ev. | 1.7e-4 (1.73) | 1.9e-4 (0.37) | 2.0e-4 (0.20) | 5 of 12 | (fixed) |
| `FE@GT` field + strain from the ground truth, 600 ev. | 3.1e-4 (3.05) | 3.4e-4 (0.69) | 3.8e-4 (0.38) | 9 of 12 | 0.010 deg |
| `FE@P` field 200 ev. + joint 600 ev., painted start | 6.7e-4 (6.7) | 6.3e-4 (1.26) | 5.3e-4 (0.53) | 10 of 12 | 0.012-0.016 deg |

Verdicts by the rule written before the runs (recovered at most 0.5 x level, partial 0.5 to 1):

- `E@GT`: 1e-4 not recovered, 5e-4 and 1e-3 recovered. The error sits near 2e-4 at every
  level: the resolution of the coverage cost on clean images, about half of pixel/distance
  (2.7e-4 to 4.4e-4). 1e-4 is below it: two grains end at a lower cost than the ground truth with
  a large error (objective shortfall, the cost is flat there, often already 0).
- `FE@GT`: 1e-4 not, 5e-4 partial, 1e-3 recovered. `FE@P`: 1e-4 and 5e-4 not, 1e-3 partial.
  `FE@P` is worse than `FE@GT` at 5e-4 (1.8 x, 3 of 4 grains), not at 1e-3 (1.4 x).
- In the joint fits 19 of 24 grains stop above the ground truth's cost: 15-number CMA with 600
  evaluations is optimizer-limited, not objective-limited, at 5e-4 and 1e-3. The orientation
  stays good (median 0.010 to 0.016 deg from the ground truth); strain error is what is left.
- Start bug found and fixed: with true-grain membership, 2 or 3 painted voxels that are off by
  many degrees pulled the least-squares start field (start cost 0.79 and 0.99 on two 1e-3
  grains). A robust start (medoid voxel; members over 1.5 deg from it set to it) fixed both
  (1e-3 median 9.3e-4 to 5.3e-4); the first run is kept as `FE_P_naive_start.json`.
- Next (GO5c): reach the cost optimum in the joint fit: staged legs (field, then strain only,
  then joint), CMA restarts or more evaluations; target `E@GT`'s 2e-4 floor from the painted
  start. Below 2e-4 needs a sharper cost (sub-pixel, e.g. intensity-weighted coverage).

## GO5c: joint strain fit that reaches the optimum (idea I39)

Question: GO5b's joint fits stopped above the ground truth's cost; does a longer or staged search
reach the optimum and bring the painted-start strain error to the `E@GT` floor (1.7e-4 /
1.9e-4 / 2.0e-4)? Same S12 grains, coverage cost and robust painted start; rows in
`benchmarks/grain_opt/go5c/*.json`. Contended timing, 5 to 12 minutes per arm.

Reading line: median strain error over 4 grains per level, multiple of the `E@GT` floor in
brackets; "short" = grains ending above the ground truth's cost.

| arm | evaluations | 1e-4 | 5e-4 | 1e-3 | short |
|---|---|---|---|---|---|
| GO5b `FE@P` (field 200, joint 600) | 826 | 6.7e-4 | 6.3e-4 (3.3) | 5.3e-4 (2.6) | 10 of 12 |
| `J@P` field 200, joint 2400 | 2626 | 6.5e-4 | 5.7e-4 (3.0) | 4.5e-4 (2.2) | 2 of 12 |
| `S@P` field 200, strain 300, joint 600 (small steps) | 1133 | 7.1e-4 | 3.4e-4 (1.8) | 3.0e-4 (1.45) | 7 of 12 |
| `SR@P` `S@P` + 2 joint legs of 300 | 1759 | 6.9e-4 | 3.7e-4 (2.0) | 3.4e-4 (1.6) | 5 of 12 |
| `S@GT` `S@P` legs from the ground-truth field | 922 | 2.0e-4 | 1.3e-4 (0.7) | 3.4e-4 (1.7) | 5 of 12 |

Verdicts by the rule written before the runs:

- No arm answers I39 (at most 4 short and the floor at 5e-4 and 1e-3). `S@P` reaches the floor
  at 1e-3 only (1.45 x); `S@P` and `SR@P` beat `FE@P` at 1e-3 (3 of 4 grains), not at 5e-4.
- The main finding corrects GO5b. `J@P` reaches the optimum (2 of 12 short; in 10 of 12 grains
  its cost is below the ground truth's, median 10% lower) yet keeps strain errors of 4.5e-4 to
  6.5e-4. With the orientation field free, the coverage cost's minimum is not at the ground-truth
  strain: the fit is objective-limited, not optimizer-limited. A 0.01 deg rotation (1.7e-4 rad)
  moves spots about as much as 2e-4 to 4e-4 strain, and the fitted fields sit 0.008 to 0.018 deg
  from the ground truth, so orientation and strain trade against each other.
- Staging helps because it stops earlier near a better point (strain first at the fitted field),
  not because it optimizes better: `S@P` ends above the ground truth's cost more often than
  `J@P` but with lower errors. From the ground-truth field (`S@GT`) the error is at the floor at
  1e-4 and 5e-4.
- Next: break the trade-off rather than search harder: fix the orientation from a cost the
  strain hardly moves (e.g. spot positions summed over all reflections, or the field from the
  unstrained fit) and fit strain alone (I41); or a sharper sub-pixel cost (I40) that separates
  rotation from strain by their different patterns across reflections.

## GO5d: fix the orientation, then strain alone (idea I41)

Question: if the grain field is fixed first from a fit that does not fit strain, and strain
alone is fitted after, does the painted-start strain error reach the `E@GT` floor? Same S12
grains, coverage cost, robust painted start; the last leg of every arm is strain only (300
evaluations) at a fixed field. Rows in `benchmarks/grain_opt/go5d/*.json`. 2 to 4 minutes per
arm.

Reading line: median strain error over 4 grains per level, median orientation error of the fixed
field (deg) in brackets.

| arm | evaluations | 1e-4 | 5e-4 | 1e-3 | short |
|---|---|---|---|---|---|
| GO5b `E@GT` (ground-truth orientations) | 309 | 1.7e-4 | 1.9e-4 | 2.0e-4 | 5 of 12 |
| `E@F` field fitted to the ground truth | 309 | 1.6e-4 (0.00002) | 1.3e-4 (0.00002) | 3.3e-4 (0.00002) | 5 of 12 |
| `O@P` field on coverage at zero strain, 600 | 920 | 1.5e-4 (0.009) | 5.3e-4 (0.022) | 6.0e-4 (0.040) | 8 of 12 |
| `Ost@P` field on the strict cost, 200 | 309 + strict | 4.7e-4 (0.015) | 5.4e-4 (0.021) | 6.3e-4 (0.039) | 11 of 12 |
| GO5c `S@P` (staged joint) | 1133 | 7.1e-4 (0.015) | 3.4e-4 (0.012) | 3.0e-4 (0.012) | 7 of 12 |

Verdicts by the rule written before the runs:

- I41 is answered no. Neither painted arm reaches the floor at 5e-4 or 1e-3 (needed at most
  2.8e-4 / 3.1e-4 against `E@GT`, or 1.9e-4 / 5.0e-4 against `E@F`), and neither beats `S@P` at
  those levels (`O@P` 1 and 0 of 4 grains better; `Ost@P` 1 and 1).
- A fit without a strain model is not strain-insensitive: the field absorbs the unmodelled
  strain as a rotation. The fixed field's orientation error grows with the strain level, 0.009 /
  0.022 / 0.040 deg on the coverage cost and 0.015 / 0.021 / 0.039 deg on the strict cost
  (about 0.04 deg, 7e-4 rad, per 1e-3 strain), against 0.012 deg for the staged joint fit.
  Strain fitted at that biased field is then worse than the joint fit's. At 1e-4 the bias is
  small and `O@P` sits at the floor (1.5e-4, 4 of 4 grains better than `S@P`; reported, not judged).
- The rule's `E@F` fallback does not apply as written: the `strain` sample's grains have one
  orientation each, so the order-1 field reproduces the ground truth to 2e-5 deg. `E@F` and
  `E@GT` are the same problem, and their difference (median 1.3e-4 vs 1.9e-4 at 5e-4, 3.3e-4 vs
  2.0e-4 at 1e-3; single grains differ up to 4x, e.g. grain 47 6.4e-5 vs 2.4e-4) is run-to-run
  scatter. The floor itself is a band of about 1e-4 to 3e-4 where the coverage cost is flat.
- Consequence: neither a staged search nor a strain-free orientation fit gets below the band;
  the cost has to resolve sub-pixel spot shifts. Next is I40 (sharper sub-pixel cost:
  intensity-weighted or distance to lit pixels), judged first at the ground-truth orientation
  (`E@GT`) and then in the staged joint fit from the painted start.

## GO5e: a sharper sub-pixel strain cost (idea I40)

Question: does a cost that sees sub-pixel shifts of the predicted coverage bring the strain-only
error at the ground-truth orientation below the coverage cost's band (reference 1.66e-4 /
1.57e-4 / 2.67e-4)? Costs on the same predicted points: `sd` (mass-weighted signed distance to
the lit set, trilinear) and `int` (1 - mass-weighted clean intensity / spot top). `E@GT`, seeds
0 and 1; rows in `benchmarks/grain_opt/go5e/*.json`. About 2 minutes per run, 15-33 ms per
evaluation (same as coverage).

| cost | 1e-4 | 5e-4 | 1e-3 | repeat scatter | short |
|---|---|---|---|---|---|
| coverage (reference, GO5b/GO5d) | 1.7e-4 | 1.6e-4 | 2.7e-4 | 3.5e-5 | 5 of 12 |
| `sd` (mean of 2 seeds) | 1.8e-3 | 1.2e-3 | 8.9e-4 | 5.8e-5 | 0 of 12 |
| `int` (mean of 2 seeds) | 2.5e-3 | 1.9e-3 | 1.5e-3 | 1.0e-4 | 0 of 12 |

Verdicts by the rule written before the runs:

- I40 is answered no for both costs: neither is sharper at any level, so the painted-start
  staged fit was not run (rule).
- Both costs are sharp but biased. Their repeat scatter is as small as coverage's, every fit
  ends below the ground truth's cost (0 of 12 short; `sd` at grain 111 -0.263 against -0.229),
  and the errors are 3e-4 to 3.6e-3. A probe along one random strain direction (grains 47, 237,
  111) showed `sd` smooth with its minimum within 5e-5 of the ground truth, so the bias lies in
  other directions: both costs reward strains that shrink or move the predicted footprint into
  the interior of the observed spot (deepest signed distance, brightest pixels) instead of
  matching it. A one-sided "predicted onto observed" cost is minimised by a footprint smaller
  than the spot, and strain changes footprint size.
- The coverage cost has the same one-sided form but saturates at 1 once all mass is lit, which
  is why it is flat (GO5d) instead of biased. A sub-pixel cost has to be two-sided (observed onto
  predicted as well) or compare positions, not depths.
- Next (I42): a two-sided cost: per-spot centroid match (predicted mass centroid against the
  observed intensity centroid, in pixel and frame units, the far-field HEXRD approach), or the
  symmetric signed distance (`sd` plus the observed pixels' distance to the predicted
  footprint).

## GO5f: a two-sided sub-pixel strain cost (idea I42)

Question: does a cost that compares predicted and observed spots from both sides beat the
coverage band at the ground-truth orientation (reference 1.66e-4 / 1.57e-4 / 2.67e-4), and then
in the painted-start staged fit? Costs: `cen` (per-spot squared distance between the predicted
mass centroid and the observed clean-intensity centroid, truncated at 4 px^2) and `ch` (chamfer:
predicted mass off the lit set plus lit pixels' distance to the nearest predicted point; the
plan's "minus its value at the ground truth" offset was left out, as a constant per spot it does
not move the minimum). Rows in `benchmarks/grain_opt/go5f/*.json`; 2 minutes per `E@GT` run,
5 minutes for `S@P`.

| cost, arm | 1e-4 | 5e-4 | 1e-3 | repeat scatter | orientation (deg) |
|---|---|---|---|---|---|
| coverage `E@GT` (reference) | 1.7e-4 | 1.6e-4 | 2.7e-4 | 3.5e-5 | (ground truth) |
| `cen` `E@GT` (mean of 2 seeds) | 1.3e-4 | 2.2e-4 | 1.7e-4 | 4.7e-6 | (ground truth) |
| `ch` `E@GT` (mean of 2 seeds) | 1.7e-3 | 9.0e-4 | 7.4e-4 | 3.9e-5 | (ground truth) |
| coverage `S@P` (GO5c) | 7.1e-4 | 3.4e-4 | 3.0e-4 | | 0.012 |
| `cen` `S@P` | 1.4e-4 | 2.2e-4 | 3.7e-4 | | 0.058 / 0.059 / 0.085 |

Verdicts by the rule written before the runs:

- I42 is answered no: `cen` is sharper only at 1e-3 (1.7e-4 against the needed 1.8e-4), not at
  1e-4 (1.3e-4 against 1.1e-4) or 5e-4; `ch` at no level. `ch` is biased like GO5e's costs.
- `cen` is the most consistent cost so far: repeat scatter 5e-6, no grain above 2.6e-4 at the
  ground-truth orientation (coverage had 7.7e-4 and 8.0e-4 on grains 111 and 435), every fit
  ends below the ground truth's cost (a small bias of about 1-2e-4, not a flat band). At the
  ground truth the per-axis centroid differences have medians near 0 (grain 47: column -0.01,
  row 0.00, frame -0.03 px); the large mean comes from spots merged with a neighbour's.
- `S@P` was run because `cen` was sharper at one level (my reading of "only if a cost is
  sharper"). From the painted start it reaches the `cen` floor at 1e-4 and 5e-4 and beats the
  coverage `S@P` there on 4 of 4 grains (1.4e-4 and 2.2e-4 against 7.1e-4 and 3.4e-4); at 1e-3
  it misses (3.7e-4, 1 of 4 better; grain 237 runs off to 0.29 deg and 1.5e-3).
- The price is the orientation: 0.06-0.09 deg median voxel error against 0.012 deg for coverage.
  A centroid fixes each spot's position, which depends on the grain's mean orientation, not on
  the field's gradient, so the gradient drifts.
- Next (I43): a hybrid, coverage for the orientation field and the centroid cost for strain:
  either the sum (centroid + weighted coverage) or staged (field on coverage, then strain and
  mean orientation on centroids), with merged spots down-weighted.

## GO5g: hybrid cost, coverage for orientation and centroids for strain (idea I43)

Question: does `cen` + lambda x coverage keep GO5f's strain from the painted start and
coverage's orientation? Benchmark S12; rows in `benchmarks/grain_opt/go5g/*.json`. 3 minutes for
`E@GT`, 9 to 11 minutes per painted-start run (about 60 ms per evaluation, contended).

Reading line: median strain error over 4 grains per level, median voxel orientation error (deg)
in brackets; the rule's limits are 1.9e-4 / 3.3e-4 / 2.5e-4 and 0.02 deg.

| arm | 1e-4 | 5e-4 | 1e-3 | short |
|---|---|---|---|---|
| GO5c coverage `S@P` | 7.1e-4 (0.015) | 3.4e-4 (0.012) | 3.0e-4 (0.012) | 7 of 12 |
| GO5f `cen` `S@P` | 1.4e-4 (0.058) | 2.2e-4 (0.059) | 3.7e-4 (0.085) | 0 of 12 |
| `hyb` lambda 2 `E@GT` (ground-truth orientation) | 9.9e-5 | 1.4e-4 | 1.1e-4 | 1 of 12 |
| `hyb` lambda 2 `S@P` | 1.3e-4 (0.007) | 1.2e-4 (0.009) | 2.5e-4 (0.021) | 3 of 12 |
| `hyb` lambda 8 `S@P` | 1.9e-4 (0.009) | 1.5e-4 (0.008) | 2.2e-4 (0.012) | 8 of 12 |
| `Sx@P` (coverage field, `cen` strain, `hyb` 2 joint) | 1.4e-4 (0.006) | 1.5e-4 (0.009) | 2.5e-4 (0.021) | 2 of 12 |

Verdicts by the rule written before the runs:

- I43 is answered yes by `hyb` lambda 8 `S@P`: strain 1.5e-4 and 2.2e-4 at 5e-4 and 1e-3 (limits
  3.3e-4, 2.5e-4), orientation 0.008-0.012 deg at every level, and better than `cen` `S@P` at
  1e-3 on 4 of 4 grains (2.2e-4 against 3.7e-4). It is the best arm meeting the rule (the only
  one).
- lambda 2 (`S@P` and `Sx@P`) misses at 1e-3 by a hair (2.5e-4 and 0.021 deg; grains 237 and 54
  keep 0.025-0.046 deg): more coverage weight holds the field in place on the strongly strained
  grains. The by-role staging adds nothing over the plain hybrid.
- Against GO5c, from the same painted start: strain error 1.9 / 1.5 / 2.2e-4 instead of 7.1 /
  3.4 / 3.0e-4, with the orientation as good or better (0.008-0.012 deg). At the ground-truth
  orientation the hybrid is below the coverage band at every level (9.9e-5 / 1.4e-4 / 1.1e-4).
  The two costs see different things: coverage pins the orientation field and the spot shapes,
  centroids pin the sub-pixel spot positions that strain moves.
- Caveats: 4 grains per level, one seed; the rule margin at 1e-3 is small (2.2e-4 against 2.5e-4).
  Next (I44): confirm on more grains (12 per level) with lambda 8 and 16, then a faster solver
  on the centroid residuals (I7 / I19, Gauss-Newton) in place of CMA.

## GO5h: confirm the hybrid strain fit on more grains (idea I44, first part)

Question: does GO5g's hybrid (`cen` + lambda x coverage, staged from the painted start) hold on
12 grains per level (S36), and is lambda 16 better than 8? Rows in
`benchmarks/grain_opt/go5h/*_S36.json`; 16 minutes for coverage, 32 minutes per hybrid arm (3
processes, contended).

Reading line: median strain error over 12 grains (75th percentile in brackets), grains better
than coverage on the same grain, median voxel orientation error (deg).

| arm | 1e-4 | 5e-4 | 1e-3 | orientation |
|---|---|---|---|---|
| coverage `S@P` | 4.6e-4 (6.9e-4) | 2.8e-4 (4.3e-4) | 6.3e-4 (8.7e-4) | 0.018 / 0.009 / 0.020 |
| `hyb` lambda 8 | 2.3e-4 (2.5e-4), 10 of 12 | 1.5e-4 (2.1e-4), 9 of 12 | 2.1e-4 (5.7e-4), 9 of 12 | 0.011 / 0.009 / 0.017 |
| `hyb` lambda 16 | 2.3e-4 (3.2e-4), 10 of 12 | 1.6e-4 (1.8e-4), 9 of 12 | 2.3e-4 (5.0e-4), 11 of 12 | 0.011 / 0.009 / 0.019 |

Verdicts by the rule written before the runs:

- Both hybrid arms confirm (limits: coverage / 1.5 = 1.9e-4 at 5e-4 and 4.2e-4 at 1e-3, at
  least 8 of 12 better, orientation at most 0.02 deg): I44's first part is yes. The better
  lambda by the rule is 8 (mean of the 5e-4 and 1e-3 medians 1.8e-4 against 1.9e-4); the two are
  within run scatter, and lambda 16 has the shorter tail at 1e-3 (75th percentile 5.0e-4 against
  5.7e-4, 11 of 12 better).
- Grains whose orientation ends above 0.05 deg: coverage 1 (485, 0.155 deg); lambda 8 two (485
  and 0, 0.059 deg); lambda 16 two (485 0.062, 432 0.055). Grain 485 is hard for every cost.
- On 36 grains the coverage fit is worse than on S12 at 1e-3 (6.3e-4 against 3.0e-4), the hybrid
  about the same (2.1e-4 against 2.2e-4): the S12 coverage numbers were on the lucky side.
- The hybrid's remaining weak spot is a tail of strongly strained (1e-3) grains at 4-9e-4.
- Next (I44 second part): replace CMA (about 50 s per grain, 1,100 evaluations of about 60 ms)
  by Gauss-Newton on the centroid residuals for strain and mean rotation, with the orientation
  field's gradient from the coverage fit; judge on S36 against `hyb` lambda 8 for accuracy and
  wall time.

## GO5i: Gauss-Newton on centroid residuals (idea I44, second part)

Question: does Levenberg-Marquardt on per-spot centroid residuals (mean rotation + strain, the
field gradient fixed from the coverage field leg) match the CMA hybrid of GO5h at a fraction of
its time? S36, painted start; rows in `benchmarks/grain_opt/go5i/*_S36.json`. Plan amended
before the judged runs: Tukey weights with the cut-off annealed 3 / 1.5 / 1 px instead of Huber
at 2 px (a one-grain smoke test showed Huber pulled by spots merged with a neighbour's: centroid
cost below the ground truth's, 6e-4 strain error). Load 1.5-3.4 (GO5h: about 3.5).

Reading line: median strain error over 12 grains (75th percentile), median voxel orientation
error (deg), mean wall time per grain.

| arm | 1e-4 | 5e-4 | 1e-3 | orientation | s per grain |
|---|---|---|---|---|---|
| GO5h `hyb` lambda 8 (CMA, reference) | 2.3e-4 (2.5e-4) | 1.5e-4 (2.1e-4) | 2.1e-4 (5.7e-4) | 0.011 / 0.009 / 0.017 | 53 |
| `GN@P` coverage field leg + Gauss-Newton | 1.2e-4 (1.5e-4) | 1.3e-4 (1.5e-4) | 2.9e-4 (7.6e-4) | 0.010 / 0.021 / 0.039 | 10 |
| `GNH@P` + hybrid CMA polish 300 | 1.6e-4 (2.1e-4) | 1.9e-4 (2.2e-4) | 2.3e-4 (4.5e-4) | 0.011 / 0.012 / 0.025 | 23 |

Verdicts by the rule written before the runs:

- I44's second part is answered no: neither arm matches. `GN@P` fails orientation at 5e-4 and
  1e-3 (0.021, 0.039 deg) and strain at 1e-3 (2.9e-4 against the limit 2.7e-4); `GNH@P` meets
  the strain limits and the time (23 s against 27 s) but its orientation at 1e-3 is 0.025 deg.
- Gauss-Newton is the most accurate strain fit at low strain so far, at a fifth of the time:
  1.2e-4 / 1.3e-4 at 1e-4 / 5e-4 against the CMA hybrid's 2.3e-4 / 1.5e-4.
- Its weakness is the field gradient: it is fitted once, on coverage at zero strain, and GO5d
  showed that fit takes up strain as rotation in proportion to the strain. Gauss-Newton refits
  only the mean rotation, so at 1e-3 the voxel orientations keep the gradient error (grains 452
  0.25 deg, 485 0.16 deg).
- Next (I45): alternate, refit the field gradient on coverage at the Gauss-Newton strain, then
  Gauss-Newton again (two rounds); judged like GO5i.

## GO5j: alternate field-gradient refit and Gauss-Newton (idea I45)

Question: does refitting the field gradient on coverage at the Gauss-Newton strain (CMA 100 over
`theta[3:]`), then Gauss-Newton again, fix GO5i's orientation and match the CMA hybrid at half
its time? S36, painted start; rows in `benchmarks/grain_opt/go5j/*_S36.json`; load 2.3-2.8.

| arm | 1e-4 | 5e-4 | 1e-3 | orientation | s per grain |
|---|---|---|---|---|---|
| GO5h CMA hybrid (reference) | 2.3e-4 | 1.5e-4 | 2.1e-4 (p75 5.7e-4) | 0.011 / 0.009 / 0.017 | 53 |
| GO5i `GN@P` | 1.2e-4 | 1.3e-4 | 2.9e-4 (p75 7.6e-4) | 0.010 / 0.021 / 0.039 | 10 |
| `GA1@P` one round | 1.2e-4 | 1.4e-4 | 3.1e-4 (p75 7.6e-4) | 0.011 / 0.013 / 0.032 | 17 |
| `GA2@P` two rounds | 1.2e-4 | 1.2e-4 | 3.1e-4 (p75 7.5e-4) | 0.009 / 0.014 / 0.034 | 24 |

Verdicts by the rule written before the runs:

- I45 is answered no: both arms meet the time and the 5e-4 limits (orientation 0.013-0.014 deg,
  fixed by the gradient refit) but miss at 1e-3 (strain 3.1e-4 against 2.7e-4, orientation
  0.032-0.034 deg against 0.02).
- At 1e-3 the failures are a few grains, and their pattern points at the Gauss-Newton loss, not
  the gradient: grains 452, 0, 485 and 349 end with strain errors 6e-4 to 4.5e-3 while the
  weighted centroid cost falls to about 0. The Tukey weights are recomputed every iteration, so
  Gauss-Newton can "escape" by pushing residuals past the 1 px cut-off, where they no longer
  count. These grains also have few usable spots (18 to 35, against 50 to 90 elsewhere): strain
  of 1e-3 moves many spots out of the windows built at the unstrained start.
- On the grains with enough spots, Gauss-Newton with the gradient refit is as good as the CMA
  hybrid or better (e.g. 439 8e-5 against 1.5e-4, 167 8e-5 against 2.0e-4) at less than half
  the time.
- Next (I46): stop the escape. Freeze each stage's inlier set and weights at the stage start,
  require at least half the spots as inliers (else keep the previous stage), and fall back to
  the CMA hybrid for grains with too few spots; also rebuild the spot windows at the fitted
  strain.

## GO5k: stop the Gauss-Newton escape (idea I46)

Question: with per-stage frozen Tukey weights (a stage skipped below half inliers) and a CMA
hybrid fallback for grains with fewer than 40 spots, does the fast fit match the CMA hybrid at
every level? S36, painted start; rows in `benchmarks/grain_opt/go5k/*_S36.json`; load 2.5-3.8.

| arm | 1e-4 | 5e-4 | 1e-3 | orientation | s per grain |
|---|---|---|---|---|---|
| GO5h CMA hybrid (reference) | 2.3e-4 | 1.5e-4 | 2.1e-4 (p75 5.7e-4) | 0.011 / 0.009 / 0.017 | 53 |
| GO5j `GA2@P` | 1.2e-4 | 1.2e-4 | 3.1e-4 (p75 7.5e-4) | 0.009 / 0.014 / 0.034 | 24 |
| `GF2@P` frozen weights | 1.3e-4 | 1.3e-4 | 5.1e-4 (p75 1.0e-3) | 0.009 / 0.015 / 0.041 | 17 |
| `GFF@P` frozen weights + fallback (5 grains) | 1.3e-4 | 1.3e-4 | 2.4e-4 (p75 6.9e-4) | 0.009 / 0.015 / 0.028 | 19 |

Verdicts by the rule written before the runs:

- I46 is answered no: `GFF@P` meets the strain limits (2.4e-4 at 1e-3 against 2.7e-4) and the
  time (19 s against 27 s) but its orientation at 1e-3 is 0.028 deg (limit 0.02). `GF2@P` fails
  strain and orientation at 1e-3.
- Frozen weights stop the escape, but where too few spots stay inliers every stage is skipped
  and the strain stays at 0 (error = 1e-3): grains 452 and 0 in `GF2@P`, and grains 273 and 232
  (48 and 56 spots, above the fallback threshold) in `GFF@P`. The fallback grains (485, 479, 452,
  349, 0) end close to the CMA reference, and are quick (10-20 s: they are small).
- On the Gauss-Newton grains at 1e-3 the orientation stays 0.03-0.04 deg (432, 351, 273, 232)
  where CMA reaches 0.006-0.021: the gradient refit (CMA 100 on coverage) and the mean rotation
  from centroids do not pin the field as well as the joint hybrid.
- Next (I47): fall back also when Gauss-Newton skips its stages, and polish with a short joint
  CMA leg on the hybrid only for grains whose fitted strain is large (RMS above 6e-4), so the
  time stays near half the reference.

## GO5l: fallback on stalls, polish large strains (idea I47)

Question: GO5k's `GFF@P` plus a CMA-hybrid fallback when Gauss-Newton skips every first-round
stage, and a joint hybrid CMA polish (400 or 600 evaluations, 0.02 deg and 1e-4) only for grains
whose Gauss-Newton strain RMS exceeds 6e-4: does it match the CMA hybrid within half its time?
S36, painted start; rows in `benchmarks/grain_opt/go5l/*_S36.json`; load 2.4-3.4.

| arm | 1e-4 | 5e-4 | 1e-3 | orientation | s per grain |
|---|---|---|---|---|---|
| GO5h CMA hybrid (reference) | 2.3e-4 | 1.5e-4 | 2.1e-4 (p75 5.7e-4) | 0.011 / 0.009 / 0.017 | 53 |
| GO5k `GFF@P` | 1.3e-4 | 1.3e-4 | 2.4e-4 (p75 6.9e-4) | 0.009 / 0.015 / 0.028 | 19 |
| `GP4@P` (polish 400) | 1.3e-4 | 1.3e-4 | 2.1e-4 (p75 3.5e-4) | 0.009 / 0.015 / 0.020 | 25 |
| `GP6@P` (polish 600) | 1.3e-4 | 1.3e-4 | 2.1e-4 (p75 3.5e-4) | 0.009 / 0.015 / 0.019 | 26 |

Verdicts by the rule written before the runs:

- I47 is answered yes: both arms match (strain at most 1.25 x the reference at 5e-4 and 1e-3,
  orientation at most 0.02 deg at every level, at most 26.7 s per grain). The faster, `GP4@P`,
  is preferred. Its orientation at 1e-3 (0.0197 deg) is just inside the limit.
- Grains: 9 fall back to the CMA hybrid (5 with fewer than 40 spots; 4 stalls: 140, 376, 273,
  232), 5 are polished (432, 439, 21, 351, 167, all at 1e-3), 22 use Gauss-Newton alone. Grain
  485 (0.175 deg) is the one left above 0.05 deg.
- Against the CMA hybrid: strain lower at 1e-4 (1.3e-4 against 2.3e-4) and 5e-4, equal at 1e-3
  with a shorter tail (75th percentile 3.5e-4 against 5.7e-4), orientation a little worse at
  5e-4 (0.015 against 0.009 deg), in half the time.
- The fast strain fit is now: robust start field, coverage field leg 200, Gauss-Newton on
  centroid residuals with frozen Tukey weights (3 / 1.5 / 1 px) and two coverage gradient
  refits, CMA hybrid where Gauss-Newton has too few spots or stalls, hybrid polish for strain
  above 6e-4 (`fit_strain.py --arm GP4@P --cost hyb --lam 8`).

## GO5m: strain fits on noisy images (idea I48)

Question: do the CMA hybrid and the fast fit keep their accuracy on images with the project's
detector noise model (`realistic_nodrop`: per spot, 5% edge-pixel flips, a hot pixel with
probability 0.05, a 2-4 px blob with probability 0.1; about 5% of lit pixels dropped and 4-5%
added per frame)? S36, painted start, `hyb` lambda 8; noisy stack built in this worktree
(`build_store.py --out_tag strainrn`); rows in `benchmarks/grain_opt/go5m/*_S36.json`.

Reading line: median strain error over 12 grains, clean -> noisy; median voxel orientation error
(deg) on noisy images.

| fit | 1e-4 | 5e-4 | 1e-3 | orientation (noisy) | s per grain |
|---|---|---|---|---|---|
| CMA hybrid | 2.3e-4 -> 2.5e-4 | 1.5e-4 -> 1.9e-4 | 2.1e-4 -> 2.8e-4 | 0.012 / 0.012 / 0.028 | 51 |
| fast fit `GP4@P` | 1.3e-4 -> 1.8e-4 | 1.3e-4 -> 2.7e-4 | 2.1e-4 -> 3.7e-4 | 0.014 / 0.017 / 0.025 | 27 |

Verdicts by the rule written before the runs:

- The CMA hybrid is robust to noise: its noisy strain is within 1.5 x its clean strain at 5e-4
  and 1e-3 (1.25 x and 1.3 x), orientation at most 0.03 deg (0.028 at 1e-3).
- The fast fit is not: its strain doubles at 5e-4 (2.7e-4 against the limit 2.0e-4) and grows
  1.7 x at 1e-3 (3.7e-4 against 3.2e-4); it no longer holds against the noisy CMA hybrid (2.7e-4
  and 3.7e-4 against 1.25 x 1.9e-4 and 2.8e-4). I48 is answered no. Its orientation stays within
  0.03 deg. Fallbacks 7, polished 8 (clean: 9 and 5).
- Reading: Gauss-Newton works at sub-pixel precision (1e-4 strain is about 0.03 px), so it feels
  the centroid noise. The noise model copies a neighbour's intensity into each added edge pixel
  and adds hot pixels and blobs near the spot; with intensity weights these move the observed
  centroid by a few hundredths of a pixel. The hybrid's coverage term does not care about a few
  pixels at the edge.
- Next (I49): a noise-robust observed centroid for Gauss-Newton: binary weights, and isolated
  lit pixels removed (fewer than 2 lit 8-neighbours in the frame), on the same noisy images.

## GO5n: noise-robust observed centroids (idea I49)

Question: on noisy images, does a binary-weight observed centroid (and one with isolated lit
pixels removed) restore the fast fit? S36, painted start, `GP4@P`, noisy `strainrn` images; rows
in `benchmarks/grain_opt/go5n/*.json`; load 2.3-3.7.

| centroid weights | 1e-4 | 5e-4 | 1e-3 | orientation | polished | s per grain |
|---|---|---|---|---|---|---|
| GO5m noisy CMA hybrid (reference) | 2.5e-4 | 1.9e-4 | 2.8e-4 | 0.012 / 0.012 / 0.028 | | 51 |
| GO5m intensity (before) | 1.8e-4 | 2.7e-4 | 3.7e-4 | 0.014 / 0.017 / 0.025 | 8 | 27 |
| `bin` binary | 4.9e-4 | 4.6e-4 | 6.7e-4 | 0.030 / 0.018 / 0.023 | 23 | 35 |
| `binclean` binary, isolated removed | 5.1e-4 | 5.0e-4 | 5.9e-4 | 0.022 / 0.020 / 0.028 | 22 | 35 |

Verdicts by the rule written before the runs:

- I49 is answered no: both variants are worse than intensity weights at every level (strain
  4.6-6.7e-4 against the limits 2.4e-4 and 3.5e-4) and slower (more grains cross the 6e-4 polish
  threshold, 22-23 against 8).
- The clean intensities carry the sub-pixel information: a pixel's intensity is its partial
  coverage by the spot, so the intensity-weighted centroid locates the spot well below a pixel,
  while the binary centroid is quantised to the lit-pixel set. Binary weights cost far more than
  the noise does.
- Next (I50): keep intensity weights and remove only isolated lit pixels (`intclean`), checked
  on both noisy and clean images.

## GO5o: intensity-weighted centroids without isolated pixels (idea I50)

Question: with intensity weights kept and only isolated lit pixels removed from the observed
centroid (`--cen_w intclean`), does the fast fit hold on noisy images without losing its clean
accuracy? S36, painted start, `GP4@P`; rows in `benchmarks/grain_opt/go5o/` (noisy) and
`go5o_clean/`; load 2.3-3.8 (higher than GO5m's 2.3-2.8).

| arm | 1e-4 | 5e-4 | 1e-3 | orientation | s per grain |
|---|---|---|---|---|---|
| noisy CMA hybrid (GO5m, reference) | 2.5e-4 | 1.9e-4 | 2.8e-4 | 0.012 / 0.012 / 0.028 | 51 |
| noisy fast fit, intensity centroids (GO5m) | 1.8e-4 | 2.7e-4 | 3.7e-4 | 0.014 / 0.017 / 0.025 | 27 |
| noisy fast fit, `intclean` | 2.2e-4 | 2.3e-4 | 3.3e-4 | 0.016 / 0.022 / 0.025 | 26.4 |
| clean fast fit (GO5l) | 1.3e-4 | 1.3e-4 | 2.1e-4 | 0.009 / 0.015 / 0.020 | 25 |
| clean fast fit, `intclean` | 1.3e-4 | 1.4e-4 | 2.3e-4 | 0.010 / 0.014 / 0.018 | 24.9 |

Verdicts by the rule written before the runs:

- Noisy: strain within the limits (2.3e-4 and 3.3e-4 against 2.4e-4 and 3.5e-4), orientation
  within 0.03 deg, but the mean time per grain is 26.4 s against the limit of 25.3 s (half the
  noisy CMA hybrid's 50.5 s): the arm misses by 1.1 s (4%), at a higher machine load than the
  reference run. Clean: no loss (1.4e-4 and 2.3e-4 against 1.7e-4 and 2.7e-4; orientation at
  most 0.018 deg). I50 is answered no, on time alone.
- Removing isolated pixels is the right direction: it takes the noisy strain at 5e-4 from 2.7e-4
  to 2.3e-4 and at 1e-3 from 3.7e-4 to 3.3e-4 at no clean cost, but the fast fit on noisy images
  is only about 1.2 x the noisy CMA hybrid's error at half its time, not better. On noisy images
  the CMA hybrid (2.5 / 1.9 / 2.8e-4, 51 s) and the fast fit `intclean` (2.2 / 2.3 / 3.3e-4,
  26 s) sit on the same frontier.
- The strain line is time-boxed here (program balance, PLAN "GO6 onwards"). Open strain ideas
  for later: a noise-aware residual weight per spot, windows rebuilt at the fitted strain.

## GO4d: cached ray index and the gated global search (goal G1b, idea I36)

Question: with the paired observed rays cached per sample, does the reconstructor's ray-indexing
seed path, run on voxels whose strict cost stays high after the ladder, recover the deformed
grains the ladder leaves wrong at 2 and 5 deg?

Part A (`ray_cache.py`): the one-time pairing took 1,152 s on `corr3` (426,061 rays from
475,861 + 530,208 spots) and 2,873 s on `rand` (623,313 rays); the saved arrays reload equal and
load in 0.01 s, so `index_seed` no longer pays the pairing per process.

Part B: `run_pipeline.py --start centre --group true`, 1,043 voxels per sample; tables
`benchmarks/grain_opt/go4d/{corr3,rand}.md`. Reading line: share of voxels more than 0.1 deg from
the ground truth (wrong, above 1 deg) per level; ms per voxel for the whole pipeline.

| sample, arm | 0.5 deg | 1 deg | 2 deg | 5 deg | gated / improved | ms / voxel |
|---|---|---|---|---|---|---|
| corr3, `L/366/3` (GO4c) | 3% (0%) | 9% (0%) | 31% (10%) | 72% (69%) | | 35 |
| corr3, `L/366/3,G/15/0.15` | 3% (0.4%) | 5% (0.4%) | 9% (4%) | 36% (28%) | 290 / 262 | 285 |
| corr3, `L/366/3,G/10/0.10` | 2% (0.4%) | 4% (0.4%) | 9% (4%) | 36% (28%) | 420 / 321 | 373 |
| rand, `L/366/3` (GO4c) | 48% (3%) | 53% (36%) | 77% (75%) | 83% (83%) | | 36 |
| rand, `L/366/3,G/15/0.15` | 26% (2%) | 15% (7%) | 52% (47%) | 44% (42%) | 637 / 541 | 592 |
| rand, `L/366/3,G/10/0.10` | 20% (2%) | 15% (7%) | 52% (47%) | 44% (42%) | 751 / 594 | 686 |

Verdicts by the rule written before the runs:

- I36 is answered yes for both arms: at 5 deg the share above 0.1 deg falls from 72% to 36% on
  corr3 (limit 48%) and from 83% to 44% on rand (limit 55%), and every level improves (rand at 1
  deg 53% -> 15%; corr3 at 2 deg 31% -> 9%). The cheaper arm, threshold 0.15, is preferred: the
  lower threshold gates 30% more voxels and changes only rand at 0.5 deg (26% -> 20%).
- Cost: about 0.3 s (corr3) and 0.6 s (rand) per voxel of the sample, against 0.035 s for the
  ladder alone; the seed path is about 1 s per gated voxel. The objective-evaluation count per
  voxel is large (1,100 and 3,500: the finisher runs FindOptimal), so this arm sits far right on
  the Pareto chart; wall time is the fairer measure.
- What remains: corr3 5 deg 28% wrong, rand 2 deg 47% and 5 deg 42% wrong. The 95th percentile
  error rises to 45 deg: some gated voxels take a seed from another grain (a lower strict cost
  at a wrong orientation, often a voxel the neighbours' spots dominate). rand at 2 deg is the
  weakest level (52%, worse than its 5 deg), where the ladder's end is near but outside the basin
  and the strict cost there is often below the gate.
- Next for G1b (later): a second pass with the improved neighbours as candidates (N stage after
  G), and a gate on the neighbour disagreement as well as the strict cost.

## GO7: rocking curve on the images (goal G2b, idea I23, part 1)

Question: how much does the orientation fit lose when the images have a rocking curve (spots
spread over omega frames by alpha / |sin eta|) but the forward model does not? Images emulated
from the clean `grad` intensities (`rock_stack.py`: overlap-weighted spread over frames, window
centred on the frame centre, capped at 8 frames, threshold at the 5th percentile of the clean
intensities). G20, painted start, `Astrict/100`; tables `benchmarks/grain_opt/go7/g20.md`.

| images | lit pixels | pixels spread past 1 frame | above 0.1 deg (853 C voxels) | median deg | grains rejected by the gate |
|---|---|---|---|---|---|
| clean | 9.44M | 0 | 3 (0.1 to 1.0%) | 0.0146 | 3 (37, 493, 64) |
| alpha 0.3 deg | 10.11M | 11% | 20 (1.5 to 3.6%) | 0.0149 | 4 (37, 493, 64, 197) |
| alpha 1.0 deg | 12.61M | 100% | 20 (1.5 to 3.6%) | 0.0134 | 2 (64, 197) |

Verdicts by the rule written before the runs:

- By the rule the rocking curve harms the fit at both alphas: the share above 0.1 deg has a
  Wilson interval entirely above the clean arm's (20 of 853 against 3). The median does not move
  (0.0149 and 0.0134 against 0.0146).
- The whole difference is one grain: grain 197 (54 voxels, 20 of them above 0.1 deg at the
  painted start) has its fitted stage rejected by the acceptance gate (mean strict cost must not
  get worse) under both rocking images, so it stays at the painted start. Every accepted grain
  fits as well as on clean images; at alpha 1.0 two grains the gate rejected on clean images
  (37, 493) are accepted.
- The mean strict cost rises with alpha (0.084 clean, 0.107, 0.144 after the fit): spreading
  weak spots over frames drops them below the threshold, so more predicted pixels land on unlit
  pixels. The single-frame model still finds the right orientation; the gate, which compares
  strict costs, is what becomes less reliable.
- Part 2 (a forward model with a rocking curve, sub-frame omega from frame ratios) changes
  `icenine/` and needs the owner's approval; given how little part 1 lost, it is low priority.

## GO7b: detector point-spread function on the images (goal G2b, idea I24)

Question: how much do the orientation fit and the fast strain fit lose when the images are
blurred by a Gaussian PSF before the threshold, while the forward model has none? Images from
the clean intensities (`psf_stack.py`; sigma 0.5 px adds 3% lit pixels, 1.0 px adds 25%).
Orientation: G20 `Astrict/100`, table `benchmarks/grain_opt/go7b/g20.md`. Strain: `GP4@P` on S36,
rows in `benchmarks/grain_opt/go7b/strain_psf*/`.

| images | orientation: above 0.1 deg (853), median | strain median 1e-4 / 5e-4 / 1e-3 | strain-fit orientation (deg) | s per grain |
|---|---|---|---|---|
| clean | 3, 0.0146 | 1.3 / 1.3 / 2.1e-4 | 0.009 / 0.015 / 0.020 | 25 |
| PSF sigma 0.5 px | 20, 0.0151 | 1.4 / 1.3 / 2.8e-4 | 0.013 / 0.020 / 0.026 | 21 |
| PSF sigma 1.0 px | 22, 0.0192 | 2.3 / 2.7 / 5.7e-4 | 0.029 / 0.039 / 0.061 | 18 |

Verdicts by the rule written before the runs:

- Orientation: harm at both sigmas by the share above 0.1 deg (Wilson intervals above the clean
  arm's); as in GO7 it comes from grain 197, which the strict-cost gate rejects (rejected grains
  37, 493, 64, 197 at 0.5 px; 37, 64, 197 at 1.0 px). At 1.0 px the median also rises 1.3 x
  (0.0192 against 0.0146; below the 1.5 x limit) and the 95th percentile from 0.035 to 0.063 deg.
- Strain: no harm at 0.5 px (within 1.5 x clean, orientation within 0.03 deg). Harm at 1.0 px:
  strain 2.7e-4 and 5.7e-4 at 5e-4 and 1e-3 (limits 2.0e-4, 3.2e-4) and orientation 0.039 and
  0.061 deg. 11 grains fall back to CMA (clean 9).
- Reading: a 1 px PSF is the first image nuisance in this program that clearly breaks the strain
  fit. Blur-then-threshold grows each observed spot by about a pixel while the predicted coverage
  keeps the sharp footprint, so the coverage term of the hybrid no longer has its optimum at the
  ground truth; the centroid term stays nearly unbiased (symmetric blur). The fit gets faster
  because more grains go to Gauss-Newton fallback and polish thresholds differently.
- Next (I52, inside the scripts, no library change): a PSF-aware predicted coverage in the strain
  costs (the predicted coverage points spread by the same Gaussian before scoring) on the 1.0 px
  images.

## GO7c: PSF-aware predicted coverage in the strain costs (goal G2b, idea I52)

Question: on the 1 px PSF images, does blurring the predicted coverage by the same PSF (4
Gaussian copies of each coverage point) restore the fast strain fit? S36, `GP4@P`; rows in
`benchmarks/grain_opt/go7c/`.

| arm | 1e-4 | 5e-4 | 1e-3 | orientation | s per grain |
|---|---|---|---|---|---|
| clean images (GO5l) | 1.3e-4 | 1.3e-4 | 2.1e-4 | 0.009 / 0.015 / 0.020 | 25 |
| 1 px PSF images, sharp prediction (GO7b) | 2.3e-4 | 2.7e-4 | 5.7e-4 | 0.029 / 0.039 / 0.061 | 18 |
| 1 px PSF images, PSF-aware coverage | 2.3e-4 | 2.6e-4 | 6.8e-4 | 0.033 / 0.046 / 0.094 | 17 |

Split by path (strain median, orientation median): Gauss-Newton-only grains clean 1.3e-4 /
0.011 deg, PSF 2.4e-4 / 0.032-0.035 deg with or without the PSF-aware coverage; fallback (CMA)
grains clean 2.2e-4 / 0.023 deg, PSF 6.4-7.5e-4 / 0.058-0.096 deg.

Verdicts by the rule written before the runs:

- I52 is answered no: no change in strain (gain 1.0 x at 1e-4 and 5e-4, worse at 1e-3) and
  orientation worse; the no-harm limits are missed at every level.
- The premise was wrong: the Gauss-Newton-only grains degrade just as much, and Gauss-Newton
  scores centroids, which a symmetric blur leaves unbiased. The likely cause is spot merging: a
  1 px blur before the threshold joins nearby spots into one connected component, so the
  observed centroid (and the lit set the coverage term scores) mixes neighbours. The fallback
  grains (fewer spots, more overlap) suffer most.
- Next for this line (not now; G2b is time-boxed like strain): split merged components by the
  predicted spot positions before taking centroids.

## GO9: Poisson likelihood on pixel counts as an orientation objective (goal G2a, idea I18)

Question: does a Poisson likelihood on the observed counts (per-reflection closed-form scale,
background b) have its optimum at the ground truth, where GO1's smooth objectives did not? GO1
protocol on `grad` grains 381 and 302 (`probe_poisson.py`; rows in
`benchmarks/grain_opt/go9/poisson.json`; gain 1,244 intensity units per count, so the weakest lit
pixels hold about 50 counts).

| objective | offset 381 (deg) | offset 302 (deg) | drop at the optimum | strict at the optimum (ground truth 0.081 / 0.052) |
|---|---|---|---|---|
| strict cost (GO1) | 0.0060 | 0.0082 | (none) | 0.082 / 0.053 |
| least squares, per-spot scale (GO1) | 0.0117 | 0.0423 | 0.6% / 3.8% | 0.103 / 0.073 |
| intensity-weighted Sinkhorn (GO1) | 0.0074 | 0.0263 | 0.1% / 0.03% | 0.102 / 0.069 |
| Poisson deviance, b = 0.5 | 0.0101 | 0.0060 | 2.3% / 1.7% | 0.092 / 0.057 |
| Poisson deviance, b = 5 | 0.0115 | 0.0146 | 2.3% / 2.1% | 0.086 / 0.058 |

Verdicts by the rule written before the runs:

- I18 is answered no: neither variant has its optimum below 0.01 deg on both grains; b = 0.5
  misses by 0.0001 deg on grain 381 (0.0101) and is the best smooth objective so far on grain 302
  (0.0060, below the strict cost's 0.0082).
- The Poisson deviance is a smooth objective that comes close to the strict cost's accuracy at
  the ground truth on both grains (worst 0.010 deg against least squares 0.042 and Sinkhorn
  0.026), but its value drops by about 2% at the optimum: a small real bias, from pixel counts the
  sharp coverage model cannot reproduce (partial-pixel intensities, the free per-spot scale).
- Next (I22): constrain the per-reflection scales by the predicted relative intensities
  (structure factor x Lorentz-polarisation) so a reflection cannot absorb misfit with its own
  scale, then the same bias test.

## GO9b: Poisson likelihood with structure-factor scales (goal G2a, idea I22)

Question: with the per-reflection scales tied to the renderer's intensity law, s_j = c F_j /
(|sin eta_j| sin 2theta_j) (one constant per grain; eta and 2 theta from each spot's observed
centroid direction), does the Poisson optimum move to the ground truth? GO1 protocol, grains 381
and 302; rows in `benchmarks/grain_opt/go9b/poisson.json` (the `pois` control reproduces GO9
exactly).

| objective | offset 381 (deg) | offset 302 (deg) | deviance at the ground truth (381 / 302) |
|---|---|---|---|
| Poisson, free scale per reflection (GO9) | 0.0101 | 0.0060 | 4.45M / 11.41M |
| Poisson, structure-factor scales | 0.0092 | 0.0303 | 6.65M / 12.92M |

Verdicts by the rule written before the runs:

- I22 is answered no: grain 381 moves inside the 0.01 deg limit (0.0092) but grain 302 moves out
  to 0.030 deg.
- The tied scales fit the observed relative intensities worse (deviance at the ground truth 1.5 x
  and 1.13 x the free-scale value): the predicted coverage mass is an area, while a rendered
  pixel's value is the reflection intensity times the number of voxel triangles whose
  rasterisation covers the pixel centre, and the centroid-based eta and 2 theta are approximate,
  so per-spot misfit is left for the orientation to absorb.
- G2a (I18, I22) closes with the Poisson deviance with free scales as the best smooth objective
  (`COSTS.md`), still short of the strict cost's rule.

## GO4e: neighbour pass after the gated seed path (goal G1b, idea I51)

Question: after GO4d's ladder + gated seed path, do a neighbour continuation (`N/30/0.3`) and a
second seed pass recover more voxels, above all rand at 2 deg (52% above 0.1 deg)? Tables
`benchmarks/grain_opt/go4e/{corr3,rand}.md`. Share above 0.1 deg (wrong, above 1 deg).

| sample, arm | 0.5 deg | 1 deg | 2 deg | 5 deg | ms / voxel |
|---|---|---|---|---|---|
| corr3, GO4d `L/366/3,G/15/0.15` | 3% (0.4%) | 5% (0.4%) | 9% (4%) | 36% (28%) | 285 |
| corr3, + `N/30/0.3` | 2% (0.4%) | 4% (0%) | 9% (4%) | 33% (27%) | 289 |
| corr3, + `N/30/0.3,G/15/0.15` | 2% (0.4%) | 4% (0%) | 9% (4%) | 33% (27%) | 394 |
| rand, GO4d `L/366/3,G/15/0.15` | 26% (2%) | 15% (7%) | 52% (47%) | 44% (42%) | 592 |
| rand, + `N/30/0.3` | 25% (2%) | 13% (7%) | 48% (43%) | 41% (39%) | 596 |
| rand, + `N/30/0.3,G/15/0.15` | 25% (2%) | 13% (7%) | 48% (43%) | 41% (39%) | 851 |

Verdicts by the rule written before the runs:

- I51 is answered no: rand at 2 deg falls from 52% to 48%, not to 35%. No level gets worse.
- The neighbour pass is a small, cheap gain everywhere (1-4 points, +1% time); it is worth
  keeping after the seed path. The second seed pass changes nothing (the same voxels are gated
  and none finds a lower strict cost) at 30-45% more time.
- rand at 2 deg keeps 43% wrong (above 1 deg) after all of this. Next (G1b, a diagnosis step):
  are those voxels gated, and is their strict cost at the ground truth lower than at their end
  (is it the search or the objective that fails)?

## GO4f: why voxels stay wrong after the G1b pipeline (diagnosis)

For the voxels still more than 1 deg from the ground truth after `L/366/3,G/15/0.15,N/30/0.3`
(`go4f_diag.py`; `benchmarks/grain_opt/go4f/`): gated by the seed path, lowered by it, and
strict cost at the ground truth against the end (medians).

| sample, level | wrong / n | gated | strict lowered by the seed path | ground truth lower than the end | strict at the ground truth / end |
|---|---|---|---|---|---|
| rand, 0.5 deg | 5 / 269 | 4 | 3 | 5 | 0.087 / 0.178 |
| rand, 1 deg | 20 / 277 | 20 | 15 | 20 | 0.083 / 0.277 |
| rand, 2 deg | 106 / 244 | 106 | 57 | 106 | 0.079 / 0.393 |
| rand, 5 deg | 98 / 253 | 98 | 81 | 98 | 0.065 / 0.405 |
| corr3, 2 deg | 9 / 244 | 9 | 8 | 9 | 0.068 / 0.398 |
| corr3, 5 deg | 68 / 253 | 68 | 56 | 68 | 0.066 / 0.399 |

- Every remaining wrong voxel (but one at rand 0.5 deg) was gated, and for every one the strict
  cost at the ground truth is far lower than at its end (about 0.07-0.09 against 0.18-0.41): the
  failures are search failures, not objective failures. The strict cost would accept the ground
  truth; the seed path's candidates (top 20 from indexing, strict cost on 5, finisher on 2) do not
  contain it, often because a voxel's attached rays are dominated by its neighbours' spots.
- The seed path lowered the strict cost of about half (rand 2 deg) to most (5 deg) of these
  voxels without reaching the ground truth: it lands on a better wrong orientation.
- Next (I54): widen the search for gated voxels the seed path leaves above the gate: more
  indexing candidates through the strict cost (top 20 all costed, finisher on 5), then the
  reconstructor's grid-search fallback (`grid_seed`) for those still above the gate.

## GO4g: wider search for the voxels the seed path misses (goal G1b, idea I54)

Question: GO4f showed the remaining wrong voxels are search failures. Does a wider search on the
voxels still above the gate recover them? Stages `W` (seed path with all 20 indexing candidates
costed and the finisher on 5) and `X` (the reconstructor's per-voxel grid search), each kept per
voxel only where it lowers the strict cost. `X` cost 16 s per searched voxel on the smoke grain,
so by the plan it ran on rand only. A first `W`/`X` smoke run was rejected by the stage gate (the
gate compared NaN objectives); fixed before the judged runs. Tables
`benchmarks/grain_opt/go4g/{corr3,rand}.md`. Share above 0.1 deg (wrong, above 1 deg).

| sample, pipeline | 0.5 deg | 1 deg | 2 deg | 5 deg | searched / improved | s / voxel |
|---|---|---|---|---|---|---|
| corr3, GO4e `L/366/3,G/15/0.15,N/30/0.3` | 2% (0.4%) | 4% (0%) | 9% (4%) | 33% (27%) | | 0.29 |
| corr3, + `W/20/0.15` | 2% (0.4%) | 3% (0%) | 3% (1%) | 21% (13%) | W 153 / 75 | 0.39 |
| rand, GO4e | 25% (2%) | 13% (7%) | 48% (43%) | 41% (39%) | | 0.60 |
| rand, + `W/20/0.15` | 25% (2%) | 5% (2%) | 21% (17%) | 14% (9%) | W 359 / 219 | 0.90 |
| rand, + `W/20/0.15,X/1/0.15` | 17% (0.4%) | 2% (1%) | 7% (5%) | 6% (3%) | X 210 / 155 | 4.0 |

Verdicts by the rule written before the runs:

- I54 is answered yes for `+W`: rand at 2 deg 48% -> 21% (limit 32%), at 5 deg 41% -> 14% (limit
  27%); corr3 at 5 deg 33% -> 21% (limit 22%); no level worse. `+W,X` passes on rand by a wide
  margin but was not run on corr3 (cost rule), so the rule's preferred arm is `+W`.
- Wrong voxels (above 1 deg) over all levels: rand 229 -> 75 (`W`) -> 25 (`W,X`) of 1,043;
  corr3 78 -> 37. The grid search finds almost all of what the indexing misses, at about 16 s
  per searched voxel; with it, rand (random voxel-to-voxel misorientation up to 5 deg RMS) is
  recovered to 6-7% above 0.1 deg at 2 and 5 deg, from 83% after the ladder alone (GO4c).
- What remains is precision more than search on rand at 0.5 deg (17%, nearly all between 0.1
  and 1 deg): a final per-voxel strict polish is the next lever. Next (later, G1b): `X` on corr3
  and a per-voxel polish (`V`) after `W,X`.

## GO4h: grid search on corr3 and a final per-voxel polish (goal G1b, idea I55)

Question: does a final per-voxel strict polish (`V/61/0.02`) after the full search fix rand's
near misses at 0.5 deg, and what does the grid search `X` do on corr3? Pipeline
`L/366/3,G/15/0.15,N/30/0.3,W/20/0.15,X/1/0.15,V/61/0.02`; tables
`benchmarks/grain_opt/go4h/{corr3,rand}.md` (runs finished before the machine restart of
2026-11-17; `X` on rand reproduces GO4g exactly). Share above 0.1 deg (wrong, above 1 deg);
median deg.

| sample, after | 0.5 deg | 1 deg | 2 deg | 5 deg | s / voxel |
|---|---|---|---|---|---|
| corr3, `W` (GO4g) | 2% (0.4%), 0.012 | 3% (0%), 0.009 | 3% (1%), 0.009 | 21% (13%), 0.016 | 0.39 |
| corr3, `W,X` | 1.5% (0%), 0.012 | 2% (0%), 0.009 | 1% (0%), 0.008 | 10% (4%), 0.011 | 1.98 |
| corr3, `W,X,V` | 1% (0%), 0.009 | 1% (0%), 0.006 | 0.4% (0%), 0.006 | 7% (4%), 0.006 | 1.98 |
| rand, `W,X` (GO4g) | 17% (0.4%), 0.018 | 2% (1%), 0.010 | 7% (5%), 0.007 | 6% (3%), 0.013 | 4.0 |
| rand, `W,X,V` | 16% (0.4%), 0.008 | 1.4% (1%), 0.006 | 5% (5%), 0.005 | 4% (3%), 0.007 | 4.1 |

Verdicts by the rule written before the runs:

- I55 is answered no: rand at 0.5 deg goes from 17% to 16% (limit 11%). No level gets worse.
- The polish halves the median error everywhere (to 0.005-0.009 deg) at about 6 s for the whole
  sample, but leaves rand's 0.5 deg near misses (0.1-0.75 deg) in place: their strict cost is
  below the 0.15 gate, so no search stage touches them, and a 0.02 deg polish cannot move 0.1-0.7
  deg.
- corr3 with the grid search: 5 deg from 21% to 10%, wrong voxels 37 -> 10; with the polish
  corr3 ends at 1 / 1 / 0.4 / 7% above 0.1 deg. With GO4c's ladder alone it was 3 / 9 / 31 / 72%.
- Follow-up (later, G1b): a per-level gate (the 0.5 deg near misses need a wider local step, or a
  gate on misfit to neighbours, not on the strict cost alone).

## GO7d: split merged spots before taking centroids (goal G2b, idea I53)

Question: on the 1 px PSF images, does assigning each observed lit pixel to the grain or a
neighbour (nearest start-predicted point, grain against the fixed shell prediction) before taking
the observed centroid restore the fast strain fit? S36, `GP4@P --cen_split`; rows in
`benchmarks/grain_opt/go7d/` (PSF) and `go7d_clean/` (clean control).

| arm | 1e-4 | 5e-4 | 1e-3 | orientation | fallbacks | s / grain |
|---|---|---|---|---|---|---|
| clean (GO5l) | 1.3e-4 | 1.3e-4 | 2.1e-4 | 0.009 / 0.015 / 0.020 | 9 | 25 |
| clean, `--cen_split` | 1.3e-4 | 1.4e-4 | 2.1e-4 | 0.012 / 0.019 / 0.015 | 5 | 24 |
| 1 px PSF (GO7b) | 2.3e-4 | 2.7e-4 | 5.7e-4 | 0.029 / 0.039 / 0.061 | 11 | 18 |
| 1 px PSF, `--cen_split` | 2.4e-4 | 2.7e-4 | 4.9e-4 | 0.029 / 0.035 / 0.061 | 11 | 18 |

Verdicts by the rule written before the runs:

- I53 is answered no: on the PSF images the strain medians change by 0.95-1.16 x (limits not
  met; no 1.5 x gain). The clean control loses nothing (within 1.25 x; fewer fallbacks, 5 against
  9).
- The real cause, found while closing: the region build's area-ratio filter (observed lit pixels
  / pixels touched by the start prediction, kept in 0.8-1.6) rejects the blur-enlarged spots. With
  a 1 px PSF a grain keeps a median 58% of its clean spots (45 against 84; as few as 34%); with 0.5
  px, 89%. Fewer spots, more fallbacks (11 against 9), worse strain. Neither merging (GO7d) nor a
  flat coverage term (GO7c) was the main effect.
- Next (I56): widen the ratio filter for blurred images (`region_kw` ratio upper bound 2.5, no
  library change) on the 1 px PSF images, with a clean control.

## GO7e: widen the region area-ratio filter (goal G2b, idea I56)

Question: with the region build's area-ratio filter widened from 0.8-1.6 to 0.8-2.5
(`--ratio_hi 2.5`), does the fast strain fit recover on the 1 px PSF images, and what happens on
clean images? S36, `GP4@P`; rows in `benchmarks/grain_opt/go7e/` (PSF) and `go7e_clean/`.

| arm | spots / GO5l clean | 1e-4 | 5e-4 | 1e-3 | orientation | s / grain |
|---|---|---|---|---|---|---|
| clean, filter 1.6 (GO5l) | 1.00 | 1.3e-4 | 1.3e-4 | 2.1e-4 | 0.009 / 0.015 / 0.020 | 25 |
| clean, filter 2.5 | 1.30 | 1.3e-4 | 1.2e-4 | 1.4e-4 | 0.010 / 0.013 / 0.011 | 35 |
| 1 px PSF, filter 1.6 (GO7b) | 0.58 | 2.3e-4 | 2.7e-4 | 5.7e-4 | 0.029 / 0.039 / 0.061 | 18 |
| 1 px PSF, filter 2.5 | 1.09 | 2.1e-4 | 2.1e-4 | 3.5e-4 | 0.027 / 0.027 / 0.047 | 29 |

Verdicts by the rule written before the runs:

- I56 is answered no: on the PSF images the strain medians improve (gain 1.27 x at 5e-4, 1.62 x
  at 1e-3) but miss the no-harm limits (2.0e-4 and 3.2e-4) and the 1.5 x gain at 5e-4;
  orientation at 1e-3 is 0.047 deg. The spot count is restored (1.09 x clean).
- The clean control is the finding: with the wider filter the clean fast fit improves at 1e-3
  from 2.1e-4 to 1.4e-4 strain and from 0.020 to 0.011 deg orientation (5e-4: 1.3e-4 -> 1.2e-4,
  0.015 -> 0.013 deg), with 30% more spots and 40% more time (35 s per grain, still below the CMA
  hybrid's 53 s). The 1.6 bound was also dropping useful spots on clean images, likely spots whose
  observed area differs from the unstrained start prediction because of the strain itself.
- The remaining PSF damage (2.1e-4 at 5e-4, orientation 0.03-0.05 deg) needs a forward model with
  the PSF (library work for the owner to approve) or more images-side work; G2b stops here.
- Recommendation for the fast strain fit: `--ratio_hi 2.5` by default (clean 1.3 / 1.2 / 1.4e-4,
  35 s per grain); confirm on the noisy images before relying on it.

## GO5p: the wider ratio filter on noisy images (idea I57)

Question: does the region ratio filter widened to 2.5 (GO7e, better on clean images) also hold on
the noisy `strainrn` images, so it can become the fast fit's default? S36, `GP4@P`; rows in
`benchmarks/grain_opt/go5p/`.

| arm (noisy images) | 1e-4 | 5e-4 | 1e-3 | orientation | s / grain |
|---|---|---|---|---|---|
| CMA hybrid (GO5m, reference) | 2.5e-4 | 1.9e-4 | 2.8e-4 | 0.012 / 0.012 / 0.028 | 51 |
| fast fit, filter 1.6 (GO5m) | 1.8e-4 | 2.7e-4 | 3.7e-4 | 0.014 / 0.017 / 0.025 | 27 |
| fast fit, filter 2.5 | 1.6e-4 | 2.6e-4 | 3.1e-4 | 0.013 / 0.017 / 0.011 | 38 |
| fast fit, filter 2.5, `intclean` | 2.0e-4 | 2.5e-4 | 2.7e-4 | 0.015 / 0.011 / 0.013 | 40 |

Verdicts by the rule written before the runs:

- I57 is answered no: both arms are better than GO5m's fast fit at 5e-4 and 1e-3 and within 1.25 x
  the noisy CMA hybrid at 1e-3, but miss it at 5e-4 (2.6e-4 and 2.5e-4 against 2.4e-4).
  Orientation is within 0.017 deg everywhere (better than GO5m at 1e-3: 0.011 against 0.025).
- So the wider filter is the better choice on clean images (GO7e) and helps on noisy images too,
  but on noisy images the fast fit still does not match the CMA hybrid at 5e-4, at 75% of its
  time. The fast fit's documented default stays as GO5l, with `--ratio_hi 2.5` recommended for
  clean data.
- The strain line closes here (program balance).

## GO4i: a lower gate for rand's near misses (goal G1b, idea I58)

Question: does a second wide-search pass (`W`, `X`) at gate 0.11, after GO4h's pipeline, recover
rand's 0.5 deg near misses (strict cost 0.11-0.20, below the 0.15 gate)? Pipeline
`L/366/3,G/15/0.15,N/30/0.3,W/20/0.15,X/1/0.15,W/20/0.11,X/1/0.11,V/61/0.02` on rand; table
`benchmarks/grain_opt/go4i/rand.md`.

| rand | 0.5 deg | 1 deg | 2 deg | 5 deg | wrong (all) | s / voxel |
|---|---|---|---|---|---|---|
| GO4h (gate 0.15 only) | 16% | 1.4% | 5% | 4% | 25 | 4.1 |
| + second pass at 0.11 | 4.5% | 0.4% | 3.3% | 3.2% | 18 | 8.1 |

Verdicts by the rule written before the runs:

- I58 is answered yes: rand at 0.5 deg falls from 16% to 4.5% (limit 11%), and every other level
  improves; median error 0.005-0.006 deg at every level, 95th percentile 0.04 deg (was 0.42).
- The second pass searched 296 voxels with `W` (45 improved) and 270 with `X` (168 improved): the
  grid search does the work.
- Over budget: the run took 140 minutes against the plan's 2-hour stop. The session that watched
  it ended in a restart before the 2-hour mark, so it was not stopped; the result is reported as
  run, and the overrun noted.
- G1b now stands at rand 4.5 / 0.4 / 3.3 / 3.2% and corr3 1 / 1 / 0.4 / 7% above 0.1 deg at 0.5 /
  1 / 2 / 5 deg RMS voxel-to-voxel misorientation, at 2-8 s per voxel (GO4c's ladder alone: rand
  48 / 53 / 77 / 83%, corr3 3 / 9 / 31 / 72%).

## GO5q: intragranular strain gradient (goal G1c, idea I59; owner request)

Question: on `straing` (the `strain` grain means plus an in-plane linear strain gradient per grain,
gradient RMS over the grain equal to its level), do the fits recover the mean and the gradient?
`fit_strain.py --sample straing --grad` (6 + 12 strain numbers per grain), `hyb` lambda 8, ratio
filter 2.5, S36; rows in `benchmarks/grain_opt/go5q/`. Relative gradient error = RMS over voxels of
the fitted minus ground-truth intragranular part, over the ground-truth part (1 = nothing
explained, 0 = all).

| arm | mean error 1e-4 / 5e-4 / 1e-3 | relative gradient error | voxel strain error | orientation (deg) | s / grain |
|---|---|---|---|---|---|
| uniform strain at the ground-truth orientation | 1.9 / 3.6 / 10.4e-4 | 1 / 1 / 1 | 1.7 / 5.4 / 10.3e-4 | 0 | 23 |
| + gradient at the ground-truth orientation | 1.4 / 1.3 / 1.7e-4 | 4.65 / 0.77 / 0.45 | 4.9 / 4.0 / 4.8e-4 | 0 | 63 |
| CMA hybrid + gradient, painted start | 2.0 / 1.7 / 2.4e-4 | 4.12 / 0.97 / 0.58 | 4.8 / 5.3 / 6.2e-4 | 0.011 / 0.013 / 0.023 | 123 |
| fast fit (filter 2.5) + gradient leg | 1.9 / 1.4 / 2.9e-4 | 2.71 / 0.79 / 0.66 | 3.4 / 4.4 / 7.6e-4 | 0.010 / 0.019 / 0.036 | 76 |

Verdicts by the rule written before the runs:

- I59 is answered no: no arm recovers the gradient (relative error at most 0.5) at both 5e-4 and
  1e-3; only the strain-only fit at the ground-truth orientation reaches it at 1e-3 (0.45). The
  mean-strain limits (1.5 x the uniform-sample results) are met by the CMA hybrid (1.7e-4, 2.4e-4
  against 2.2e-4, 3.2e-4) and by the fast fit at 5e-4, not at 1e-3 (2.9e-4 against 2.1e-4).
- Fitting the gradient is what makes the mean right: a uniform fit on this sample misses the mean
  by 3.6e-4 and 1.0e-3 at 5e-4 and 1e-3 (the gradient biases it), the gradient fits by 1.3-2.9e-4.
- The gradient itself is weakly observed: a linear gradient across a grain changes each spot's
  shape, not its centroid (the voxel shifts average out), and the coverage term resolves shape
  changes of a fraction of a pixel poorly. At 1e-4 the 12 gradient numbers fit noise (relative
  error 2.7-4.7, about 4e-4 per voxel). The voxel strain error floor is about 4-5e-4 in every
  gradient arm.
- Next (I60, queued): a Gauss-Newton fit with per-spot second-moment (spot spread) residuals as
  well as centroids, since the gradient acts on the spread; judged against these arms.

## GO5r: Gauss-Newton gradient fit with spot-spread residuals (goal G1c, idea I60)

Question: does Gauss-Newton on per-spot centroids plus second moments (row variance and
column-row covariance; the column variance was dropped before the runs, biased by -4.2 px^2 at
the ground truth from the in-frame omega smear) recover the linear strain gradient on `straing`
better and faster than GO5q's CMA fits? Rows in `benchmarks/grain_opt/go5r/`. When every
Gauss-Newton stage was skipped (9 of 36 grains from the painted start), GO5q's CMA gradient leg ran
instead (amendment before the runs).

| arm | mean error 1e-4 / 5e-4 / 1e-3 | relative gradient error | orientation | s / grain |
|---|---|---|---|---|
| CMA, strain only at the ground-truth orientation (GO5q) | 1.4 / 1.3 / 1.7e-4 | 4.65 / 0.77 / 0.45 | 0 | 63 |
| Gauss-Newton + moments, at the ground-truth orientation | 0.95 / 1.4 / 1.1e-4 | 12.3 / 1.28 / 0.81 | 0 | 15 |
| CMA hybrid + gradient, painted start (GO5q) | 2.0 / 1.7 / 2.4e-4 | 4.12 / 0.97 / 0.58 | 0.011 / 0.013 / 0.023 | 123 |
| fast fit + CMA gradient leg (GO5q) | 1.9 / 1.4 / 2.8e-4 | 2.71 / 0.79 / 0.66 | 0.010 / 0.019 / 0.036 | 76 |
| fast fit + Gauss-Newton + moments | 1.3 / 1.4 / 2.8e-4 | 12.7 / 1.38 / 0.74 | 0.010 / 0.020 / 0.038 | 50 |

Verdicts by the rule written before the runs:

- I60 is answered no: the gradient error is worse than the CMA hybrid's at 5e-4 and 1e-3 (1.38 /
  0.74 against 0.97 / 0.58), not 0.5 or better. The spread residuals pull the 12 gradient numbers
  toward noise (relative error above 1 means a gradient larger than the ground truth's in the
  wrong direction); the single smoke grain (237: 0.26) was not representative.
- Where Gauss-Newton helps is the mean and the time: at the ground-truth orientation it gives the
  best mean strain so far (0.95 / 1.4 / 1.1e-4) in 15 s per grain against CMA's 63 s, and from
  the painted start the same mean as the fast fit at 50 s against 76 s.
- The intragranular gradient stays weakly observed by every cost tried (coverage, centroids,
  spreads): the spot-shape signal of a linear gradient across a grain is a fraction of a pixel
  and comparable to the model-image mismatch in spot shape (the omega smear). A smoother
  parameterisation (I62) or a prior (I64) can only help if the signal is there; the non-linear
  sample (I61) tests the case where the field is richer.

## GO5s: non-linear strain-gradient sample and the linear fit on it (goal G1c, idea I61)

Sample `straingq`: the `strain` grain means plus, per grain, an in-plane linear term, a quadratic
term and a smooth Gaussian random field (correlation length 3 voxel spacings), each of tensor RMS
level / sqrt(3) over the grain; amended before any run so that the quadratic and random parts are
orthogonal to the linear basis per grain (without that the oracle linear share was 0.60). Oracle
linear share (least-squares linear projection of the ground-truth intragranular field, median over
grains): 0.33, so a perfect linear fit has relative error sqrt(0.67) = 0.82. Rendered clean (9.5M
lit pixels). Linear fits (6 + 12) as GO5q; rows in `benchmarks/grain_opt/go5s/`.

| arm on `straingq` | mean error 1e-4 / 5e-4 / 1e-3 | relative error, full intragranular field | voxel error | s / grain |
|---|---|---|---|---|
| linear, strain only at the ground-truth orientation | 1.5 / 1.7 / 1.7e-4 | 3.73 / 1.05 / 0.98 | 4.0 / 5.5 / 10.3e-4 | 60 |
| fast fit + linear CMA leg, painted start | 1.8 / 1.9 / 2.5e-4 | 2.54 / 1.10 / 0.98 | 2.9 / 5.9 / 10.4e-4 | 72 |
| oracle linear projection (floor) | | 0.82 | | |

Verdicts by the rule written before the runs:

- I61 is answered (sample delivered): it renders, and the oracle linear share is 0.33 (limit 0.5).
- The linear fits hold their mean on the non-linear sample (within 1.5 x of GO5q at 5e-4 and 1e-3:
  1.7e-4 and 1.7e-4 at the ground-truth orientation, 1.9e-4 and 2.5e-4 from the painted start).
- Against the full intragranular field they explain almost nothing (relative error 0.98-1.10 at
  5e-4 and 1e-3, against the linear floor of 0.82): the linear part they can see is a third of the
  field, and they recover only part of that (GO5q).
- Next (I62): the quadratic field (6 + 12 + 18 numbers) on `straingq`, against these rows; the
  oracle floor for linear + quadratic is about sqrt(1/3) = 0.58.

## GO5t: quadratic per-grain strain field (goal G1c, idea I62)

Question: on `straingq`, does a quadratic per-grain strain field (mean 6 + linear 12 + quadratic 18
= 36 numbers) recover more of the intragranular field than the linear field (GO5s)? Rerun from
scratch after the 2026-11-19 shutdown checkpoint (the partial first run was not judged). Rows in
`benchmarks/grain_opt/go5t/`.

| arm on `straingq` | relative intragranular error 1e-4 / 5e-4 / 1e-3 | mean error 1e-4 / 5e-4 / 1e-3 | voxel error at 1e-3 | s / grain |
|---|---|---|---|---|
| linear, ground-truth orientation (GO5s) | 3.73 / 1.05 / 0.98 | 1.5 / 1.7 / 1.7e-4 | 10.3e-4 | 60 |
| quadratic, ground-truth orientation | 4.90 / 1.24 / 0.72 | 1.6 / 1.2 / 2.1e-4 | 7.7e-4 | 97 |
| fast fit + linear leg (GO5s) | 2.54 / 1.10 / 0.98 | 1.8 / 1.9 / 2.5e-4 | 10.4e-4 | 72 |
| fast fit + quadratic leg | 3.89 / 1.09 / 0.74 | 1.7 / 2.0 / 2.1e-4 | 7.7e-4 | 96 |
| oracle floors: linear 0.82, linear + quadratic 0.58 | | | | |

Verdicts by the rule written before the runs:

- I62 is answered no: the quadratic field is at least 0.1 better than the linear field at 1e-3 (0.72
  and 0.74 against 0.98; below the linear oracle floor of 0.82) but not at 5e-4 (1.24 and 1.09
  against 1.05 and 1.10). Mean strains stay within 1.5 x of the linear arms.
- I63 (coordinate network) is not run: the rule runs it only if the quadratic arm passes.
- Reading: at 1e-3 the intragranular field is large enough for the extra 18 numbers to capture
  real structure (voxel error 1.03e-3 -> 0.77e-3); at 5e-4 and below the added freedom mostly fits
  noise (relative error above 1). A field model's usable order depends on the strain level, which
  is what a prior (I64) is meant to handle: let the data decide how much intragranular structure to
  admit.

## GO5u: GP prior on the intragranular strain residual (goal G1c, idea I64)

Question: on `straingq`, does a GP prior on a smooth residual (5 eigenvectors per strain component
of a squared-exponential kernel over the grain, orthogonal to mean + linear; prior weight B on the
coefficients scaled by the eigenvalues; 48 strain numbers per grain) beat the quadratic field at
5e-4 without losing at 1e-3? Rows in `benchmarks/grain_opt/go5u/`.

| arm on `straingq` | relative intragranular error 1e-4 / 5e-4 / 1e-3 | voxel error 1e-4 / 5e-4 / 1e-3 | mean error 5e-4 / 1e-3 | s / grain |
|---|---|---|---|---|
| quadratic, ground-truth orientation (GO5t) | 4.90 / 1.24 / 0.72 | 5.2 / 6.2 / 7.7e-4 | 1.2 / 2.1e-4 | 97 |
| GP basis, no prior (B = 0), ground-truth orientation | 1.00 / 2.06 / 0.91 | 1.4 / 10.2 / 8.7e-4 | 1.7 / 1.7e-4 | 101 |
| GP basis + prior (B = 0.05), ground-truth orientation | 1.00 / 1.73 / 0.78 | 1.4 / 8.3 / 7.5e-4 | 1.6 / 1.7e-4 | 101 |
| fast fit + quadratic leg (GO5t) | 3.89 / 1.09 / 0.74 | 3.9 / 5.8 / 7.7e-4 | 2.0 / 2.1e-4 | 96 |
| fast fit + GP leg + prior (B = 0.05) | 1.00 / 1.23 / 0.84 | 1.7 / 6.2 / 8.7e-4 | 2.0 / 2.5e-4 | 101 |

Verdicts by the rule written before the runs:

- I64 is answered no: the prior arms are worse than the quadratic at 5e-4 (1.73 and 1.23 against
  1.24 and 1.09) and at 1e-3 more than 0.05 worse (0.78 and 0.84 against 0.72 and 0.74). Mean
  strains are within the limits.
- The prior itself helps the same basis (B = 0 -> 0.05: 2.06 -> 1.73 at 5e-4, 0.91 -> 0.78 at
  1e-3): it pulls the coefficients toward a smooth, small residual. At 1e-4 it keeps the fit
  nearly flat (relative error 1.00, voxel error 1.4-1.7e-4 against 3.9-5.2e-4 for the quadratic),
  which is where an unregularised field fits noise.
- The GP eigenvector basis is a worse basis than the quadratic for this sample at 5e-4 and 1e-3,
  and CMA over 48 numbers is a weak optimiser for it. The I59-I64 queue ends here: across linear,
  quadratic and GP-prior fields, the intragranular strain field is recovered only at the 1e-3 level
  (best relative error 0.72) and only partly; the observable signal (spot shape below a pixel) is
  the limit, not the field model.
