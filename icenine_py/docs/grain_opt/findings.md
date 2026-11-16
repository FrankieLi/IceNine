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
