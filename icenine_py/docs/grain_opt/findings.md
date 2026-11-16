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
