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
