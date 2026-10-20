# Findings: smoothing the cost for ranking (Task E, 2026-10-15)

Run record: `benchmarks/eval_reduction/E/RUN_LOG.md`. Library: `icenine/smooth_cost.py`, config key `SmoothCostRank`
(opt-in, off by default). The full technical record, with every table, is MIGRATION_HISTORY "Task E". This doc is the
plain-English reading of it.

## What we asked

The search scores an orientation with the **strict cost**: one minus the share of predicted spot pixels that land on
lit pixels. Task B showed that this cost changes in steps. It jumps when a spot edge crosses a pixel, when a spot moves
into the next omega frame, and when a peak enters or leaves the list.

Steps matter most at the coarse levels (0, 1 and 2). At each of them, every candidate gets a short Monte Carlo polish
(the **quick MC**). Then the list is sorted by cost and only the best quarter is kept (the **quarter cut**). A
candidate 1 to 3 degrees from the truth often has no spot overlap with the truth's spots at all, so the strict cost
scores it near 1, as bad as a random orientation. Yet it is close enough that the next level's finer grid would still
reach the truth. If the cut throws it out, the voxel can end up wrong.

So we asked: can a smoother cost order the candidates better at the quarter cut? The strict cost would still decide
everything else: the quick MC moves, level 3, FindOptimal and the final answer.

## What we found, in short

1. **On the E0 test set, a smooth ranking fixes most wrong voxels and breaks none.** E0 is our set of 200 single
   voxels. Each one has its own simulated images, holding only that voxel's spots. With the smooth ranking, wrong
   voxels fell from 41 to 7 or 5 (lenient rule) and from 56 to 9 (C++ rule). No voxel that was right became wrong.

2. **This is not yet shown on the full sample.** The full-sample images also hold the neighbouring grains. Our six
   full-sample test voxels are all right with the strict cost, so they can show harm but not a gain. They showed no
   harm. Whether the gain holds on full-sample voxels that are wrong today is the next thing to test.

3. **Each smoothing removes one kind of step.** Credit for near misses removes the pixel steps. Blending neighbouring
   frames removes the frame steps. Soft peak weights remove the steps where a peak enters or leaves the list. All three
   together leave 19 steps on 120 rays, against 28,290 for the strict cost.

4. **Smoothing alone can add false dips.** Credit for near misses, used on its own, adds dips away from the truth on
   the full-sample images but not on the E0 images. That is consistent with it rewarding lit pixels of *other* grains
   near a predicted spot; we did not test that directly. The combinations that include the soft peak weights and the
   observed-ray term have the fewest dips, fewer than the strict cost.

5. **The price is time.** As shipped, the cheaper carried combination adds about 33% to the time of a full-sample
   voxel and the richer one about 116% (a model, close to a contended measurement of one voxel). Much of that is a
   16.6 s set-up step that the library repeats on every voxel. Done once per image set, it would be about 20% and
   103%; that cache is not built yet (Task G). On E0 voxels, measured on a quiet machine, the ranking more than
   doubles the time.

6. **The observed-ray term agrees with Task D.** It is not sharper than the pixel cost. It helps only inside a
   combination, but it has a gradient, which Task F can use.

7. **No default changed.** With `SmoothCostRank` off, the search is bit-identical to before. A test compares it with
   a search recorded on the commit before Task E. Also, the 200 E0 baseline runs of this task equal Task C's baseline
   runs on all 200 voxels (same error, same evaluation count).

## Terms

* **Strict cost.** The pixel-overlap cost the search uses today (above).
* **Arm.** One way of smoothing the cost. Arms are used only to *order* candidates at the quarter cut of levels 0-2.
* **Truth basin.** The candidates close enough to the truth that the next level's grid contains it: within the
  level's diameter divided by 1.5. With a 5 deg diameter at level 0, that is 5 / 1.5 = 3.33 deg at level 0, then
  2.22 deg at level 1 and 1.48 deg at level 2 (worked out from the diameters, not read from a results file).
* **Kept, fixed, broken.** The basin is *kept* when at least one basin candidate survives the quarter cut. A voxel or
  level is *fixed* when the arm keeps it and the strict cost does not. It is *broken* in the reverse case.
* **Wrong.** A final orientation more than 1 deg from the truth.
* **Lenient rule, C++ rule.** Two versions of the strict cost that differ in how edge pixels count (Task A). Both
  were tested.
* **Ray.** A straight walk out from the true orientation in one direction, 0 to 5 deg (Task B's method).
* **Isolated step.** A jump between two neighbouring points of a ray, much larger than the change on either side.
* **Dip of depth d.** A local minimum along a ray, away from the truth, at least d below the rise around it. A dip can
  trap a search.
* **Half-way rise.** How far from the truth the cost first gets half-way from its value at the truth to 1. Smaller
  means a sharper cost.

## The arms

Every arm predicts the same spots as the strict cost. They differ only in how a predicted spot is scored. Each spot is
checked at four points: its three corners and its centre.

* **S, the control.** A point scores 1 on a lit pixel, else 0. It tests whether the four-point sampling by itself
  changes anything. (It does not.)
* **E1, credit for near misses.** A point's score falls off smoothly with its distance to the nearest lit pixel,
  instead of jumping from hit to miss.
* **E2, frame blending.** A spot near a frame edge also takes credit from the neighbouring frame, weighted by how close
  it is to that edge.
* **E3, soft peak weights.** A peak's weight fades to zero near a detector edge, the eta limit or the ends of the omega
  scan. A peak entering or leaving the list then no longer causes a jump. E3 also drops the fixed half credit for spots
  that miss the far detector, so the cost at the truth becomes 0.
* **E4, observed rays.** Task D matches spots across the two detectors into observed rays. E4 measures the angle
  between each predicted peak's direction and the nearest observed ray through the voxel. It is smooth, with a
  gradient.
* **Combinations** are written with "+". When E4 is combined with an image arm, the cost is the mean of the two.

The arm widths follow the level, from the quick MC's step size. At level 0 they are 16.5 pixels for E1, 1.67 frames
for E2 and 0.42 deg for E4. The widths were set once by rule of thumb, not tuned.

## Measurement 1: the shape of each cost along the rays

We walked each cost along Task B's rays: 6 full-sample voxels x 20 directions, and 30 E0 voxels x 20 directions.

One ray shows best what happens. This is Task B's ray plot (voxel 15901, ray 0, out to 0.5 deg), with each arm in red
on its own panel against the strict cost in black. Read: a staircase is a cost that jumps; a smooth curve is one that
does not.

![Each arm against the strict cost](../benchmarks/eval_reduction/E/figures/ray_each_arm.png)

* The strict cost (black) is a staircase all the way out.
* S (the control) and E3 alone follow the same staircase: four-point sampling and soft peak weights do not remove the
  pixel and frame steps.
* E1 alone (distance credit) has fewer but larger jumps, with flat stretches between them. These are the frame steps it
  leaves behind.
* E2 alone (frame blending) is smoother but still has small pixel steps.
* E1 + E2 is smooth apart from one jump near 0.46 deg, where a peak leaves the list. E1 + E2 + E3 is smooth: E3
  removes that last kind of step.
* E4 (observed rays) is smooth, but it starts high at the truth (about 0.39 at tau 0.05 deg). That is not because
  predicted peaks lack an observed ray: almost all of them have one nearby. The nearest ray is just not close enough,
  compared with tau (see the table after the step table below).
* E3 + E4 keeps small steps from its image part, but they are smaller than the strict cost's.

The same ray with every arm on shared axes, at three ranges (0 to 0.1, 0 to 0.5 and 0 to 5 deg). Read the 5-degree
panel for the far field: the strict cost reaches 1 by about 2 deg and stays flat, while E1 at 12 px and E1 + E2 level
off well below 1. Far from the truth they separate orientations less well.

![The strict cost and every arm along Task B's ray](../benchmarks/eval_reduction/E/figures/ray_strict_vs_arms.png)

The tables below count the same thing over all 120 full-sample rays.

Steps along the full-sample rays. Read across a row: a good arm has few steps per degree and a low share of its change
in steps. The strict cost changes almost entirely in steps.

How the columns are measured (`scripts/eval_reduction/smooth_landscape_summary.py`, `isolated_steps`):

* Each ray is sampled every 0.0005 deg, so 2,000 points per degree.
* A change between two neighbouring points is a **step** when it is above 1e-6 and more than three times the change on
  both neighbouring intervals. A sudden jump counts. A steady slope never counts, however steep it is.
* **Steps per deg** is the number of steps in a distance band from the truth, divided by the degrees walked in that
  band, summed over all 120 rays. Example: 113.8 per deg under 0.05 deg means about 6 jumps in the first 0.05 deg of a
  ray (113.8 x 0.05), roughly one every 18 sample intervals.
* **Variation in steps** is the share of the cost's total change (the sum of the size of every change) that happens in
  steps. The step count treats a jump of 1e-5 the same as one of 0.05; this column weights each jump by its size. So
  read the two together. E1 at 3 px has fewer steps than the strict cost (49 against 114 per deg), but 97.7% of its
  change under 0.05 deg still comes in jumps.
* **A flat cost has no steps.** Beyond about 2 deg the strict cost sits at 1 (no overlap with the truth's spots) and
  does not move. That is part of why the C++ rule has only 25 steps per deg from 0.5 to 5 deg: it is flat there, not
  smooth. The lenient rule has 108.
* Two jumps on back-to-back intervals are compared with each other, so neither counts. At this spacing that should be
  rare, but it was not checked.

<!-- table:evalE_land_steps_full -->
| cost | cost at the truth (median) | steps per deg, under 0.05 deg | steps per deg, 0.05 to 0.5 deg | steps per deg, 0.5 to 5 deg | variation in steps, under 0.05 deg | variation in steps, 0.05 to 0.5 deg | variation in steps, 0.5 to 5 deg |
|---|---|---|---|---|---|---|---|
| strict cost, C++ rule | 0.0703 | 113.8 | 261.5 | 25.0 | 88.2% | 84.2% | 94.0% |
| strict cost, lenient rule | 0.0703 | 116.3 | 323.6 | 108.2 | 87.5% | 78.4% | 89.6% |
| S: sampled binary (control) | 0.0717 | 96.8 | 239.0 | 46.9 | 87.6% | 75.1% | 93.5% |
| E1: distance kernel 3 px | 0.0717 | 49.3 | 45.1 | 30.3 | 97.7% | 80.5% | 73.5% |
| E1: distance kernel 12 px | 0.0694 | 49.3 | 49.5 | 47.7 | 97.5% | 93.6% | 89.4% |
| E2: omega blending, 1 frame | 0.0717 | 55.2 | 236.3 | 89.8 | 54.9% | 52.6% | 65.4% |
| E3: soft peak weights | 0.0000 | 96.5 | 241.2 | 45.4 | 87.3% | 75.0% | 92.8% |
| E4: ray residual, tau 0.05 deg | 0.4136 | 0.2 | 1.4 | 8.0 | 0.0% | 1.7% | 2.6% |
| E4: ray residual, tau 0.25 deg | 0.2930 | 4.5 | 4.3 | 13.2 | 0.3% | 6.5% | 9.7% |
| E1 + E2 | 0.0717 | 0.0 | 2.2 | 2.0 | 0.0% | 7.8% | 4.5% |
| E1 + E2 + E3 | 0.0000 | 0.0 | 0.0 | 0.0 | 0.0% | 0.0% | 0.1% |
| E2 + E3 | 0.0000 | 55.8 | 237.9 | 88.1 | 56.3% | 52.9% | 65.5% |
| E3 + E4 (0.25 deg) | 0.1464 | 96.3 | 239.7 | 54.8 | 84.2% | 63.0% | 55.5% |
<!-- /table:evalE_land_steps_full -->

* E1 alone removes the pixel steps but leaves the frame steps. E2 alone does the reverse.
* E1 + E2 together leave only the steps where a peak enters or leaves the list, and E3 removes those.
* E4 has almost no steps. But its cost at the truth is high: 0.41 at tau 0.05 deg on the full-sample rays.

Why E4 is high at the truth. E4 scores each predicted peak by its distance to the nearest observed ray within one
frame. A peak with no such ray scores 1, the worst. The table splits E4's cost at the truth into the part from those
peaks and the part from peaks that have a ray (`scripts/eval_reduction/smooth_e4_truth.py`). Read: almost no peak lacks
a ray. On the full sample, peaks with no ray give at most 5% of any voxel's cost, so nearly all of it comes from peaks
whose nearest ray is a few hundredths of a degree away. On E0 it varies more: 9 of the 30 voxels have at least one
peak with no ray, and on one voxel those peaks give 30% of the cost. The columns
are medians over voxels, so the two parts need not add up to the cost exactly.

<!-- table:evalE_e4_truth -->
| voxel set | n | E4 cost at the truth, tau 0.05 deg (median) | part from peaks with no observed ray nearby (median) | part from peaks with a ray (median) | predicted peaks (median) | share of peaks with no ray nearby, median (range) | residual to the nearest ray, peaks with a ray (median over voxels of the voxel median, deg) | observed rays attached (median) |
|---|---|---|---|---|---|---|---|---|
| six full-sample voxels | 6 | 0.4136 | 0.0000 | 0.4051 | 58.0 | 0.0% (0.0%-1.7%) | 0.0256 | 217.5 |
| 30 E0 voxels | 30 | 0.0977 | 0.0000 | 0.0904 | 58.0 | 0.0% (0.0%-3.2%) | 0.0141 | 58.0 |
<!-- /table:evalE_e4_truth -->

* On both image sets the median share of predicted peaks with no ray nearby is 0%. So missing rays do not explain the
  high cost at the truth; an earlier draft of this doc said they did, and that was wrong.
* The nearest ray is farther away on the full-sample images (median 0.0256 deg) than on the E0 images (0.0141 deg).
  The full-sample images also attach many more rays (217.5 against 58). Both fit spots of neighbouring grains that
  merge with or shift the voxel's own spots; we did not test that.

Dips away from the truth, at the widths the search uses at level 0. Read: an arm passes when it has no more rays with
a dip than the strict cost (C++ rule), at both depths. The first table is the full sample, the second E0.

<!-- table:evalE_landL0_minima_full -->
| cost | rays with a minimum of depth 0.01, dense, 0 to 5 deg | rays with a minimum of depth 0.05, dense, 0 to 5 deg | rays with a minimum of depth 0.01, 0 to 0.5 deg | rays with a minimum of depth 0.05, 0 to 0.5 deg | rays where the truth is not the lowest point | lowest point when not the truth: median deg (max) | ... and its depth below the truth (median) |
|---|---|---|---|---|---|---|---|
| strict cost, C++ rule | 72 of 120 (51.1-68.3%) | 0 of 120 (0.0-3.1%) | 6 of 120 (2.3-10.5%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| strict cost, lenient rule | 97 of 120 (72.9-86.9%) | 0 of 120 (0.0-3.1%) | 5 of 120 (1.8-9.4%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| E1, level-0 width (16.5 px) | 120 of 120 (96.9-100.0%) | 46 of 120 (30.1-47.3%) | 5 of 120 (1.8-9.4%) | 0 of 120 (0.0-3.1%) | 36 of 120 | 0.0203 (0.0460) | 0.0000 |
| E2, level-0 reach (1.67 frames) | 113 of 120 (88.4-97.1%) | 6 of 120 (2.3-10.5%) | 18 of 120 (9.7-22.5%) | 0 of 120 (0.0-3.1%) | 2 of 120 | 0.2083 (0.2230) | 0.0057 |
| E4, level-0 tau (0.42 deg) | 39 of 120 (24.8-41.3%) | 0 of 120 (0.0-3.1%) | 24 of 120 (13.8-28.0%) | 0 of 120 (0.0-3.1%) | 55 of 120 | 0.0235 (0.2230) | 0.0001 |
| E1 + E2, level-0 widths | 74 of 120 (52.7-69.9%) | 0 of 120 (0.0-3.1%) | 16 of 120 (8.4-20.6%) | 0 of 120 (0.0-3.1%) | 35 of 120 | 0.0960 (0.2950) | 0.0019 |
| E1 + E2 + E3, level-0 widths | 34 of 120 (21.0-37.0%) | 1 of 120 (0.1-4.6%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| E2 + E3, level-0 reach | 115 of 120 (90.6-98.2%) | 10 of 120 (4.6-14.7%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 |  |  |
| E3 + E4, level-0 tau | 46 of 120 (30.1-47.3%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 65 of 120 | 0.0145 (0.0575) | 0.0000 |
| E1 + E2 + E3 + E4, level-0 widths | 16 of 120 (8.4-20.6%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 0 of 120 (0.0-3.1%) | 72 of 120 | 0.0255 (0.0765) | 0.0000 |
<!-- /table:evalE_landL0_minima_full -->

<!-- table:evalE_landL0_minima_e0 -->
| cost | rays with a minimum of depth 0.01, dense, 0 to 0.5 deg | rays with a minimum of depth 0.05, dense, 0 to 0.5 deg | rays where the truth is not the lowest point | lowest point when not the truth: median deg (max) | ... and its depth below the truth (median) |
|---|---|---|---|---|---|
| strict cost, C++ rule | 12 of 600 (1.1-3.5%) | 0 of 600 (0.0-0.6%) | 38 of 600 | 0.0005 (0.0010) | 0.0021 |
| strict cost, lenient rule | 13 of 600 (1.3-3.7%) | 0 of 600 (0.0-0.6%) | 38 of 600 | 0.0005 (0.0010) | 0.0021 |
| E1, level-0 width (16.5 px) | 61 of 600 (8.0-12.8%) | 0 of 600 (0.0-0.6%) | 0 of 600 |  |  |
| E2, level-0 reach (1.67 frames) | 17 of 600 (1.8-4.5%) | 0 of 600 (0.0-0.6%) | 0 of 600 |  |  |
| E4, level-0 tau (0.42 deg) | 70 of 600 (9.3-14.5%) | 0 of 600 (0.0-0.6%) | 292 of 600 | 0.0085 (0.1730) | 0.0000 |
| E1 + E2, level-0 widths | 83 of 600 (11.3-16.8%) | 0 of 600 (0.0-0.6%) | 49 of 600 | 0.1370 (0.4475) | 0.0119 |
| E1 + E2 + E3, level-0 widths | 6 of 600 (0.5-2.2%) | 0 of 600 (0.0-0.6%) | 30 of 600 | 0.0235 (0.0450) | 0.0000 |
| E2 + E3, level-0 reach | 0 of 600 (0.0-0.6%) | 0 of 600 (0.0-0.6%) | 6 of 600 | 0.0008 (0.0030) | 0.0000 |
| E3 + E4, level-0 tau | 6 of 600 (0.5-2.2%) | 0 of 600 (0.0-0.6%) | 124 of 600 | 0.0010 (0.0060) | 0.0000 |
| E1 + E2 + E3 + E4, level-0 widths | 6 of 600 (0.5-2.2%) | 0 of 600 (0.0-0.6%) | 281 of 600 | 0.0080 (0.0330) | 0.0000 |
<!-- /table:evalE_landL0_minima_e0 -->

* E1 alone adds deep dips on the full-sample rays (46 of 120 at depth 0.05; strict cost 0). On the E0 images, which
  hold no other grain, it adds none. That fits the idea that E1 rewards neighbouring grains' spots, but we did not test
  that directly.
* E4 alone adds shallow dips on the E0 rays (70 of 600, strict cost 12).
* E1 + E2 + E3 + E4 has the fewest dips: 16 of 120 full-sample rays at depth 0.01, against 72 for the strict cost.
* E4's lowest point is often a little off the truth (median about 0.01 deg). That is far below what the coarse levels
  resolve, but it is a reason to rank with E4 and never accept on it.

Basin size along the full-sample rays. Read: a larger radius means the cost falls toward the truth from farther away.

<!-- table:evalE_land_basin_full -->
| cost | half-way rise from the truth (median deg) | basin radius, eps 0.01 (median deg) | basin radius, eps 0.01 (lower quartile deg) | rays with no drop over eps out to 5 deg, eps 0.01 | basin radius, eps 0.05 (median deg) | rays with no drop over eps out to 5 deg, eps 0.05 |
|---|---|---|---|---|---|---|
| strict cost, C++ rule | 0.312 | 3.693 | 1.729 | 45 of 120 | 5.000 | 120 of 120 |
| strict cost, lenient rule | 0.333 | 2.783 | 1.971 | 15 of 120 | 5.000 | 120 of 120 |
| S: sampled binary (control) | 0.348 | 1.688 | 1.315 | 0 of 120 | 5.000 | 111 of 120 |
| E1: distance kernel 3 px | 0.667 | 1.697 | 1.287 | 0 of 120 | 5.000 | 95 of 120 |
| E1: distance kernel 12 px | 1.269 | 1.081 | 0.758 | 0 of 120 | 4.701 | 54 of 120 |
| E2: omega blending, 1 frame | 0.482 | 2.111 | 1.626 | 4 of 120 | 5.000 | 110 of 120 |
| E3: soft peak weights | 0.348 | 1.698 | 1.267 | 0 of 120 | 5.000 | 105 of 120 |
| E4: ray residual, tau 0.05 deg | 0.428 | 5.000 | 5.000 | 117 of 120 | 5.000 | 120 of 120 |
| E4: ray residual, tau 0.25 deg | 1.075 | 5.000 | 3.427 | 87 of 120 | 5.000 | 120 of 120 |
| E1 + E2 | 1.064 | 3.750 | 1.901 | 25 of 120 | 5.000 | 111 of 120 |
| E1 + E2 + E3 | 1.069 | 4.238 | 3.420 | 31 of 120 | 5.000 | 109 of 120 |
| E2 + E3 | 0.481 | 2.228 | 1.787 | 0 of 120 | 5.000 | 103 of 120 |
| E3 + E4 (0.25 deg) | 0.549 | 4.983 | 3.837 | 60 of 120 | 5.000 | 120 of 120 |
<!-- /table:evalE_land_basin_full -->

The arms at the level-0 widths, along the same ray:

![The arms at the level-0 widths](../benchmarks/eval_reduction/E/figures/ray_level0_widths.png)

## Measurement 2: ranking the strict search's own candidate lists

Here the search ran with the strict cost as usual. At levels 0, 1 and 2 we took the list after the quick MC and
applied the quarter cut twice: once ordered by the strict cost, once ordered by the arm. The search itself is the
same in both, so this measures the ranking alone.

E0, lenient rule. Read a cell as "kept of present (+fixed / -broken)". The last two columns add up levels 0 to 2.

<!-- table:evalE_rank_e0_old_next -->
| ranking key | level 0: basin kept | level 0: best rank (median) | level 1: basin kept | level 1: best rank (median) | level 2: basin kept | level 2: best rank (median) | levels 0-2: kept (95% Wilson) | levels 0-2: fixed / broken |
|---|---|---|---|---|---|---|---|---|
| strict (baseline) | 186 of 197 | 1 | 166 of 176 | 0 | 160 of 163 | 0 | 512 of 536 (93.4-97.0%) |  |
| S | 188 of 197 (+3 / -1) | 1 | 164 of 176 (+0 / -2) | 0 | 160 of 163 (+0 / -0) | 0 | 512 of 536 (93.4-97.0%) | 3 / 3 |
| E1s3 | 195 of 197 (+9 / -0) | 0 | 175 of 176 (+9 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 532 of 536 (98.1-99.7%) | 20 / 0 |
| E1s8 | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
| E1a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
| E2a | 188 of 197 (+4 / -2) | 0 | 169 of 176 (+3 / -0) | 0 | 160 of 163 (+1 / -1) | 0 | 517 of 536 (94.5-97.7%) | 8 / 3 |
| E3 | 188 of 197 (+4 / -2) | 1 | 163 of 176 (+1 / -4) | 0 | 159 of 163 (+0 / -1) | 0 | 510 of 536 (93.0-96.7%) | 5 / 7 |
| E4t05 | 195 of 197 (+10 / -1) | 0 | 174 of 176 (+8 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 531 of 536 (97.8-99.6%) | 20 / 1 |
| E4t25 | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E4a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E12a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E123a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E23a | 188 of 197 (+4 / -2) | 0 | 170 of 176 (+5 / -1) | 0 | 160 of 163 (+1 / -1) | 0 | 518 of 536 (94.8-97.9%) | 10 / 4 |
| E34a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
| E1234a | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E1s12 | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
| E2w1 | 186 of 197 (+2 / -2) | 1 | 167 of 176 (+3 / -2) | 0 | 160 of 163 (+1 / -1) | 0 | 513 of 536 (93.6-97.1%) | 6 / 5 |
| E12s3 | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 536 of 536 (99.3-100.0%) | 24 / 0 |
| E123s3 | 196 of 197 (+10 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 163 of 163 (+3 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
| E23w1 | 186 of 197 (+3 / -3) | 1 | 169 of 176 (+5 / -2) | 0 | 160 of 163 (+1 / -1) | 0 | 515 of 536 (94.1-97.4%) | 9 / 6 |
| E34t25 | 197 of 197 (+11 / -0) | 0 | 176 of 176 (+10 / -0) | 0 | 162 of 163 (+2 / -0) | 0 | 535 of 536 (99.0-100.0%) | 23 / 0 |
<!-- /table:evalE_rank_e0_old_next -->

* The strict ranking keeps the basin in 512 of 536 level cases. The arms built on E1 or E4 keep it in 535 or 536, and
  break none.
* E2 alone, E3 alone and the control S barely change anything. The gain comes from credit for near misses and from the
  observed rays.
* The basin candidates the strict cost cuts are mostly 1 to 3 deg from the truth, where their spots no longer overlap
  the truth's. That is exactly where a smooth cost helps.
* The C++ rule gives the same picture: strict 474 of 508, E1 + E2 + E3 + E4 all 508 (table `evalE_rank_e0_cpp_next`
  in `benchmarks/eval_reduction/E/rank/tables.md`).
* On the six full-sample voxels the strict ranking already keeps the basin in all 18 level cases, and so does every
  arm.

## The choice of arms for the end-to-end test

The rule was written down before the ranking data was in. An arm qualifies when it keeps the basin more often than it
loses it, and adds no dips beyond the strict cost's, on both image sets. At most three arms are carried.

* Two arms qualified: **E3 + E4** (tau 0.25 deg) and **E1 + E2 + E3 + E4**.
* The closest miss was E1 + E2 + E3: one full-sample ray out of 120 had a dip of depth 0.05.
* E1 alone, E4 alone, E2 alone, E1 + E2 and E2 + E3 add dips and were not carried. The third slot stayed empty.

## Measurement 3: the whole search with the smooth ranking on

`SmoothCostRank` on, levels 0-2 ranked by the arm, everything else strict, compared voxel by voxel with the strict
search. Every broken voxel would be traced to its cause (owner rule D-20261014-N13). None broke, so none needed
tracing.

Read: "wrong" counts carry a 95% interval. "Fixed: traced to arm" means the arm's ranking kept a basin candidate that
the strict cost would have cut on the same list. "Traced to draw" means the fix came from a later knock-on difference,
such as other random draws.

<!-- table:evalE_e2e -->
| voxel set | rule | arm (SmoothCostRank) | n | wrong, baseline | wrong, arm | fixed | broken | both wrong | broken: traced to arm / draw / untraced | fixed: traced to arm / draw | median error of voxels right in both (deg), baseline / arm | strict evaluations per voxel, baseline | strict evaluations per voxel, arm | smooth evaluations per voxel |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E0 200 | lenient | E3+E4,tau=0.25 | 200 | 41 (15.5-26.6%) | 7 (1.7-7.0%) | 34 | 0 | 7 | 0 / 0 / 0 | 17 / 17 | 0.0130 / 0.0054 | 48700 | 49215 | 252 |
| E0 200 | C++ | E3+E4,tau=0.25 | 200 | 56 (22.2-34.6%) | 9 (2.4-8.3%) | 47 | 0 | 9 | 0 / 0 / 0 | 21 / 26 | 0.0223 / 0.0082 | 47230 | 47551 | 172 |
| six full-sample | lenient | E3+E4,tau=0.25 | 6 | 0 (0.0-39.0%) | 0 (0.0-39.0%) | 0 | 0 | 0 | 0 / 0 / 0 | 0 / 0 | 0.0181 / 0.0163 | 332046 | 334486 | 15030 |
| six full-sample | C++ | E3+E4,tau=0.25 | 6 | 0 (0.0-39.0%) | 0 (0.0-39.0%) | 0 | 0 | 0 | 0 / 0 / 0 | 0 / 0 | 0.0286 / 0.0291 | 228146 | 235524 | 11084 |
| E0 200 | lenient | E1+E2+E3+E4 | 200 | 41 (15.5-26.6%) | 5 (1.1-5.7%) | 36 | 0 | 5 | 0 / 0 / 0 | 15 / 21 | 0.0130 / 0.0050 | 48700 | 49653 | 263 |
| E0 200 | C++ | E1+E2+E3+E4 | 200 | 56 (22.2-34.6%) | 9 (2.4-8.3%) | 47 | 0 | 9 | 0 / 0 / 0 | 21 / 26 | 0.0223 / 0.0112 | 47230 | 47767 | 177 |
| six full-sample | lenient | E1+E2+E3+E4 | 6 | 0 (0.0-39.0%) | 0 (0.0-39.0%) | 0 | 0 | 0 | 0 / 0 / 0 | 0 / 0 | 0.0181 / 0.0284 | 332046 | 348018 | 15755 |
<!-- /table:evalE_e2e -->

* **E0, 200 voxels.** Lenient rule: 41 wrong become 7 (E3 + E4) and 5 (E1 + E2 + E3 + E4). C++ rule: 56 become 9 with
  both arms. About half of the fixes are traced to the arm itself.
* **No breaks, and the pairing is sound.** The code review checked it. Each arm run and its baseline cover the same
  200 voxels, with the same seeds, the same images and the same rule. The baseline was run on this branch with the key
  off, and it equals Task C's baseline on all 200 voxels.
* **Why zero breaks is plausible.** With the arms, only 7 or 5 of 200 voxels stay wrong (lenient rule), and 9 (C++
  rule). Every one of them ends between 59.9 and 60.0 deg from the truth, a near-degenerate solution, and every one
  was already wrong in the baseline. In the baseline, 30 of the 41 wrong voxels (lenient) end within 1 deg of 60 deg;
  the arms fix 10 of the other 11. So the arms seem to remove the losses at the quarter cut and leave the
  near-degenerate solutions. That reading fits the counts; it was not tested.
* **The "14 to 23 flips" calibration does not apply here.** Changing only the random seed of the strict search flips
  14 to 23 E0 voxels each way (`s3_ab_paired` in `benchmarks/s3_steps1_5/tables.md`). Those pairs compare two runs
  with about the same number of wrong voxels (35 to 39 of 200 each). An arm that leaves only 5 to 9 wrong voxels has
  far fewer voxels that could break, so that figure does not set the number of breaks to expect.
* **Voxels right in both runs end closer to the truth on E0.** E0, lenient rule, the 159 voxels right in both runs:
  median error 0.0130 deg became 0.0054 deg with E3 + E4 and 0.0050 deg with E1 + E2 + E3 + E4. On the six
  full-sample voxels the medians moved both ways: E3 + E4, lenient rule, 0.0181 to 0.0163 deg (2 voxels farther, 4
  closer); E3 + E4, C++ rule, 0.0286 to 0.0291 deg (3 farther, 3 closer); E1 + E2 + E3 + E4, lenient rule, 0.0181 to
  0.0284 deg (3 farther, 3 closer). With six voxels none of these is a measured difference.
* **Six full-sample voxels:** 0 of 6 wrong before and after. E1 + E2 + E3 + E4 was run on them with the lenient rule
  only; the C++-rule run of that arm on the six voxels was not done.
* **Evaluations:** the strict cost's own evaluations rose by 0.7% to 4.8%. The smooth cost is called separately. With
  the lenient rule that is 252 to 263 times per E0 voxel and 15,030 to 15,755 times per full-sample voxel. With the
  C++ rule it is 172 to 177 times per E0 voxel and 11,084 per full-sample voxel (E3 + E4 only).

## Measurement 4: time

Two kinds of time were measured. The time per evaluation was taken in one process on a quiet machine (the timing
preflight passed), at the level-0 widths, which are the slowest. The time of whole voxels was also taken in one
process. For the four E0 voxels the preflight passed (single-worker). For the full-sample voxel it failed (a busy
machine), so that time is labelled contended, and its repeat on a quiet machine was not done (owner decision
D-20261018-N16).

Time per evaluation. Read: "x strict" is how many times slower than the strict cost.

<!-- table:evalE_timing_ms -->
| cost | full-sample: ms per evaluation, batch of 256 | full-sample: ms per evaluation, one per call | full-sample: x strict | E0: ms per evaluation, batch of 256 | E0: ms per evaluation, one per call | E0: x strict |
|---|---|---|---|---|---|---|
| strict cost (one call per orientation) | 0.381 | 0.381 | 1.000 | 0.392 | 0.392 | 1.000 |
| S | 0.470 | 0.809 | 1.233 | 0.228 | 0.601 | 0.581 |
| E1 | 1.332 | 1.649 | 3.492 | 0.240 | 0.634 | 0.613 |
| E2 | 2.080 | 2.675 | 5.453 | 0.886 | 1.590 | 2.261 |
| E3 | 0.473 | 0.844 | 1.240 | 0.227 | 0.640 | 0.580 |
| E4 | 1.186 | 2.153 | 3.109 | 0.313 | 1.248 | 0.798 |
| E1+E2 | 6.663 | 7.268 | 17.466 | 0.950 | 1.710 | 2.424 |
| E1+E2+E3 | 6.728 | 7.485 | 17.637 | 0.952 | 1.747 | 2.430 |
| E2+E3 | 2.136 | 2.792 | 5.599 | 0.888 | 1.627 | 2.265 |
| E3+E4 | 1.639 | 2.811 | 4.297 | 0.495 | 1.632 | 1.262 |
| E1+E2+E3+E4 | 7.903 | 9.499 | 20.718 | 1.223 | 2.745 | 3.120 |
<!-- /table:evalE_timing_ms -->

Whole voxels, key off against key on (`scripts/eval_reduction/smooth_voxel_timing.py`, lenient rule, the seeds of
Measurement 3). Read: "extra time" is the key-on time over the key-off time, minus one. Each run ended with the same
error as in Measurement 3, so these are the same searches.

<!-- table:evalE_timing_voxel_wall -->
| voxel set | voxel | timing | key off (s) | E3+E4,tau=0.25 (s) | E3+E4,tau=0.25: extra time | E1+E2+E3+E4 (s) | E1+E2+E3+E4: extra time | same error as the end-to-end run |
|---|---|---|---|---|---|---|---|---|
| E0 | 7336 | single-worker, median of 2 | 2.1 | 6.0 | +179% | 6.4 | +201% | yes |
| E0 | 3108 | single-worker, median of 2 | 2.4 | 6.5 | +173% | 7.1 | +195% | yes |
| E0 | 2910 | single-worker, median of 2 | 3.0 | 7.0 | +130% | 7.7 | +153% | yes |
| E0 | 16182 | single-worker, median of 2 | 2.5 | 6.5 | +164% | 6.9 | +179% | yes |
| full-sample | 15901 | contended, median of 1 | 117.4 | 159.6 | +36% | 263.9 | +125% | yes |
<!-- /table:evalE_timing_voxel_wall -->

* **On E0 the smooth ranking more than doubles the time of a voxel** (single-worker). The extra time per E0 voxel
  is close to the 3.6 s it takes to build a `SmoothCost` for one E0 image set. Every E0 voxel has its own images, so
  the build is paid on every voxel, and a cache per image set would not help there.
* **On the full sample** (voxel 15901, contended) E3 + E4 adds 36% and E1 + E2 + E3 + E4 adds 125%.

The full-sample voxel as a model. Per voxel: the strict cost's evaluations times its time per evaluation, plus the
smooth evaluations times theirs, plus the build. Today the library builds a `SmoothCost` in every `reconstruct_voxel`
call (16.6 s on the full sample). That is the "as shipped" column. A cache per image set is not built yet; Task G
builds it if it uses the smooth key (D-20261018-N17). So the last column is a projection, not a measurement.

<!-- table:evalE_timing_per_voxel -->
| voxel set | arm | modelled strict cost time per voxel, baseline (s) | modelled strict cost time per voxel, arm (s) | smooth evaluations per voxel | modelled smooth cost time per voxel (s) | SmoothCost build (s) | extra time as shipped (a build per voxel) | projection: one build per image set (cache not built) |
|---|---|---|---|---|---|---|---|---|
| six full-sample | E3+E4,tau=0.25 | 126.7 | 127.6 | 15030 | 24.6 | 16.6 | +33.3% | +20.2% |
| six full-sample | E1+E2+E3+E4 | 126.7 | 132.8 | 15755 | 124.5 | 16.6 | +116.2% | +103.1% |
<!-- /table:evalE_timing_per_voxel -->

* As shipped, the model gives +33.3% (E3 + E4) and +116.2% (E1 + E2 + E3 + E4). The contended wall time of voxel 15901
  (+36% and +125%) is close to that.
* With one build per image set, the projection is +20.2% and +103.1%.
* The model uses the time per evaluation at a batch of 256. Shorter lists cost more per evaluation: one orientation
  per call takes 2.811 ms against 1.639 ms at a batch of 256 (E3 + E4, full sample). So the model is not an upper
  estimate.
* The model is not used for E0. There most of the counted evaluations are screen and re-cost evaluations, which cost
  far less than one full evaluation, so the model overstated the E0 voxel time several times over (Task E review).

E1 is slow on the full sample. That is consistent with its nearest-lit-pixel search, which covers millions of lit
pixels there; we did not profile it. On the sparse E0 images it costs about the same as the strict cost.

## What smoothing removed and what it cost

**Removed:**
* the pixel staircase (E1);
* the frame edge (E2);
* the jumps when a peak enters or leaves the list, and the fixed half credit at the truth (E3);
* on E0, most of the truth-basin losses at the quarter cut, and with them most wrong voxels.

**Cost:**
* time: as shipped, about +33% (E3 + E4) or +116% (E1 + E2 + E3 + E4) per full-sample voxel (modelled; a contended
  check on one voxel gave +36% and +125%). With a build cache per image set, not built yet, the projection is about
  +20% and +103%. On E0, more than double (measured);
* false dips when an arm is used alone, consistent with credit given to other grains' spots (not tested);
* weaker separation far from the truth. That is why the smooth costs only rank, after the quick MC, and the strict
  cost still accepts.

**Not shown:** any gain on full-sample voxels that are wrong today.

## What this means

* The smooth ranking is a strong candidate, but only E0 evidence supports it. Before any default changes, it needs a
  full-sample set that has wrong voxels: the other 24 voxels of the 30-voxel set, which Task G runs before the smooth
  key is chosen (owner decision D-20261018-N15).
* E3 + E4 (tau 0.25 deg) is the cheaper candidate.
* For a full-sample run, the `SmoothCost` build should be cached per image set. It is not built yet; Task G builds it
  if it uses the smooth key (D-20261018-N17).
* E4's gradient is ready for Task F's outward search.
* **For Track S (S4).** All four arms could in principle be computed from a per-sample level-0 table built by
  translating spots. E1, E2 and E3 use the translated spot points or position-free quantities. E4 needs only each
  voxel's observed rays, which come from the images once.

## Limits

* E0 is the only set with wrong voxels, and its images hold one voxel each, with no neighbouring grains.
* The full-sample test was cut from 30 voxels to the six interior ones, for time. That was decided before any
  full-sample ranking result was read. The other 24 voxels run in Task G (D-20261018-N15).
* The whole-voxel time of the full-sample voxel in Measurement 4 was taken on a busy machine (contended). Its repeat
  on a quiet machine was not done (D-20261018-N16).
* The arm widths were set once, by rule of thumb.
* E2 is a tolerance across frame edges, not a model of how a real spot spreads over frames.
* E3's fade-out ramps are fixed by hand.
* E4 inherits Task D's assumptions (two detectors, binary centroids, pairing settings chosen on this sample).
* The arms rank only. Using a smooth cost to move the quick MC was not tested.
* Clean images only, one image family, no extra seeds.
