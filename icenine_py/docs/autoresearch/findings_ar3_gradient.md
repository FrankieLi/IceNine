# AR3 step 3: paint as a start, then refine, on grains with gradients

Step 3 of "AR3" in `docs/autoresearch/PLAN.md`, task `ar3c` (D-20261110-N4). New:
`icenine/gradient_sample.py` (no config key, no reconstructor change). Scripts
`scripts/autoresearch/ar3/`; results `benchmarks/autoresearch/ar3/gradient/`. This page was placed
by the main session from the text the implementer returned (the tool refused the implementer's own
write of the file); the tables are filled from the committed `tables.md`. The review fixes of
2026-11-10 are applied from the stored records; no run was added.

Terms: a **grain entry** is one accepted seed; to **paint** is to give each voxel the passing grain
of lowest strict cost; the **start** is the painted orientation; the **refinement** is one BFS
local fit from the start; a grain's **level** is its largest misorientation from its centre
orientation (0, 0.1, 0.2, 0.5, 1.0 deg); **error** is the symmetry-reduced misorientation to the
voxel's own truth, "wrong" above 1 deg; a voxel is **hurt** by a shortcut when it ends more than
0.05 deg farther from the truth than the full fit. Counts by voxel, Wilson 95% intervals, no
p-values; a difference is called a gain or loss only when intervals do not overlap.

The sample: the layer's truth (24,570 voxels, 497 grains); each grain keeps its orientation at its
centre and turns linearly with position to its level (98/81/101/101/116 grains at
0/0.1/0.2/0.5/1.0). Measured spread equals the level to 5.1e-5 deg. Clean images rendered with
`render_full.py --mic` (tag `grad`): 9,444,080 lit pixels and 104,912 spots, against 7,414,841 and
57,389 for the stored images. The packed image stack has 9,444,080 non-zero pixels, the grad
render's count; how it was built, and the sha256 of the truth map, are in `grad_stack.json`
(`scripts/autoresearch/ar3/grad_stack_note.py`). The gradient is a modelling choice, not
measured on a real sample.

Runs (Q-max 8, FusedCost, one worker, voxel order of default_rng(0)): P paint through the key; R
one fit from the start on each of the 19,013 painted voxels; R2 a second fit on every fourth
(4,768); D the BFS with indexed seeds. Trace check: first 200 fits identical with trace on and off
(200 of 200).

## The result in one picture

Reading the picture. Top left: error by level (median solid, 95th dashed) for start, after R, and
D. Top middle: share of voxels where the centre orientation hits under 0.9 times the true one
(Reading 0). Top right: evaluations skipped against share hurt for the skip rules (dotted line:
1%). Bottom left: share of a fit's evaluations after it is within 0.02 deg of its end. Bottom
middle: stop inside the fit. Bottom right: wrong voxels by level.

![Error by level, Reading 0, skip and stop rules, wrong voxels](../../benchmarks/autoresearch/ar3/gradient/ar3_gradient.png)

- The cost can tell the truth from the painted start, but the fit does not get there. On 16,591
  of 19,013 painted voxels (86.8 to 87.7%) the strict cost of the true orientation is below the
  cost of the start, and on all 2,240 voxels whose start is more than 0.1 deg off. The fit lowered
  the cost on only 36 of the 19,013.
- The fit's trial steps are many times the start's error. On the painted voxels the start is a
  median 0.0406 deg from the truth. In the fit's first block the trial steps have a median of
  0.70 deg (5th percentile 0.32 deg). After a failed block the search restarts a median 4.9 to
  16.8 deg from the start, with larger steps. 18,500 of the 19,013 fits are four failed blocks.
  This is consistent with a step too coarse for this start; a run with a finer box would test it
  (planned as step 4).
- The refinement did nothing detectable: 9 voxels came closer by more than 0.05 deg and 2 were
  hurt; the fit ended within 0.02 deg of its start on 18,977 (99.7 to 99.9%). Skipping or stopping
  the fit saves almost everything at no measured cost here. That is the fit not moving, not a
  finding about when a fit can stop.
- The start: 5,369 entries for 497 grains (21.9% of voxels are seeds); median error 0.0342 deg
  over all voxels, 170 wrong of 24,570 (0.60 to 0.80%); 78,142,641 evaluations, 74.8% of them in
  172 grid-fallback entries. A start this close is consistent with the many seeds; untested,
  since no arm with fewer seeds was run.
- Under the plan's rule P beats D on the share of voxels above 0.1 deg: 2,531 (9.9 to 10.7%) with
  R against 3,545 (14.0 to 14.9%), and at levels 0.1, 0.2 and 0.5. It also beats D on wrong voxels
  at level 0 (3 against 17), not overall (170 against 203, intervals overlap). The gain is P's (R
  adds nothing), at 90,985,770 evaluations against 118,350,225. One run of each arm.

## Reading 0: what an offset costs (truth used; analysis only)

Reading the table. One strict score of each voxel's true orientation and of its grain's centre
orientation. `ht`, `hc`: median hit ratios; `rm`: median of hc/ht; `r5` its 5th percentile;
`below`: voxels under the 0.9 test; `ct`, `cc`: median costs. The images took the gradient: level 0
is unchanged, higher levels lose more hit ratio.

<!-- table:ar3c_reading0 -->
| level | n | ht | hc | rm | r5 | below | ct | cc |
|---|---|---|---|---|---|---|---|---|
| level 0 | 4,877 | 1.000 | 1.000 | 1.000 | 1.000 | 0 of 4877 (0.00 to 0.08%) | 0.073 | 0.073 |
| level 0.1 | 4,488 | 1.000 | 0.951 | 0.957 | 0.815 | 979 of 4488 (20.6 to 23.0%) | 0.074 | 0.118 |
| level 0.2 | 5,041 | 1.000 | 0.897 | 0.901 | 0.667 | 2498 of 5041 (48.2 to 50.9%) | 0.063 | 0.159 |
| level 0.5 | 4,625 | 1.000 | 0.763 | 0.767 | 0.451 | 3618 of 4625 (77.0 to 79.4%) | 0.068 | 0.285 |
| level 1 | 5,539 | 1.000 | 0.594 | 0.596 | 0.281 | 5168 of 5539 (92.6 to 93.9%) | 0.068 | 0.449 |
<!-- /table:ar3c_reading0 -->

A flood from a centre orientation stops early at high levels: about half the voxels at level 0.2,
about 93% at level 1.0, would fail the 0.9 test against their own true orientation.

## 1. The start

Reading the table. Error of the painted orientation of every voxel by level.

<!-- table:ar3c_start -->
| group | n | median | p95 | gt01 | wrong |
|---|---|---|---|---|---|
| all voxels | 24,570 | 0.0342 | 0.1304 | 2541 of 24570 (10.0 to 10.7%) | 170 of 24570 (0.60 to 0.80%) |
| level 0 | 4,877 | 0.0161 | 0.0518 | 3 of 4877 (0.02 to 0.18%) | 3 of 4877 (0.02 to 0.18%) |
| level 0.1 | 4,488 | 0.0313 | 0.0699 | 41 of 4488 (0.7 to 1.2%) | 10 of 4488 (0.12 to 0.41%) |
| level 0.2 | 5,041 | 0.0401 | 0.0896 | 151 of 5041 (2.6 to 3.5%) | 23 of 5041 (0.30 to 0.68%) |
| level 0.5 | 4,625 | 0.0571 | 0.1375 | 773 of 4625 (15.7 to 17.8%) | 39 of 4625 (0.6 to 1.2%) |
| level 1 | 5,539 | 0.0551 | 0.1875 | 1573 of 5539 (27.2 to 29.6%) | 95 of 5539 (1.4 to 2.1%) |
<!-- /table:ar3c_start -->

Reading the table. Where voxels came from by level (`grain` painted, `fit` settled by a fit before
seeding, `seed` an entry, `failed`).

<!-- table:ar3c_sources -->
| group | n | grain | fit | seed | failed | none |
|---|---|---|---|---|---|---|
| all voxels | 24,570 | 19013 | 184 | 5369 | 4 | 0 |
| level 0 | 4,877 | 4774 | 0 | 103 | 0 | 0 |
| level 0.1 | 4,488 | 4152 | 2 | 334 | 0 | 0 |
| level 0.2 | 5,041 | 4332 | 1 | 707 | 1 | 0 |
| level 0.5 | 4,625 | 3196 | 48 | 1381 | 0 | 0 |
| level 1 | 5,539 | 2559 | 133 | 2844 | 3 | 0 |
<!-- /table:ar3c_sources -->

Reading the table. Entries per true grain (`multi`: grains with more than one entry; `none`: grains
with no entry; `grid`: grid-fallback entries).

<!-- table:ar3c_entries -->
| level | grains | entries | mean | mx | multi | none | grid |
|---|---|---|---|---|---|---|---|
| level 0 | 98 | 103 | 1.05 | 2 | 5 of 98 (2.2 to 11.4%) | 0 | 1 |
| level 0.1 | 81 | 334 | 4.12 | 7 | 81 of 81 (95.5 to 100.0%) | 0 | 1 |
| level 0.2 | 101 | 707 | 7.00 | 13 | 101 of 101 (96.3 to 100.0%) | 0 | 2 |
| level 0.5 | 101 | 1381 | 13.67 | 29 | 101 of 101 (96.3 to 100.0%) | 0 | 31 |
| level 1 | 116 | 2844 | 24.52 | 48 | 116 of 116 (96.8 to 100.0%) | 0 | 137 |
<!-- /table:ar3c_entries -->

Of 5,369 entries, 4,844 are within 1 deg of an earlier entry; 4,748 of them lie in the same true
grain as that entry and 96 in another. The median entry floods to 3 voxels (its `flood_passes`).
The fit-settled branch was tried 5,263 times and accepted 184 times (all at levels 0.1 and above).

Reading the tables. Evaluations by stage, then per entry. The 172 grid-fallback entries take
74.8% of P's evaluations (76.6% with the 4 rejected grid searches of the failed voxels). The
fallback is part of the key's design and ran in D too: 142 times, 48,273,080 evaluations.

<!-- table:ar3c_evals_stage -->
| stage | evals |
|---|---|
| seed_index | 14,106,247 |
| score | 291,783 |
| fit_rejected | 3,459,521 |
| fit_accepted | 124,976 |
| index_rejected | 336,721 |
| seed_grid | 58,461,276 |
| grid_rejected | 1,362,117 |
| total | 78,142,641 |
<!-- /table:ar3c_evals_stage -->

<!-- table:ar3c_seed_cost -->
| what | n | per |
|---|---|---|
| entries from the seed path | 5,197 | 2,714 |
| entries from the grid search (fallback) | 172 | 339,891 |
<!-- /table:ar3c_seed_cost -->

The whole-layer paint on step 2's sample (one orientation per grain) cost 1,316,345 evaluations
with 528 entries; here 78,142,641 with 5,369. That is 59 times the evaluations with the grid
fallback and 13.9 times without it (18,319,248), for 10.2 times the entries. The key-on BFS of
that sample (33,040,126) is not a like-for-like comparison.

## 2. The refinement

Reading the table. The 19,013 painted voxels only. Error at the start, after R, and in D.

<!-- table:ar3c_refine_all -->
| group | n | median | p95 | gt01 | wrong |
|---|---|---|---|---|---|
| start, painted voxels | 19,013 | 0.0406 | 0.1329 | 2240 of 19013 (11.3 to 12.2%) | 87 of 19013 (0.37 to 0.56%) |
| after R | 19,013 | 0.0406 | 0.1325 | 2230 of 19013 (11.3 to 12.2%) | 87 of 19013 (0.37 to 0.56%) |
| D (BFS) | 19,013 | 0.0440 | 0.1397 | 2601 of 19013 (13.2 to 14.2%) | 105 of 19013 (0.46 to 0.67%) |
<!-- /table:ar3c_refine_all -->

Reading the table. By level (`closer`, `hurt`: moved more than 0.05 deg toward or away from the
truth; `med`, `p95`, `w`: median, 95th percentile, wrong count of start (s), R (r), D (d)).

<!-- table:ar3c_refine_level -->
| level | n | closer | hurt | med_s | med_r | med_d | p95_s | p95_r | p95_d | w_s | w_r | w_d |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| level 0 | 4,774 | 0 of 4774 (0.00 to 0.08%) | 0 of 4774 (0.00 to 0.08%) | 0.0161 | 0.0161 | 0.0170 | 0.0518 | 0.0518 | 0.0554 | 3 | 3 | 16 |
| level 0.1 | 4,152 | 0 of 4152 (0.00 to 0.09%) | 0 of 4152 (0.00 to 0.09%) | 0.0324 | 0.0324 | 0.0448 | 0.0704 | 0.0704 | 0.0968 | 9 | 9 | 20 |
| level 0.2 | 4,332 | 0 of 4332 (0.00 to 0.09%) | 1 of 4332 (0.00 to 0.13%) | 0.0433 | 0.0433 | 0.0559 | 0.0917 | 0.0917 | 0.1185 | 20 | 20 | 27 |
| level 0.5 | 3,196 | 2 of 3196 (0.02 to 0.23%) | 0 of 3196 (0.00 to 0.12%) | 0.0719 | 0.0718 | 0.0765 | 0.1440 | 0.1432 | 0.1647 | 25 | 25 | 26 |
| level 1 | 2,559 | 7 of 2559 (0.13 to 0.56%) | 1 of 2559 (0.01 to 0.22%) | 0.1022 | 0.1018 | 0.0723 | 0.1989 | 0.1986 | 0.1914 | 30 | 30 | 16 |
| all | 19,013 | 9 of 19013 (0.02 to 0.09%) | 2 of 19013 (0.00 to 0.04%) | 0.0406 | 0.0406 | 0.0440 | 0.1329 | 0.1325 | 0.1397 | 87 | 87 | 105 |
<!-- /table:ar3c_refine_level -->

Reading the table. By the size of the start's error (`moved`: more than 0.02 deg).

<!-- table:ar3c_refine_bins -->
| err | n | moved | closer | hurt |
|---|---|---|---|---|
| 0 to 0.05 | 11,423 | 3 of 11423 (0.01 to 0.08%) | 0 of 11423 (0.00 to 0.03%) | 0 of 11423 (0.00 to 0.03%) |
| 0.05 to 0.1 | 5,350 | 11 of 5350 (0.11 to 0.37%) | 1 of 5350 (0.00 to 0.11%) | 1 of 5350 (0.00 to 0.11%) |
| 0.1 to 0.2 | 2,042 | 20 of 2042 (0.6 to 1.5%) | 6 of 2042 (0.13 to 0.64%) | 1 of 2042 (0.01 to 0.28%) |
| 0.2 to 0.5 | 111 | 2 of 111 (0.5 to 6.3%) | 2 of 111 (0.5 to 6.3%) | 0 of 111 (0.0 to 3.3%) |
| 0.5 and above | 87 | 0 of 87 (0.0 to 4.2%) | 0 of 87 (0.0 to 4.2%) | 0 of 87 (0.0 to 4.2%) |
<!-- /table:ar3c_refine_bins -->

Reading the table. Evaluations of each run.

<!-- table:ar3c_evals -->
| arm | evals |
|---|---|
| P (paint) | 78,142,641 |
| R (one fit on each painted voxel) | 12,843,129 |
| P plus R | 90,985,770 |
| R2 (second fit, every fourth voxel) | 3,217,669 |
| D (BFS, indexed seeds) | 118,350,225 |
<!-- /table:ar3c_evals -->

### Why the fit moves almost nothing

Two questions, answered from the stored records (the truth is used here for analysis only). Can
the strict cost tell the start from the truth? And can the fit's trial steps reach the truth?

Reading the table. Each painted voxel's strict cost at its true orientation (Reading 0) against
the cost at its painted start (run R's record), by level and by the start's error in degrees.
`truth_lower`: voxels where the truth's cost is strictly below the start's; `equal`,
`truth_higher`: the other two cases; `fit_lowered`: voxels where R's fit ended below the start's
cost; `med_gap`: median of the start's cost minus the truth's cost.

<!-- table:ar3c_truth_vs_start -->
| group | n | truth_lower | equal | truth_higher | fit_lowered | med_gap |
|---|---|---|---|---|---|---|
| all painted voxels | 19,013 | 16591 of 19013 (86.8 to 87.7%) | 2219 | 203 | 36 of 19013 (0.14 to 0.26%) | 0.0207 |
| level 0 | 4,774 | 2938 of 4774 (60.2 to 62.9%) | 1713 | 123 | 0 of 4774 (0.00 to 0.08%) | 0.0029 |
| level 0.1 | 4,152 | 3760 of 4152 (89.6 to 91.4%) | 358 | 34 | 0 of 4152 (0.00 to 0.09%) | 0.0161 |
| level 0.2 | 4,332 | 4178 of 4332 (95.9 to 97.0%) | 131 | 23 | 7 of 4332 (0.08 to 0.33%) | 0.0231 |
| level 0.5 | 3,196 | 3167 of 3196 (98.7 to 99.4%) | 16 | 13 | 10 of 3196 (0.17 to 0.58%) | 0.0346 |
| level 1 | 2,559 | 2548 of 2559 (99.2 to 99.8%) | 1 | 10 | 19 of 2559 (0.5 to 1.2%) | 0.0506 |
| start error 0 to 0.05 | 11,423 | 9058 of 11423 (78.5 to 80.0%) | 2180 | 185 | 3 of 11423 (0.01 to 0.08%) | 0.0100 |
| start error 0.05 to 0.1 | 5,350 | 5293 of 5350 (98.6 to 99.2%) | 39 | 18 | 11 of 5350 (0.11 to 0.37%) | 0.0345 |
| start error 0.1 to 0.2 | 2,042 | 2042 of 2042 (99.8 to 100.0%) | 0 | 0 | 20 of 2042 (0.6 to 1.5%) | 0.0588 |
| start error 0.2 to 0.5 | 111 | 111 of 111 (96.7 to 100.0%) | 0 | 0 | 2 of 111 (0.5 to 6.3%) | 0.0815 |
| start error 0.5 and above | 87 | 87 of 87 (95.8 to 100.0%) | 0 | 0 | 0 of 87 (0.0 to 4.2%) | 0.0738 |
<!-- /table:ar3c_truth_vs_start -->

So the cost can tell them apart for most voxels: the truth's cost is lower on 16,591 of 19,013,
and on all 2,240 voxels whose start is more than 0.1 deg off. The truth's hit ratio is also above
the start's on 16,953 voxels. A flat or quantised cost at this scale does not explain the
inertness. Only at level 0, where the start is closest, are ties common.

Reading the table. The trial steps of the fit at the run's settings (LocalOrientationGridRadius
5, MinLocalResolution 0, so a box of 5 deg; MaxMCSteps 200), as in a fit whose blocks all fail.
`radius_deg`: the block's subregion radius; `steps`: its trial steps; `evals_at_end`: evaluations
when it ends; `step_med`, `step_p05`: median and 5th percentile of a trial step's angle;
`restart_med`: median angle from the start of the point the block starts from (block 1 starts at
the start). The angles are measured with the library's `QuaternionGrid.get_near_identity_point`
on 10,000 draws per row, in `gradient_summary.py`.

<!-- table:ar3c_steps -->
| block | radius_deg | steps | evals_at_end | step_med | step_p05 | restart_med |
|---|---|---|---|---|---|---|
| 1 | 0.724 | 10 | 12 | 0.70 | 0.32 | none |
| 2 | 1.447 | 19 | 32 | 1.41 | 0.64 | 4.9 |
| 3 | 2.894 | 126 | 159 | 2.82 | 1.31 | 9.7 |
| 4 | 5.000 | 501 | 661 | 4.86 | 2.29 | 16.8 |
<!-- /table:ar3c_steps -->

With the final evaluation, these four blocks make 662 evaluations, the median at every level.
18,500 of the 19,013 R fits are exactly this pattern of four failed blocks. Only block 1 searches
near the start, with 10 trials whose steps are many times the start's error (median 0.0406 deg).
A trial counts only if its cost is strictly lower, so a trial almost never lands in the small
region near the truth where the cost is lower. After block 1 fails, the search moves degrees away.
This is consistent with a step too coarse for this start. It is untested: a run of R with a finer
box, paired on voxels against R, would test it (planned as step 4).

The C++ source has the same logic: `ContinuousSearch.h` line 70 sets the box width to
LocalGridRadius / 2^nLocalResolution, `SearchDetails.cpp` line 50 sets nLocalResolution to
MinLocalResolution, and `AdaptiveSamplingZeroTemp` (`OrientationSearch.cpp` lines 206 to 287)
has the same radius, step count, doubling and restart. The C++ program itself was not run here, so
the same inertness in C++ at these settings is expected, not measured.

## 3. Stop before the fit

Skip when the start's hit ratio is at least t times the entry's best (t 0.90 to 1.00) or when its
strict cost is at most c (deciles of A's start costs). Picked on pick set A (even grain ids): the
rule skipping the most evaluations while hurting at most 1% of A's painted voxels; reported on test
set B (odd ids), threshold fixed. Both picks sit at the grid's edge: t = 0.90, the lowest t, and
c = 0.3054, the top of the grid (the 100th percentile, that is A's largest start cost, so every A
voxel is skipped). Every threshold was within the limit because the fit moved nothing worth the
name. The rules are uninformative here.

Reading the tables. The picked rules on B, then by level. `A_share`: share of A's fit evaluations
skipped; `B_skipped`: B's painted voxels whose fit is skipped; `B_evals`: B's fit evaluations
skipped; `B_hurt`: B's painted voxels hurt by the skip; `B_med`, `B_p95`, `B_gt01`, `B_wrong`:
median error, 95th percentile, share above 0.1 deg and share wrong of the resulting map over all
of B's voxels (12,184).

<!-- table:ar3c_skip -->
| rule | threshold | A_share | B_skipped | B_evals | B_hurt | B_med | B_p95 | B_gt01 | B_wrong |
|---|---|---|---|---|---|---|---|---|---|
| start hit ratio >= t x entry best | 0.9 | 0.998 | 9299 of 9305 | 6,307,982 of 6,311,954 | 3 of 9305 (0.01 to 0.09%) | 0.0365 | 0.1378 | 1467 of 12184 (11.5 to 12.6%) | 86 of 12184 (0.57 to 0.87%) |
| start cost <= c | 0.3054 | 1.000 | 9305 of 9305 | 6,311,954 of 6,311,954 | 3 of 9305 (0.01 to 0.09%) | 0.0365 | 0.1378 | 1467 of 12184 (11.5 to 12.6%) | 86 of 12184 (0.57 to 0.87%) |
<!-- /table:ar3c_skip -->

<!-- table:ar3c_skip_level -->
| rule | level | n | skipped | hurt | skipped_evals |
|---|---|---|---|---|---|
| hit | level 0 | 2170 | 2170 | 0 | 1,451,098 |
| hit | level 0.1 | 1921 | 1920 | 0 | 1,287,104 |
| hit | level 0.2 | 2097 | 2097 | 0 | 1,410,390 |
| hit | level 0.5 | 1459 | 1456 | 1 | 990,108 |
| hit | level 1 | 1658 | 1656 | 2 | 1,169,282 |
| cost | level 0 | 2170 | 2170 | 0 | 1,451,098 |
| cost | level 0.1 | 1921 | 1921 | 0 | 1,287,766 |
| cost | level 0.2 | 2097 | 2097 | 0 | 1,410,390 |
| cost | level 0.5 | 1459 | 1459 | 1 | 992,094 |
| cost | level 1 | 1658 | 1658 | 2 | 1,170,606 |
<!-- /table:ar3c_skip_level -->

## 4. Stop inside the fit

Reading the table. From the traces of 19,013 fits. `share`: share of a level's fit evaluations
after the fit is within 0.02 deg of its end; `start_within`: fits whose start already is. A fit has
a median 662 evaluations at every level.

<!-- table:ar3c_converged -->
| level | n | med_evals | share | med_share | start_within |
|---|---|---|---|---|---|
| level 0 | 4,774 | 662 | 0.999 | 0.998 | 4774 of 4774 (99.9 to 100.0%) |
| level 0.1 | 4,152 | 662 | 0.999 | 0.998 | 4152 of 4152 (99.9 to 100.0%) |
| level 0.2 | 4,332 | 662 | 0.998 | 0.998 | 4325 of 4332 (99.7 to 99.9%) |
| level 0.5 | 3,196 | 662 | 0.998 | 0.998 | 3186 of 3196 (99.4 to 99.8%) |
| level 1 | 2,559 | 662 | 0.998 | 0.998 | 2540 of 2559 (98.8 to 99.5%) |
| all | 19,013 | 662 | 0.998 | 0.998 | 18977 of 19013 (99.7 to 99.9%) |
<!-- /table:ar3c_converged -->

Two simulated rules at block ends: stop at the first failed block; stop when the best cost is at
most c (deciles of A's end costs). Picked on A, reported on B. A simulated stop is exact for its
own voxel; in a real run the shared generator would differ for later voxels (not simulated).

Both rules are uninformative here. The cost pick c = 0.3054 is the top of its grid (the 100th
percentile, A's largest end cost). "Stop at the first failed block" saves 98.1% because 18,500 of
the 19,013 fits never improve: all four blocks fail, so the rule ends each of them after block 1.

Reading the tables. The picked rules on B, then by level. `A_saved_share`: share of A's fit
evaluations saved; `B_saved`, `B_saved_share`: B's fit evaluations saved; `B_hurt`: B's fits hurt
by the stop; `B_med`, `B_p95`, `B_gt01`, `B_wrong`: the four error numbers of the resulting map
over all of B's voxels (12,184).

<!-- table:ar3c_stop -->
| rule | threshold | A_saved_share | B_saved | B_saved_share | B_hurt | B_med | B_p95 | B_gt01 | B_wrong |
|---|---|---|---|---|---|---|---|---|---|
| stop at the first failed block | none | 0.981 | 6,190,769 of 6,311,954 | 0.981 | 0 of 9305 (0.00 to 0.04%) | 0.0364 | 0.1376 | 1462 of 12184 (11.4 to 12.6%) | 86 of 12184 (0.57 to 0.87%) |
| stop when best cost <= c | 0.3054 | 0.981 | 6,190,989 of 6,311,954 | 0.981 | 0 of 9305 (0.00 to 0.04%) | 0.0364 | 0.1376 | 1462 of 12184 (11.4 to 12.6%) | 86 of 12184 (0.57 to 0.87%) |
<!-- /table:ar3c_stop -->

<!-- table:ar3c_stop_level -->
| rule | level | n | saved | total | hurt |
|---|---|---|---|---|---|
| first_fail | level 0 | 2170 | 1,422,888 | 1,451,098 | 0 |
| first_fail | level 0.1 | 1921 | 1,262,793 | 1,287,766 | 0 |
| first_fail | level 0.2 | 2097 | 1,383,085 | 1,410,390 | 0 |
| first_fail | level 0.5 | 1459 | 973,061 | 992,094 | 0 |
| first_fail | level 1 | 1658 | 1,148,942 | 1,170,606 | 0 |
| cost | level 0 | 2170 | 1,422,888 | 1,451,098 | 0 |
| cost | level 0.1 | 1921 | 1,262,793 | 1,287,766 | 0 |
| cost | level 0.2 | 2097 | 1,383,129 | 1,410,390 | 0 |
| cost | level 0.5 | 1459 | 973,127 | 992,094 | 0 |
| cost | level 1 | 1658 | 1,149,052 | 1,170,606 | 0 |
<!-- /table:ar3c_stop_level -->

## 5. Does a second fit help?

R2 on 4,768 voxels moved 9 by more than 0.05 deg and brought 4 closer; none pushed away.

Reading the tables. The counts with intervals, then by level.

<!-- table:ar3c_second_wil -->
| what | v |
|---|---|
| moved by more than 0.05 deg | 9 of 4768 (0.10 to 0.36%) |
| brought closer by more than 0.05 deg | 4 of 4768 (0.03 to 0.22%) |
| pushed away by more than 0.05 deg | 0 of 4768 (0.00 to 0.08%) |
<!-- /table:ar3c_second_wil -->

<!-- table:ar3c_second -->
| level | n | moved | helped | hurt |
|---|---|---|---|---|
| level 0 | 1206 | 0 | 0 | 0 |
| level 0.1 | 1052 | 0 | 0 | 0 |
| level 0.2 | 1088 | 1 | 1 | 0 |
| level 0.5 | 754 | 1 | 0 | 0 |
| level 1 | 668 | 7 | 3 | 0 |
| all | 4768 | 9 | 4 | 0 |
<!-- /table:ar3c_second -->

A threshold on the first fit's end hit ratio over the entry's best ("a second fit will not help"),
picked on A (lowest t with at most 1% of predicted voxels helped), reported on B. It sits at the
grid edge and says "will not help" for almost every voxel, so its 99.9% (2,288 of 2,290) is near
the base rate (`B_base`). Both B voxels a second fit helped were predicted not to need it. The
rule is not informative.

Reading the table. The picked rule on B. `A_pred`, `B_pred`: voxels the rule says a second fit
will not help; `A_helped`: of A's predicted voxels, those a second fit helped; `B_right`: B's
predicted voxels a second fit did not help; `B_base`: B's voxels a second fit helped.

<!-- table:ar3c_second_rule -->
| threshold | A_pred | A_helped | B_pred | B_right | B_base |
|---|---|---|---|---|---|
| 0.90 | 2471 | 2 | 2290 | 2288 of 2290 (99.7 to 100.0%) | 2 of 2292 (0.02 to 0.32%) |
<!-- /table:ar3c_second_rule -->

## 6. The BFS under gradients

Reading the table. D's fit log by voxel level (`fits`: every local fit of D, all phases; `same`:
returned the start bit for bit; `not_lowered`: the end's cost is not below the start's, the fit's
own test; `mv`: distance moved; `s_med`, `e_med`: median error of the start and end of the fit;
`e_wrong`: fits whose end is wrong; a wrong fit end is not a wrong voxel, since later visits
change it).

<!-- table:ar3c_dfits -->
| level | fits | evals | med_evals | mv_med | mv_p95 | same | not_lowered | s_med | e_med | e_p95 | e_wrong |
|---|---|---|---|---|---|---|---|---|---|---|---|
| level 0 | 8,617 | 5,695,834 | 662 | 0.0000 | 0.3226 | 4908 of 8617 (55.9 to 58.0%) | 7866 of 8617 (90.7 to 91.9%) | 0.0484 | 0.0484 | 53.2871 | 3846 of 8617 (43.6 to 45.7%) |
| level 0.1 | 11,893 | 7,875,695 | 662 | 0.0000 | 0.3705 | 6865 of 11893 (56.8 to 58.6%) | 10621 of 11893 (88.7 to 89.8%) | 0.1953 | 0.1954 | 54.0878 | 5930 of 11893 (49.0 to 50.8%) |
| level 0.2 | 18,608 | 12,369,093 | 662 | 0.0000 | 0.3983 | 10077 of 18608 (53.4 to 54.9%) | 16382 of 18608 (87.6 to 88.5%) | 20.7688 | 20.7151 | 54.3637 | 9843 of 18608 (52.2 to 53.6%) |
| level 0.5 | 22,516 | 15,059,166 | 662 | 0.0000 | 0.4481 | 11119 of 22516 (48.7 to 50.0%) | 18821 of 22516 (83.1 to 84.1%) | 0.4467 | 0.3998 | 53.4110 | 10847 of 22516 (47.5 to 48.8%) |
| level 1 | 28,558 | 19,320,043 | 662 | 0.0000 | 0.6240 | 12484 of 28558 (43.1 to 44.3%) | 21049 of 28558 (73.2 to 74.2%) | 0.6528 | 0.5059 | 53.8450 | 12963 of 28558 (44.8 to 46.0%) |
| all | 90,192 | 60,319,831 | 662 | 0.0000 | 0.5093 | 45453 of 90192 (50.1 to 50.7%) | 74739 of 90192 (82.6 to 83.1%) | 0.6252 | 0.4857 | 53.9183 | 43429 of 90192 (47.8 to 48.5%) |
<!-- /table:ar3c_dfits -->

`same` is a bitwise test, so it understates how often the fit does nothing. The fit did not lower
the cost on 74,739 of 90,192 fits (82.6 to 83.1%). Of these, 45,453 return the start bit for bit
and 29,286 return it up to the quaternion round trip (median 3.6e-15 deg).

Reading the table. D's fits by the error of their start: `lowered`, fits that ended below the
start's cost; `truth_lower`, fits whose voxel's true orientation has a lower strict cost than the
start (Reading 0).

<!-- table:ar3c_dfit_bins -->
| start_error | fits | lowered | truth_lower |
|---|---|---|---|
| 0 to 0.05 | 8,973 | 6 of 8973 (0.03 to 0.15%) | 7240 of 8973 (79.9 to 81.5%) |
| 0.05 to 0.1 | 8,240 | 62 of 8240 (0.59 to 0.96%) | 8213 of 8240 (99.5 to 99.8%) |
| 0.1 to 0.2 | 11,121 | 518 of 11121 (4.3 to 5.1%) | 11121 of 11121 (100.0 to 100.0%) |
| 0.2 to 0.5 | 14,762 | 3931 of 14762 (25.9 to 27.3%) | 14762 of 14762 (100.0 to 100.0%) |
| 0.5 to 1 | 3,493 | 2212 of 3493 (61.7 to 64.9%) | 3493 of 3493 (99.9 to 100.0%) |
| 1 to 5 | 167 | 91 of 167 (46.9 to 61.9%) | 167 of 167 (97.8 to 100.0%) |
| 5 and above | 43,436 | 8633 of 43436 (19.5 to 20.3%) | 43436 of 43436 (100.0 to 100.0%) |
<!-- /table:ar3c_dfit_bins -->

Up to 1 deg, the share of fits that lower the cost rises with the start's error, while the truth's
cost is lower than the start's in almost every fit above 0.05 deg. The fit succeeds once the start
error nears the step scale of the table above. This is consistent with the same step-size cause
as in R. The 43,436 fits that start 5 deg or more from the truth (31,522 of them revisits) are a
different case: they start from another grain's orientation.

Reading the tables. P plus R against D on the maps by level, then the four error numbers.
`P_wrong`, `PR_wrong`, `D_wrong`: wrong voxels of P, of P plus R and of D; `fixed`: wrong in D,
right in P plus R; `broken`: right in D, wrong in P plus R; `both`: wrong in both; `iv`: whether
the intervals of the wrong share of P plus R and of D overlap; `iv_gt01`: the same for the share
above 0.1 deg.

<!-- table:ar3c_map_vs_D -->
| level | n | P_wrong | PR_wrong | D_wrong | fixed | broken | both | iv | iv_gt01 |
|---|---|---|---|---|---|---|---|---|---|
| level 0 | 4,877 | 3 of 4877 (0.02 to 0.18%) | 3 of 4877 (0.02 to 0.18%) | 17 of 4877 (0.22 to 0.56%) | 16 | 2 | 1 | do not overlap | do not overlap |
| level 0.1 | 4,488 | 10 of 4488 (0.12 to 0.41%) | 10 of 4488 (0.12 to 0.41%) | 23 of 4488 (0.34 to 0.77%) | 16 | 3 | 7 | overlap | do not overlap |
| level 0.2 | 5,041 | 23 of 5041 (0.30 to 0.68%) | 23 of 5041 (0.30 to 0.68%) | 38 of 5041 (0.5 to 1.0%) | 23 | 8 | 15 | overlap | do not overlap |
| level 0.5 | 4,625 | 39 of 4625 (0.6 to 1.2%) | 39 of 4625 (0.6 to 1.2%) | 47 of 4625 (0.8 to 1.3%) | 22 | 14 | 25 | overlap | do not overlap |
| level 1 | 5,539 | 95 of 5539 (1.4 to 2.1%) | 95 of 5539 (1.4 to 2.1%) | 78 of 5539 (1.1 to 1.8%) | 14 | 31 | 64 | overlap | overlap |
| all | 24,570 | 170 of 24570 (0.60 to 0.80%) | 170 of 24570 (0.60 to 0.80%) | 203 of 24570 (0.72 to 0.95%) | 91 | 58 | 112 | overlap | do not overlap |
<!-- /table:ar3c_map_vs_D -->

<!-- table:ar3c_map_four -->
| group | n | median | p95 | gt01 | wrong |
|---|---|---|---|---|---|
| P plus R, all | 24,570 | 0.0342 | 0.1302 | 2531 of 24570 (9.9 to 10.7%) | 170 of 24570 (0.60 to 0.80%) |
| D, all | 24,570 | 0.0427 | 0.1459 | 3545 of 24570 (14.0 to 14.9%) | 203 of 24570 (0.72 to 0.95%) |
| P (start), all | 24,570 | 0.0342 | 0.1304 | 2541 of 24570 (10.0 to 10.7%) | 170 of 24570 (0.60 to 0.80%) |
| P plus R, level 0 | 4,877 | 0.0161 | 0.0518 | 3 of 4877 (0.02 to 0.18%) | 3 of 4877 (0.02 to 0.18%) |
| D, level 0 | 4,877 | 0.0169 | 0.0554 | 17 of 4877 (0.22 to 0.56%) | 17 of 4877 (0.22 to 0.56%) |
| P plus R, level 0.1 | 4,488 | 0.0313 | 0.0699 | 41 of 4488 (0.7 to 1.2%) | 10 of 4488 (0.12 to 0.41%) |
| D, level 0.1 | 4,488 | 0.0436 | 0.0964 | 179 of 4488 (3.5 to 4.6%) | 23 of 4488 (0.34 to 0.77%) |
| P plus R, level 0.2 | 5,041 | 0.0401 | 0.0897 | 152 of 5041 (2.6 to 3.5%) | 23 of 5041 (0.30 to 0.68%) |
| D, level 0.2 | 5,041 | 0.0533 | 0.1167 | 554 of 5041 (10.2 to 11.9%) | 38 of 5041 (0.5 to 1.0%) |
| P plus R, level 0.5 | 4,625 | 0.0570 | 0.1374 | 770 of 4625 (15.6 to 17.7%) | 39 of 4625 (0.6 to 1.2%) |
| D, level 0.5 | 4,625 | 0.0665 | 0.1618 | 1251 of 4625 (25.8 to 28.3%) | 47 of 4625 (0.8 to 1.3%) |
| P plus R, level 1 | 5,539 | 0.0550 | 0.1873 | 1565 of 5539 (27.1 to 29.5%) | 95 of 5539 (1.4 to 2.1%) |
| D, level 1 | 5,539 | 0.0535 | 0.1960 | 1544 of 5539 (26.7 to 29.1%) | 78 of 5539 (1.1 to 1.8%) |
<!-- /table:ar3c_map_four -->

Under the plan's rule (a gain only when the intervals do not overlap), P plus R beats D on the
share of voxels above 0.1 deg: 2,531 (9.9 to 10.7%) against 3,545 (14.0 to 14.9%) overall, and at
levels 0.1, 0.2 and 0.5. On wrong voxels it wins only at level 0 (3 against 17); overall the
intervals overlap (170 against 203). The gain is P's: P alone has 2,541 voxels above 0.1 deg and
the same wrong counts, so R adds nothing. P plus R costs 90,985,770 evaluations against D's
118,350,225. This is one run of each arm.

Scale: two key-off pilot runs of step 2 differ by 20 voxels one way and 7 the other, of 2,000; the
run-to-run difference on the whole layer is not measured.

## 7. The control

Reading the table. Level-0 grains (4,877 voxels) in P and D here against the same voxels in step
2's C and D. Neighbours differ, so this is a reading, not a repeat.

<!-- table:ar3c_control -->
| group | n | median | p95 | gt01 | wrong |
|---|---|---|---|---|---|
| this run P (paint) | 4,877 | 0.0161 | 0.0518 | 3 of 4877 (0.02 to 0.18%) | 3 of 4877 (0.02 to 0.18%) |
| this run P plus R | 4,877 | 0.0161 | 0.0518 | 3 of 4877 (0.02 to 0.18%) | 3 of 4877 (0.02 to 0.18%) |
| this run D (BFS) | 4,877 | 0.0169 | 0.0554 | 17 of 4877 (0.22 to 0.56%) | 17 of 4877 (0.22 to 0.56%) |
| step 2 C (paint) | 4,877 | 0.0177 | 0.0486 | 1 of 4877 (0.00 to 0.12%) | 1 of 4877 (0.00 to 0.12%) |
| step 2 D (BFS) | 4,877 | 0.0194 | 0.0486 | 37 of 4877 (0.6 to 1.0%) | 37 of 4877 (0.6 to 1.0%) |
<!-- /table:ar3c_control -->

## 8. Wall time

Reading the table. P and D ran side by side with a busy editor helper process. R and R2 ran
later in a separate refine-stage process, with no preflight, while D was still running (R
throughout, R2 at its start). Every wall time is contended and no timing comparison is claimed.

<!-- table:ar3c_wall -->
| run | wall | label |
|---|---|---|
| P (paint, with set-up of pairing) | 7390 | contended: timing preflight refused (cpptools at 100% CPU, Python at 99% CPU); single run, one worker |
| R (one fit per painted voxel) | 841 | contended: separate refine-stage process, no preflight; run D running |
| R2 (second fit, every fourth voxel) | 219 | contended: separate refine-stage process, no preflight; run D running at start |
| D (BFS) | 9369 | contended: timing preflight refused (cpptools at 114% CPU, Python at 100% CPU); single run, one worker |
<!-- /table:ar3c_wall -->

## What is and is not claimed

Claimed: on this gradient sample, at these settings, the painted start is close to the truth
(consistent with the many seeds, 21.9% of voxels; untested, no arm with fewer seeds), and the BFS
local fit moved almost none of the 19,013 painted voxels. The stored records show that the strict
cost does tell the truth from the start (lower at the truth on 16,591 of 19,013, and on all 2,240
starts more than 0.1 deg off), yet the fit lowered the cost on only 36. The fit's first-block
trial steps (median 0.70 deg, 5th percentile 0.32 deg) are many times the start's error, and its
restarts land 4.9 to 16.8 deg (median) from the start. Consistent with a step too coarse for this
start; a run with a finer box would test it (step 4). Step 2's sample showed the same inertness.
P (with R) beats D on the share above 0.1 deg, in one run of each arm.

Not claimed: that the gradient model matches a real sample (the owner is asked which sizes and
shapes are typical); sharp sub-grain boundaries; realistic images, another Q-max, layer or voxel
order; a tested stop rule (the rules are simulated from records, at grid edges); a tested wall
time; another local optimizer (the CMA key exists, not run); the run-to-run difference of one arm
on the layer.

Deviations: the 90-minute stop on P and the 150-minute stop on D were not coded (the constants
were defined and unused, now removed). P ran 7,390 s (output is written only at the end, so it was
left to finish, under twice the budget); D 9,369 s. The total compute was about 4.95 h against the
plan's 4-hour budget: P 2.05 h, D 2.60 h, R and R2 0.29 h. R and R2 ran as a stage loading P's
stored paint and generator state; the P map outside painted voxels is read from `recon.mic`
(Euler angles to 1e-4 deg). D's map is stored gzipped (`recon.mic.gz`).
