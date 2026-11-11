# AR4 step 1: a learned verifier, measured offline against the hit ratio

Step 1 of "AR4" in `docs/autoresearch/PLAN.md`, task `ar4a`. No file under `icenine/` changed. No
key was added and no default changed. Scripts are in `scripts/autoresearch/ar4/` (`verify_*.py`).
Results are in `benchmarks/autoresearch/ar4/`. Every table below is filled from `tables.md`
there, and every number in the text comes from `summary.json` or `thresholds.json` there.

The question: can a small model, fitted once on two samples, tell a right orientation from a wrong
one at a voxel better than the hit ratio does? The hit ratio is the only test today. It sees one
number per evaluation: the share of predicted pixels that fall on observed intensity. The model
sees the same evaluation split by predicted peak. Below, `N1` is the learned set network and
`B0` the hit ratio alone; the section "Terms" defines these and the other terms.

**Short answer.** On `mid`, no model passes the plan's rule for use 1 (flag the unfinished
voxels). `N1` flags 23 of the 31 wrong voxels; the rule needs as many as the hit ratio flags
(31).

Use 2 is to accept a seed without a grid search. On the 15 grid-fallback cases of `mid`, `N1`
accepts 7 of the 10 right ends and 0 of the 5 wrong ends. The hit ratio `B0` accepts none of
them, but it cannot do otherwise: a case is a fallback because its hit ratios are under 0.80,
and the cut of `B0` is 0.926. So this comparison with `B0` says nothing about `N1`.

Use 2 can be read in two ways, and the reading decides whether it pays:

- **OR rule** (`N1` is added to the seed test: a seed is accepted when the hit ratio or `N1`
  says so). Nothing is refused. The 7 accepted right ends skip 2,325,301 of the 5,002,750
  grid evaluations of the 15 cases.
- **Replace rule** (`N1` replaces the seed test). `N1` also refuses 20 right passing seeds. At
  the mean grid search of the 15 cases that is about 6,670,333 extra evaluations, more than the
  2,325,301 saved.

Every evaluation count in these two lines is arithmetic from stored numbers, not a run. Only the
OR rule is worth a test. Predictions 1 and 3 agree with the numbers; predictions 2 and 4 do not.

The plot below leads because it shows the effect that decides use 1. Reading it: each colour is
one score. Each dot is one cut that flags a given number of the 31 wrong voxels of `mid`
(x axis) at the smallest cost in flagged voxels; the y axis is the share of all 24,570 voxels
that this cut flags. Only measured cuts are drawn, with no line between them. Dots near the
bottom are a cheap flag. The open star is the cut that was
picked on the train samples (before `mid` was touched); it is the point the pass rule judges.
The star of `N1` is at 23 of 31 wrong voxels, left of where it would have to be (31).

![Flag curve on mid](../../benchmarks/autoresearch/ar4/ar4_flag_curve.png)

## Terms

- A **pair** is one voxel with one orientation.
- A pair is **right** if its orientation is within 1 degree of the voxel's truth, by the library's
  symmetry-reduced angle. Nothing here is specific to one crystal system. The peak-family input is
  a rank by the length of the reciprocal vector, not Miller indices.
- A **verifier** reads numbers from one evaluation of a pair (no truth) and returns a score from
  0 to 1 that the pair is right. A pair is **flagged** when its score is at or under a **cut**.
- **Local optimization** is the refinement that follows a candidate orientation.
- A **seed step** is one call of the library's `index_seed`: index the voxel from its rays, cost
  20 candidates, refine up to 5. If the best end fails the seed test (hit ratio 0.80), the voxel
  goes to a full **grid search**. These are the **grid-fallback** cases.
- The **hit ratio** is the share of predicted pixels on observed intensity. `S8` is the stored run
  with the boundary pass; "after the pass" is its final map.
- `B0` is the hit ratio alone (today's rule). `B1` is a logistic regression on six summary
  numbers: hit ratio, cost, confidence, the share of peaks with any lit pixel, the lowest quarter
  of the per-peak shares, and the number of attached rays the orientation explains. `N1` is a
  set network: one shared layer for each predicted peak, then the mean and the maximum over the
  peaks, then two layers, with the summary numbers added. It has 21,953 parameters. `lit` is
  the share of predicted peaks with any lit pixel, with no training.
- `N1` reads five numbers per predicted peak: the peak-family rank; "dist"; the rotation angle
  (omega); the share of the peak's pixels on observed intensity; and the peak's pixel count.
  Two of them differ from the plan's list. "dist" is the radial position of the spot on
  detector 0, as a fraction of the image size, not the detector distance the plan names. The
  pixel count was not in the plan's list; it was added.
- The **OR rule** and the **replace rule** are the two readings of use 2 (section "Use 2").

## What was run

Train samples: `big` and the first gradient sample (`grad`). Test sample: `mid`. No model was
trained and no cut was picked on `mid`. The pairs are of eight kinds (plan, "Pairs"): every voxel
after the pass and, where it differs, before it; every distinct neighbour orientation at each
boundary voxel (`try_start`); the refined ends of every stored seed step (`seed_end`); the costed
candidates of every 10th seed step (`seed_cand`); the truth turned by 0.1, 0.3 and 0.5 degree
about a random axis at every 20th voxel.

Reading the next table: one row per sample, the number of pairs of each kind. Each pair costs one
evaluation of the strict cost function. The three samples hold 229,651 pairs and so 229,651
evaluations.

<!-- table:ar4_pairs -->
| sample | pairs | after | before | try_start | seed_end | seed_cand | turn_01 | turn_03 | turn_05 |
|---|---|---|---|---|---|---|---|---|---|
| big | 49103 | 24570 | 1056 | 18665 | 465 | 660 | 1229 | 1229 | 1229 |
| grad | 94306 | 24570 | 752 | 60381 | 2356 | 2560 | 1229 | 1229 | 1229 |
| mid | 86242 | 24570 | 1429 | 50099 | 3217 | 3240 | 1229 | 1229 | 1229 |
<!-- /table:ar4_pairs -->

The library's cost function returns only the totals of one evaluation. The script takes the
per-peak pixel counts from the same overlap data, and the library is not changed. Each pair
takes two library calls. The first is the fused evaluation that the stored runs use; it gives
the score of the pair. The second runs the non-fused path (torch, then the stage-D C function
once per predicted peak); it gives the per-peak counts. Checks on all 229,651 pairs: the
per-peak counts sum exactly to the totals of the non-fused path (largest difference 0 pixels).
C++ also counts the pixels of each predicted peak and adds them (`Src/OverlapInfo.h`), so the
split itself matches C++.

Against the fused evaluation the per-peak sums differ on 2,271 of the 86,242 pairs of `mid`, by
at most 6 pixels. The hit ratios differ by at most 0.0047 (0.0039 on `big`, 0.0043 on the first
gradient sample). This is consistent with the known pixel-edge differences between the fused and
the torch paths (parity rule D-20261027-I2); it was not traced pair by pair. One flipped pixel
moves a peak's share by 1/17 at the median peak size of 17 pixels, and by 1/6 at the smallest
peak (6 pixels).

Reading the next table: per sample, how many right and wrong pairs differ. The shares are about
equal for right and wrong pairs, so no difference by label is detectable, and train and test
samples are alike. So the difference cannot favour one model over another.

<!-- table:ar4_sumback -->
| sample | right | wrong |
|---|---|---|
| big | 932 of 37196 (2.51%) | 326 of 11907 (2.74%) |
| grad | 956 of 36607 (2.61%) | 1618 of 57699 (2.80%) |
| mid | 1373 of 52188 (2.63%) | 898 of 34054 (2.64%) |
<!-- /table:ar4_sumback -->

The replay of the seed steps reproduced the stored local optimization evaluations and hit ratios
in all 329, 1,279 and 1,611 steps of `big`, the first gradient sample and `mid` (0 mismatches).

One planned summary number is missing: `best_votes`, the best vote count among the rays an
orientation does not explain. A pilot on 40 voxels of `big` put it at 50 to 200 ms a pair. That
is an unrecorded estimate: the pilot's output was not kept. At that cost it is hours for this
many pairs, over the 45-minute budget, so it was not computed and no model uses it. AR5 step 1
stored vote counts, but only for 2,017, 3,803 and 4,420 voxels of the 24,570 per sample (the
voxels its own flag picked plus every 20th other boundary voxel), and only for the orientation
after the pass. The other pairs have none, so the stored counts could not feed a model trained
on all pairs. The number of attached rays that the orientation explains was computed and is
used. The inputs took about 17 minutes for the two train samples (run side by side) and 31
minutes for `mid`, with two processes at most and other jobs running (contended; not timing
results).

## Models and cuts (train samples only)

Each model was trained once on all 143,409 pairs of the two train samples, with generator seed 0.
`N1` ran 20 epochs in 125 s. The layer sizes, loss and epochs were committed before the run. The
cut of each score is the one that flags the fewest voxels while flagging at least 95% of the 87
wrong voxels (24 and 63) of the maps after the pass, the two samples pooled. Ties at the cut are
all flagged. The thresholds file was committed (`669b9a3`) before any script read `mid`. The
cuts are picked on the same voxels the models trained on, so the train numbers are in-sample.

Reading the next table: one row per score. `cut` is the picked cut; `flagged` is the number of the
49,140 train voxels flagged; `wrong` is how many of the 87 wrong voxels that includes.

<!-- table:ar4_cuts -->
| score | cut | flagged | wrong | share |
|---|---|---|---|---|
| B0 (hit ratio) | 0.926 | 1741 | 83 | 3.54% |
| B1 (logistic, summary numbers) | 0.976 | 25009 | 83 | 50.89% |
| N1 (set network) | 0.749 | 110 | 83 | 0.22% |
| share of peaks with a lit pixel | 1.000 | 49140 | 87 | 100.00% |
<!-- /table:ar4_cuts -->

The `lit` cut is 1.000: 20 of the 87 wrong voxels have a lit share under 1, so no cut under 1.0
flags 83 of the 87, and the cut flags every voxel.

## Use 1: flag the unfinished voxels on `mid`

Reading the next table: one row per score, on the 24,570 voxels of the `S8` map after the pass.
There are 31 wrong voxels. `wrong` is how many of them the train cut flags, with a Wilson 95%
interval (a range of shares that fits the count). `flagged` is all voxels flagged. `auc` (area
under the curve) is the share of (right, wrong) voxel pairs in which the right voxel scores
higher (a tie counts one half); it is a count of ordered pairs, not a test. `rule` is the
plan's pass rule: flag at least as many of the 31 as `B0` does, and at most half as many voxels
as `B0` does (880). The rule compares with `B0`, so `B0`'s cell reads
"baseline".

<!-- table:ar4_use1_mid -->
| score | wrong | flagged | share | auc | rule |
|---|---|---|---|---|---|
| B0 (hit ratio) | 31 of 31 (0.89 to 1.00) | 1760 of 24570 | 7.16% | 0.981 | baseline |
| B1 (logistic, summary numbers) | 31 of 31 (0.89 to 1.00) | 13653 of 24570 | 55.57% | 0.875 | no |
| N1 (set network) | 23 of 31 (0.57 to 0.86) | 190 of 24570 | 0.77% | 0.994 | no |
| share of peaks with a lit pixel | 31 of 31 (0.89 to 1.00) | 24570 of 24570 | 100.00% | 0.506 | no |
<!-- /table:ar4_use1_mid -->

`B0` flags all 31 and 1,760 voxels. `N1` flags 190 voxels, far under half, but only 23 of the 31
wrong voxels: it fails the first half of the rule. `B1` flags all 31 but 13,653 voxels. No model
passes. The same cuts on the train samples (in-sample) flag 21 and 22 of the 24 wrong voxels of
`big` for `B0` and `N1`, and 62 and 61 of the 63 of the first gradient sample.

`N1` misses 8 of the 31 wrong voxels of `mid`. On the train samples its cut misses 4 of the 87
wrong voxels, by construction: the cut was picked there to flag at least 95% of them. The gap is
consistent with the cut being picked on the voxels `N1` trained on; that cause is untested.

The 8 missed voxels are boundary voxels, but that does not single them out: all 31 wrong voxels
of `mid` are boundary voxels. Their errors are 34 to 56 degrees, their hit ratios 0.889 to 0.910
and their `N1` scores 0.83 to 0.99.

Reading the next table: one row per missed voxel; `err` in degrees, `hit` the hit ratio, `n1` the
`N1` score (the cut is 0.749), `bnd` whether the voxel has a try in the boundary pass.

<!-- table:ar4_missed -->
| voxel | err | hit | n1 | bnd |
|---|---|---|---|---|
| 97 | 55.06 | 0.903 | 0.986 | yes |
| 2995 | 34.40 | 0.908 | 0.959 | yes |
| 3302 | 54.20 | 0.889 | 0.912 | yes |
| 4758 | 52.10 | 0.904 | 0.982 | yes |
| 5807 | 38.94 | 0.905 | 0.828 | yes |
| 6419 | 52.50 | 0.910 | 0.899 | yes |
| 8427 | 55.99 | 0.899 | 0.973 | yes |
| 21630 | 54.73 | 0.900 | 0.926 | yes |
<!-- /table:ar4_missed -->

A what-if, not a result: the same scores with the cut picked on `mid` itself. Reading the next
table: `p95` is the cut that flags 30 of the 31 wrong voxels (95%), `p100` the cut that flags all
31; each cell is the voxels flagged and their share of all voxels. This is hindsight; the cut
moves between samples (the cut of `B0` picked on `mid` flags fewer voxels than its train cut).

<!-- table:ar4_whatif -->
| score | p95 | p100 |
|---|---|---|
| B0 (hit ratio) | 1016 (4.14%) | 1029 (4.19%) |
| B1 (logistic, summary numbers) | 9548 (38.86%) | 11910 (48.47%) |
| N1 (set network) | 779 (3.17%) | 854 (3.48%) |
| share of peaks with a lit pixel | 24570 (100.00%) | 24570 (100.00%) |
<!-- /table:ar4_whatif -->

## Boundary voxels and the share of lit peaks (prediction 4)

A boundary voxel is a voxel with at least one try in the boundary pass: 17,859 of the voxels of
`mid`, holding all 31 wrong ones. Reading the next table: the same scores on those voxels alone.
`auc` as above; `wrong` and `flagged` at the train cuts.

<!-- table:ar4_boundary -->
| score | auc | wrong | flagged |
|---|---|---|---|
| B0 (hit ratio) | 0.974 | 31 of 31 | 1735 of 17859 |
| B1 (logistic, summary numbers) | 0.882 | 31 of 31 | 9487 of 17859 |
| N1 (set network) | 0.991 | 23 of 31 | 189 of 17859 |
| share of peaks with a lit pixel | 0.496 | 31 of 31 | 17859 of 17859 |
<!-- /table:ar4_boundary -->

The share of lit peaks does not separate right from wrong here (`auc` 0.496 against 0.974 for the
hit ratio). Reading the next table: after the pass, the share of voxels whose predicted peaks are
all lit, the number of wrong voxels with all peaks lit, and the lowest and median lit share of the
wrong voxels.

<!-- table:ar4_lit -->
| sample | all_lit | wrong_all_lit | lowest | median |
|---|---|---|---|---|
| big | 77.1% | 19 of 24 | 0.972 | 1.000 |
| grad | 69.8% | 48 of 63 | 0.900 | 1.000 |
| mid | 79.1% | 24 of 31 | 0.984 | 1.000 |
<!-- /table:ar4_lit -->

On these clean images most predicted peaks of a wrong orientation still touch some observed
intensity (consistent with dense clean images; not tested with a sparser image). At the boundary
voxels of `mid`, every peak is lit for 77.1% of the right voxels and 77.4% of the wrong ones.
The lowest lit share is 0.906 for a right voxel and 0.984 for a wrong one. With a median of 65
peaks a pair, one unlit peak moves the share by about 0.015, so the share takes few values: 76
distinct values for the right voxels and 5 for the wrong ones. So the `auc` of 0.496 is a
measured result, not a constant input. The count of peaks with a lit pixel carries almost no
signal here. `N1` scores higher than the hit ratio on `auc` (0.991 against 0.974), by a small
margin.

## Use 2: accept a seed without a grid search

The seed steps whose best refined end failed the seed test are the 31 grid-fallback cases of the
plan (16 on the first gradient sample, 15 on `mid`; the voxels are the same as in AR5 step 1). The
"end" of a step is the one the library takes (the first end that reaches hit ratio 0.95, else the
lowest final cost). A seed is accepted when its score is at or above the cut of use 1. As a
control, the same scores are applied to the seed steps that passed (1,596 on `mid`).

Reading the plot: each dot is the score of the end of one of the 15 fallback seed steps of `mid`.
Green dots are right ends (10), red dots wrong ends (5). The dashed line is the cut. A dot above
the line is accepted.

![Scores of the fallback ends](../../benchmarks/autoresearch/ar4/ar4_use2.png)

The grid evaluations below are those of the grid searches alone. AR5 step 1's table gives
larger totals because, for a reseed case, its field is the whole fit from the seed: the seed
step plus the grid search. The seed step runs before any accept decision, so it is never saved.
Reading the next table: per sample, the number of fallback cases and of reseed cases, AR5's
total, the seed-step evaluations of the reseed cases (from the stored log of `S8`), the
grid-only total, and the library's own counter of grid-fallback evaluations in `S8`'s run
file. The grid-only total equals the counter on both samples.

<!-- table:ar4_grid_evals -->
| sample | cases | reseed | ar5_field | seed_steps | grid_only | counter |
|---|---|---|---|---|---|---|
| grad | 16 | 15 | 5,518,587 | 69,086 | 5,449,501 | 5,449,501 |
| mid | 15 | 14 | 5,087,764 | 85,014 | 5,002,750 | 5,002,750 |
<!-- /table:ar4_grid_evals -->

Reading the next table: `right` and `wrong` are the ends accepted, of the 10 right and 5 wrong
ends, with a Wilson 95% interval. `saved` is a what-if, not a result: the stored grid-only
evaluations of the grid searches that accepting the right ends would have skipped. It does not
count the scoring call (one more library call per failed seed, 15 on `mid`; not measured) or the
cost of refused passing seeds. `refused` is the passing seeds that the cut would refuse, and
how many of those ends are right.

<!-- table:ar4_use2_mid -->
| score | right | wrong | saved | refused |
|---|---|---|---|---|
| B0 (hit ratio) | 0 of 10 (0.00 to 0.28) | 0 of 5 (0.00 to 0.43) | 0 | 141 of 1596 (130 right) |
| B1 (logistic, summary numbers) | 3 of 10 (0.11 to 0.60) | 0 of 5 (0.00 to 0.43) | 999,081 | 582 of 1596 (571 right) |
| N1 (set network) | 7 of 10 (0.40 to 0.89) | 0 of 5 (0.00 to 0.43) | 2,325,301 | 29 of 1596 (20 right) |
| share of peaks with a lit pixel | 2 of 10 (0.06 to 0.51) | 4 of 5 (0.38 to 0.96) | 659,324 | 476 of 1596 (472 right) |
<!-- /table:ar4_use2_mid -->

These are small counts and are read as counts. On `mid`, `N1` accepts 7 of the 10 right ends and
none of the 5 wrong ones. Its 7 skip stored grid searches of 2,325,301 evaluations, out of
5,002,750 for the 15 cases. It refuses 29 of the 1,596 passing seeds, 20 of them right and 9
wrong (11 of the 1,596 passing ends are wrong). `B0` accepts none of the 15 fallback ends by
construction: their hit ratios are under 0.80, so none can reach its cut of 0.926. At that cut
it would refuse 141 passing seeds, 130 of them right. `lit` accepts 4 of the 5 wrong ends.

**Which rule use 2 means.** The plan says "a seed is accepted when its score is at or above the
threshold" and counts the refused passing seeds as a control, but does not say whether the
score replaces the seed test or is added to it. The two readings give opposite answers:

- **OR rule:** a seed is accepted when the seed test passes (hit ratio 0.80) or the score is at
  or above the cut. No passing seed is refused. For `N1` on `mid`: 7 of 10 right ends and 0 of
  5 wrong ends accepted; 2,325,301 grid evaluations skipped, 46% of the 5,002,750.
- **Replace rule:** the score alone decides. The 29 passing seeds that `N1` refuses then each
  go to a grid search. The 20 right ones are pure cost: at the mean grid search of the 15 cases
  (333,517 evaluations) that is about 6,670,333 extra evaluations, against 2,325,301 saved. The
  9 wrong ones would cost about 3,001,650 more; a grid search might fix them, which is untested.

Both are arithmetic from stored numbers, not a run, and the mean grid search is an assumed
price for a search that was not run. Use 2 is usable only as an OR rule.

Reading the next table: one row per score on `mid`, both readings. `or_saved` is the grid
evaluations the accepted right ends skip under the OR rule, and `or_wrong` the wrong ends
accepted. `rep_refused` is the passing seeds refused under the replace rule, right and wrong.
`rep_extra_right` and `rep_extra_wrong` price those refusals at 333,517 evaluations each. All
evaluation counts in it are arithmetic, not a run.

<!-- table:ar4_use2_rules -->
| score | or_saved | or_wrong | rep_refused | rep_extra_right | rep_extra_wrong |
|---|---|---|---|---|---|
| B0 (hit ratio) | 0 | 0 of 5 | 130 right, 11 wrong | 43,357,167 | 3,668,683 |
| B1 (logistic, summary numbers) | 999,081 | 0 of 5 | 571 right, 11 wrong | 190,438,017 | 3,668,683 |
| N1 (set network) | 2,325,301 | 0 of 5 | 20 right, 9 wrong | 6,670,333 | 3,001,650 |
| share of peaks with a lit pixel | 659,324 | 4 of 5 | 472 right, 4 wrong | 157,419,867 | 1,334,067 |
<!-- /table:ar4_use2_rules -->

The same table for the first gradient sample, where `N1`'s training pairs included these ends
(in-sample, so not evidence):

<!-- table:ar4_use2_grad -->
| score | right | wrong | saved | refused |
|---|---|---|---|---|
| B0 (hit ratio) | 0 of 13 (0.00 to 0.23) | 0 of 3 (0.00 to 0.56) | 0 | 116 of 1263 (94 right) |
| B1 (logistic, summary numbers) | 6 of 13 (0.23 to 0.71) | 0 of 3 (0.00 to 0.56) | 2,041,128 | 356 of 1263 (337 right) |
| N1 (set network) | 11 of 13 (0.58 to 0.96) | 0 of 3 (0.00 to 0.56) | 3,744,760 | 31 of 1263 (10 right) |
| share of peaks with a lit pixel | 4 of 13 (0.13 to 0.58) | 1 of 3 (0.06 to 0.79) | 1,363,405 | 347 of 1263 (340 right) |
<!-- /table:ar4_use2_grad -->

## Predictions

1. `B1` does not pass use 1: **agrees.** It flags all 31 wrong voxels but 13,653 voxels (55.57%).
   This tests a reduced `B1`: it lacks `best_votes`, one of the plan's summary numbers (see
   "What was run" for why).
2. `N1` passes use 1 on `mid`: **does not agree.** It flags 23 of 31 (0.57 to 0.86).
3. `N1` accepts at least two thirds of the 10 right fallback ends and at most one wrong end:
   **agrees.** 7 of 10 (0.40 to 0.89), 0 of 5. Seven is just over the two-thirds line (6.7).
   The prediction does not say which rule use 2 means; it pays only as an OR rule.
4. The share of lit peaks separates right from wrong better than the hit ratio: **does not
   agree.** `auc` 0.496 against 0.974.

## Scope and limits

- Simulated samples with clean images only. **Realistic images (overlap and detector noise) are
  untested**, so nothing here says how a verifier behaves on them.
- One run per model, generator seed 0, no other seed. Counts and Wilson intervals; no p-values.
  Voxels in a grain are not independent, so the intervals are too narrow across grains.
- The numbers 31 (wrong voxels of `mid`) and 15 (its fallbacks) were known from AR5 step 1 and the
  plan before the run. The order of work for `mid` was: model definitions, training on the train
  samples and the thresholds file committed (`669b9a3`, 12:43 on 11 Nov 2026); the first
  script that read `mid` started at 12:44. Before that, no script touched `mid`, not even to fill
  a cache or test a loader. (Links to `mid`'s stored files were made just before the run.)
- The cuts were picked in-sample (same voxels as the training pairs); the test on `mid` is not.
  `N1` misses 4 of 87 wrong train voxels and 8 of 31 on `mid` (see "Use 1").
- `best_votes` was not computed (see above). A model with it is untested.
- The result is for one set of 31 wrong voxels, 8 of which `N1` misses. I describe the misses; I
  did not test why `N1` scores them high.
- The use 2 savings and costs are arithmetic from stored numbers. No whole-layer run tested
  either rule.
- Wall times were measured with other jobs running and are not timing results.

## What it means for the next step

By the plan's rule step 2 (wire a verifier in as a script layer and run the three samples) is run
only for a model that passes use 1. None did. Use 2 alone is a separate question for the owner.
As an OR rule, `N1` accepted 7 of 10 right fallback ends and 0 of 5 wrong ones on `mid`, and by
arithmetic would skip 2,325,301 grid evaluations. That rests on 15 cases, with wide intervals,
and with a cut picked for use 1. As a replace rule it would cost more than it saves.
