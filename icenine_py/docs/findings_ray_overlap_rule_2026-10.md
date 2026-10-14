# Ray-consistency overlap rule: results (2026-10-11 to 2026-10-12)

Task branch `feature/finisher-mc-study-ray-overlap` (parent `feature/finisher-mc-study`). Commits 8d35ff0, 362799e
(rule, steps 1-2), f6ad076 (case f rejection dropped, ray check in C), 5f0cd9c (window decision), 387f3b6 (30-voxel
coarse test and end-to-end test). Raw tables: `benchmarks/phase_d_seed_diag/ray_overlap/` (`tables.md`,
`tables30.md`, `e2e_tables.md`); run details in its `RUN_LOG.md`.

**Update 2026-10-13 (centroid task, branch `feature/finisher-mc-study-centroid`).** The lenient rule's screen used to
project the first triangle vertex instead of the voxel centroid that C++ and the other rules use. Every rule now
starts from the centroid. The clean-image numbers of the lenient rule in this write-up (sections 3 to 6, the summary
table and the decisions) were re-run with that start point and replaced; the previous values are named in section
10, which also holds the before / after tables. The realistic-image numbers of the lenient rule were not re-run
(owner, clean images only) and still come from the vertex-0 screen. The C++ rule and both ray windows are
unchanged (re-run: same results).

## Summary

On clean images (the images to rely on; see the caveat below), the ray rule recovers nearly all the accuracy that the
strict C++ overlap rule loses: 40 and 41 wrong of 200 for the two windows, against 56 for the C++ rule and 41 for the
original lenient Python rule (39 before the lenient rule's start point was changed to the centroid, section 10).
Paired by voxel, the 2-pixel window against the lenient rule is 10 fixed, 9 broken and 31 both wrong, so the two are
tied on these 200 voxels. The ray rule costs 3 to 4% more cost evaluations than the C++
rule per E0 voxel, and on full-sample voxels 6% (2-pixel window) and 21% (geometric window) more than the C++ rule,
against 46% more for the lenient rule (re-run in the current build, section 6).

Caveat on the realistic-image columns: the realistic images use a noise model the owner has ruled unphysical (every
spot on each detector is dropped at random with probability 0.1; open question 8 of the parent findings). On those
images the lenient rule is more accurate (18 wrong of 200 against 29 to 32 for the ray rule), but this is reported for
completeness and should not drive decisions.

Update 2026-10-14 (Q8, `MIGRATION_HISTORY.md` "Q8: realistic images without the random spot drop"): the same E0 runs
were repeated on images without the random spot drop. With `min_sin_eta` 0 the wrong counts of 200 are, clean / old
realistic / no drop: lenient 41 / 25 / 22, C++ rule 56 / 56 / 38, 2-pixel window 40 / 32 / 14, geometric window 41 / 29 / 14
(the earlier 18 for the lenient rule was the vertex-0 start; it is 25 in this build; the lenient rule with the vertex-0 screen
reproduces the 18, `MIGRATION_HISTORY.md` "Q8" Step 8). On the no-drop images the ray rule is
ahead of the lenient rule (2-pixel window against lenient: 15 fixed, 7 broken, 7 both wrong), so the lenient lead above
comes with the drop; it is one sample of 200 voxels and one seed. Mean evaluations per voxel on the no-drop images: lenient 51,901,
2-pixel window 51,652, geometric window 51,677, C++ rule 49,281 (clean: 48,700, 48,641, 48,649, 47,230).

| Overlap rule | Wrong of 200 (clean images) | Wrong of 200 (realistic images, see caveat) | Time per full-sample voxel (one process, mixed quiet and contended, three sessions) |
|---|---|---|---|
| Original lenient Python rule | 41 | 18 (vertex-0 start point, not re-run) | 180 s (six processes at once; 176 s before) |
| Strict C++ rule | 56 | 56 | 116 s |
| Ray rule, 2-pixel window | 40 | 32 | 131 s |
| Ray rule, geometric window | 41 | 29 | 152 s |

Each arm is one run per voxel (seed 0); the randomness only drives the search's exploration, so the comparison rests on the paired per-voxel counts. The two ray-rule windows cannot be told apart on these numbers. The times are means over
six full-sample voxels; the preflight was quiet for 0 of 6 runs (lenient, re-run 2026-10-13 as six processes at once), 5 of 6 (C++), 2 of 6 (2-pixel window) and
5 of 6 (geometric window), so compare the rules on evaluations (section 6) rather than on these seconds.

## 1. What the rule does

When a predicted diffraction spot lands on lit pixels on the near detector (detector 0), the cost function must
decide whether that counts as a matched peak. The study compares three ways to decide:

| Overlap rule | How it decides whether a lit spot on detector 0 counts |
|---|---|
| Original lenient Python rule | Counts almost every overlap. |
| Strict C++ rule | The C++ `CountQualifiedPeaks` rule: the spot counts only if the same orientation's predicted spot on the far detector (detector 1) is also lit. A predicted spot that misses detector 1 altogether is forgiven (counted as valid). |
| Ray rule (owner's design) | Traces the observed ray from the voxel through the lit pixels on detector 0 out to detector 1, and counts the spot only if detector 1 is lit near where that ray lands. |

"Near where the ray lands" needs a tolerance window on detector 1. Two windows were carried forward:

| Window | Definition |
|---|---|
| 2-pixel window | A square of 2 pixels on each side of the pixel the ray hits on detector 1. |
| Geometric window | A window sized from the sources of error, each term a half-width (distance either side of the hit point). Isotropic part: voxel size times (L1/L0 − 1) plus detector-0 pixel size times L1/L0 plus the detector registration error, where L0 and L1 are the distances of detectors 0 and 1 from the sample. Plus the arc a voxel sweeps during one 2° rotation step, r·dω·(L1/L0 − 1) with r the voxel's distance from the rotation axis (larger for voxels far from the axis), applied along the arc direction. |

Decisions on record (owner, 2026-10-11 and 2026-10-12):

| Decision | Detail |
|---|---|
| Where the ray rule is used | Only in the coarse search, levels 0, 1 and 2. Level 3, the final fit (FindOptimal), the variance-minimizing step and the final reported cost all use the strict C++ rule. The two costs differ by design. |
| Window sizes | Each geometric term is a half-width. Both windows are left as set, untuned. Search the code for "WINDOW DECISION". |
| Case f | When the orientation predicts its spot misses detector 1 but the observed ray lands on dark detector-1 pixels, the spot is always counted (forgiven), as in C++. Rejecting it was tested and dropped: at the true orientation it threw away about 3% of real peaks (96.6% kept instead of 99.4%) with no benefit seen. May revisit. |
| Default | The ray rule is opt-in (config key `RayOverlapRule`), off by default in the code. Reference ("golden") results are unchanged. The owner accepted the ray rule as the future default (2026-10-12, section 9); changing the code default is a follow-up, not part of this task. Setting `RayOverlapRule` without `CppOverlapRule` gives a warning (not studied). |

## 2. Test data

Two different image sets are used:

| Image set | Used in | Clean images | Realistic images |
|---|---|---|---|
| Full-sample images | Sections 3, 4 and 6 (peak level, coarse search, full-sample cost) | The whole 500-grain sample simulated together, so neighbouring grains light pixels physically; no noise. | The same images plus detector noise (`scripts/phase_d/noise.py`): each spot on each detector is dropped entirely with probability 0.1, edge pixels flip with probability 0.05, plus stray hot pixels (0.05) and blobs (0.1). |
| Per-voxel images (E0 set) | Section 5 (end to end) | Images of the voxel alone. | Simulated at Q-max 16 with overlap from neighbouring and twin grains plus detector noise; reconstructed at Q-max 8. |

| Term | Meaning |
|---|---|
| E0 set | 200 test voxels, each with its own images (one random seed, seed 0), run with `min_sin_eta` 0 so that the C++ comparison is like for like. |
| Full-sample voxels | Voxels of the 500-grain full sample: 30 in sections 3 and 4, 6 in section 6. |
| Wrong voxel | Final orientation more than 1° from the truth (after symmetry). |

## 3. Peak level: what each rule keeps (30 voxels)

Share of detector-0 overlaps that each rule counts as a valid match, at three kinds of orientation (1.000 = keeps
everything). The two tables below are the screen numbers (initial screen: lit pixels within a pixel radius of 3 of the
predicted spot, spot placed at the voxel centroid); the local numbers (radius 0) differ, for example 0.590 against
0.656 for the C++ rule at 1° on clean images. Ideally a rule keeps nearly all real peaks at the true orientation, keeps most peaks 1° away (so the
coarse search still sees a signal near the answer), and rejects most chance overlaps at random orientations.

Clean images:

| Overlap rule | At the true orientation | 1° from the truth | 5° from the truth | Random orientation |
|---|---|---|---|---|
| Original lenient Python rule | 1.000 | 1.000 | 1.000 | 1.000 |
| Strict C++ rule | 0.994 | 0.656 | 0.279 | 0.209 |
| Ray rule, 2-pixel window | 0.994 | 0.907 | 0.371 | 0.258 |
| Ray rule, geometric window | 0.994 | 0.927 | 0.501 | 0.435 |

Realistic images:

| Overlap rule | At the true orientation | 1° from the truth | 5° from the truth | Random orientation |
|---|---|---|---|---|
| Original lenient Python rule | 0.998 | 1.005 | 1.004 | 1.000 |
| Strict C++ rule | 0.909 | 0.624 | 0.281 | 0.206 |
| Ray rule, 2-pixel window | 0.909 | 0.831 | 0.362 | 0.253 |
| Ray rule, geometric window | 0.918 | 0.856 | 0.478 | 0.409 |

How the share is computed: the denominator is the number of detector-0 overlaps found by the detector-0 test that
all ray-rule and C++ variants share, which places the predicted spot at the voxel centroid. Since the centroid task
(2026-10-13) the lenient rule's screen places the spot at the centroid too, so on the clean images it finds exactly
the same detector-0 overlaps and its share is exactly 1.000 (before: 0.9995 at 1° from the truth, 1.006 at 5° and
0.9994 at random orientations, because its screen then used the first vertex). The lenient rule never rejects a
detector-0 overlap. The realistic-image table above is the earlier vertex-0 run, so its lenient row is still slightly
off 1.000 (0.998, 1.005, 1.004, 1.000) for that reason. (With no search radius, all rules use the same test and the
lenient rule is exactly 1.000.)

Reading: the ray rule keeps far more peaks near the answer than the C++ rule (about 0.83 to 0.93 against 0.62 to
0.66 at 1°), while still rejecting most chance overlaps at random orientations. The 2-pixel window rejects chance
overlaps almost as well as the C++ rule; the geometric window is looser (keeps 0.41 to 0.44 at random
orientations).

### 3.1 Real peaks rejected and chance overlaps accepted (counts)

The same data as counts, read as two error types:

- **Real peaks rejected:** at the true orientation the detector-0 overlaps are real peaks, so every one a rule
  rejects is a rejected real peak.
- **Chance overlaps accepted:** at a random orientation the detector-0 overlaps are chance hits on light from other
  grains, so every one a rule accepts is an accepted chance overlap.

Counts are pooled over the 30 voxels (1 true orientation and 200 random orientations per voxel). The denominator is
the number of detector-0 overlaps found by the centroid test of section 3 that the C++ and ray variants share. The
lenient rule's count now equals it on clean images (it was 8,937 of 8,942 with the vertex-0 screen; the realistic images below are still that earlier run).

Initial screen (lit pixels within 3 pixels of the predicted spot), clean images:

| Overlap rule | Real peaks rejected | Chance overlaps accepted |
|---|---|---|
| Original lenient Python rule | 0 of 524 (0%) | 8,942 of 8,942 (100%) |
| Strict C++ rule | 3 of 524 (0.6%) | 1,870 of 8,942 (21%) |
| Ray rule, 2-pixel window | 3 of 524 (0.6%) | 2,305 of 8,942 (26%) |
| Ray rule, geometric window | 3 of 524 (0.6%) | 3,889 of 8,942 (43%) |
| Ray rule, geometric window, rejecting in case f (dropped) | 18 of 524 (3.4%) | 3,738 of 8,942 (42%) |

Initial screen, realistic images:

| Overlap rule | Real peaks rejected | Chance overlaps accepted |
|---|---|---|
| Original lenient Python rule | 1 of 475 (0.2%) | about all of 8,590 |
| Strict C++ rule | 43 of 475 (9.1%) | 1,772 of 8,590 (21%) |
| Ray rule, 2-pixel window | 43 of 475 (9.1%) | 2,174 of 8,590 (25%) |
| Ray rule, geometric window | 39 of 475 (8.2%) | 3,512 of 8,590 (41%) |

Cost function with no search radius (the setting that scores the final fit; only the predicted spot's own pixels
count as an overlap):

| Overlap rule | Clean: real peaks rejected | Clean: chance overlaps accepted | Realistic: real peaks rejected | Realistic: chance overlaps accepted |
|---|---|---|---|---|
| Original lenient Python rule | 0 of 521 (0%) | 2,894 of 2,894 (100%) | 0 of 471 (0%) | 3,070 of 3,070 (100%) |
| Strict C++ rule | 0 of 521 (0%) | 414 of 2,894 (14%) | 41 of 471 (8.7%) | 479 of 3,070 (16%) |
| Ray rule, 2-pixel window | 0 of 521 (0%) | 583 of 2,894 (20%) | 40 of 471 (8.5%) | 636 of 3,070 (21%) |
| Ray rule, geometric window | 0 of 521 (0%) | 1,183 of 2,894 (41%) | 36 of 471 (7.6%) | 1,195 of 3,070 (39%) |

Reading:

- The lenient rule's "1 of 475" and "0 of 524" are not rejections: the lenient rule never rejects. The 1 on
  realistic images is the counting difference explained above section 3.1 (first vertex against centroid; the clean
  row no longer has it).
- On real peaks the strict C++ rule and both ray windows behave almost the same. On realistic images all three reject
  about 8 to 9% of real peaks. This is consistent with the noise model: each spot on each detector is dropped with
  probability 0.1, so about 10% of real peaks lit on detector 0 would have no spot on detector 1, and any rule that
  asks detector 1 to confirm would reject them. Not tested: why the rate is a little under 10% (merging of spots
  with a neighbour's spot is one possible reason, not checked). The noise model itself is ruled unphysical (open
  question 8 of the parent findings), so this rejection rate is an artefact of the test data, not a property of real
  images. Q8 (2026-10-14): on images without the random drop the same three rules reject 3 of 524 (initial screen) and 0 of 521
  real peaks, the same as on clean images, and 0 of 12,439 on the E0 images (against 974 of 11,174 for the C++ rule with the
  drop at `min_sin_eta` 0.3), which fits the drop being the cause of the 8 to 9%.
- On clean images the C++ rule and both ray windows reject 3 real peaks of 524: the peak is lit on detector 0
  but both the predicted spot and the observed ray land on dark pixels of detector 1. The cause is not yet checked.
- The windows differ in how many chance overlaps they accept: the 2-pixel window a few more than the C++ rule, the
  geometric window about twice as many.
- So the ray rule's gain is not in these two error types. It is in the orientations near the answer (1° from the
  truth, table above), where the C++ rule keeps only 0.62 to 0.66 of the peaks and the ray rule 0.83 to 0.93.

Limit: these labels come from the orientation, not from each peak. On realistic images some overlaps at the true
orientation may be a neighbour's light sitting where the voxel's own peak is missing, and these cannot be separated.
True per-peak labels would need the simulator to record which voxel lit each pixel (not done).

## 4. Coarse search: is a candidate near the truth kept? (30 voxels, levels 0 to 2)

The coarse search keeps the best-scoring quarter of its candidates at each level. The question is whether a
candidate near the true orientation is among those kept.

### 4.1 Candidates at level 0 and how highly the truth ranks

"Rank of closest candidate" is the median, over 30 voxels, of the rank of the candidate closest to the truth before
the quick Monte Carlo refinement moves it (rank 0 = best score). Lower is better.

| Overlap rule | Level-0 candidates (clean) | Kept (clean) | Rank of closest candidate (clean) | Level-0 candidates (realistic) | Kept (realistic) | Rank of closest candidate (realistic) |
|---|---|---|---|---|---|---|
| Original lenient Python rule | 10,252 | 2,563 | 215.5 | 10,236 | 2,558 | 453 |
| Strict C++ rule | 7,940 | 1,985 | 108 | 7,614 | 1,903 | 382 |
| Ray rule, 2-pixel window | 8,422 | 2,105 | 22 | 8,113 | 2,028 | 102 |
| Ray rule, geometric window | 9,359 | 2,340 | 72 | 9,140 | 2,285 | 186 |

Paired by voxel (level 0, 30 voxels; ties = equal rank), the closest candidate ranks higher (smaller rank) under the
ray rule than under the strict C++ rule in this many voxels:

| Images | Rule | Higher than C++ rule | Lower than C++ rule | Tied |
|---|---|---|---|---|
| Clean | Ray rule, 2-pixel window | 17 | 10 | 3 |
| Clean | Ray rule, geometric window | 12 | 15 | 3 |
| Realistic | Ray rule, 2-pixel window | 20 | 7 | 3 |
| Realistic | Ray rule, geometric window | 20 | 7 | 3 |

(from `review_extras.json`, `scripts/phase_d/ray_overlap_review_extras.py`.) Reading: the median ranks favour the ray
rule, but the paired counts are mixed on clean images (the 2-pixel window ranks higher in 17 of 30 voxels, the
geometric window in 12 of 30 and lower in 15), and clearer on realistic images (20 against 7 for both windows). Because
about 2,000 candidates are kept, every rule keeps the closest candidate anyway, so this does not show up as survival
below.

### 4.2 Survival of a candidate within 3° of the truth (voxels of 30)

| Overlap rule | Level 0 clean | Level 1 clean | Level 2 clean | Level 0 realistic | Level 1 realistic | Level 2 realistic |
|---|---|---|---|---|---|---|
| Original lenient Python rule | 30 | 30 | 30 | 30 | 30 | 30 |
| Strict C++ rule | 30 | 30 | 30 | 30 | 29 | 29 |
| Ray rule, 2-pixel window | 30 | 30 | 30 | 30 | 30 | 30 |
| Ray rule, geometric window | 30 | 30 | 30 | 29 | 29 | 29 |

The one realistic voxel lost at levels 1 and 2 by the geometric window is the same voxel the strict C++ rule loses.
The geometric window also loses one realistic voxel at level 0 (29 of 30), which the 2-pixel window does not.

### 4.3 Survival within 2° and 1° (voxels of 30)

| Overlap rule | Within 2°, level 0 clean | Within 2°, level 0 realistic | Within 1°, level 2 clean | Within 1°, level 2 realistic |
|---|---|---|---|---|
| Original lenient Python rule | 27 | 27 | 30 | 30 |
| Strict C++ rule | 27 | 24 | 30 | 29 |
| Ray rule, 2-pixel window | 27 | 26 | 30 | 30 |
| Ray rule, geometric window | 28 | 28 | 29 | 29 |

The 1° check at level 0 is not shown: only 5 of 30 voxels have any level-0 candidate within 1°, so it measures the
grid spacing, not the rule.

Conclusion for the coarse search: with 30 voxels the rules cannot be separated by survival (each interval spans
roughly 0.83 to 1.00). The rank measure (4.1) is the only clear difference, and it favours the ray rule.

## 5. End to end: wrong voxels after the full search (E0, 200 voxels)

### 5.1 Wrong counts

Intervals are 95% Wilson intervals on the wrong fraction. "Evaluations" is the mean number of cost evaluations per
voxel. "Time" is the mean wall time per voxel with 10 workers sharing the machine, with stored runs from different
sessions: it is contended and not a timing claim (the run metadata says so). Compare cost on the evaluations column
(the ray windows add 3 to 4% over the C++ rule on clean images) and on the single-worker coarse-stage timings in
section 7. "Median error of right voxels" is a median over each rule's own right voxels (159, 144, 160 and 159 clean;
182, 144, 168 and 171 realistic, in table order), so the subsets differ.

| Images | Overlap rule | Wrong of 200 | Interval (wrong fraction) | Evaluations per voxel | Time per voxel | Median error of right voxels |
|---|---|---|---|---|---|---|
| Clean | Original lenient Python rule | 41 | 0.15 to 0.27 | 48,700 | 26.0 s | 0.0130° |
| Clean | Strict C++ rule | 56 | 0.22 to 0.35 | 47,230 | 24.1 s | 0.0223° |
| Clean | Ray rule, 2-pixel window | 40 | 0.15 to 0.26 | 48,641 | 26.6 s | 0.0097° |
| Clean | Ray rule, geometric window | 41 | 0.16 to 0.27 | 48,649 | 27.2 s | 0.0114° |
| Realistic | Original lenient Python rule | 18 | 0.06 to 0.14 | 51,163 | 26.1 s | 0.0107° |
| Realistic | Strict C++ rule | 56 | 0.22 to 0.35 | 48,769 | 25.2 s | 0.0132° |
| Realistic | Ray rule, 2-pixel window | 32 | 0.12 to 0.22 | 50,915 | 28.3 s | 0.0087° |
| Realistic | Ray rule, geometric window | 29 | 0.10 to 0.20 | 50,967 | 29.0 s | 0.0100° |

### 5.2 Voxel by voxel

Each row compares a rule with a reference rule on the same 200 voxels. "Fixed" = wrong under the reference, right
under the rule; "broken" = right under the reference, wrong under the rule; "both wrong" = wrong under both.

| Images | Rule | Compared with | Fixed | Broken | Both wrong |
|---|---|---|---|---|---|
| Clean | Ray rule, 2-pixel window | Strict C++ rule | 32 | 16 | 24 |
| Clean | Ray rule, geometric window | Strict C++ rule | 31 | 16 | 25 |
| Realistic | Ray rule, 2-pixel window | Strict C++ rule | 38 | 14 | 18 |
| Realistic | Ray rule, geometric window | Strict C++ rule | 37 | 10 | 19 |
| Clean | Ray rule, 2-pixel window | Original lenient Python rule | 10 | 9 | 31 |
| Clean | Ray rule, geometric window | Original lenient Python rule | 11 | 11 | 30 |
| Realistic | Ray rule, 2-pixel window | Original lenient Python rule | 7 | 21 | 11 |
| Realistic | Ray rule, geometric window | Original lenient Python rule | 9 | 20 | 9 |
| Clean | Strict C++ rule | Original lenient Python rule | 15 | 30 | 26 |
| Realistic | Strict C++ rule | Original lenient Python rule | 6 | 44 | 12 |

Reading:

- Against the strict C++ rule, the ray rule fixes two to nearly four voxels for every one it breaks, on both kinds
  of image.
- Against the original lenient rule, the ray rule ties on clean images (10 fixed, 9 broken for the 2-pixel window; 11
  and 11 for the geometric window) and is worse on realistic images (7 to 9 fixed, about 20 broken; the lenient
  realistic run is the earlier vertex-0 one).
- Over the voxels right under all four rules (121 clean, 125 realistic), the median error is 0.0121° (lenient),
  0.0232° (C++), 0.0097° (2-pixel window) and 0.0098° (geometric window) on clean images, and 0.0096°, 0.0133°,
  0.0087° and 0.0083° on realistic images. So on the same voxels the ray rule's right answers are more precise than
  the C++ rule's on clean images, and on clean images also a little more precise than the lenient rule's (0.0097° and
  0.0098° against 0.0121°; with the earlier vertex-0 start point the lenient rule was the most precise, 0.0082° on 113
  voxels, which changes the ordering - a difference of about 0.004°, one run per voxel, so not firm). On realistic
  images they are about equal to or better.

## 6. Cost on full-sample voxels (6 voxels)

All four rules got all six voxels right (error under 1°), so this measures cost only. One process per run, clean
images, seed 1. The lenient-rule column was re-run on 2026-10-13 with the centroid screen (tag `base_centroid`, six processes at
once; `benchmarks/phase_d_seed_diag/centroid/cache/run_all.sh`); the 2026-10-12 run with the vertex-0 screen (tag
`lenientnow`, after the C++-faithful Monte Carlo port and the variance-stage fix; `scripts/phase_d/cache/ray_lenient_rerun.sh`) is in section 10. The earlier
lenient runs (`base_full_*`, 2026-10-08) predate both changes, ran about 1.38 million variance-stage evaluations on
voxel 15901 alone, and are no longer used. The C++ and ray columns are the runs of 2026-10-10 and 2026-10-12.

| Voxel | Lenient rule: error / evaluations / time | C++ rule: error / evaluations / time | Ray, 2-pixel window: error / evaluations / time | Ray, geometric window: error / evaluations / time |
|---|---|---|---|---|
| 15901 | 0.038° / 327,619 / 179 s | 0.029° / 237,739 / 125 s | 0.014° / 250,495 / 136 s | 0.010° / 279,319 / 158 s |
| 17034 | 0.022° / 340,422 / 185 s | 0.064° / 246,127 / 125 s | 0.025° / 256,771 / 136 s | 0.037° / 284,960 / 157 s |
| 18558 | 0.005° / 329,915 / 179 s | 0.003° / 209,728 / 106 s | 0.005° / 223,744 / 119 s | 0.004° / 268,489 / 146 s |
| 19889 | 0.015° / 329,724 / 179 s | 0.021° / 226,453 / 115 s | 0.033° / 243,577 / 135 s | 0.018° / 276,634 / 150 s |
| 22867 | 0.012° / 333,046 / 182 s | 0.053° / 244,186 / 125 s | 0.061° / 257,829 / 141 s | 0.037° / 285,600 / 156 s |
| 5242 | 0.060° / 331,552 / 179 s | 0.028° / 204,646 / 102 s | 0.023° / 218,577 / 117 s | 0.059° / 265,129 / 143 s |
| **Mean** | 332,000 evaluations, 180 s | 228,000 evaluations, 116 s | 242,000 evaluations, 131 s | 277,000 evaluations, 152 s |

Timing labels ("single process"): the six lenient runs of 2026-10-13 were started together (six processes at once,
load 9.7 at the start), so none was quiet and their wall times are contended; the preflight was quiet for 5 of 6
C++ runs (17034 not quiet), 2 of 6 for the 2-pixel window (15901, 17034, 18558 and 22867 not quiet; load 1.54 to
2.14) and 5 of 6 for the geometric window (18558 not quiet). The C++ and ray runs come from two sessions (10-10 and
10-12). With the vertex-0 screen the lenient wall times were 172 to 179 s whether the start was quiet or not. The
other columns are not like for like, so rely on the evaluations.

Relative to the strict C++ rule, the 2-pixel window costs 6% more evaluations and the geometric window 21% more; the
re-run lenient rule costs 46% more (1.46 times; wall 1.55 times, but its wall times are from six processes at once, so
rely on the evaluations). Against the ray rule, the lenient rule
costs 1.37 times (2-pixel window) and 1.20 times (geometric window) as many evaluations. The lenient rule's extra
evaluations are in the coarse search (about 181,000 quick-search evaluations per voxel against 138,000 for the C++ rule
on voxel 15901). This is consistent with the lenient rule passing more candidates to the next level (30-voxel test,
section 4.1: about 10,300 against 7,900 level-0 candidates); that mechanism was not tested on these six voxels.
The earlier statement that the lenient rule costs about four to four-and-a-half times as much was a measurement of the
old variance-stage budget, not of the overlap rule, and is withdrawn. On the per-voxel E0 images (section 5.1) all
rules cost about the same.

## 7. Speed work

The coarse search scores about 233,000 orientations per voxel, one at a time, at about 0.5 ms each. Almost all of that
is Python and PyTorch overhead; the overlap counting in C is about 5% of it. The ray check was moved from a Python
loop into C, with identical results on 5,784 checked evaluations. Timing label: one process, clean images, voxels
22867, 15901 and 17034, levels 0 to 2, 3 repeats each. The preflight (`preflight_speed.json`) showed a 1-minute load
of 1.52, just above the 1.5 threshold used for the other preflights, and 8 torch threads (the full-sample runs in
section 6 used 1). So this is "single process, load 1.52, 8 threads", not a strictly quiet single-worker timing.
Relative to the C++ rule the C version of the 2-pixel window costs about 12% more and the geometric window about 27%
more:

| Overlap rule | Coarse search time per voxel, ray check in Python | Coarse search time per voxel, ray check in C |
|---|---|---|
| Original lenient Python rule | 162 to 168 s | unchanged |
| Strict C++ rule | 116 to 118 s | unchanged |
| Ray rule, 2-pixel window | 289 to 294 s | 131 to 132 s |
| Ray rule, geometric window | 335 to 340 s | 148 to 149 s |

Scoring many orientations per call (batching) is the remaining large speed-up and is recorded as a separate task:
open question 5 in `findings_finisher_mc_study_2026-10.md`.

## 8. Limits of these results

| Limit | Effect |
|---|---|
| Noise model of the realistic images | The realistic images use a noise model the owner has ruled unphysical (every spot on each detector is dropped at random with probability 0.1; open question 8 of the parent findings, deferred). The realistic-image numbers in this write-up (sections 3 to 5) are therefore reported but should not drive decisions; the clean-image numbers come first. Q8 (2026-10-14) re-ran the E0 rows on images without the drop (`MIGRATION_HISTORY.md` "Q8"); sections 3 to 5 still show the old realistic numbers. |
| One run per voxel for E0 (seed 0) | The differences between the two ray-rule windows (32 against 29 wrong on realistic images) are a handful of voxels either way. More seeds are not needed (owner, 2026-10-12): the search uses randomness to explore, not as independent samples. |
| Realistic-image rows of the lenient and C++ rules are stored earlier runs | On clean images all four rules were re-run in this build (section 10). The realistic columns were not: the lenient rule's realistic numbers still use the vertex-0 start point. |
| Full-sample voxels all right | The six full-sample voxels measure cost only, not accuracy. |
| Not separated | Whether the gain comes from the ray rule in the initial screen, the re-scoring, or the quick refinement (all three use it at levels 0 to 2). |
| Not tested | More than two detectors. |

## 9. Decisions for the owner

Clean-image results first: on the clean E0 images the 2-pixel ray window (40 wrong of 200) ties the lenient rule (41
wrong; paired 10 fixed, 9 broken, 31 both wrong) and is better than the C++ rule (56 wrong). The realistic-image gap in
the lenient rule's favour (18 against 29 to 32) rests on a noise model the owner has ruled unphysical (open question 8
of the parent findings), so it is not used for the decisions below. On the no-drop images (Q8, 2026-10-14) the gap reverses:
14 wrong for each ray window against 22 for the lenient rule (2-pixel window: 15 fixed, 7 broken, 7 both wrong), against
38 for the C++ rule, at 51,652 and 51,677 evaluations per voxel for the ray windows against 51,901 (lenient) and 49,281 (C++ rule);
whether that changes a decision below is the owner's.

| Question | Suggestion or decision |
|---|---|
| Which window goes forward: 2-pixel or geometric? | 2-pixel window: about equally accurate, cheaper (6% against 21% more evaluations than the C++ rule on the six full-sample voxels), and it ranks the closest candidate higher than the C++ rule in 17 of 30 clean voxels (lower in 10; the geometric window 12 against 15). |
| Make the ray rule the default? | Decided (owner, 2026-10-12): the ray rule as the default is fine. It ties the lenient rule on clean images and costs 6% more evaluations than the C++ rule on full-sample voxels, against 46% for the lenient rule (re-run in the current build). The code default is NOT changed in this task (still off); changing it is a follow-up. |
| `RayOverlapRule` without `CppOverlapRule`? | Decided (owner, 2026-10-12): warn only (a `UserWarning` in `SearchParameters.from_config`); the setting is not forced. |
| Cases b and e pixel counting in triangle mode (+1 / -1 pixel) | Decided (owner, 2026-10-12): keep as is, documented in `ray_overlap.py` and `_rasterize.c`; open question 9 in `findings_finisher_mc_study_2026-10.md`, fix later. |
| More seeds? | Decided (owner, 2026-10-12): no. The algorithm does not rely on independence of the random streams; it uses randomness to explore the search space, so repeat seeds are not a measure of the rule. |

## 10. Centroid start point for every rule (2026-10-13)

What changed: the pixel-radius screen asks "is a lit pixel within 3 pixels of where this reflection lands?". C++
places the reflection by projecting the voxel centroid. The Python lenient rule used to project the first corner of
the voxel's triangle instead; the C++ rule and the ray windows already used the centroid. Now every rule uses the
centroid (commit d96c413). The change touches the lenient rule only; the other three rules are a check that nothing
else moved. All numbers below are clean images, one seed, the same voxels and seeds as before. "Before" = lenient rule
with the first-corner screen (stored runs); "after" = centroid screen. Run record:
`benchmarks/phase_d_seed_diag/centroid/RUN_LOG.md`.

Arm names in the tables: `old` = original lenient Python rule, `cpp` = strict C++ rule, `f0_R2` = ray rule with the
2-pixel window, `f0_A` = ray rule with the geometric window. In the paired tables a row "X vs old (after)" counts
fixed = wrong under old and right under X, broken = the reverse.

### 10.1 Peak level (30 voxels)

Columns: `kept_*` = peaks kept / detector-0 overlaps; `share_*` = their ratio; `pass_*` = evaluations with at least one
kept peak / evaluations. Lenient rule: after the change it keeps exactly all detector-0 overlaps in every row (524 of
524 at the truth, 2,067 of 2,067 at 1°, 8,942 of 8,942 at random orientations); before it showed 8,937 of 8,942 at
random orientations and 2,066 of 2,067 at 1°. The C++ and ray rows are identical before and after. The old statement
that the lenient rule "rejected 1" came from the first-corner screen finding a slightly different set of overlaps than
the centroid test used as the denominator, not from a rejection.

<!-- table:centroid_peak_clean -->
| cost | cls | arm | kept_before | kept_after | share_before | share_after | pass_before | pass_after |
|---|---|---|---|---|---|---|---|---|
| screen | truth | old | 524/524 | 524/524 | 1.0000 | 1.0000 | 30/30 | 30/30 |
| screen | truth | cpp | 521/524 | 521/524 | 0.9943 | 0.9943 | 30/30 | 30/30 |
| screen | truth | f0_R2 | 521/524 | 521/524 | 0.9943 | 0.9943 | 30/30 | 30/30 |
| screen | truth | f0_A | 521/524 | 521/524 | 0.9943 | 0.9943 | 30/30 | 30/30 |
| screen | near1.0 | old | 2066/2067 | 2067/2067 | 0.9995 | 1.0000 | 299/300 | 299/300 |
| screen | near1.0 | cpp | 1356/2067 | 1356/2067 | 0.6560 | 0.6560 | 292/300 | 292/300 |
| screen | near1.0 | f0_R2 | 1875/2067 | 1875/2067 | 0.9071 | 0.9071 | 292/300 | 292/300 |
| screen | near1.0 | f0_A | 1916/2067 | 1916/2067 | 0.9269 | 0.9269 | 294/300 | 294/300 |
| screen | random | old | 8937/8942 | 8942/8942 | 0.9994 | 1.0000 | 4686/6000 | 4696/6000 |
| screen | random | cpp | 1870/8942 | 1870/8942 | 0.2091 | 0.2091 | 1598/6000 | 1598/6000 |
| screen | random | f0_R2 | 2305/8942 | 2305/8942 | 0.2578 | 0.2578 | 1874/6000 | 1874/6000 |
| screen | random | f0_A | 3889/8942 | 3889/8942 | 0.4349 | 0.4349 | 2843/6000 | 2843/6000 |
| local | truth | old | 521/521 | 521/521 | 1.0000 | 1.0000 | 30/30 | 30/30 |
| local | truth | cpp | 521/521 | 521/521 | 1.0000 | 1.0000 | 30/30 | 30/30 |
| local | truth | f0_R2 | 521/521 | 521/521 | 1.0000 | 1.0000 | 30/30 | 30/30 |
| local | truth | f0_A | 521/521 | 521/521 | 1.0000 | 1.0000 | 30/30 | 30/30 |
| local | near1.0 | old | 742/742 | 742/742 | 1.0000 | 1.0000 | 282/300 | 282/300 |
| local | near1.0 | cpp | 438/742 | 438/742 | 0.5903 | 0.5903 | 244/300 | 244/300 |
| local | near1.0 | f0_R2 | 651/742 | 651/742 | 0.8774 | 0.8774 | 271/300 | 271/300 |
| local | near1.0 | f0_A | 689/742 | 689/742 | 0.9286 | 0.9286 | 275/300 | 275/300 |
| local | random | old | 2894/2894 | 2894/2894 | 1.0000 | 1.0000 | 2307/6000 | 2307/6000 |
| local | random | cpp | 414/2894 | 414/2894 | 0.1431 | 0.1431 | 399/6000 | 399/6000 |
| local | random | f0_R2 | 583/2894 | 583/2894 | 0.2015 | 0.2015 | 555/6000 | 555/6000 |
| local | random | f0_A | 1183/2894 | 1183/2894 | 0.4088 | 0.4088 | 1070/6000 | 1070/6000 |
<!-- /table:centroid_peak_clean -->

### 10.2 Coarse search (30 voxels, levels 0 to 2)

Columns: `surv3` / `surv2` / `surv1` = voxels (of 30, with Wilson interval for `surv3`) that keep a candidate within 3°
/ 2° / 1° of the truth after the quick refinement; `in_disc3` = voxels whose level list has a candidate within 3° before
the refinement; `keep_median` = candidates kept at the ranking cut; `rank_disc_median` = rank of the closest
candidate before the refinement (smaller is better); `cand_median` and `evals_median` = candidates and cumulative
evaluations per voxel (medians). The C++ and both ray rows are identical before and after (identical medians, including the
evaluation counts). The lenient rule: level 0 survival at 3° goes from 29 to 30 of 30 (one voxel fixed, none broken;
levels 1 and 2 unchanged at 30), survival at 1° at level 1 from 25 to 29, level-0 candidates from 10,256 to 10,252
and median evaluations at level 2 from 318,258 to 318,468 (+0.07%). These are small shifts; with 30 voxels the
intervals still span about 0.89 to 1.00.

<!-- table:centroid_coarse_clean -->
| arm | level | when | surv3 | surv2 | surv1 | keep_median | rank_disc_median | in_disc3 | cand_median | evals_median |
|---|---|---|---|---|---|---|---|---|---|---|
| old | 0 | before | 29/30 [0.83, 0.99] | 26 | 18 | 2.56e+03 | 214 | 30 | 10256 | 202008 |
| old | 0 | after | 30/30 [0.89, 1.00] | 27 | 20 | 2.56e+03 | 216 | 30 | 10252 | 202028 |
| old | 1 | before | 30/30 [0.89, 1.00] | 29 | 25 | 881 | 1 | 30 | 3524 | 289357 |
| old | 1 | after | 30/30 [0.89, 1.00] | 30 | 29 | 881 | 0 | 30 | 3526 | 289421 |
| old | 2 | before | 30/30 [0.89, 1.00] | 30 | 30 | 279 | 0 | 30 | 1116 | 318258 |
| old | 2 | after | 30/30 [0.89, 1.00] | 30 | 30 | 278 | 0 | 30 | 1114 | 318468 |
| cpp | 0 | before | 30/30 [0.89, 1.00] | 27 | 20 | 1.98e+03 | 108 | 30 | 7940 | 151172 |
| cpp | 0 | after | 30/30 [0.89, 1.00] | 27 | 20 | 1.98e+03 | 108 | 30 | 7940 | 151172 |
| cpp | 1 | before | 30/30 [0.89, 1.00] | 30 | 29 | 604 | 0 | 30 | 2416 | 207912 |
| cpp | 1 | after | 30/30 [0.89, 1.00] | 30 | 29 | 604 | 0 | 30 | 2416 | 207912 |
| cpp | 2 | before | 30/30 [0.89, 1.00] | 30 | 30 | 164 | 0 | 30 | 658 | 225284 |
| cpp | 2 | after | 30/30 [0.89, 1.00] | 30 | 30 | 164 | 0 | 30 | 658 | 225284 |
| f0_R2 | 0 | before | 30/30 [0.89, 1.00] | 27 | 19 | 2.1e+03 | 21.5 | 30 | 8422 | 159068 |
| f0_R2 | 0 | after | 30/30 [0.89, 1.00] | 27 | 19 | 2.1e+03 | 21.5 | 30 | 8422 | 159068 |
| f0_R2 | 1 | before | 30/30 [0.89, 1.00] | 30 | 29 | 642 | 0.5 | 30 | 2568 | 220919 |
| f0_R2 | 1 | after | 30/30 [0.89, 1.00] | 30 | 29 | 642 | 0.5 | 30 | 2568 | 220919 |
| f0_R2 | 2 | before | 30/30 [0.89, 1.00] | 30 | 30 | 177 | 0 | 30 | 710 | 239760 |
| f0_R2 | 2 | after | 30/30 [0.89, 1.00] | 30 | 30 | 177 | 0 | 30 | 710 | 239760 |
| f0_A | 0 | before | 30/30 [0.89, 1.00] | 28 | 21 | 2.34e+03 | 71.5 | 30 | 9359 | 177942 |
| f0_A | 0 | after | 30/30 [0.89, 1.00] | 28 | 21 | 2.34e+03 | 71.5 | 30 | 9359 | 177942 |
| f0_A | 1 | before | 30/30 [0.89, 1.00] | 30 | 28 | 698 | 0 | 30 | 2792 | 249730 |
| f0_A | 1 | after | 30/30 [0.89, 1.00] | 30 | 28 | 698 | 0 | 30 | 2792 | 249730 |
| f0_A | 2 | before | 30/30 [0.89, 1.00] | 30 | 29 | 191 | 0 | 30 | 764 | 270908 |
| f0_A | 2 | after | 30/30 [0.89, 1.00] | 30 | 29 | 191 | 0 | 30 | 764 | 270908 |
<!-- /table:centroid_coarse_clean -->

Paired by voxel, after against before (fixed = fails before and survives after; broken = the reverse; both fail =
fails both times):

<!-- table:centroid_coarse_paired -->
| arm | level | fixed | broken | both_fail |
|---|---|---|---|---|
| old | 0 | 1 | 0 | 0 |
| old | 1 | 0 | 0 | 0 |
| old | 2 | 0 | 0 | 0 |
| cpp | 0 | 0 | 0 | 0 |
| cpp | 1 | 0 | 0 | 0 |
| cpp | 2 | 0 | 0 | 0 |
| f0_R2 | 0 | 0 | 0 | 0 |
| f0_R2 | 1 | 0 | 0 | 0 |
| f0_R2 | 2 | 0 | 0 | 0 |
| f0_A | 0 | 0 | 0 | 0 |
| f0_A | 1 | 0 | 0 | 0 |
| f0_A | 2 | 0 | 0 | 0 |
<!-- /table:centroid_coarse_paired -->

### 10.3 End to end (E0, 200 voxels, min_sin_eta 0, clean)

Wrong = error above the usual threshold. The C++, 2-pixel and geometric rows are identical before and after: the same
wrong voxels, the same mean evaluations (so the stored runs reproduce in this build). The lenient rule gives 41 wrong
(interval 0.15 to 0.27) against 39 before. Voxel by voxel the start-point change alone fixed 14, broke 16 and left 25
wrong, so 30 of 200 voxels (15%) change outcome in either direction and the net is +2 wrong. With one run per voxel,
this does not show that either start point is more accurate; it shows that a small change in the candidate list
(a few of 10,256 level-0 candidates) moves single runs about as much as the rule differences of section 5 (not
tested: that the flips come from the stochastic search path rather than from a systematic effect). Evaluation cost
per voxel: 48,700 against 48,501, 0.4% more.

<!-- table:centroid_e0_clean -->
| arm | when | wrong | mean_evals | median_err_right |
|---|---|---|---|---|
| old | before | 39/200 [0.146, 0.255] | 48501 | 0.0086 |
| old | after | 41/200 [0.155, 0.266] | 48700 | 0.0130 |
| cpp | before | 56/200 [0.222, 0.346] | 47230 | 0.0223 |
| cpp | after | 56/200 [0.222, 0.346] | 47230 | 0.0223 |
| f0_R2 | before | 40/200 [0.150, 0.261] | 48641 | 0.0097 |
| f0_R2 | after | 40/200 [0.150, 0.261] | 48641 | 0.0097 |
| f0_A | before | 41/200 [0.155, 0.266] | 48649 | 0.0114 |
| f0_A | after | 41/200 [0.155, 0.266] | 48649 | 0.0114 |
<!-- /table:centroid_e0_clean -->

Paired, after against before (rows for the other arms are the check that nothing else moved), then the other rules
against the new lenient result (the fixed / broken split of section 5.2 with the new lenient column):

<!-- table:centroid_e0_paired -->
| arm | fixed | broken | both_wrong |
|---|---|---|---|
| old | 14 | 16 | 25 |
| cpp | 0 | 0 | 56 |
| f0_R2 | 0 | 0 | 40 |
| f0_A | 0 | 0 | 41 |
| cpp vs old (after) | 15 | 30 | 26 |
| f0_R2 vs old (after) | 10 | 9 | 31 |
| f0_A vs old (after) | 11 | 11 | 30 |
<!-- /table:centroid_e0_paired -->

### 10.4 Six full-sample voxels, lenient rule

All six are right (error below 1°) before and after. Errors move in both directions (up on 15901, 18558 and 5242, down on
the other three), as one run per voxel would. Mean evaluations 332,320 before and 332,046 after (-0.1%). The "after"
wall times are from six processes at once (contended); "before" is the earlier run (`lenientnow`, mixed labels), so do
not read the 3% difference in seconds as an effect.

<!-- table:centroid_six_old -->
| voxel | err_before | err_after | evals_before | evals_after | wall_before | wall_after |
|---|---|---|---|---|---|---|
| 15901 | 0.0170 | 0.0378 | 328519 | 327619 | 172 | 179 |
| 17034 | 0.0284 | 0.0217 | 339561 | 340422 | 176 | 185 |
| 18558 | 0.0048 | 0.0050 | 329367 | 329915 | 175 | 179 |
| 19889 | 0.0210 | 0.0145 | 330209 | 329724 | 175 | 179 |
| 22867 | 0.0263 | 0.0119 | 335687 | 333046 | 176 | 182 |
| 5242 | 0.0493 | 0.0604 | 330577 | 331552 | 179 | 179 |
<!-- /table:centroid_six_old -->
