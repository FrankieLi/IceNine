# Ray-consistency overlap rule: results (2026-10-11 to 2026-10-12)

Task branch `feature/finisher-mc-study-ray-overlap` (parent `feature/finisher-mc-study`). Commits 8d35ff0, 362799e
(rule, steps 1-2), f6ad076 (case f rejection dropped, ray check in C), 5f0cd9c (window decision), 387f3b6 (30-voxel
coarse test and end-to-end test). Raw tables: `benchmarks/phase_d_seed_diag/ray_overlap/` (`tables.md`,
`tables30.md`, `e2e_tables.md`); run details in its `RUN_LOG.md`.

## Summary

The ray rule recovers most of the accuracy that the strict C++ overlap rule loses, at nearly the C++ rule's speed.
It does not reach the accuracy of the original lenient Python rule on realistic images, but the original rule
needs about four times as many cost evaluations on full-sample voxels.

| Overlap rule | Wrong of 200 (clean images) | Wrong of 200 (realistic images) | Time per full-sample voxel |
|---|---|---|---|
| Original lenient Python rule | 39 | 18 | 554 s |
| Strict C++ rule | 56 | 56 | 116 s |
| Ray rule, 2-pixel window | 40 | 32 | 131 s |
| Ray rule, geometric window | 41 | 29 | 152 s |

This is a single random seed. The two ray-rule windows cannot be told apart on these numbers.

## 1. What the rule does

When a predicted diffraction spot lands on lit pixels on the near detector (detector 0), the cost function must
decide whether that counts as a matched peak. The study compares three ways to decide:

| Overlap rule | How it decides whether a lit spot on detector 0 counts |
|---|---|
| Original lenient Python rule | Counts almost every overlap. |
| Strict C++ rule | The C++ `CountQualifiedPeaks` rule: the spot counts only if the same orientation's predicted spot on the far detector (detector 1) is also lit. |
| Ray rule (owner's design) | Traces the observed ray from the voxel through the lit pixels on detector 0 out to detector 1, and counts the spot only if detector 1 is lit near where that ray lands. |

"Near where the ray lands" needs a tolerance window on detector 1. Two windows were carried forward:

| Window | Definition |
|---|---|
| 2-pixel window | A square of 2 pixels on each side of the pixel the ray hits on detector 1. |
| Geometric window | A window sized from the sources of error, each term a half-width (distance either side of the hit point): voxel size and detector-0 pixel size magnified to detector 1, the detector registration error, plus the arc a voxel sweeps during one 2° rotation step (larger for voxels far from the rotation axis), applied along the arc direction. |

Decisions on record (owner, 2026-10-12):

| Decision | Detail |
|---|---|
| Where the ray rule is used | Only in the coarse search, levels 0, 1 and 2. Level 3, the final fit (FindOptimal), the variance-minimizing step and the final reported cost all use the strict C++ rule. The two costs differ by design. |
| Window sizes | Each geometric term is a half-width. Both windows are left as set, untuned. Search the code for "WINDOW DECISION". |
| Case f | When the orientation predicts its spot misses detector 1 but the observed ray lands on dark detector-1 pixels, the spot is always counted (forgiven), as in C++. Rejecting it was tested and dropped: at the true orientation it threw away about 3% of real peaks (96.6% kept instead of 99.4%) with no benefit seen. May revisit. |
| Default | The ray rule is opt-in (config key `RayOverlapRule`), off by default. Reference ("golden") results are unchanged. |

## 2. Test data

| Term | Meaning |
|---|---|
| Clean images | Simulated detector images of the voxel alone. |
| Realistic images | Simulated at Q-max 16 with overlap from neighbouring and twin grains plus detector noise; reconstructed at Q-max 8. |
| E0 set | 200 test voxels, each with its own images (one random seed, seed 0), run with `min_sin_eta` 0 so that the C++ comparison is like for like. |
| Full-sample voxels | 6 voxels from the 500-grain full sample, where neighbouring grains light many pixels. |
| Wrong voxel | Final orientation more than 1° from the truth (after symmetry). |

## 3. Peak level: what each rule keeps (30 voxels)

Share of detector-0 overlaps that each rule counts as a valid match, at three kinds of orientation (1.000 = keeps
everything). Ideally a rule keeps nearly all real peaks at the true orientation, keeps most peaks 1° away (so the
coarse search still sees a signal near the answer), and rejects most chance overlaps at random orientations.

Clean images:

| Overlap rule | At the true orientation | 1° from the truth | 5° from the truth | Random orientation |
|---|---|---|---|---|
| Original lenient Python rule | 1.000 | 1.000 | 1.006 | 0.999 |
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

(Values slightly above 1.000 for the lenient rule come from how overlaps are counted; read them as "keeps
everything".)

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
the number of detector-0 overlaps found by the original lenient rule.

Initial screen (lit pixels within 3 pixels of the predicted spot), clean images:

| Overlap rule | Real peaks rejected | Chance overlaps accepted |
|---|---|---|
| Original lenient Python rule | 0 of 524 (0%) | 8,937 of 8,942 (about 100%) |
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

- On real peaks the strict C++ rule and both ray windows behave almost the same. On realistic images all three reject
  about 8 to 9% of real peaks, nearly the same count. The likely cause is the data (a real peak whose detector-1 spot
  is missing or faint), not the rule; these data cannot confirm it.
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
| Original lenient Python rule | 10,256 | 2,564 | 214 | 10,236 | 2,558 | 453 |
| Strict C++ rule | 7,940 | 1,985 | 108 | 7,614 | 1,903 | 382 |
| Ray rule, 2-pixel window | 8,422 | 2,105 | 22 | 8,113 | 2,028 | 102 |
| Ray rule, geometric window | 9,359 | 2,340 | 72 | 9,140 | 2,285 | 186 |

Reading: the ray rule, especially the 2-pixel window, ranks the true orientation much nearer the top. Because about
2,000 candidates are kept, every rule keeps it anyway, so this improvement does not show up as survival below.

### 4.2 Survival of a candidate within 3° of the truth (voxels of 30)

| Overlap rule | Level 0 clean | Level 1 clean | Level 2 clean | Level 0 realistic | Level 1 realistic | Level 2 realistic |
|---|---|---|---|---|---|---|
| Original lenient Python rule | 29 | 30 | 30 | 30 | 30 | 30 |
| Strict C++ rule | 30 | 30 | 30 | 30 | 29 | 29 |
| Ray rule, 2-pixel window | 30 | 30 | 30 | 30 | 30 | 30 |
| Ray rule, geometric window | 30 | 30 | 30 | 29 | 29 | 29 |

The one realistic voxel lost at levels 1 and 2 by the geometric window is the same voxel the strict C++ rule loses.

### 4.3 Survival within 2° and 1° (voxels of 30)

| Overlap rule | Within 2°, level 0 clean | Within 2°, level 0 realistic | Within 1°, level 2 clean | Within 1°, level 2 realistic |
|---|---|---|---|---|
| Original lenient Python rule | 26 | 27 | 30 | 30 |
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
voxel; "time" is the mean wall time per voxel (10 workers sharing the machine).

| Images | Overlap rule | Wrong of 200 | Interval (wrong fraction) | Evaluations per voxel | Time per voxel | Median error of right voxels |
|---|---|---|---|---|---|---|
| Clean | Original lenient Python rule | 39 | 0.15 to 0.26 | 48,501 | 24.2 s | 0.0086° |
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
| Clean | Ray rule, 2-pixel window | Original lenient Python rule | 19 | 20 | 20 |
| Clean | Ray rule, geometric window | Original lenient Python rule | 18 | 20 | 21 |
| Realistic | Ray rule, 2-pixel window | Original lenient Python rule | 7 | 21 | 11 |
| Realistic | Ray rule, geometric window | Original lenient Python rule | 9 | 20 | 9 |
| Clean | Strict C++ rule | Original lenient Python rule | 17 | 34 | 22 |
| Realistic | Strict C++ rule | Original lenient Python rule | 6 | 44 | 12 |

Reading:

- Against the strict C++ rule, the ray rule fixes two to nearly four voxels for every one it breaks, on both kinds
  of image.
- Against the original lenient rule, the ray rule ties on clean images (about 19 fixed, 20 broken) and is worse on
  realistic images (7 to 9 fixed, about 20 broken).
- The ray rule's right voxels are also more precise than the C++ rule's (median error 0.009° to 0.011° against
  0.013° to 0.022°).

## 6. Cost on full-sample voxels (6 voxels)

All four rules got all six voxels right, so this measures cost only. Single process per run.

| Voxel | Lenient rule: error / evaluations / time | C++ rule: error / evaluations / time | Ray, 2-pixel window: error / evaluations / time | Ray, geometric window: error / evaluations / time |
|---|---|---|---|---|
| 15901 | 0.042° / 1,687,675 / 897 s | 0.029° / 237,739 / 125 s | 0.014° / 250,495 / 136 s | 0.010° / 279,319 / 158 s |
| 17034 | 0.030° / 952,819 / 509 s | 0.064° / 246,127 / 125 s | 0.025° / 256,771 / 136 s | 0.037° / 284,960 / 157 s |
| 18558 | 0.006° / 463,626 / 245 s | 0.003° / 209,728 / 106 s | 0.005° / 223,744 / 119 s | 0.004° / 268,489 / 146 s |
| 19889 | 0.024° / 626,184 / 333 s | 0.021° / 226,453 / 115 s | 0.033° / 243,577 / 135 s | 0.018° / 276,634 / 150 s |
| 22867 | 0.024° / 1,339,563 / 753 s | 0.053° / 244,186 / 125 s | 0.061° / 257,829 / 141 s | 0.037° / 285,600 / 156 s |
| 5242 | 0.028° / 1,063,389 / 589 s | 0.028° / 204,646 / 102 s | 0.023° / 218,577 / 117 s | 0.059° / 265,129 / 143 s |
| **Mean** | 1.02 million evaluations, 554 s | 228,000 evaluations, 116 s | 242,000 evaluations, 131 s | 277,000 evaluations, 152 s |

Relative to the strict C++ rule, the 2-pixel window costs 6% more evaluations and the geometric window 21% more. The
original lenient rule costs about 4.5 times as much. On the per-voxel E0 images (section 5.1) all rules cost about
the same; the difference appears only on full-sample images, where stray light from neighbouring grains makes the
lenient rule accept many more candidates.

## 7. Speed work

The coarse search scores about 233,000 orientations per voxel, one at a time, at about 0.5 ms each. Almost all of that
is Python and PyTorch overhead; the overlap counting in C is about 5% of it. The ray check was moved from a Python
loop into C, with identical results on 5,784 checked evaluations:

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
| One random seed for E0 | The differences between the two ray-rule windows (32 against 29 wrong on realistic images) are within noise. |
| Lenient and C++ rows of the end-to-end test are stored earlier runs | Their code path is unchanged in this build, but they were not re-run alongside the ray rule. |
| Full-sample voxels all right | The six full-sample voxels measure cost only, not accuracy. |
| Not separated | Whether the gain comes from the ray rule in the initial screen, the re-scoring, or the quick refinement (all three use it at levels 0 to 2). |
| Not tested | More than two detectors. |

## 9. Decisions for the owner

| Question | Suggestion |
|---|---|
| Which window goes forward: 2-pixel or geometric? | 2-pixel window: about equally accurate, cheaper (6% against 21% more evaluations than the C++ rule), and it ranks the truth highest in the coarse search. |
| Make the ray rule the default? | Not yet. The original lenient rule is still more accurate on realistic images (18 against 29 to 32 wrong); its 4.5 times higher cost on full-sample voxels is the trade-off to weigh before the 500-grain run. |
| More seeds? | Run two more E0 seeds on realistic images with the lenient rule, the C++ rule and the 2-pixel ray rule (about 1 hour) before a final decision. |
