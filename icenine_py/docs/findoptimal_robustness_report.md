---
title: "FindOptimal robustness: why the multi-level search returns CSL relatives, and what fixes it"
subtitle: "Diagnosis, cheap fixes, and a feasibility test of a low-cost classifier"
date: "2026-10-06"
geometry: margin=1in
fontsize: 11pt
---

# 1. Question

Full multi-level reconstruction from scratch (`AdaptiveVoxelReconstructor.reconstruct_voxel`:
coarse levels, then FindOptimal, VarianceMinimizing, final overlap) returns a wrong orientation
(more than 1 degree from the truth) in 21 of 50 clean and 17 of 50 realistic single-voxel runs
(`MIGRATION_HISTORY.md`, "Comparison with multi-level reconstruction"). The wrong answers are
mostly coincidence-site-lattice (CSL) relatives of the truth, and in every wrong case the local cost
at the truth is lower than at the returned answer: the search fails, not the cost. Four questions:

1. Where in the search is the truth lost?
2. Which cheap fixes (no machine learning) remove the failures, and at what cost?
3. (2a) Would a low-cost classifier make the search more robust?
4. (2b) What training data would such a classifier need?

Terminology: "realistic" is the detector image with neighbours, grain mates, noise and edits
(`all` in the code); "clean" is the target voxel's own spots only.

# 2. Setup

- **Voxels.** 200 random voxels of `Examples/Example2.ManyGrains/SimInput/rand_500grains_1mm_inFZ.mic`:
  the perturbation sweep's selection (seed 0, the 30 NN-dataset voxels excluded) extended from 50 to
  200 voxels. The first 50 are the sweep's 50 voxels (asserted). The voxel list and their properties
  are in `benchmarks/findoptimal_robustness/voxels.npz` (the 30 excluded NN-dataset voxels are in
  `perturbation_sweep_raw.npz`, `excluded_indices`).
- **Images.** The per-voxel detector images of experiment A of `scripts/findoptimal_sweep.py`
  (target spots at the truth; realistic adds distractor sources over the whole detector and the
  sweep's window-level realism edits). Nothing simulates the whole sample. No voxel failed to build.
- **Search settings.** `ReconstructQ8.config` (5 deg grid radius, 4 levels, 200 MC steps, 2 restarts,
  30 candidates), Q_max 8, |sin eta| >= 0.3, as in the sweep.
- **Runs.** 3 reconstruction seeds per (voxel, variant): 1200 instrumented runs. Seed 0 uses
  experiment A's random stream and reproduces it exactly (21/50 and 17/50 wrong on the first 50
  voxels, bit-identical orientations), which also shows that the recorder hook changes nothing.
- **Error and labels.** Cubic-symmetry-reduced misorientation to the truth; "wrong" is > 1 degree.
  A candidate "in the basin" is within 1 or 3 degrees of the truth. CSL label by the Brandon
  criterion (15 deg / sqrt(Sigma)), Sigma <= 29 (21 entries, 13a/13b ... 29a/29b).
- **Intervals.** Wilson 95% intervals on rates. Runs of one voxel are not independent, so intervals
  over runs are somewhat optimistic; the per-voxel comparison in Section 4 shows how strong the
  voxel dependence is.
- **Compute.** E0 took about 57 minutes on 10 CPU workers (median 25 s per run). The whole study,
  including the re-run fixes and the classifier runs, took about 4 hours of wall time.

The recorder hook (`AdaptiveVoxelReconstructor.recorder`) records, per level, every candidate
(discrete score, post-quick-MC cost, rank, kept flag) and each FindOptimal result; it never touches
the random stream. `tests/test_findoptimal_refactor.py::test_recorder_hook_leaves_reconstruct_voxel_bit_identical`
checks it against the golden numbers. The search knobs used for the fixes (`keep_fraction`,
`keep_union_discrete`, `n_q_start_offset`, `global_pixel_radius`, `extra_candidates`, `rank_key`)
default to the existing behaviour.

# 3. Port against C++

I compared the candidate selection of `Src/DiscreteAdaptive.tmpl.cpp` (`ReconstructVoxel`,
`RunDiscreteSearch`) and `Src/DiscreteSearch.h` (`GetSpacedCandidates`) with the Python:
per level a discrete search over the FZ x local grid with the pixel-radius-3 cost and Q_max
5 + level, the local-cost re-evaluation, a quick MC (10 steps, 5 restarts), the spacing filter, sort,
keep a quarter, diameter / 1.5, FindOptimal on the top `nMaxDiscreteCandidates`. I found no difference
that could produce this failure. One small difference: C++ keeps `floor(n / 4)` candidates, which is
0 when fewer than 4 candidates exist, the Python keeps at least 1; it matters only in degenerate runs.
I did not compare the MC internals line by line. Together with the TODO note that C++ also fails on
ThreeVoxels voxel 2, this points to the algorithm, not to a port bug.

**Structural fact that matters below.** The candidate set shrinks by a factor of four per level
(median clean 172, 45, 11, 2 candidates after the quick MC at levels 0..3; realistic 248, 67, 17, 4).
FindOptimal therefore receives 2 to 4 candidates, and `max_discrete_candidates = 30` never binds.
At every pruning step a single ranking, the post-quick-MC local cost, decides what survives.

# 4. Diagnosis (E0)

Wrong rates (1200 runs; 200 voxels x 3 seeds):

| | clean | realistic |
|---|---|---|
| wrong, all runs | 204/600 = 34.0% (CI 30.3-37.9) | 143/600 = 23.8% (CI 20.6-27.4) |
| per seed | 33.5%, 35.5%, 33.0% | 26.0%, 23.0%, 22.5% |
| first 50 voxels, seed 0 | 21/50 | 17/50 |
| the 150 new voxels, seed 0 | 46/150 = 30.7% | 35/150 = 23.3% |

The cost at the truth is below the cost of the returned answer in 347 of 347 wrong runs.

**Where the truth is lost** (wrong runs; plot `where_lost.png`). A basin candidate is one within
3 degrees after the quick MC.

| Class | Meaning | clean (204) | realistic (143) |
|---|---|---|---|
| S1 | no basin candidate at level 0 | 12 (5.9%) | 5 (3.5%) |
| S2 | basin candidate present, then pruned | 179 (87.7%) | 118 (82.5%) |
|  | of which pruned at level 0 / 1 / 2 | 70 / 78 / 31 | 52 / 46 / 20 |
| S2b | kept at level L, absent at level L + 1 | 4 | 1 |
| S3 | in the hand-off to FindOptimal, but it returns a trap | 9 (4.4%) | 19 (13.3%) |

- A basin candidate exists at level 0 in 98% (clean) and 99% (realistic) of all runs; it survives
  the pruning at level 0 in 86% and 90%, and reaches the hand-off in 68% and 79% of runs.
- **S2 is the failure.** At the pruning step the best basin candidate has median rank 28 of 51
  (clean; the cut is at 12) and 38 of 105 (realistic; cut at 26). 70% and 75% of them rank within
  twice the cut, so a moderately larger keep would catch many (F2a below).
- **Mechanism (`s2_mechanism.txt`).** The pruning ranks by the local cost (pixel radius 0, all
  reflections). A basin candidate after the quick MC is still about 1-2 degrees off the truth
  (median 0.9 to 2.2 degrees), and the cost is extremely sharp: its cost is 0.94-0.99, close to "no
  overlap", although the cost at the exact truth is 0.06 (clean) or 0.27 (realistic). The candidate
  that ranks first is a wrong one with cost 0.82-0.93, usually 38 or 60 degrees from the truth (a
  Sigma7/Sigma5/Sigma3 relative or something unrelated). The basin candidate's cost is higher in 100%
  of these cases: it loses a ranking among near-1 costs, decided by a handful of chance hits. The wrong
  candidate that survives then becomes polished by the later levels. Ranking by the discrete-stage
  score (radius-3 tolerance) alone would have kept the basin candidate in 58% (clean) and 65%
  (realistic) of the S2 cases.
- **S3.** Of the 9 clean and 19 realistic S3 runs, FindOptimal refined a basin candidate to under
  1 degree in 5 and 17, but that refined result had a higher cost than the trap's: it stalled short of
  the optimum (FindOptimal's search box is 0.33 degrees and fixed, it cannot undo more than about one
  degree), while the trap's cost is lower than that of a candidate stuck near the truth.
- **CSL class of the wrong answers.** Clean: Sigma3 126, Sigma5 31, Sigma7 26, Sigma11 3, Sigma13b 2,
  Sigma19b 2, Sigma21a, 25a, 29a 1 each, none within Brandon 11. Realistic: Sigma3 94, Sigma5 17,
  Sigma7 11, Sigma9 2, Sigma17a 2, Sigma19a 2, Sigma11 1, Sigma25a 1, none 13. Sigma3 is 62%
  and 66% of them.
- **Seed luck against systematic.** In clean data 28 voxels (14%) are wrong in all three seeds, 87 in
  some but not all, 85 never; under independent seeds with the pooled rate one expects 7.9, 134.6,
  57.5. Realistic: 13, 80, 107 against 2.7, 108.9, 88.4. So the failure is strongly voxel-dependent
  (systematic), and a second seed also helps a lot (F4).
- **Dependences.** The wrong rate does not depend on r_perp (clean quartiles 0.34, 0.33, 0.35, 0.34;
  realistic 0.31, 0.21, 0.21, 0.22), and not on proximity to a grain boundary (clean 31.8% near
  against 41.1% not near; realistic 23.1% against 26.2%; if anything the opposite sign). With the
  fraction of the truth's eligible reflections that are shared with at least one Sigma relative
  (`e0_dependence.txt`), the dependence is weak: Spearman +0.09 (Sigma3) and +0.10 (Sigma7) in clean,
  +0.11 and +0.12 in realistic data (p about 0.01; over runs, not independent), wrong rate by tertile
  of the Sigma3 shared fraction 0.27, 0.36, 0.38 (clean) and 0.21, 0.19, 0.31 (realistic); no
  dependence for Sigma5, 9, 11 or the number of eligible spots. The shared fraction varies little
  between voxels (0.60-0.76 for Sigma3), so this is not a strong predictor.

# 5. Cheap fixes (E1)

All rows: wrong rate with Wilson 95% interval; cost = cost-function evaluations per answer
(global + local; the baseline is about 49,000 clean and 52,000 realistic, mostly the global
radius-3 evaluations at level 0) and measured wall time. Plot: `fixes_wrong_rate_vs_cost.png`. Table:
`fixes_summary.txt`.

Which sets: baseline, F1 and F4 use the E0 runs (F1 on all 3 seeds, 600 per variant; "seed 0" rows
200); F1b and the classifier rows are re-runs on all 200 voxels, seed 0 (the same rng as E0); F2 and
F3 are re-runs of seed 0 on a subset (all baseline-wrong seed-0 runs plus an equal number of
randomly drawn baseline-right runs, 121 clean and 117 realistic cases) to save compute, and their
overall rate is implied: (still-wrong fraction of the baseline-wrong runs x number of baseline-wrong
runs + broken fraction of the baseline-right runs x number of baseline-right runs) / 200, with a
conservative interval (both components at their interval ends at once).

| Fix | clean wrong rate | realistic wrong rate | extra evaluations (clean / realistic) | extra time (s) | notes |
|---|---|---|---|---|---|
| baseline (seed 0, n=200) | 33.5% [27.3, 40.3] | 26.0% [20.4, 32.5] | 0 | 0 | median error of right answers 0.030 / 0.028 deg |
| F1 CSL check after the search, Sigma<=11 (n=600) | 3.8% [2.6, 5.7] | 3.5% [2.3, 5.3] | +2.7k / +2.9k | +1.3 / +1.5 | fixes 181 / 122, breaks 0 / 0 |
| F1, Sigma<=29 (n=600) | 3.5% [2.3, 5.3] | 3.3% [2.2, 5.1] | +5.3k / +5.4k | +2.6 / +2.7 | fixes 183 / 123, breaks 0 / 0 |
| F1, Sigma<=29, seed 0 (n=200) | 1.5% [0.5, 4.3] | 4.0% [2.0, 7.7] | +5.3k / +5.4k | +2.7 / +2.5 | |
| F1b CSL expansion at each level, seed 0 | 0.0% [0.0, 1.9] | 1.5% [0.5, 4.3] | +13.8k / +13.5k | +10.9 / +10.9 | fixes 67 / 50, breaks 0 / 1 |
| F2a keep 1/2 per level (cap 60) | 11.0% implied (still wrong 22/67) | 5.5% implied (11/52) | +7.9k / +12.0k | +5.1 / +7.6 | breaks 0 / 0 |
| F2b keep union of top 1/4 by post-MC cost and top 1/4 by discrete score | 11.0% implied (22/67) | 9.6% implied (17/52) | +2.7k / +4.1k | +1.9 / +2.9 | breaks 0 / 1 |
| F3a coarse levels at Q_max 8 | 15.2% implied (18/67) | 11.8% implied (10/52) | +15.2k / +23.6k | +18.7 / +23.3 | breaks 5/54 and 6/65 right runs |
| F3b coarse pixel radius 1 | 56% implied (61/67) | 40% implied (40/52) | -2.5k / -3.9k | -1.2 / -2.2 | breaks 21/54 and 18/65: worse than baseline |
| F4 best of 2 seeds (by local cost) | 20.5% [15.5, 26.6] | 11.0% [7.4, 16.1] | +49k / +52k | +24.5 / +26.5 | one answer per voxel, n=200 |
| F4 best of 3 seeds | 14.0% [9.9, 19.5] | 6.5% [3.8, 10.8] | +98k / +105k | +49 / +53 | |
| F4+F1: best of 2 / 3 F1 answers | 0.0% [0.0, 1.9] / 0.0% | 1.5% [0.5, 4.3] / 0.5% [0.1, 2.8] | +60k / +114k (clean) | +30 / +57 | |

(The extra time of F1 and F4+F1 rows is estimated from the evaluation counts at the measured time
per evaluation; the re-run fixes' times are measured.)

**What the results say.**

- **F1 is the clear winner on cost.** For a wrong answer that is a CSL relative of the truth, generating
  all its distinct Sigma relatives (Sigma<=29: 264 orientations), running the same quick MC on each,
  and refining the best three per Sigma subset recovers the truth: 90% (clean) and 86% (realistic) of
  the wrong runs fixed, 0 of 1200 right runs broken, for +11% evaluations and +10% time. Sigma<=11
  (38 relatives) gets almost all of it. In 6 of the 21 clean and 2 of the 20 realistic runs still wrong
  afterwards a relative's quick-MC result was within 3 degrees of the truth but did not rank first by
  cost. About half of the original wrong answers of the remaining runs (11 of 21 clean, 13 of 20
  realistic) have no Sigma <= 29 relation to the truth; those are not relatives of the truth, so F1
  cannot reach them (`f1_remaining.txt`; their E0 classes: all S2 in clean data, S2 12, S3 5, S1 3
  in realistic data).
- **F1b (put the relatives into the candidate set at every level) is the best single search change**
  (0/200 and 3/200 wrong) at +28% evaluations. Its advantage over F1 (1.5% and 4.0% on the same seed)
  is within the intervals at n=200.
- **F2: keeping more helps**, but less than F1: keep 1/2 per level cuts the rate to about 11% / 5.5%
  at +15-23% evaluations; the cheaper dual-ranking union (F2b) is as good in clean data and worse in
  realistic data. These address the cause (S2), not the symptom, but a wider cut still loses
  basin candidates that rank very low.
- **F3: the coarse-cost hypothesis (TODO hypothesis 1) is only partly right.** Starting at the full
  Q_max (F3a) helps (about 15% / 12%) at +30-45% evaluations and breaks some right runs (9% of the
  sampled right runs); a smaller pixel tolerance (F3b) is much worse: the radius-3 tolerance is what
  lets the coarse grid find the basin at all. The coarse-cost settings are not the cause; the
  pruning rank is.
- **F4 (best of N seeds by local cost)** helps a lot but costs N times; it is dominated by F1.
  With F1 on top, 2 seeds reach 0.0% / 1.5%.
- Median error of right answers: 0.030 to 0.026 degrees in all rows that change nothing in the final
  stage; F1b improves it slightly (0.020 / 0.018) because it adds candidates near the truth.

# 6. Classifier feasibility (2a)

**Data (E2).** For each of the 400 (voxel, variant) cases, features and labels of (i) every candidate
of E0 seeds 0 and 1 (post-quick-MC candidates of all levels and FindOptimal's results: 242,396
candidates) and (ii) synthetic hard examples: the truth perturbed by residuals (24 per case) and every
distinct Sigma<=29 relative of the truth perturbed by a residual (270 per case), residual angles drawn
from the E0 error of the nearest candidate to the truth at levels 0-2 (n=486, median 0.8, 90th
percentile 3.0 degrees, random axis). Total 359,996 candidates; positives (< 3 degrees) are 4.9%.

**Features (62), from one forward pass per candidate** (about 6 ms): number of eligible (peak,
detector) pairs; hit fractions with the centre pixel, centre or any triangle vertex, +-1 and +-3 pixel
boxes (pixel radius 0, 1 and 3); the same by |q| family (8 families), by detector (2); and CSL-aware
features for Sigma in {3, 5, 7, 9, 11}: fraction of the candidate's reflections shared with at least
one relative, hit fraction on the shared reflections, and over the relatives the minimum and mean hit
fractions on the reflections NOT shared with that relative (a trap matches the shared subset and
misses the rest). Plus the two cost features (local cost, radius-3 cost). Models: logistic
regression and histogram gradient-boosted trees (scikit-learn), trained to predict "within 3
degrees", balanced class weights, fixed hyperparameters.

**Protocol.** Four folds of about 50 voxels split by GRAIN (both variants of a voxel and all voxels of
a grain together: the 200 voxels come from 158 grains, and 28 grains have voxels in more than one
voxel-disjoint fold, which could leak orientations and CSL relatives; the first analysis used
voxel-disjoint folds, kept in `*_voxel_disjoint_folds.*` and within 0.004 AUC and 0.004 pruning recall
of the grain-disjoint numbers below). A model is evaluated only on grains it never saw. Evaluation sets: A "contested" (candidates that matter: ranking below twice the
typical cut, the last level, FindOptimal results; basin < 1 degree against trap > 3 degrees; n=71,388),
B synthetic basin against synthetic CSL relatives (n=117,600), C basin against candidates within 3
degrees of an exact CSL relative (n=96,823), D pruning recall (per level group with a basin candidate:
is the best one within the n_keep best by the score; 2134 groups), E final-choice precision (per run,
among FindOptimal's results with at least one basin and one trap candidate: is the top-scored one a
basin candidate; 514 groups). Baseline: the local cost alone (`e2_summary.txt`).

| Score | AUC A | AUC B | AUC C | pruning recall D | final precision E |
|---|---|---|---|---|---|
| local cost (no training) | 0.963 | 0.851 | 0.969 | 0.900 | 0.965 |
| logistic regression | 0.991 | 0.930 | 0.994 | 0.980 | 0.928 |
| gradient-boosted trees | 0.997 | 0.937 | 0.998 | 0.985 | 0.986 |
| gradient-boosted trees, no cost features | 0.997 | 0.936 | 0.998 | 0.985 | 0.992 |

Roc curves: `e2_roc.png`. By trap class (basin against candidates near an exact relative), the cost's
AUC is lowest for the classes that dominate the failures: Sigma3 0.896, Sigma5 0.923, Sigma7 0.958
against gradient-boosted trees 0.986, 1.000, 1.000 (all other Sigma 0.97-1.00 for the cost and 0.99-1.00 for
the trees).

**Ablation (`e2_ablation.txt`, trees, grain-disjoint).** All 62 features: AUC A 0.9967, pruning recall
0.985. Without the cost features: 0.9970, 0.985. Without the CSL-aware features: 0.9936, 0.985 (AUC B
0.927 against 0.937). Without any radius-1 or radius-3 feature: 0.9936 and 0.965 (AUC B drops to
0.853). Only the five overall hit fractions (centre, any vertex, +-1 and +-3 boxes, and the number of pairs):
0.9915, 0.982. So the gain over the cost comes mainly from combining the exact hit rate with the
tolerant hit rates, not from the CSL-specific features (they add about 0.003 AUC) and not from the
cost features.

**Cross-variant.** Trained on clean voxels and tested on realistic held-out grains: trees AUC A 0.989,
pruning recall 0.993, but final precision 0.831 (cost 0.949) on the realistic test data; trained on
realistic and tested on clean: AUC A 1.000, pruning 0.981, final precision 1.000. A classifier trained
only on clean data is fine at pruning and worse than the cost at the final choice.

**Learning curve (`e2_learning_curve.png`)**, trees (and logistic regression), held-out grains, mean
over folds and random subsets of the training voxels: pruning recall 0.977 (0.976) at 25 training
voxels, 0.982 (0.980) at 50, 0.985 (0.979) at 100, 0.985 (0.980) at 150; AUC A 0.988 (0.989), 0.989
(0.991), 0.995 (0.991), 0.997 (0.991); final precision (trees) 0.929, 0.965, 0.976, 0.988. The baseline
cost has 0.900 and 0.963. Twenty-five voxels (50 cases, both variants) already give most of the gain
at pruning; AUC and the final-choice precision keep improving up to 100-150 voxels (spread across
folds and subsets about +-0.01 to +-0.03, so the late differences are small).

**End to end (`e2_endtoend.txt`).** The model of the fold not containing the voxel's grain; seed 0, the
same random stream as E0; each candidate scored costs 3 evaluation equivalents (two cost evaluations
and a geometry pass; measured time +2-3 s).

| | clean | realistic |
|---|---|---|
| baseline, seed 0 | 67/200 = 33.5% | 52/200 = 26.0% |
| (a) classifier rerank at each pruning step | 8/200 = 4.0% [2.0, 7.7], +2.0k evals, +2.1 s; fixes 59, breaks 0 | 9/200 = 4.5% [2.4, 8.3], +3.2k, +3.3 s; fixes 46, breaks 3 |
| (a) + F1 | 4/200 = 2.0% [0.8, 5.0] | 7/200 = 3.5% [1.7, 7.0] |
| F1, seed 0 (for comparison) | 3/200 = 1.5% [0.5, 4.3] | 8/200 = 4.0% [2.0, 7.7] |
| F1b, seed 0 (for comparison) | 0/200 [0.0, 1.9] | 3/200 = 1.5% [0.5, 4.3] |
| (b) choose the final answer by the classifier among the answer and F1's refined relatives (600 each) | 20/600, identical to F1 | 25/600 = 4.2% against F1's 20/600 = 3.3%; breaks 5 |

**Answer to 2a.** A one-pass classifier does make the search more robust where it matters: used
as the pruning ranking it removes almost 90% of the failures (33.5% to 4.0%, 26.0% to 4.5%), at
+4-6% evaluations, because it separates the basin candidate (1-2 degrees off, cost near 1) from
traps far better than the local cost (AUC 0.997 against 0.963; pruning recall 0.985 against
0.900). But on these data it does not beat the non-learned fix: F1 reaches 1.5% / 4.0% on seed 0
(3.5% / 3.3% over 3 seeds) at +11% evaluations, F1b 0.0% / 1.5%, and the combination of classifier and
F1 (2.0%, 3.5%) is not distinguishable from F1 alone at n=200. For the final choice the cost already
separates basin from trap (final precision 0.965 against 0.986 for the best classifier on 514 groups;
end to end no gain, and slightly worse on realistic data), so the classifier adds nothing there. Its
value is that it addresses the cause (S2), not only CSL relatives; it would matter more if real data
produce non-CSL traps. That could not be tested here (see caveats).

# 7. Training data requirements (2b)

Backed by the learning curve of Section 6 (25 voxels already give most of the pruning gain;
100 give the full AUC and final-choice precision) for this one sample and geometry.

- **Labels.** The truth orientation of the voxel, hence simulation or a fully characterised sample: each
  training example is a candidate orientation plus the images of its voxel, labelled by the
  symmetry-reduced angle to the truth ("basin" under 1-3 degrees, "trap" otherwise; the Sigma class is
  secondary information, useful for stratifying but not needed). Experimental data without ground truth
  can be used only through a surrogate truth (a reconstruction confirmed by an independent measurement or by
  a much more expensive search), which carries the surrogate's own errors.
- **Candidate distribution.** Search-generated candidates are the right negatives: 239,434 level candidates
  (1.6% within 1 degree) and 2,962 FindOptimal results (37% within 1 degree) from the instrumented search,
  because they are what the ranking sees. Synthetic additions are cheap and fill the Sigma classes the
  search happens not to visit (108,000 CSL relatives with coarse-level residuals, 9,600 perturbed truths);
  on the synthetic examples (set B) the classifier reaches AUC 0.94 against the cost's 0.85. I did not train
  on harvested-only or synthetic-only data, so the relative value of the two sources is not measured.
- **Class balance.** 4.9% positives overall; balanced class weights were used; no resampling was tried.
- **Volume.** 200 voxels (400 cases with two image variants, 360,000 candidates, 6 ms of forward model
  each, about 40 CPU minutes to compute the features, plus the E0 runs: about 10 CPU hours for the
  instrumented reconstructions) were used; the curve says 25-50 voxels give most of the gain at the pruning decision (recall 0.977-0.982 against 0.985 with
  150), and 100-150 are needed for the full AUC and the final-choice precision. This is one
  sample type; do not read it as a requirement for a different geometry.
- **Diversity needed** (not varied here, so an inference): orientations across the fundamental zone
  (random grains: covered), positions and r_perp (covered, no dependence found), voxel sizes (all one
  size: not covered), Q_max (8 only), detector geometry (one 2-detector layout; the features
  by detector and |q| family are geometry-specific, the hit fractions are not), image realism (neighbour and
  grain-mate spots, noise, dropped/grown/hot pixels: both clean and realistic were mixed, and a model
  trained on clean data only lost on the final choice), and full-sample renders with real overlap of spots
  between voxels (not covered: images here are per voxel).
- **Held-out protocol.** By voxel at a minimum, but by grain (a grain's voxels share the
  orientation and all its CSL relatives, so they must stay together: done in the final analysis, the 200
  voxels come from 158 grains; the voxel-disjoint split gave nearly the same numbers), and by sample for any claim about
  a new material or geometry.
- **Domain-shift risks for real data and mitigations.** (1) Different noise, background and
  detector response than the synthetic images: train on measured-noise-like realism and calibrate the hit
  tests to the data's pixel statistics; prefer normalised hit fractions over raw counts. (2) Different Q_max
  or number of reflections: retrain, or drop the |q|-family features (the ablation shows the overall
  tolerant hit fractions carry most of the gain). (3) Voxel size and grain-boundary voxels: the voxel
  footprint changes the hit statistics. (4) Systematic positional errors (detector distance, tilt) shift
  all spots at once: the tolerant radii help, and the classifier should be trained with the geometry
  uncertainty simulated. (5) A model trained on one search configuration is evaluated inside a different
  candidate distribution when it changes the search (as the rerank does): the end-to-end test here shows
  it works for this search; re-check after any change of the search settings. Mitigation for all:
  ship the non-learned F1 check as the safety net, since the cost always separates basin from trap at the
  final choice.

# 8. Recommendation

1. **Add F1 (CSL check after the search) now.** About 150 lines, no training data, +11% evaluations,
   removes 86-90% of the failures (34% to 3.5% clean, 24% to 3.3% realistic, on 600 runs each), 0 right runs
   broken. Sigma<=11 gives almost all of it for half the cost.
2. **Consider F1b where robustness matters more than 30% extra time** (0.0% clean and 1.5% realistic wrong
   on 200 voxels). Whether it is significantly better than F1 needs more voxels.
3. **Cheap structural improvement: keep more candidates per level** (F2a, keep 1/2: to about 11% / 5.5%
   at +15-23% evaluations), or at least rank the survivors by something that tolerates the 1-2 degree
   offset of a coarse candidate. Do not reduce the global pixel tolerance (F3b is much worse). Do
   not run extra seeds as the main remedy (F4: 2 to 3 times the cost for less than F1).
4. **The classifier is not needed for the CSL failure.** It is the better fix of the cause (S2) and is
   cheap to evaluate, but on this data it is no better than F1; it is worth a trial only if real data
   show failures that are not CSL relatives (the remaining failures after F1 are non-CSL: 21/600 and 20/600) and only
   with training data that includes realistic data. Never use it as the final chooser; the cost does
   that job.
5. **A change to test next (inferred, not tested):** the ablation says five overall hit fractions
   (exact and tolerant) carry most of the trained model's gain (pruning recall 0.982 against 0.985 for
   all features, with a trained tree on those five); a hand-made pruning key that combines the exact and the
   tolerant (radius 1 and 3) hit rates might give much of the gain without a trained model.

# 9. Caveats

- One sample (500 grains, one 2-detector geometry, copper, Q_max 8, |sin eta| >= 0.3); per-voxel
  images, not a full-sample render; neighbour interference in realistic data is the sweep's
  distractor model, not real overlap. Real data may differ in the failure mix.
- The 3 seeds per voxel are not independent samples of the problem; intervals on rates over runs are
  optimistic; the F4 rows have one answer per voxel (n=200).
- F2 and F3 were run on a wrong+right subset (seed 0), and their overall rates are implied estimates with
  conservative intervals. F3c (F3a and F3b combined) was not run; F2b/F2a/F3 on seeds 1-2 and combinations
  of F2 with F1 were not run.
- The classifier was trained on candidates from E0 (the baseline search); the end-to-end rerank
  runs show it works inside the search, but only for seed 0 and 200 voxels; the fold models were
  trained on the other grains (about 150 voxels). Hyperparameters were fixed, not tuned. The CSL-aware features helped little.
- Cost units: one "evaluation" is one cost-function call; global (radius-3, partial Q) and
  local calls cost about the same in wall time; classifier scoring is counted as 3 evaluations
  per candidate, measured time gives +2-3 s per run. The runtime of F1 and F4+F1 rows is estimated from evaluation counts.
- The S-classes use the 3-degree basin and the post-quick-MC orientation; other thresholds were not
  scanned. Candidates within 3 degrees but beyond 1 degree count as "present".
- The analysis assumes the data images are those of the voxel itself; in a BFS reconstruction
  FindOptimal starts from neighbours' orientations, which this study does not exercise.

# 10. Reproduction

From `icenine_py/` (details: `scripts/findoptimal_robustness/README.md`):

```bash
uv run python scripts/findoptimal_robustness/e0_run.py select
uv run python scripts/findoptimal_robustness/e0_run.py run --workers 10 --n-voxels 204   # ~1 h
uv run python scripts/findoptimal_robustness/analyze_e0.py
uv run python scripts/findoptimal_robustness/f1_run.py run --workers 10                   # ~12 min
uv run python scripts/findoptimal_robustness/e2_dataset.py run --workers 10               # ~4 min
uv run python scripts/findoptimal_robustness/analyze_dependence.py
uv run python scripts/findoptimal_robustness/s2_mechanism.py
uv run python scripts/findoptimal_robustness/fixes_run.py run --fix F1b --seeds 0 --workers 10
for F in F2a F2b F3a F3b; do uv run python scripts/findoptimal_robustness/fixes_run.py run --fix $F --seeds 0 --workers 10 --subset-wr; done
OMP_NUM_THREADS=2 uv run python scripts/findoptimal_robustness/e2_models.py run --save-models   # ~12 min
OMP_NUM_THREADS=2 uv run python scripts/findoptimal_robustness/e2_ablation.py
uv run python scripts/findoptimal_robustness/e2_endtoend.py final --model GBT --workers 10
uv run python scripts/findoptimal_robustness/e2_endtoend.py rerank --model GBT --workers 10    # ~25 min
uv run python scripts/findoptimal_robustness/f1_run.py run --workers 10 --source e2_rerank_GBT
uv run python scripts/findoptimal_robustness/summarize_fixes.py
uv run python scripts/findoptimal_robustness/f1_remaining.py
uv run python scripts/findoptimal_robustness/e2_summary.py
uv run python scripts/findoptimal_robustness/plots.py
```

Raw per-run caches (about 300 MB) are in `scripts/findoptimal_robustness/cache/` (gitignored). Committed
results are in `benchmarks/findoptimal_robustness/`: `e0_summary.{txt,json}`, `e0_runs.npz`,
`e0_dependence.txt`, `s2_mechanism.txt`, `fixes_summary.{txt,json}`, `f1_remaining.txt`, `e2_summary.txt`,
`e2_results.json`, `e2_ablation.txt`, `e2_endtoend.{txt,json}`, `*_voxel_disjoint_folds.*` (first, voxel-disjoint run of E2), `e2_roc.npz`, `voxels.npz`,
`residual_pool.npy` and the four plots `where_lost.png`, `fixes_wrong_rate_vs_cost.png`, `e2_roc.png`,
`e2_learning_curve.png`.
