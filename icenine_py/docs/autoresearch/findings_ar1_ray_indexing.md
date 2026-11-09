# AR1: can the rays be indexed instead of grid-searched?

Step AR1 of `docs/autoresearch/PLAN.md` (2026-11-09), hypothesis H1. Code:
`icenine/ray_indexing.py`. Scripts: `scripts/autoresearch/ar1/`. Results:
`benchmarks/autoresearch/ar1/`. The code, the runs and the numbers are the AR1 implementer's
(commits 3a35b30, 780dd97, d6dbbcd). This page was placed by the main session from the
implementer's draft, because the implementer's write of this file was refused by the tool.

## The result in one picture

Reading the picture. Left: for each K, the share of voxels whose true orientation is within
2 degrees of one of the first K candidates. Candidates are ranked by the number of rays they
explain. Right: evaluations per voxel for one full search, on a log axis. Grey is the stored grid
search. Blue is the indexing path set in advance: the strict cost on 5 candidates, then the usual
finisher on the best one. E0 is the easy case: its images hold only each voxel's own spots. The
pilot is the realistic test.

![Recall at K and evaluations per voxel](../../benchmarks/autoresearch/ar1/ar1_recall.png)

- On the pilot's 69 full-search voxels the truth is among the first 3 candidates for every voxel.
- The first candidate alone is right for 41 of the 69. In the other 28 it is the true orientation
  of another grain. In 27 of those 28 that grain touches the voxel.
- The indexing path set in advance takes a median of 1,138 evaluations per voxel. The stored grid
  search took 332,154 on the same voxels. The path gets 3 voxels wrong against 4 stored: it fixes
  4 and breaks 3. The 3 breaks come from the arm's own rule, so they count against it (step 5).
- The evaluation count leaves out three costs: the one-time pairing of all spots (17.069 s),
  attaching rays to a voxel (median 0.002 s) and indexing it (median 0.177 s). Indexing one voxel
  takes about as long as 400 strict evaluations. All times are contended (timed while other
  processes ran).
- The plan's gate for step AR2 was a recall at K = 20 of at least 90% on these voxels. It is met:
  69 of 69.

## Terms

- **Ray**: one diffracted beam. It is known from a spot on detector 0 paired with its spot on
  detector 1.
- **Attached rays** of a voxel: the rays whose line through the sample passes the voxel. Most of
  them belong to other grains.
- **Scattering vector** of a ray: the ray's direction minus the beam's direction, scaled by the
  wave number, turned into the sample frame with the rotation angle (omega) of the ray's frame.
  For a crystal with orientation R it equals R h for one reciprocal-lattice vector h.
- **Family**: the reciprocal-lattice vectors that the crystal's own symmetry rotations map onto
  each other. The code takes the rotations from the crystal and turns them into the Cartesian
  frame of its lattice; nothing in it is specific to cubic crystals. All measurements on this page
  are on one cubic copper sample. The hexagonal case is covered only by a unit test of the
  construction (hcp titanium: the vectors close under the rotations, and exact vectors index to
  the right orientation).
- **Vote**: two rays of different families give one orientation when the angle between their
  scattering vectors matches the angle between two crystal vectors of those families.
- **Explained ray**: a ray within 0.7 degrees and 1% in length of some R h.
- **Candidate**: a cluster of votes, refined by least squares on the rays it explains. Candidates
  are ranked by the number of rays they explain.
- **Recall at K**: the share of voxels whose true orientation is within 2 degrees of one of the
  first K candidates.
- **Wrong**: more than 1 degree from the truth.
- **Evaluation**: one call of the strict cost function for one voxel and one orientation.
- **E0**: 200 hard voxels. Its clean images hold only each voxel's own spots, so E0 is the easy
  case for indexing.
- **Pilot full searches**: the 69 voxels of the 2,000-voxel BFS pilot (key-off rerun) that took
  more than 100,000 evaluations. The images hold the spots of all grains, so this is the realistic
  test. By their final record 62 are seeds. The other 7 are seeds that the seed test rejected
  (`seed_rejected`): 5 were later revisited and 2 stayed unresolved. For those 7 the stored error
  is the later revisit or the failed fit, not the full search's own result.
- **Stored full search**: for E0 the arm `e0_old_base`; for the pilot the key-off rerun.
- **Touching grain**: a grain one of whose voxels shares a vertex with the voxel (from the
  simulation's grain map).
- **Boundary crossing** (decision D-20261109-N6): a boundary voxel that takes the orientation of a
  touching grain.
- **Contended**: timed while other processes ran on the machine.

All images are clean (no noise model). The pilot images hold the spots of all grains; the E0
images hold only each voxel's own spots. Reconstruction uses Q-max 8. Counts carry a Wilson 95%
interval in percent.

## Step 1: are the scattering vectors right?

This was first checked on a few voxels before any voting; the table below is the later full
rerun (see the end of this step).

Reading the table: one row per test set. `attached` is the number of attached rays over all
voxels. `t05` counts the rays within 0.5 degrees and 1% in length of a vector of the true
orientation. `c05` counts the same for a control orientation (the truth turned by 37 degrees).
`t1` and `c1` are the same within 1 degree. `ang50`, `ang95` and `angmax` are the median, the 95th
percentile and the largest angular error in degrees of the rays within 1 degree of the truth.
`len95` is the 95th percentile of their relative length error.

<!-- table:ar1_sanity -->
| set | voxels | attached | t05 | c05 | t1 | c1 | ang50 | ang95 | angmax | len95 |
|---|---|---|---|---|---|---|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | 200 | 11312 | 11270 (99.6%) | 3 (0.03%) | 11312 | 14 | 0.203 | 0.439 | 0.559 | 0.0026 |
| Pilot full searches (69 voxels) | 69 | 27164 | 3065 (11.3%) | 4 (0.01%) | 3115 | 28 | 0.206 | 0.451 | 0.985 | 0.0056 |
<!-- /table:ar1_sanity -->

The convention is right. On E0 nearly all attached rays match the truth. On the pilot 11% match,
as expected, since most attached rays belong to other grains. The control matches 4 of 27,164
rays within 0.5 degrees. The angular errors stay below 0.6 degrees on E0. That is consistent with
the 1-degree frame width (omega is known to half a frame); the share due to spot centroids was not
tested. The tolerances are 0.7 degrees to explain a ray, 1 degree for the angle between two
vectors, 1% in length and 1 degree for clustering. They were set from a sanity run on 10 E0 and 6
pilot voxels, before recall was measured. The table above is the full rerun, made after recall
and after the end-to-end runs, with the tolerances unchanged.

## Step 2: recall

Reading the table: one row per test set. `K1`, `K5` and `K20` give the voxels whose truth is
within 2 degrees of one of the first 1, 5 or 20 candidates. E0 is the easy case (own spots only);
the pilot row is the realistic test.

<!-- table:ar1_recall -->
| set | K1 | K5 | K20 |
|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | 200/200 (100.0%, 98.1-100.0) | 200/200 (100.0%, 98.1-100.0) | 200/200 (100.0%, 98.1-100.0) |
| Pilot full searches (69 voxels) | 41/69 (59.4%, 47.6-70.2) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
<!-- /table:ar1_recall -->

Reading the table: the same counts when "within 2 degrees" is replaced by 0.5, 2 and 5 degrees.

<!-- table:ar1_recall_thr -->
| set | within | K1 | K5 | K20 |
|---|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | 0.5 deg | 200/200 | 200/200 | 200/200 |
| E0 (200 voxels, easy: own spots only) | 2 deg | 200/200 | 200/200 | 200/200 |
| E0 (200 voxels, easy: own spots only) | 5 deg | 200/200 | 200/200 | 200/200 |
| Pilot full searches (69 voxels) | 0.5 deg | 41/69 | 69/69 | 69/69 |
| Pilot full searches (69 voxels) | 2 deg | 41/69 | 69/69 | 69/69 |
| Pilot full searches (69 voxels) | 5 deg | 41/69 | 69/69 | 69/69 |
<!-- /table:ar1_recall_thr -->

The counts do not change between 0.5 and 5 degrees. A candidate is at the truth or far from it.

Reading the table: `best_med`, `best_p90` and `best_max` are the median, the 90th percentile and
the largest distance in degrees from the truth to the nearest of the 20 candidates. `top1_med` is
the median distance of the first candidate where it is within 2 degrees. `n_cand_med` is the
median number of candidates.

<!-- table:ar1_dist -->
| set | best_med | best_p90 | best_max | top1_med | n_cand_med |
|---|---|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | 0.038 | 0.076 | 0.125 | 0.038 | 20 |
| Pilot full searches (69 voxels) | 0.054 | 0.082 | 0.137 | 0.048 | 20 |
<!-- /table:ar1_dist -->

Reading the table: `med_ms` and `p90_ms` are the indexing time per voxel in milliseconds
(contended; no cost function is called). The E0 times were measured with its 3 recall parts
running at once, and for part of that time the pilot end-to-end run as well. `med_votes` is the median number of votes. `med_att` is the
median number of attached rays. `med_true` is the median number of rays the truth explains.

<!-- table:ar1_index_time -->
| set | med_ms | p90_ms | med_votes | med_att | med_true |
|---|---|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | 91 | 101 | 5022 | 57 | 57 |
| Pilot full searches (69 voxels) | 179 | 234 | 28248 | 413 | 44 |
<!-- /table:ar1_index_time -->

A pilot voxel has about ten times more attached rays than its own grain explains. The voting
still finds the grain.

## Step 3: why the first candidate is sometimes not the truth

The forward model records which grain drew each spot. That gives the true owner of every ray.

Reading the table: medians over the 69 pilot voxels of the attached rays, by who drew them.

<!-- table:ar1_ray_kinds -->
| quantity | median |
|---|---|
| attached rays | 413 |
| drawn only by the voxel's own grain | 45 |
| drawn by the own grain and another grain | 1 |
| drawn by other grains only | 366 |
| explained by the true orientation | 44 |
<!-- /table:ar1_ray_kinds -->

Reading the table: the 28 voxels whose first candidate is more than 2 degrees from the truth, by
reason; whether the grain of the first candidate touches the voxel (shares a vertex, or a whole
edge); and the rank at which the truth appears.

<!-- table:ar1_misses -->
| what | n |
|---|---|
| rank-1 candidate is not the truth | 28/69 |
|   reason: too few rays (truth explains < 8) | 0 |
|   reason: another real grain explains more rays | 28 |
|   reason: matches no grain (coincidence of rays) | 0 |
| rank-1 grain touches the voxel (shares a vertex / an edge) | 27 / 16 |
| truth is candidate 1 / 2 / 3 (of 20) | 41 / 22 / 6 |
<!-- /table:ar1_misses -->

Every miss has the same reason. The first candidate is the true orientation of another grain
(within 1 degree), and it explains more of the attached rays than the voxel's own grain does: a
median of 5 more. Both counts are taken after the same least-squares refinement. One case is a
tie at 41 rays each; candidates are ranked by explained rays, then by votes, so the votes broke
it. In 27 of the 28 misses the other grain touches the voxel (16 share an edge with it). In the
remaining one (voxel 11294) it is one voxel side away. No miss comes from too few rays. No miss
comes from a chance alignment of rays. Mixed spots are rare: 65 of the 3,095 rays that the true
orientations explain.

So the count of explained rays did not tell the two grains apart in 28 of the 69 voxels. In the
3 voxels where the end-to-end arm then went wrong, the strict cost after the finisher was lower
for the true grain (step 5). The higher count of the other grain is consistent with its rays
being attached because their lines pass the voxel; the attachment rule was not varied.

## Step 4: fewer rays, and other tolerances

Two stress tests on the pilot voxels. First, keep only every m-th attached ray (a fixed stride,
no random draw). Second, halve and double each tolerance.

Reading the table: one row per case. `rays` is the median number of rays the truth explains. The
other columns are as in the recall table. `stride m` keeps every m-th ray. `ang`, `pair` and
`cluster` name the tolerance that was changed and its new value in degrees.

<!-- table:ar1_sparse -->
| case | rays | K1 | K5 | K20 |
|---|---|---|---|---|
| stride 1 | 44 | 41/69 (59.4%, 47.6-70.2) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| stride 2 | 23 | 42/69 (60.9%, 49.1-71.5) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| stride 3 | 15 | 37/69 (53.6%, 42.0-64.9) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| stride 4 | 11 | 33/69 (47.8%, 36.5-59.4) | 67/69 (97.1%, 90.0-99.2) | 68/69 (98.6%, 92.2-99.7) |
| stride 6 | 7 | 32/69 (46.4%, 35.1-58.0) | 62/69 (89.9%, 80.5-95.0) | 67/69 (97.1%, 90.0-99.2) |
| stride 8 | 5 | 30/69 (43.5%, 32.4-55.2) | 54/69 (78.3%, 67.2-86.4) | 66/69 (95.7%, 88.0-98.5) |
| ang 0.35 | 36 | 44/69 (63.8%, 52.0-74.1) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| ang 1.4 | 46 | 41/69 (59.4%, 47.6-70.2) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| pair 0.5 | 44 | 42/69 (60.9%, 49.1-71.5) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| pair 2 | 44 | 42/69 (60.9%, 49.1-71.5) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| cluster 0.5 | 44 | 42/69 (60.9%, 49.1-71.5) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
| cluster 2 | 44 | 42/69 (60.9%, 49.1-71.5) | 69/69 (100.0%, 94.7-100.0) | 69/69 (100.0%, 94.7-100.0) |
<!-- /table:ar1_sparse -->

With 5 true rays left (stride 8) the truth is still among 20 candidates for 66 of 69 voxels, and
among 5 for 54. The tolerance changes move the first-candidate count by at most 3 voxels (41 to
44) and leave K = 5 and K = 20 at 69 of 69.

## Step 5: end to end

The arms:

- the strict cost on the first 1, 5 or 20 candidates, with no finisher;
- "AR1, best 1": the strict cost on the first 5 candidates, then the usual finisher on the best
  one. This arm was set in advance;
- the same with the finisher on the best 2 and on the best 3, keeping the lowest final cost. Both
  were added after the best-1 result was seen. Neither was set in advance, so their counts are not
  a test.

Reading the table: `wrong` is the voxels more than 1 degree from the truth. `fixed` is wrong in
the stored search and right here. `broken` is right in the stored search and wrong here. `both` is
wrong in both. `mean_ev` and `med_ev` are the evaluations per voxel (strict cost on the candidates
plus the finisher; the stored count includes the screen's calls). Pairing, attachment and indexing
are not in these counts (see the end of this step). `med_err` is the median error in degrees.

<!-- table:ar1_e2e -->
| set | arm | wrong | fixed | broken | both | mean_ev | med_ev | med_err |
|---|---|---|---|---|---|---|---|---|
| E0 (200 voxels, easy: own spots only) | stored full search | 41/200 (20.5%, 15.5-26.6) | 0 | 0 | 41 | 48,700 | 48,724 | 0.021 |
| E0 (200 voxels, easy: own spots only) | strict cost on first 1, no finisher | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 1 | 1 | 0.038 |
| E0 (200 voxels, easy: own spots only) | strict cost on first 5, no finisher | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 5 | 5 | 0.038 |
| E0 (200 voxels, easy: own spots only) | strict cost on first 20, no finisher | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 20 | 20 | 0.038 |
| E0 (200 voxels, easy: own spots only) | AR1: strict on 5, finisher on best 1 (set in advance) | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 519 | 530 | 0.017 |
| E0 (200 voxels, easy: own spots only) | AR1: strict on 5, finisher on best 2 (chosen after seeing best 1) | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 950 | 951 | 0.017 |
| E0 (200 voxels, easy: own spots only) | AR1: strict on 5, finisher on best 3 (chosen after seeing best 1) | 0/200 (0.0%, 0.0-1.9) | 41 | 0 | 0 | 1,393 | 1,392 | 0.017 |
| Pilot full searches (69 voxels) | stored full search | 4/69 (5.8%, 2.3-14.0) | 0 | 0 | 4 | 332,380 | 332,154 | 0.015 |
| Pilot full searches (69 voxels) | strict cost on first 1, no finisher | 28/69 (40.6%, 29.8-52.4) | 2 | 26 | 2 | 1 | 1 | 0.073 |
| Pilot full searches (69 voxels) | strict cost on first 5, no finisher | 3/69 (4.3%, 1.5-12.0) | 4 | 3 | 0 | 5 | 5 | 0.054 |
| Pilot full searches (69 voxels) | strict cost on first 20, no finisher | 3/69 (4.3%, 1.5-12.0) | 4 | 3 | 0 | 20 | 20 | 0.054 |
| Pilot full searches (69 voxels) | AR1: strict on 5, finisher on best 1 (set in advance) | 3/69 (4.3%, 1.5-12.0) | 4 | 3 | 0 | 1,125 | 1,138 | 0.020 |
| Pilot full searches (69 voxels) | AR1: strict on 5, finisher on best 2 (chosen after seeing best 1) | 0/69 (0.0%, 0.0-5.3) | 4 | 0 | 0 | 1,926 | 1,907 | 0.019 |
| Pilot full searches (69 voxels) | AR1: strict on 5, finisher on best 3 (chosen after seeing best 1) | 0/69 (0.0%, 0.0-5.3) | 4 | 0 | 0 | 2,546 | 2,512 | 0.019 |
<!-- /table:ar1_e2e -->

Reading the table: the pilot rows of the table above for the stored search and the three
finisher arms, split by the stored record. `seeds` are the 62 voxels whose stored result is the
full search's own. `rejected seeds` are the 7 seeds that the seed test rejected; their stored
result is a later revisit or a failed fit, not the full search.

<!-- table:ar1_e2e_split -->
| subset | arm | wrong | fixed | broken | both |
|---|---|---|---|---|---|
| 62 seeds | stored full search | 2/62 (3.2%, 0.9-11.0) | 0 | 0 | 2 |
| 62 seeds | AR1: strict on 5, finisher on best 1 (set in advance) | 3/62 (4.8%, 1.7-13.3) | 2 | 3 | 0 |
| 62 seeds | AR1: strict on 5, finisher on best 2 (chosen after seeing best 1) | 0/62 (0.0%, 0.0-5.8) | 2 | 0 | 0 |
| 62 seeds | AR1: strict on 5, finisher on best 3 (chosen after seeing best 1) | 0/62 (0.0%, 0.0-5.8) | 2 | 0 | 0 |
| 7 rejected seeds | stored full search | 2/7 (28.6%, 8.2-64.1) | 0 | 0 | 2 |
| 7 rejected seeds | AR1: strict on 5, finisher on best 1 (set in advance) | 0/7 (0.0%, 0.0-35.4) | 2 | 0 | 0 |
| 7 rejected seeds | AR1: strict on 5, finisher on best 2 (chosen after seeing best 1) | 0/7 (0.0%, 0.0-35.4) | 2 | 0 | 0 |
| 7 rejected seeds | AR1: strict on 5, finisher on best 3 (chosen after seeing best 1) | 0/7 (0.0%, 0.0-35.4) | 2 | 0 | 0 |
<!-- /table:ar1_e2e_split -->

Reading the table: the 3 pilot voxels that the best-1 arm breaks. `taken` is the grain whose
orientation the arm returned. `boundary` says whether that grain touches the voxel, which makes
the break a boundary crossing. `cost_taken` and `cost_truth` are the strict costs of the two
candidates before refinement; `fin_taken` and `fin_truth` are their strict costs after the
finisher. `start_deg` is the distance in degrees from the truth to the finisher's start.
`against_best1` says whether the break counts against the arm.

<!-- table:ar1_broken -->
| voxel | taken | boundary | cost_taken | cost_truth | fin_taken | fin_truth | start_deg | against_best1 |
|---|---|---|---|---|---|---|---|---|
| 7141 | grain 23 | yes (shares an edge) | 0.183 | 0.201 | 0.177 | 0.080 | 49.5 | yes |
| 8848 | grain 395 | yes (shares an edge) | 0.187 | 0.190 | 0.139 | 0.028 | 43.8 | yes |
| 23942 | grain 286 | yes (shares an edge) | 0.271 | 0.335 | 0.209 | 0.111 | 54.4 | yes |
<!-- /table:ar1_broken -->

- E0, the easy case: no voxel is wrong in any arm. The stored search has 41 wrong; all 41 are
  fixed and none is broken. The median is 530 evaluations against 48,724.
- Pilot, best 1 (set in advance): 3 wrong against 4 stored. It fixes the 4 stored wrong voxels
  (5876, 8866, 15749, 24334) and breaks 3 others (7141, 8848, 23942). Two of the 4 fixed voxels
  (5876, 24334) are rejected seeds whose stored result is an unresolved fit, not a full search.
  Against the 62 seeds alone the arm fixes 2 and breaks 3.
- The 3 broken voxels. The arm picks one candidate by its strict cost before refinement, then
  refines only that one. In all 3 the other grain's candidate has the lower cost before
  refinement (0.183 against 0.201 in voxel 7141). After the finisher the true grain has the
  clearly lower cost in all 3 (0.080 against 0.177 in voxel 7141). So these breaks trace to the
  arm's rule (pick before refining), not to a preference of the strict cost for the other grain.
  Under the broken-voxel rule (D-20261014-N13) a break caused by the arm counts against it, so
  all 3 count against the best-1 arm (D-20261109-N7). The grain taken shares an edge with the
  voxel in all 3, so they are also boundary crossings (D-20261109-N6). That column does not excuse
  them here.
- Pilot, best 2 and best 3: all 69 right, at a mean of 1,926 and 2,546 evaluations against 1,125.
  Both were chosen after the fact, so they need a test on voxels they have not seen.
- The finisher's start. In the 66 voxels where the picked candidate is the truth, the finisher
  starts within 0.120 degrees of it and takes a median of 1,137.5 evaluations. In the 3 broken
  voxels it starts 43.8 to 54.4 degrees away, at the other grain, and takes 835 to 877. Its search
  box is fixed, whatever the start. Starts between 1 degree and the 43.8 degrees of the broken
  voxels were not tested.

Costs outside the evaluation count. The one-time pairing of all spots (17.069 s), attaching rays
to a voxel (median 0.002 s) and indexing it (median 0.177 s: voting, clustering, refining) are not
in the counts above. The finisher takes a median of 0.44 ms per evaluation. At that rate, indexing
one voxel takes as long as about 400 strict evaluations, about a third of the finisher's median
time (0.49 s). All these times are contended.

## Time

Reading the table: seconds on the 69 pilot voxels, one process. The timing is labelled contended:
the timing preflight refused the machine as busy (load 2.59, one other process at 97% CPU; see
`preflight_timing.json`).

<!-- table:ar1_timing -->
| quantity | seconds |
|---|---|
| one-time pairing of all spots (26,019 rays) | 17.069 |
| attach rays to a voxel (median) | 0.002 |
| index a voxel (median) | 0.177 |
| index a voxel (90th percentile) | 0.230 |
|   of which voting (median) | 0.021 |
|   of which clustering (median) | 0.085 |
|   of which refining (median) | 0.069 |
<!-- /table:ar1_timing -->

Clustering is the largest part. It is a Python loop over seeds; a faster version was not tried.
The times of the strict cost and the finisher are in the end-to-end JSON files, also contended.

## What this means for the plan

- Indexing passes the AR2 gate on this sample. Whether it can replace the grid search for seeds
  is AR2's test.
- Best 2 and best 3 were both chosen after the best-1 result was seen. AR2 will use the finisher
  on the best 2 as its main arm (D-20261109-N7). Its test is a holdout on voxels outside the pilot
  region, with the arms fixed in advance.
- A boundary voxel sees two real grains among its candidates: in 27 of the 28 misses the other
  grain touches the voxel. That is not an error of the indexing; both grains are there. It is the
  reason to cost more than one candidate, and it is the starting point of AR3: the candidates of a
  voxel are the grains around it.
- The comparison is in evaluations. The wall-time comparison of the whole path against the grid
  search, on a quiet machine, is still to do (AR2).

## Not tested

- Realistic images. The noise model is an open question; these are clean images.
- Another sample, material or geometry. This is one cubic copper sample and one region of 2,000
  of 24,570 voxels with 60 grains. Non-cubic crystals are covered only by a unit test of the
  construction (hexagonal).
- Voxels in small grains with few paired spots, except by thinning the rays.
- Local fits. Only full searches were indexed.
- The cause of the other grain's higher ray count. It is consistent with the attachment rule, not
  tested.
- The 3 broken of 69 are counts on one region, not a rate.
- The full test suite was not run (no existing library file changed).
