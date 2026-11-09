# Autoresearch stage (Track A): plan

Owner mandate: D-20261109-N1 (2026-11-09). Run by Newt Agent on task branches of
`feature/eval-reduction`. The saved state is `docs/autoresearch/STATE.md`; read it first after any
restart. Ideas set aside earlier are in `docs/todo_autoresearch.md` (Ilium's file, idea A1).

## The goals, in the owner's words

1. Reduce both the total number of operations and the total wall-clock time of the reconstruction
   of a strain-free sample, under the current physical assumptions.
2. Look for deep-learning approximations that accelerate the grid search, using the repetitive
   nature of the problem.
3. The voxel-by-voxel assumption may be dropped.
4. A lower-quality approximation followed by a final optimization is fine.
5. General early stopping: predict whether further optimization will improve a voxel.

## Where the operations go today

Terms. An **evaluation** is one call of the cost function for one voxel and one orientation. A
**full search** finds a voxel's orientation with no starting guess. A **local fit** starts from a
neighbour's orientation. A **ray** is one diffracted beam, known from a spot seen on both detectors.
**Wrong** means more than 1 degree from the truth.

Reading the table: one row per part of the 2,000-voxel BFS pilot (key-off rerun, 60 grains), with
its count, its evaluations and its share of the wall time.

| part | count | evaluations | share of time |
|---|---|---|---|
| full searches (40 seeds, 29 reseeds) | 69 | 22.88 M | 87% |
| local fits (neighbours and revisits) | 4,049 | 2.66 M | 13% |
| total | | 25.55 M | |

So one full search per grain is nine tenths of the work. A full search costs about 332,000
evaluations: 61% at grid level 0 (screen 44,000, re-cost 35,000, quick MC 123,000), 26% at level 1,
9% at level 2, 3% at level 3 and 1% in the finishers (six voxels of the G2 ledger). A local fit
costs about 660. The whole layer has 24,570 voxels and 497 grains. Sources:
`docs/findings_eval_reduction_2026-10.md` lines 109-146 and 425-430; `G/g2/ledger/runs/`.

Three earlier results make a different route possible:

- **Rays are known before any search.** Task D pairs a spot on detector 0 with its spot on
  detector 1 (`icenine/spot_pairing.py`): 98.7% of one-grain pairs are found and 88.7% of the
  matches are right. A pair gives the beam's direction, so its scattering vector, without knowing
  the voxel. It also gives a line in the sample that the beam came from (2.7 µm median error).
  A spot belongs to a grain, not to a voxel.
- **A start within 2 degrees is enough.** Task F's Gauss-Newton step on the rays (`RayGaussNewton`)
  reaches the truth from 1 and 2 degrees in 100% of cases, in a few hundred evaluations.
- **One orientation can be costed over the whole layer at once.** S4 showed a spot's pixel position
  is linear in the voxel position, so one orientation's spots are computed once and shifted per
  voxel (6.58 s for the whole level-0 grid on 10 threads).

## Hypotheses

H1. **Index the rays; do not grid-search.** Two scattering vectors of one grain fix its orientation.
Take the rays whose lines pass near a seed voxel, assign each to a lattice-plane family by its ring
angle, and let every pair with an allowed angle between planes vote for an orientation. The
orientations with the most votes are the candidates. Check the best few with the strict cost and
finish with the usual finisher. Expected: a full search of about 1,000 evaluations in place of
332,000. Task D named this ("seeding orientations by indexing") and never tested it.

H2. **Index the layer, then paint grains (no voxel-by-voxel search).** Vote with all rays of the
layer at once. Each cluster of votes is a grain: its orientation, its rays, and from the rays' lines
its rough place and extent. Then cost each grain's orientation over all voxels near it in one
batched pass and give each voxel the best grain. Only voxels with no clear winner get a local fit.
Expected: about 5 batched evaluations per voxel in place of today's roughly 12,800 (25.55 M over
2,000 voxels).

H3. **A stop rule from what is left unexplained.** For a voxel and its current orientation, two
numbers say whether more search can help: the hit ratio (how many predicted spots are lit), and the
best vote count of any other orientation among the rays near the voxel that the current
orientation does not explain. If nothing else gets votes, more search cannot find a better
orientation. H3 has a part that needs no new run: on the stored E0 results, how well do the final
cost and the per-level costs separate right from wrong voxels, and how many evaluations would a
stop at each level have saved?

H4. **A learned proposer (deep learning).** Where the voting fails (grains with few or mixed
spots), train a network that reads the rays near a voxel as a set of points (family, direction in
the sample frame) and returns likely orientations as a ranking over the 4,886 grid points of the
fundamental zone. Training data is free: the forward model at random orientations, with spots of
neighbouring grains added as clutter. The network is fitted once per setup and material and reused
for every voxel, layer and sample; that is where the repetition pays. It plugs into the same
check-then-finish path as H1. It is measured against H1 on the same voxels, so we learn whether the
network adds anything over plain voting. Earlier networks here re-ranked or refined a start (the
toy orientation net, the cost proxy); none proposed orientations from the data.

H5. **Order of the layer-wide screen and the BFS** (idea A1). Kept as a fallback arm for the voxels
where H1, H2 and H4 give no candidate.

## How every arm is judged

- Test beds: E0 (200 hard voxels), the 2,000-voxel BFS pilot region, and the whole layer for the
  arms that pass. Clean full-sample images. Reconstruct at Q-max 8.
- Report, side by side: wrong counts (paired: fixed and broken against the baseline), evaluations
  in strict-equivalent units, and wall time with `FusedCost` on, single-worker under the timing
  preflight. No p-values. No extra random seeds.
- For a proposer (H1, H4): recall at K, the share of voxels whose true orientation is within
  2 degrees of one of the first K candidates, for K = 1, 5, 20.
- Baselines: the key-off BFS pilot rerun (37 wrong, 25.55 M evaluations); the same with `FusedCost`
  on (review item 14); the full-sample run (item 5); and a measured C++ time (the owner's question
  on finding 6).
- New switches are opt-in and off by default. No default flips and no golden regeneration.

**Library changes:** new modules only (`icenine/ray_indexing.py`, later `icenine/grain_paint.py`
and `icenine/ray_proposer_nn.py`), and opt-in keys in `reconstructor.py`, `config_file.py` and
`orientation_search.py` whose off path is bit-identical. Steps AR0 and AR1 change nothing under
`icenine_py/icenine/` except adding `ray_indexing.py`.

**Stop rule:** stop and report if a step's run passes twice its stated time budget, if a baseline
does not reproduce (wrong count or evaluation count differs from the stored one on unchanged
code), or if the plan cannot be followed as written.

# Steps

## AR0: what the stored results say about stopping

No new search runs. Read the stored E0 results (`benchmarks/eval_reduction/E/cache/rank/e0_*/v*.npz`
in the main checkout; fields `err_deg`, `cost_final`, `evals`, and per level `L{k}_cost`,
`L{k}_dist`, `L{k}_screen_recost`, `L{k}_quick`) and the pilot's `records.json` and `recon.mic`.
Answer with a plot first, then a table: (a) does the final cost separate right from wrong voxels,
and at which threshold; (b) after each level, is the true basin already the best candidate, and how
many evaluations come after that point; (c) for a rule "stop after level k when the best cost is
at most c", the evaluations saved against the voxels broken, for a grid of k and c. Scripts in
`scripts/autoresearch/ar0/`, results in `benchmarks/autoresearch/ar0/`, findings in
`docs/autoresearch/findings_ar0_stopping.md`. Budget: 30 minutes of compute.

## AR1: recall of the ray-indexing proposer

Write `icenine/ray_indexing.py`: from `SmoothCost.observed_rays()` and
`attached_rays(voxel_vertices)`, compute each ray's scattering vector in the sample frame, assign
families by ring angle, vote over pairs, and return the top K orientations (reduced to the
fundamental zone, duplicates merged). Symmetry comes from the crystal's own operators, nothing
cubic-specific. Measure on E0 and on the pilot's 69 full-search voxels: recall at K = 1, 5, 20
within 2 degrees, the candidates' distance to the truth, the time per voxel, and for the misses the
reason (too few rays, mixed spots, wrong pairs). Use the forward model's `peak_hook` for the true
grain of each spot. Then check candidates end to end: strict cost on the top K, finisher on the
best, wrong count and evaluations against the stored full-search results on the same voxels.
Scripts in `scripts/autoresearch/ar1/`, results in `benchmarks/autoresearch/ar1/`, findings in
`docs/autoresearch/findings_ar1_ray_indexing.md`, unit tests in `tests/test_ray_indexing.py`.
Budget: 2 hours of compute.

## AR2: indexing as the BFS seed search

Only if AR1's recall at 20 is at least 90% on the pilot's full-search voxels. An opt-in key
(`SeedByIndexing`) replaces the grid search of seeds and reseeds with the AR1 path and falls back to
the grid search when no candidate passes the seed test. Run the pilot against the key-off rerun,
paired. Budget: 2 hours.

Added 2026-11-09 after AR1 and its review (D-20261109-N7). The gate is met (69 of 69). AR2 has two
steps, in this order.

**Step 1, holdout of the finisher arms.** AR1 chose "finisher on the best 2" and "on the best 3"
after it saw the result of "on the best 1", on the same 69 voxels. So their 0 of 69 is not a test.
The test: voxels of the layer outside the pilot region, picked by a fixed stride over the voxel
index (no random draw), in two sets fixed before the run: 150 boundary voxels (a voxel with a
neighbour of another grain in the true map) and 150 interior voxels. Index each voxel, take the
strict cost on the first 5 candidates, then run the finisher on the best 1, the best 2 and the
best 3. Report, per set and arm: wrong voxels against the truth with intervals, recall at 1, 5 and
20, and evaluations. The main arm for step 2 is best 2. It is kept unless it has more wrong voxels
than best 3 on the holdout; then best 3 is used. No stored grid search exists for these voxels, so
this step gives wrong counts and evaluations of the indexing path only. Budget: 30 minutes.

**Step 2, the key.** `SeedByIndexing` (off by default; the off path stays bit-identical) in
`config_file.py` and `reconstructor.py`. The rays are paired once per run. A seed or reseed is
indexed, its first 5 candidates get the strict cost, the finisher runs on the best 2 (or the arm
step 1 picks), and the result goes through the same seed test as today. If no candidate passes,
the grid search runs as today. The pilot region is run with the key on and compared with the
key-off rerun on the same voxels: wrong counts (fixed, broken, each break traced, boundary
crossings in their own column), evaluations, the number of fallbacks, and wall time. For wall time
both arms run with `FusedCost` on, one worker, under the timing preflight; the key-off arm of that
pair is the pilot run of review item 14. Indexing time and the one-time pairing are reported beside
the evaluations, since they are outside the evaluation count. Budget: 2 hours.

## AR3: index the layer and paint grains

Only if AR1 passes. Layer-level voting, grain list, batched painting, local fits for unclear
voxels. Run on the pilot region, then the whole layer. Compare with the full-sample baseline.

## AR4: the learned proposer

Train and measure H4 on the same voxels as AR1. It starts after AR1 so that it can be compared
with the voting and aimed at the voxels the voting misses.

## AR5: the stop rule

Combine AR0's thresholds with the unexplained-ray vote count of H3; test as a stop rule inside the
arms that survive.

## Machine-time runs from the review page

Run when the machine is otherwise quiet: item 14 (both keys on E0 and the pilot, with the timing
A/B of item 3), item 5 (the full-sample baseline), item 18 (does C++ fuse the clip), and a measured
C++ time for the pilot region.
