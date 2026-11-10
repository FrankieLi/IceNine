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

Added 2026-11-10 after AR2 step 2 and its review (D-20261110-N1). The four lines above are the
first draft and are kept for the record. The plan below replaces them.

**Why this step.** With the key `SeedByIndexing` on, the grid searches are gone. Of the 2,819,906
evaluations of the pilot, 2,697,708 are local fits: 4,095 fits at about 658 each. So the target
is now the local fits. And every wrong voxel so far is a boundary voxel that took the orientation
of a grain it touches. The BFS gives a voxel to the first grain that reaches it and passes the
test. It never asks whether the other grain fits better.

**The idea.** The sample is strain-free, so one grain has one orientation. Once a grain's
orientation is known from its seed, a voxel of that grain needs no fit. It needs one evaluation,
to check that the orientation fits there. Terms:

- A **grain entry** is one accepted seed: its voxel, its orientation and its hit ratio. The hit
  ratio is the share of predicted spots that are lit (the BFS's own acceptance number).
- To **score** a grain on a voxel is one strict evaluation of the grain's orientation on that
  voxel, with no fit. It gives a cost and a hit ratio.
- A grain **passes** on a voxel by the BFS's own neighbour test, unchanged: hit ratio above 0.9
  times the best hit ratio seen so far in that grain's flood.
- To **flood** a grain is to score it on the neighbours of its seed, then on the neighbours of
  every voxel where it passes, and so on. The neighbourhood is the BFS's (radius 2 side lengths).
  A flood runs through every voxel where the grain passes, also where another grain passes. So
  one grain's flood does not depend on any other grain.
- A voxel is **contested** when two or more grains pass on it. It is **unpainted** when no grain
  passes on it.
- To **paint** is to give each voxel the passing grain with the lowest strict cost.

AR3 has two steps, in this order.

**Task split.** Task `ar3a` is step 1 only. Step 2 is a later task (`ar3b`); it is not started
before step 1 is reviewed and merged, and no part of it belongs in `ar3a`. For the implementer of
either task: the tool refuses a helper agent's write of a new `.md` file. So do not write the
findings page or the hand-back report as files. Return each as text in the final message, the
findings page between the lines `=== FINDINGS PAGE ===` and `=== END FINDINGS PAGE ===`, the
report between `=== REPORT ===` and `=== END REPORT ===`. The main session places both and says
so in them. Everything else (code, tests, scripts, result files, plots) is committed on the task
branch as usual.

**Step 1, painting on the pilot region (a study; no key yet).**

What is built. A new module `icenine/grain_paint.py` with unit tests in
`tests/test_grain_paint.py` (a fake cost on a small grid with two grains and a known boundary).
Nothing else under `icenine_py/icenine/` changes. The module uses the pieces the key already
uses: the seed path `ray_indexing.seed_by_indexing` (index, strict cost on 5 candidates, finisher
on the best 2), the reconstructor's strict cost built once per run, and its local fit. The loop:

1. Take the voxels in the same shuffled order as the BFS runs of AR2 (the order the config's seed
   gives; no new random draw).
2. The next voxel that is unpainted and has not been a seed becomes a seed. If at least one grain
   was scored on it, it first gets one local fit from the lowest-cost of those grains, as the BFS
   would. If that fit passes the neighbour test, the voxel is settled "by a fit" and the loop goes
   on. Otherwise it is indexed by the seed path. If the result passes the seed test (the same as
   today), it becomes a grain entry and is flooded. If not, the grid search runs as with the key
   (the fallback), and a voxel whose search is rejected is failed, as today.
3. Repeat until every voxel is painted, settled by a fit, a seed, or failed.
4. Paint. Record for every voxel each grain scored on it (cost, hit ratio, passed or not), the
   stage that settled it, and the evaluations by stage.

A new seed whose orientation is within 1 degree of an earlier grain entry is still its own entry.
Such pairs are counted and reported. Distances use the crystal's own symmetry operators.

Before the paint run, two short measurements:

- **Where the 213 s go.** Rerun the key-on pilot (`FusedCost` on, one worker, timing preflight)
  under a profiler. Report the time by part: cost evaluations inside local fits, the optimizer's
  own work, indexing, pairing, set-up, the rest. This says what a batched pass could save.
- **Is the premise true on this sample.** From the true map: the number of distinct orientations
  inside each true grain of the region. Painting assumes one.

The arms. All on the 2,000-voxel pilot region, clean images, Q-max 8.

| arm | what | source |
|---|---|---|
| key off | BFS, grid search on seeds | stored run of the key-off rerun (37 wrong, 25,547,935 evaluations) |
| key on | BFS, seeds indexed, a local fit on every other voxel | stored run of AR2 step 2 (22 wrong, 2,819,906 evaluations) |
| paint | the loop above; a contested voxel takes the lowest cost | new, one run |
| paint + fit the contested (side arm) | as paint, then on every contested voxel one local fit from each of its two lowest-cost grains; the lower final cost is kept | new, same run, its evaluations counted apart |

Reading the table: one row per arm. The two stored arms are not rerun, except for the profile.

How each is judged.

- **Evaluations** (the main number): paint against key on and key off, split by stage (seed
  paths, scoring, local fits, fallbacks). One score is one evaluation. Indexing time and the
  one-time pairing are reported beside them, as in AR2.
- **Wrong voxels**, against the truth, with intervals, and paired against the key-on run: fixed,
  broken, both wrong. Each break is traced. Boundary crossings are in their own column, split
  into edge-sharing and corner-only (N6, N8, N9). Two runs of the same key-off arm differ by 20
  voxels one way and 7 the other. So a difference of that size, either way, is not called a gain
  or a loss.
- **The direct test of the boundary failure.** On every contested voxel: is the grain with the
  lowest strict cost the true grain? Counts with an interval, split into edge-sharing and
  corner-only boundary voxels and interior voxels. This asks whether the cost itself can tell the
  two touching grains apart when both are scored. It does not depend on the BFS route. Also: how
  often the lowest cost and the highest hit ratio name different grains.
- **How far a wrong grain floods.** For each grain, the voxels where it passes outside its true
  grain, by distance to the true boundary in voxels.
- **Accuracy without a fit.** The distance to the truth of painted voxels against the same voxels
  in the key-on run (which fitted each one): the plot of both distributions first, then the
  table. This says whether a final fit per voxel, or one per grain, is needed (goal 4).
- **The side arm.** For a margin m, a rule "fit a contested voxel when its two lowest costs are
  closer than m" has a cost (the fits) and a gain (the voxels fixed). The side arm gives both for
  every m. The curve is read on this run, so it picks an m; it does not test it. The test is the
  rest of the layer, in step 2.
- **Wall time**: the paint run with `FusedCost` on, one worker, under the timing preflight,
  beside the profiled key-on run of the same session. The paint loop scores one voxel at a time
  in Python. So the wall time is a first reading of an unbatched loop, labelled as such.

The gate for step 2, fixed now: go on unless paint breaks more than 20 voxels against the key-on
run, or needs more than half of the key-on run's evaluations.

What is not claimed in step 1: nothing about the whole layer (the pilot pairs rays with the hull
of its own region); nothing about samples whose grains hold more than one orientation; nothing
about realistic images or another Q-max; no margin rule (it is only picked here); no tested wall
time; nothing about another seed order (the loop is run with one fixed order); nothing about the
fallback if it runs on no voxel again. A layer-wide vote over all rays (the first draft of H2) is
not built: grains are listed by indexing unpainted voxels one at a time. If the seed paths turn
out to be the largest share of the paint arm's evaluations, a cheaper seed (fewer finisher runs,
or one finisher run per grain over several voxels) is the next thing to try, in step 2 or AR5.

Scripts in `scripts/autoresearch/ar3/`, results in `benchmarks/autoresearch/ar3/`, findings in
`docs/autoresearch/findings_ar3_paint.md`. Budget: 1 hour of compute.

**Step 2, the key and the whole layer.** Only if the gate is met. Its details are fixed after
step 1 and its review, as AR2's were. In outline: an opt-in key `PaintGrains` (off by default; the
off path stays bit-identical) in `config_file.py` and `reconstructor.py`, which runs the loop of
`grain_paint.py` in place of the BFS. The pilot through the key must give step 1's numbers again.
Then the whole layer (24,570 voxels, 497 grains), rays paired with the layer's own hull, with the
key-on BFS (`SeedByIndexing`) on the whole layer as the arm to compare with. The margin rule
picked in step 1 is tested on the voxels outside the pilot region. Budget: 3 hours.

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
