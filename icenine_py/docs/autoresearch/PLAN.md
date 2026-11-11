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

**Library changes:** new modules only (`icenine/ray_indexing.py`, later `icenine/grain_paint.py`,
`icenine/gradient_sample.py` (AR3 step 3, D-20261110-N4)
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

**The idea.** Once a grain's orientation is known from its seed, a voxel of that grain can
start from it with no search. One evaluation checks that the orientation fits there. On the
simulated test sample, which has exactly one orientation per grain, that start is already the
answer and no fit is needed. On a real sample it is a starting approximation, and a refinement
of each voxel follows (step 3).

Corrected 2026-11-10 on the owner's word (D-20261110-N3). This paragraph first read "The sample
is strain-free, so one grain has one orientation." That is not true of a real sample: a
strain-free grain can hold micro-gradients of orientation, so every voxel in a grain may have
its own orientation. Steps 1 and 2 were planned and run under the first wording. Their results
hold for a sample with one orientation per grain and for nothing else.

Terms:

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
before step 1 is reviewed and merged, and no part of it belongs in `ar3a`. Added 2026-11-10:
step 1 is merged (8e3098f), so task `ar3b` is step 2 as written under "Step 2, the details"
below; step 1 is not redone. Task `ar3c` is step 3 as written under "Step 3" below; it starts
after `ar3b` is merged and redoes nothing of steps 1 and 2. Task `ar3d` is step 4 as written
under "Step 4" below; it starts after `ar3c` is merged and reruns nothing of step 3. Task
`ar3e` is step 5 as written under "Step 5" below; it starts after `ar3d` is merged and reruns
none of step 4's whole-layer runs. For the implementer of
any of these tasks: the tool refuses a helper agent's write of a new `.md` file. So do not write the
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

**Step 2, the details.** Added 2026-11-10 after step 1 and its review (D-20261110-N2). The gate
is met: paint broke 0 voxels against the key-on run (limit 20) and used 5.6% of its evaluations
(limit 50%). Step 1's page is `docs/autoresearch/findings_ar3_paint.md`. Three things from step 1
shape this step:

- The side arm fixed 0 and broke 0 voxels at every margin. So step 1 picked no margin. The rule
  it leaves is "no fit on a contested voxel". Step 2 tests that rule.
- No local fit was accepted (0 of 48). So the rule "a voxel settled by a fit keeps its fit" was
  never used.
- The key-on BFS map holds only 64 distinct orientations over 2,000 voxels. So its local fits
  nearly always seem to return their start. They are 2,697,708 of its 2,819,906 evaluations.

What is built.

1. **The key.** `PaintGrains` (0 or 1, default 0) in `config_file.py`, carried as `paint_grains`
   in the search parameters of `orientation_search.py`, beside `seed_by_indexing`. With the key
   on, `BFSReconstruction.reconstruct_sample` runs the loop `grain_paint.paint_grains` through
   `ReconstructorHooks` in place of the BFS. It uses the voxel order and the random generator
   that the BFS would use. It writes what the BFS writes: the map, one record per voxel (source
   `seed`, `paint`, `fit` or `failed`) and the evaluations by stage in `stats`. The paint result
   stays on the object so that a script can save it. The loop in `grain_paint.py` is not changed:
   the rule for a fit-settled voxel stays as in step 1.
2. **What the key refuses.** `PaintGrains 1` needs `SeedByIndexing 1` (the seed path is the
   indexing). With `SeedByIndexing 0`, or with `OutwardSearch 1`, the set-up raises a
   `ValueError` that names the keys. The BFS's revisits, restart pass and reseeding do not run
   under the key; the README says so. Only the single-process path is built.
3. **Tests** in `tests/test_paint_key.py`: the key parses, default off; with the key off the
   paint loop is never called (the loop is patched to raise) and the result of a small BFS run
   is unchanged; the two refused settings raise; with the key on, a small real case gives the
   same map and the same evaluations by stage as a direct call of `paint_grains` with the same
   order and generator, and the stages sum to the cost function's own count.
4. **Study records, no change to the loop.** The hooks keep, per seed, what the seed path already
   returns (`SeedResult`: the strict costs of the costed candidates, the finisher's final costs,
   its evaluations, the candidate taken). A script wraps the BFS's local fit to log each fit
   (see "The fit log" below). Both live in the scripts or in the hooks, not in the BFS.

The runs. Clean images, Q-max 8, `FusedCost` on, one worker, one fixed voxel order. Runs longer
than a few minutes are started in the background with a log (as
`scripts/autoresearch/ar2/key_run_bg.sh` does) and watched with `scripts/dev/job_status.sh`. The
two whole-layer runs go one after the other, not at the same time, so that their wall times can
be read side by side.

| run | sample | what | compared with |
|---|---|---|---|
| A | pilot, 2,000 voxels | paint through the key | step 1's stored paint run |
| B | pilot, 2,000 voxels | key-on BFS with the fit log | the stored key-on run (`key_on_fused`) |
| C | whole layer, 24,570 voxels | paint through the key, then the side arm | run D and the truth |
| D | whole layer, 24,570 voxels | key-on BFS with the fit log | run C and the truth |

Reading the table: one row per run. A and B are checks that must pass before C and D start. C
and D are the result. The whole layer pairs rays with the layer's own hull, so the pilot voxels
inside C and D do not see the same rays as in A and B.

The checks, in this order. If one fails, stop and report; do not go on to C and D.

- **A gives step 1 again**: 156,766 evaluations with the same split by stage, 62 entries, 0 wrong
  voxels, and the same orientation on every voxel as `runs/paint/recon.mic`.
- **B is the stored run**: 2,816,897 evaluations and the same orientation on every voxel as the
  stored `key_on_fused` map. This shows two things: the off path is unchanged by the key, and
  the fit log does not disturb the run.
- The full test suite passes.

How C and D are judged. Every count is given for the whole layer and split into the pilot region
(2,000 voxels, where the method was worked out) and the rest (22,570 voxels, the test).

- **Evaluations** (the main number): paint against the key-on BFS, by stage, as in step 1.
  Indexing time and the one-time pairing beside them.
- **Wrong voxels** against the truth, with intervals, and paired against run D: fixed, broken,
  both wrong. Boundary crossings in their own column, edge-sharing and corner-only apart. Each
  break is traced, or a sample of 20 if there are more.
- **The wording rule for a paired count** (fixed here, for this page and later ones). A paired
  count is given as numbers: fixed, broken, both wrong. It is called a gain or a loss only when
  the two error rates' intervals do not overlap. Next to it, each fixed voxel is put in one of
  two bins by the scores of run C: the grain that run D gave it passes in paint and loses on
  strict cost, or it does not pass at all. The run-to-run difference of one arm on the whole
  layer is not measured. The pilot's (20 and 7 voxels of 2,000) is quoted beside the counts as
  the scale. No extra seeds are run.
- **The direct test.** On every contested voxel whose passing grains are more than 1 degree
  apart: is the lowest strict cost the true grain? Count and interval, split into edge-sharing,
  corner-only and interior. Contested voxels between two entries of one grain are counted apart;
  they cannot fail.
- **The rule "no fit on a contested voxel".** The side arm of step 1, run on the contested voxels
  whose passing grains are more than 1 degree apart: one local fit from each of the two
  lowest-cost grains, the lower final cost kept. Fixed and broken against run C, evaluations
  counted apart. The rule holds if the side arm fixes no more than it breaks.
- **The grain list.** Entries against the 497 true grains: grains with no entry, and entries
  within 1 degree of an earlier one (duplicates) with what they cost (seed path and flood). For
  the two duplicates of the pilot and for those of the layer: which test let the voxel become a
  seed (read from the loop's log: the hit ratio of the rejected fit against the grain's best,
  and the seed's own hit ratio). No merge rule is built here.
- **Voxels no grain paints**: unpainted, failed, and every grid-search fallback that ran, each
  listed with its outcome. This is the first test of the fallback on real voxels, if it runs.
- **Fit-settled voxels.** How many fits are accepted. For each: does a grain pass there with a
  lower strict cost than the fit's? The wrong count is given under both rules (the fit is kept;
  the lowest cost wins). Both come from the scores of run C, with no new evaluations. The key
  keeps step 1's rule; the page says which rule the numbers favour.
- **Accuracy without a fit.** Distance to the truth, run C against run D on the same voxels: the
  plot first, then the table.
- **Wall time.** C and D with the timing preflight, labelled as it comes out (single-worker, or
  contended with the reason). The loop still scores one voxel at a time.

The fit log (the first reading for an early stop, goal 5 of the owner's mandate). In runs B and
D, each local fit of the BFS is logged: the start and the end orientation, the angle between
them (with the crystal's symmetry operators), the final cost, the evaluations used, and the
strict cost and hit ratio of the start. The start's score is one extra evaluation on a second
cost function; it is counted apart and does not touch the run's generator. Reported, plot first:

- the share of fits by how far they moved (not at all; under 0.01 degrees; 0.01 to 0.1; over
  0.1), with the evaluations spent in each bin;
- the cost fall against the start's cost;
- for the fits that moved more than 0.1 degrees: was the voxel right at the end?

This is a reading on one sample whose grains each hold exactly one orientation. It says how much
of the BFS's work could have been skipped here. It picks a stop rule for AR5; it does not test
one. Paint is the simplest such rule run for real (score once, keep the start if it passes).

The seed record (a reading for a cheaper seed). Seed paths were 76% of paint's evaluations in
step 1, at 1,920 per seed. From the per-seed record of run C: how often the candidate with the
lowest strict cost is also the one the finisher ends on; how far the finisher moved the taken
candidate; how much the cost fell. This says what "finisher on one candidate" or "a shorter
finisher" would have changed. It is read on this run only and is tested later (AR5).

Stop rules. If check A or B fails, stop. If run D is not done after 120 minutes of wall time,
stop it, say so, and judge run C against the truth alone. If the whole-layer pairing does not
fit in memory, report that and stop. Budget: 3 hours of compute.

What is not claimed in step 2: nothing about realistic images or another Q-max; nothing about
another voxel order (one fixed order); nothing about samples whose grains hold more than one
orientation; nothing about other layers; no tested stop rule and no tested cheaper seed (both
are readings); no tested wall time unless the preflight passes; nothing about a parallel run.

Files. Scripts in `scripts/autoresearch/ar3/`, results in `benchmarks/autoresearch/ar3/layer/`,
findings in `docs/autoresearch/findings_ar3_layer.md` (returned as text, see "Task split"). The
README gets the key; `MIGRATION_HISTORY.md` gets the summary.

### Step 3: paint as a start, then refine, on grains with gradients

Added 2026-11-10. Task `ar3c`. It starts after `ar3b` (step 2) is merged.

**Why.** The owner corrected the premise on 2026-11-10 (D-20261110-N3). A strain-free grain can
hold micro-gradients of orientation, so every voxel in a grain may have its own orientation.
The test sample of steps 1 and 2 has exactly one orientation per grain. So their two main
readings, "no fit is needed" and "the BFS fits rarely move", were never able to fail. The
owner's design: painting gives each voxel a start, and a refinement of each voxel follows.
Step 3 asks three things on a sample where grains hold gradients. How good is the painted
start? What does the refinement cost and give? When can the refinement stop (goal 5 of the
mandate: predict whether more optimization will give a better orientation)?

**Terms.**

- The **gradient sample** is a copy of the layer's truth map with the same voxels and the same
  voxel-to-grain map. Each grain keeps its old orientation at its centre. Inside the grain the
  orientation turns smoothly with position.
- The **spread** of a grain is the misorientation between the orientation at its centre and the
  orientation of the voxel farthest from the centre.
- The **level** of a grain is the spread chosen for it: 0, 0.1, 0.2, 0.5 or 1.0 degrees.
- The **start** of a voxel is the orientation painting gives it.
- The **refinement** of a voxel is one local fit of the BFS (`local_optimization`, unchanged
  settings) from its start.
- The **error** of a voxel is the misorientation between its orientation and its own true
  orientation, with the crystal's symmetry operators. "Wrong" stays: error above 1 degree.
- A voxel is **hurt** by a shortcut when its error with the shortcut is more than 0.05 degrees
  above its error without it.

**The sample.** New module `icenine/gradient_sample.py` with one pure function. Inputs: voxel
centres, the grain id of each voxel, each grain's orientation, each grain's level, a numpy
generator. For grain `g` with centre `c` (the mean of its voxel centres) it draws a 3 by 2
matrix `M` of standard normal numbers. The rotation vector at position `x` is `M (x - c)`,
scaled so that its largest length over the grain's voxels equals the level. The voxel's
orientation is that small rotation applied to the grain's orientation; the docstring states on
which side. A grain of one voxel, and a grain of level 0, keep their orientation exactly.
Levels are drawn per grain, uniformly from the five values, in grain-id order, generator seed
0; one draw, no second seed. Tests in `tests/test_gradient_sample.py`: the measured spread
equals the level; level 0 is unchanged bit for bit; the centre orientation is kept; two calls
with the same seed agree; positions and grain ids are untouched.

A script writes the gradient truth map, a grid-only copy for the reconstruction (orientation
columns zeroed, as for the stored sample) and a small JSON (seed, level of each grain, measured
spread, sha256 of the truth map). The two `.mic` files are committed only if the pre-commit
size check lets them through; no override is used. The JSON and the script are committed.

**The images.** Rendered with `scripts/phase_d/render_full.py` (the Python forward model that
made the stored images): clean only, same geometry and Q-max as the images of runs C and D of
step 2, a new tag `grad`. If the script cannot take another truth map, add an option; without
the option it must do what it did. The rendering took about 48 s for the stored sample.

**Reading 0: what an offset costs** (uses the truth; for analysis only, in no decision). For
each voxel, one strict score of its true orientation and one of its grain's centre
orientation. By level: the hit ratio and the cost of both, and the ratio of the two hit ratios
against the neighbour test's 0.9. This says how far a flood can reach before it must stop.

**The runs.** Clean images, Q-max 8, `FusedCost 1`, one worker, the voxel order of step 2.

- **P**: paint through the key (`SeedByIndexing 1`, `PaintGrains 1`), as run C of step 2.
- **R**: the refinement of P's map, a script on `ReconstructorHooks`. In the run's voxel order,
  every voxel that a grain painted gets one `local_fit` from its start. Seeds, voxels settled
  by a fit and failed voxels keep what P gave them. Per voxel: the start's cost and hit ratio,
  the entry's best hit ratio, the end orientation, cost and hit ratio, the evaluations, and the
  trace.
- **R2**: a second fit from R's end, on every fourth voxel of the voxel order (fixed here). Same
  record.
- **D**: the BFS with indexed seeds (`SeedByIndexing 1`, `PaintGrains 0`) with step 2's fit log.

The **trace** of a fit: after each block of the Monte Carlo fit, the evaluations so far, the
best cost, the best orientation, and whether the block failed. `MCOptimizer.optimize` already
takes a `trajectory` list; use it if it gives this, else a wrapper in the script as for step
2's fit log. Check: the first 200 fits give the same result with the trace on and off (same
generator state). No change to the library's fit is planned. If one is needed, it is off by
default, the off path is bit-identical, and a test shows it.

P and D may run side by side. R and R2 follow P. Every wall time is labelled as the timing
preflight says.

**The split for every rule that is picked.** Grains with an even id are the pick set A. Grains
with an odd id are the test set B. A rule's threshold is chosen on A and then reported on B
with the threshold fixed. Fixed here, before any run.

**Judging.** Every count is given by level. The plot comes before the table. A difference
between two arms is called a gain or a loss only when the two rates' intervals do not overlap.

1. **The start.** Error of the painted orientation: median, 95th percentile, the share above
   0.1 degrees and the share wrong. Where the voxels came from (a grain, a fit, a seed,
   failed). Entries per true grain. Evaluations by stage. This is the first run in which the
   fit-settled branch and the extra seeds are expected to carry real weight.
2. **The refinement.** Error after R against the start and against D, paired on voxels: the
   same four numbers, and how many voxels the fit brings closer by more than 0.05 degrees and
   how many it hurts. Evaluations of P, of P plus R, and of D.
3. **Stop before the fit** (goal 5). Two rule families, both from run-time numbers only. "Skip
   the fit when the start's hit ratio is at least t times the entry's best hit ratio" (t from
   0.90 to 1.00 in steps of 0.01). "Skip the fit when the start's strict cost is at most c" (c
   on the deciles of A's start costs). On A: the rule that skips the most evaluations while it
   hurts at most 1% of A's painted voxels. On B, with that threshold: the evaluations skipped,
   the share hurt with its interval, and the four error numbers of the resulting map.
4. **Stop inside the fit** (goal 5). From the trace: how many of a fit's evaluations come
   after it is within 0.02 degrees of its own end. Then two rule families, simulated on the
   trace at block ends. "Stop at the first block that fails." "Stop when the best cost is at
   most c" (the fit's own `max_convergence_cost`; c on the deciles of A's end costs). Picked on
   A and reported on B as in point 3. A simulated stop is exact for its own voxel. In a real
   run the shared generator would then differ for later voxels; say so.
5. **Does a second fit help?** (goal 5 in the owner's words). From R2: the share of voxels a
   second fit moves by more than 0.05 degrees, and the share it brings closer to the truth by
   more than 0.05 degrees. Then one threshold on the first fit's end hit ratio over the entry's
   best, picked on A, that says "a second fit will not help"; on B: how often it is right.
6. **The BFS under gradients.** D's fit log by level: how far the fits move, the error, the
   evaluations. P plus R against D on the four error numbers and on evaluations.
7. **The control.** Level-0 grains inside P and D against the same voxels in step 2's C and D.
   Their neighbours differ, so this is a reading, not a repeat.
8. **Wall time**, labelled.

**Stop rules of the task.** If the sample tests fail, or Reading 0 shows no difference between
the levels (the images did not take the gradient), stop. P is stopped after 90 minutes, R
after 60, D after 150; a stopped run is reported as far as it got. Budget: 4 hours of compute.

**Not claimed in step 3.** The gradient model is a choice (linear in position, five levels),
not measured on a real sample; the owner is asked which sizes and shapes are typical. Nothing
about sharp sub-grain boundaries. Clean images, one layer, one voxel order. The stop rules of
points 3 to 5 are simulated from records; a later step runs the chosen one inside the loop.
No tested wall time unless the preflight passes.

**What comes after** (not part of `ar3c`). Step 4 is planned from these readings: the chosen
stop rule run for real inside the loop; a grain-wide model of the orientation field fitted to
a few refined voxels and scored once on the others (this drops the voxel-by-voxel assumption,
goal 3), on a sample that also has sharp sub-grain boundaries so that the model can fail; the
cheaper seed of step 2's seed record. Note from step 2's fixes (2026-11-10): on the 15 seeds
where the finisher's end was not the lowest-cost candidate, that candidate's own finisher run
ends 14 to 56 degrees from the truth. So "finisher on the lowest-cost candidate only" is not
safe as it stands; it needs the seed test's verdict on those ends, which was not recorded.

**Practical notes for the implementer.** Step 2 is merged (741e8fc): reuse its scripts
(`scripts/autoresearch/ar3/run_layer.py`, `run_layer_bg.sh`, `ar3_common.py`, the fit log of
`install_fit_log`) and its run configs (`benchmarks/autoresearch/ar3/layer/runs/C` and `D`) as
the templates for P and D. The image data is gitignored and lives in the main checkout: a
worktree needs a symlink `Examples/Example2.ManyGrains/ScatteringData_PhaseD` to the main
checkout's directory of that name; the new images go below it under the tag `grad`. Do not
commit the symlink, the images or any cache. The spot pairing of the new images is new work
(about 20 s); do not reuse step 2's pairing cache for them.

**Files.** `icenine/gradient_sample.py`, `tests/test_gradient_sample.py`; scripts in
`scripts/autoresearch/ar3/`; results in `benchmarks/autoresearch/ar3/gradient/`; findings in
`docs/autoresearch/findings_ar3_gradient.md` (returned as text, see "Task split"). The README
gets the new module; `MIGRATION_HISTORY.md` gets the summary. Only the result files that the
summary script reads are committed (D-20261110-N4); a per-voxel record above 2 MB is stored
gzipped and read gzipped.

### Step 4: a local fit whose step matches the error

Added 2026-11-10 after the review of step 3. Task `ar3d`.

**Why.** Step 3 found that the local fit moves almost none of the painted voxels. The review
found the reason in the stored records and in the code. The cost does tell the truth from the
painted start: the true orientation costs less than the start on 16,591 of 19,013 painted
voxels, and on every voxel that starts more than 0.1 degrees off. But the fit searches a box of
5 degrees. Its first trial steps are about 0.7 degrees, and after each failed block it restarts
5 to 17 degrees away. The painted starts are 0.04 degrees off (median). A trial of that size
almost never lands in the small region near the truth. So the fit spends about 660 evaluations
and returns its start. The C++ source has the same logic, and all 43 config files in this
repository set the same box (LocalOrientationGridRadius 5, MinLocalResolution 0).

**The idea.** Start the fit with a box that matches the error of its start, and widen it only
when that fails. A box of 5/8 degree gives a first trial step of about 0.09 degrees. By the
step rule of the code, a fit in a small box runs about 46 evaluations where the stock fit runs
662. So a finer box may be both closer to the truth and about 14 times cheaper. Both are to be
measured; neither is claimed here.

**Terms.**

- The **box** is the side of the region that the local fit searches (`search_box_side` of
  `MCOptimizer.variance_minimizing_optimize`). The **stock** box is 5 degrees.
- A **fine fit** is the same fit with a smaller box. Nothing else changes.
- A **ladder** runs a fine fit first and a wider one only when the fine one failed.
- **Start**, **error**, **wrong**, **level**, **hurt**, pick set **A** and test set **B** are as
  in step 3.

**Part 1: fits from the stored painted starts.** No whole-layer run. Reuse step 3's refine
stage: P's stored paint and generator state, the 19,013 painted voxels of the gradient sample,
the same images. Every arm starts every voxel from the same painted start and from the same
generator state. One run per arm.

| arm | the fit |
|---|---|
| `B5` | the stock fit, box 5 degrees. This is step 3's run R. It is read from the stored record and not rerun. |
| `B06` | box 5/8 degree |
| `B03` | box 5/16 degree |
| `B015` | box 5/32 degree |
| `L` | a ladder: box 5/32; when the cost did not fall, box 5/8; when it still did not fall, the stock box. The evaluations of every rung count. |
| `C` | the existing CMA local fit (`local_optimizer` "cma"), first step 0.05 degrees, at most 60 evaluations. Run it only if it can be called from a script with no library change; otherwise leave it out and say so. |

The box is passed to the optimizer by the script. The config key MinLocalResolution is not
used for this, because it also changes the Q range of the adaptive search and the trial box.

Readings of part 1, for each arm, on all voxels, by level and by the size of the start's error:

1. How many voxels came closer to the truth by more than 0.05 degrees, and how many were hurt.
2. Median and 95th percentile of the error, the share above 0.1 degrees, the wrong count.
3. The share of fits that lowered the cost. The share whose end costs no more than the true
   orientation (the truth is used here for analysis only).
4. Evaluations per fit (median and total), with the final evaluation counted as in step 3.

**The pick.** On A: the arm with the lowest share of voxels above 0.1 degrees after the fit.
Among arms whose intervals overlap with that one, the arm with the fewest evaluations. The pick
is then reported on B. It is called a gain over the start, or over `B5`, only when the
intervals on B do not overlap.

**Gate.** If no arm is a gain over the start on B, stop after part 1 and report. Parts 2 and 3
are then not run.

**Part 2: the stop readings again, on a fit that moves.** With the picked arm, redo step 3's
readings 3, 4 and 5 with step 3's summary code: skip the fit from what is known before it; stop
inside the fit from its trace; a second fit on every fourth voxel and whether it helps. Rules
are picked on A and reported on B. A pick that sits at the edge of its grid is called
uninformative, in the section where it appears.

**Part 3: the whole layer with the picked fit.** Still no library change. A script wrapper
replaces the fit inside the reconstructor's `local_optimization` and nowhere else: not in the
seed path, not in the grid search, not in the finisher. The implementer checks this by counting
the calls and says so on the page. The wrapped fit is "fine first, stock as fallback": the
picked fit runs first; when its end does not pass the acceptance test that the caller applies
next, the stock fit runs from the same start, and the end with the lower cost is returned.
Read that test in `reconstructor.py` and `grain_paint.py` and name it on the page. If it cannot
be evaluated inside a script wrapper, use "the cost did not fall" as the trigger and say so.

| run | what | compared with | cap |
|---|---|---|---|
| `D4` | BFS with indexed seeds on the gradient sample, wrapped fit | step 3's D (118,350,225 evaluations, 203 wrong) | 150 min |
| `P4` | paint through the key on the gradient sample, wrapped fit in the fit-settled branch; then `R4`, one picked fit from each painted start | step 3's P plus R (90,985,770 evaluations, 170 wrong) | 120 min |
| `E4` | `D4`'s protocol on step 2's sample (one orientation per grain, stored images) | step 2's D (33,040,126 evaluations, 153 wrong) | 90 min |

The three runs may run side by side. Each cap is coded: the background script kills the run at
its cap, and a killed run is reported as not finished. Where it is cheap, a run writes what it
has at intervals so that a killed run still leaves a record.

Readings of part 3: evaluations by stage; the number of grain entries and of grid fallbacks;
wrong voxels, the share above 0.1 degrees, median and 95th percentile of the error, all by
level; paired counts by voxel against the compared run (fixed, broken, wrong in both). A
difference is a gain or a loss only when the intervals do not overlap. Each arm is one run with
one voxel order, and the page says so next to each comparison.

**Stop rules of the task.** Part 1 and part 2 together: 60 minutes. Part 3: the caps above.
Total compute about 4.5 hours. The general stop rule holds as well (twice a budget, a baseline
that does not reproduce, a plan that cannot be followed).

**Not claimed.** That C++ behaves the same (the source has the same logic; the program is not
run here). That the gradient model matches a real sample. A change of any default. A tested
wall time unless the preflight passes.

**What comes after** (not part of `ar3d`). If the wrapped fit is a gain on the whole layer: an
opt-in key for it in the library, with the off path bit-identical. Then the stop rule inside
the loop, the grain-wide model of the orientation field, and the cheaper seed, as listed under
step 3.

**Practical notes for the implementer.** Step 3 is merged: reuse `run_gradient.py` (its refine
stage and its trace wrapper), `gradient_common.py` and `gradient_summary.py`. The untracked run
files of step 3 (P's `recon.mic` and `records.json`, the pairing cache, the image symlink) are
in the worktree named in `docs/autoresearch/STATE.md`; link them, do not commit them. Check
first that the trace wrapper does not change a fine fit: the first 200 fits identical with the
trace on and off.

**Files.** No file under `icenine_py/icenine/` changes. Scripts in `scripts/autoresearch/ar3/`;
results in `benchmarks/autoresearch/ar3/fine/`; findings in
`docs/autoresearch/findings_ar3_fine.md` (returned as text, see "Task split").
`MIGRATION_HISTORY.md` gets the summary. Only the result files that the summary script reads
are committed; a per-voxel record above 2 MB is stored gzipped and read gzipped.

### Step 5: let neighbouring grains compete for the boundary voxels

Added 2026-11-10 after the review of step 4. Task `ar3e`.

**Why.** Step 4 found a local fit that moves (the existing CMA fit, first step 0.05 degrees, 61
evaluations). In the whole-layer search it cut the evaluations from 118,350,225 to 35,129,496.
But two of the three whole-layer runs ended with more wrong voxels (281 against 170; 305
against 153). The review read the stored records. In the paint run, all 134 newly wrong voxels
were already wrong after the paint, 107 of them were settled by a fit, and 104 of those 107
fits started from a grain that was seeded in another true grain. All 107 passed the caller's
0.9 test. On step 2's sample, the 170 newly wrong voxels were accepted at the first neighbour
fit (114) or at a revisit (56). So the loss sits in fits that the caller accepts: a voxel next
to a grain boundary takes the orientation of the wrong side, the small fit polishes it, and
the 0.9 test lets it pass. The stock fit used to fail that test, and the voxel was then
settled another way. A fallback for rejected voxels does not reach an accepted fit.

**The idea.** A voxel that passed from one side of a boundary was never tried from the other
side. Try it. After a run, every voxel that has a neighbour with a clearly different
orientation gets one small fit from that neighbour's orientation, and keeps the better of the
two ends. This costs 61 evaluations per try and is run only at boundaries. It uses nothing but
the reconstructed map, and it does not depend on the crystal symmetry. The docstring of
`_expand` in `reconstructor.py` says that C++ refits a fitted voxel when another neighbour
reaches it and keeps the end with the higher confidence, and that the Python port does not.
The implementer reads `BreadthFirstReconstructor.tmpl.cpp` and reports what C++ does there.
Whether the pass removes the extra wrong voxels, and what it costs, are to be measured;
neither is claimed here.

**Terms.**

- Two voxels are **neighbours** when the reconstructor's own neighbour relation says so.
- A **boundary voxel** is a voxel with at least one neighbour whose orientation in the map
  differs from its own by more than 0.5 degrees (symmetry-reduced). The 0.5 is set, not tuned:
  inside a grain of this sample, neighbours differ by far less; across a boundary by far more.
- The **candidates** of a boundary voxel are its neighbours' orientations that differ from its
  own by more than 0.5 degrees, with candidates within 0.5 degrees of each other counted once.
- A **try** is one CMA fit (first step 0.05 degrees, at most 60 evaluations plus the final
  one) at the voxel, started from a candidate.
- The **pass** makes every try of every boundary voxel. The candidates are read from the map
  as it was before the pass, so the order of the voxels does not matter. A **keep rule** then
  decides, per voxel, between its old end and the ends of its tries.
- **Error**, **wrong**, **level**, pick set **A** and test set **B** are as in step 3.

**Part 0: can the pass work at all? (records and single scores; the truth is used for analysis
only.)** For every wrong voxel of six stored maps (step 4's `D4`, `P4` plus `R4`, `E4`; and
the runs they were compared with: step 3's D, step 3's P plus R, step 2's D):

1. Is the voxel a boundary voxel? Is there a candidate within 1 degree of its true orientation?
2. One strict score of the true orientation and of the voxel's end: which costs less, which
   has the higher hit ratio, which has the higher confidence.

Also, for all voxels of each map: the number of boundary voxels and of candidates, so that the
cost of the pass is known before it runs.

**Gate.** Take the newly wrong voxels of the three step 4 runs (134, 170, and the 108 of
`D4`). If, on fewer than half of them, a candidate lies within 1 degree of the truth and the
truth scores better than the end by at least one of the three measures, stop after part 0 and
report. The pass would then have little to offer.

**Part 1: the pass on the stored maps.** No whole-layer run. Run the pass once on each of the
six maps, with the images and Q-max of the run that made the map, one worker, one generator
state. The tries do not depend on the keep rule, so each map needs one set of tries; the rules
are applied to the stored tries.

| keep rule | the voxel takes the end with |
|---|---|
| `cost` | the lowest strict cost |
| `hit` | the highest hit ratio (a tie: the lowest cost) |
| `conf` | the highest confidence (the rule the docstring gives for C++) |

**The pick.** On A, over the three step 4 maps together: the rule with the fewest wrong voxels
after the pass. Among rules whose intervals overlap with that one, the rule that changes the
fewest voxels. The pick is reported on B.

Readings of part 1, per map, for the picked rule (the other two rules in one short table):

1. Boundary voxels, tries, evaluations. The evaluations of the pass as a share of the run's.
2. Voxels changed. Paired counts against the map before the pass: fixed, broken, wrong in both.
3. Wrong voxels, the share above 0.1 degrees, median and 95th percentile of the error, all
   voxels and by level, before and after.
4. The whole pipeline: run plus pass against the compared run, and against the compared run
   plus pass (the same pass on both sides is the fair comparison). Evaluations and wrong
   voxels, with paired counts. A gain or a loss only when the intervals do not overlap.
5. A second pass on the map the first pass made: how many voxels change again.
6. When is a try needed? From what is known before the pass (the voxel's hit ratio over the
   highest hit ratio among its neighbours; its cost), a rule that skips the tries of a voxel.
   Thresholds on a fixed grid, picked on A (the rule that skips the most tries while losing at
   most 1% of the voxels the pass fixed), reported on B. A pick at the edge of its grid is
   called uninformative.

**Part 2: the same saving with the keys that exist.** Step 4 used a script wrapper. The
library already has the keys `LocalOptimizer cma`, `CMASigma0`, `CMANeighborMaxEvals` and
`CMARetrySigma0` (a wider second fit when the caller rejects the first; this is close to the
fallback step 4's plan asked for). Two whole-layer runs with these keys and no wrapper:

| run | what | compared with | cap |
|---|---|---|---|
| `K` | BFS with indexed seeds on the gradient sample; neighbour and revisit fits are the CMA fit of step 4 (first step 0.05 degrees, 60 evaluations), set through the keys; then the pass with the picked rule | `D4` plus pass; step 3's D plus pass | 60 min |
| `KE` | the same on step 2's sample (one orientation per grain, stored images) | `E4` plus pass; step 2's D plus pass | 60 min |

Leave every other key as in the compared runs. List every key set, and say what the seed's own
fit does in this mode (read it in the code; do not assume it equals `D4`'s). The two runs may
run side by side. Each cap is coded, as in step 4. Readings: as part 3 of step 4 (evaluations
by stage, entries, grid fallbacks, error by level, paired counts), before and after the pass.
The paint run is not repeated here; whether the paint needs the pass is read from part 1.

**Stop rules of the task.** Part 0 and part 1 together: 90 minutes. Part 2: the caps above.
Total compute about 3.5 hours. The general stop rule holds as well.

**Not claimed.** That C++ behaves the same (the source is read; the program is not run). That
the gradient model matches a real sample. A change of any default. A tested wall time unless
the preflight passes. That a rule picked here holds on another sample.

**What comes after** (not part of `ar3e`). If the pass is a gain: an opt-in key for it in the
library (a new module, the off path bit-identical). Then the seed cost, which is now the
largest part: in `D4` the grid search behind 76 failed index seeds took 25,833,494 of the
35,129,496 evaluations. Then the stop rule inside the grid search (AR5), and the grain-wide
model of the orientation field.

**Practical notes for the implementer.** Step 4 is merged: reuse `run_fine.py` (how a script
calls the CMA fit and counts its evaluations), `wrapped_summary.py` (paired counts, tables by
level) and `ar3_common.py`. The stored maps are `recon.mic.gz` (or `recon.mic`) in each run's
directory; the untracked ones are in the worktrees named in `docs/autoresearch/STATE.md`
(link them, do not commit them). Check first that one try reproduces a fit of step 4: for 200
voxels of run `C`, a try from the painted start gives the stored end and evaluation count.

**Files.** No file under `icenine_py/icenine/` changes. Scripts in `scripts/autoresearch/ar3/`;
results in `benchmarks/autoresearch/ar3/boundary/`; findings in
`docs/autoresearch/findings_ar3_boundary.md` (returned as text if the write is refused, see
"Task split"). `MIGRATION_HISTORY.md` gets the summary. Only the result files that the summary
script reads are committed; a per-voxel record above 2 MB is stored gzipped and read gzipped.

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
