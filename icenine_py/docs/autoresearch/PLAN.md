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
`icenine/gradient_sample.py` (AR3 step 3, D-20261110-N4),
`icenine/boundary_pass.py` (AR3 step 10, D-20261111-N8)
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

### Step 6: the seed cost

Added 2026-11-10 after the review of step 5 (D-20261110-N7). Task `ar3f`.

In this step a **local optimization** is a refinement from a starting orientation, with no
grid search. Earlier steps of this plan call it "finisher" or "local fit"; it is the same
thing. A **fit** is one run of a local optimization.

**Why.** With the small CMA fit of step 4 and the boundary pass of step 5, the neighbour and
revisit fits are no longer the largest cost. In step 4's run `D4` (35,129,496 evaluations) the
grid search took 25,833,494. It ran 76 times, each time because the index seed failed the seed
test. Run `K` of step 5 had spent 31,067,562 of its 38,149,476 evaluations the same way when
its cap stopped it. The index seeds themselves took 3,389,254 evaluations in `D4` (1,298
seeds, about 2,600 each); that is more than the neighbour fits (2,817,382) or the revisit fits
(3,089,366). So the seed step is now where the operations go.

**What the stored records of `D4` already show** (read by the main session; the implementer
derives these again and the page quotes the implementer's numbers, not these):

- Each grid search costs about 340,000 evaluations, whatever the voxel.
- 74 of the 76 are at reseeds. 2 are at first seeds (2 of 344). So the index fails mostly
  on voxels that the BFS had already failed to settle.
- After the grid search the voxel's hit ratio is 0.88 or higher on 73 of the 76, and 1.0 on
  24. So an orientation that passes the seed test exists there, and the index seed did not
  reach it. The records do not say why.
- 36 of the 74 reseeds were for a cluster of one voxel and 11 for a cluster of two. A grid
  search of 340,000 evaluations then settles one or two voxels.
- Of all 954 reseeds, 260 were for a cluster of one voxel and 101 for a cluster of two. Each
  pays an index seed of about 2,600 evaluations. A try from a neighbour costs 61.
- Step 2's sample (`E4`) had no failed index seed. The failures belong to the gradient sample.

**The questions.** (1) Why does the index seed fail where the grid search succeeds? (2) Which
cheaper step reaches an end as good as the grid search's on those voxels? (3) Do small
clusters need a seed at all? Nothing is claimed before it is measured.

**Terms.**

- The **seed step** of a voxel is what `seed_by_indexing` does: index the voxel (the
  candidates), take the strict cost of the first 5, run a local optimization from the best 2,
  keep the lowest cost. The **seed test** is the caller's: the end must have a hit ratio of at
  least the configured threshold (0.80 here; read the value from the run's config).
- A **grid fallback** is the full grid search that the caller runs when the seed step returns
  no candidate or its end fails the seed test.
- A **fallback voxel** is a voxel where a stored run made a grid fallback.
- A **cluster** is the set of unsettled voxels that one reseed is made for.
- **Error**, **wrong** (error above 1 degree), **level**, pick set **A** and test set **B**
  are as in step 3. A voxel belongs to A or B by its true grain.

**Part 0: why the index seed fails (stored records and a replay; no whole-layer run; the
truth is used for analysis only).** The seed step does not depend on the state of the BFS: it
reads the voxel's observed rays and uses its own generator. So it can be run again at any
voxel.

1. List the fallback voxels of every stored run on the gradient sample that has per-voxel
   records (`D4`, `P4`, step 3's D and P). Say how each is identified. Give the counts by
   seed or reseed, by cluster size, by level, and by whether the voxel is at a true grain
   boundary. Give the evaluations by cluster size for all reseeds.
2. Check: the replayed seed step fails the seed test on these voxels, and passes on 200
   index seeds that the run used (a fixed draw, `default_rng(0)`). If more than 5% of either
   group disagree with the stored run, stop and report.
3. Replay the seed step on each fallback voxel with the indexer returning 50 candidates, and
   log: the number of candidates; the rank of the first candidate within 1 degree of the
   truth; the strict cost and the hit ratio of the first 20 candidates and of the truth;
   which candidates the seed step costed and which it refined; the end of each local
   optimization (error, cost, hit ratio, evaluations).
4. Put each fallback voxel in one class, in this order:
   (a) the truth itself fails the seed test;
   (b) no candidate within 1 degree of the truth;
   (c) such a candidate exists but beyond the first 5 (not costed);
   (d) costed but not among the 2 refined;
   (e) refined, and the local optimization ends away from the truth or below the seed test;
   (f) the replay passes (should be empty after check 2).
5. The grid search's end at each fallback voxel: error against the truth, hit ratio. Is the
   grid's end right? Count the settled neighbours (in the final map) whose orientation is
   within 1 degree of the voxel's truth: could a try from a neighbour have started near it?
6. The same log (items 3 to 5 without the grid) for the 200 control seeds, so that the cost
   of each remedy on ordinary seeds is known.

**Gate.** If class (a) holds more than half of the fallback voxels, no change of the seed
step can pass the test there. Then part 1 tests only remedies 4 and 5 and the page says that
the seed test, not the index, is the question.

**Part 1: remedies, by replay on the fallback voxels and the controls.** No whole-layer run.
Each remedy replaces or follows the seed step at one voxel. All use the images and Q-max of
the stored run, one worker, a fixed generator.

| remedy | what it does at the voxel |
|---|---|
| 1 `more` | strict cost on the first 20 candidates in place of 5; local optimization from the best 2, and from the best 5 |
| 2 `cma-seed` | the seed step of the library's `LocalOptimizer cma` mode, as run `KE` used it (read it in the code; list what it does) |
| 3 `wide` | when the seed step fails: one CMA fit from its best end with a first step of 1.5 degrees and at most 1,000 evaluations |
| 4 `neighbours` | before the seed step: one try (CMA, 0.05 degrees, 61 evaluations) from each distinct orientation among the voxel's settled neighbours in the final map; keep rule `hit` of step 5 |
| 5 `defer` | no seed step and no grid search for a cluster of one or two voxels: only the tries of remedy 4 |
| 6 `outward` | when the seed step fails: the existing outward search (the policy behind the `OutwardSearch` key) from its best end, if a script can call it without a library change; say so if it cannot |

Remedy 4 reads the final map, which the run did not have at that moment. Say so as a limit:
it shows what a try from a neighbour can reach, not what the run would have had.

Readings of part 1, per remedy, on the fallback voxels and on the controls:

1. Voxels that pass the seed test. Voxels that end right. Paired counts against the grid
   search's end at the same voxel (right in both, fixed, broken, wrong in both).
2. Evaluations per voxel (median, mean, largest), against the grid search's.
3. On the controls: does the remedy change an end that was right? What does it add in cost?
4. Small clusters (question 3): for every reseed of `D4` with a cluster of one or two voxels,
   does remedy 5 end within 0.1 degrees of the stored end, with a hit ratio that is no lower?
   The evaluations it would have saved.

**The pick.** On the fallback voxels of A: the remedy (or the chain of two remedies, cheapest
first) with the fewest voxels that still need the grid search, among those that break no
more right ends than they fix. A tie goes to the fewer evaluations. Reported on B. The counts
are small (about 40 fallback voxels per set from `D4`); give them as counts with intervals and
do not call a difference of a few voxels a gain.

**Part 2: whole-layer runs with the pick.** A script wrapper around the seed step, as step 4
wrapped the fit; no library change. Two runs on the gradient sample, side by side:

| run | what | compared with | cap |
|---|---|---|---|
| `S` | `D4`'s setting with the picked remedy in the seed step; then the boundary pass with rule `hit` | `D4` plus pass; step 3's D plus pass | 60 min |
| `K2` | run `K` of step 5 again (the library's cma keys and no fit wrapper) with the picked remedy; then the pass | `S`; `D4` plus pass | 60 min |

Start both on a machine with no other heavy job of this track. Each cap is coded, as in step
4. Readings: evaluations by stage (seed step, grid fallbacks, neighbour, revisit, pass),
seeds, reseeds, grid fallbacks, wrong voxels and error by level before and after the pass,
paired counts. The comparison that counts has the pass on both sides. If part 1 finds no
remedy that removes at least half of the grid fallbacks on A, part 2 is not run and the page
says so.

**Stop rules of the task.** Parts 0 and 1 together: 90 minutes of compute. Part 2: the caps
above. The general stop rule holds as well.

**Not claimed.** That C++ behaves the same. That the gradient model matches a real sample. A
change of any default. A tested wall time unless the preflight passes. That a remedy picked
on about 40 voxels holds on another sample.

**What comes after** (not part of `ar3f`). If a remedy removes the grid fallbacks: the same
check on a sample with larger grains, then opt-in keys for the remedy and for the boundary
pass. The learned proposer (AR4) is a candidate for the voxels that every remedy here leaves
to the grid search. Then the stop rule inside the grid search (AR5).

**Practical notes for the implementer.** Step 5 is merged: reuse `boundary_lib.py` (a try,
the keep rules, the neighbour relation), `boundary_common.py` and `run_keys.py`. `D4`'s
`run.json` has no per-seed log; a fallback voxel is found from its record (a seed record with
far more evaluations than an index seed) and from `stats.reseeds`. Check that against
`index_fallback_seed_test`. The stored maps and images are where step 5's page says.

**Files.** No file under `icenine_py/icenine/` changes. Scripts in `scripts/autoresearch/ar3/`;
results in `benchmarks/autoresearch/ar3/seed/`; findings in
`docs/autoresearch/findings_ar3_seed.md` (returned as text if the write is refused, see "Task
split"). `MIGRATION_HISTORY.md` gets the summary. The page leads with a plot, says "local
optimization", and defines each term at first use. Only the result files that the summary
script reads are committed, and they hold only the fields it reads: per-voxel logs as gzipped
JSON, no per-try arrays of unused fields. The committed results of this step stay under 5 MB.

### Step 7: the same pipelines on larger grains

Task `ar3g`. It starts after step 6 is merged. No file under `icenine_py/icenine/` changes.

**Why.** Step 6 left this split of run `S5` (24,570 voxels, with the boundary pass, 20.6
million evaluations): the seed step 6.0 million, 15 grid fallbacks 5.1 million, the local
optimizations 5.9 million, the pass 3.7 million. A seed and a grid fallback happen once per
reseed, so their cost follows the number of grains. The gradient sample has 497 grains, about
49 voxels each. A sample with larger grains has fewer reseeds per voxel. This step measures
what the pipelines cost there and whether step 6's seed step still does no harm. It is the
check that D-20261111-N1 puts before any opt-in key.

**Terms.** *Reseed*: the BFS has no neighbour to start from and runs the seed step. *Seed
step*: ray indexing proposes candidates and the best are refined. *Grid fallback*: the seed
step fails and the full grid search runs. *Boundary pass*: step 5's pass, rule `hit`.

**The sample `big`.** 32 grains on the grid of the gradient sample (24,570 voxels, the same
voxel side), about 770 voxels per grain.

- Orientations come from `Examples/Example2.ManyGrains/SimInput/random_cubic_32grains.mic`
  (6,144 voxels, one size coarser). Each voxel of `rand_500grains_1mm_neworient_s0_grid.mic`
  takes the orientation of the coarse voxel that holds its centre. Check: every fine voxel
  finds exactly one coarse voxel, and there are 32 distinct orientations.
- A grain map is written from the 32 orientations (connected voxels of one orientation).
- Gradients: `gradient_orientations` as in step 3, seed 0, one level per grain, with the
  levels 0, 0.2, 0.5, 1, 2 and 4 degrees of spread across a grain. A grain here is about 4
  times wider than in the gradient sample. So a spread of 0.4 to 4 degrees here gives the same
  change from voxel to voxel as 0.1 to 1 degree there, and 0.2 to 1 degree here gives a
  gentler one. If the function takes no list of levels, the script sets the levels; the
  library file is not edited.
- Images: clean, rendered with `scripts/phase_d/render_full.py --mic`, tag `big`, packed into
  one stack as in step 3. The stack is not committed.
- Check before any run: at the truth, 200 voxels drawn with seed 0 pass the seed test at the
  same share as on the gradient sample (step 3's check).

**Runs.** Four whole-layer runs on `big`, side by side, each with a coded cap and each
followed by the boundary pass. The sample and stack are chosen by an argument or an
environment variable in the scripts; the default stays the gradient sample.

| run | what | cap |
|---|---|---|
| `D` | the stock BFS with the grid search (step 3's setting D) | 180 min |
| `D4` | step 5's run: index seeds, CMA local optimization, grid fallback | 90 min |
| `S5` | `D4` with step 6's seed step (refine the best 5 of 20 candidates) | 90 min |
| `P4R4` | step 4's run: index the layer, paint grains, CMA local optimization | 90 min |

**Written down before the runs.** (1) Reseeds fall to the order of the number of grains, and
the seed step and the grid fallbacks together take less than a tenth of the evaluations of
`S5`. (2) The cost of `S5` with the pass comes close to its local optimizations and pass
alone, about 390 evaluations per voxel on the gradient sample. (3) At the levels with the same
voxel-to-voxel change, the share of wrong voxels is no higher than on the gradient sample.
The page says for each whether the run agrees.

**Readings.** Evaluations by stage (seed step, grid fallbacks, neighbour, revisit, pass);
seeds, reseeds, grid fallbacks; wrong voxels and error by level, before and after the pass,
with 95% intervals; paired counts (fixed, broken) for `S5` against `D4`, `S5` against `D` and
`P4R4` against `S5`, all with the pass. One table puts the two samples side by side for each
pipeline. The page leads with a plot: evaluations per voxel by stage, per pipeline, on the two
samples; then the share wrong. Wall times are labelled as contended. A level holds only about
five grains: say so, and give its result as a result on those grains.

**What the result decides.** If `S5` with the pass has no more wrong voxels than `D4` with
the pass (intervals) at fewer or equal evaluations, the opt-in keys step uses step 6's seed
step. If it has more, the keys step uses `D4`'s seed step and the page says so. If the seed
step and the grid fallbacks are a small part on `big`, the next target is what is left: the
local optimization per voxel and the pass. The split is reported for that.

**Stop rules of the task.** Sample and images: 30 minutes of compute. Runs: the caps above. A
run that hits its cap is reported with what it reached, not repeated. The general stop rule
holds as well.

**Not claimed.** That C++ behaves the same. That the gradient model or these levels match a
real sample. A change of any default. A tested wall time unless the preflight passes. That 32
grains stand for all large-grained samples: it is one sample, one run each.

**Files.** Scripts in `scripts/autoresearch/ar3/` (`big_*.py`, and a sample switch in the
common files); results in `benchmarks/autoresearch/ar3/big/`, slim, under 5 MB; findings in
`docs/autoresearch/findings_ar3_big.md` (returned as text if the write is refused).
`MIGRATION_HISTORY.md` gets the summary. The page says "local optimization" and defines each
term at first use. The stored maps and images of earlier steps are where step 6's page says.

### Step 8: why the seed step costs so much in grains with strong gradients

Task `ar3h`. It starts after step 7 is merged. No file under `icenine_py/icenine/` changes.
No key is added and no default changes.

**Why.** Step 7 left the seed step as the largest stage on the sample `big`: 7.47 of the 13.00
million evaluations of run `S5` with the boundary pass (57.5%), and, with the grid fallbacks,
6.47 of the 11.94 million of run `D4` (54%). Two numbers make this cost, and step 7 measured
neither cause.

- *How many seed steps.* `S5` ran 349 on 32 grains. 267 of them are in the 6 grains with a 4
  degree gradient and 54 in the 6 grains with 2 degrees. The grains without a gradient have 5
  seed steps and no reseed.
- *What one seed step costs.* About 21,400 evaluations in `S5` on `big` (7,469,574 over 349)
  against about 4,600 on the first gradient sample (5,990,794 over 1,293). `D4`'s seed step
  shows the same factor: 12,208 against 2,611.

If a seed step on `big` cost what it costs on the first gradient sample, `S5` would need about
7 million evaluations there, not 13. This step finds out where the evaluations go, tries
remedies on stored cases first, and only then spends whole-layer runs.

**Terms.** *BFS*: the breadth-first search that spreads from voxel to voxel; each voxel gets a
local optimization started from a fitted neighbour's orientation. *Relative test*: the BFS
accepts that fit if its hit ratio is above 0.9 times the best hit ratio seen so far in the
same flood. *Reseed*: after the BFS, the rejected voxels form connected clusters, and the seed
step runs on the voxel nearest the centre of a cluster (`_reseed_clusters`). *Seed step*: ray
indexing proposes candidate orientations, the candidates are costed, and the best are refined
by a local optimization (`D4`: 5 costed, 2 refined; `S5`: 20 costed, 5 refined). *Seed test*:
a hit ratio of at least 0.80. *Grid fallback*: the seed step fails the seed test and the full
grid search runs (a global search: no starting guess). *Boundary pass*: step 5's pass, rule
`hit`. *Right*: within the angle from the truth that step 7's page uses for a right voxel.

**Part 1: where the evaluations go.** From the stored runs of step 7 (`S5` and `D4` on `big`)
and of step 6 (the same two on the first gradient sample). If a stored log lacks a field, one
logged rerun of `S5` on `big` is allowed (cap 90 minutes, the same voxel order); the page
then shows that its counts match the stored run.

1. *One seed step, split.* For every seed step: the evaluations spent costing candidates, and
   the evaluations of each refinement in the order they ran; which refined candidate was
   kept; the hit ratio of each refined candidate; whether a candidate refined earlier already
   passed the seed test at a right orientation. Reported by gradient level and by sample.
2. *Why the BFS rejected.* For every voxel in a reseed cluster, the rejected fit from a
   neighbour falls in one class, judged with the truth:
   - A, reached and rejected: the fit ended right, and the relative test rejected it;
   - B, not reached: the start came from the voxel's own true grain, and the fit did not end
     right (split by the start's angle from the truth: under 0.1, 0.1 to 0.2, 0.2 to 0.5,
     over 0.5 degrees);
   - C, other grain: the start came from another true grain.
   Reported with the hit ratio, the flood's best hit ratio at that moment, and the voxel's
   distance to the nearest true grain boundary, by gradient level.
3. *The clusters.* Size of each reseed cluster, its level, its distance to the boundary, and
   how many voxels the flood from its reseed then settles.

**Part 2: remedies, on stored cases first.** The seed steps and clusters are split by true
grain: grains with an even index are the pick set, grains with an odd index the test set. A
remedy's settings are chosen on the pick set and written into the results file before the
test set is read. The list is fixed here; a remedy is tried only if Part 1 shows its class or
its cost share is at least a fifth of its total.

| remedy | aimed at | what it does |
|---|---|---|
| `E1` | one seed step | refine in order of cost; stop at the first that passes the seed test |
| `E2` | class A | keep a rejected fit that passes a second test, chosen on the pick set |
| `E3` | class B | try a rejected fit again from a better start |

For `E2` the second test is, for example, the seed test, or 0.9 times the hit ratio of the
neighbour the fit started from. For `E3` the better start is a wider first step, or a start
extrapolated from two fitted neighbours in a line.

Step 6 showed the risk of `E2` and `E3`: tries from the neighbours settled 602 voxels and 59
of them were wrong, because a boundary voxel accepted the orientation of the grain next to it.
So each remedy is reported with the count of cases it settles wrong, paired with what the
stored run did on the same cases (fixed, broken). A remedy goes on to Part 3 if, on the test
set, it cuts the evaluations of its target by at least a third and breaks no more cases than
it fixes.

**Part 3: whole-layer runs.** Only if a remedy passes. Run `S8` is `S5` with every remedy
that passed, on `big` and on the first gradient sample, each followed by the boundary pass
(cap 90 minutes each). If more than one remedy passed, one more run on `big` with only the
remedy that saved most in Part 2. No other whole-layer run. Same accounting as step 7:
evaluations by stage, seeds, reseeds, grid fallbacks, wrong voxels before and after the pass
with 95% intervals, paired counts against `S5` and `D4` with the pass.

**Written down before the runs.** (1) Most rejected fits in the clusters of the grains at 2
and 4 degrees are of class A or B, not C. (2) A seed step is dearer on `big` because of the
refinements, not because of the costing of candidates. (3) `S8` with the pass needs under 9.0
million evaluations on `big` (`D4`: 11.94) and no more than `S5`'s 20.6 million on the first
gradient sample, with intervals of wrong voxels that overlap `D4`'s on `big` and `S5`'s on
the first gradient sample. The page says for each whether the result agrees.

**What the result decides.** If prediction 3 holds, `S8`'s seed step is my recommendation for
the opt-in key, and the owner's choice between the old and the new seed step gets simpler. If
no remedy passes Part 2, the page says which class or cost share is left and why the listed
remedies do not reach it; the next step is then AR4 or AR5 on those cases. The choice of the
seed step for the key stays the owner's in both cases.

**Rules of the study.** One run per arm, the stored voxel order; no extra random seeds. No
p-values: counts, paired fixed and broken, 95% intervals. Intervals treat voxels as
independent although wrong voxels cluster in grains; the page says so. A gradient level holds
1 to 9 grains: each result by level is a result on those grains. Nothing here is specific to
cubic symmetry or to special misorientations.

**Stop rules of the task.** Part 1: 60 minutes of compute besides the one allowed rerun.
Part 2: 60 minutes. Part 3: the caps above. A run that hits its cap is reported with what it
reached, not repeated. The general stop rule holds as well.

**Not claimed.** That C++ behaves the same. That the gradient model or these levels match a
real sample. A change of any default. A tested wall time unless the preflight passes. That
two samples stand for all samples.

**Files.** Scripts in `scripts/autoresearch/ar3/` (`reseed_*.py`); results in
`benchmarks/autoresearch/ar3/reseed/`, slim, under 5 MB; findings in
`docs/autoresearch/findings_ar3_reseed.md` (returned as text if the write is refused).
`MIGRATION_HISTORY.md` gets the summary. The page says "local optimization", defines each
term at first use, and leads with a plot: evaluations of one seed step by part, per level, on
the two samples; then the classes of rejected fits per level. The stored maps, logs and image
stacks of steps 6 and 7 are in the worktrees that `STATE.md` lists; they are linked, not
copied into git.

### Step 9: a start test for the neighbour fits, the rejected fits that end right, a third sample

Task `ar3i`. It starts after step 8 is merged. No file under `icenine_py/icenine/` changes.
No key is added and no default changes.

**Why.** Step 8 leaves three things open.

- *Run `S8` is not confirmed.* Its stop test (hit ratio 0.95) was added after results
  (D-20261111-N5). It has to be an arm of a whole-layer run on a sample it was not chosen on.
- *The neighbour fits are now the largest stage on `big`.* In `S8` with the boundary pass they
  are 3.49 of 8.34 million evaluations, and the revisits 0.85 million more. A fit costs up to
  61 evaluations. Nothing stops it when its start is already good.
- *Class A.* 1,069 of the 1,471 voxels in reseed clusters of `S5` on `big` had a fit from a
  neighbour that ended right and was rejected by the relative test. The 292 reseeds of `S8`
  cost most of its 2.86 million seed-step evaluations. Step 8 named two explanations and
  tested neither: the fit stops a few tenths of a degree short of the best orientation, or the
  true orientation itself scores under 0.9 times the flood's best.

**Terms.** As in step 8 (BFS, relative test, reseed, seed step, seed test, grid fallback,
boundary pass, classes A, B and C). *Neighbour fit*: the local optimization of a voxel started
from a fitted neighbour's orientation. *Revisit*: a later fit of a voxel the BFS marked for
another try. *Start test*: one evaluation of the start before a fit; if its hit ratio is at
least a threshold, the start is kept and the fit is skipped. *Right*: within 1 degree of the
truth, as on step 8's page. *Error*: the angle between a voxel's orientation and the truth.

**Part 1: what the stored runs say.** From step 8's logged run of `S5` on `big` and its
stored cases on the first gradient sample. One more logged rerun is allowed if a field is
missing (cap 90 minutes, same voxel order; the page shows its counts match).

1. *The start of every neighbour fit and revisit.* The hit ratio and the error at the start
   and at the end, and the evaluations used. Reported by gradient level: the share of fits
   whose start already has a hit ratio of at least 0.90, 0.95 and 0.98, and how far those
   fits moved.
2. *Class A, the check.* For every class A voxel: the hit ratio at the true orientation,
   against 0.9 times the flood's best at that moment. Then the same after one more local
   optimization from the rejected fit's end (first step 0.1 degrees, at most 60 evaluations).
   The page says which explanation holds for how many voxels, and whether the stored rejected
   fit already includes the BFS's own retry.

**Part 2: remedies, on stored cases first.** The split is step 8's: true grains with an even
index are the pick set, with an odd index the test set. The remedies and their settings are
fixed here. No remedy, setting or test that is not in this table is added after a result. The
pick is written to a results file and committed before the script that reads the test set is
run, so that git shows the order.

| remedy | aimed at | what it does | settings |
|---|---|---|---|
| `F1` | neighbour fits, revisits | the start test: keep the start, skip the fit | threshold 0.90, 0.95, 0.98 |
| `F2` | class A | before a fit is rejected, one more fit from its end | first step 0.05, 0.1, 0.2 degrees |

`F2` uses at most 60 evaluations and is then judged by the relative test. It runs on every
rejected fit, because the BFS does not know the class; its cost on classes B and C counts. It
is tried only if Part 1 shows that the fit stopped short in at least half of class A.

*The pick.* Among a remedy's settings, the one that saves most evaluations on the pick set
among those that break no more cases than they fix. If none qualifies, the remedy fails. *The
pass rule.* On the test set the remedy cuts the evaluations of its target by at least a third
and breaks no more cases than it fixes. For `F1` the target is the neighbour fits and
revisits; a case is broken if the stored fit ended right and the kept start is wrong. For
`F2` the target is the seed steps of the reseed clusters that vanish, less `F2`'s own cost.

Stored cases cannot show one risk of `F1`: a kept start is passed on to the next voxel, so
the error can add up across a grain with a gradient. Only the whole-layer run measures that.
The page therefore reports the error (median, 90th and 99th percentile), not only the count
of wrong voxels.

**Part 3: whole-layer runs.** A third sample, `mid`: step 7's generator with its stored
settings and seed value, about 120 grains on the same 24,570-voxel grid (about 205 voxels per
grain, between the two samples so far), the same six gradient levels in step 7's proportions,
clean images. No setting was chosen on it.

- On `mid`: runs `D4`, `S5` and `S8`, each with the boundary pass. This is the confirmation
  that D-20261111-N5 asks for. It runs whatever Part 2 gives.
- Run `S9` is `S8` with every remedy that passed: on `mid`, on `big` and on the first
  gradient sample, each with the boundary pass. Only if a remedy passed.

Cap 90 minutes per run, 60 minutes for the image stack of `mid`. No other whole-layer run.
Accounting as in step 8: evaluations by stage, seeds, reseeds, grid fallbacks, wrong voxels
before and after the pass with 95% intervals, the error percentiles, paired counts.

**Written down before the runs.** (1) In most class A voxels the truth passes the relative
test: the fit stopped short. (2) In `S5` on `big`, more than half of the neighbour fits in the
grains at 1 degree or less start with a hit ratio of at least 0.95, and fewer than a quarter
in the grains at 4 degrees. (3) On `mid`, `S8` with the pass needs fewer evaluations than
`S5` and than `D4` with the pass, and its interval of wrong voxels overlaps that of the one
with fewer wrong voxels. (4) `S9` with the pass needs under 6.5 million evaluations on `big`
(`S8`: 8.34), with an interval of wrong voxels that overlaps `S8`'s and a median error no
more than 0.05 degrees above `S8`'s. The page says for each whether the result agrees.

**What the result decides.** If prediction 3 holds, `S8`'s seed step counts as confirmed on a
sample it was not chosen on, and it is my recommendation for the opt-in key. If it does not
hold, the page says so and the recommendation falls back to "both numbers as settings, no
stop test". If `F1` passes and prediction 4 holds, the start test is a candidate for a second
opt-in key. The choice of what a key recommends stays the owner's in every case. AR5 proper
(a stop rule inside the grid search, with the vote count of unexplained rays) is not part of
this step: it aims at the grid fallbacks, 5.45 of `S8`'s 18.00 million evaluations on the
first gradient sample.

**Rules of the study.** One run per arm; no extra random seeds (`mid` is a new grain size,
not a repeat). No p-values: counts, paired fixed and broken, 95% intervals. Intervals treat
voxels, fits and seed steps as independent although they cluster in grains; the page says
so. A gradient level holds few grains: each result by level is a result on those grains.
Every number chosen after a result is listed as such. Nothing here is specific to cubic
symmetry or to special misorientations.

**Stop rules of the task.** Part 1: 60 minutes of compute besides the one allowed rerun.
Part 2: 60 minutes. Part 3: the caps above. A run that hits its cap is reported with what it
reached, not repeated. The general stop rule holds as well.

**Not claimed.** That C++ behaves the same. That the gradient model or these levels match a
real sample. A change of any default. A tested wall time unless the preflight passes. That
three samples stand for all samples.

**Files.** Scripts in `scripts/autoresearch/ar3/` (`start_*.py`); results in
`benchmarks/autoresearch/ar3/start/`, slim, under 5 MB; findings in
`docs/autoresearch/findings_ar3_start.md` (returned as text if the write is refused).
`MIGRATION_HISTORY.md` gets the summary. The page says "local optimization", defines each
term at first use, and leads with a plot: evaluations by stage for `D4`, `S5`, `S8` and `S9`
on the three samples; then the share of fits whose start passes, per level. The stored maps,
logs and image stacks of steps 6 to 8 are in the worktrees that `STATE.md` lists; they are
linked, not copied into git.

### Step 10: opt-in keys for the indexed seed step and the boundary pass

Task `ar3j`. It starts after step 9 and its fixes round are merged. This is the first AR3 step
that changes files under `icenine_py/icenine/`. It adds opt-in keys only. No default changes.

**Goal.** Steps 6 to 9 ran their seed steps and the boundary pass as script wrappers around the
library BFS. A user cannot switch them on. This step moves that code into the library behind
opt-in keys, so that a config file alone gives runs `D4`, `S5` and `S8`, each with the pass.
With every new key absent, the library does exactly what it does today. The step recommends no
setting. It documents what each key does and what was measured; the choice stays the owner's.

**Decided with the plan (D-20261111-N8; the owner can overrule each).**

- A third key, `BFSNeighborCMAFirst`, is part of the step. Without it no stored run can be
  reproduced, and a user who switched on the seed keys alone would get a pipeline nobody ran.
- `icenine/boundary_pass.py` is added to the list of new modules in "How every arm is judged".
  `ray_indexing.py`, a module of this track, gets new arguments whose defaults are today's code.
- The pass is "exact on the stored map, a prediction in the run" (see the validation runs). The
  library does not round its map to imitate the script.
- The pass draws from its own `default_rng(0)`, as `seed_by_indexing` does. With the run's
  generator no stored pass could be reproduced.
- The pass uses every voxel that has a fitted orientation, as a target and as a source of
  candidates. A voxel that was never fitted is skipped in both roles. In the stored runs every
  voxel was fitted, so this changes no stored number.
- The fixes round of this step may change the library files the step itself touched, with the
  tests of this step and the bit-identity checks (a) and (b) run again. This is an exception
  to the runbook's rule for fixes rounds, for this step only.
- The implementer checks two things before writing code and reports them: whether any test or
  golden compares a whole `CMAResult` or `SeedResult` object (the new fields must not break
  it), and whether `cma_sigma0_deg` or `cma_neighbor_max_evals` is read on any path of
  `LocalOptimizer mc` (the grid search, `refine_from_candidates`). If either is so, the step
  stops and reports before the validation runs.
- If a stored map, record or `before_S8.npy` is missing from the worktrees that `STATE.md`
  lists, the exact check that needs it cannot run. The step says so; it does not rebuild it.

**Terms.** *BFS*: `BFSReconstruction` in `reconstructor.py`; it spreads from voxel to voxel.
*Seed step*: ray indexing proposes candidates, the first ones get the strict cost, and the
cheapest are refined by a local optimization (`ray_indexing.seed_by_indexing`, called through
`AdaptiveVoxelReconstructor.index_seed` and `BFSReconstruction._index_seed`). *Seed test*: a hit
ratio of at least `MinAccelerationThreshold` (0.80 here). *Early stop*: the seed step stops at
the first refinement whose end has a hit ratio of at least a threshold (`S8`: 0.95).
*Neighbour fit*: the local optimization of a voxel from a fitted neighbour's orientation
(`_local_fit`); a *revisit* is a later fit of a rejected voxel. *Boundary pass*: step 5's pass,
rule `hit`. *Off path*: a run with every new key absent or written at its default.

**What the scripts do today, and what the library lacks** (read in the code):

| run | script code | library code it calls | what the script adds |
|---|---|---|---|
| all | `run_wrapped.Wrap` | `local_optimization` | CMA first, MC if no gain |
| `S5` | `seed_pick.SeedPick` | `index_seed(..., 20, 5)` | 20 and 5, not 5 and 2 |
| `S8` | `reseed_s8.E1Pick` | `index_voxel`, `refine_from_candidates` | the early stop loop |
| pass | `boundary_pass.py`, `boundary_*` | `local_optimization` | all of the pass |

Reading the table: one row per piece of script code. The first row matters most. `D4`, `S5` and
`S8` all ran with step 4's fit wrapper. It replaces the instance method `local_optimization`:
a CMA fit (first step 0.05 degrees, at most 60 evaluations) runs first; if its end cost is not
below the start cost, the stock MC fit runs from the same start, and the lower cost is kept
(the CMA end on a tie). The library has no such mode. Its `LocalOptimizer cma` with
`CMARetrySigma0` (run `K`) is different: it also changes the seed path and retries on rejection.
So the seed keys alone cannot reproduce any stored run. This step therefore needs a third key.

**The keys.** Parsed in `config_file.py` (the parser table in `ConfigFile`, beside
`SeedByIndexing`), carried by `SearchParameters` in `orientation_search.py` (`from_config` and
`__post_init__` checks), and checked together in a new `check_bfs_keys(params)` in
`reconstructor.py`, called from `BFSReconstruction.__init__` beside `check_paint_grains_keys`.

| key (new) | type, default | field (new) | where it acts |
|---|---|---|---|
| `BFSNeighborCMAFirst` | 0 or 1, 0 | `bfs_neighbor_cma_first` | `_local_fit` |
| `SeedIndexCandidates` | two ints, `5 2` | `seed_index_n_cost`, `_n_refine` | `_index_seed` |
| `SeedIndexStopHitRatio` | float, 0 (off) | `seed_index_stop_hit` | `seed_by_indexing` |
| `BoundaryPass` | 0 or 1, 0 | `boundary_pass` | `reconstruct_sample` |

Reading the table: one row per key; the defaults equal today's behaviour. The second field is
`seed_index_n_refine`. `_local_fit`, `_index_seed` and `reconstruct_sample` are methods of
`BFSReconstruction`; `seed_by_indexing` is in `ray_indexing.py`.

- `BFSNeighborCMAFirst 1`: each neighbour and revisit fit is the wrapper's rule. The CMA fit
  uses `CMASigma0` and `CMANeighborMaxEvals` (the runs set 0.05 and 60). Needs `LocalOptimizer
  mc`. The seed path and the grid search stay as they are. `CMARetrySigma0` has no effect here.
- `SeedIndexCandidates n_cost n_refine`: the seed step costs `n_cost` candidates and refines the
  `n_refine` cheapest. `S5` and `S8` use `20 5`. Valid: 1 <= `n_refine` <= `n_cost` <=
  `ray_indexing.SEED_TOP_K` (20, the indexer's list length).
- `SeedIndexStopHitRatio t`: with 0 < t <= 1, the early stop at hit ratio t. If no refinement
  reaches t, all are refined and the lowest cost is kept, as today. The hit ratio of each end is
  one overlap evaluation, not counted in the totals (as in `S8` and in today's seed test).
- `BoundaryPass 1`: after the reseeds (and the restart pass, if on), before the `.mic` is
  written, the pass runs once on the whole map with rule `hit`. Its settings are constants of
  the new module, not keys: the 0.5 degree threshold, CMA tries of 0.05 degrees and 60
  evaluations, its own generator `default_rng(0)`, voxels in index order, candidates in
  neighbour-index order. These are the measured settings.

Invalid combinations raise `ValueError` before any search: either seed key off its default
without `SeedByIndexing 1`; any new key with `PaintGrains 1` (the paint path is not part of
this step); `BFSNeighborCMAFirst 1` with `LocalOptimizer cma`, with `OutwardSearch 1`, or with
the field `use_hybrid_optimizer` (no config key) set; these were never measured together. Bad
values raise at parse time: `SeedIndexCandidates 30 5`, `5 6` or one number;
`SeedIndexStopHitRatio 1.5` or negative; `BoundaryPass 2`.

**The design.** Code moves; nothing that reads the truth, logs to files, caches runs, counts by
monkeypatching or watches a time cap moves (`Env.err`, `G.Truth`, `CostRecorder`, progress
threads, watchdogs, `seed_log`, the pick and test sets stay in the scripts).

1. `orientation_search.py`: `CMAResult` gets `start_cost` (new, `Optional[float] = None`), set
   in `CMAOptimizer.optimize` from the start's own evaluation. This replaces the wrapper's
   `CostRecorder` without an extra evaluation.
2. `reconstructor.py`: `local_optimization` gets `optimizer: Optional[str] = None` (None: as
   today, from `params.local_optimizer`). `_local_fit`, when `bfs_neighbor_cma_first`, calls it
   with `optimizer="cma"` and the CMA overrides, then with `optimizer="mc"` if needed, and sets
   `last_local_optimization_evals` to the sum, as `Wrap.wrapped` does. Counters
   `cma_first_fallback_runs` and `cma_first_fallback_won` are made by `_count` (key on only).
3. `ray_indexing.seed_by_indexing` gets `stop_hit: float = 0.0` and `hit_fn=None`; `SeedResult`
   gets `stopped_at: int = -1` (new). With `stop_hit` 0 the loop is the current code.
   `AdaptiveVoxelReconstructor.index_seed` builds `hit_fn` from `evaluate_overlap` and
   `BFSReconstruction._ratios`. `_index_seed` passes the two counts and the threshold, and adds
   `stopped_at` to `index_log` and `index_stopped_early` to the stats only when the key is on.
4. New module `icenine/boundary_pass.py` (the scope list gets this name):
   `boundary_candidates(R, neighbours, sym)`, `keep_hit(old, tries)` and
   `run_boundary_pass(bfs, mic)`. It ports `Env.boundary`, `Env._merge`, `try_fit` and
   `apply_rule(..., "hit")`; neighbours from `MicFile.get_neighbors` (radius 2 side lengths,
   sorted), angles from `ray_indexing.misorientation_deg` and `proper_rotations`. Only boundary
   voxels get the strict score of their old end (the script scored all voxels; the score was
   never counted). A changed voxel takes the end's orientation, cost, confidence and hit ratio.
   Stats (key on only): boundary voxels, tries, voxels changed, `n_evals_boundary_pass`,
   `n_evals_boundary_score`, wall time; `bfs.boundary_log` keeps each changed voxel's old and
   new orientation, so the map before the pass can be rebuilt.

**How the off path stays bit-identical.** Every new branch is behind a new field at its
default, and the default calls are the current calls (`index_seed(vertices, phase, 5, 2)` is
what `_index_seed` makes today). Checks: (a) the whole suite and every golden unchanged, no
golden file touched; (b) the fingerprint tests below; (c) one validation run, outside the suite:
the BFS with `scripts/autoresearch/ar2/key_on_fused.config` on the 2,000-voxel pilot region, the
new keys absent, once on the parent commit and once on the task branch; orientations, records,
stats (wall times left out) and the written `.mic` bytes must be equal.

**Tests first** (written and committed before the implementation; each runs in seconds; new
file `tests/test_bfs_ar3_keys.py`, reusing `test_bfs_refit.FakeReconstructor`, `_run` and the
`_fingerprint` of `test_seed_by_indexing.py`):

1. Defaults and parsing of the four keys; `SearchParameters()` holds today's values.
2. Bad values and the invalid combinations above raise `ValueError`.
3. Keys absent equals keys written at their defaults, with `SeedByIndexing` and `BFSReseed` on
   and off: same orientations, records, stats; no new stats key; no new `index_log` field.
4. `test_key_off_scripted_grid_still_matches_the_golden`, again with the new keys written at
   their defaults (`bfs_grid_golden.npz`).
5. `_index_seed` passes `n_cost`, `n_refine` and the threshold to `index_seed` (a fake records
   the call).
6. The early stop on a fake reconstructor with scripted refinement costs and hit ratios: it
   stops at the first end at or above the threshold and keeps it although a later end would cost
   less; with none above, the lowest cost as today; the counted evaluations are the costing plus
   the refinements made; with threshold 0 the `SeedResult` equals today's field by field.
7. On the ThreeVoxels toy (real images, `needs_c`): `index_seed` with `20 5` and 0.95 equals
   the loop of `reseed_s8.E1Pick.e1_step` written out by hand, as
   `test_library_path_equals_the_script_path_on_the_toy` does for the old seed step.
8. `BFSNeighborCMAFirst` on the toy: from a start where the CMA fit lowers the cost, the CMA end
   and its evaluations; from one where it does not (the truth itself), the MC fit also runs, the
   lower cost is kept, the evaluations add up; equal to `Wrap.wrapped` written out by hand.
   `refine_from_candidates` makes no CMA call with the key on.
9. `start_cost` equals the first evaluation of a CMA run, and CMA results are otherwise
   unchanged.
10. `boundary_candidates` on a hand-made 6 x 6 grid with two orientations: the boundary voxels,
    the merge of candidates within 0.5 degrees (the first in index order kept), no truth used.
11. `keep_hit`: highest hit ratio, then lowest cost, then order; the old end wins a full tie.
12. The pass in the scripted BFS: it runs after the reseeds and before the write; candidates
    come from the map before the pass; with the key off it is never called.
13. On the toy with one voxel turned 5 degrees: `run_boundary_pass` equals `Env.boundary`,
    `try_fit` and `apply_rule(..., "hit")` written out by hand (tries, ends, evaluations).

**Validation runs** (not in the suite; scripts `scripts/autoresearch/ar3/libkeys_*.py`, new;
no truth inside the library; the scripts score against the truth afterwards). All on `big`,
clean images, Q-max 8, `key_on_fused.config` plus the keys, `default_rng(0)` for the voxel
order, one worker each. The configs:

| run | keys added to `key_on_fused.config` |
|---|---|
| `D4L` | `BFSNeighborCMAFirst 1`, `CMASigma0 0.05`, `CMANeighborMaxEvals 60`, `BoundaryPass 1` |
| `S5L` | as `D4L`, plus `SeedIndexCandidates 20 5` |
| `S8L` | as `S5L`, plus `SeedIndexStopHitRatio 0.95` |

Reading the table: one row per library run; each should give the stored run of the same name.

**Written down before the runs.** Exact (every number equal):

| run | seeds | reseeds | grid | evals before the pass | wrong before | source |
|---|---|---|---|---|---|---|
| `D4L` | 38 | 291 | 8 | 10,811,715 | 68 | `big/slim/D4_run.json` |
| `S5L` | 35 | 314 | 0 | 11,867,509 | 71 | `big/slim/S5_run.json` |
| `S8L` | 37 | 292 | 0 | 7,201,863 | 67 | `reseed/slim/s8_big.json` |

Reading the table: counts of the BFS before the pass. Paths are under
`benchmarks/autoresearch/ar3/`; the wrong counts are in `reseed/tables.md` (`s8_wrong`). Also
exact: every stats counter except wall times; for `S8L`, 322 early stops in 329 seed steps and
2,856,712 seed-step evaluations; for `S5L`, 3,235 MC fallbacks in 32,673 neighbour fits (the
wrapper counts). The map before the pass, rebuilt from `boundary_log`, written as a `.mic`,
equals the stored `recon.mic.gz` byte for byte (stored in the worktrees that `STATE.md` lists).

The pass is checked twice. (1) Exact: `run_boundary_pass` applied to the stored map that the
script read (`before_S8.npy`, `reseed/cache/passes_big/` in the step 8 worktree) gives the stored
`tries.npz` (18,665 tries, 1,138,565 evaluations) and the same map after rule `hit`. (2) In the
run: the library pass reads the map in memory, but the script read the `.mic` file, whose angles
are rounded (4 and 6 decimals). So the tries can start a tiny angle apart. Prediction, not
guaranteed: the same tries and evaluations as stored (`D4L` 18,557 and 1,131,977; `S5L` 18,563
and 1,132,343; `S8L` 18,665 and 1,138,565, from `big/passes/*_pass.json` and `s8_big.json`), and
the same wrong counts after the pass (21, 30 and 24). A difference here is reported with the
voxels it touches. It is not tuned away.

If an exact number differs, the step stops and reports the first differing voxel and counter.
No setting is changed to make a number match.

**Time budgets.** Stored wall times, all contended: `D4` 970 s, `S5` 1,027 s, `S8` 721 s; each
pass about 122 s of tries (4 passes side by side). Budget 30 minutes per whole-layer run with its
pass, 10 minutes for the pilot pair of check (c), 10 minutes for pass check (1), 15 minutes for
scoring. Total about 2 hours of compute; the three runs go side by side (3 heavy processes; at
most 6). The general stop rule holds: a run past twice its budget is stopped and reported.

**Optional, separable: a wall-clock pair.** Every AR3 run so far was contended, so no stored
wall time is a timing. If, and only if, `preflight.require_quiet()` passes: `D4L` and `S8L`
without the pass, one after the other, one worker, `FusedCost` on, preflight JSON saved beside
them, labelled "single-worker". If the preflight refuses, this part is not run here; it goes to
"Machine-time runs from the review page". Nothing else in the step depends on it.

**What the step must not do.** No default flip. No golden change. No new dependency. No
recommendation of a seed step. No key for the start test (`F1`, run `S9`). Nothing learned. No
change to C++. No grid-search change (the 16 grid fallbacks of `grad` stay). No new run on
`grad` or `mid`.

**Documentation.** `docs/autoresearch/findings_ar3_keys.md` (returned as text if the write is
refused): each key, what it does, the measured numbers of `D4`, `S5` and `S8` with the pass on
the three samples (`big`: 11.94, 13.00, 8.34 million evaluations and 21, 30, 24 wrong; first
gradient sample: 38.81, 20.63, 18.00 million and 74, 62, 63; `mid`: `D4` 40.56 million and 48,
`S8` 26.44 million and 31, `S5` from step 9's page), and the validation tables. It leads with a
plot: evaluations by stage, script run against library run. `icenine_py/README.md`: the four
keys in "Key Config Parameters for Reconstruction" and the new module in "Package Modules".
`MIGRATION_HISTORY.md`: a section for step 10. Results in `benchmarks/autoresearch/ar3/keys/`,
slim, under 5 MB tracked; maps and records are linked from the worktrees, not committed.

### AR3 step 11: a held-out sample with its `S8` baseline

Task `ar3k`. Scripts and results only: no file under `icenine_py/icenine/` changes. It starts
after AR5 step 2 is merged, because the baseline uses that step's script layer. Planned on
2026-11-11 (D-20261111-N15).

**Why.** Three samples exist (`big`, the first gradient sample, `mid`) and every one has now
been read for a choice. The mixed seed refinement of AR5 step 3 was found in a replay that
included `mid`. The cut of the learned verifier of AR4 step 1 was set on the other two, and
its use for accepting a seed was then read on `mid`. A rule picked on all the samples can fit
them. The next whole-layer tests therefore need one sample that no choice has seen. This is a
test of whether a result carries over to a new sample. It is not a repeat with a new random
seed to make a comparison look firmer.

**The sample `held`.** Built like `mid` (`scripts/autoresearch/ar3/start_make_mid.py`), with
three differences: 250 Voronoi cells instead of 120, generator seed 1 instead of 0 for the
centres, the orientations and the gradient levels, and the tag `held`. Same 24,570 voxel
centres, same six gradient levels (0, 0.2, 0.5, 1, 2 and 4 degrees), same rule that redraws an
orientation closer than 5 degrees to a neighbour grain's. Images: clean, rendered and packed as
for `mid` (`scripts/phase_d/render_full.py`, the `clean` variant). The step records the grain
sizes, the grains per level, the measured spread per level, the lit pixels and the spots, as
`mid_check.json` does for `mid`.

**The baseline.** Run `S8` on `held` with the boundary pass, with the script layer of AR5
step 2 and both of its changes off (the layer that reproduced stored `S8` on `big` exactly).
First repeat that check on `big` (7,201,863 evaluations, 67 wrong before the pass); if it
differs, stop and report. At most 120 minutes for the run on `held`; a run at its cap is
reported as unfinished. Keep the per-voxel log, the seed-step log and the maps before and
after the pass in the cache, in the layout of step 9, so that later steps can replay seed
steps and pair voxels.

**What is reported.** Evaluations by stage, grid fallbacks, wrong voxels before and after the
pass with Wilson 95% intervals, the median error over all voxels, the boundary share of the
voxels. One table with `held` next to the three other samples. One row pair in
`docs/autoresearch/pareto/pareto_rows.json` (before and after the pass).

**What must not happen.** No other run on `held`: no other arm, no replay with another
refinement, no scoring by a model. The sample stays unread for choices until the plans of the
steps that use it are committed with their predictions.

**Predictions, written before the run.**

1. `S8` on `held` with the pass needs between 20 and 40 million evaluations (more grains mean
   more seed steps; `mid` with 120 grains needed 26.44, the first gradient sample with 497
   grains 18.00).
2. Wrong voxels after the pass: between 30 and 120 (0.12% to 0.49%).
3. Grid fallbacks: between 5 and 40.

**Budgets.** Sample and images: 15 minutes. Baseline: at most 120 minutes, at most 2 heavy
processes. All wall times are contended and are not timing results.

**Amended on 2026-11-11 (D-20261111-N17).** The run reached the 120-minute cap unfinished (all
24,570 voxels had a record; no map was written). The cap was a budget taken from prediction 1,
which missed. The same run goes again with a cap of 480 minutes and nothing else changed.

**Files.** `scripts/autoresearch/ar3/held_*.py`, a test of the sample builder's pure parts in
`tests/test_ar3_held_sample.py`, results in `benchmarks/autoresearch/ar3/held/` (under 5 MB
tracked; the image stack and the caches stay untracked), a short page
`docs/autoresearch/findings_ar3_held.md`.

**Not in this step.** Any comparison of methods on `held`. A timing. Realistic images.

## AR4: learning where the search is uncertain

The first plan put the learned part at the proposer (H4): a network that proposes orientations
where the voting of ray indexing fails. The results since then say the voting does not fail on
these samples. In AR1 the truth was among the first 5 candidates for 69 of 69 and 200 of 200
voxels. In AR5 step 1, in all 31 grid-fallback cases of run `S8` a candidate within 2 degrees
of the truth was costed and refined. A network that proposes has nothing to add there. Ray
indexing on realistic images (overlap and detector noise) has not been tested; if it fails
there, H4 comes back.

What does fail is the test of a result. The only test is the hit ratio, and the hit ratios of
right and wrong results overlap. A right seed that ends at a hit ratio of 0.66 to 0.80 fails
the seed test and costs a full grid search (5.45 and 5.00 million evaluations on the two
samples that have fallbacks, 30% and 19% of `S8` with the pass). A flag at 0.95 catches every
wrong voxel but marks 8.8% and 15.1% of all voxels; the check `V1` on the flagged voxels of
`mid` then spent 27,717,299 of its 27,882,533 evaluations on right voxels. A seed test
lowered to 0.78 would accept 14 right ends and 1 wrong end of the 31 (both from the review
of AR5 step 1). So the learned part moves from proposing
to verifying: "is this orientation right for this voxel, or can more search find a better
one?" This is the owner's goal 5 (predict whether more optimization can help) with the tool
of goal 2 (a model fitted once per setup and material, reused for every voxel).

### AR4 step 1: a learned verifier, measured offline against the hit ratio

Task `ar4a`. Scripts and results only: no file under `icenine_py/icenine/` changes. It starts
after AR5 step 1 is merged, because it reads that step's cases and vote counts.

**Terms.** A **pair** is one voxel with one orientation. A pair is **right** if the
orientation is within 1 degree of the voxel's truth (the symmetry-reduced angle of the
library; nothing here is specific to cubic crystals). A **verifier** reads numbers that one
evaluation of the pair gives, with no truth, and returns a score from 0 to 1 that the pair is
right. A pair is **flagged** when its score is under a threshold. Train samples: `big` and the
first gradient sample. Test sample: `mid`; no model is trained and no threshold is picked on
it, and the file of thresholds is committed before any script scores a pair of `mid`.

**The idea to test.** At a boundary voxel the footprint lies partly in a neighbouring grain.
The right orientation then lights part of every predicted peak. A wrong orientation (most
often the neighbour grain's) lights some peaks fully and misses others. The hit ratio adds
all pixels together and cannot tell these apart; the per-peak pattern may.

**Inputs of a pair** (from one evaluation and the stored ray lists):

- per predicted peak: its lattice-plane family as a rank by length of the reciprocal vector,
  the detector distance, the rotation angle, and the share of its pixels that fall on observed
  intensity;
- summary numbers: hit ratio, cost, confidence, the share of peaks with any lit pixel, the
  lowest quarter of the per-peak shares, the number of attached rays the orientation explains,
  and the best vote count among the rays it does not explain (AR5 step 1's number).

If the cost function does not give the per-peak shares, a script computes them from the same
overlap data; the library is not changed.

**Models.**

- `B0`: the hit ratio alone (today's rule).
- `B1`: a logistic regression on the summary numbers.
- `N1`: a small set network on the per-peak inputs (a shared layer per peak, mean and max
  over the peaks, then two layers; under 50,000 parameters), plus the summary numbers.
  PyTorch, CPU, one training run with generator seed 0, at most 10 minutes. The layer sizes,
  the loss and the number of epochs are written in the script and committed before training.

**Pairs.** For each sample: (i) every voxel with its `S8` orientation before the boundary
pass and after it; (ii) at every boundary voxel, each distinct neighbour orientation (the
rival-grain case); (iii) the refined ends of every stored seed step, and the costed candidates
of every 10th seed step; (iv) the truth turned by 0.1, 0.3 and 0.5 degree about a random axis
at every 20th voxel (right, but not exact). Labels come from the truth. One evaluation per
pair; the count of pairs and of evaluations is reported.

**Two uses, each with its own test on `mid`.**

- *Use 1: the flag of unfinished voxels.* On the `S8` maps after the pass. For each model the
  threshold is the one that flags the fewest voxels while flagging at least 95% of the wrong
  voxels of the train samples, pooled. On `mid`: wrong voxels flagged (of 31, with a Wilson
  95% interval) and the share of all voxels flagged. Pass rule for a model: it flags at least
  as many of the 31 as `B0` at `B0`'s own threshold, and at most half as many voxels.
- *Use 2: accept a seed without a grid search.* On the stored seed steps whose refined end
  failed the seed test (the 31 fallback cases; 16 on a train sample, 15 on `mid`) and, as a
  control, on the seed steps that passed. A seed is accepted when its score is at or above the
  threshold of use 1. Count on `mid`: right ends accepted (each one is a grid search saved;
  give the stored evaluations of those searches as a what-if, not as a result), wrong ends
  accepted, and passing seeds that would now be refused. These are small counts and are read
  as counts.

**Predictions, written before any training.**

1. `B1` does not pass use 1 (the summary numbers add little to the hit ratio).
2. `N1` passes use 1 on `mid`.
3. Use 2 on `mid`: `N1` accepts at least two thirds of the right ends among the 15 fallback
   cases and at most one wrong end.
4. The share of peaks with any lit pixel alone (one number, no training) separates right from
   wrong at boundary voxels better than the hit ratio does. If it does as well as `N1`, the
   page says so and the network is not needed.

**Budgets.** Inputs of all pairs: 45 minutes, at most 2 heavy processes. Training: 10
minutes a model. No whole-layer run.

**Files.** `scripts/autoresearch/ar4/verify_*.py`, tests of the pure functions (inputs,
pooling, threshold pick) in `tests/test_ar4_verifier.py`, results in
`benchmarks/autoresearch/ar4/` (under 5 MB tracked; the pair inputs stay untracked in a
cache), page `docs/autoresearch/findings_ar4_verifier.md`.

**Not in this step.** The verifier inside the BFS, a whole-layer run, a key. If a model
passes, step 2 wires it in as a script layer and runs the three samples. Realistic images.

## AR5: the stop rule

The first line of this section in the first plan read: "Combine AR0's thresholds with the
unexplained-ray vote count of H3; test as a stop rule inside the arms that survive." The arm
that survives is run `S8` with the boundary pass. Step 1 below tests stop rules at the three
places where that run still spends evaluations without knowing whether they can help.

### AR5 step 1: stop rules for the boundary pass, the grid fallbacks and the finished map

Task `ar5a`. No file under `icenine_py/icenine/` changes. No key is added and no default
changes. Script code only, as in AR3 steps 4 to 9.

**Why.** The owner asked for a way to tell whether more optimization can still improve a
voxel, so that the search can stop earlier. Reading the table: one row per sample; the
evaluations of run `S8` with the boundary pass, by stage (from `findings_ar3_start.md`, table
`p3_stages`).

| sample | seed steps | grid fallbacks | neighbour fits | revisits | boundary pass | total |
|---|---|---|---|---|---|---|
| `big` | 2,856,712 | 0 | 3,491,437 | 853,714 | 1,138,565 | 8,340,428 |
| first gradient sample | 3,008,757 | 5,449,501 | 2,798,001 | 3,060,730 | 3,683,241 | 18,000,230 |
| `mid` | 12,007,352 | 5,002,750 | 3,859,661 | 2,517,272 | 3,056,039 | 26,443,074 |

- *The boundary pass* is 12% to 20% of the total. It makes one try of up to 61 evaluations
  from every distinct neighbour orientation. In step 5 only 739 of 60,290 tries changed a
  voxel (run `D4`, first gradient sample). Nothing tells a try in advance that it cannot win.
- *The grid fallbacks* are 19% and 30% of the total on the two samples that have them: 15 and
  16 full grid searches of about 333,000 evaluations each. AR0 found that a stop between the
  grid levels can save at most 39% of one search, and that the best cost after a level does
  not tell right from wrong. So the saving has to come from cheaper tries before the grid.
- *The finished map.* After the pass 24, 63 and 31 voxels are wrong. Nothing in the run says
  which voxels are finished and which would gain from more search. AR0 found on the pilot
  that a hit ratio under 0.93 flags all 37 wrong voxels and 15 right ones. That was a small
  region without gradients.

The start test for the neighbour fits (run `S9`, step 9) is not repeated: it broke more
voxels than it fixed.

**Terms.** As in AR3 step 9 (BFS, neighbour fit, revisit, reseed, seed step, seed test, grid
fallback, boundary pass, hit ratio, right, wrong, error). *Try*: one local optimization of a
boundary voxel from one neighbour's orientation, in the boundary pass. *Start score*: the hit
ratio and cost of a try's start, one evaluation. *Try test*: a rule that reads the start
score and skips the try. *Net gain of the pass*: wrong voxels before the pass minus wrong
voxels after it. *Rung*: a cheaper search tried before the full grid search. *Ladder*: the
rungs in order, then the grid search if all fail. *Flag*: a rule without the truth that marks
a voxel as possibly not finished. *Unexplained rays*: the rays attached to a voxel that the
voxel's orientation does not explain (`explained_mask` in `icenine/ray_indexing.py`).

**Samples and stored data.** `big`, the first gradient sample and `mid`, clean images, Q-max
8, the stored voxel order. The stored maps, tries (`tries.npz`), seed-step logs
(`e1_log.json.gz`), records and errors of `S8` are in the worktrees of steps 8 to 10 that
`STATE.md` lists; they are linked, not copied into git. The truth is read only by scoring
code. Settings are picked on `big` and the first gradient sample and tested on `mid`, with
nothing changed after the pick. The pick is written to a results file and committed before
the script that reads `mid` is run, so that git shows the order.

**Part 1: a try test for the boundary pass.** For every stored try of `S8`'s pass on the
three samples, compute the start score (one evaluation each). The voxel's score before the
pass is stored. The rules and their settings are fixed here; none is added after a result.

| rule | what it does | settings |
|---|---|---|
| `P1` | skip a try whose start hit ratio is more than `d` under the voxel's hit ratio | `d` = 0, 0.02, 0.05, 0.10, 0.20 |
| `P2` | per voxel, keep only the `m` tries with the highest start hit ratio | `m` = 1, 2 |
| `P3` | `P1` and `P2` together | the picked `d` with `m` = 1, 2 |

A rule's result is computed from the stored tries: the voxel keeps the best end among its
kept tries and its old orientation (rule `hit`, as in the pass). Its cost is one evaluation
per try plus the evaluations of the kept tries. Reported per rule and sample: tries kept,
evaluations, wrong voxels after the pass, and paired against the full pass the voxels fixed
and broken. *The pick*: the setting with the fewest evaluations among those that lose at most
5% of the net gain of the pass, pooled over the two pick samples. *The pass rule*: on `mid`
the picked setting cuts the evaluations of the pass by at least half and loses at most 5% of
its net gain.

**Part 2: a ladder before the grid search.** The cases are the seed steps of `S8` that ended
in a grid fallback on the first gradient sample and on `mid` (16 and 15; an entry of
`e1_log.json.gz` with `stop` = -1 and no refined hit ratio of 0.80). These are few cases; the
page says so, lists every case, and picks no setting on them.

1. *What failed.* For each case: the number of candidates the indexing gave; the rank of the
   first candidate within 2 degrees of the truth; whether it was among the 20 costed and
   among the 5 refined; the end of the grid search (hit ratio, error). A helper's reading of
   the caches gave 15 of 16 and 13 of 15 grid ends right; the step recomputes this.
2. *The rungs*, each run on every case, then as a ladder:
   - `R1`, more candidates: cost the first 100 candidates, refine the 5 cheapest with the
     seed step's local optimization, stop at the first hit ratio of 0.95.
   - `R3`, grid level 0: run only level 0 of the grid search (its screen, re-cost and quick
     MC, as the library does them), refine its 5 best orientations the same way, stop at the
     first hit ratio of 0.95. If the library cannot continue a search after level 0, the
     ladder's last rung is a full grid search and its level 0 is counted twice; the page
     says which.
   - A rung passes a case if its end passes the seed test (hit ratio 0.80).
   - A rung that moves the seed to another voxel of the cluster is not in this step: it
     needs the BFS's cluster list.
3. *Reported* per case and rung: evaluations, seed test passed or not, right or wrong. For
   the ladder: the total against the stored grid fallbacks of the same cases.

*The pass rule*: over all cases the ladder needs at most two thirds of the evaluations of
the grid fallbacks it replaces, and ends right in at most one case fewer than the grid
search did. Every case the grid search got right and the ladder got wrong is listed.

**Part 3: is a voxel finished?** On the maps of `S8` after the pass, three samples.

1. *A flag from the hit ratio.* The hit ratio of each voxel after the pass comes from the
   stored scores and tries. For thresholds 0.50 to 0.95 in steps of 0.05: voxels under the
   threshold, the wrong ones among them, the wrong ones above it. Also by the voxel's source
   (seed, neighbour fit, revisit, reseed) and for boundary voxels alone. *The pick*: the
   lowest threshold that flags at least 80% of the wrong voxels pooled over the pick samples.
2. *The vote count of H3.* For every flagged voxel and for every 20th unflagged boundary
   voxel by index: index the unexplained rays (`index_rays`) and take the highest vote
   count. Shown as two strips of dots, wrong and right; then the counts at the cut that
   separates best, labelled as chosen after the result.
3. *A check that costs evaluations*, `V1`: at each flagged voxel run `S8`'s seed step (20
   candidates costed, up to 5 refined, stop at hit ratio 0.95) and keep its end if its hit
   ratio is higher than the voxel's. Reported: evaluations, fixed and broken against the map
   after the pass. *The pass rule*: on `mid`, at the picked threshold, `V1` fixes more voxels
   than it breaks, fixes at least a quarter of the flagged wrong voxels, and costs under 5%
   of the run's total.
4. *What "finished" would mean.* The page gives, for the picked threshold, the wrong voxels
   among the unflagged ones with a 95% interval. That count is what a user who stops at the
   flag leaves behind.

**Part 4: whole-layer run `T1`.** `S8` with every rule that passed: the ladder at each grid
fallback, the try test in the pass, `V1` after the pass. On the three samples, one run each,
cap 90 minutes per run. Only if a rule passed. `S8` is run by the stored script runner of
steps 8 and 9, or by the keys of step 10 if that step is merged when the run starts. On `big`
the run has no grid fallback, so before the pass it must equal the stored `S8` exactly; if it
does not, stop and report. Accounting as in step 9: evaluations by stage, seeds, reseeds,
grid fallbacks, wrong voxels before and after the pass with 95% intervals, the error
percentiles, counts paired against `S8`. The runs go into
`docs/autoresearch/pareto/pareto_rows.json` as rows of a new run `T1`.

**Written down before the runs.** (1) Some setting of the try test removes at least 70% of
the pass's evaluations on the pick samples and loses at most 5% of its net gain; it passes
on `mid`. (2) In at least half of the fallback cases a candidate within 2 degrees of the
truth is among the first 100 and not among the 20 costed. (3) The ladder passes its rule.
(4) A hit-ratio threshold flags at least 80% of the wrong voxels and at most 5% of all
voxels, on the pick samples and on `mid`. (5) The vote count of the unexplained rays does not
separate wrong from right among flagged voxels: a ray's line crosses many voxels, so rays of
the neighbouring grain are attached to every boundary voxel. (6) `T1` with the pass needs
under 14 million evaluations on the first gradient sample (`S8`: 18.00) and under 22 million
on `mid` (`S8`: 26.44), with intervals of wrong voxels that overlap `S8`'s. The page says for
each whether the result agrees.

**What the result decides.** A rule that passes on stored data and in `T1` is a candidate
for an opt-in key in a later step; what a key recommends stays the owner's choice. A rule
that fails is reported as failed, with its numbers. If prediction 5 holds, the vote count
of H3 is dropped as a stop rule and the page says what replaced it.

**Rules of the study.** One run per arm; no extra random seeds. No p-values: counts, paired
fixed and broken, 95% intervals. Intervals treat voxels and tries as independent although
they cluster in grains; the page says so. Every number chosen after a result is listed as
such. Nothing here is specific to cubic symmetry or to special misorientations.

**Stop rules of the task.** Part 1: 45 minutes of compute. Part 2: 90 minutes. Part 3: 60
minutes. Part 4: the caps above. A run that hits its cap is reported with what it reached,
not repeated. The general stop rule holds as well.

**Not claimed.** That C++ behaves the same. That the gradient model matches a real sample.
A change of any default. A tested wall time unless the timing preflight passes. That three
samples stand for all samples. That realistic images give the same thresholds.

**Files.** Scripts in `scripts/autoresearch/ar5/` (`stop_*.py`); results in
`benchmarks/autoresearch/ar5/`, slim, under 5 MB; findings in
`docs/autoresearch/findings_ar5_stop.md` (returned as text if the write is refused).
`MIGRATION_HISTORY.md` gets the summary. The page says "local optimization", defines each
term at first use, and leads with a plot: evaluations by stage for `S8` and `T1` on the
three samples (or, if no rule passed, the tries kept against the net gain lost for each
setting of the try test); then the tables.

### AR5 step 2: two long local optimizations that rarely change the result

Task `ar5b`. Scripts and results only: no file under `icenine_py/icenine/` changes, no key is
added, no default changes. Planned on 2026-11-11 from stored numbers; step 1's result (no rule
passed) does not change it.

**Why.** Two stages of run `S8` run the stock MC optimizer (`FindOptimal`, then the
variance-minimizing stage) from a start that is already close. Stored numbers, not yet put side
by side on a page:

- The seed step refines each indexed candidate with the stock MC optimizer. The candidate is a
  least-squares (Procrustes) output; in AR1 its median distance from the truth was 0.038 and
  0.054 degree. One refinement costs 6,129 evaluations on average on `big`, 1,266 on the first
  gradient sample and 3,723 on `mid`. The seed step is 2,856,712, 3,008,757 and 12,007,352
  evaluations; costing the candidates is a tiny part of that (6,580 on `big`).
- A neighbour fit runs CMA first (61 evaluations). If the CMA end cost is not below the start
  cost, the stock MC fit runs from the same start (the "MC fallback"). In `S8` it ran in 3,197,
  3,903 and 4,568 fits (9.8%, 7.3%, 8.8%) at a mean of 738, 668 and 706 evaluations: about
  2.36, 2.61 and 3.23 million. In the runs where its wins were counted it ended lower than the
  CMA end in 6 of 3,974 (`D4`, first gradient sample), 4 of 3,194 (`D4`, `big`) and 3 of 3,235
  (`S5`, `big`).

Together that is about 62%, 31% and 58% of the evaluations of `S8` with the boundary pass
(8,340,428, 18,000,230 and 26,443,074). These shares are my arithmetic from the stored
numbers; Part 0 recomputes them from the files.

**Terms.** Evaluation, seed step, seed test (hit ratio at least 0.80), early stop (0.95), grid
fallback, neighbour fit, revisit, boundary pass, wrong (more than 1 degree from the truth): as
in `findings_ar3_keys.md`. Pick samples: `big` and the first gradient sample. Test sample:
`mid`. A pick file is committed before any script reads `mid`. No script runs on `mid` before
that commit, not even to fill a cache (the review of step 1 found such a run; the page must
be able to say there was none).

**Part 0: the stored numbers.** From the stored `S8` runs of the three samples: per sample,
the number of seed refinements and their evaluations, the number of MC fallbacks, their
evaluations and how many ended lower than the CMA end. If a stored run does not hold the last
count, say so and take it from Part 2's base run. One table.

**Part 1: a seed refinement whose step matches the error of its start (replay).** For every
stored seed step of `S8` (329, 1,279 and 1,611 steps), index the voxel again and check what
the log holds (the number of candidates, the number costed). Then refine with each arm, in `S8`'s
order, with `S8`'s early stop (the first end with a hit ratio of at least 0.95 is kept; if none,
the lowest cost of up to 5 ends):

- `M`: the stored MC refinement. Its numbers come from the log. Rerun it on every 10th step
  and check that evaluations and hit ratio equal the stored ones; if they do not, stop and
  report.
- `C60`: a CMA fit, first step 0.05 degree, at most 60 evaluations (the neighbour fit's
  settings).
- `C200`: the same with at most 200 evaluations.
- `C0`: no refinement; the candidate as costed.

Per step and arm: evaluations, end hit ratio, error against the truth, seed test passed or
not. Pick rule (on the pick samples, pooled, paired by step): the cheapest arm for which (a)
the steps that pass the seed test are at most 1% of all steps fewer than with `M`, and (b) the
steps that end within 1 degree of the truth are at most 1% of all steps fewer than with `M`.
Pass rule on `mid`: (a) and (b) hold there, and the seed-step evaluations fall by at least
80%. Report the paired counts (passes with `M` only, with the arm only), no p-values.

**Part 2: the whole layer.** One base check first: the script layer with both changes off must
give stored `S8` exactly on `big` before the pass (7,201,863 evaluations, 67 wrong). Then run
`S10` = `S8` with (i) the picked seed refinement, if Part 1's pick passed on `mid`, and (ii)
no MC fallback: a neighbour fit or revisit returns what the CMA fit returned. Three samples,
with the boundary pass, at most 90 minutes a run; a run at its cap is reported as unfinished.
Evaluations by stage, grid fallbacks, wrong voxels before and after the pass with Wilson 95%
intervals, fixed and broken against `S8` paired by voxel, median error.

Pass rule for `S10`, on each sample: at least 20% fewer evaluations than `S8` with the pass;
broken minus fixed at most 5 voxels; median error at most 0.005 degree above `S8`'s. If it
fails on a sample, run the two changes alone on that sample (`S10a`: seed refinement only;
`S10b`: no MC fallback only) to see which one fails. The rows of every finished run go into
`docs/autoresearch/pareto/pareto_rows.json`.

**Predictions, written before any run.**

1. Part 0: in `S8` the MC fallback ended lower than the CMA end in under 1% of its runs on
   every sample that holds the count.
2. `C60` is the pick and passes on `mid`; it cuts the seed-step evaluations by at least 90% on
   all three samples.
3. `C0` fails rule (a): without a refinement too few candidates pass the seed test.
4. `S10` with the pass needs under 4.0 million evaluations on `big` (`S8`: 8.34), under 14
   million on the first gradient sample (18.00) and under 13.5 million on `mid` (26.44), and
   passes its rule on all three.
5. The number of grid fallbacks of `S10` is within 2 of `S8`'s on every sample (0, 16, 15).

**Budgets.** Part 0 minutes; Part 1 30 minutes; Part 2 at most 90 minutes a run, at most 3
heavy processes at a time. All wall times are contended and are not timing results.

**Files.** `scripts/autoresearch/ar5/fit_*.py`, tests of the pure functions in
`tests/test_ar5_fit_rules.py`, results in `benchmarks/autoresearch/ar5/fit/` (under 5 MB
tracked), page `docs/autoresearch/findings_ar5_fits.md`.

**Not in this step.** A library key (a later step, if `S10` passes). A timing. A real sample
or realistic images. The try test `P2 m=2` of step 1 (it lost no net gain on any sample and
cut 16.6%, 30.2% and 38.9% of the pass, under the plan's 50% rule): whether it becomes an
opt-in key is the owner's choice.

### AR5 step 3: a mixed seed refinement, tested on whole layers

Task `ar5c`. Scripts and results only: no file under `icenine_py/icenine/` changes, no key is
added, no default changes. It starts after AR5 step 2 is merged, because it uses that step's
script layer (`fit_*.py`), its replay cache and its base check. Planned on 2026-11-11
(D-20261111-N14).

**Why.** In step 2 no CMA seed refinement met the 1% rule alone: the 200-evaluation arm
(`C200`) lost 60 seed-test passes on the 1,608 seed steps of the pick samples (limit 16). The
review of step 2 then read the stored replay step by step and found that the losses are the
steps where the CMA end stays under a hit ratio of 0.95. A what-if from those stored numbers
(arithmetic, not a run): refine with `C200` first, and run the stock MC optimizer from the same
start only when the CMA end hit ratio is below 0.95.

| what-if, from the stored replay | MC (`M`) | mixed |
|---|---|---|
| seed-test passes, pick samples (1,608 steps) | 1,592 | 1,596 |
| ends within 1 degree, pick samples | 1,583 | 1,589 |
| seed-step evaluations, pick samples | 5,865,469 | 2,501,991 |
| seed-test passes, `mid` (1,611 steps) | 1,596 | 1,600 |
| seed-step evaluations, `mid` | 12,007,352 | 5,643,374 |

Reading the table: one row per count, the stored MC refinement against the mixed rule. The
mixed rule keeps the passes and needs 57.3% and 53.0% fewer seed-step evaluations.

**What this what-if is not.** It was read after all three samples were seen, `mid` included,
so no sample is held out for the seed rule. It treats each seed step alone: in a whole layer a
different seed end changes what the neighbours start from. The test of this step is therefore
the whole-layer run, with the predictions below written before it.

**Terms.** As in step 2 and in `findings_ar5_fits.md`. **Mixed seed refinement** (`X60`,
`X200`): for each indexed candidate, in `S8`'s order, a CMA fit with a first step of 0.05
degree and at most 60 or 200 evaluations; if its end hit ratio is at least 0.95 that end is
kept and the seed step stops (`S8`'s early stop). If not, the stock MC optimizer runs from the
same candidate, and the candidate's end is the one of the two with the lower cost. If no end
reaches 0.95, the lowest cost of up to 5 ends is kept, as in `S8`. Evaluations count both.

**Part 1: the mixed rule from the replay cache (no new evaluation).** A committed script and
tests replace the review's hand count. For `X60` and `X200` and each sample: seed-test passes,
ends within 1 degree, seed-step evaluations, the paired counts against `M` (passes with `M`
only, with the arm only), and the share of candidates where the MC optimizer still ran. Pick:
the arm with the fewer seed-step evaluations on the pick samples, pooled, among those that
meet step 2's rules (a) and (b) there. The pick is committed before any whole-layer run. If
neither arm meets (a) and (b), stop and report.

**Part 2: the whole layer.** Run `S11` = `S8` with the picked mixed seed refinement and with
no MC fallback after CMA in a neighbour local optimization or revisit (step 2's run `S10` plus
the seed rule). Three samples, with the boundary pass, at most 90 minutes a run; a run at its
cap is reported as unfinished. Before the runs, check on 20 seed steps of `big` that the
script layer's mixed refinement gives the ends and evaluation counts that Part 1 computed from
the cache; if not, stop and report. Report evaluations by stage, grid fallbacks, wrong voxels
before and after the pass with Wilson 95% intervals, fixed and broken against `S8` and against
`S10` paired by voxel, and the median error over all 24,570 voxels.

Pass rule for `S11`, on each sample, the same as step 2's and not re-tuned: at least 20% fewer
evaluations than `S8` with the pass; broken minus fixed against `S8` at most 5 voxels; median
error at most 0.005 degree above `S8`'s. If it fails on a sample, run `S11a` there (the mixed
seed refinement with the MC fallback kept) to see which change fails. The rows of every
finished run go into `docs/autoresearch/pareto/pareto_rows.json`.

**Predictions, written before any run.**

1. `X200` is the pick: `X60` reaches 0.95 too rarely, so the MC optimizer runs on most of its
   candidates.
2. `S11` with the pass needs under 4.4 million evaluations on `big` (`S8`: 8.34, `S10`: 5.90),
   under 14.2 million on the first gradient sample (18.00, 15.17) and under 17.5 million on
   `mid` (26.44, 22.67), and passes its rule on all three.
3. The number of grid fallbacks of `S11` is within 2 of `S8`'s on every sample (0, 16, 15).
4. The seed-step evaluations of `S11` in the whole layer are within 15% of Part 1's count for
   the same sample (the path changes little).

**Budgets.** Part 1 minutes; Part 2 at most 90 minutes a run, at most 3 heavy processes at a
time. All wall times are contended and are not timing results.

**Files.** `scripts/autoresearch/ar5/mix_*.py` (reusing `fit_common.py` and `fit_run.py`; a
change to a step 2 script must leave step 2's results reproducible), tests of the pure
functions in `tests/test_ar5_mix_rules.py`, results in `benchmarks/autoresearch/ar5/mix/`
(under 5 MB tracked), page `docs/autoresearch/findings_ar5_mix.md`.

**Not in this step.** A library key (the next step, for both changes together, if `S11`
passes). A timing. A real sample or realistic images. Other thresholds than 0.95 and other
budgets than 60 and 200: the two values come from keys that already exist, and tuning them on
three samples that have all been read would fit the samples.

### AR5 step 4: the mixed seed refinement and the learned verifier on the held-out sample

Task `ar5d`. Scripts and results only: no file under `icenine_py/icenine/` changes, no key is
added, no default changes. Drafted on 2026-11-11 (commit 17a27c1), while the `S8` baseline of
AR3 step 11 was still running: the predictions below were written before any number of `held`
was read (D-20261111-N16).

**Amended on 2026-11-11 (D-20261111-N18).** The baseline then reached its cap of 120 minutes
with 48.3 million evaluations and 61 grid fallbacks (both lower bounds) and runs again with a
cap of 480 minutes (D-20261111-N17). Two things changed here and nothing else: the cap of a
run (150 to 480 minutes), and the start. The step no longer waits for step 11 to merge: its
task branch merges the branch of step 11 first, its runs go next to the baseline run, and it
merges after step 11. No comparison is read before the baseline has finished. The
predictions were not touched.

**Amended on 2026-11-11 (D-20261111-N20), after the review.** The runs `S10` and `S11` stand.
The run `S11v` did not test the planned rule: each of the 24 seeds that `N1` accepted was then
refused by the library's own seed check (the same hit-ratio threshold as the seed test), so
its grain did not grow. The run only skipped grid fallbacks. That is a defect of the script
layer, not a result about the verifier. The fixes round repairs the hook in the script layer
(no file under `icenine_py/icenine/` changes) so that an accepted seed is kept and its grain
grows. The repair is checked on the stored grid-fallback cases of `mid`, not on `held`. Then
`S11v` runs once more on `held` as `S11v2`: the same frozen model and cut, the same rule, the
same pass rule, a cap of 480 minutes. Nothing is tuned and no setting is chosen from the first
run. The first run stays on the page, labelled as a defective arm, and its rows leave the
plots. The rule part of prediction 5 and prediction 6 are judged on `S11v2`; the predictions
are not rewritten. The page says that `S11v2` ran after the numbers of the defective run on
`held` had been seen. This is a repeat because of a bug. It is not a repeat with a new random
seed to firm up a comparison.

**Why.** Run `S11` (AR5 step 3) needs 32% to 52% fewer evaluations than `S8` on three samples,
but its seed rule was found in a replay that included one of them (`mid`). The learned
verifier `N1` (AR4 step 1) accepted 7 of 10 right and 0 of 5 wrong seeds that had failed the
seed test on `mid`, but that is 15 cases and arithmetic, not a run. `held` (250 grains, a new
generator seed) has been read for no choice. This step tests whether both results carry over.

**Terms.** As in `findings_ar5_mix.md` and `findings_ar4_verifier.md`. **Seed test**: the
hit-ratio test that a refined seed must pass before its grain is painted. **Grid fallback**:
the full adaptive grid search that runs when the seed step fails that test. **Added
seed-accept rule** (the OR rule of AR4 step 1): a seed that fails the seed test is still
accepted when the score of `N1` for its end is at or above the stored cut; the grid fallback
is then skipped.

**Frozen inputs.** Nothing is tuned on `held`. The mixed seed refinement is `X200` as in step
3. `N1` is the stored model of AR4 step 1 (`models.npz`) with its stored cut. The runs use the
script layer of AR5 steps 2 and 3.

**Runs, all on `held`, clean images, with the boundary pass, one run each.**

| run | what it is |
|---|---|
| `S10` | `S8` with no MC fallback after CMA in a neighbour local optimization or revisit |
| `S11` | `S10` plus the mixed seed refinement `X200` |
| `S11v` | `S11` plus `N1` as an added seed-accept rule |

Reading the table: one row per run; each row adds one change to the row above it.

At most 480 minutes a run, at most 3 heavy processes next to the baseline run; a run at its
cap is reported as unfinished. All wall times are contended and are not timing results.

**Checks before the runs.** (1) On 20 seed steps of `big`, the script layer still gives the
ends and counts of step 3. (2) On the 15 stored grid-fallback cases of `mid`, the score of
`N1` computed inside the run layer equals the stored offline score to 1e-6. If either check
differs, stop and report. The evaluations that `N1` needs for its inputs are counted and
reported; they are part of the total of `S11v`.

**What is reported.** For each run, next to the `S8` baseline: evaluations by stage (seed
step, grid fallback, neighbour local optimization, revisit, boundary pass); grid fallbacks;
wrong voxels before and after the pass with Wilson 95% intervals; fixed and broken against
`S8` paired by voxel (and `S11v` against `S11`); the median error over all voxels. The grid
fallback cases are listed one row per seed step, paired across the runs by the seed voxel:
which runs fell back there, and for `S11v` the score, whether the seed was accepted, and
whether the accepted seed was within 1 degree of the truth. The rows of every finished run go
into `docs/autoresearch/pareto/pareto_rows.json`.

**Pass rule, the same as steps 2 and 3, not re-tuned.** Against `S8` with the pass: at least
20% fewer evaluations; broken minus fixed at most 5 voxels; median error at most 0.005 degree
higher. For `S11v` the same three against `S8`, and in addition against `S11`: fewer
evaluations, and broken minus fixed at most 5 voxels.

**Predictions, written before the baseline of `held` was read.**

1. `S10` needs 12% to 20% fewer evaluations than `S8` and so misses the 20% rule (it gave
   29%, 16% and 14% on the three earlier samples).
2. `S11` needs 30% to 45% fewer evaluations than `S8` and meets all three parts of the rule.
3. `S11` has no more grid fallbacks than `S8`.
4. Wrong voxels of `S11` after the pass are within 10 of `S8`'s.
5. `S11v` accepts at least half of the seeds that fail the seed test, at most one of the
   accepted seeds is more than 1 degree from the truth, and it meets its rule against `S11`.
6. The saving of `S11v` against `S11` is under 15% of the total of `S11` (few cases).

**What follows.** If `S11` passes: the next step adds one opt-in key pair to the library for
the mixed seed refinement and "no MC fallback", with the tests of the library path. If it
fails: the failing stage is read from the split, and no key is added. `S11v` decides only
whether the verifier gets a whole-layer step of its own. After this step `held` is no longer
unread for these rules.

**Files.** `scripts/autoresearch/ar5/held_*.py` (reusing `fit_*.py`, `mix_*.py` and
`ar4/verify_lib.py`; a change to an earlier script must leave its results reproducible),
tests of the pure functions in `tests/test_ar5_held_rules.py`, results in
`benchmarks/autoresearch/ar5/held/` (under 5 MB tracked), page
`docs/autoresearch/findings_ar5_held.md`.

**Not in this step.** Any tuning on `held`. A library key. A timing. Realistic images. A
second run of any arm with another random seed.


### AR5 step 5: opt-in keys for the mixed seed refinement and for no MC fallback

Task `ar5e`. A library step: two opt-in keys, both off by default. No default changes, no
golden file is regenerated, C++ is not touched. Planned on 2026-11-11 (D-20261111-N21), after
the review of AR5 step 4 (D-20261111-N20).

**Why.** Run `S11` is `S8` with two changes that so far exist only in scripts. With the
boundary pass it needs 51.9%, 31.6% and 39.2% fewer evaluations than `S8` on `big`, the first
gradient sample and `mid` (AR5 step 3), and 23.3% fewer on the held-out sample `held` (AR5
step 4, reviewed). It met the pass rule on all four. Decisions N14 and N16 said: no key before
the test on `held`. That test is done. A user of the library cannot switch the two changes on
today. This step moves them into the library, as AR3 step 10 did for the seed step of `S8`.

**Terms.** As in `findings_ar5_mix.md`. **Mixed seed refinement**: for each indexed seed
candidate, a CMA fit first; its end is kept when its hit ratio reaches the stop ratio; if not,
the stock MC optimizer runs from the same candidate and the end with the lower cost is kept.
**No MC fallback**: in a neighbour local optimization or revisit that starts with CMA, the CMA
end is returned and the MC optimizer does not run after it. **Parity**: a library run with the
keys gives the same evaluation counts, the same per-voxel records and the same map file, byte
for byte, as the stored script run.

**The two keys.**

- `SeedIndexCMAFirst n`. Default 0 (off). A value of 2 or more switches the mixed seed
  refinement on, with `n` as the largest number of CMA evaluations for a candidate (`S11` uses
  200). The first step of the fit comes from `CMASigma0` and the stop ratio from
  `SeedIndexStopHitRatio`; no new constant is added. The key is refused with a clear message
  unless `SeedByIndexing 1` is set and the stop ratio is above 0, and it is refused together
  with `PaintGrains 1`, `LocalOptimizer cma` or the hybrid optimizer. A value of 1 is refused.
- `BFSNeighborCMANoFallback 0|1`. Default 0. With 1 the CMA end of a neighbour local
  optimization or revisit is returned as it is. It is refused unless `BFSNeighborCMAFirst 1`
  is set.

Counters that exist only when a key is on: how often the CMA end did not fall below its start
cost, and per seed candidate whether the CMA end reached the stop ratio and whether the MC
optimizer ran. Evaluations are counted as in the scripts (a CMA fit of 200 counts 201: the
fit and the final overlap).

**Where.** `icenine/reconstructor.py` (the CMA-first fit; the seed step that builds the
refinement), `icenine/ray_indexing.py` (the loop over seed candidates), `icenine/config_file.py`
and `icenine/orientation_search.py` (the keys, their checks), following the pattern of step 10.
The stock paths (`local_optimization`, `refine_from_candidates`, `_reconstruct_seed`, the
branch without a stop ratio) are not edited.

**Order of work.** Tests first: `tests/test_bfs_ar5_keys.py` is written and committed before
the library code (parsing, defaults, every refused combination, the off path, the rule itself
on small synthetic cases with a fixed generator). Then the code. Then the parity runs.

**Parity checks.** All against stored runs; the stored files are read only.

1. Keys off: the full suite passes with no golden change, and the library run of `S8` on `big`
   still gives 7,201,863 evaluations and 67 wrong voxels before the pass.
2. `BFSNeighborCMANoFallback 1` alone, on top of the keys of `S8`, against the stored script
   run `S10` on `big` (4,767,135 evaluations before the pass, 36 seeds, 297 reseeds): parity.
3. Both keys (`SeedIndexCMAFirst 200`) against the stored script run `S11` on `big`
   (2,879,705 evaluations before the pass, 35 seeds, 295 reseeds): parity.
4. Both keys against the stored `S11` on the first gradient sample: parity. `big` has no grid
   fallback, so this run is the one that tests the mixed seed refinement next to grid
   fallbacks (11 in the stored run).
5. The boundary pass after the library runs of checks 3 and 4, started from the stored maps
   before the pass: the totals with the pass equal the stored ones (4,015,342 and 12,307,994).
   The pass inside a library run reads the unrounded map and may differ by single tries, as
   step 10 found; that difference is reported, not repaired.

If a parity check fails, the step stops and reports the first voxel that differs. Nothing is
tuned to make a count match. No run on `mid` or `held` is needed; none is made.

**Predictions, written before any code.**

1. Checks 1, 2 and 3 are exact at the first full run of each.
2. Check 4 is exact. If one check fails, it is this one (the order of evaluations around a
   grid fallback is the part that the small tests cover least).
3. Check 5 is exact from the stored maps; the pass inside the library run differs from the
   stored pass by under 0.1% of the pass evaluations on each sample.
4. The step needs between 25 and 50 new tests, and no existing test changes.

**What is reported.** The page `docs/autoresearch/findings_ar5_keys.md`: the two keys with an
example config block, each parity check with its counts, the predictions against the results,
and the limits. The keys go into the key table of `icenine_py/README.md` with a plain statement
that they are off by default and that no setting is recommended yet. A section in
`MIGRATION_HISTORY.md`. No new row on the plots of cost against accuracy: the runs are the
stored `S10` and `S11` again.

**Limits to state.** Parity is shown on two samples with clean images, one voxel order each.
Evaluation counts are costs; no timing was made. The keys were not run on a real sample or on
realistic images. Whether the documentation recommends them is the owner's choice.

**Budgets.** Each run on `big` about 15 minutes, the run on the first gradient sample at most
90 minutes; at most 2 heavy processes at a time (the rerun of AR5 step 4 may share the
machine). All wall times are contended and are not timing results.

**Not in this step.** A key for the learned verifier (AR5 step 4 decides whether it gets a
further test). A default change. A recommended setting. A timing. C++.

### AR5 step 6: the frozen seed-accept rule on two more whole layers, against a plain rule

Task `ar5f`. Direction: D-20261111-N23; plan approved by D-20261111-N24. Starts after AR5
step 4 (`ar5d`) merges (it reuses that step's script layer).

**Library changes:** none. No file under `icenine_py/icenine/` changes; no key, no default.

**Why.** On `held`, `S11v2` (`S11` plus the learned verifier `N1` as an added seed-accept rule)
needed 18.4% fewer evaluations than `S11`, broken minus fixed 4 (limit 5): one run, one sample.
The owner's question: can a cheap test say "more search will not help here"? The accepted seeds
sat just under the seed test, so a plain rule with no network is tested next to `N1`.

**Terms.** **Seed test**: a refined seed is accepted when its hit ratio (the share of predicted
pixels on observed intensity) is at least 0.80. **Grid fallback**: the global search (a search
with no starting guess) on the adaptive grid that runs when a seed fails the seed test.
**Added seed-accept rule**: a seed that fails the seed test is still kept when a score of its
end is at or above a cut; its grid fallback is skipped and its grain grows by the library's
neighbour rule. **Kept seed**: an accepted seed that the layer marked fitted. An end is
**right** when it is within 1 degree of the truth (the library's symmetry-reduced angle).

**Which samples are clean (question 1).** Reading line: one row per sample, what has touched it.

| sample | grains, gradient levels | `N1` trained | cut picked | rule run before | stored `S8`, `S11` with pass |
|---|---|---|---|---|---|
| `big` | 32, 0 to 4 degrees | yes | yes | no | yes, yes (0 grid fallbacks) |
| `grad` | 497, 0 to 1 degree | yes | yes | no | yes, yes (11 fallbacks in `S11`) |
| `mid` | 120, 0 to 4 degrees | no | no | offline only | yes, yes (12 fallbacks in `S11`) |
| `held` | 250, 0 to 4 degrees | no | no | `S11v`, `S11v2` | yes, yes |

Sources: `findings_ar3_big.md`, `findings_ar3_start.md`, `findings_ar3_gradient.md` (levels);
`findings_ar4_verifier.md` (trained on all 143,409 pairs of `big` and `grad`,
including the ends of `grad`'s 16 fallbacks; cut picked on the same voxels), `findings_ar5_mix.md`
(tables "stages", "fallback_cost"), step 4's page. No sample is fully clean. `mid` is the
cleanest, but AR4 step 1 scored its 15 `S8` fallback ends (7 of 10 right, 0 of 5 wrong accepted),
the OR reading of the rule came from that arithmetic, and step 4 used `mid` for two checks (the
score on 15 cases; a partial run that saw 3 failed seeds). Every sample has grains with
gradients; `grad`'s scores are in-sample. On `big` the arm would equal `S11` (no fallback).

**Choice.** Test sample `mid`. `grad` runs too, labelled in-sample for the scores: only its
whole-layer effects (growth of kept seeds, the pass, fixed and broken) are new. No new sample.

**The plain control (question 2).** AR4 step 1 compared `N1` with the hit ratio at its use-1 cut
(0.926; it accepts no failed seed by construction), a logistic rule and the share of lit peaks.
It did not compare a lowered seed test. The review of AR5 step 1 noted that 0.78 would accept 14
right and 1 wrong of 31 ends, but that pool included `mid`. The C++ rule
(`Src/BreadthFirstReconstructor.tmpl.cpp` lines 141-148) keeps a failed seed when its search
converged, and convergence there means cost under `MaxConvergenceCost` (`Src/CostFunctions.h`,
`ClassifyConvergence`). The stored configs set it to 0.0001. The lowest end cost of a failed seed
in the stored `S11` run on `held` is 0.227. So the C++ rule is the seed test itself here: its
run is `S11`, already stored. The control is therefore a lowered hit ratio:

- `S11h`: a failed seed is kept when its end's hit ratio is at least `h_c`; 0 extra evaluations.
- `h_c` is fixed in phase A, before any run, from a pool without `mid`: the 16 `S8` fallback ends
  of `grad` (`ar4/summary.json`, `use2.grad.rows_fallback`) and the 63 failed ends of `S11v2` on
  `held` (`n1_log.json`, errors from `summary.json`). Rule: the lowest cut that accepts no wrong
  end of the pool, midway between the highest wrong hit ratio and the next right one above it
  (near 0.791 by a look at those files). This cut is picked on accept data, `N1`'s for another
  use, so the comparison favours the control.

**Arms, one run each, clean images, Q-max 8, with the boundary pass.**

| run | samples | what it is |
|---|---|---|
| `S11n` | `mid`, `grad` | `S11` plus `N1` (stored `models.npz`, cut 0.74914), repaired hook `held_keep.py` |
| `S11h` | `mid`, `grad` | `S11` plus the plain rule with `h_c`, through the same hook |

Reading the table: `S11n` is step 4's `S11v2` renamed; both change only the accept decision of a
failed seed, in the script layer of `S11v2` (`held_run.RUNS`).

**Frozen before any run** (`ar5/ver/frozen.json`, committed): sha256 of `models.npz`, the cut,
`h_c` and its pool counts, `held_keep.LIB_SHA256`, `X200`, configs (stored `S11` plus the hook).

**Pass rule, unchanged.** Against `S8` with the pass: at least 20% fewer evaluations; broken
minus fixed at most 5 voxels; median error at most 0.005 degree higher. Against `S11`: fewer
evaluations, and broken minus fixed at most 5. Each arm on each sample is judged alone.
**"Plain does as well"** is fixed now: on `mid`, `S11h` passes and saves at least 80% of
`S11n`'s saving against `S11`. Then the page says the verifier is not needed on this evidence.

**Validity, asserted after each run and before any result is read.** Every accepted seed has
`seed_rejected` false and its "Seed kept" log line; accepted equals kept; no accepted seed is
refused downstream; `S11n` spends 2 evaluations per scored seed, `S11h` none; the run ended
before its cap; the hook's library guard passed. If any part fails, the stats script stops.

**What is measured**, per arm and sample, next to `S8` and `S11`: evaluations by stage (seed
step, grid fallback, neighbour local optimization, revisit, pass); grid fallbacks; failed seeds
scored, accepted and kept, with Wilson 95% intervals; the error of each accepted end; voxels
grown from each kept seed (from the library's line "k of n cluster voxels fitted", not the hook
field, which misses reseed growth) and how many are wrong before and after the pass; fixed and
broken by voxel against `S8`, `S11`, and `S11n` against `S11h`; wrong voxels with Wilson
intervals; median error over all voxels. For goal 5, each refused seed's log row also stores
its grid end, so the page can say how often the global search moved a failed end by more than 1
degree, split by whether the seed end was right.

**Predictions, written before any run** (`mid` per-fallback cost 333,715, `grad` 340,128).

1. `S11n` on `mid` keeps 4 to 8 failed seeds; at most 1 accepted end is over 1 degree.
2. `S11n` on `mid` needs 8% to 18% fewer evaluations than `S11` (44% to 50% fewer than `S8`)
   and passes both rules.
3. `S11h` on `mid` keeps at most 3 failed seeds and saves under half of `S11n`'s saving, so
   "plain does as well" is not met.
4. `S11n` on `grad` keeps 6 to 10 failed seeds and needs 15% to 28% fewer evaluations than
   `S11`; it passes both rules.
5. On both samples the median kept seed grows at most 3 voxels, seed included.
6. Wrong voxels after the pass of each `S11n` run are within 5 of `S11`'s.
7. Of the refused seeds whose end was right, the grid fallback moves the end by more than 1
   degree in at most 1 case per sample.

**Cost.** At most 56.8 million evaluations (each run as costly as `S11`), about 50 million
expected. Stored `S11` fits ran at 19 to 22 million an hour per process (`grad` 8,622,191 in
1,485 s, `mid` 13,034,526 in 2,144 s): about 3 process-hours, about 1.5 hours with four runs side
by side. Cap 240 minutes a run (at the cap: unfinished). At most 4 heavy processes. Wall times
are contended and are not timing results.

**Stop rule:** stop and report if the `held_keep` guard refuses the current library (step 5
edited `reconstructor.py`), a phase A check differs, or a validity assertion fails. Nothing is
tuned to pass a check; no arm reruns with another random seed.

**Phase A: prepare and check, stop at "READY TO RUN".**
1. Merge the parent; confirm `held_keep.install` accepts the current `_fit_from_seed`.
2. Link the stored `S8` and `S11` runs and passes of `mid` and `grad`; check their totals
   against `findings_ar5_mix.md` (26,443,074; 18,000,230; 16,082,696; 12,307,994).
3. Tests first: `tests/test_ar5_ver_rules.py` (the `h_c` pick on a synthetic pool, the plain
   accept, the validity assertion, the growth-line parser), committed before the scripts.
4. `ver_pick.py` computes `h_c`; `frozen.json` is committed before any run on `mid` or `grad`.
5. Checks: 20 seed steps of `mid` give the stored `S11` ends and counts; the in-run `N1` score
   equals the offline score on `mid`'s 15 stored cases to 1e-6; two forced pilots of `S11h` on
   `big` voxels (accept path kept and grown; refuse path runs the grid). No pilot on `grad`.
6. `ver_run_all.sh --dry-run` passes its guards. Write "READY TO RUN" and stop.

The main session starts `ver_run_all.sh` (four runs, each then its pass).

**Phase B.** Validity, `ver_stats.py`, `ver_tables.py`; the plot (evaluations by stage and wrong
voxels per run and sample) before the tables; predictions, pareto rows, page, history entry.
**Files:** `scripts/autoresearch/ar5/ver_*.py` and `ver_run_all.sh` (`held_*.py` unchanged),
`tests/test_ar5_ver_rules.py`, `benchmarks/autoresearch/ar5/ver/` (under 5 MB tracked), page
`docs/autoresearch/findings_ar5_ver.md`. **Limit:** about 12 failed seeds per sample (63 on
`held`), so counts are small. **What follows:** `S11n` passes on `mid` and plain does not do as
well: a decision on an opt-in key; plain does as well: any key is the plain cut; else no key.

**Not in this step.** Realistic images. Another material or crystal system. A new sample. Any
tuning of `N1`, its cut or `h_c` on `mid` or `grad`. `N1` as a replacement for the seed test, or
as the flag of unfinished voxels. A library key. A timing. A second run of any arm.

**Choices made.**
- Test sample `mid`, with `grad` labelled in-sample. Alternatives: `mid` alone (cheaper, about
  12 cases); a new sample as a different test condition (for example grains with gradients
  above 4 degrees), which would also need `S8` and `S11` runs and its own decision.
- `big` not run (no grid fallback). Control `h_c` with no wrong end in the pool; alternatives: a
  cut matching `N1`'s accept rate on the pool; the C++ convergence rule (equals `S11` here).
- The step-4 script layer, not the step-5 keys plus a hook (key parity not shown on `mid`).
- "Plain does as well" at 80% of the saving. Alternative: any saving within one grid fallback.

### AR5 step 7: skip the boundary pass at voxels that already fit

Task `ar5g`. Approved: D-20261112-N4. Starts from `feature/eval-reduction`; it links the
stored passes of steps 3, 4 and 6 from their kept worktrees and does not need steps 4 and 6
merged. Its scripts read those stored files themselves and import nothing from the two
unmerged task branches. One implementer session runs phases A, B and C in order and stops
only on the stop rule; no run needs a quiet machine.

**Why.** After `N1`, the boundary pass is the largest stage of `S11n` on `grad` and the second
on `mid`. It tries every neighbour orientation at every boundary voxel, also where the voxel's
end already lights nearly all its predicted pixels. Step 1 tested skips that read each try's
start; they failed, because winning tries start low. This step reads only the voxel's own end
before any try: "can a nearby start still help this voxel?" That is goal 5, as a plain rule.

**Where the evaluations go, by stage** (written with the draft of this plan, 2026-11-12).
Files read (all stored, nothing rerun): step 6 `ar5/ver/summary_mid.json`, `summary_grad.json`,
`tables.md` (fixes worktree `agent-a9420a98460078e81`); the run logs `mix_log.json` of `S11n`
(`ar5/ver/cache/`, worktree `agent-a89458addd73332d9`) and of `S11` (`ar5/mix/cache/`, worktree
`agent-a5a3601b746680640`). Clean images, Q-max 8, every total with the boundary pass.

Reading line: one row per run and sample; each cell is millions of evaluations, then in brackets
the share of that run's total in percent. "MC part" is the part of the seed step spent in
refinements where the MC optimizer ran after the CMA fit ended under hit ratio 0.95 (their CMA
part included); `S8` has no CMA fit there. "Nbr" is the neighbour local optimizations.

| sample | run | seed step | MC part | grid fallback | nbr | revisit | pass | total |
|---|---|---|---|---|---|---|---|---|
| mid | S8 | 12.01(45.4) | n/a | 5.00(18.9) | 3.86(14.6) | 2.52(9.5) | 3.06(11.6) | 26.44 |
| mid | S11 | 5.83(36.2) | 5.57(34.6) | 4.00(24.9) | 1.49(9.3) | 1.72(10.7) | 3.05(19) | 16.08 |
| mid | S11n | 5.52(39.2) | 5.27(37.3) | 2.33(16.6) | 1.49(10.6) | 1.70(12) | 3.06(21.7) | 14.11 |
| grad | S8 | 3.01(16.7) | n/a | 5.45(30.3) | 2.80(15.5) | 3.06(17) | 3.68(20.5) | 18.00 |
| grad | S11 | 1.63(13.2) | 1.38(11.2) | 3.74(30.4) | 1.48(12) | 1.78(14.4) | 3.69(29.9) | 12.31 |
| grad | S11n | 1.57(16.5) | 1.32(13.9) | 1.02(10.7) | 1.48(15.5) | 1.78(18.7) | 3.68(38.6) | 9.54 |

What the rest of the seed step is: costing 20 candidates per seed step (33,040 on `mid`), the
CMA-only refinements (223,110) and the `N1` score (26). Grid fallbacks: 7 on `mid` and 3 on
`grad` in `S11n`, about 334,000 and 340,000 each. The pass has 50,151 and 60,408 tries of exactly
61 evaluations, over 17,877 and 23,108 boundary voxels; it fixed 43 and broke 3 voxels on `mid`.

So two pieces are left. On `mid` the MC part of the seed step is largest (37.3%). On `grad`
the boundary pass is (38.6%), and it is second on `mid` (21.7%). The pass is the same size in
every run (3.0 to 3.7 million per map), so its share grows as the other stages shrink.

A what-if from the stored tries (not a result; computed while writing this draft). A voxel's
pass tries are skipped when the hit ratio of its own end before the pass is at least 0.95. The
replay is exact for the stored tries: a skipped voxel keeps its end from before the pass.
Reading line: per stored pass, its evaluations, the share the skip removes, and wrong voxels.

| sample, run | pass evaluations | cut | wrong, full | wrong, skip | fixed / broken | median error |
|---|---|---|---|---|---|---|
| mid, S11n | 3,059,211 | 68.9% | 35 | 35 | 0 / 0 | 0.02913 to 0.02913 |
| grad, S11n | 3,684,888 | 83.3% | 59 | 59 | 0 / 0 | 0.02118 to 0.02122 |
| big, S11 | 1,135,637 | 75.9% | 18 | 18 | 0 / 0 | 0.03077 to 0.03084 |
| grad, S11 | 3,685,803 | 83.5% | 62 | 62 | 0 / 0 | 0.02104 to 0.02105 |
| mid, S11 | 3,048,170 | 69.2% | 36 | 36 | 0 / 0 | 0.02913 to 0.02913 |

The highest hit ratio before the pass of any wrong voxel is 0.944 (`mid`, `big`) and 0.936
(`grad`), so the margin under 0.95 is small. At 0.92 the skip breaks 6 voxels on `mid`.

**Terms.** **Boundary pass, try**: as in `icenine/boundary_pass.py`; a try is one CMA local
optimization of 61 evaluations from one neighbour's orientation. **Own hit ratio**: the hit
ratio of the voxel's end before the pass; the pass computes it already (one uncounted overlap
evaluation per voxel, 24,570 on `mid`). **Skip cut** `h_s`: a voxel with own hit ratio at least
`h_s` gets no try. **Net gain**: wrong voxels before the pass minus after it. **Replay**: the
stored tries with the skipped ones removed. **Live pass**: the library pass run with the key.
**Wrong**: more than 1 degree from the truth (the library's symmetry-reduced angle).

**Samples, and what has touched them.** `big`, `grad`, `mid`: this draft computed the skip at
0.95 on their stored passes (table 2 below); disclosed, not clean. `held` (250 grains, gradients
0 to 4 degrees): its `S8`, `S11`, `S11v2` passes and wrong counts were read in step 4, but no
skip has been computed on it. `held` is the test sample. `S11v2` on `held` is `S11n`'s rule.

**Frozen now.** `h_s` = 0.95, the library's existing stop ratio (`SeedIndexStopHitRatio`) and
step 1's flag pick on `big` and `grad`. No other value is tried as a decision.

**What is measured**, per sample and run: pass evaluations with and without the skip, the cut,
voxels and tries skipped; wrong voxels after the pass (Wilson 95%); fixed and broken against the
full pass and against `S8`; median error; run totals with the skipped pass, against `S8`; the own
hit ratio of every voxel the full pass fixed. A curve over `h_s` from 0.90 to 0.99 is a what-if.

**Predictions, written before any run.**
1. On `held` the replay at 0.95 cuts 55% to 85% of the pass and fixed / broken against the full
   pass is at most 0 / 2.
2. On `held`, `S11v2` with the skipped pass needs 5% to 12% fewer evaluations than with the full
   pass (its pass is a smaller share: an estimate of 3.5 of about 36 million, from step 4's
   totals; this draft did not open `held`'s pass files).
3. The live pass spends exactly the replay's evaluations (61 per kept try) on every sample.
4. Live wrong voxels are within 2 of the replay's on every sample (one generator for the pass,
   so skipped tries shift the random draws of later tries).
5. Median error rises by at most 0.0005 degree on every sample.

**Decision rules, fixed now.** The skip passes on a sample if: the pass is cut by at least 50%;
broken minus fixed against the full pass is at most 5% of the pass's net gain, rounded up, at
least 1 voxel; median error is at most 0.001 degree higher; and the run with the skipped pass
meets the house rules against `S8` (20% fewer, broken minus fixed at most 5, median at most
0.005 degree higher). It passes the step if it passes on `held` live and on every disclosed
sample. Then the owner decides whether the documentation recommends the key. If it fails on
`held`: the key stays, off, the page says so; `h_s` is not tuned on `held`.

**Machine.** Phase A reads stored files only, any time. Phase B runs passes alone on stored maps
(no BFS): four skipped passes of 0.27 to 1.0 million evaluations and one key-off parity pass of
3.06 million, under 10 minutes each, at most 2 side by side. Wall times are not timing results.

**Stop rule:** stop and report if a stored total does not reproduce, if the replay at 0.95 fails
the rule on a disclosed sample (before any library change), if live evaluations differ from the
replay, or if the key-off pass differs from the stored one in any try.

**Phase A (offline).** 1. Link the stored passes and maps; check the totals 3,059,211,
3,684,888, 1,135,637, 3,685,803, 3,048,170. 2. Tests first: `tests/test_ar5_gate.py` (replay
rule on a synthetic pass, the fixed/broken count). 3. `gate_replay.py` on the disclosed passes;
commit `ar5/gate/frozen.json` (`h_s`, rule, file hashes). 4. Then `held`'s passes; per-voxel
errors from the maps and the truth, no evaluation. If `held` has no stored `old_score`, it is
tested live only. Write "READY FOR KEY".
**Phase B (key and live passes).** Add the key; unit tests; key-off parity pass on `mid` `S11n`;
live skipped passes on `mid`, `grad` (`S11n`), `big` (`S11`), `held` (`S11v2`).
**Phase C.** Stats, the plot (cut against net gain lost, one line per sample, 0.95 marked)
before the tables, predictions, page `findings_ar5_gate.md`, history entry, pareto rows.
**Files:** `scripts/autoresearch/ar5/gate_*.py`, `tests/test_ar5_gate.py`,
`benchmarks/autoresearch/ar5/gate/` (under 2 MB tracked), the page.

**Library changes:** one opt-in key, `BoundaryPassSkipHitRatio h` (float; default 0, off). In
`icenine/boundary_pass.py` a boundary voxel whose own hit ratio is at least `h` gets no try; its
score is still recorded; `bfs.stats` gains `boundary_skipped`. `icenine/config_file.py` parses
the key; `reconstructor.py` passes it with `BoundaryPass`. The off path is bit-identical; no
default changes and no golden regeneration.

**Not in this step.** Realistic images (lower hit ratios: more tries run, less saving). Another
crystal system. `N1` as the skip score. The skip stacked with step 1's `P2 m=2`. Any change
of the seed step or the grid fallback. A timing. A second run of any arm.

**Choices made.**
- This step is candidate (b) aimed at the pass: a saving already seen offline on five stored
  passes with no voxel lost, so low risk. Not (a), the `N1` key: it adds no saving beyond the
  12.3% measured; it makes that usable, moving the feature code of `held_keep.py` and `ver_*.py`
  into the library with a parity test on the stored `S11n` runs. Other files; it can run beside.
- Not (b) aimed at the MC part of the seed step (37.3% on `mid`, 13.9% on `grad`): the log keeps
  only the CMA end, not the MC end or trace, so a first phase needs a replay of about 5 million
  evaluations, and steps 2 and 3 already tried replacing that MC. It is the next target.
- Not (c) `N1` as the stop test of a grown voxel's local optimization: those fits are short (about
  660 evaluations in the pilot), no per-evaluation traces are stored, and `N1` costs 2
  evaluations a call.
- `h_s` = 0.95, not 0.97 (cut 53.7% on `mid`, safer margin) or a cut fit on the pass data.
- One generator for the pass kept; alternative: a generator per voxel, which would make the live
  pass equal the replay but would change the key-off tries.

### AR5 step 8: fit every candidate with CMA before any MC refinement

Task `ar5h`. Approved: D-20261112-N7. Starts from `feature/eval-reduction`; links the
stored runs of steps 2, 3 and 6 from their kept worktrees, imports nothing from unmerged branches.
Phase A stops at "READY FOR KEY"; phase B starts only after a recorded decision.

**Arms (amended, D-20261112-N8).** Both arms are the run `S11` (the merged keys
`SeedIndexCMAFirst 200` and `BFSNeighborCMANoFallback 1`, script layer of step 3), without the
`N1` seed-accept rule. That rule's code is on the unmerged branches of steps 4 and 6, and it acts
only on failed seeds (13 on `mid`), after the seed step, so it does not touch the order tested
here. The step then needs no unmerged code. Where the text below says `S11n` for the arms, read
`S11`: step 1 checks the stored `S11` totals of step 3 instead of 14,105,137 and 9,539,809, and
the runs are a little larger (about 16.1 and 12.3 million on `mid` and `grad`). The predictions
and the decision rules are unchanged. The sweep order with `N1` is a later combined run.

**Why.** After step 7 the largest stage left in `S11n` is the MC part of the seed step: the seed
refinements where the stock MC optimizer runs because the CMA end stayed under hit ratio 0.95.
The mixed seed refinement does this candidate by candidate: CMA on candidate 1, then MC on
candidate 1 if needed, then CMA on candidate 2, and so on. Often the MC end of candidate 1 is
thrown away, because a later candidate's plain CMA end reaches 0.95. This step changes only the
order: CMA on every candidate first, MC only if no CMA end reached 0.95. It is the owner's goal 5
as a plain rule: run the cheap fit everywhere before deciding that the costly one is needed.

**What the stored logs say about each MC refinement** (computed while drafting; files read, nothing
rerun: `mix_log.json`, `index_log.json.gz`, `e1_log.json.gz`, `n1_log.json` of `S11n` in
`ar5/ver/cache/`, worktree `agent-a89458addd73332d9`; voxel errors from `errors_<sample>.npz`).
Stored per refinement: the CMA end hit ratio, whether MC ran, evaluations, the hit ratio and cost
of the kept end, the candidate's place in the order, and whether the step took it. Per failed seed:
the `N1` score, accepted or not, the seed end and its error. Not stored: the MC end apart from the
kept one, any trace, or an error per refinement. Reading line: one row per count; "steps" are
seed steps; evaluations include the 201 of the CMA fit. (The `S11` run folder has no index log.)

| per MC refinement, run `S11n` | mid | grad |
|---|---|---|
| MC refinements (seed steps with one) | 1,653 (811 of 1,652) | 996 (491 of 1,343) |
| evaluations, total (share of the run with the pass) | 5,266,711 (37.3%) | 1,322,010 (13.9%) |
| evaluations each: mean, median, largest | 3,186; 1,973; 33,499 | 1,327; 1,126; 9,893 |
| MC end not below the CMA end in cost (MC gained nothing) | 481 | 520 |
| end reached 0.95 and stopped the step: refinements; evaluations | 406; 2,352,519 | 132; 297,908 |
| taken as the lowest cost, under 0.95 | 136; 387,459 | 94; 174,071 |
| discarded, a later candidate was taken | 1,111; 2,526,733 | 770; 850,031 |
| step accepted, voxel right before the pass: steps; refinements; evaluations | 786; 1,528; 5,086,143 | 465; 866; 1,175,685 |
| step accepted, voxel wrong before the pass | 12; 60; 89,588 | 17; 85; 96,813 |
| step failed the seed test, seed kept by `N1` | 6; 30; 45,948 | 6; 30; 30,683 |
| step failed the seed test, grid fallback | 7; 35; 45,032 | 3; 15; 18,829 |
| MC-only evaluations of steps that a later CMA end stopped | 832,456 (269 steps) | 350,959 (265 steps) |

"Voxel right" uses the voxel's error before the pass, which is the seed end unless a later revisit
changed it. Most of the MC part (96.6% and 88.9%) is in steps that end with a right seed; 48.0%
and 64.3% of it buys an end that the step then discards.

**What steps 2 and 3 tried.** Step 2 replaced the MC refinement by CMA alone (`C60`, `C200`, `C0`);
`C200` lost 60 seed-test passes in 1,608 steps (limit 16), so MC is sometimes needed. Step 3 kept
MC but ran it only after that candidate's CMA end missed 0.95 (`X200`). This step keeps both the
MC and the 0.95 test and changes only when MC runs. No new threshold, budget or model.

**A what-if on step 2's stored replay** (computed while drafting; arithmetic, not a run). The replay
(`ar5/mix/cache/replay_*.json.gz`, worktree `agent-a5a3601b746680640`) holds, for every `S8` seed
step, the `C200` end and the MC end of each candidate up to where each arm stopped. Every refinement
starts a fresh generator (`default_rng(0)`), so its end depends only on its start, and the replay
gives both orders exactly. No step needed an end the replay lacked.
Reading line: one row per sample; seed-step evaluations of the two orders; then counts.

| sample | seed steps | `X200` | sweep order | cut | passes, both | right, both | paired split | taken end differs | median error |
|---|---|---|---|---|---|---|---|---|---|
| big | 329 | 938,875 | 230,604 | 75.4% | 329 | 329 | 0 / 0 | 37 | 0.01232 to 0.01232 |
| grad | 1,279 | 1,563,116 | 1,301,274 | 16.8% | 1,267 | 1,260 | 0 / 0 | 7 | 0.00898 to 0.00906 |
| mid | 1,611 | 5,643,575 | 3,268,675 | 42.1% | 1,600 | 1,599 | 0 / 0 | 188 | 0.01135 to 0.01121 |

(`X200` on `mid` is 201 above step 3's 5,643,374: the one-step variant explained in `mix_part1.py`.)
On `mid` the saving splits into steps stopped later by a CMA end (888,498) and steps where `X200`
stopped on an MC end but a later CMA end also reaches 0.95 (188 steps, 1,616,248), less extra CMA
fits where MC is still needed (205 steps, 129,846). On `grad` only 7 steps are of the second kind.

**Terms.** **Candidate**: one of the indexer's orientations, refined in the order of strict cost,
at most 5. **CMA end**: the end of the 200-evaluation CMA fit. **MC refinement**: the stock MC
optimizer from the same candidate. **Interleaved order** (`X200`, today's `SeedIndexCMAFirst 200`):
per candidate CMA, then MC if the CMA end is under 0.95. **Sweep order** (`Y200`): CMA on candidates
in order until a CMA end reaches 0.95; if none does, MC on candidates in order, each candidate's end
the lower cost of its two, stopping at the first end that reaches 0.95; else the lowest cost of all
ends. **Seed test**: hit ratio at least 0.80. **Wrong**: over 1 degree from the truth (the library's
symmetry-reduced angle). **Aligned draws**: every local optimization, grid fallback and pass try in
the layer draws from its own generator, seeded by the voxel, the kind of call and its count there.

**Samples, and what has touched them.** `big`, `grad`, `mid`: the what-if above was computed on
them while drafting; disclosed, not clean. `grad` was also `N1`'s fitting sample (that part is
unchanged). `held` (250 grains, gradients 0 to 4 degrees): its whole-layer runs were read in steps
4, 6 and 7, but no seed-order count was computed on it and this draft opened none of its seed logs.
The rule has no value fit to any sample. `held` is the test sample; `S11v2` there is `S11n`'s rule.

**Frozen now.** The sweep order as defined; 0.95, 200, 5 candidates and the `N1` cut unchanged.
**What is measured**, per sample and arm (with the pass): evaluations by stage; seed-step
evaluations and MC refinements; grid fallbacks; seeds kept by `N1`; wrong voxels (Wilson 95%);
fixed and broken by voxel against the other arm and against `S8`; median error; the seed steps
whose taken end differs, and how many of those are wrong. Each MC refinement's trace (best cost
after each improvement, with the evaluation count) is written to the cache and not read here.

**Predictions, written before any run.** Arms `A` (interleaved) and `Y` (sweep), aligned draws.
1. `held`: `Y` needs 20% to 45% fewer seed-step evaluations and 6% to 16% fewer in total than `A`.
2. `mid`: 30% to 45% fewer in the seed step, 11% to 18% in total; `grad`: 10% to 25% and 1% to
   5%; `big`: 60% to 80% and 12% to 22%.
3. On every sample, grid fallbacks of `Y` are within 2 of `A`'s, and broken minus fixed against
   `A` is between -3 and 3 voxels.
4. Median error moves by at most 0.0005 degree on every sample.
5. Where both arms refine the same candidate at the same voxel, the ends are identical.

**Decision rules, fixed now.** `Y` passes on a sample if, against `A` on that sample: (1) seed-step
evaluations at least 10% fewer and the total fewer; (2) broken minus fixed at most 5; (3) median
error at most 0.001 degree higher; (4) at most 2 more grid fallbacks; and (5) the house rule against
`S8` with the pass: at least 20% fewer evaluations, broken minus fixed at most 5, median error at
most 0.005 degree higher. If `A` itself fails (5) on a sample, the page says so and `Y` is judged
there by (1) to (4). The step passes if `Y` passes on `held` and on all three disclosed samples.
Then phase B adds the key, and the owner decides on any recommendation. If it fails on `held`: no
key, the page says so, nothing is tuned.

**Machine.** Phase A: checks under 10 million evaluations; eight whole-layer runs (`A`, `Y` on four
samples) of about 37.5 + 33 (`held`), 14 + 12 (`mid`), 9.5 + 9.3 (`grad`), 4 + 3.3 (`big`)
million, and one run `Yl` on `big` (3.3 million): about 126 million. At 19 to 22 million an hour
per process, about 6.5 process-hours, about 2.5 hours wall with 4 processes side by side. Cap 300
minutes for `held`, 240 for the rest; at most 4 heavy processes; wall times are not timing results.

**Stop rule:** stop and report if a stored total or a stored refinement end does not reproduce, if
two runs of `A` on `big` differ, if a validity check fails (a step with a CMA end at 0.95 ran MC in
`Y`; a refinement shared by both arms differs), or if the key-on library run differs from `Yl`.

**Phase A (script layer only).** 1. Link the stored replay and `S11n` runs; check the totals
14,105,137 and 9,539,809 and the replay counts above. 2. Tests first: `tests/test_ar5_order.py`
(the sweep rule on synthetic steps, the aligned-draw seeding). 3. `order_whatif.py` writes the
what-if table; commit `ar5/order/frozen.json` (rule, file hashes). 4. Checks: 20 stored seed steps
of `mid` give the stored per-candidate ends in both orders; `A` twice on `big` is identical.
5. The eight runs, `held` first, then `Yl` (sweep order with the library's own draws). Validity,
stats, the plot (seed-step evaluations against wrong voxels, one point per arm and sample) before
the tables, predictions. Write "READY FOR KEY" and stop.
**Phase B (key).** Add the key; unit tests; key-off suite with no golden change and the library
`S11` on `big` at 2,879,705 before the pass (step 5, check 3); key on, library run on `big` equal to
`Yl` byte for byte. Page `findings_ar5_order.md`, history entry, pareto rows.
**Files:** `scripts/autoresearch/ar5/order_*.py`, `tests/test_ar5_order.py`,
`benchmarks/autoresearch/ar5/order/` (under 3 MB tracked; traces untracked), the page.

**Library changes:** phase A none. Phase B, only if the step passes: one opt-in key,
`SeedIndexCMASweep 0|1` (default 0, off), refused unless `SeedIndexCMAFirst` is 2 or more. With 1,
`ray_indexing.seed_by_indexing` runs the sweep order with the existing CMA and MC functions of
`reconstructor.index_seed`; `config_file.py` and `orientation_search.py` parse and check it. The
off path is bit-identical; no default changes, no golden regeneration. The aligned draws stay in
the script layer as a test device; the library keeps its own draws.

**Not in this step.** An early stop inside the MC run (traces only recorded). A skip of MC from the
CMA end's hit ratio. Other values of 0.95, 200 or 5. Realistic images. Another crystal system. A
timing. A second run of any arm with another random seed.

**Choices made.**
- Sweep order, not (a) an early stop read from the MC run's own trace: that needs a replay of about
  5 million evaluations, a fitted predictor and a sample to fit it on; it is the natural next step
  on what is left, and phase A records the traces for it at no extra evaluation.
- Not (b) skipping MC when the CMA end is low: on `mid`, CMA ends under 0.5 hold 462 MC runs
  (772,418 evaluations, 5.5% of `S11n`), of which 37 still reached 0.95; it needs a new cut.
- Test sample `held`, not a new sample from `scripts/autoresearch/ar3/held_make_sample.py`
  (another generator seed): `held` saw no count of this rule. A new sample would also need `S8`.
- Both arms rerun with aligned draws (step 7's lesson), not `Y` against the stored `S11n`.

### AR5 step 9: rescue a failed seed with a wide local optimization before the grid search

Task `ar5i`. Approved: D-20261112-N10. Starts from `feature/eval-reduction` after `ar5h` merges
as a study (D-20261112-N9); it reuses `order_run.py`, `order_pass.py`, `order_draws.py`. Phase A
ends at "READY FOR KEY" or at its gate; phase B starts only after a recorded decision.

**Why.** When a seed fails the seed test (hit ratio 0.80), the run does a grid fallback: a global
search on the adaptive grid, about 335,000 evaluations each. On `held` it is the largest stage
(`S8`: 75 fallbacks, 43% of its evaluations; `S11`: about 63 to 65, about 21 of an estimated 37.5
million). The owner asks for a cheap approximation of that search, with a final local optimization
after it. The stored ends say the grid search mostly lands next to an end the seed step already had.

**What was read** (stored files only, nothing rerun). Arm `A` of step 8 (run `S11`, aligned
draws), `ar5/order/cache/A_mid`, `A_grad`: `e1_log.json.gz`, `records.json.gz`, `recon.mic.gz`,
`stats/part2_A_<sample>_err.npz`; `ar5/order/stats_<sample>.json` for big, grad, mid; the truth
(`gradient_truth.mic`, `sample.json`) of `mid` and `grad`; step 1's page `findings_ar5_stop.md`
(Part 2: the 31 `S8` fallbacks) and AR0's `findings_ar0_stopping.md` (the split by level). The
stored `S11` runs of step 3 were not read (arm `A` is the same run with aligned draws). No `held`
file was opened or listed. `big` has no fallback in arm `A`.

**What a fallback's log holds.** `records.json.gz`: one line per voxel (evaluations, hit ratio of
the end). `e1_log.json.gz`: the failed seed step's 10 fits (CMA and MC on 5 candidates: hit ratio,
cost, evaluations, end orientation). The map gives the fallback's end, the truth its error. Not
stored: anything per level (candidates, costs, evaluations), the 20 indexed candidates, the cost
of the fallback's end. AR0 measured the split by level on six full-sample voxels: level 0 screen
23.7%, level 0 quick MC 37.1%, level 1 26.3%, level 2 8.8%, level 3 2.8%, final local
optimization 1.3%. A per-level record of the 23 disclosed fallbacks costs about 7.9 million
evaluations (about 25 process-minutes); phase A step 5 makes it, as data for a later learned step.

Reading line: one column per sample; each row counts arm `A`'s fallbacks (failed seed steps).

| arm `A` fallbacks | grad | mid |
|---|---|---|
| fallbacks; evaluations (share of the run with the pass) | 11; 3,726,723 (30.6%) | 12; 4,007,241 (25.6%) |
| a kept seed end within 0.6 degree of the fallback's end | 11 | 11 |
| the first such end's place by cost: 1st / 2nd / 3rd | 9 / 1 / 1 | 8 / 2 / 1 |
| nearest kept end to the fallback's end, degrees | 0.01 to 0.58 | 0.08 to 0.56 |
| lowest-cost kept end is wrong (26 to 55 degrees off) | 2 | 4 |
| fallback's end right (at most 1 degree); hit ratio at least 0.95 | 11; 10 | 12; 11 |
| lowest-cost kept end: hit ratio; cost | 0.66 to 0.80; 0.25 to 0.38 | 0.71 to 0.80; 0.21 to 0.32 |

The one case with no near end is `mid` voxel 16979 (nearest 15.0 degrees). So in 22 of 23 cases
(Wilson 95% 79% to 99%) the grid search ends within 0.6 degree of a kept seed end, and its end is
right in 23 of 23 (86% to 100%), with errors 0.001 to 0.09 degree (one 0.55).
Ends with hit ratio at least 0.95 cost 0.072 (median, `mid`), the failed ends 0.21 to 0.38. That
is consistent with the seed refinements being trapped 0.07 to 0.58 degree from a much better
orientation, not with a missing candidate (step 1 found the same on `S8`: the right candidate was
1st to 3rd in all 31). All 12 `mid` fallbacks are in grains with gradients of 2 or 4 degrees (34%
of its voxels); all 11 on `grad` at its top level, 1 degree (23%). `held` has gradients 0 to 4.

**Cases.** (a) Near an end the seed step had: 22 of 23. (b) Already best after level 0: not
stored; AR0 found the level-0 best within 1 degree for only 71 of 159 right voxels, and a stop
after level 0 saves at most 39.2% of a search. (c) No good answer: 2 of 23 fallback ends stay under
0.95 (0.70, 0.87), both right; none is off the sample.

**What-ifs** (arithmetic; a replay is a cost estimate only, step 8's lesson). Reading line: per
idea, the evaluations it would save as a share of arm `A`'s run, and what blocks it.

| idea | grad | mid | what blocks it |
|---|---|---|---|
| `G1`: wide CMA from the kept ends, take a 0.95 end, else grid | 27.7% | 21.2% | untested whether CMA escapes the trap |
| stop the grid after level 0, finish the best 5 with CMA | at most 12% | at most 10% | level-0 best often over 1 degree off |
| a wider indexed list before the grid (step 1's `R1`) | 0 | 0 | passed the seed test in 0 of 31 |
| skip the fallback, leave the voxel to growth and the pass | up to 30.6% | up to 25.6% | ends were right in 23 of 23; no data on growth |

`G1`'s numbers assume it reaches 0.95 exactly where a kept end is within 0.6 degree of a grid end
that reached 0.95 (10 of 11 on `grad`, 10 of 12 on `mid`), at 1,000 evaluations a try (24 and 18
tries). They are an upper estimate until the replay in phase A.

**Terms.** **Kept end**: per candidate, the lower-cost end of its CMA fit and its MC refinement.
**Rescue**: the CMA local optimization (the library's `CMAOptimizer`) from one kept end, initial
step 0.25 degree, 1,000 evaluations, its own generator (`default_rng(0)`). **Rule `G1`**: after a
failed seed test, rescue the kept ends in order of cost and take the first end with hit ratio at
least 0.95; if none, run the grid fallback as today. **Arm `G`**: arm `A` with `G1`. **Wrong**:
over 1 degree from the truth (the library's symmetry-reduced angle). **Aligned draws**: as step 8.

**Samples, and what has touched them.** `grad`, `mid`: this draft read their fallback cases and
chose the values on them; disclosed. `big`: no fallback, so `G` must equal `A` there exactly.
`held` (250 grains, gradients 0 to 4 degrees): totals and fallback counts read in steps 4 and 8;
none of its seed logs or fallback ends opened. `held` is the test sample. Clean images, Q-max 8.

**Frozen now.** `G1` with 0.25 degree, 1,000 evaluations, at most 5 tries, accept at 0.95 (the
library's `SeedIndexStopHitRatio`). Chosen from the table: the farthest near end is 0.58 degree,
about 2.3 initial steps; 1,000 is five times the seed step's CMA budget. Not tuned after any run.

**What is measured**, per sample: evaluations by stage; failed seed steps, rescues tried, rescued,
grid fallbacks left; error and hit ratio of each rescued end and of the grid end it replaced (in
the replay); wrong rescued seeds; seed steps whose taken end differs; wrong voxels (Wilson 95%);
fixed and broken against `A` and `S8`; median error.

**Predictions, written before any run.**
1. Cost check: in at least 21 of 23 cases the grid's end costs less than every kept end.
2. Replay (23 disclosed cases): `G1` reaches 0.95 in 12 to 21; every rescued end is right; its
   median error is at most 0.05 degree; at most 2.5 tries a case on average.
3. `held`: grid fallbacks 40% to 85% fewer than `A`'s; total 20% to 45% fewer.
4. `mid` 12% to 21% fewer in total; `grad` 15% to 28%; `big` identical to `A`.
5. Broken minus fixed against `A` between -3 and 5 voxels; median error moves at most 0.001 degree.

**Decision rules, fixed now.** Gate before any whole-layer run: the replay must reach 0.95 in at
least 12 of 23 cases with no wrong rescued end; else stop, report, no tuning. Then `G` passes on a
sample if, against `A`: (1) on `held` and `mid` fewer grid fallbacks and the total at least 10%
fewer, on `grad` the total fewer, on `big` identical; (2) broken minus fixed at most 5; (3) median
error at most 0.001 degree higher; (4) at most 1 wrong rescued seed; (5) the house rule against
`S8` with the pass: at least 20% fewer evaluations, broken minus fixed at most 5, median error at
most 0.005 degree higher. If `A` fails (5) on a sample, `G` is judged there by (1) to (4). The step
passes if `G` passes on `held` and on all disclosed samples. If it fails on `held`: no key.

**Machine.** Cost check: 46 evaluations. Replay: at most 115,000. Per-level record: about 7.9
million. Arm `G` whole-layer: about 18 to 30 (`held`), 12.3 (`mid`), 8.8 (`grad`), 3.8 (`big`)
million. About 51 to 63 million, about 3 process-hours: about 1.5 to 2 hours wall with 4
processes side by side. Caps: 300 minutes `held`, 240 the rest; at most 4 heavy processes.
Arm `A` is the stored step 8 run (`G` on `big` doubles as its parity check).

**Stop rule:** stop and report if a stored total does not reproduce (4,007,241 and 3,726,723
fallback evaluations; 15,645,472, 12,189,777, 3,807,747 totals), if the gate fails, if `G` on
`big` differs from `A_big`, if a rescue on a disclosed sample ends wrong, or if the key-on library
run differs from the script arm with library draws.

**Phase A (script layer).** 1. Link step 8's `A` runs; check totals; the cost check. 2. Tests
first: `tests/test_ar5_rescue.py` (rule on synthetic steps; aligned draws unchanged). 3. Commit
`ar5/rescue/frozen.json` (rule, values, hashes). 4. Replay `G1` on the 23 disclosed cases; gate.
5. Per-level record of the 23 fallbacks (beside step 6). 6. `G` on `held` first, then `mid`,
`grad`, `big`; passes; validity, stats, the plot (run evaluations against wrong voxels) before
the tables; predictions. "READY FOR KEY".
**Phase B (key).** The key; unit tests; key-off suite with no golden change; key on, a library
run on `mid` equal to a script arm `Gl` (library draws). Page `findings_ar5_rescue.md`, history.
**Files:** `scripts/autoresearch/ar5/rescue_*.py`, `tests/test_ar5_rescue.py`,
`benchmarks/autoresearch/ar5/rescue/` (under 3 MB tracked; records untracked), the page.

**Library changes:** phase A none. Phase B, only if the step passes: one opt-in key,
`SeedRescueCMA sigma evals` (default off), used by `reconstructor.index_seed` after a failed seed
test and before the grid fallback; `config_file.py`, `orientation_search.py` parse and check it.
Off path bit-identical; no default changes, no golden regeneration.

**Not in this step.** A learned approximation of the grid search (the per-level record is its
data). Skipping fallbacks. Other values of 0.80, 0.95. Realistic images. Another crystal system.
A timing. A second run of any arm with another random seed.

**Choices made.**
- `G1`, not a stop after level 0 (at most 39.2% of a search, and its best is often off).
- All kept ends in order of cost, not only the lowest: in 6 of 23 the lowest is wrong.
- Accept at 0.95, not 0.80: failed ends right and wrong overlap in hit ratio (step 1).
- This step, not an early stop of an MC refinement after 2,000 evaluations with no improvement
  (drafted too: it saves 12.6% on `mid`, 8.0% on `big`, 0.4% on `grad`, and arm `A` is 1 voxel
  inside the `S8` limit on `mid`). This one aims at `held`'s largest stage and has a cheap gate
  (under 0.2 million evaluations) before any whole-layer run.
- Not a library key for the seed-accept rule of `N1` now: its code is on two branches that wait
  for the owner.

### AR5 step 10: a cheaper level 0 for the grid search that is left after the rescue

Task `ar5j`. Approved: D-20261113-N3. Starts from `feature/eval-reduction` at `a5a473d` (step 9
merged, `28b75fa`); reuses `order_run.py`, `order_pass.py`, `order_draws.py`, `rescue_rules.py`.
Phase A ends at "READY FOR KEY" or at its gate; phase B starts only after a recorded decision.

**Why.** After rule `G1` the grid search still runs 19 times on `held` (6.56 million evaluations,
21.3% of arm `G`), 4 on `mid`, 5 on `grad`. The owner asks for a cheap approximation of that
search, with a final local optimization after it. Level 0 is about 61% of a grid search (AR0's
split), mostly a quick MC refinement of each of about 10,200 listed candidates. The per-level
record says a candidate near the right orientation is always near the top of that list.

**What was read** (stored files only, nothing rerun). Arm `G` of step 9: `ar5/rescue/stats_*.json`
and `tables.md` (all four samples, evaluations by stage); in the untracked cache (worktree
`agent-a8d20b3821621fd41`) `G_mid`, `G_grad`, `G_big`: `records.json.gz`, `rescue_log.json.gz`,
`mc_traces.json.gz`, `run.json`; the per-level record of arm `A`'s 23 grid searches on `mid` and
`grad` (`levels_<sample>.json`, `levels/<sample>_<voxel>.npz`); step 9's MC-stop draft. Angles
to the final end use the 24 cubic rotations, only to read the record; the rule uses no symmetry.

**Where arm `G`'s evaluations go.** Reading line: one row per sample; each stage's evaluations
as a share of the run with the boundary pass. "Seed step" includes its MC refinements.

| sample | total | seed step | rescue | grid search | neighbour | revisit | pass |
|---|---|---|---|---|---|---|---|
| held | 30,768,867 | 46.0% | 0.5% | 21.3% (19) | 4.8% | 10.4% | 16.9% |
| mid | 13,004,457 | 41.7% | 0.2% | 10.3% (4) | 11.4% | 12.9% | 23.5% |
| grad | 10,196,692 | 15.0% | 0.3% | 16.6% (5) | 14.5% | 17.5% | 36.1% |
| big | 3,807,747 | 18.1% | 0 | 0 (0) | 39.3% | 12.8% | 29.8% |

Seed-step MC refinements (arm `A`'s traces) are 37.3% (`mid`), 10.7% (`grad`), 15.7% (`big`) of
`G`'s run. A neighbour fit is a fixed 61 evaluations; a revisit is 122 to 1,220.

**What the per-level record says** (23 grid searches of arm `A`; the 9 that `G1` did not rescue
are exactly arm `G`'s grid searches on `mid` and `grad`). Reading line: one row per finding;
"place" counts from 0 in the level-0 list sorted by the discrete score.

| finding, 23 grid searches | mid (12) | grad (11) |
|---|---|---|
| evaluations each; of which discrete stages | 329,721 to 338,416; 23% | 333,524 to 342,487; 23% |
| level-0 list length | 10,092 to 10,249 | 10,300 to 10,439 |
| place of the first candidate within 5 degrees of the end: 0 to 4; largest | 11; 92 | 8; 165 |
| after the quick MC, first within 3 degrees: 0 to 4; largest | 12; 4 | 10; 88 |
| level-0 best right (at most 1 degree); wrong by 17 to 60 degrees | 5; 7 | 4; 5 (2 at 1.07) |
| level-3 best right; final end right | 12; 12 | 9; 11 |
| `G1` failed: rescue end nearest the grid end, degrees | 0.08, 0.19, 0.29; 15.0 | 0.004 to 0.58 |
| `G1` failed: hit ratio of that nearest end | 0.72, 0.77, 0.75; 0.20 | 0.59 to 0.81 |

So `G1` failed in 8 of 9 cases with an end already 0.004 to 0.58 degree from the grid's end that
stayed under 0.95; in one (`mid` 16979) no end was near. A stop after a level does not work: the
level-0 best is wrong in 14 of 23, the level-3 best in 2, and level bests cost 0.13 to 0.84
against about 0.07 for a 0.95 end. All 23 are in grains with orientation gradients.

**Terms.** **Grid search**: the global search on the adaptive grid, levels 0 to 3. Each level
has a discrete stage (screen of the start set times a 9-point local grid, strict re-cost, spacing
filter; a list sorted by the **discrete score**), a quick MC refinement of every listed candidate
(10 steps, 5 restarts), and keeps the best quarter as the next start set. Then the final local
optimization (FindOptimal on the best 30, the variance stage). **Level-0 cut**: run the quick MC
of level 0 only on the first N of its list; the quarter kept is counted from the cut list. **Arm
`C`**: arm `G` with the level-0 cut, N = 1,000. **Right**/**wrong**: at most / over 1 degree from
the truth (the library's symmetry-reduced angle). **Aligned draws**: as step 8; each grid search
keeps its generator (voxel, "grid", 0), so the cut moves no other draw in the layer.

**Samples, and what has touched them.** `mid`, `grad`: this draft read their records and chose
N on them; disclosed. `big`: no grid search, so `C` must equal `G` there exactly. `held` (250
grains, gradients 0 to 4 degrees): this draft read arm `G`'s and `A`'s stage totals and counts in
`stats_held.json` and `tables.md`, and listed the file names in `cache/G_held`; no seed log, trace,
rescue log, map or level record of `held` was opened. `held` is the test sample. Clean, Q-max 8.

**Frozen now.** N = 1,000 at level 0 only; levels 1 to 3, the final local optimization, `G1`,
0.80, 0.95 unchanged. Chosen from the record: 1,000 is about 6 times the largest place (165) and
about 10% of the list. Not tuned after any run.

**What is measured**, per sample: evaluations by stage; grid searches and evaluations each;
error, hit ratio and cost of each cut search's end and of the full search's end it replaces (in
the replay); wrong grid ends; voxels whose end changed against `G`; wrong voxels (Wilson 95%);
fixed and broken against `G` and `S8`; median error. The runs also write each neighbour and
revisit fit's best cost per generation to the cache (no evaluation added), not read here.

**Predictions, written before any run.**
1. Replay (23 disclosed cases): the cut search's end is right in 23; hit ratio at least 0.95 in
   20 or 21 (the full search: 21); within 0.1 degree of the full search's end in at least 18;
   85,000 to 135,000 evaluations a search (AR0's split: level 0's screen and re-cost stay).
2. `held`: grid evaluations 60% to 75% fewer than `G`'s; total 12% to 17% fewer.
3. `mid` 5% to 8% fewer in total; `grad` 9% to 13%; `big` identical to `G`.
4. Grid searches within 2 of `G`'s on each sample; no wrong grid end on `mid` or `grad`.
5. Broken minus fixed against `G` between -3 and 3 voxels; median error moves at most 0.001 degree.

**Decision rules, fixed now.** Gate before any whole-layer run: in the replay, the cut search
is right in at least 22 of 23, reaches 0.95 in at least 20, and the 23 need at least 50% fewer
evaluations than the full searches (7,733,964); else stop, report, no tuning. The replay's count
is a cost estimate only. Then `C` passes on a sample if, against `G`: (1) on `held`, `mid`, `grad`
grid evaluations at least 50% fewer and the total fewer, on `held` at least 8% fewer; on `big`
identical; (2) broken minus fixed at most 5; (3) median error at most 0.001 degree higher; (4) at
most 1 more wrong grid end; (5) the house rule against `S8` with the pass: at least 20% fewer
evaluations, broken minus fixed at most 5, median error at most 0.005 degree higher. The step
passes if `C` passes on `held` and on all disclosed samples. If it fails on `held`: no key.

**Machine.** Reproduction of 2 stored searches: about 0.7 million evaluations. Replay: about 2.5
million. Arm `C` whole-layer: about 26 (`held`), 12.1 (`mid`), 9.0 (`grad`), 3.8 (`big`) million.
About 54 million, 2.5 to 2.9 process-hours at 19 to 22 million an hour: about 1.5 hours wall with
4 processes side by side (`held` is the long pole). Phase B: two runs on `mid` before the pass,
about 18 million, about 30 minutes wall. Caps: 300 minutes `held`, 240 the rest; at most 4 heavy
processes. Arm `G` is the stored step 9 run; `C` on `big` doubles as its parity check.

**Stop rule:** stop and report if a stored total does not reproduce (the 2 searches: 335,077 for
`mid` 7167 and 341,304 for `grad` 19058; the record's totals 4,007,241 and 3,726,723), if the gate
fails, if `C` on `big` differs from `G_big`, if a cut search on a disclosed sample ends wrong in a
whole-layer run, or if the key-on library run differs from the script arm before the pass.

**Phase A (script layer).** 1. Link step 9's `G` runs and records; check totals; rerun the 2
full searches with the recorder. 2. Tests first: `tests/test_ar5_cut.py` (the cut on a synthetic
list: first N kept, level 0 only, quarter from the cut list; draws outside the grid unchanged).
3. Commit `ar5/cut/frozen.json` (rule, N, hashes) before any `held` run. 4. Replay the cut on the
23 cases (wrapping `run_discrete_search_spaced` in the script layer); gate. 5. `C` on `held`
first, then `mid`, `grad`, `big`; passes; validity, stats, the plot (run evaluations against wrong
voxels) before the tables; predictions. "READY FOR KEY".
**Phase B (key).** The key; unit tests; key-off suite with no golden change; key on, a library
run on `mid` equal to a script arm `Cl` (library draws) before the pass. Page
`findings_ar5_cut.md`, history, pareto rows.
**Files:** `scripts/autoresearch/ar5/cut_*.py`, `tests/test_ar5_cut.py`,
`benchmarks/autoresearch/ar5/cut/` (under 3 MB tracked; records and traces untracked), the page.

**Library changes:** phase A none. Phase B, only if the step passes: one opt-in key,
`GridLevel0Top N` (default 0, off), used by `AdaptiveVoxelReconstructor.reconstruct_voxel` after
level 0's discrete stage and before its quick MC; refused with `CandidateDedup` or
`SmoothCostRank`; `config_file.py`, `orientation_search.py` parse and check it. Off path
bit-identical; no default changes, no golden regeneration.

**Not in this step.** A cut of the level-0 screen or re-cost (about 79,000 of what is left). A
learned ranking of level-0 candidates (the 23 records are its data). A stronger local optimization
from the near ends. An early stop of MC or neighbour fits. The pass-skip key in the base arm.
Other values of N, 0.80, 0.95. Realistic images. Another crystal system. A timing. A second run
of any arm with another random seed.

**Choices made.**
- (a) The level-0 cut, not a stop after a level (the level bests above are often wrong and never
  near a 0.95 cost) and not a learned ranker yet: the discrete score already puts a candidate
  within 5 degrees of the end at place 165 or better in 23 of 23 (Wilson 95% 86% to 100%).
- N = 1,000, not 500 (3 times the largest place) or 2,000 (at least 12,000 more a search).
- Not a stronger local optimization from `G1`'s near ends first: 8 of 9 failures have one, but
  both the seed step's MC refinement and `G1`'s CMA stopped there under 0.95, so what would reach
  0.95 is unknown, and it cannot help a case like `mid` 16979. A next step; the cut search would
  be its fallback.
- Not (b) the MC stop `P(2000)`: about 15% of `G` on `mid`, 8% on `big`, 0.5% on `grad`; its
  `held` share is unknown (no `held` trace read); it changes 39 taken seed ends on `mid`, where
  `G` against `S8` is already at net 4 broken (limit 5); its replay costs about 6.5 million.
- Not (c) an early stop of neighbour fits: they are fixed at 61 evaluations and no trace is
  stored; this step records the traces at no cost.
- Not (d) the pass: the skip key of step 7 is measured (69% to 83% of the pass); stacking it is
  a combined run, not a new test.

### AR5 step 11: stop the seed step's MC refinement after a long run with no improvement

Task `ar5k`. Approved: D-20261113-N5. Starts from `feature/eval-reduction` at `8ffc5d4` (step 10
merged as a study, D-20261113-N4); reuses `order_run.py`, `order_pass.py`, `order_draws.py`,
`rescue_rules.py`. Phase A ends at "READY FOR KEY" or at its gate; phase B only after a decision.

**Why.** The owner asks for a way to tell whether more local optimization will still give a better
orientation, so that it can stop earlier. The seed step is arm `G`'s largest stage on `held`
(46.0%) and `mid` (41.7%); its MC refinements alone are 37.3% of `G` on `mid`, and about 76% of
their evaluations come after their last improvement (the "tail").

**What was read** (stored files only, nothing rerun). Step 9's cache (worktree `agent-a8d20b38...`),
`G_mid`, `G_grad`, `G_big`: `e1_log.json.gz` (seed steps' fits), `mc_traces.json.gz` (best cost at
each improvement; no hit ratio, no orientation), `records.json.gz` (seed rows), `rescue.json`; arm
`A`'s two logs only to check they equal `G`'s (yes on `mid`, `big`; `G_grad` has 2 more seed steps,
so `G`'s are used). Code: `refine_from_candidates`. Step 9's MC-stop draft, steps 9 and 10 here,
D-20261113-N4 and its review.

**How a seed step takes its end.** Per candidate (up to 5, in order): a CMA fit (201 evaluations);
if under hit ratio 0.95, an MC refinement; the kept end is the lower-cost of the two. The first
kept end at 0.95 stops the step and is taken; else the lowest-cost kept end is taken and the seed
test (0.80) chooses seed or rule `G1`'s rescue. Seed rows agree with this in 118 of 130, 87 of 97,
4 of 4 such steps (the rest were rescued). Each refinement has its own generator, so a stop after
its last improvement keeps its end exactly and moves no other draw.

**Where the MC evaluations go.** Reading line: one row per role of a refinement in its step; each
cell is refinements, then evaluations as a share of `G`'s run with the pass (tail in brackets).

| role of the MC refinement | mid | grad | big |
|---|---|---|---|
| reached 0.95, taken | 410; 17.6% (14.5%) | 137; 2.5% (2.1%) | 45; 8.7% (6.9%) |
| taken in a step with no 0.95 end | 103; 2.2% (1.3%) | 66; 1.1% (0.7%) | 3; 0.3% (0.2%) |
| not taken (late gains idle) | 1,092; 17.5% (12.4%) | 788; 7.1% (5.9%) | 105; 6.6% (4.9%) |

**Three kinds of stop rule, walked through whole seed steps** (arithmetic on the traces; a cost
estimate only, step 8's lesson): kept ends, whether the step still stops, which end is taken.
Assumed: a cut end of a refinement that never reached 0.95 stays under 0.95. Reading line: per
rule, the saving as a share of `G`'s run, then seed steps whose taken end changes.

| rule (fitted on `mid` and `grad`) | mid | grad |
|---|---|---|
| (i) `P(1000)`: stop after 1,000 evaluations with no improvement | 22.5%; 85 | 1.9%; 10 |
| (i) `P(2000)` | 15.1%; 37 | 0.4%; 2 |
| (i) `P(3000)` | 11.4%; 24 | 0.2%; 1 |
| (i) `P(5000)` | 7.2%; 14 | 0.0%; 0 |
| (ii) best of 120 cost-gated rules (1,000 if cost over 0.3, else 2,000) | 16.2%; 42 | 0.6%; 3 |
| (ii) cut at a checkpoint if cost tops every 0.95 refinement's (in-sample) | 1.0%; - | 2.1%; - |
| (iii) stop after 500 with no gain while worse than an end the step has | 9.0%; 54 | 2.2%; 7 |

(ii) is the "will it still help" predictor. At a stall of 1,000 evaluations no trace number ranks
refinements that improve again above the rest much better than chance (score 0.5 chance, 1.0
perfect): 0.38 to 0.71 for best cost, evaluations so far, number and size of gains, rate of fall.
So no learned model. (iii) fails: refinements worse than a kept end often come back to 0.95.

**Terms.** **MC refinement**: the seed step's stock MC local optimization from one candidate (MC,
then the variance stage). **`P(k)`**: stop after k evaluations with no strict gain of the best
cost, over the whole refinement; keep the first lowest-cost orientation; one final evaluation.
**Whole changed path**: the seed step with the cut, any candidates it now goes on to, the seed
test, `G1`'s rescue and the grid search. **Arm `S`**: `G` with `P(3000)`. **Right**/**wrong**: at
most / over 1 degree from the truth (symmetry-reduced). **Aligned draws**: as step 8.

**Samples.** `mid`, `grad`: traces and logs read, rule and k chosen on them; disclosed. `big`:
traces read, nothing chosen; disclosed. `held` (250 grains, gradients 0 to 4 degrees): this draft
read only arm `G`'s stage totals as printed in step 10's plan and saw the names `A_held`, `G_held`
in a cache listing; no `held` trace, log, record or map opened. `held` is the test sample.

**Frozen now.** `P(3000)`, seed-step MC refinements only; 0.80, 0.95, 5 candidates, `G1`, the grid
search unchanged. On `mid` it changes 24 of 1,607 taken ends (1.5%, Wilson 95% 1.0% to 2.2%): 22
where an MC end stopped the step (cut cost 0.118 median, at most 0.155, against 0.095 at the end),
2 with no 0.95 end; `grad` 1 of 1,280; `big` 4 of 317 (6.1% saved). If each lost 0.95 end made
its step run all remaining candidates, `mid` saves 9.4%. Clean images, Q-max 8.

**Measured**, per sample: evaluations by stage; cuts; taken ends changed (hit ratio, cost, error,
angle to `G`'s); failed seed tests, rescues, grid searches; wrong voxels (Wilson 95%); fixed and
broken against `G`, `S8`; median error. The replay records the hit ratio at each improvement.

**Predictions, written before any run.**
1. Replay: every uncut step and every refinement cut after its last improvement equals `G`'s;
   taken ends changed 20 to 28 (`mid`), 0 to 2 (`grad`), 3 to 5 (`big`); 12 to 20 of `mid`'s 22
   stop-type changes still reach 0.95; changed ends 0.01 to 0.3 degree from `G`'s (median under
   0.1), none wrong; seed-test failures within 1 of `G`'s on each sample.
2. Replay saving, extra evaluations counted: `mid` 9% to 11.4%; `big` 4% to 6.1%; `grad` 0 to 0.3%.
3. `held`: seed-step evaluations 10% to 30% fewer than `G`'s; total 4% to 12% fewer. 4. In total
   `mid` 8% to 12% fewer; `big` 4% to 7%; `grad` -1% to 1%.
5. Broken minus fixed against `G` between -4 and 4; median error moves at most 0.0005 degree.

**Decision rules, fixed now.** Gate before any whole-layer run (step 10's lesson: replay the whole
changed path, not one stage): replay every seed step with an MC refinement on `mid`, `grad`, `big`
along its whole changed path. It passes if (g1) every uncut step, and every refinement cut after
its last improvement, equals `G`'s stored one (orientation, cost); (g2) taken ends changed at most
30 (`mid`), 3 (`grad`), 6 (`big`); (g3) no changed taken end is wrong where `G`'s was right, and at
most 1 more failed seed test than `G` per sample; (g4) the replay saves at least 8% of `G`'s run on
`mid`, every extra evaluation counted. Else stop, report, no tuning, no other k. Then `S` passes on
a sample if, against `G` (aligned draws): (1) on `held` the total at least 4% fewer and seed-step
evaluations at least 10% fewer; the total at least 7% fewer on `mid`, 3% on `big`, at most 1% more
on `grad`; (2) broken minus fixed at most 5; (3) median error at most 0.001 degree higher; (4) at
most 1 more wrong seed (taken seed end wrong where `G`'s was right) and grid searches within 1 of
`G`'s; (5) the house rule against `S8` with the pass: at least 20% fewer evaluations, broken minus
fixed at most 5, median error at most 0.005 degree higher. If `G` fails (5) on a sample, `S` is
judged there by (1) to (4). The step passes if `S` passes on `held` and on all disclosed samples.
If it fails on `held`: no key.

**The margin on `mid`, said plainly.** `G` against `S8` on `mid` is at broken minus fixed 4 (limit
5). `S` changes about 24 taken ends there, each changing later growth; step 8 lost rule 5 by 2 to 3
voxels this way (within draw noise; `S8` is not draw-aligned). So rule 5 on `mid` is about even
odds, whatever the rule's merit. The plan does not move the limit, drop `mid`, or rerun with other
draws. It lowers the risk only through k (3000, not 2000: 24 changed ends, not 37), states the odds
now, and accepts the result: a fail on `mid` alone means no key and a study merge.

**Machine.** Gate replay: about 5.5 million evaluations plus at most 1 million for new rescues or
grid searches (about 20 process-minutes). Arm `S`: about 27 to 30 (`held`), 11.5, 10.2, 3.6
million; about 60 million, 2.8 to 3.2 process-hours, about 1.5 hours wall with 4 processes. Phase B:
two `mid` runs before the pass, about 18 million, 30 minutes. Caps 300 minutes `held`, 240 the rest.

**Stop rule:** stop and report if a stored count does not reproduce (totals 13,004,457, 10,196,692,
3,807,747; seed steps 5,417,024, 1,528,284, 687,902; traces 1,605, 991, 153), if the gate fails, if
a cut after the last improvement moves an end, or if the key-on run differs before the pass.

**Phase A (script layer).** 1. Link step 9's `G` runs; check counts; `patience_whatif.py` writes
the tables above. 2. Tests first: `tests/test_ar5_patience.py` (`P(k)` on synthetic traces: stop
point, first lowest cost kept, final evaluation; the whole-step walk; other draws unchanged).
3. Commit `ar5/patience/frozen.json` (rule, k, hashes) before any `held` run. 4. The gate replay
(`order_check.py`'s step replay, extended with the seed test, `G1` and the grid search). 5. `S` on
`held` first, then `mid`, `grad`, `big`; passes; validity, stats, the plot (run evaluations against
wrong voxels, `G` and `S`) before the tables; predictions. "READY FOR KEY".
**Phase B (key).** The key; unit tests; key-off suite with no golden change; key on, a library run
on `mid` equal to a script arm `Sl` (library draws) before the pass. Page
`findings_ar5_patience.md`, history, pareto rows. **Files:** `scripts/autoresearch/ar5/patience_*`,
`tests/test_ar5_patience.py`, `benchmarks/autoresearch/ar5/patience/` (under 3 MB tracked).

**Library changes:** phase A none. Phase B, only if the step passes: one opt-in key,
`SeedIndexMCPatience k` (default 0, off), only in the seed step's MC refinement (`index_seed`): a
counting wrapper keeps the first lowest-cost orientation, ends after k evaluations with no strict
gain, keeps the final evaluation; parsed in `config_file.py`, `orientation_search.py`. Off path
bit-identical; no default or golden change.

**Not in this step.** A stop that reads the running best's hit ratio ("stalled, still under 0.95"):
the most direct "will it still help" test, but no trace holds a hit ratio (the replay records it).
Patience for neighbour fits (fixed 61), revisits (no trace), the grid search's final local
optimization (about 1.3% of a search). A learned predictor. Other k, 0.80, 0.95. Realistic images.
Another crystal system. A timing. A rerun with another random seed.

**Choices made.**
- Plain patience, not a cost-based predictor: no trace number ranks well, and the best fitted
  rule buys 1.1 points over `P(2000)` on `mid` for 5 more changed ends. Not (iii): 54 changed ends.
- k = 3000, not 2000 (15.1%, 37 changed) or 5000 (7.2%, 14): the largest k whose pessimistic
  `mid` saving stays near 10% (9.4%), with fewer changed ends for `mid`'s margin.
- Not a stop when the running best first reaches 0.95: 410 taken ends change on `mid`, and with
  orientation gradients in grains a less polished seed shifts its growth.
- `grad`-like layers gain almost nothing (90% of gains end by evaluation 278); `held` decides.

### AR5 step 12: stop a stalled MC refinement only when its best already reaches 0.95

Task `ar5l`. Approved: D-20261113-N7. Starts from `feature/eval-reduction` after step 11's study merge
(D-20261113-N6); reuses `patience_rules.py`, `patience_replay.py`, `patience_run.py`, `order_*.py`,
`rescue_rules.py`. Phase A ends at "READY FOR KEY" or at its gate; phase B only after a decision.

**Why.** The owner asks for a general early stop: predict whether more local optimization will
still improve a voxel. Step 11's `P(3000)` saved 10.6% on `mid` but changed 20 taken seed ends
there and missed house rule 5 by one voxel. Most of its saving came from refinements already at
hit ratio 0.95 when they stalled; those are taken anyway, so a cut there keeps the step's path.

**What was read** (stored files only, nothing rerun). Step 11 (worktree `agent-a69be862...`):
`tables.md`, `replay_*.json` (row keys), and for `mid`, `grad`, `big` only `cache/{G,S}_<s>/`
`e1_log.json.gz` and `mc_traces.json.gz`; `patience_rules.cut_point`. Steps 10 and 11 here,
D-20261113-N4 to N6, the review `20261113-1000-review-ar5k.md`. Scratch scripts pair each S fit
with G's by (voxel, occurrence, candidate). A finding: step 11's replay did not store the hit ratio
at each improvement (its rows hold step totals only). The running best's hit ratio is on record at
one point only: the first 3,000-evaluation stall, as S's cut end. So k stays 3000.

**What the stored record says.** Reading line: per rule, the share of `G`'s run (with the pass)
saved, then seed steps whose taken end changes. "Gated high" is this step's rule; "gated low"
is its mirror (abandon a stalled refinement still under 0.95). From traces, a cost estimate only.

| rule (k = 3000) | mid | grad | big |
|---|---|---|---|
| `P(3000)`, step 11, measured whole run (changed: matched steps) | 10.60%; 20 | 0.14%; 1 | 5.01%; 4 |
| gated high: cut only if the best is at 0.95 or more | 7.5% to 8.3%; 7 to 11 | 0.16%; 0 | 4.5%; 1 |
| gated low: cut only if the best is under 0.95 | 3.3%; 13 | 0.04%; 1 | 1.7%; 3 |

How the `mid` row of "gated high" is built. 178 cuts come after the last gain with a G end at 0.95
or more: end unchanged, 893,957 evaluations (6.87%). 7 more come at a second stall, after a first
under 0.95 (52,886, 0.41%). 7 come before a later gain with the stall best at 0.95 (S's cut end):
the step stops at the same candidate, 0.011 to 0.078 degree from G's end (32,874, 0.25%). 4 come
before a later gain with an unknown stall hit ratio and a taken G end at 0.95: up to 4 more changed
ends and 92,765 (0.71%). Known: 7 changed of 1,607 steps, 0.4% (Wilson 95% 0.2% to 0.9%), against
20 to 25 for `P(3000)`. `big`: 1 (0.032 degree). In all 8 known cases a stall best at 0.95 belongs
to a refinement whose G end is at 0.95.

Other k, gated high. Reading line: exact saving (cuts after the last gain); bound on changed ends.

| k | mid exact; bound | grad exact; bound | big exact; bound |
|---|---|---|---|
| 2000 | 8.40%; 32 | 0.33%; 0 | 3.97%; 5 |
| 3000 | 6.87%; 22 (7 known, 11 under 0.95 so not cut, 4 unknown) | 0.16%; 0 | 3.81%; 4 |
| 5000 | 4.78%; 14 | 0.04%; 0 | 2.99%; 3 |

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels; a seed step ends
at the first kept end with 0.95 or more. **Stall**: k evaluations of an MC refinement with no
strict gain of its best cost. **Stall check**: at a stall, one evaluation of the running best's
hit ratio. **Rule `PH(3000)`**: in the seed step's MC refinement, at each stall check, stop if the
best's hit ratio is at least 0.95 (keep the first lowest-cost orientation; the check is the final
evaluation); else go on, and check again 3,000 evaluations after the next gain. **Arm `H`**: arm
`G` with `PH(3000)`. **Whole changed path**: as step 11. **Right**/**wrong**: at most / over 1
degree from the truth (symmetry-reduced). **Aligned draws**: as step 8. **Frozen now.**
`PH(3000)` in seed-step MC refinements only; 0.80, 0.95, 5 candidates, CMA fits, `G1`, the grid
search unchanged. Not tuned after any run.

**Samples, and what has touched them.** `mid`, `grad`, `big`: logs and traces read, the rule
chosen on them; disclosed. `held`: run under `G` and `S` in step 11, whose implementer and review
opened its seed logs, traces and maps. This draft opened only step 11's `tables.md` and the review,
which print `held`'s stage totals, its 763 cuts (420 under 0.95, 343 at 0.95 or more), 41 changed
taken ends and its paired counts; no `held` log, trace, record or map. `H`'s cuts are a subset of
`S`'s, so its seed-step outcome on `held` follows almost exactly from the stored `S_held` and
`G_held` logs. So `held` is no longer an untouched test sample for this rule: it becomes disclosed.
Honest test: a new sample `held2`, built by `held_make_sample.py` exactly as `held` (250 Voronoi
grains, gradients 0 to 4 degrees) with generator seed 2: a new sample for a carry-over test (AR3
step 11's reasoning), not an extra random seed for a comparison. Cost: sample and images about 15
minutes; `S8` about 60 million evaluations (from step 11's 51.6%), `G` about 31, `H` about 29.
Clean images, Q-max 8.

**Measured**, per sample: evaluations by stage; stall checks, cuts; taken ends changed (hit ratio,
cost, error, angle to G's); failed seed tests, rescues, grid searches; wrong voxels (Wilson 95%);
fixed and broken against `G`, `S8`; median error. Matched and unmatched seed steps both counted
(review lesson). New: the replay and runs store the hit ratio at every stall check.

**Predictions, written before any run.**
1. Replay: every uncut refinement, every cut after the last gain and every refinement whose stall
   best is under 0.95 equals G's; changed taken ends 7 to 11 (`mid`), 0 (`grad`), 1 to 2 (`big`),
   all at 0.95 or more at the same candidate, 0.005 to 0.15 degree from G's; failed seed tests,
   rescues and grid searches equal to G's.
2. Replay saving, extra evaluations counted: `mid` 7.0% to 8.3%; `big` 4.0% to 4.8%; `grad` 0.1%
   to 0.2%. 3. Whole layer: `mid` 6% to 9% fewer; `big` 3.5% to 5.5%; `grad` -0.5% to 0.5%;
   `held` 2% to 6%; `held2` total 2% to 8% fewer, seed step 6% to 18% fewer. 4. Broken minus
   fixed against `G` in -3 to 3; median error moves at most 0.0005 degree. 5. `mid` passes rule 5.

**Decision rules, fixed now.** Gate before any whole-layer run: replay every seed step with an MC
refinement on `mid`, `grad`, `big` along its whole changed path. It passes if (g1) as prediction
1's first clause; (g2) changed taken ends at most 12 (`mid`), 1 (`grad`), 3 (`big`); (g3) no
changed taken end wrong where G's was right, failed seed tests equal to G's; (g4) the replay saves
at least 6% of G's run on `mid`, every extra evaluation counted. Else stop, report, no tuning, no
other k. Then `H` passes on a sample if, against `G` (aligned draws): (1) on `held2` the total at
least 2% fewer and the seed step at least 5% fewer; at least 5% fewer on `mid`, 2.5% on `big`,
fewer on `held`, at most 1% more on `grad`; (2) broken minus fixed at most 5; (3) median error at
most 0.001 degree higher; (4) at most 1 more wrong seed (matched or unmatched) and grid searches
within 1; (5) the house rule against `S8` with the pass: at least 20% fewer evaluations, broken
minus fixed at most 5, median error at most 0.005 degree higher. If `G` fails (5) on a sample, `H`
is judged there by (1) to (4). The step passes if `H` passes on `held2` and on `held`, `mid`,
`grad`, `big`. If it fails on `held2`: no key.

**The margin on `mid`, said plainly.** `G` is at net 4 of limit 5 against `S8` there. `H` changes
7 to 11 taken seed ends on `mid` (step 11: 20), each by under 0.1 degree, and any change can move
later growth: step 11's 2 broken voxels came through changed growth paths, not changed ends. So
rule 5 on `mid` stays a risk: lower than step 11's even odds, roughly one in three (a judgement,
not a measurement). No stop can change zero ends: 7 known cuts come before a gain unseen at the
stall. The limit is not moved, `mid` not dropped, no rerun: a fail on `mid` alone means no key.

**Machine.** Gate replay about 6 million evaluations (about 20 process-minutes). Then, after the
freeze commit: `held2` built; `S8`, `G`, `H` on `held2` and `H` on `held` side by side (about 60,
31, 29, 29 million); then `H` on `mid`, `grad`, `big` (about 12, 10.2, 3.6 million). About 181
million, 8.2 to 9.5 process-hours, about 3.5 hours wall with 4 processes. Caps: 480 minutes for
`S8` on `held2` (`held`'s needed over 120), 300 the rest. Phase B: two `mid` runs, 18 million.

**Stop rule:** stop and report if a stored count does not reproduce (totals 13,004,457, 10,196,692,
3,807,747; seed steps 5,417,024, 1,528,284, 687,902), if the gate fails, if a cut after the last
gain or an uncut refinement moves an end, if a `held2` run hits its cap, or if the key-on run
differs from the script arm before the pass.

**Phase A (script layer).** 1. Link step 9's `G` runs and step 11's caches; check counts;
`stallhit_whatif.py` writes the tables above. 2. Tests first: `tests/test_ar5_stallhit.py` (stop
at 0.95, go on under it, re-check after a gain, first lowest cost kept, one evaluation per check,
whole-step walk, other draws unchanged). 3. Commit `ar5/stallhit/frozen.json` (rule, k, hashes of
every run, replay and helper script) before building `held2`. 4. Gate replay. 5. Build and render
`held2` (seed 2); nobody reads it before its runs end. 6. Runs; passes; validity, stats, the plot
(run evaluations against wrong voxels, `G` and `H`) before the tables; predictions. "READY FOR KEY".
**Phase B (key).** The key; unit tests; key-off suite with no golden change; key on, a library run
on `mid` equal to a script arm `Hl` (library draws) before the pass. Page
`findings_ar5_stallhit.md`, history, pareto rows with a status field. **Files:** `stallhit_*`
scripts and test, `benchmarks/autoresearch/ar5/stallhit/` (under 3 MB tracked).

**Library changes:** phase A none. Phase B, only if the step passes: one opt-in key,
`SeedIndexMCStallHit k` (default 0, off), only in the seed step's MC refinement (`index_seed`): a
counting wrapper keeps the first lowest-cost orientation and, after k evaluations with no strict
gain, evaluates its hit ratio once; ends at 0.95 or more, else goes on. Parsed in `config_file.py`,
`orientation_search.py`. Off path bit-identical; no default or golden change.

**Not in this step.** Gated low. Other k, 0.80, 0.95. A stronger local optimization from near ends.
Neighbour, revisit or pass stops. Learned predictors. Realistic images. Other crystals. A timing.
Another random seed.

**Choices made.**
- (a) gated high, not gated low: low saves 3.3% on `mid` for 13 changed ends, 11 of them a later
  0.95 end lost (the step then runs more candidates). High saves 7.5% or more for 7 to 11, each
  ending at 0.95 at the same candidate; 185 of its 192 known `mid` cuts keep the end exactly.
- k = 3000, not 2000 (up to 32 changed ends) or 5000 (4.78%): the only k with the stall hit ratio
  on record, so the changed-end count is known, not guessed.
- Not (b) a stronger local optimization from `G1`'s near ends: 8 of 9 failures had an end within
  0.58 degree, but what would reach 0.95 from there is not on record (MC and CMA both stopped under
  it), 9 cases from `mid` and `grad` is thin, and it changes every rescued seed end. Later.
- Not (c) neighbour or revisit stops: no trace is stored (step 10's recorder never ran a layer;
  step 11's caches hold none). `held2`, not `held`, as the test: `held`'s result is nearly fixed.

### AR5 step 13: arm H with the pass skip, as the best current recipe and the skip's clean test

Task `ar5m`. Approved: D-20261113-N9. Starts from `feature/eval-reduction` after step 12's merge
(D-20261113-N8); phase A needs no new key; phase B needs step 12's key `SeedIndexMCStallHit`.
Reuses `gate_live.py`, `gate_rules.py`, step 12's cache. Phase B only after a recorded decision.

**Why.** The owner asks for fewer evaluations and a general early stop: predict whether more
local optimization can still help a voxel. Under arm `H` the boundary pass is 16.8% to 36.1% of
the run. Step 7's key `BoundaryPassSkipHitRatio` is such a stop (a voxel whose own end already
reaches hit ratio 0.95 gets no try), but it was never run with `H`, and its only test on an unread
sample was not clean: variant 2 was chosen after `held` had been read live (D-20261112-N5).

**What was read** (stored files only, nothing rerun). PLAN steps 7, 10, 11, 12; D-20261113-N4 to
N8; reviews `20261113-1000-review-ar5k.md`, `20261113-1600-review-ar5l.md`. Step 12 (worktree
`agent-a0ffe667...`, `ar5/stallhit/`): `tables.md`; `stats_<s>.json` (totals of `H`, `S8`) and,
for `mid`, `grad`, `big`, `held` only: `cache/passes_<s>/D4_H/tries.npz`, `cache/stats/`
`part2_H_<s>_err.npz`, `G_<s>/rescue_log.json.gz`, `passes_<s>/before_G.npy`; one row of
`levels_mid.json`. Step 7's `live_held_S11v2_h0.95_v2.json`. Code: `run_boundary_pass`.

**Where arm `H`'s evaluations go.** Reading line: one row per sample; each stage as a share of
`H`'s run with the pass (step 12's `tables.md`); grid searches in brackets.

| sample | total | seed step | rescue | grid search | neighbour fits | revisits | boundary pass |
|---|---|---|---|---|---|---|---|
| held2 | 31,143,763 | 44.5% | 0.6% | 23.3% (21) | 4.8% | 10.1% | 16.8% |
| held | 29,588,170 | 43.8% | 0.6% | 22.2% (19) | 5.0% | 10.8% | 17.6% |
| mid | 11,950,397 | 36.5% | 0.3% | 11.2% (4) | 12.5% | 14.1% | 25.5% |
| grad | 10,180,191 | 14.9% | 0.3% | 16.7% (5) | 14.5% | 17.5% | 36.1% |
| big | 3,635,256 | 14.2% | 0 | 0 (0) | 41.2% | 13.4% | 31.2% |

After the seed step: the pass leads on `mid`, `grad`, `big`; the grid search on `held`, `held2`.

**The skip on `H`'s stored passes** (this draft, arithmetic). The replay is exact for the library
pass: candidate lists are built from the map before the pass, before any try, and a skipped voxel
draws its tries' seeds and discards them, so every kept try equals the full pass's (step 7, `held`
variant 2: 23,557 kept tries, 0 differences). Reading line: per disclosed sample, the pass, the
share skipped, then `HK` against `H` and against `S8` (after the pass).

| sample | pass evaluations | cut | HK total | fewer than H | fixed / broken vs H | vs S8: fewer; broken minus fixed | median shift (deg) | top own hit of a voxel wrong before |
|---|---|---|---|---|---|---|---|---|
| mid | 3,052,257 | 68.7% | 9,854,315 | 17.5% | 0 / 0 | 62.7%; +4 | +0.00016 | 0.9388 |
| grad | 3,679,520 | 82.8% | 7,134,949 | 29.9% | 0 / 0 | 60.4%; -7 | +0.00002 | 0.9411 |
| big | 1,135,332 | 76.5% | 2,767,226 | 23.9% | 0 / 0 | 66.8%; -3 | +0.00008 | 0.9486 |
| held | 5,195,553 | 72.8% | 25,806,780 | 12.8% | 0 / 0 | 56.9%; -1 | +0.00011 | 0.9218 |

Skipped voxels that the full pass would have moved: 277, 113, 307, 191; each by at most 0.27
degree, none across the 1-degree line. The margin under 0.95 is thin on `big` (0.9486).

**Terms.** **Boundary pass, try, own hit ratio**: as step 7 (a try is one 61-evaluation CMA local
optimization from a neighbour's orientation; own hit ratio is that of the voxel's end before the
pass, already scored by the pass). **Skip** (`h` = 0.95): a boundary voxel with own hit ratio at
least 0.95 gets no try. **Arm `HK`**: arm `H` with the skip. **Replay**: `H`'s stored tries with the
skipped ones removed. **Live pass**: the library pass with the key, on `H`'s map before the pass.
**Whole changed path**: here the pass alone: nothing before the pass reads the key. **Right** /
**wrong**: at most / over 1 degree from the truth (symmetry-reduced). **Aligned draws**: as step 8.

**Samples, and what has touched them.** `mid`, `grad`, `big`: step 7 chose 0.95 on their passes;
this draft replayed the skip on `H`'s; disclosed. `held`: step 7 ran and read a live skipped pass
on it; this draft replayed `H`'s; disclosed. `held2` (step 12, built after its freeze): step 12's
implementer and review read its seed logs, traces, maps and pass summaries. This draft opened only
`tables.md` (stage totals, wrong counts, paired counts) and listed `passes_held2/`. The skip's
outcome needs each voxel's own hit ratio, stored only in `D4_H/tries.npz`: its access time equals
its write time (13:14:07), so nobody read it, and step 12's scripts never load `old_score`. `h`
was fixed in step 7, before `held2` existed. So `held2` is an honest, untouched test of the skip.
Knowing that 97 voxels were wrong before `H`'s pass and 67 after does not say where the fixes sit
in own hit ratio. This draft did not open `held2`'s rescue logs. Clean images, Q-max 8.

**Frozen now.** `h` = 0.95 (step 7); arm `H` unchanged. Nothing is tuned after any run.

**Measured**, per sample: pass evaluations with and without the skip; skipped voxels, tries; run
total; wrong voxels after the pass (Wilson 95%); fixed and broken against `H` and `S8` (not
draw-aligned: voxel-paired only); median error; own hit ratio of each voxel the full pass moved
across the 1-degree line; live against replay, try by try.

**Predictions, written before any run.**
1. The replay reproduces the table above exactly; every live pass equals its replay try by try,
   with the same skipped voxels and evaluations (61 per kept try).
2. `held2`: the skip cuts 68% to 78% of `H`'s pass (5,227,517); `HK` needs 11% to 13.5% fewer
   evaluations than `H`; fixed / broken against `H` 0 / 0 (at most 0 / 1); median error up by at
   most 0.0003 degree; the highest own hit ratio of a voxel wrong before the pass is under 0.95.
3. `HK` against `S8` on `held2`: 54% to 56.5% fewer; broken minus fixed -5 to -3 (`H`: -4).

**Decision rules, fixed now.** Gate before `held2`'s tries are read: (g1) the committed replay
script reproduces the table (cuts 2,096,082, 3,045,242, 868,030, 3,781,390; 0 fixed, 0 broken);
(g2) the live pass with the key on `mid`'s map before `H`'s pass equals the replay in every kept
try (voxel, candidate, end, cost, hit ratio, evaluations) and in the skipped set; (g3) `mid`'s
broken minus fixed against `S8` stays 4. This replays the whole changed path. Else stop, report,
no tuning. Then `HK` passes on a sample if, against `H`: (1) the pass is cut by at least 50% and the
total is fewer (on `held2` at least 8% fewer); (2) broken minus fixed is at most 5% of `H`'s pass
net gain (wrong before minus after the pass), rounded up, at least 1 (`held2`: 30, so 2);
(3) median error at most 0.001 degree higher; (4) the house rule against `S8` with the pass: at
least 20% fewer evaluations, broken minus fixed at most 5, median error at most 0.005 degree
higher. The step passes if `HK` passes on `held2` and on `held`, `mid`, `grad`, `big`. Then `HK`
is the best current recipe (pareto row, a recipe config); a default is the owner's choice. If it
fails on `held2`: the skip stays opt-in, the page says its clean test failed, no other `h`.

**The margin on `mid`, said plainly.** `H` is at net 4 of limit 5 against `S8` there. The skip
changes no wrong voxel on `mid` in the exact replay (0 fixed, 0 broken; net stays 4), and the live
pass equals the replay by construction. So rule 5 on `mid` is not at risk in this step except
through a bug, which (g2) stops. The step does not widen the margin: any later rule that moves
growth on `mid` keeps the one-voxel risk. The limit is not moved, `mid` not dropped.

**Machine.** Replay: seconds. Live passes with the key on `H`'s stored maps (kept tries plus one
score per boundary voxel): about 0.98 (`mid`), 0.66 (`grad`), 0.29 (`big`), 1.44 (`held`), 1.2 to
1.7 (`held2`) million evaluations; about 5 million, 20 to 30 process-minutes, at most 2 processes
side by side. Phase B: one recipe run on `mid`, about 9.9 million, 35 minutes. About 15 million in
all, under 1.5 hours wall. Caps 60 minutes a pass, 120 the run.

**Stop rule:** stop and report if a stored count does not reproduce (`H`'s passes 3,052,257,
3,679,520, 1,135,332, 5,195,553, 5,227,517; totals as above), if the gate fails, if a live pass
differs from its replay in a kept try or the skipped set, if `held2`'s tries are read before the
freeze commit, or if the recipe run differs from the script arm (before the pass, and in it).

**Phase A (script layer).** 1. Link step 12's cache (`H` runs, passes, error files, `S8`
errors); check counts. 2. Tests first: `tests/test_ar5_recipe.py` (replay on a synthetic pass:
skipped voxels keep their end, kept tries unchanged, counts; the judge; a guard that refuses
`held2` before the freeze). 3. Commit `ar5/recipe/frozen.json` (`h`, rules, hashes of `H`'s runs
and passes, every script) before any `held2` tries are read. 4. Replay on `mid`, `grad`, `big`,
`held`; gate (live `mid` pass). 5. Live passes on `grad`, `big`, `held`, then `held2`; `held2`'s
replay from the same tries as a cross-check. 6. Stats; the plot first (run evaluations against
wrong voxels after the pass, `S8`, `G`, `H`, `HK`, one panel per sample); tables; predictions.

**Phase B (recipe).** After step 12's phase B is merged: a recipe config (`SeedRescueCMA`,
`SeedIndexMCStallHit 3000`, `BoundaryPass 1`, `BoundaryPassSkipHitRatio 0.95`); one library run on
`mid`. Before the pass it equals step 12's key-on library run `Hl` (fits, map, counts); its pass
equals the replay of the skip on `Hl`'s pass, try by try. Page `findings_ar5_recipe.md` (plot
first), history, README line, pareto rows with a status field ("best current recipe, opt-in keys,
not a default"). **Files:** `recipe_*` scripts and test, `benchmarks/autoresearch/ar5/recipe/`
(under 2 MB tracked).

**Library changes:** none. All three keys exist (steps 7, 9, 12); off paths untouched.
**Not in this step.** A longer local optimization or local grid search from the rescue's ends.
Neighbour, revisit or grid stops. Other `h`. A learned predictor. A timing (the machine is shared;
`HK` is the arm for a later quiet timing against `S8`). Realistic images. Other crystals. Another
random seed. A new sample.

**Choices made.**
- (c), the combined run, because the data supports it best: on four disclosed samples its effect
  is exact arithmetic, 12.8% to 29.9% fewer than `H` with 0 fixed and 0 broken, and `held2` gives
  the skip its first clean test at a cost of about 15 million evaluations.
- Not (a) yet. Opened: the failed rescues of `mid` (4), `grad` (5), `held` (19), with the map
  before the pass. A rescue end lies within 0.6 degree of the grid end in 20 of 28 and within
  1.06 in 27, but it is the lowest-cost end in only 13; on `held` 2 grid ends are wrong (33, 40
  degrees) with a rescue end 0.1 to 0.5 degree from them. What local optimization reaches 0.95
  from 0.1 to 0.6 degree is not on record (the MC refinement and the CMA both stopped at 0.59 to
  0.81), and step 10 showed such readings fail later stages. Next: a probe replay of a polish from
  every rescue end on the 28 (under 1 million evaluations), `held2`'s rescue logs as the test. If
  it worked it could remove much of the grid search, 23.3% of `H` on `held2`.
- Not (b) revisits or neighbour fits: no per-fit trace in any cache (step 10's recorder never ran
  a layer); neighbour fits are a fixed 61 evaluations. A recording run would come first.

### AR5 step 14: a second local optimization from a failed rescue's ends, before the grid search

Task `ar5n`. Approved: D-20261114-N4. Starts from `feature/eval-reduction` (step 12 merged at fe3d277);
reuses `rescue_common.cma_rescue`, `rescue_rules.py`, `stallhit_*` helpers and step 12's cache
(read only). Phase A ends at "READY FOR KEY" or at a gate; phase B only after a recorded decision.

**Why.** The owner asks for fewer evaluations, and for a general test of "will more local
optimization help?". Under arm `H` the grid search is 23.3% of the run on `held2` and 22.2% on
`held`. It runs only when rule `G1` fails. The stored record says the grid search mostly ends next
to an end `G1` already reached, at a hit ratio of 0.64 to 0.81. A second, different local
optimization from those ends may reach 0.95 for a few thousand evaluations, not 345,000.

**What was read.** PLAN steps 9, 12, 13; D-20261113-N1 to N10, D-20261114-N1 to N3; `STATE.md`;
step 12's reports (grep); `rescue_rules.py`, `rescue_common.py`, `rescue_replay.py`,
`stallhit_stats.py`, `sampling.get_structured_local_grid`. Stored, step 12 cache: `H_<s>/`
`rescue_log.json.gz` and `records.json.gz` of `mid`, `grad`, `held` only; for `held2` file times
only (`stat`). Nothing under `passes_held2/` was opened or listed with sizes.

**What the stored record says.** Reading line: one row per sample; arm `H`'s run with the pass,
its grid searches (one per failed rescue), their evaluations, how many grid ends reach hit ratio
0.95, and the best hit ratio among the rescue ends of a failed case.

| sample | H total | failed rescues | grid evaluations (share of H) | grid end at 0.95 | best rescue end, hit ratio |
|---|---|---|---|---|---|
| held2 | 31,143,763 | 21 | 7,252,187 (23.3%) | not opened | not opened |
| held | 29,588,170 | 19 | 6,563,817 (22.2%) | 13 | 0.64 to 0.80 |
| mid | 11,950,397 | 4 | 1,334,072 (11.2%) | 4 | 0.72 to 0.77 |
| grad | 10,180,191 | 5 | 1,697,698 (16.7%) | 4 | 0.66 to 0.81 |
| big | 3,635,256 | 0 | 0 | - | - |

On the 28 disclosed cases (step 13's draft, map before the pass): a rescue end lies within
0.6 degree of the grid end in 20 (Wilson 95% 53% to 85%), within 1.06 in 27. That near end is the
lowest-cost one in only 13 (30% to 64%), so picking an end by cost does not work. On `held` 2 grid
ends are wrong (33, 40 degrees), each with a rescue end 0.1 to 0.5 degree away. The grid end reaches
0.95 in 21 of 28 (57% to 87%); those 21 hold 7,181,436 of the 9,595,587 grid evaluations (74.8%).
That is the ceiling: a rule that skips a grid search only when it finds a 0.95 end nearby can save
at most about 75% of them, about 17% of `H` on `held2`, less its own cost. So the idea can pay.
What reaches 0.95 from 0.1 to 0.6 degree is not on record: that is what the probe measures.

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels. **Seed test**: a
seed step passes at 0.80. **Kept end**: per seed candidate, the lower-cost end of its CMA fit and
MC refinement. **Rescue (`G1`)**: after a failed seed test, a CMA local optimization
(first step 0.25 degree, 1,000 evaluations) from each kept end in order of cost, at most 5; the
first end at 0.95 or more is taken. **Rescue end**: the end of one `G1` try (stored).
**Failed rescue**: no rescue end at 0.95. **Grid search**: the global search (no starting guess)
that runs after a failed rescue, about 345,000 evaluations. **Rule `G2`**: after a failed rescue,
a second local optimization from every rescue end in `G1`'s order; the first end at 0.95 or more is
taken; if none, the grid search runs unchanged. Two forms, both frozen now. `G2a` (more of the
same): CMA from the rescue end, first step 0.5 degree, 3,000 evaluations. `G2b` (short starts
nearby): the library's fixed local grid about the rescue end (`get_structured_local_grid`, side
1.5 degrees, level 2: 64 points about 0.375 degree apart), each scored once, then CMA (0.1 degree,
500 evaluations) from its 2 lowest-cost points. Own generator (`default_rng(1)`); every evaluation
counts, one hit ratio check per CMA end included: at most 15,005 (`G2a`) or 5,330 (`G2b`) a case.
**Arm `J`**: arm `H` with `G2`. **Right** / **wrong**: at most / over 1 degree from the truth
(symmetry-reduced). **Aligned draws**: as step 8; `G2` takes no draw from the book, so a grid
search after a failed `G2` equals `H`'s. **Probe replay**: `G2` from the stored rescue ends, grid
evaluations from the log. It replays the whole changed path at the voxel (`G2` acts only after a
failed rescue; all other seed steps stay as in `H`), but growth after the voxel can change, which
only whole layers show. So it is a cost estimate and a voxel check, not the result.

**Samples, and what has touched them.** `mid`, `grad`, `held`: their 28 failed rescues were read by
step 13's draft and this one, and `G2`'s values are chosen on them; disclosed. `big`: no failed
rescue, so `J` must equal `H` there (a parity check). `held2` is the test; no new sample and no new
random seed. Disclosed limit: step 12's review opened `held2`'s seed logs (`e1_log.json.gz`, read
16:01 and 16:03, 13 Nov), its map (15:55) and `rescue.json` (counts: 64 failed seed tests, 43
rescued, 21 grid searches). So it could know which 21 voxels failed and how good their grid ends
are. `held2`'s rescue logs were never opened: in `G_held2` and `H_held2` their access time equals
their write time (13:05:49, 13:04:31). `G2`'s outcome depends on those rescue ends and a new
computation, so nobody has seen it. This draft did not open them. `passes_held2/` holds step 13's
unopened test: step 14 opens nothing under it until step 13's `held2` result is on record. Clean
images, Q-max 8.

**Measured**, per probe case: each `G2` try (start, end, hit ratio, cost, evaluations, error,
angle to the grid end); accepted or not; evaluations saved. Per whole layer: evaluations by stage;
failed seed tests, rescues, `G2` accepted, grid searches left; wrong voxels (Wilson 95%); fixed
and broken against `H` and `S8`; median error; wrong seeds, matched and unmatched; error of each
accepted end and its angle to `H`'s grid end.

**Predictions, written before any run.**
1. Re-scoring each stored rescue end gives its stored cost and hit ratio (to 1e-9).
2. Disclosed probe: the better form accepts 12 to 20 of 28, `G2b` at least as many as `G2a`; none
   wrong; 0.002 to 0.2 degree from the grid end; net saving 40% to 70% of 9,595,587.
3. `held2` probe: 7 to 15 of 21 accepted, none wrong; net saving 30% to 70% of 7,252,187.
4. `J` against `H`, with the pass: `held2` 6% to 16% fewer; `held` 6% to 15%; `mid` 0% to 8%;
   `grad` 0% to 10%; `big` identical. Grid searches: `held2` 21 to 6 to 14; `held` 19 to 5 to 12.
5. Broken minus fixed against `H` in -2 to 3; median error moves at most 0.0005 degree.

**Decision rules, fixed now.** Gate 1, disclosed probe, for each form: (g0) prediction 1; (g1) at
least 13 of 28 accepted; (g2) no accepted end wrong where the grid end is right; (g3) median error
of accepted ends at most 0.05 degree, none over 0.3; (g4) grid evaluations of accepted cases minus
all `G2` evaluations on the 28 at least 40% of 9,595,587. If both forms pass, take the larger
(g4); a tie goes to `G2b`. If neither passes: stop, report, no tuning, no third form, no whole
layer. Gate 2, `held2` probe of the chosen form (after the freeze commit): (t1) at least 7 of 21
accepted; (t2) as (g2); (t3) net saving at least 25% of 7,252,187. Else stop, no whole layer.
Then `J` passes on a sample if, against `H` (aligned draws): (1) total with the pass at least 6%
fewer on `held2`, at least 5% fewer on `held`, at most 0.5% more on `mid` and `grad`, identical on
`big`; (2) broken minus fixed at most 5; (3) median error at most 0.001 degree higher; (4) at most
1 more wrong seed (matched or unmatched); (5) the house rule against `S8` with the pass: at least
20% fewer evaluations, broken minus fixed at most 5, median error at most 0.005 degree higher. If
`H` fails (5) on a sample, `J` is judged there by (1) to (4). The step passes if `J` passes on
`held2` and on `held`, `mid`, `grad`, `big`. If it fails on `held2`: no key.

**The margin on `mid`, said plainly.** `H` is at net +4 of limit 5 against `S8` there. `J` can
change at most 4 seed ends on `mid` (its 4 failed rescues), and growth from them can move later
voxels. So rule 5 on `mid` is a risk: about one in four (a judgement, not a measurement). The limit
is not moved, `mid` is not dropped, no rerun: a fail on `mid` alone means no key.

**Machine.** Gate 1: at most 28 x 20,335 evaluations plus re-scoring, about 0.6 million, under 10
process-minutes, 2 processes; gate 2 at most 0.32 million. Whole layers (`J` with its pass): about
27, 26, 11.5, 9.7, 3.6 million (`held2`, `held`, `mid`, `grad`, `big`), 78 million; at about 20
million an hour (step 12) about 4 process-hours, 1.5 hours wall with 4 processes; caps 300 minutes
a run. Phase B: one `mid` run, 11.5 million. In all about 90 million, about 4.6 process-hours.

**Stop rule:** stop and report if a stored count does not reproduce (`H` totals and grid
evaluations above), if prediction 1 fails, if a gate fails, if `held2`'s rescue logs are read
before the freeze commit, if anything under `passes_held2/` is opened before step 13's `held2`
result is recorded, if a run hits its cap, if `J` on `big` differs from `H`, or if the key-on run
differs from the script arm before the pass.

**Phase A (script layer).** 1. Link step 12's cache read only (not `passes_held2`); check counts;
write only into step 14's own cache (D-20261113-N10). 2. Tests first: `tests/test_ar5_polish.py`
(order of tries, first 0.95 taken, grid call unchanged after a failed `G2`, own generator, aligned
draws untouched, every evaluation counted, guards on `held2`'s rescue logs and `passes_held2`).
3. Gate 1; choose the form. 4. Commit `ar5/polish/frozen.json` (form, values, hashes). 5. Gate 2.
6. `J` on all five samples (4 at a time), then its passes; checks before the pass. 7. After step
13's `held2` result is on record (else the owner decides): `held2` after the pass. Stats; the plot
first (run evaluations against wrong voxels after the pass: `S8`, `H`, `J`, a panel per sample);
tables; predictions. "READY FOR KEY".

**Phase B (key).** One opt-in key, `SeedRescuePolish` (the frozen form's values; default off),
used by `reconstructor.index_seed` after `G1` fails and before the grid search; parsed in
`config_file.py`, `orientation_search.py`. Unit tests; key-off suite with no golden change; key
on, a library run on `mid` equal to a script arm `Jl` (library draws) before the pass. Page
`findings_ar5_polish.md` (plot first), history, pareto rows with a status field.
**Files:** `ar5/polish_*.py`, `tests/test_ar5_polish.py`, `benchmarks/autoresearch/ar5/polish/`
(under 2 MB tracked). **Library changes:** phase A none; phase B the key, only if the step passes.

**Not in this step.** A choice of rescue end by cost (the near end is lowest in 13 of 28). A
change to `G1` or to rescued seeds. Other values of 0.80, 0.95 or `G2`. Arm `HK` (step 13,
blocked): `J` builds on `H`, the accepted recipe. A learned ranking of level 0. Neighbour or
revisit stops. Realistic images. Other crystals. A timing. Another random seed. A new sample.

**Choices made.**
- This idea, not (c) or (d): it can pay (ceiling about 75% of a stage that is 22% to 23% of `H`
  on `held`, `held2`); a case costs at most 15,005 evaluations against 345,000; a gate under 1
  million comes first. (c) needs a recording run first; (d) repeats step 10's risk.
- After a failed rescue only, not a stronger `G1`, which would change every rescued end (step 12).
- Two forms, chosen by a fixed rule on disclosed cases: `G2a` asks if more of the same helps; `G2b`
  does what the grid search's finer levels do (short starts close by), out to 0.75 degree.
- `held2` as the test, with its disclosed limit, not a new sample: its rescue logs are unread.

### AR5 step 15: a two-sided stall check in the seed step's MC refinement

Task `ar5o`. Approved: D-20261114-N6. Starts from `feature/eval-reduction` after the merge of
step 14's phase A (scripts only; its key still waits). Reuses `stallhit_rules.py`,
`stallhit_replay.py`, `polish_*` and step 14's cache (read only). Phase A ends at "READY FOR KEY"
or at a gate; phase B only after a recorded decision.

**Why.** The owner asks for fewer evaluations and for general early stops. Under arm `J` the seed
step's MC refinements are the largest part left: 46.5% of the run on `held`, 33.6% on `mid`, 11.8%
on `grad`. Of their evaluations, 67.6%, 71.1% and 81.4% come after the refinement's last gain. Rule
`PH(3000)` stops only at 3,000 evaluations with no gain and only at hit ratio 0.95. But the median
tail after the last gain is 1,058, 1,363 and 753 evaluations, so most refinements end on their own
budget first. Most MC evaluations (75%, 67%, 78%) go to refinements that end under 0.95.

**What was read.** PLAN steps 12 to 14; `STATE.md`; step 14's `ar5/polish/stats_*.json`,
`tables.md`; its cache (worktree `agent-a1bd2e7a7dbfe1efd`): `J_<s>/` `e1_log`, `mc_traces`,
`run`, `records` of `held`, `mid`, `grad`, `big`. Nothing of `held2`, no `passes_held2/`.

**What the stored record says.** Reading line: one row per sample; the parts of `J`'s run with the
pass. Full table and checks: `plan_step15_breakdown.md`.

| sample | J total | MC (share) | MC after the last gain | revisits | boundary pass | grid left |
|---|---|---|---|---|---|---|
| held | 24,457,773 | 11,381,609 (46.5%) | 7,694,387 (31.5%) | 13.1% | 21.3% | 5.7% |
| mid | 11,295,158 | 3,795,130 (33.6%) | 2,699,510 (23.9%) | 14.9% | 27.0% | 5.9% |
| grad | 9,175,867 | 1,078,188 (11.8%) | 877,460 (9.6%) | 19.4% | 40.1% | 7.4% |

The pass is step 13's (blocked). So the MC tail is the largest open part. What-if from traces.
"Exact": the first k-stall comes after the last gain, so the running best at the check is the
stored end, its hit ratio known; a stop there keeps the end. "Unsure": a gain comes after the
first k-stall; the hit ratio at that check is not stored. Reading line: per k, the exact saving of
rule `PB(k)` (below) as a share of `J`'s run; then the unsure ones (taken at 0.95 in brackets).

| k | held exact | mid exact | grad exact | big exact | unsure held; mid |
|---|---|---|---|---|---|
| 500 | 2,443,984 (9.99%) | 1,048,943 (9.29%) | 294,038 (3.20%) | 149,731 (4.12%) | 1,173 (250); 312 (83) |
| 1000 | 1,467,261 (6.00%) | 696,861 (6.17%) | 104,873 (1.14%) | 104,350 (2.87%) | 651 (138); 193 (58) |

Unsure ones hold 14.4% more on `held` (k = 500); the bracketed ones are the risk (a lost 0.95 end).

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels. **Seed step**:
index the voxel, cost 20 candidates, then for the best 5 in turn a CMA fit and, under 0.95, an MC
refinement; the step ends at the first end with 0.95 or more; else the lowest-cost end is taken.
**Seed test**: the taken end passes at 0.80; else rescue `G1`, then `G2b`, then the grid search
(the global search, about 345,000 evaluations). **Gain**: a strict drop of the running best cost.
**Stall**: k evaluations with no gain. **Stall check**: one evaluation (counted) of the running
best's hit ratio. **Rule `PB(k)`, frozen now**: in the seed step's MC refinement, at each stall
check, stop if the hit ratio is at least 0.95 or under 0.80; keep the first lowest-cost
orientation as the end; else go on and check again k evaluations after the next gain. 0.95 is the
step's stop value and 0.80 the seed test value: no new number. It replaces `PH(3000)`; a stopped
end is a kept end like any other. Two forms, both frozen now: `PB(1000)` and `PB(500)`. CMA fits,
5 candidates, 0.80, 0.95, `G1`, `G2b` and the grid search stay as in `J`. **Arm `JB`**: arm `J`
with the chosen form. **Right**/**wrong**: at most / over 1 degree from the truth
(symmetry-reduced). **Aligned draws**: as step 8. **Changed taken end**: a seed step whose taken
end differs from `J`'s. **Lost 0.95 end**: `J` took this candidate's MC end at 0.95; `JB` stopped
it earlier under 0.80 and went on.

**Samples, and what has touched them.** `held`, `mid`, `grad`, `big`: traces read by this draft;
the forms and gate are set on them; disclosed. `held2`: not used. Step 12's review opened its seed
logs, so it is not clean for a seed-step rule, and its pass waits for step 13. Clean test: a new
sample `held3`, built by `held_make_sample.py` exactly as `held` (250 Voronoi grains, gradients 0
to 4 degrees) with generator seed 3. It is a new sample for a carry-over test (step 12's reasoning),
not an extra random seed for a comparison. Nobody reads its logs before its runs end. Clean
images, Q-max 8.

**Measured.** Probe: the running best's hit ratio, cost and orientation at every 500 and 1000
stall of every seed-step MC refinement; per form, stops, changed taken ends (error, angle to `J`'s,
hit ratio), lost 0.95 ends, failed seed tests, rescues, `G2b`, grid searches, net saving. Whole
layer: evaluations by stage; stall checks, stops high and low; wrong voxels (Wilson 95%); fixed and
broken against `J` and `S8`; median error; wrong seeds, matched and unmatched.

**Predictions, written before any run.**
1. The recording replay gives every `J` refinement's evaluations and end cost (to 1e-9), and the
   stored `PH(3000)` check hit ratios.
2. Probe: of the bracketed refinements, under a third have a stall best under 0.80 at k = 1000,
   under half at k = 500. Net saving `PB(1000)`: `mid` 6% to 9%, `held` 6% to 10%; `PB(500)`:
   `mid` 8% to 13%, `held` 9% to 15%. Changed taken ends on `mid` 5 to 20 (`PB(1000)`), 10 to 35
   (`PB(500)`); none wrong where `J` is right; no extra grid search on `mid`, `grad`, `big`.
3. `JB` against `J`, with the pass: `held3` 6% to 14% fewer; `held` 6% to 15%; `mid` 5% to 12%;
   `grad` 1% to 4%; `big` 2% to 5%.
4. Broken minus fixed against `J` in -3 to 4; median error moves at most 0.0005 degree.
5. `mid` passes rule 5 (about even odds; see the margin paragraph).

**Decision rules, fixed now.** Gate, for each form, on the four disclosed samples: (g0) prediction
1; (g1) no changed taken end wrong where `J`'s is right; (g2) grid searches equal to `J`'s on
`mid`, `grad`, `big`, at most 1 more on `held`; (g3) changed taken ends at most 30 on `mid`, 80 on
`held`, 5 on `grad`, 5 on `big`; (g4) net saving, every extra evaluation counted (checks, extra
candidates, rescues, `G2b`, grid searches), at least 6% of `J`'s run on `mid` and on `held`. If
both pass, take the larger summed net saving on `mid` and `held`; a tie goes to `PB(1000)`. If
neither passes: stop, report, no tuning, no third form, no whole layer. Then `JB` passes on a
sample if, against `J` (aligned draws): (1) total with the pass at least 5% fewer on `held3` and
`held`, at least 4% fewer on `mid`, at least 1% fewer on `big`, at most 0.5% more on `grad`;
(2) broken minus fixed at most 5; (3) median error at most 0.001 degree higher; (4) at most 1 more
wrong seed (matched or unmatched), grid searches within 1; (5) the house rule against `S8` with
the pass: at least 20% fewer evaluations, broken minus fixed at most 5, median error at most
0.005 degree higher. If `J` fails (5) on a sample, `JB` is judged there by (1) to (4). The step
passes if `JB` passes on `held3` and on `held`, `mid`, `grad`, `big`. Fails on `held3`: no key.

**The margin on `mid`, said plainly.** `J` is at net +4 of limit 5 against `S8` there (13 broken,
9 fixed). `JB` changes more taken ends than `H` did (step 12: 8 on `mid`); step 11 changed 20 and
missed rule 5 by one voxel through changed growth. So rule 5 on `mid` is about even odds (a
judgement, not a measurement). The limit is not moved, `mid` is not dropped, no rerun: a fail on
`mid` alone means no key and a study merge.

**Machine.** Probe: the recording replay re-runs every seed step's CMA fits and MC refinements on
the four samples, about 19 million evaluations, about 1 process-hour, 4 processes (`held` the
longest, about 40 minutes). Changed paths per form, under 3 million. Whole layers with the pass:
`JB` on `held`, `mid`, `grad`, `big` (about 21, 10, 9, 3.5 million); `held3` built (about 15
minutes), then `S8`, `J`, `JB` on it (about 60, 25, 22 million). In all about 172 million, about
8.6 process-hours at 20 million an hour, about 3.5 hours wall with at most 4 heavy processes (`S8`
on `held3` the long pole). Caps: 480 minutes for `S8` on `held3`, 300 the rest. Phase B: one `mid`
run, about 10 million.

**Stop rule:** stop and report if a stored count does not reproduce (`J` totals 24,457,773,
11,295,158, 9,175,867, 3,635,256; MC 11,381,609, 3,795,130, 1,078,188, 424,450), if prediction 1
fails, if the gate fails, if `held3`'s logs are read before its runs end, if anything of `held2`
is opened, if a run hits its cap, or if the key-on run differs from the script arm before the pass.

**Phase A (script layer).** 1. Link step 14's cache read only; check the counts above; write only
into step 15's own cache. 2. Tests first: `tests/test_ar5_band.py` (stop at 0.95 and under 0.80,
go on between, re-check after a gain, first lowest cost kept, one evaluation per check, stopped
end kept, whole-step walk, other draws untouched, guards on `held2`). 3. Recording replay;
prediction 1. 4. Commit `ar5/band/frozen.json` (both forms, hashes) before reading the stall hit
ratios. 5. Gate; choose the form; commit it. 6. Build `held3`. 7. Runs (4 at a time), passes.
Stats; the plot first (run evaluations against wrong voxels after the pass: `S8`, `J`, `JB`, a
panel per sample); tables; predictions. "READY FOR KEY".

**Phase B (key).** One opt-in key, `SeedIndexMCStallLow h` (default 0, off), read beside
`SeedIndexMCStallHit k`: at a stall check, also stop under h. Used only in `index_seed`'s MC
refinement; parsed in `config_file.py`, `orientation_search.py`. Unit tests; key-off suite with
no golden change; key on, a library run on `mid` equal to a script arm `JBl` (library draws)
before the pass. Page `findings_ar5_band.md` (plot first), history, pareto rows with a status
field. **Files:** `ar5/band_*.py`, `tests/test_ar5_band.py`, `benchmarks/autoresearch/ar5/band/`
(under 3 MB tracked). **Library changes:** phase A none; phase B the key, only if the step passes.

**Not in this step.** Other k or low values. The pass skip (step 13). A learned stop. Revisit or
neighbour stops. CMA, candidates, 0.80, 0.95. Realistic images. Other crystals. A timing. Another
random seed. `held2`.

**Choices made.**
- This idea: it attacks the largest open part (MC, 36.2% of `held`+`mid`+`grad`), its exact part
  alone saves 6% to 10% on `held` and `mid`, and a 1-process-hour probe decides it before any
  whole layer. Its weak side: little on `grad` (MC 11.8%) and `big`.
- Two-sided, not low only: step 12's gated-low at k = 3000 lost 11 later 0.95 ends on `mid`;
  this check stops low only under 0.80, and the high side at k under 3000 keeps the step's path.
- Runner-up (a), a learned stop: a small model on the trace so far, the CMA end and the stall hit
  ratio predicts a later gain to 0.95. It needs this step's recorded stall hit ratios to train.
- Runner-up (b), a revisit start test: skip a revisit whose start is near one already rejected
  there. Revisits are 13% to 19% of the run, 53% to 81% rejected; needs a recording run first.
- `held3`, not `held2`: `held2`'s seed logs were opened by step 12's review.

### AR5 step 16: the high stop alone at a 500-evaluation stall in the seed step's MC refinement

Task `ar5p`. Approved: D-20261114-N8. Starts from `feature/eval-reduction` after step 15's study merge
(053fd66). Reuses `band_rules.py` (its walk with no low side), `band_record.py`, `band_paths.py`,
`stallhit_*`, step 15's recording rows and step 14's cache (read only). Phase A ends at "READY FOR
KEY" or at a gate; phase B only after a recorded decision.

**Why.** The owner asks for fewer evaluations and for general early stops. Under arm `J` the seed
step's MC refinements are 46.5% of the run on `held`, 33.6% on `mid`, 11.8% on `grad`. Step 15's
two-sided check failed its gate, and all its damage came from the low side: a stall best under 0.80
does not predict a failing refinement. Its high side broke nothing. Rule `PH(3000)` (step 12, in
`J`) already stops high, but only after 3,000 evaluations with no gain; the median tail after the
last gain is 1,058 / 1,363 / 753. This step asks: does the same high stop after only 500
evaluations with no gain save about 6% on `held` and `mid`, with the same wrong voxels?

**What was read.** PLAN steps 12, 14, 15; step 15's review; `ar5/band/gate.json`,
`gate_tables.md`; step 15's recording rows (`band/cache/record/*.jsonl`, worktree
`agent-a38701500191fde36`) and step 14's `J_<s>/` `e1_log`, `mc_traces` (worktree
`agent-a1bd2e7a7dbfe1efd`) of `held`, `mid`, `grad`, `big`. Nothing of `held2`.

**What the stored record says.** The drafter walked `PH(500)` and `PH(1000)` over step 15's recorded
stall events (`plan_step16_candidates.md`). Reading line: one row per sample; per form, steps whose
taken end moves, steps whose path changes, and the net saving as a share of `J`'s run with the pass
(every check counted; growth after a changed voxel not replayed).

| sample | `PH(500)` moved / path change | `PH(500)` net saving | `PH(1000)` moved / path | `PH(1000)` saving |
|---|---|---|---|---|
| held | 8 / 0 | 1,584,152 (6.48%) | 3 / 0 | 4.58% |
| mid | 6 / 1 | 740,633 (6.56%) | 6 / 1 | 4.89% |
| grad | 1 / 0 | 134,974 (1.47%) | 0 / 0 | 0.83% |
| big | 1 / 0 | 88,144 (2.42%) | 1 / 0 | 1.87% |

A moved step keeps its path: the stopped end and `J`'s end are both at 0.95, so the step ends at
the same candidate; only the seed's orientation differs a little. The one path change (`mid` step
466: stall best at 0.95, `J`'s later end under 0.95) needs a live replay. These numbers are on
disclosed samples and were seen before this plan: the gate below checks them; it does not test.

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels. **Seed step**:
index the voxel, cost 20 candidates, then for the best 5 in turn a CMA fit and, under 0.95, an MC
refinement; the step ends at the first end at 0.95 or more; else the lowest-cost end is taken.
**Seed test**: the taken end passes at 0.80; else rescue `G1`, then `G2b`, then the grid search
(the global search, about 345,000 evaluations). **Gain**: a strict drop of the running best cost.
**Stall**: k evaluations with no gain. **Stall check**: one counted evaluation of the running
best's hit ratio. **Rule `PH(k)`** (step 12's rule, frozen there): at a stall check, stop if the
hit ratio is at least 0.95 and keep the first lowest-cost orientation as the end; else go on and
check again k evaluations after the next gain. **This step's form, frozen now: `PH(500)`**, in place
of `PH(3000)`. No low side, no other k. CMA fits, 5 candidates, 0.80, 0.95, `G1`, `G2b` and the
grid search stay as in `J`. **Arm `J5`**: arm `J` with `PH(500)`. **Right**/**wrong**: at most /
over 1 degree from the truth (symmetry-reduced). **Aligned draws**: as step 8; each MC refinement
has its own generator, so a stop moves no other draw. **Changed taken end**: a seed step whose taken
end differs from `J`'s.

**Samples, and what has touched them.** `held`, `mid`, `grad`, `big`: disclosed; the form and the
gate numbers were seen on them. `held2`: not used, not opened. Clean test: the new sample `held3`,
built by `held_make_sample.py` exactly as `held` (250 Voronoi grains, gradients 0 to 4 degrees)
with generator seed 3 (step 15's `band_make_held3.py`, untested end to end). It is a new sample for
a carry-over test, not an extra random seed. Nobody reads its logs before its runs end. Clean
images, Q-max 8.

**Measured.** Gate: per sample, high stops, checks that went on, moved steps, path changes, net
saving; for the path-change steps, the live path (taken end, error, seed test, rescue, grid
search). Whole layer: evaluations by stage; stall checks and stops; wrong voxels (Wilson 95%);
fixed and broken against `J` and `S8`; median error; wrong seeds, matched and unmatched.

**Predictions, written before any run.**
1. The gate reproduces the table above to the evaluation on the arithmetic steps (moved 8 / 6 / 1
   / 1, path changes 0 / 1 / 0 / 0) and the stored `J` counts.
2. `mid` step 466 live: the taken end stays right; no extra failed seed test or grid search.
3. `J5` against `J`, with the pass: `held3` 5% to 8% fewer; `held` 5.5% to 7%; `mid` 5.5% to 7.5%;
   `grad` 1% to 2%; `big` 1.5% to 3%.
4. Broken minus fixed against `J` in -2 to 2 on every sample; median error moves at most 0.0005
   degree; wrong seeds and grid searches equal to `J`'s.
5. `mid` passes rule 5 (about 2 in 3; see the margin paragraph).

**Decision rules, fixed now.** Gate, on the four disclosed samples, before any whole layer: (g0)
prediction 1; (g1) moved steps at most 10 on `held` and `mid`, at most 2 on `grad` and `big`; path
changes at most 2 per sample; (g2) on the live path-change steps, no taken end wrong where `J`'s is
right, no extra failed seed test, no extra grid search; (g3) net saving, every check counted, at
least 5.5% of `J`'s run on `held` and on `mid`. Gate fails: stop, report, no tuning, no other k, no
whole layer. Then `J5` passes on a sample if, against `J` (aligned draws): (1) total with the pass
at least 4% fewer on `held3`, 5% on `held` and `mid`, 1% on `big`, at most 0.5% more on `grad`;
(2) broken minus fixed at most 5; (3) median error at most 0.001 degree higher; (4) at most 1 more
wrong seed (matched or unmatched), grid searches within 1; (5) the house rule against `S8` with
the pass: at least 20% fewer evaluations, broken minus fixed at most 5, median error at most 0.005
degree higher. If `J` fails (5) on a sample, `J5` is judged there by (1) to (4). The step passes if
`J5` passes on `held3` and on `held`, `mid`, `grad`, `big`. Fails on `held3`: no key value.

**The margin on `mid`, said plainly.** `J` is at net +4 of limit 5 against `S8` there (13 broken,
9 fixed). So `J5` can afford at most one more net broken voxel on `mid` than `J`; two fail the
step. `J5` changes 6 taken ends on `mid`, all small moves at 0.95, fewer than step 12's `H` (8,
same wrong voxels as its baseline) and far fewer than step 11 (20, missed rule 5 by one voxel). So
rule 5 on `mid` is about 2 in 3 (a judgement, not a measurement). The limit is not moved, `mid` is
not dropped, no rerun: a fail on `mid` alone means no key value and a study merge.

**Machine.** Gate: arithmetic on stored rows (minutes) and one live seed step on `mid`. Whole layers
with the pass: `J5` on `held`, `mid`, `grad`, `big` (about 23, 10.5, 9, 3.5 million); `held3`
built (about 15 minutes), then `S8`, `J`, `J5` on it (about 60, 25, 23.5 million). In all about
155 million, about 7.7 process-hours at 20 million an hour, about 3.5 hours wall with at most 4
heavy processes (`S8` on `held3` the long pole). Caps: 480 minutes for `S8` on `held3`, 300 the
rest. Phase B: one `mid` run, about 10.5 million. The `held3` runs of `S8` and `J` stay as the clean
baselines for later steps.

**Stop rule:** stop and report if a stored count does not reproduce (`J` totals 24,457,773,
11,295,158, 9,175,867, 3,635,256; MC 11,381,609, 3,795,130, 1,078,188, 424,450), if prediction 1
fails, if the gate fails, if `held3`'s logs are read before its runs end, if anything of `held2`
is opened or listed, if a run hits its cap, or if the key-on run differs from the script arm before
the pass.

**Phase A (script layer).** 1. Link step 15's recording rows and step 14's cache read only; check
the counts above; write only into step 16's own cache. Keep step 14's freeze check active (step
15's review). 2. Tests first: `tests/test_ar5_high500.py` (stop at 0.95 only, never low, go on
under 0.95, re-check after a gain, first lowest cost kept, one evaluation per check, `PH(3000)`
absent, other draws untouched, guards on `held2` and on `held3` logs). 3. Commit
`ar5/high500/frozen.json` (the form, gate, decision rules, hashes) before any script runs on the
rows. 4. Gate; live replay of the path-change steps. 5. Pilot the `held3` build and one short `J`
run on it (setup lines only, no result read); then build `held3`. 6. Runs (4 at a time), passes.
Stats; the plot first (run evaluations against wrong voxels after the pass: `S8`, `J`, `J5`, a
panel per sample); tables; predictions. "READY FOR KEY".

**Phase B (key).** No new code: the key `SeedIndexMCStallHit k` (step 12, default 0, off) already
takes k = 500. Unit test that k = 500 parses and stops as the script rule; key-off suite with no
golden change; key on (with `J`'s other keys), a library run on `mid` equal to a script arm `J5l`
(library draws) before the pass. Page `findings_ar5_high500.md` (plot first), history, pareto rows
with a status field. **Files:** `ar5/high500_*.py`, `tests/test_ar5_high500.py`,
`benchmarks/autoresearch/ar5/high500/` (under 3 MB tracked). **Library changes:** none in phase A;
phase B a docstring note that k = 500 is tested, only if the step passes.

**Not in this step.** A low side of any kind. Other k (250, 1000). A stop at the first gain whose
best reaches 0.95. A learned stop. Revisit or neighbour start tests. CMA, candidates, 0.80, 0.95.
The pass skip (step 13). Realistic images. Other crystals. A timing (MC evaluations may cost more
wall time than others; pointer only). Another random seed. `held2`.

**Choices made.**
- This idea: the only part of step 15 that did no harm, with a saving known before any run (6.48%
  / 6.56% on `held` / `mid`), 8 / 6 / 1 / 1 changed taken ends, one live step for the gate, and an
  existing key. Weak side: little on `grad` (1.47%) and `big` (2.42%), and the `mid` margin.
- `PH(500)`, not `PH(1000)`: about 2 points more on `held` and `mid` with the same moved count on
  `mid`. Not `PH(250)`: its 250-stall hit ratios are not all recorded (855 unsure on `held`).
- Runner-up (a), a learned stop: an exploratory look (train on three disclosed samples, test on
  the fourth) adds about 1 point over this rule on `held` and loses 2 to 4 0.95 ends; a plain hit
  threshold does as well. A perfect stop would save 26.8% / 20.8%: the stall summary lacks the
  information. Needs richer inputs (per-peak) first.
- Runner-up (c), a reject-start test for revisits and neighbour fits: skip a fit whose start hit
  is very low. About 5% on `grad` in S8-era data, unknown on `held` and `mid`; needs a recording
  run of `J` first (about 2.4 process-hours). Wrong skips change growth.
- Runner-up (d1), stop at the first gain that reaches 0.95: at least 8.4% / 8.2% on `held` /
  `mid`, but it moves most seed ends at 0.95 (hundreds); try it once `mid` has margin.

### AR5 step 17: record the BFS local optimizations of arm `J5`, then an early skip chosen on it

Task `ar5q`. Approved: D-20261115-N2. Starts after step 16's merge (its review approved `J5`
unchanged, D-20261115-N1). Reuses step 8's aligned draws (`order_draws.py`),
step 16's `high500_*` and cache (read only). Phase A ends at "READY FOR KEY" or at a gate; phase B
only after a recorded decision. The pass skip (step 13's idea) stays an owner option.

**Why.** The owner asks for fewer evaluations and for general early stops. Under arm `J5` the BFS
local optimizations (neighbour fits and revisits, 61 evaluations each) are 18.8% of the run on
`held3`, 20.5% on `held`, 30.0% on `mid`, 36.0% on `grad`, 55.6% on `big`. Of them 56% to 73% are
rejected (25% on `big`), and a rejected fit almost never decides a final answer. Nothing stored
says which fits will lose before they run. This step records that, then tests one rule.

**What was read.** PLAN steps 14 to 16; `STATE.md`; `reconstructor.py` (`_expand`, reseed);
step 16's `J5_<s>/run.json`, `records.json.gz` (worktree `agent-a8d3953c88c593c42`); AR3 step
9's `ar3/start/feat_{grad,big}.npz` (`S8`-era fit records, worktree `agent-a9c88b610fb10e1b6`).
Scripts `s17c_rec.py`, `s17c_feat.py` (scratchpad). Nothing of `held2` was opened or listed.

**What the stored record says.** Reading line: per sample, `J5`'s run with the pass; rejected
revisits and neighbour fits not accepted at their first fit (counts); their evaluations as a share
of `J5` (the ceiling: every losing fit skipped at no cost); and revisits spent on voxels that in
the end were fitted by a reseed (a full search), as a share of `J5`.

| sample | `J5` total | rejected revisits | first-rejected neighbour fits | ceiling | revisits on reseeded voxels |
|---|---|---|---|---|---|
| held3 | 23,503,385 | 38,314 of 48,181 | 12,953 of 24,277 | 13.30% | 25,735 (6.68%) |
| held | 22,939,812 | 42,756 of 52,638 | 13,295 of 24,273 | 14.91% | 28,904 (7.69%) |
| mid | 10,589,493 | 18,778 of 27,618 | 10,297 of 24,409 | 16.75% | 9,923 (5.72%) |
| grad | 9,041,686 | 15,416 of 29,176 | 14,693 of 24,222 | 20.31% | 3,617 (2.44%) |
| big | 3,567,098 | 2,873 of 8,010 | 5,424 of 24,533 | 14.19% | 1,058 (1.81%) |

A rejected end is the final answer only at voxels that end unresolved (5 / 3 / 0 / 0 / 0); all
others are later replaced. `S8`-era records (other arm): a start test (skip if the start's hit
ratio is under t times the best) skips 26.2% of fits on `grad` at t = 0.3 with 14 of them
accepted, 35.0% with 119 at 0.4; `big` 7.4% with 7, 10.9% with 48. At 61 evaluations a fit and
one check per fit, t = 0.3 saves about 8.8% of `J5` on `grad`, 3.2% on `big`; others: no record.

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels. **Fit**: one BFS
local optimization (CMA, 60 evaluations plus one final hit evaluation) from an inherited start: a
**neighbour fit** (a voxel's first) or a **revisit** (a refit of a rejected voxel from a later
neighbour). **Best**: the expansion's running best hit ratio; a fit is accepted if its hit ratio
is over 0.9 times the best. **Wrong skip**: a skipped or stopped fit that `J5` accepted (it changes
growth). **Arm `J5`**: arm `J` with `PH(500)` (step 16). **Arm `J5C`**: `J5` with the chosen rule.
**Right**/**wrong**: at most / over 1 degree from the truth (symmetry-reduced). **Aligned draws**:
as step 8; a skipped fit still takes and discards its draw, so no other draw moves.

**Recording run (`J5r`).** `J5` re-run on the five disclosed samples, aligned draws, the pass
unchanged. Per fit it stores: voxel, kind, expanding voxel, order, the best, start orientation,
start hit ratio and cost, the running best's cost at every evaluation and its hit ratio after 20
evaluations, evaluations, end orientation, hit ratio and cost, accepted or not, kept by the push
rule or not, and whether the end is the voxel's final answer before and after the pass. The two
hit evaluations per fit are recording-only and listed apart; `J5r`'s counted run must equal `J5`.

**Rule family, frozen now, before any recording is read.** A skipped or stopped fit is a rejected
fit (the voxel becomes or stays a refit voxel; a skipped neighbour fit holds its start, a skipped
revisit keeps the old end). Every check is one counted evaluation.
- `S(t, scope)`: before the fit, check the start; skip if its hit ratio is under t times the best.
- `C(t, scope)`: after 20 evaluations, check the running best; stop if under t times the best.
- `D(d)`: skip a revisit whose start lies within d degrees of a start already rejected at that
  voxel since its last reset (no evaluation).
t in {0.3, 0.4, 0.5, 0.6}; scope revisits, neighbour fits, or both; d in {0.1, 0.25, 0.5}. 27 rules.
**Selection, frozen now.** A rule is admissible if on each of the five samples (a) wrong skips are
at most 0.1% of that sample's fits and (b) at most 2 final answers before the pass change. Among
admissible rules take the largest sum over the five samples of the net saving as a share of `J5`
(saved evaluations of skipped fits that `J5` rejected, minus every check; wrong skips earn nothing).
Ties: fewer wrong skips in all, then `D` before `S` before `C`, then the smaller t or d.

**Samples, and what has touched them.** `held3`, `held`, `mid`, `grad`, `big`: disclosed; the rule
is chosen on them. `held2`: not used, not opened. Clean test: a new sample `held4`, built by
`held_make_sample.py` exactly as `held` (250 Voronoi grains on the 24,570-voxel grid, gradients
0 to 4 degrees) with generator seed 4, as step 16 built `held3`: a new sample for a carry-over
test, not an extra random seed. Nobody reads its logs or maps before its runs end. Clean, Q-max 8.

**Predictions, written before any run.**
1. `J5r` gives `J5`'s map before the pass at every voxel, its fit counts (table) and evaluations.
2. The chosen rule is `S(0.3, both)` or a `D` rule (about even odds); no `C` rule is chosen.
3. Net saving of the chosen rule on the recording: `held` 5% to 8%, `held3` 4.5% to 7.5%, `mid`
   5% to 9%, `grad` 7% to 10%, `big` 2% to 4% of `J5`; wrong skips 5 to 40 per sample.
4. `J5C` against `J5`, with the pass: `held4` 4% to 9% fewer; the disclosed samples within 1.5
   points of their recording saving; broken minus fixed -2 to 3; median error up at most 0.0005.
5. `mid` passes rule 5 (about 3 in 5; see the margin paragraph).

**Decision rules, fixed now.** Gate, on the disclosed samples, before `held4` is built: (g0)
prediction 1; (g1) at least one admissible rule; (g2) the chosen rule's net saving at least 4% of
`J5` on `held` and `mid`, at least 3% on `held3` and `grad`, at least 1% on `big`; (g3) exact
replay: on `mid` the live `J5C` run, up to its first wrong skip in run order, equals the replay
fit by fit (skips, ends, evaluations), and over the whole run if it has none. The replay is exact
for skips of fits that `J5` rejected: such a skip leaves the best, the refit set, the revisit
counts and the reseed choice (geometry only) as they were. Gate fails: stop, report, no tuning, no
second pick, no whole layer. Then `J5C` passes on a sample if, against `J5` (aligned draws):
(1) total with the pass at least 3% fewer on `held4`, `held3`, `held`, `mid`, `grad`, at least 1%
on `big`; (2) broken minus fixed at most 5; (3) median error at most 0.001 degree higher; (4) at
most 1 more wrong seed (matched or unmatched), grid searches within 1, reseeds at most 2% more;
(5) the house rule against `S8` with the pass: at least 20% fewer evaluations, broken minus fixed
at most 5, median error at most 0.005 degree higher. If `J5` fails (5) on a sample, `J5C` is
judged there by (1) to (4). The step passes if `J5C` passes on `held4` and on the five disclosed
samples. Fails on `held4`: no key value; study merge.

**The margin on `mid`, said plainly.** `J5` is at net +4 of limit 5 against `S8` there (13 broken,
9 fixed). `J5C` changes growth only through wrong skips; a skipped voxel can still be fitted by a
later revisit or a reseed, but each can also cost a voxel. So `J5C` can afford at most one more
net broken voxel on `mid` than `J5`. My judgement: about 3 in 5 (not a measurement). No limit is
moved, `mid` is not dropped, no rerun: a fail on `mid` alone means no key value.

**Machine.** `J5r`: `J5` before the pass is about 18.6, 17.7, 7.5, 5.4, 2.4 million (`held3`,
`held`, `mid`, `grad`, `big`) plus 0.6 million recording evaluations, about 52 million, 2.6
process-hours at 20 million an hour, about 1.7 hours wall on 4 processes (`J5` took 5,780 s on
`held3`, contended). Gate: arithmetic (minutes), then `J5C` on `mid` with its pass (about 10
million, 0.7 hours). Then `held4` built (about 15 minutes); `J5C` on `held3`, `held`, `grad`,
`big` (about 56 million); on `held4` `S8`, `J5`, `J5C` (about 67, 23.5, 22 million). In all about
230 million, about 11.5 process-hours, about 6 hours wall with at most 4 heavy processes. Caps:
480 minutes for `S8`, 300 the rest. `S8` and `J5` on `held4` stay as clean baselines.

**Stop rule:** stop and report if a stored count does not reproduce (`J5` totals 23,503,385,
22,939,812, 10,589,493, 9,041,686, 3,567,098; fit counts in the table), if prediction 1 fails,
if the gate fails, if anyone reads the recording before `ar5/fitskip/frozen.json` is committed,
if `held4`'s logs or maps are read before its runs end, if anything of `held2` is opened or
listed, if a run hits its cap, or if step 16's review changes `J5`.

**Phase A (script layer).** 1. Merge step 16; link its caches read only; check the counts above;
write only into this step's cache. 2. Tests first: `tests/test_ar5_fitskip.py` (each form skips
or stops only below its line; checks counted; a skipped fit is rejected and keeps the stated end;
discarded draws leave other fits unchanged; recording evaluations kept apart; guards on `held2`
and on `held4`). 3. Commit `ar5/fitskip/frozen.json` (family, selection, gate, rules, hashes).
4. `J5r` on the five samples; prediction 1. 5. Selection; commit the pick. 6. Gate. 7. Build
`held4` (pilot first), then the runs and passes, 4 at a time. Stats; the plot first (run
evaluations against wrong voxels after the pass: `S8`, `J5`, `J5C`, six panels); tables;
predictions. "READY FOR KEY".

**Phase B (key).** One opt-in key for the chosen rule (default off), read in `_expand`; parsed in
`config_file.py`, `orientation_search.py`. Unit tests; key-off suite with no golden change; key on,
a library run on `mid` equal to a script arm `J5Cl` (library draws) before the pass. Page
`findings_ar5_fitskip.md` (plot first), history, pareto rows with a status field.
**Files:** `ar5/fitskip_*.py`, `tests/test_ar5_fitskip.py`, `benchmarks/autoresearch/ar5/fitskip/`
(under 3 MB tracked; the recording stays in the cache). **Library changes:** none in phase A;
phase B the key, only if the step passes.

**Not in this step.** Two rules together. Accepting a start without a fit (AR3 step 9's `F1`
broke voxels through growth). Other fit budgets. The pass skip. Seed-step MC stops. A learned
rule (the recording can train one later). Realistic images. Other crystals. A timing. Another
random seed. `held2`.

**Choices made.**
- This idea: a ceiling of 13% to 20% of `J5`, never recorded on `held` or `mid`. The rule only
  rejects, so it cannot grow a wrong orientation; its risk is a missed accept. Weak side: no size
  on `held`, `held3`, `mid` before the recording; the `S8`-era numbers are another arm's; about
  11.5 process-hours, twice step 16. The family and selection are fixed first, so the pick is a
  procedure, not a look; `held4` tests it once.
- Not chosen: the pass skip (owner option, step 13); a seed MC stop at the first gain at 0.95
  (about 1.6 to 1.9 points over `PH(500)` on `held`, `mid`; moves hundreds of seed ends).

### AR5 step 18: a learned skip of BFS local optimizations, trained on the `J5r` recording

Task `ar5r`. Approved: D-20261115-N4. Starts from `feature/eval-reduction` after step 17's study merge
(D-20261115-N3). Reuses step 8's aligned draws, step 17's `fitskip_*` (its live skip mode and its
`held4` guards exist but are untested) and the recording `J5r` (worktree `agent-aeba15b7c3850b007`;
manifest `ar5/fitskip/manifest.json`), read only. Phase A ends at "READY FOR KEY" or at a gate;
phase B only after a recorded decision. Rule `F95` (below) stays an option.

**Why.** The owner asks for fewer evaluations, learned approximations of a repetitive search, and
early stops. BFS local optimizations are 18.8% to 55.6% of `J5`'s run; 56% to 73% are rejected (25%
on `big`); skipping every loser free would save 13.3% to 20.3%. Step 17's hand rule `S(0.3, both)`
saved 3.20% to 8.74% and failed on `mid`. This step trains a model to predict a rejection first.

**What was read.** PLAN steps 16, 17; step 17's reports; `J5r`'s rows. Scripts `s18_*.py` (scratch).

**What the stored record says.** Leave one sample out: train on four disclosed samples, score the
fifth; the threshold comes from the four only (procedure below). Reading line: per held-out sample,
net saving as a share of `J5` with the pass; wrong skips (limit: 0.1% of the sample's fits);
changed final answers (limit 2). Each check is one counted evaluation.

| sample | fits | `S(0.3, both)` | learned skip (this plan's form) | ceiling |
|---|---|---|---|---|
| held3 | 72,458 | 3.98%, 21, 1 | 5.46%, 30 of 72 allowed, 1 | 13.31% |
| held | 76,911 | 4.21%, 34, 0 | 6.16%, 70 of 76, 0 | 14.90% |
| mid | 52,027 | 3.62%, 26, 0 | 5.79%, 46 of 52, 0 | 16.75% |
| grad | 53,398 | 8.74%, 24, 0 | 11.83%, 43 of 53, 0 | 20.31% |
| big | 32,543 | 3.20%, 6, 0 | 5.38%, 17 of 32, 0 | 14.19% |

All five are inside the limits. Said plainly: the first sizing had a 17th feature (fit index over
the run's fit count), not known during a run; the table is without it. One other variant was seen
and not chosen (candidates, section 4). Seen before this plan: the gate checks; it does not test.

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels. **Fit**: one BFS
local optimization (CMA, 60 evaluations plus one final hit evaluation): a **neighbour fit** (a
voxel's first) or a **revisit** (a refit of a rejected voxel from a later neighbour). **Best**: the
expansion's running best hit ratio; a fit is accepted if its hit ratio is over 0.9 times the best.
**Kept end**: the end the push rule stores at the voxel. **Wrong skip**: a skipped fit that `J5`
accepted. **Changed final answer** (step 17's review): a skipped fit that `J5` rejected and whose
end was the voxel's final answer before the pass. It is narrow: rejected ends are final only at
unresolved voxels. Wrong skips, counted apart, hit fits whose end is final; the replay cannot tell
if a later revisit or reseed recovers them. **Arm `J5L`**: `J5` with the learned skip. **Right** /
**wrong**: at most / over 1 degree from the truth. **Aligned draws**: as step 8; a skipped fit
takes and discards its draw.

**Model, frozen now.** Label: `J5` accepted the fit. One row per fit of `J5r`, unweighted. Library:
scikit-learn 1.6.1 (already in the `dev` extra and `uv.lock`; training only).
`HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, random_state=0)`, all else at the
1.6.1 defaults: 31 leaves, 20 rows a leaf, no L2 term, 255 bins, early stopping on (over 10,000
rows: a 10% validation split drawn with seed 0, stop after 10 rounds with no loss gain). The trees
are exported to JSON (feature, threshold, children, leaf value, the baseline) and evaluated in plain
Python at run time: no new runtime dependency. Score: the predicted chance that the fit is accepted.

**The 16 features, frozen now** (all from before the fit; 4 to 6 come from the check):
1. kind: 1 for a revisit, 0 for a neighbour fit. 2. revisit count of the voxel (0 for a neighbour
fit). 3. best. 4. start hit ratio. 5. start cost. 6. start confidence. 7. start hit ratio over best.
8. distance to rejected starts: smallest misorientation (degrees, crystal symmetry from the config)
from the start to starts rejected at this voxel since its last reset (a revisit count that did not
go up); capped at 10, and 10 if none. 9. parent unknown: 1 if the expanding voxel has no kept fit
end yet (a seed, say). 10. parent hit: the hit ratio of that kept end (best if unknown). 11. parent
cost: its cost (-1 if unknown). 12. start hit ratio over parent hit (1 if unknown). 13. own hit:
the hit ratio of this voxel's kept fit end (-1 if none). 14. fits so far at the voxel. 15.
rejected fits so far at the voxel. 16. place in the expansion (the fit's index in it, from 0).
Only fit ends count as kept ends; seed, reseed and pass ends do not.

**Threshold procedure, frozen now.** Skip a fit if its score is under the threshold. Each training
sample is scored by a model trained on the other training samples (out of fold). Candidate
thresholds: the quantiles at 0, 0.1%, ..., 80% of all out-of-fold scores. Take the one with the
largest sum of net shares over the training samples among those where every training sample has
wrong skips at most 0.1% of its fits and changed final answers at most 2; none: no skip. The gate
model for `mid` is trained on `held3`, `held`, `grad`, `big` (threshold 0.00957 in the sizing). The
final model for `held4` is trained on all five (246 rounds, threshold 0.01112, 0.63 MB of JSON).

**Accounting.** Every fit is checked: one counted evaluation of the start gives its hit ratio,
confidence and cost (as rule `S`). A skipped fit costs no further evaluation. Features 1 to 3 and
7 to 16 need none. Net saving: evaluations of skipped fits that `J5` rejected, minus all checks.
A model call is not a cost-function evaluation: counted as zero, its wall time and the feature
time apart in `run.json`. Sized: 0.08 ms a fit, about 6 s a run; a cost evaluation is about 0.3 ms.

**Samples.** `held3`, `held`, `mid`, `grad`, `big`: disclosed, trained on; their whole layers are
in-sample consistency checks. `held2`: not used. Clean test, the only out-of-sample whole layer:
`held4`, built by `held_make_sample.py` exactly as `held` with generator seed 4 (a new sample, not
an extra random seed), unread until its runs end. Clean images, Q-max 8.

**Predictions, written before any run.**
1. The tracked script reproduces the table above to the evaluation.
2. Live `mid` (gate model): matches its replay to its first wrong skip; 4.5% to 7% under `J5`.
3. `J5L` (final model) against `J5`, with the pass: `held4` 4% to 7% fewer; `held3`, `held`,
   `mid`, `big` 5% to 7.5%; `grad` 10% to 13.5%; broken minus fixed -2 to 4; median error up at
   most 0.0005 degree; reseeds at most 2% more.
4. `mid` passes rule 5 (about 1 in 2; see the margin paragraph).

**Decision rules, fixed now.** Gate, before `held4` is built: (g0) prediction 1, and `J5r` has
72,458 / 76,911 / 52,027 / 53,398 / 32,543 fits; (g1) on each held-out sample, wrong skips at most
0.1% of its fits, changed final answers at most 2, net saving at least 4.5% (9% on `grad`); (g2)
the live `J5L` run on `mid` with the gate model, up to its first wrong skip in run order, equals
the replay fit by fit (features to 1e-9, scores, skips, ends, evaluations), and over the whole run
if it has none; (g3) that run's whole `mid` layer with its pass passes rules (1) to (5) below.
Gate fails: stop, report, no tuning, no other model or threshold, no whole layer. Then `J5L`
passes on a sample if, against `J5` (aligned draws): (1) total with the pass at least 3% fewer;
(2) broken minus fixed at most 5; (3) median error at most 0.001 degree higher; (4) at most 1 more
wrong seed (matched or unmatched), grid searches within 1, reseeds at most 2% more; (5) the house
rule against `S8` with the pass: at least 20% fewer evaluations, broken minus fixed at most 5,
median error at most 0.005 degree higher. If `J5` fails (5) on a sample, `J5L` is judged there by
(1) to (4). The step passes if `J5L` passes on `held4` and on the five disclosed samples. Fails on
`held4`: no key value; study merge.

**The margin on `mid`, said plainly.** `J5` is at net +4 of limit 5 against `S8` there (13 broken,
9 fixed). `J5L` can afford at most one more net broken voxel on `mid` than `J5`; two fail the step.
It makes about 46 wrong skips on `mid` (`S(0.3, both)`: 26), and the final threshold sits at
`mid`'s wrong-skip limit out of fold (52 of 52). My judgement: about 1 in 2 (not a measurement).
No limit is moved, `mid` is not dropped, no rerun: a fail on `mid` alone means no key value.

**Machine.** Training: about 31 tree fits, about 10 CPU-minutes on 2 processes. Gate: the live
`mid` run with its pass (about 10 million, 0.6 hours). Then `held4` built (about 15 minutes);
`J5L` on `held3`, `held`, `grad`, `big` (about 55 million); on `held4` `S8`, `J5`, `J5L` (about 67,
23.5, 22 million). In all about 178 million, about 9 process-hours at 20 million an hour, about 4.5
hours wall, at most 4 heavy processes. No new recording. Caps: 480 minutes for `S8`, 300 the rest.

**Stop rule:** stop and report if a stored count does not reproduce (`J5` totals 23,503,385,
22,939,812, 10,589,493, 9,041,686, 3,567,098; the fit counts above), if prediction 1 fails, if the
gate fails, if training runs before `ar5/learnskip/frozen.json` is committed, if `held4` is read
before its runs end, if anything of `held2` is opened or listed, or if a run hits its cap.

**Phase A (script layer).** 1. Merge step 17; link `J5r` read only by its manifest; check counts.
2. Tests first, `tests/test_ar5_learnskip.py`: each feature on hand-made rows; JSON evaluator equals
scikit-learn; skip only under the threshold; one counted check a fit; a skipped fit is rejected and
keeps the stated end; discarded draws; the `held2`, `held4` guards. 3. Commit `frozen.json`. 4.
Leave-one-out table (g0, g1). 5. Live `mid` (g2, g3). 6. Final model and JSON. 7. `held4` (pilot
first), runs, passes, 4 at a time. Plot first (evaluations against wrong voxels after the pass,
`S8`, `J5`, `J5L`, six panels); tables; predictions. "READY FOR KEY".

**Phase B (key).** One opt-in key `BFSLearnedSkip <weights file>` (default empty, off), read in
`_expand`; parsed in `config_file.py`, `orientation_search.py`. Weights: the final JSON (trees,
threshold, features, training hash), tracked, under 3 MB (0.63 MB sized), and a plain-Python
evaluator. Unit tests; key-off suite with no golden change; key on, a library run on `mid` equal
to a script arm `J5Ll` (library draws) before the pass. Page `findings_ar5_learnskip.md`, history,
pareto rows. **Files:** `ar5/learnskip_*.py`, `tests/test_ar5_learnskip.py`,
`benchmarks/autoresearch/ar5/learnskip/` (under 3 MB). **Library changes:** none in phase A.

**Not in this step.** A stop at 20 evaluations (needs a CMA stop hook); `F95`; rules `S`, `D`;
retraining on `held4`; other models, features; realistic images; other crystals; a timing; seeds.

**Choices made.**
- This idea: exact and out of sample on all five samples, inside the limits, 5.4% to 6.2% (11.8%
  on `grad`), larger than `F95`'s upper bound everywhere. Weak side: `mid`'s margin and its tight
  threshold; whole layers on disclosed samples are in-sample; wrong skips' effect on growth is not
  in the replay. A small MLP sized about the same but broke the wrong-skip limit on `held`.
- Runners-up: `F95` (seed MC stop at the first gain reaching 0.95), 1.34-3.31 / 1.73-4.59 /
  1.63-4.54 / 0.58-0.90 / 0.65-2.72%; `S(0.3, both)`, failed step 17's gate (3.62% on `mid`); a
  stop from the 20-evaluation trace, 5.61 / 5.66 / 4.65 / 9.48 / 4.05%, no live form.

### AR5 step 19: after the pass, refit contested voxels with the MC local optimization

Task `ar5s`. Approved: D-20261115-N6. Starts from `feature/eval-reduction` after step 18's study merge
(D-20261115-N5). Reuses step 16's cache (worktree `agent-a8d3953c88c593c42`, `ar5/high500/cache`:
`J5` runs, maps, pass tries), the stored `S8` runs, step 18's live `J5L` run of `mid` (worktree
`agent-adedc789a656fa5ed`), all read only. Phase A ends at "READY FOR KEY" or at a gate; phase B
only after a recorded decision.

**Why.** The owner asks for fewer evaluations; a cheaper approximation followed by a final
optimization is welcome. Arm `J5` uses about 60% fewer evaluations than `S8`, but on `mid` it sits
at broken minus fixed +4 against `S8` (limit 5). Step 17 fell short of its saving on `mid`; step
18's learned skip saved 5.65% there and broke 4 voxels (a replay cannot see what a wrong skip does
to BFS growth). This step saves no evaluations itself. It asks: can a cheap final repair after the
pass fix the kind of voxel that `J5` and step 18 get wrong, and so give `mid` room for savings?

**What was read.** PLAN steps 16 to 18; both gate reports; `J5`'s maps before and after the pass,
its pass tries and `stats_<s>.json`; `S8`'s maps and pass tries; the truths (scoring only); `J5L`'s
map, pass tries and `broken_mid.json`. Scripts `s19_*.py` (scratch). Nothing of `held2` opened.

**What the stored record says.** Every wrong voxel of `J5` after the pass (233 on five samples) is a
boundary voxel holding a neighbouring grain's orientation. The pass tried a start within 1 degree
of the truth at 219 of them, and at 218 a try ended right, but always at a lower hit ratio than the
wrong end, so the pass kept the wrong one. Where `S8` is right and `J5` wrong (82 voxels), `S8`'s
end has the higher hit at 79. `S8` fits with the MC local optimization (about 92 to 158 evaluations
a fit on `mid`); `J5`'s fits and pass tries are 60-evaluation CMA runs. Step 18's breaks: alike.
Reading line: per sample, the flag below: voxels flagged; wrong ones among them; repair fits; cost
as a share of `J5` (160 evaluations a fit, assumed); a proxy result (at a flagged voxel, `S8`'s end
if its hit is higher): fixed / broken against `J5`, and broken minus fixed against `S8`.

| sample | flagged | wrong caught | fits | cost | proxy vs `J5` | vs `S8`: `J5`, proxy |
|---|---|---|---|---|---|---|
| held3 | 193 | 28 of 54 | 1,217 | 0.83% | 14 / 0 | -7, -21 |
| held | 183 | 30 of 67 | 1,035 | 0.72% | 16 / 0 | -3, -19 |
| mid | 143 | 22 of 35 | 643 | 0.97% | 10 / 0 | +4, -6 |
| grad | 238 | 47 of 56 | 967 | 1.71% | 23 / 0 | -7, -30 |
| big | 56 | 18 of 21 | 161 | 0.72% | 4 / 0 | -3, -7 |
| mid, `J5L` | 146 | 26 of 39 | 658 | 1.05% | 15 / 0 | +8, -7 |

The proxy is another arm's ends, not a replay: it cannot show a new break, and a live MC fit may
fall short of `S8`'s end. The flag was chosen among 12 forms looking at these samples: seen before
this plan. The flag catches all four of step 18's extra breaks.

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels. **Boundary pass**:
step 5's pass in `J5` (a voxel with a neighbour more than 0.5 degree away gets one 60-evaluation
CMA try from each candidate; rule `hit` keeps the highest hit ratio). **Candidates**: a voxel's
neighbour orientations (radius 2 side lengths) more than 0.5 degree from its own, merged within
0.5 degree. **Ratio**: a voxel's hit ratio after the pass over the median of its neighbours'.
**Margin**: its hit ratio minus the best hit among its other pass ends (the old end and try ends
more than 1 degree from the kept end); a voxel with none is never flagged. **Contested voxel**:
flagged by `M(0.95, 0.15)`: ratio under 0.95 and margin under 0.15. No truth is used. **MC local
optimization**: the run's classic variance-minimizing MC fit (`optimizer="mc"`), as `S8`'s BFS
fits. **Arm `J5R`**: `J5` with the repair. **Right**/**wrong**: at most / over 1 degree from the
truth (symmetry-reduced). **Fixed**/**broken**: wrong in the older arm and right in the newer / the
reverse, counted by voxel.

**The repair, frozen now.** After the pass, on a snapshot of the map: flag the contested voxels.
For each, in index order: one MC local optimization from its own end and one from each candidate
of the snapshot; keep, among the current end and these ends, the highest hit ratio, then the lowest
cost, then the current end. Its own generator, `default_rng(0)`; no other draw moves. A voxel's hit
ratio is the one the pass or the BFS record stores; a missing one costs one counted evaluation (the
sizing used the stored score of every voxel). Every evaluation is counted. BFS, seeds and the pass
are untouched, so `J5R`'s run is `J5`'s stored run plus the repair: exact, no rerun.

**Samples, and what has touched them.** `held3`, `held`, `mid`, `grad`, `big`: disclosed; the flag
was chosen on them. `J5L` (`mid`): seen. `held2`: not used, not opened. Clean test: a new sample
`held4`, built by `held_make_sample.py` exactly as `held` (250 Voronoi grains, 24,570 voxels,
gradients 0 to 4 degrees) with generator seed 4, as step 16 built `held3`: a new sample, not an
extra random seed. Nobody reads its logs or maps before its runs and repair end. Clean, Q-max 8.

**Measured.** Flagged voxels, fits, evaluations (and per fit), voxels changed; wrong voxels (Wilson
95%), fixed and broken against `J5` and `S8`, median error; `J5L` the same, and against `J5R`.

**Predictions, written before any run.**
1. The scripts reproduce the table's flagged voxels and fits, and the stored counts.
2. MC fits at contested voxels take 100 to 250 evaluations (median); repair cost 0.5% to 2% of
   `J5` on `held3`, `held`, `mid`, `big`, 1% to 3.5% on `grad`.
3. Fixed against `J5`: `held3` 6 to 14, `held` 7 to 16, `mid` 4 to 10, `grad` 10 to 23, `big` 2
   to 4; broken 0 to 2 per sample; median error moves at most 0.0005 degree.
4. `mid` after the repair: broken minus fixed against `S8` from -6 to +1.
5. `held4`: repair cost under 2% of `J5`; fixed minus broken against `J5` at least 3.
6. `J5L` repaired, on `mid`: -7 to +2 against `S8`.

**Decision rules, fixed now.** Gate, on the disclosed samples, before `held4` is built: (g0)
prediction 1 (flagged 193 / 183 / 143 / 238 / 56, fits 1,217 / 1,035 / 643 / 967 / 161; `J5`
wrong 54 / 67 / 35 / 56 / 21); (g1) `J5R` on `mid`: fixed minus broken against `J5` at least 3,
repair evaluations at most 2.5% of `J5`; (g2) `J5R` on `held3`, `held`, `grad`, `big`: broken
minus fixed against `J5` at most 0, repair evaluations at most 2.5% (3.5% on `grad`). Gate fails:
stop, report, no tuning, no other form, no `held4`. Then `J5R` passes on a sample if, against `J5`
(the same run, so draws align by construction): (1) repair evaluations at most 2.5% of `J5`'s run
with the pass (3.5% on `grad`); (2) broken minus fixed at most 0; (3) median error at most 0.001
degree higher; (4) the house rule against `S8` with the pass: at least 20% fewer evaluations,
broken minus fixed at most 5, median error at most 0.005 degree higher. The step passes if `J5R`
passes on `held4` and on the five disclosed samples, and on `mid` broken minus fixed against `S8`
is at most +1. Fails on `held4`: no key value; study merge. `J5L` repaired is a measurement only:
it informs the owner; it does not reopen step 18.

**The margin on `mid`, said plainly.** `J5` is at +4 of limit 5 against `S8` there (13 broken, 9
fixed). This step aims to move that to +1 or lower, so a later step may cost up to 4 voxels there.
The repair can also break a voxel (an MC fit from a neighbour's orientation may beat a right end);
rule (2) catches that. My judgement: `mid` reaches +1 or lower about 3 in 4 (not a measurement).

**Machine.** Gate: the repair on five stored maps and `J5L`, about 0.75 million evaluations plus
setup, under 1 process-hour, 2 processes. Then `held4` built (about 15 minutes); `S8` and `J5` on
it with their passes (about 67 and 23.5 million); the repair (about 0.2 million). In all about 92
million, about 4.6 process-hours at 20 million an hour, about 3.6 hours wall (`S8` the long pole),
at most 4 heavy processes. Caps: 480 minutes for `S8`, 300 the rest. `S8` and `J5` on `held4` stay
as clean baselines for later steps.

**Stop rule:** stop and report if a stored count does not reproduce (`J5` totals 23,503,385,
22,939,812, 10,589,493, 9,041,686, 3,567,098; flagged and fits above), if prediction 1 fails, if
the gate fails, if the repair runs before `ar5/repair/frozen.json` is committed, if `held4`'s logs
or maps are read before its runs and repair end, if anything of `held2` is opened or listed, or if
a run hits its cap.

**Phase A (script layer).** 1. Link step 16's cache, the `S8` dirs and `J5L` read only; check the
counts; write only into this step's cache. 2. Tests first, `tests/test_ar5_repair.py`: ratio and
margin on hand-made maps; a voxel with no other pass end is never flagged; candidates as the pass;
rule `hit` with the current end first on ties; one MC fit per start; every evaluation counted; its
own generator leaves other draws alone; guards on `held2` and `held4`. 3. Commit `frozen.json`
(rule, gate, decision rules, hashes). 4. Flags (g0) by arithmetic. 5. Repair on the five maps and
`J5L` (g1, g2). 6. Build `held4` (pilot first), `S8` and `J5` with passes, the repair, 4 at a time.
Plot first (run evaluations against wrong voxels after the pass: `S8`, `J5`, `J5R`, six panels;
`mid` with `J5L` and `J5L` repaired); tables; predictions. "READY FOR KEY".

**Phase B (key).** One opt-in key `BoundaryPassRepair` (default 0, off), run at the end of
`run_boundary_pass`; its settings (0.95, 0.15, MC, generator) are module constants, as the pass's
are; parsed in `config_file.py`. Unit tests; key-off suite with no golden change; key on, a library
run on `mid` equal to `J5R`. Page `findings_ar5_repair.md` (plot first), history, pareto rows.
**Files:** `ar5/repair_*.py`, `tests/test_ar5_repair.py`, `benchmarks/autoresearch/ar5/repair/`
(under 3 MB tracked). **Library changes:** none in phase A; phase B the key, if the step passes.

**Not in this step.** Retrying steps 17 or 18, or `F95`, on `J5R` (each needs its own plan and, for
step 18, an owner decision). Other flag forms; two fits a voxel; a second repair round; a longer
CMA in place of MC; any change to the pass or to its skip (step 13's idea stays blocked); a change
to the house rule; realistic images; other crystals; a timing; another random seed; `held2`.

**Choices made.**
- This idea: it targets the one kind of voxel that costs `mid` its margin, with a test that needs
  no truth and the run's own MC fit; it leaves BFS growth alone, so the whole layer is exact and
  cheap. Sized by proxy at 0.7% to 1.7% of `J5`, 4 to 23 fixes a sample, `mid` from +4 to -6.
  Weak side: the proxy uses `S8`'s ends; the MC cost at contested voxels is assumed; the flag was
  picked looking at disclosed data; it saves no evaluations itself.
- `M(0.95, 0.15)` over `M(0.94, 0.10)` (0.32% to 1.10%, `mid` to -5): the wider form catches all
  four of step 18's breaks, the narrower catches one. Ratio alone costs 4.0% to 6.7%.
- Runner-up `F95` (seed MC stop at the first gain reaching 0.95): 1.34-3.31 / 1.73-4.59 / 1.63-4.54
  / 0.58-0.90 / 0.65-2.72%, but it moves hundreds of seed ends; best tried once `mid` has room.
  Inside the pass, an early stop of a try needs traces not stored; nothing else sizable.
- Owner option, not proposed: the house limit of 5 voxels is the same on every sample; `mid`'s 13
  broken and 9 fixed are 0.05% and 0.04% of its voxels. Keeping it means every step that touches
  growth must hold `mid` within one voxel; this step works around it instead.

**Amendment (D-20261115-N7), before the freeze and before any repair ran.** The sizing script
used its own neighbour relation and angle, not the pass's. The flag is frozen with the pass's own
`get_neighbors` and angle ("candidates as the pass"). The numbers of g0 and prediction 1 become:
flagged 194 / 183 / 142 / 238 / 56, fits 1,214 / 1,026 / 635 / 965 / 162 (`held3`, `held`, `mid`,
`grad`, `big`); `J5L` on `mid` 145 / 650. Nothing else changes.

### AR5 step 20: a cheaper repair: a 120-evaluation CMA fit in place of the MC fit

Task `ar5t`. Approved: D-20261115-N9. Starts from `feature/eval-reduction` after step 19's study merge
(D-20261115-N8). Reuses step 19's cache (worktree `agent-af5e5060b23027220`, `ar5/repair/cache`:
links to step 16's `J5` maps and pass tries, the `S8` maps, step 18's `J5L`), all read only.
Phase A ends at "READY FOR KEY" or at a gate; phase B only after a recorded decision.

**Why.** Step 19's repair worked on accuracy (fixed / broken against `J5`: 3/0, 9/0, 6/0, 8/0,
5/1; `mid` from +4 to -2 against `S8`) but failed its gate on cost: 3.14% to 7.09% of `J5`
against a cap of 2.5% (3.5% on `grad`), because an MC fit took a median of 662 evaluations, not
160. This step asks for the same repair at a fraction of the cost, keeping its fixes, so that
`mid` gains room against `S8` for later savings. It saves no evaluations itself.

**What was read and run.** PLAN step 19 with its amendment (D-20261115-N7); the gate report; step
19's repair scripts and cache (read only); the library's MC and CMA local optimizations. Run in
scratch: step 19's fits again with a trace of every gain, on the five `J5` maps and `J5L`, with
step 19's generator, with one generator per fit, and as 240-evaluation CMA runs; about 8 million
evaluations, under 10 minutes, 2 processes. Nothing of `held2` opened; `held4` not built.

**What the record says.** The trace reproduces step 19's fits exactly (evaluations, cost, hit
ratio per fit). Every MC fit runs its schedule to the end (its stop needs a cost under 0.0001,
never reached here). Only 33% of the 4,002 MC fits beat their start at all, and 96% make their
last gain by evaluation 100. A CMA fit beats its start in 99% to 100% of fits. MC draws its random
numbers from the one generator of the repair, three per trial step, so a fit stopped early shifts
every later fit's draws; CMA draws one seed per fit. The generator alone moves step 19's fixed
counts by up to 3 a sample. Reading line: per sample, `J5` wrong after the pass; the repair's cost
as a share of `J5` and fixed / broken against `J5`, for step 19 as run, for full MC with one
generator per fit, for the cheapest MC stop (stall 50) and for the chosen form; broken minus fixed
against `S8` for `J5`, step 19 and the chosen form.

| sample | `J5` wrong | step 19 | MC, own generators | MC stall 50 | CMA 120 | vs `S8`: `J5`, 19, CMA 120 |
|---|---|---|---|---|---|---|
| held3 | 54 | 3.59%, 3/0 | 3.63%, 5/0 | 0.35%, 5/0 | 0.62%, 14/0 | -7, -10, -21 |
| held | 67 | 3.14%, 9/0 | 3.17%, 6/0 | 0.30%, 6/0 | 0.54%, 18/0 | -3, -12, -21 |
| mid | 35 | 4.26%, 6/0 | 4.19%, 7/0 | 0.40%, 7/0 | 0.73%, 15/0 | +4, -2, -11 |
| grad | 56 | 7.09%, 8/0 | 7.08%, 6/0 | 0.63%, 5/0 | 1.29%, 31/0 | -7, -15, -38 |
| big | 21 | 3.47%, 5/1 | 3.39%, 5/0 | 0.29%, 5/0 | 0.55%, 14/0 | -3, -7, -17 |
| mid, `J5L` | 39 | 4.70%, 8/0 | 4.57%, 8/0 | 0.43%, 8/0 | 0.79%, 18/0 | +8, +0, -10 |

The chosen form was picked among about 50 forms (budgets, stall stops, win and lose stops, fewer
fits, longer CMA, tighter flags) looking at these samples, by a rule written after the first
table: see "Choices made". Median error falls by 0.00009 to 0.00037 degree with it.

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels. **Boundary pass**,
**candidates**, **ratio**, **margin**, **contested voxel** (`M(0.95, 0.15)`): as step 19, with the
pass's own `get_neighbors` (radius 2 side lengths) and angle (D-20261115-N7). **CMA fit**: one
local CMA run of the library (`optimizer="cma"`), first step 0.05 degree (the pass's try step), at
most 120 evaluations with the start's own, plus the final evaluation (121 in all). **Per-fit
generator**: `default_rng([0, v, c + 1])` for the fit at voxel `v` from start voxel `c` (-1 for
the voxel's own end). **Arm `J5C`**: `J5` with this repair. **Right**/**wrong**: at most / over 1
degree from the truth (symmetry-reduced). **Fixed**/**broken**: wrong in the older arm and right
in the newer / the reverse, counted by voxel.

**The repair, frozen now.** After the pass, on a snapshot of the map: flag the contested voxels
exactly as step 19's frozen `repair_rules` does (imported, not copied). For each, in index order:
one CMA fit from its own end and one from each candidate of the snapshot, each with its per-fit
generator; keep, among the current end and these ends, the highest hit ratio, then the lowest
cost, then the current end. A voxel's hit ratio is the stored one; a missing one costs one counted
evaluation. Every evaluation is counted. No stop, budget or fit depends on another fit, so no
draw moves when a fit is added, dropped or cut. BFS, seeds and the pass are untouched, so `J5C`'s
run is `J5`'s stored run plus the repair: exact, no rerun.

**Samples, and what has touched them.** `held3`, `held`, `mid`, `grad`, `big`: disclosed; step 19's
flag and this form were chosen on them. `J5L` (`mid`): seen. `held2`: not used, not opened. Clean
test: a new sample `held4`, built by `held_make_sample.py` exactly as `held` (250 Voronoi grains,
24,570 voxels, gradients 0 to 4 degrees) with generator seed 4, as step 16 built `held3`: a new
sample, not an extra random seed. It is still unbuilt and unread; nobody reads its logs or maps
before its runs and repair end. Clean images, Q-max 8.

**Measured.** Flagged voxels, fits, evaluations (and per fit), voxels changed; wrong voxels (Wilson
95%), fixed and broken against `J5` and `S8`, median error; `J5L` the same, and against `S8`.

**Predictions, written before any run.**
1. The tracked scripts reproduce the sizing exactly: flagged / fits 194 / 1,214, 183 / 1,026,
   142 / 635, 238 / 965, 56 / 162; repair evaluations 146,894 / 124,146 / 76,835 / 116,765 /
   19,602; fixed / broken 14/0, 18/0, 15/0, 31/0, 14/0; `J5L` 145 / 650, 78,650, 18/0.
2. `held4`: flagged 150 to 250; repair cost 0.4% to 1.0% of `J5`; fixed against `J5` 8 to 25,
   broken 0 to 2; median error not higher by more than 0.0001 degree.
3. `held4`, against `S8`: `J5` broken minus fixed -10 to +5; `J5C` -30 to -5.

**Decision rules, fixed now.** Gate, on the disclosed samples, before `held4` is built: (g0) step
19's flags and `J5` totals reproduce (above; `J5` totals 23,503,385, 22,939,812, 10,589,493,
9,041,686, 3,567,098); (g1) prediction 1, exactly; (g2) on each disclosed sample: repair
evaluations at most 1.5% of `J5` (2.0% on `grad`), fixed minus broken against `J5` at least 3,
and on `mid` broken minus fixed against `S8` at most +1. Gate fails: stop, report, no tuning, no
other form, no `held4`. Then `J5C` passes on a sample if, against `J5` (the same run, so draws
align by construction): (1) repair evaluations at most 1.5% of `J5`'s run with the pass (2.0% on
`grad`); (2) broken minus fixed at most 0; (3) median error at most 0.001 degree higher; (4) the
house rule against `S8` with the pass: at least 20% fewer evaluations, broken minus fixed at most
5, median error at most 0.005 degree higher. The step passes if `J5C` passes on `held4` and on the
five disclosed samples, and on `mid` broken minus fixed against `S8` is at most +1. Fails on
`held4`: no key value; study merge. `J5L` repaired is a measurement only.

**Machine.** Gate: the repair on five stored maps and `J5L`, about 0.56 million evaluations plus
setup, minutes, 2 processes. Then `held4` built (about 15 minutes); `S8` and `J5` on it with their
passes (about 67 and 23.5 million); the repair (about 0.15 million). In all about 91 million,
about 4.6 process-hours at 20 million an hour, about 3.6 hours wall (`S8` the long pole), at most
4 heavy processes. Caps: 480 minutes for `S8`, 300 the rest. `S8` and `J5` on `held4` stay as
clean baselines for later steps.

**Stop rule:** stop and report if a stored count does not reproduce, if prediction 1 fails, if the
gate fails, if the repair runs before `ar5/cmarepair/frozen.json` is committed, if `held4`'s logs
or maps are read before its runs and repair end, if anything of `held2` is opened or listed, or if
a run hits its cap.

**Phase A (script layer).** 1. Link step 19's cache read only; check the counts; write only into
this step's cache. 2. Tests first, `tests/test_ar5_cmarepair.py`: a fit's result does not depend
on which other fits ran or how long (per-fit generators); a CMA fit counts 121 evaluations; the
keep rule with the current end first on ties; the flag is step 19's frozen function; every
evaluation counted; guards on `held2` and `held4`. 3. Commit `frozen.json` (rule, gate, decision
rules, hashes). 4. Gate (g0 to g2); `J5L`. 5. Build `held4` (pilot first), `S8` and `J5` with
passes, the repair, 4 at a time. Plot first (run evaluations against wrong voxels after the pass:
`S8`, `J5`, `J5R` where stored, `J5C`, six panels; `mid` with `J5L`); tables; predictions.
"READY FOR KEY".

**Phase B (key).** One opt-in key `BoundaryPassRepair` (default 0, off), run at the end of
`run_boundary_pass`; its settings (0.95, 0.15, CMA, 0.05 degree, 120 evaluations, per-fit
generator) are module constants, as the pass's are; parsed in `config_file.py`. Unit tests;
key-off suite with no golden change; key on, a library run on `mid` equal to `J5C` (the same map
and evaluation count). Page `findings_ar5_cmarepair.md` (plot first), history, pareto rows.
**Files:** `ar5/cmarepair_*.py`, `tests/test_ar5_cmarepair.py`,
`benchmarks/autoresearch/ar5/cmarepair/` (under 3 MB tracked). **Library changes:** none in phase
A; phase B the key, if the step passes.

**Not in this step.** The runners-up below; an MC stop; retrying steps 17 or 18, or `F95`, on
`J5C` (each needs its own plan and, for step 18, an owner decision); a second repair round; any
change to the pass or to its skip (step 13's idea stays blocked); a change to the house rule;
realistic images; other crystals; a timing; another random seed; `held2`.

**Choices made.**
- This form, by the rule: on every disclosed sample, fixed at least step 19's minus one and broken
  at most step 19's; cost at most step 19's cap minus 0.5 point; median error at most 0.001
  degree higher; `mid` at most +1 against `S8`; among forms that change one thing in step 19's
  repair, the cheapest. Only CMA 120 and CMA 240 meet it; CMA 120 is half the cost. Weak side: it
  was chosen looking at the disclosed samples, by a rule written after the first table; its fixes
  are 2 to 5 times step 19's, which is consistent with MC's trial steps being coarse for these
  starts (not tested); it moves nearly every flagged voxel a little (the own-end fit refines it),
  so a break on `held4` would come from a candidate end beating a right end (none so far).
- Per-fit generators, not one stream: then a stop or a dropped fit moves no other draw, and every
  replay is exact. Cost: `J5C` does not use step 19's draws, so `J5C` against `J5R` mixes the
  form with the draws (up to 3 voxels a sample, as above).
- Runners-up: CMA 120 with only the top two candidates by pass-try hit ratio (0.28% to 0.92%;
  14/0, 15/0, 13/0, 24/0, 14/0; `mid` -9) and CMA 120 with flag `M(0.95, 0.10)` (0.29% to 1.00%;
  12/0, 14/0, 14/0, 29/0, 12/0; `mid` -10): about 40% cheaper again, but two changes, picked after
  seeing CMA 120. MC stall 50 (0.29% to 0.63%; 5/0, 6/0, 7/0, 5/0, 5/0; `mid` -3) keeps the MC fit
  and full MC's fixes but fails the rule on `held` (6 against 9 minus one). CMA 240 (1.08% to
  2.57%) buys 2 to 5 more fixes a sample for twice the cost.
- Owner option, not proposed: accept step 19's repair as it is (MC, 3.14% to 7.09% of `J5`, caps
  of about 4.3% on `mid` and 7.1% on `grad`); it fixes fewer than CMA 120 at six times the cost.

### AR5 step 21: score each pass candidate once, and let a wider repair do the fitting

Task `ar5u`. Approved: D-20261115-N11. Starts from `feature/eval-reduction` after step 20's key phase.
Reuses step 20's `cmarepair_*` (flag, repair, per-fit generators) and cache (`J5`'s maps before
and after the pass, pass tries, `J5C` maps, `S8` maps), read only. Phase A ends at "READY FOR KEY"
or at a gate; phase B only after a recorded decision.

**Why.** The owner asks for fewer evaluations at equal accuracy, and for early stops. The boundary
pass is 20% to 41% of `J5`'s run (`held4` 20.2%, `grad` 40.7%): one 60-evaluation CMA fit per
candidate per boundary voxel. Step 20's repair (a 120-evaluation CMA fit at contested voxels) now
sits behind it. This step stops every pass try at its first evaluation (a screen) and lets the
repair do the fitting, on a wider set of voxels. The pass runs after BFS, so growth is untouched:
a whole layer is `J5`'s stored run plus the new pass and repair. Every in-sample number is exact.

**What was read and run.** PLAN steps 17 to 20; DECISIONS N3 to N10; step 20's scripts and cache.
In scratch: all 385,943 stored pass tries of the six samples re-run with a trace of every gain
(exact: 0 differ), about 23.5 million evaluations, 7 to 9 minutes a run, 2 processes; then about
25 pass forms, each followed by step 20's repair run live. Form "full tries" reproduces `J5C`
exactly on all six samples. Nothing of `held2` opened. `held5` does not exist.

**What the record says.** Cutting the tries alone breaks voxels: at 30 evaluations 5 to 10 a
sample against `J5C`, at one evaluation 22 to 66. They are voxels that `J5`'s pass got wrong and
`J5C`'s repair fixed; with short tries the flag misses them, because its margin is measured
against lower try ends. A wider flag (below) gives them back, and then the try length hardly
matters. Reading line: per sample, the share of `J5C`'s total saved; fixed / broken against `J5C`
(by voxel); broken minus fixed against `S8` for `J5C` and for the chosen form (house limit 5).

| sample | tries cut to 30 | race (best on) | one evaluation | chosen | vs S8: J5C, chosen |
|---|---|---|---|---|---|
| held4 | 10.0%, 0/6 | 11.4%, 0/3 | 19.7%, 0/22 | 14.1%, 3/2 | -27, -28 |
| held3 | 10.3%, 0/8 | 11.4%, 0/2 | 20.5%, 0/27 | 14.4%, 4/2 | -21, -23 |
| held | 11.2%, 0/10 | 12.5%, 0/4 | 22.1%, 0/29 | 15.6%, 9/3 | -21, -27 |
| mid | 14.2%, 0/8 | 13.6%, 0/3 | 28.1%, 0/23 | 21.3%, 3/3 | -11, -11 |
| grad | 19.7%, 0/5 | 18.3%, 0/3 | 39.1%, 0/66 | 31.7%, 2/2 | -38, -38 |
| big | 15.6%, 0/8 | 9.5%, 0/2 | 30.8%, 0/29 | 26.6%, 3/1 | -17, -19 |

The first three columns keep step 20's flag. With the chosen form the median error falls by
0.0003 to 0.0031 degree; the repair grows to 4.6% to 8.6% of `J5`, the pass shrinks to 0.7% to
1.3%.

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels. **Boundary pass**,
**candidate**, **try**, **ratio**, **margin**, **contested** (`M(0.95, 0.15)`), **repair**: as step
20. **Screen try**: a try with a budget of one evaluation (`cma_max_evals=1`): the candidate start
is scored, with the library's final evaluation 2 counted; its end is the start. **Wide flag**: a
voxel is flagged if it is contested, or if its ratio is under 0.95 and one of its pass candidates
lies more than 1 degree from its kept end (a poor fit next to another orientation). **Arm `J5S`**:
`J5` with screen tries and the repair on the wide flag. **Right** / **wrong**: at most / over 1
degree from the truth (symmetry-reduced). **Fixed** / **broken**: wrong in the older arm and right
in the newer / the reverse, by voxel.

**The rule, frozen now.** BFS and seeds as `J5`. The pass as the library's, with every try a
screen try; the keep rule `hit` unchanged; each try's generator as `J5`'s pass (`make_rng(v,
"pass", k)`), though a one-evaluation try draws nothing that matters. Then step 20's repair,
imported unchanged (CMA 120, first step 0.05 degree, own end and every candidate of the snapshot,
per-fit generators, keep rule `hit`), on the wide flag. No voxel is skipped by its own hit ratio:
every boundary voxel scores all its candidates. This is not step 13's skip.

**Samples, and what has touched them.** `held4`, `held3`, `held`, `mid`, `grad`, `big`: disclosed;
the form was chosen on them. `held2`: not used, not opened. Clean test: a new sample `held5`,
built by `held_make_sample.py` exactly as `held` (250 Voronoi grains, 24,570 voxels, gradients 0
to 4 degrees) with generator seed 5, as steps 16 and 20 built `held3` and `held4`: a new sample,
not an extra random seed. Nobody reads its logs or maps before its runs and repairs end. Clean
images, Q-max 8.

**Measured.** Pass and repair evaluations, flagged voxels (contested and wide apart), fits, voxels
changed; wrong voxels (Wilson 95%), fixed and broken against `J5C`, `J5` and `S8`, median error.

**Predictions, written before any run.**
1. The tracked scripts reproduce the sizing exactly on (held4, held3, held, mid, grad, big): pass
   evaluations 181,534 / 161,864 / 170,478 / 100,058 / 120,710 / 37,242; flagged 1,869 / 1,787 /
   1,832 / 1,187 / 1,445 / 366; repair evaluations 1,617,286 / 1,512,258 / 1,547,106 / 759,275 /
   776,336 / 162,987; totals 23,612,288 / 20,240,655 / 19,457,817 / 8,397,057 / 6,257,077 /
   2,631,446; fixed / broken against `J5C` 3/2, 4/2, 9/3, 3/3, 2/2, 3/1.
2. `held5`: 11% to 18% fewer evaluations than `J5C`; flagged 1,500 to 2,200; repair 5% to 8% of
   `J5`; against `J5C` fixed 2 to 12, broken 0 to 4; median error 0.001 to 0.004 degree lower.
3. `held5`, against `S8`: `J5C` broken minus fixed -35 to -5; `J5S` within 3 of `J5C`.

**Decision rules, fixed now.** Gate, on the disclosed samples, before `held5` is built: (g0) the
traced re-run of `J5`'s tries and the full-tries replay reproduce `J5`'s tries and `J5C` exactly
(wrong voxels, repair evaluations 137,214 / 146,894 / 124,146 / 76,835 / 116,765 / 19,602); (g1)
prediction 1, exactly; (g2) on each disclosed sample: at least 10% fewer evaluations than `J5C`,
broken minus fixed against `J5C` at most 0, median error not higher, and on `mid` broken minus
fixed against `S8` at most +1. Gate fails: stop, report, no tuning, no other form, no `held5`.
Then `J5S` passes on a sample if, against `J5C` (the same `J5` run, so BFS is identical): (1) at
least 8% fewer evaluations; (2) broken minus fixed at most 2; (3) median error at most 0.001 degree
higher; (4) the house rule against `S8` with the pass: at least 20% fewer evaluations, broken minus
fixed at most 5, median error at most 0.005 degree higher. The step passes if `J5S` passes on
`held5` and on the six disclosed samples. Fails on `held5`: no key value; study merge. Rule (2)
allows 2, not 0, because `J5S` refits ten times more voxels than `J5C` and a clean sample can lose
one or two by chance; set now, before any `held5` number.

**Validity checks.** A screen try run by the library (`cma_max_evals=1`) equals the trace cut at
one evaluation, try by try, on `mid`; the arm reads no truth (the flag uses the grid's neighbours);
every evaluation counted, the repair's per-fit evaluations kept apart from the pass's; the wide
flag's extra voxels listed apart; guards on `held2` and on `held5` before its runs end.

**Machine.** Gate: the pass and repair on six stored maps, about 6 million evaluations, about 15
minutes on 2 processes; the library check on `mid` (about 0.9 million). Then `held5` built (about
15 minutes); `S8` and `J5` on it with their passes (about 68 and 27 million); then from `J5`'s
stored map `J5C` and `J5S` (about 0.15 and 1.8 million). In all about 104 million, about 5.2
process-hours at 20 million an hour, about 3.8 hours wall (`S8` the long pole), at most 3 heavy
processes. Caps: 480 minutes for `S8`, 300 the rest. `S8`, `J5`, `J5C` on `held5` stay as clean
baselines for later steps.

**Stop rule:** stop and report if a stored count does not reproduce, if prediction 1 fails, if the
gate fails, if `J5S` runs on a disclosed sample before `ar5/screen/frozen.json` is committed, if
`held5`'s logs or maps are read before its runs and repairs end, if anything of `held2` is opened
or listed, or if a run hits its cap.

**Phase A (script layer).** 1. Link step 20's cache read only; check the counts; write only into
this step's cache. 2. Tests first, `tests/test_ar5_screen.py`: a screen try costs 2 evaluations and
ends at its start; the wide flag on hand-made maps (ratio, 1 degree, the contested part unchanged);
the repair is step 20's (imported); per-fit generators; every evaluation counted; guards on `held2`
and `held5`. 3. Commit `frozen.json` (rule, gate, decision rules, hashes). 4. Gate (g0 to g2). 5.
Build `held5` (pilot first), `S8` and `J5` with passes, then `J5C` and `J5S`. Plot first (run
evaluations against wrong voxels after the repair: `S8`, `J5`, `J5C`, `J5S`, seven panels); tables;
predictions. "READY FOR KEY".

**Phase B (key).** One opt-in key `BoundaryPassScreen` (default 0, off), valid only with
`BoundaryPassRepair 1`: the pass's try budget becomes 1 and the repair takes the wide flag; its
settings (0.95, 1 degree) are module constants; parsed in `config_file.py`. Unit tests; key-off
suite with no golden change; key on, a library run on `mid` equal to `J5S` (the same map and
evaluation count). Page `findings_ar5_screen.md` (plot first), history, pareto rows.
**Files:** `ar5/screen_*.py`, `tests/test_ar5_screen.py`, `benchmarks/autoresearch/ar5/screen/`
(under 3 MB tracked; traces stay in the cache). **Library changes:** none in phase A; phase B the
key, if the step passes.

**What would make us stop the line.** The gate fails; `held5` saves under 8% or breaks more than
2 net against `J5C`; or the repair's cost on `held5` is above 10% of `J5` (then the wide flag does
not carry over, and the screen is not worth it).

**Not in this step.** Cutting tries to 8 or 15 evaluations; races; the top-2 repair; a learned
flag; a BFS skip or seed stop on top (later steps, with `J5S` as their base); counting a screen try
as 1 evaluation (about 0.3% to 0.7% of `J5C` more); the pass skip (step 13); the house rule;
realistic images; other crystals; a timing; another random seed; `held2`.

**Choices made.**
- This idea: the largest exact saving found (14.1% to 31.7% of `J5C`), no growth risk, no model,
  in-sample no net loss on any sample and `mid` at -11 as with `J5C`. Weak side: picked among
  about 25 forms looking at all six disclosed samples, so `held5` is the only clean test; the
  repair becomes the main cost and scales with boundary area; 4 of 13 broken voxels are not
  flagged at all; the accuracy gain is consistent with more voxels getting a 120-evaluation fit
  from their own end (not tested).
- One evaluation, not 8 or 15: equal fixed / broken on five samples (`grad` 2/2 against 2/1 at
  15), 4.6 to 9.2 points more saving, and the simplest rule (a budget of 1 in the library).
- Runners-up: the race (best try runs on; a general early stop): 9.5% to 18.3%, but 2 to 4 broken
  and none fixed against `J5C`. Step 17's BFS skip behind the repair (3.2% to 8.7% of `J5` in
  replay; needs live runs; stacks with this step later). `F95` (0.6% to 4.6%; needs a recording).
  A learned grid search (ceiling 0% to 14.5%; a wrong seed grows a grain the repair cannot fix).
- Not proposed: step 18's learned skip again (needs an owner decision, step 20's plan text).

### AR5 step 22: defer a stalled seed refinement while untried candidates remain

Task `ar5v`. Approved: D-20261116-N2. Starts from `feature/eval-reduction` after step 21's key phase.
Reuses step 16's `high500_*` and step 8's aligned draws (`order_draws.py`), step 21's `screen_*`
and cache (`J5` runs and logs of `held5`, the `J5S` maps), step 16's and step 20's `J5` runs, all
read only. Phase A ends at "READY FOR KEY" or at a gate; phase B only after a recorded decision.

**Why.** The owner asks for fewer evaluations at equal accuracy, and for general early stops that
predict whether more optimization will help. After step 21 the seed step is the largest part of
`J5S` on six of seven samples: 43% to 59% (`grad` 22%, `big` 17%), almost all of it MC
refinements. Most of that goes to refinements that never reach 0.95. Under rule `PH(500)` such a
refinement, once stalled, keeps going to its budget while the step still has untried candidates.
This step stops paying for that: a stalled refinement waits until the other candidates are tried.

**What was read and run.** PLAN steps 17 to 21; `plan_step21_candidates.md`; step 21's tables and
`held5` report; `J5`'s `run.json`, `polish.json`, `e1_log` and `mc_traces` on all seven samples.
Arithmetic only, in scratch; no run. Nothing of `held2`.

**What the record says.** Seed-step MC refinements that end at 0.95 or more hold only 1.6% to 6.7%
of `J5S`; the tail after the last gain in refinements `PH(500)` did not stop holds 7.5% to 29.5%.
Reading line: per sample, `J5S`'s total; the seed step's share (MC in brackets); steps of type A
(exact, below) and the exact saving as a share of `J5S`; steps of type C (path changes); steps of
type N (no candidate at 0.95); the bounds on the net saving (C's extra cost from 201 evaluations a
step up to (4 - k) mean paused candidates).

| sample | J5S total | seed step (MC) | A steps, exact saving | C steps | N steps | net |
|---|---|---|---|---|---|---|
| held5 | 22,706,889 | 48.2% (41.3%) | 998, 10.54% | 191 | 546 | 9.1-11.7% |
| held4 | 23,612,288 | 58.5% (50.9%) | 1,106, 13.14% | 251 | 593 | 11.7-15.1% |
| held3 | 20,240,655 | 53.1% (45.5%) | 869, 9.75% | 172 | 599 | 8.4-11.0% |
| held | 19,457,817 | 58.8% (50.6%) | 1,010, 12.28% | 238 | 512 | 10.2-14.1% |
| mid | 8,397,057 | 43.5% (36.7%) | 358, 10.65% | 79 | 127 | 8.8-11.9% |
| grad | 6,257,077 | 22.0% (15.1%) | 269, 2.40% | 10 | 98 | 2.0-2.4% |
| big | 2,631,446 | 17.0% (13.5%) | 62, 4.95% | 9 | 4 | 4.5-5.6% |

The rest of `J5S`: grid searches 0% to 22.8% (0 to 15 searches), BFS local optimizations 20% to
75%, the screen pass 0.75% to 1.93%, the repair 6.2% to 12.4%. A pure stop at the check saves
23% to 31% on the `held` types but changes about 400 seed steps a sample.

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels; here always the
voxel's **own** hit ratio, never a ratio to its neighbours. **Seed step**: index the voxel, cost
20 candidates, then for the best 5 in turn a CMA fit and, under 0.95, an MC refinement; it stops
at the first end at 0.95 or more, else takes the lowest-cost end. **Gain**: a strict drop of the
running best cost. **Check**: `J5`'s `PH(500)` check, one counted evaluation of the running best's
own hit ratio after 500 evaluations with no gain; at 0.95 the refinement stops, else it **goes
on**. **Pause** / **resume**: keep a refinement's state and generator / continue it exactly where
it stopped. **Step types**: **A**, found at a later candidate than every paused refinement; **N**,
no end at 0.95; **C**, the found refinement itself would be paused. **Arm `J5D`**: `J5` with rule
`DEF`; **`J5SD`**: `J5D` with step 21's screen pass and wide repair. **Right** / **wrong**: at
most / over 1 degree from the truth (symmetry-reduced). **Fixed** / **broken**: by voxel.

**The rule, frozen now.** Rule `DEF`: at a check that goes on, if the seed step has an untried
candidate, the refinement is paused, and the next candidate is tried. The step stops at the first
end at 0.95 or more, as now. If no later candidate reaches 0.95, the paused refinements are resumed
in candidate order under `J5`'s rules until one reaches 0.95; else the lowest-cost end is taken,
as now. The last candidate is never paused. No new check, no threshold, no neighbour, no model.
CMA fits, 5 candidates, 0.80, 0.95, `PH(500)`, `G1`, `G2b`, the grid search, BFS, the screen pass
and the wide repair are unchanged. Type A and N steps are exactly `J5`'s; only C steps change.

**Samples, and what has touched them.** `held5`, `held4`, `held3`, `held`, `mid`, `grad`, `big`:
disclosed; the rule was chosen on them. `held2`: not used, not opened. Clean test: a new sample
`held6`, built by `held_make_sample.py` exactly as `held` (250 Voronoi grains, 24,570 voxels,
gradients 0 to 4 degrees) with generator seed 6, as earlier steps built `held3` to `held5`: a new
sample, not an extra random seed. Nobody reads its logs or maps before its runs and repairs end.
Clean images, Q-max 8.

**Measured.** Seed-step evaluations, pauses, resumes, steps A / N / C, C steps whose taken end
moves; grid searches, rescues; evaluations by stage; wrong voxels (Wilson 95%); fixed and broken
against `J5S`, `J5` and `S8`; median error over all 24,570 voxels.

**Predictions, written before any run.**
1. The tracked replay reproduces the table exactly on the seven disclosed samples: A steps, exact
   savings, C steps and N steps as above.
2. Live `J5SD` against `J5S`, disclosed samples: savings inside the net bounds widened by 1.5
   points; 30% to 70% of C steps end at another candidate; broken minus fixed -3 to +3; grid
   searches within 2 of `J5`'s.
3. `held6`: 8% to 14% fewer evaluations than `J5S`; broken minus fixed against `J5S` -3 to +3;
   median error within 0.0005 degree; the house rule against `S8` passes.

**Decision rules, fixed now.** Gate, on the disclosed samples, before `held6` is built: (g0) with
`DEF` off the scripts reproduce `J5`'s seed steps on `mid` (ends, evaluations); a resumed refinement
equals `J5`'s unpaused one, evaluation by evaluation; (g1) prediction 1, exactly; (g2) the live
`J5D` run on `mid` equals the replay in every seed step up to the first C step in run order; (g3)
`J5SD` on `mid` and `held`: at least 7% fewer evaluations than `J5S`, broken minus fixed at most
3, median error at most 0.001 degree higher, grid searches at most 2 more than `J5`'s. Gate fails:
stop, report, no tuning, no threshold, no `held6`. Then `J5SD` passes on a sample if, against
`J5S`: (1) at least 7% fewer evaluations on the `held` types and `mid`, 1.5% on `grad`, 3% on
`big`; (2) broken minus fixed at most 3; (3) median error at most 0.001 degree higher; (4) the
house rule against `S8`: at least 20% fewer evaluations, broken minus fixed at most 5, median
error at most 0.005 degree higher. The step passes if `J5SD` passes on `held6` and on the seven
disclosed samples. Fails on `held6`: no key value; study merge. Rule (2) allows 3, not 0: C steps
move seed ends, and a moved end changes BFS growth (route noise); set now, before any run.

**Validity checks.** Each refinement's generator is keyed by voxel, kind and call count; a resume
continues the same generator and is not a new call (checked on `mid`). Phase A may resume by an
exact re-run of the paused part: those evaluations are checked identical to the first run and not
counted, and are listed apart. Every other evaluation counted, checks included. The arm reads no
truth. C steps listed with old and new taken ends. Guards on `held2`, and on `held6`.

**Machine.** Gate: replay (minutes); live `J5D` on `mid` and `held` (about 24 million), each then
the screen pass and wide repair (about 2.5 million). Then `J5SD` live on `held5`, `held4`,
`held3`, `grad`, `big` (about 68 million). Then `held6` built (about 15 minutes); `S8` and `J5`
with passes (about 70 and 27 million); `J5S` from `J5`'s map (about 1.7 million); `J5SD` live
(about 22 million). In all about 215 million, about 10.8 process-hours at 20 million an hour,
about 5 hours wall at 3 heavy processes. Caps: 480 minutes for `S8`, 300 the rest. `S8`, `J5`,
`J5S` on `held6` stay as clean baselines.

**Stop rule:** stop and report if a stored count does not reproduce, if prediction 1 fails, if the
gate fails, if `J5D` runs on any sample before `ar5/defer/frozen.json` is committed, if `held6`'s
logs or maps are read before its runs and repairs end, if anything of `held2` is opened or listed,
or if a run hits its cap.

**Phase A (script layer).** 1. Link step 21's and step 16's caches read only; check the counts;
write only into this step's cache. 2. Tests first, `tests/test_ar5_defer.py`: a pause only at a
check that goes on and never at the last candidate; resume order and exactness; type A and N steps
equal `J5`'s; generators unchanged; every evaluation counted; guards on `held2` and `held6`.
3. Commit `frozen.json` (rule, gate, decision rules, hashes). 4. Replay (g0, g1). 5. Gate (g2,
g3). 6. `J5SD` on the other five. 7. Build `held6` (pilot first), `S8` and `J5` with passes,
`J5S`, `J5SD`. Plot first (run evaluations against wrong voxels: `S8`, `J5S`, `J5SD`, eight panels);
tables; predictions. "READY FOR KEY".

**Phase B (key).** One opt-in key `SeedIndexMCDefer` (default 0, off), valid only with
`SeedIndexMCStallHit` on: the library's MC refinement gets a true pause and resume (its state and
generator kept); parsed in `config_file.py`. Unit tests; key-off suite with no golden change; key
on, a library run on `mid` equal to `J5D` (the same map and evaluation count). Page
`findings_ar5_defer.md` (plot first), history, pareto rows. **Files:** `ar5/defer_*.py`,
`tests/test_ar5_defer.py`, `benchmarks/autoresearch/ar5/defer/` (under 3 MB tracked).
**Library changes:** none in phase A; phase B the key, if the step passes.

**What would make us stop the line.** The gate fails; `held6` saves under 7% against `J5S` or
breaks more than 3 net; or extra grid searches (about 345,000 evaluations each) eat the saving.

**Not in this step.** A pure stop at the check; a go-on threshold on the own hit ratio; `F95`; a
stall stop in the repair's CMA fits; step 17's BFS skip; learned proposals for seeds or the grid
search; the boundary-pass skip (step 13); the house rule; realistic images; other crystals; a
timing; another random seed; `held2`.

**Choices made.**
- This idea: the largest saving with an exact core (9.75% to 13.14% of `J5S` on the `held` types,
  10.65% on `mid`), a general early stop (do not polish a stalled answer while untried ones
  remain), no model and no threshold. Weak side: 172 to 251 C steps a `held` sample move seed
  ends, so growth changes and only live runs judge it; small on `grad` (2.0% to 2.4%) and `big`;
  chosen looking at seven disclosed samples, so `held6` is the only clean test; about 10.8
  process-hours.
- Runners-up: a pure stop at the check (23% to 31%, about 400 changed seed steps a sample); a skip
  of the MC by the CMA end's own hit ratio (3.0% to 3.7% at 0.3 with 25 to 46 wrong skips); step
  17's BFS skip behind the repair (4.3% to 5.0%, `grad` 12.6%; live); `F95` (0.8% to 5.4%; needs
  a recording); a stall stop in the repair (1.4% to 2.9% at k = 40; the own hit ratio still rises
  late in the fits); learned grid proposals (ceiling 0% to 22.8%; step 10's cut failed). Not
  proposed: step 18's learned skip again; the pass skip by a voxel's own hit ratio (step 13).

### Owner task of 16 Nov, part 1: check the C++ pilot result for a bug

Task `cppchk`. Owner task given 2026-11-16 (part 1 of 2; part 2, task `fine1`, goes ahead whatever
this part finds). Starts from `feature/eval-reduction` after AR5 step 22. A check, not a study:
it names the cause of what the stored C++ run shows, or reports that no bug was found.

**Why.** The C++ time study (`cpptime`, merged at 003cc5d, D-20261114-N3) left 416 of 2,000 pilot
voxels wrong in C++ against 24 in Python. Two things in its record look like a bug: 393 voxels
marked for a refit were never fitted again, and all 1,921 C++ local optimizations print an end
orientation equal to the start. C++ is the reference program, so the owner wants this checked
before C++ numbers are used again (part 2 runs C++ on a new layer).

**What was read.** PLAN section "C++ time on the pilot"; `cpptime/tables.md`, `trace_C1.json`;
D-20261114-N3; the C++ sources named below; Python `config_file.py`, `orientation_search.py`,
`reconstructor.py`. The run cache (`cpptime/cache/runs/C1/`: client logs, stdout, `.mic.LBFS`)
no longer exists: `finish_task.sh` removed the worktree that held it (record of 2026-11-14 10:38).
Only the tracked summaries remain. So phase A starts with a rerun of the stored command.

**What the record says.** Reading line: each row is one count from `trace_C1.json` or
`tables.md` of run C1 (one fitting client, stored binary of 2026-03-23, Release).

| count in run C1 | value |
|---|---|
| wrong voxels (over 1 degree from the truth) | 416 of 2,000 |
| voxels marked REFIT (hit ratio over the running best at most 0.9), all wrong | 393 |
| voxels accepted above 0.9 that are wrong | 23 |
| local-optimization records whose printed end equals the printed "Original" | 1,921 of 1,921 |
| distinct orientations in the output map | 80 |
| wrong voxels within 1 degree of another grain's truth | 416 |
| paired with Python A1: wrong only in C++ / only in Python / both / neither | 404 / 12 / 12 / 1,572 |
| seconds in global searches / in local optimizations (C1) | 1,137 / 85 |

**Terms.** **Global search**: a search with no starting guess (C++ `ReconstructVoxel`).
**Local optimization**: the MC fit from a neighbour's orientation (C++ `LocalOptimization`).
**Cost**: 1 minus the overlap quality; quality is the pixel hit ratio times the detector ratio
(`OverlapInfo.h:138,164`). **REFIT**: the mark C++ puts on a voxel whose hit ratio is at most 0.9
of the running best of its BFS. **Server** / **client**: MPI rank 0 (hands out seeds) and rank 1
(fits). **Wrong**: over 1 degree from the truth, cubic symmetry. **Bug**: code that does not do
what its own comments, names or config keys say, or two programs given different inputs, or a
scoring script that misreads a file. A behaviour the code was written to have is a design
difference, reported as such, not a bug.

**Hypotheses, to test, not concluded.** Reading line: per check, the candidate cause, where it
lives, and the observation that would confirm or refute it.

| # | candidate cause | code | confirms / refutes |
|---|---|---|---|
| R1 | client grid marks never cleared between seeds: `Reset` empties the queue only, though its comment says visited points "have to be reset to NON_VISITED"; `InsertSeed` queues only NOT_VISITED voxels | `Src/ReconstructionStrategies.tmpl.cpp:314-326`, `:335-356`; called at `BreadthFirstReconstructor.tmpl.cpp:115` | a REFIT voxel borders a voxel fitted by a later BFS and is never queued / no REFIT voxel ever borders a later BFS |
| R2 | server never re-serves REFIT voxels: `Pop` skips all but NOT_VISITED; `Refit` runs only on RESTART_FIT, sent only with `PartialResultFilename` | `ReconstructionStrategies.tmpl.cpp:124-145`, `:151-196`; `BreadthFirstReconstructor.tmpl.cpp:84-105`, `:223-228`, `:361-366`; `ConfigFile.cpp:689-693` | server log: queue empty while REFIT voxels remain, no RESTART_FIT sent / any REFIT voxel served again |
| L1 | start already at the lowest reachable cost; MC accepts only a strictly lower cost | `OrientationSearch.cpp:174`, `:241-261` | stdout `C_i` equals `C_f` and is near 0 at accepted voxels / `C_f` below `C_i` |
| L2 | stop after the first 10-step block: cost under `MaxConvergenceCost` (0.0001) and variance under 0.02 squared | `OrientationSearch.cpp:228-231`, `:263-264`; `ContinuousSearch.h:170-186` | stdout `Step` = 10 on most records / larger step counts |
| L3 | log reading: "Original" is the BFS seed (`oResult`), not the voxel's own start; Euler angles print with 6 significant digits | `BreadthFirstReconstructor.tmpl.cpp:174-184` | stdout misorientation to the start nonzero while the log shows equal / exactly 0 |
| L4 | settings differ from Python: box width = `LocalOrientationGridRadius` / 2^`MinLocalResolution`; 0.5 degree resolution fixed; `MCRadiusScaleFactor` unused; `SuccessiveRestarts` passed but unused | `ContinuousSearch.h:68-73`, `:122-133`, `:177`; `SearchDetails.cpp:46-57`; `ConfigFile.cpp:489-520`; Python `config_file.py:460-470`, `orientation_search.py:362`, `:1015-1110`, `reconstructor.py:1510-1521` | parsed values or the use of a key differ between the two / equal |
| L5 | the C++ local cost differs from Python's at the same orientation (`ShapePixelOverlapCounter`, `DiscreteAdaptive.h:119-127`) | `DiscreteAdaptive.tmpl.cpp:280-317` | per-orientation cost differs by over 0.006 (the old parity limit) / within it |
| I1 | inputs differ (config, detector, images, grid, structure) | `cpptime/make_config.py`; `pilot_bfs.py:76` | any hash or parsed value differs / all equal |
| S1 | the scoring script misreads the C++ map (columns, degrees, Euler convention, order, symmetry) | `cpptime/score_cpp.py:35-62`, `:76-86`; `XDMServer.cpp:84-133`, `:219-222` | an independent reading disagrees on any voxel / agrees on all |

**Samples.** The pilot only (2,000 voxels, `benchmarks/phase_d_pilot/region_grid.mic`, clean
images, Q-max 8). Python side: item 14's stored run `item14/runs/A1` (`recon.mic`,
`records.json`), read only, not rerun. The 10 voxels: the first 10 in pilot-grid row order that are
wrong in C++ and right in Python (from C3's map, scored by `score_cpp.py`). Fixed now, before any
cost is looked at. The traced local optimization: the first record in the C3 client log, in log
order, whose voxel is one of the 10.

**Measured.** (1) Per REFIT voxel: the BFS that marked it, the later BFS runs that reached a
neighbour of it, whether it was queued again. (2) Per local optimization (stdout of C3): start
cost `C_i`, end cost `C_f`, steps, restarts, misorientation end to start. (3) Per chosen voxel:
the cost in C++ and in Python at four orientations (truth, Python end, C++ end, C++ start).
(4) The step-by-step trace: every MC sample (offset, orientation, cost, accepted or not), block
radius, block variance. (5) Inputs: hashes and parsed values. (6) Scoring: per-voxel errors by
two independent readings.

**Predictions, written before any run.**
1. C3 (rerun of the stored command) reproduces every count of the table above exactly.
2. R1 and R2 both hold: every one of the 393 REFIT voxels is adjacent to at least one voxel
   fitted later, and none is queued or served again.
3. On at least 90% of accepted records `C_i` equals `C_f`; L3 explains part of "all equal".
4. L4: the parsed settings agree; the two ports share the same MC rule.
5. L5: costs agree within 0.006 on all 40 voxel-orientation pairs.
6. I1 and S1: no difference found.

**Decision rules, fixed now.** A hypothesis is confirmed only by its listed observation, on the
rerun's own logs. Each check ends as "bug found" (with file and line), "design difference" or
"no cause found". If a bug is found: a fix on its own task branch `feature/eval-reduction-cppfix`,
with a test (a C++ run on a small sub-grid, read by a pytest), an Opus review, and then a rerun of
the pilot (C4, untimed, same command); both C1/C3 and C4 are kept and reported side by side. One
fix per bug; no tuning to match Python. If no bug is found: the page says which hypotheses were
refuted, which design differences explain the 416 (with counts), and that C++ stays as it is.

**Validity checks.** C3 equals the stored counts, else stop. A Release build from current `Src/`
(with the evaluation counter of b9b566b) gives a `.mic.LBFS` byte-equal to C3's, else the stored
binary stays the reference. The trace build's map is byte-equal to C3's (prints change no draw).
C++ costs at given orientations come from a restart run (`PartialResultFilename` on a 10-voxel
sub-grid: `Refit` prints the start cost before any move); if that old path fails, from a small
driver copied from `Src/CostFunctionBenchmark.cpp`. Scoring: Python A1's map read through
`score_cpp.py` gives its stored 24 wrong; a second reader (raw text, `scipy` rotations) agrees
voxel by voxel. Counts and paired splits only.

**Machine.** C3: about 21 minutes, 2 busy processes, untimed. Builds: two Release builds, a few
minutes each. Trace run on the pilot: about 21 minutes. Restart runs: 4 short runs, seconds each.
Scripts: minutes. If a fix: C4, about 21 minutes. In all under 2 hours of machine.
**Process limit exception:** the C++ build command is fixed (`make -j8`, 8 compile jobs); it runs
only when no other heavy process runs, for its few minutes.

**Stop rule:** stop and report if C3 does not reproduce the stored counts, if a run passes 1 hour
or its log does not grow for 10 minutes, if stderr holds `RUNTIME_ASSERT`, if any file of an
earlier step would need an edit, or if a fix would change Python. Part 2 starts in any case.

**Phase A (scripts, task branch `feature/eval-reduction-cppchk`).** 1. `read_on_fetch.sh --save`.
The cache lives outside the worktree (a gitignored folder of the main checkout, linked), so
`finish_task.sh` cannot delete it again. 2. Tests first (`tests/test_cppchk.py`): log parsers on
hand-made lines, the 10-voxel rule, the second reader on a hand-made map. 3. C3 with the stored
binary and `cpptime/make_config.py`; stdout kept. 4. Checks I1 and S1. 5. R1, R2, L1 to L3 from
C3's logs. 6. L4: both parsers on the same config; the C++ use of each key read from the source.
**Phase B (C++, task branch `feature/eval-reduction-cpptrace`).** Trace prints behind an
environment variable naming one voxel location (`OrientationSearch.cpp:143-198`, `:207-287`) and a
one-line evaluation count at client exit (near `XDMClient.cpp:103`, also needed by part 2).
Release build only: `cmake -DCMAKE_BUILD_TYPE=Release . && make -j8`. Opus review before use.
Trace run; restart runs for L5; Python costs at the same orientations. **Phase C.** Verdict per
check; if a bug, the fix branch, test, review, C4. Plot first (per check, counts); tables; page.

**Files.** `scripts/autoresearch/cppchk/`, `tests/test_cppchk.py`,
`benchmarks/autoresearch/cppchk/` (summaries under 1 MB tracked; logs and maps in the cache).
**Library changes:** none in Python. C++: the trace and count prints on `cpptrace`; a fix only on
`cppfix`. No frozen file of an earlier step is edited. Nothing merges to develop.

**Not in this step.** C++ with several clients; a Python arm with C++ rules; timing; another
sample; changing the C++ search design (a revisit pass) unless a bug is found; another seed.

**Choices made.**
- Rerun first, because the stored logs are gone; the run is deterministic (one client), so C3
  stands in for C1 if the counts match.
- Costs from the restart path before a new driver: no code change for that part.
- One count print added at client exit: the open choice of N3 ("a C++ evaluation count by a print
  and a rebuild"); part 2 needs it. Owner question: allow it, or report C++ evaluations as missing.
- R1 is the lead candidate for the dropped refits (comment and code disagree); it may still turn
  out to be intended for one client, so the decision rule above decides, not this text.

### Owner task of 16 Nov, part 2: a whole layer one generation finer

Task `fine1`. Owner task given 2026-11-16 (part 2 of 2). Goes ahead whatever part 1 (`cppchk`)
finds. Starts from `feature/eval-reduction` after AR5 step 22. A measurement on a new sample: no
method is changed or tuned here.

**Why.** Every whole-layer number so far is on 24,570-voxel layers. The owner wants the three
methods on a layer with half the voxel side (four times the voxels), recorded as evaluations,
wrong voxels and wall-clock time, on the cost-against-accuracy chart of the review page.

**What was read.** `ar3/held_make_sample.py`, `ar3/held_sample_lib.py`, `gradient_sample.py`,
`ar5/screen_build_held5.py`; `phase_d/render_full.py` (arguments); the `cpptime` plan and scripts;
step 21's `screen/lib_check.json` and `summary.json`; item 14's `runs/A1/run.json`;
`pareto/pareto_rows.json`; RUNBOOK "The review page"; C++ `ReconstructionStrategies.h:66-80`.

**Terms.** **Generation**: the triangle level of a `.mic` grid; side = 0.6 mm / 2^generation.
`held` is generation 6 (side 0.009375 mm, 24,570 voxels). **Fine layer**: generation 7, each
triangle split into four (side 0.0046875 mm, 98,280 voxels, the same area). **Neighbour pair**:
two voxels of one grain that share an edge. **Spread** of a grain: the largest misorientation
between one of its voxels and its centre orientation. **Global search**: a search with no starting
guess. **Local optimization**: a fit from a neighbour's orientation. **Wrong**: over 1 degree from
the truth, cubic symmetry; a voxel with no output orientation counts as wrong and is listed apart.
**Quiet**: the timing preflight passes (`timing_preflight.py --require`) and the `quiet.py`
sampler (`item14/quiet.py`, rule R4) finds no other busy process. **Slice**: the 2,000 voxels of
the fine layer nearest the layer centre (as `render_full.py --n-voxels`), fixed now.

**The arms.** Reading line: each row is one run on the whole fine layer and its exact settings.

| arm | what it is | settings |
|---|---|---|
| `PS` stock Python | grid search for every seed, stock MC local optimization: the Python run `cpptime` set beside C++ (item 14 arm A) | `ar2/key_off_fused.config` (with `BFSRevisitRefit 1`, `BFSReseed 1`, `FusedCost 1`) |
| `CS` stock C++ | mode `p`, `mpirun -n 2` (one server, one fitting client), as `cpptime` | `cpptime/make_config.py` on the same config (drops the three Python-only lines) |
| `J5S` latest recipe | AR5 step 21's arm, run fully from library keys (`lib_check.json`: library equals the script layer on `mid`) | `ar2/key_on_fused.config` plus `BFSNeighborCMAFirst 1`, `CMASigma0 0.05`, `CMANeighborMaxEvals 60`, `BoundaryPass 1`, `SeedIndexCandidates 20 5`, `SeedIndexStopHitRatio 0.95`, `BFSNeighborCMANoFallback 1`, `SeedIndexCMAFirst 200`, `SeedRescueCMA 0.25 1000`, `SeedIndexMCStallHit 500`, `BoundaryPassRepair 1`, `BoundaryPassScreen 1` |
| `CF` (only if part 1 makes a fix) | `CS` with part 1's fixed build, an extra labelled run | as `CS` |

Every config changes only: `SampleFilename` (fine grid), `InfileBasename` (fine images),
`MinSideLength` and `MaxInitSideLength` 0.0046875 (C++ splits any voxel above `MinSideLength`;
Python does not read them), `OutStructureBasename`. Clean images, Q-max 8 as before.

**Sample.** Name `fine7`. New script `fine1_make_grid.py` splits every triangle of
`SimInput/rand_500grains_1mm_neworient_s0_grid.mic` into four (three corner triangles and the
flipped middle one). New builder `fine1_make_sample.py`: `held_make_sample.py`'s steps with its
helpers imported unchanged (250 Voronoi cells, connected grains, orientations redrawn if under
5 degrees from a neighbour grain, gradient levels 0, 0.2, 0.5, 1, 2, 4 degrees, linear model),
on the fine grid, generator seed 7. A new seed is a new sample, as `held3` to `held6` were built;
it is not an extra random seed for a comparison. One addition, the gradient cap: after a grain's
gradient is drawn, if its largest neighbour-pair misorientation is 0.45 degree or more, its
rotation vectors are scaled down so that it is 0.45 degree; each capped grain is listed with its
drawn and kept level. A crude estimate from `held`'s grain sizes (not a measurement): about 2 of
250 grains need the cap at the fine grid, about 45 would at the coarse one. Images: `render_full.py
--tag fine7 --variants clean --mic <truth>` (the render config of `held`), then the stack.

**Gradient check, a hard gate before the freeze.** `fine1_gradient_check.py` reads the written
truth back from the `.mic` and reports, over all neighbour pairs: the count, the maximum, the
median, the 90th, 99th and 99.9th percentiles and a histogram of neighbour misorientation; per
grain: spread, level, voxels, largest neighbour step. Gate: the maximum is under 0.5 degree, and
the largest spread is at least 3 degrees (the "several degrees" of the spec). Fails: stop.

**Freeze.** Committed before any method runs: `fine1/frozen.json` (seed, grid rule, cap rule,
arms, slice, projection formula, decision rules, sha256 of the grid, truth, slice grid, every
image frame and the stack), the generator scripts, `sample.json`, `grainmap.npz`, the gradient
check JSON. The truth `.mic` (about 10 MB) and images are not tracked, as `held`'s are not; their
hashes are. Tracked files under 3 MB.

**Measured.** Per run: evaluations (all, and global searches apart); wrong voxels with a Wilson
95% interval (`common/stats.py`); missing voxels; median error; wall-clock time, labelled with the
machine, the process count, what else ran and the quiet verdict. C++ evaluations come from part
1's count print; without it they are reported as missing.

**Budget: slices first, then a projection.** Each arm runs once on the slice, timed alone. The
projection: T = t_g x G + t_r x (98,280 / 2,000), with t_g the slice's mean seconds per global
search, G the projected global searches and t_r the slice's seconds outside global searches. Low
G = slice searches x (grains in the layer / grains touching the slice), from grain-map counts
only; high G = slice searches x 49.14. Weak points: a slice holds only about 5 to 10 grains, so
G is poorly known; seeds and reseeds per voxel and grid searches do not scale linearly; smaller
spots may change hit ratios and so rejected seeds; Python per-voxel overheads may grow with the
layer. From the slice runs only counts and seconds are read, never maps or wrong counts.

**Stored-number estimates (not the projection).** Reading line: per arm, the stored run used,
evaluations on the fine layer (low: global part fixed per grain; high: all times four), the
rate used, and hours on one worker.

| arm | stored basis | evaluations, low-high | rate (per process-hour) | hours |
|---|---|---|---|---|
| `PS` | pilot A1: 69 global searches at about 332,000 each, 1,349 per voxel otherwise | 228 M - 1.26 G | 57 M (A1, near-quiet) | 4 - 22 |
| `CS` | C1: 85 searches at 13.4 s on 60 grains; 1,921 local opt. in 85 s | not measured | seconds | 2.5 - 17 |
| `J5S` | `held5`: 22,706,889, seed step 48.2% | 58 M - 91 M | 10.6 M (`J5` on `held5`) | 5.5 - 8.6 |

Not verified: the global-search count on the fine layer (scaled from 60 pilot grains to 250); the
`PS` high end assumes reseeds grow with voxels; `CS` local optimizations stay at 44 ms only while
they return their start (a part 1 fix may change that); rates were measured contended. For
comparison `S8` on `held5`: 75,879,690 evaluations, 10,135 s contended for its search part.

**Decision rules, fixed now.** If any arm's high projection passes 48 hours, stop and ask the owner
before starting it. Caps: 1.5 times each arm's high projection. No arm is changed after the slice.
Reported: the three (or four) points as measured; no winner rule, this is a measurement.

**Machine.** At most 6 heavy processes. Render with up to 4 workers (about 10 to 40 minutes);
slices alone, in the order `CS`, `J5S`, `PS` (each 5 to 15 minutes, under 1 hour in all, `CS` = 2
processes). Whole-layer runs timed alone, one after another, in the order `CS`, `J5S`, `PS`, then
`CF` if it exists: about 12 to 48 hours serial. Each run detached (`nohup caffeinate -i`), with
`--resume`, a `--dry-run`, a log and an end marker `<arm>.done`; a launcher starts the next only
after the previous marker. The cache lives outside the worktree (lesson of `cpptime`).
**Process limit exception:** none; the C++ build is part 1's.

**Validity checks.** The grid has 98,280 voxels, the same area, every child inside its parent, and
equals C++'s own split (`SetMinResolution`, `ReconstructionStrategies.h:66-80`) on a small sub-grid
within 1e-6 mm. The stack equals the frames. Both programs get the same inputs (part 1's input
check reused). Freeze committed before any output exists (file times checked, as step 21). The
runs read no truth. `score_cpp.py` matching works at the fine side (neighbour centres 0.0027 mm
apart, tolerance 1e-3 mm, up/down also matched). Every evaluation counted.

**Stop rule:** stop and report if the gradient gate fails, if the grid check fails, if a run starts
before `fine1/frozen.json` is committed, if any log or map of `fine7` is read before all end markers
exist, if a projection passes 48 hours (ask the owner), if a run hits its cap or its log stalls for
30 minutes, or if a timed run's quiet check fails (relabel as contended and report; no rerun).

**Phases.** A: tests first (`tests/test_fine1.py`: grid split on hand-made triangles, the cap, the
neighbour-pair measure, guards); grid; sample; gradient gate; render; stack; freeze. B: slices,
projection file committed; ask if over 48 hours. C: whole-layer runs. D: after all markers, score;
plot first (a new panel `fine7_98280` in `pareto_rows.json`: evaluations against wrong voxels with
Wilson bars, wall time in each label; `pareto_build.py`, `pareto_insert.py` on a fresh read of the
live page); tables; page `findings_fine1.md`; history.

**Files.** `scripts/autoresearch/fine1/`, `tests/test_fine1.py`, `benchmarks/autoresearch/fine1/`.
**Library changes:** none. C++: only part 1's builds.

**Known limits.** `J5S` was chosen on 24,570-voxel samples; this layer is a new clean sample for it.
One sample, one layer, clean images, cubic crystal. Nobody reads its logs or maps before the end
markers. Times on one machine; `CS` uses one client, not its usual several.

**Not in this step.** Realistic images; other Q-max; several C++ clients; `J5SD` (step 22) unless
the owner says; tuning any arm to the fine grid; another random seed; `held2`.

**Choices made.**
- "Stock Python" read as `PS` (unmodified grid-search search and stock MC local optimization, the
  run `cpptime` compared with C++). Alternatives: `S8` (the AR5 house baseline, a tuned seed step);
  or no opt-in key at all (no revisit, no reseed: 740 wrong on the pilot). Owner question.
- The gradient cap at 0.45 degree, so the bound holds by construction and is still measured.
- Missing voxels count as wrong; timed runs alone (slower, but quiet walls); seed 7 (next unused).
- Part 1's fix, if any, adds `CF`; `CS` stays the build of the branch as it is.

### AR5 step 23: end the seed refinement's second MC stage when it has converged

Task `ar5w`. Approved: D-20261119-N1. Starts from `feature/eval-reduction` after step 22's key
phase. Reuses step 22's cache (worktree `agent-aba6e1bdb95dff72a`: the `J5D` runs, their
`mc_traces`, `e1_log`, `defer_steps`, `run.json`, and the `J5SD` passes), read only. Phase A ends
at "READY FOR KEY" or at a gate; phase B only after a recorded decision.

**Why.** The owner asks for fewer evaluations and less wall-clock time at equal accuracy, and for
general early stops. The restart study (`.claude/reports/20261116-restart-distance.md`) found that
the stock MC local optimization rarely refines a start that is already close: from starts 0.1,
0.3, 0.5 and 1.0 degree away it ended within 0.05 degree in 0, 5, 12 and 8 runs of 40; the same
evaluations spent in a shrinking true ball did so in 35, 33, 29 and 34. The schedule is the weak
part, not the restart width. But in `J5SD` the stock MC fit from a neighbour's orientation never
runs: the BFS fits are CMA fits with no MC fallback (0 fallback runs on all eight samples). The
same sampler does run in `J5SD`, as the second stage of the seed step's MC refinement. There it
costs 35% to 43% of `J5SD` on the `held` types, mostly after its last gain. So this step changes
that stage only.

**What was read and run.** The restart report; PLAN steps 11, 12, 15, 16 and 22; step 22's
`tables.md` and `summary.json`; `J5D`'s `run.json`, `mc_traces`, `e1_log` and `defer_steps` on all
eight disclosed samples; `refine_from_candidates`, `MCOptimizer.optimize` and
`variance_minimizing_optimize` in the library. Arithmetic only, in scratch; no run. Nothing of
`held2`; nothing of `fine7`.

**What the record says.** The MC refinements are the largest single part of `J5SD` on the `held`
types; the stock MC fit in BFS costs nothing. Reading line: per sample, `J5SD`'s total; the seed
step's share (costing, CMA fits, MC refinements, rescue, polish and grid searches, as in step 22);
the MC refinements alone; the BFS fits (all CMA); the grid searches (their count in brackets);
the screen pass and the repair.

| sample | J5SD total | seed step | MC refinements | BFS fits | grid (count) | screen, repair |
|---|---|---|---|---|---|---|
| held5 | 20,076,940 | 69.3% | 34.9% | 22.8% | 24.0% (14) | 7.9% |
| held4 | 20,685,850 | 67.2% | 43.4% | 24.0% | 13.3% (8) | 8.8% |
| held3 | 18,363,331 | 66.7% | 39.3% | 24.1% | 16.9% (9) | 9.2% |
| held | 17,114,589 | 62.7% | 43.4% | 27.4% | 8.1% (4) | 9.9% |
| mid | 7,383,368 | 45.4% | 27.7% | 43.0% | 9.0% (2) | 11.5% |
| grad | 6,134,231 | 32.3% | 13.3% | 53.1% | 11.0% (2) | 14.6% |
| big | 2,452,596 | 10.7% | 6.7% | 81.2% | 0.0% (0) | 8.1% |
| held6 | 22,690,096 | 70.0% | 36.1% | 22.0% | 24.4% (16) | 8.0% |

Most of a refinement's evaluations come after its last gain. Reading line: per sample, the MC
refinements of `J5D`'s seed steps; their mean evaluations; the tail after the last gain as a share
of `J5SD`, and the part of it in refinements that `DEF` resumed; the refinements that stopped at
0.95, with those whose last gain came after the first stage in brackets; the ceiling (below) at a
cap of 120 and of 300 evaluations for the second stage.

| sample | refinements | mean evals | tail | in resumed | stops at 0.95 | ceiling 120 | ceil. 300 |
|---|---|---|---|---|---|---|---|
| held5 | 5,926 | 1,183 | 22.5% | 12.8% | 835 (450) | 25.2% | 19.9% |
| held4 | 6,695 | 1,341 | 26.8% | 16.2% | 1,027 (568) | 32.8% | 27.0% |
| held3 | 5,797 | 1,245 | 25.3% | 15.6% | 774 (391) | 28.9% | 23.3% |
| held | 5,919 | 1,256 | 27.0% | 15.4% | 944 (506) | 32.1% | 25.9% |
| mid | 1,724 | 1,184 | 18.2% | 9.9% | 361 (165) | 20.0% | 15.8% |
| grad | 1,026 | 793 | 10.0% | 4.5% | 127 (29) | 7.8% | 4.8% |
| big | 172 | 953 | 4.5% | 1.0% | 44 (15) | 4.4% | 3.1% |
| held6 | 6,382 | 1,283 | 23.5% | 14.2% | 950 (479) | 26.8% | 21.8% |

The **ceiling** is what `J5SD` would save if every stored refinement ended at its own count or at
208 + cap + 1 evaluations, whichever is smaller, with the same refinements and no path change. It
is arithmetic, not a prediction: a new second stage changes paths, and a refinement that no longer
reaches 0.95 sends its step to the next candidate, to the rescue or to a grid search (about
345,000 evaluations each). Three more facts. The first stage holds most early gains: the median
last gain of a 0.95 stop is at 200 to 223 evaluations on the `held` types and `mid`. The second
stage reaches far: on the `held` types the MC ends at 0.95 lie a median 0.44 to 0.49 degree from
the CMA end of the same candidate, and 40% to 47% lie over 0.5 degree away. And late finds are
real: on `held` step 15's low stop lost 68 (k = 1000) and 161 (k = 500) later 0.95 ends, and
step 22 found 79 to 115 seeds a `held` sample only in a resumed refinement.

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels; always the voxel's
**own** hit ratio. **Seed step**: index the voxel, cost 20 candidates, then for the best 5 in turn
a CMA fit and, under 0.95, an MC refinement from the same candidate; it stops at the first end at
0.95 or more, else takes the lowest-cost end. **MC refinement**: the library's
`refine_from_candidates` on one candidate. Its **first stage** (`FindOptimal`) runs blocks of 31
samples in a small cube about the current state, step 0.13 degree halved on each gain, at most
200 trial steps, and stops after three failed blocks: at most 208 evaluations. Its **second
stage** (the variance stage) runs blocks of 10 samples from the first stage's best, doubles the
radius after a failed block up to the box, restarts about that start with the raw radius (median
0.32, 0.64 and then 1.11 degrees out), and grows its budget of 200 while a block's cost variance
exceeds 0.02 squared.
**Box**: `max(diameter/3, 0.2 deg) / 2^MinLocalResolution`, 0.329 degree in these runs. **Gain**:
a strict drop of the running best cost. **Tail**: the evaluations after the last gain. **Ball of
radius r**: a rotation about a uniformly random axis by the angle r times the cube root of a
uniform number, so the points fill a ball of true angle r evenly. **Check**, **`PH(500)`**,
**`DEF`**: as in step 22. **Right** / **wrong**: at most / over 1 degree from the truth
(symmetry-reduced). **Fixed** / **broken**: by voxel. **Arm `J5DC`**: `J5D` with rule `CB`;
**`J5SDC`**: `J5DC` with step 21's screen pass and wide repair. **Replay**: each stored seed step of
`J5D` run again alone, in the script layer, from the same voxel and image.

**The rule, frozen now.** Rule `CB` (converging ball). In the seed step's MC refinement the first
stage runs as now. The second stage is replaced. It starts at the first stage's best with radius r
equal to the box. Each proposal is the current best turned by a rotation drawn in a ball of
radius r. A proposal with a strictly lower cost becomes the best and the new centre. After 10
proposals in a row with no gain, r halves. The stage ends when r falls under 0.01 degree, or after
300 proposals. Then the final evaluation, as now. The refinement's own generator goes on from the
first stage. No hit ratio, no neighbour, no truth and no crystal symmetry inside the stage; each
voxel is fitted on its own, so grains with orientation gradients need nothing special.
Unchanged: the 20 costed and 5 refined candidates, CMA fits, 0.80, 0.95, `PH(500)` and `DEF` (inert
under `CB`, see validity), `G1`, `G2b`, the grid search with its own final refinement, BFS, the
screen pass and the wide repair. The constants 10, 0.01 degree and 300 are set now and not tuned:
10 failures and halving come from the restart study's illustration; 0.01 degree is under half the
median error (0.019 to 0.030 degree); 300 is 1.5 times the stock stage's starting budget.

**Samples, and what has touched them.** `held5`, `held4`, `held3`, `held`, `mid`, `grad`, `big`
and `held6`: disclosed; the rule was chosen looking at their records. `held2`: not used, not
opened. `fine7`: not read. Clean test: a new sample `held7`, built by `held_make_sample.py` exactly
as `held` (250 Voronoi grains, 24,570 voxels, gradients 0 to 4 degrees) with generator seed 7, as
earlier steps built `held3` to `held6`: a new sample, not an extra random seed. Nobody reads its
logs or maps before its runs and passes end. Clean images, Q-max 8.

**Measured.** Seed-step evaluations, split into costing, CMA fits, first and second stage;
second-stage proposals, halvings and ends (floor or cap); candidates tried per step; steps with an
end at 0.95; rescues (`G1`, `G2b`) and grid searches; evaluations by stage; wrong voxels (Wilson
95%); fixed and broken against `J5SD` and `S8`; median error over all 24,570 voxels; wall seconds
of the `mid` runs under the timing preflight (reported, not gated).

**Predictions, written before any run.**
1. With `CB` the first stage of every refinement equals `J5D`'s, evaluation by evaluation (`mid`);
   with `CB` off the replay reproduces `J5D`'s seed steps on `mid` (ends, evaluations).
2. Replay, eight disclosed samples: the seed step's index stage 40% to 60% fewer evaluations than
   `J5D`'s on the `held` types and `mid`, 20% to 45% on `grad` and `big`; 250 to 450 evaluations
   per refinement (`J5D`: 793 to 1,341); the second stage ends at the floor in at least 90% of
   refinements; steps with an end at 0.95 within 3% of `J5D`'s count.
3. Live `J5SDC` against `J5SD`: 15% to 28% fewer evaluations on the `held` types, 10% to 20% on
   `mid`, 3% to 8% on `grad`, 1% to 4% on `big`; broken minus fixed -3 to +4; grid searches within
   2 of `J5D`'s; median error within 0.0005 degree.
4. `held7`: 15% to 28% fewer evaluations than `J5SD`; broken minus fixed -3 to +4; median error
   within 0.0005 degree; the house rule against `S8` passes.

**Decision rules, fixed now.** Gate, on the disclosed samples, before `held7` is built: (g0)
prediction 1, exactly; (g1) the replay on all eight: the index stage at least 30% fewer
evaluations than `J5D`'s on the `held` types and `mid`, and steps that lose their 0.95 end (`J5D`
had one, `CB` has none) at most 3% of the steps; (g2) live `J5SDC` on `mid` and `held`: at least
10% fewer evaluations than `J5SD`, broken minus fixed at most 4, median error at most 0.001 degree
higher, grid searches at most 2 more than `J5D`'s. Gate fails: stop, report, no tuning (no change
of 10, 0.01 degree, 300 or the start radius), no `held7`. Then `J5SDC` passes on a sample if,
against `J5SD`: (1) at least 10% fewer evaluations on the `held` types and `mid`, 2% on `grad`, 1%
on `big`; (2) broken minus fixed at most 4; (3) median error at most 0.001 degree higher; (4) the
house rule against `S8`: at least 20% fewer evaluations, broken minus fixed at most 5, median error
at most 0.005 degree higher. The step passes if `J5SDC` passes on `held7` and on the five disclosed
samples run live (`mid`, `held`, `held6`, `grad`, `big`). Fails on `held7`: no key value; study
merge. Rule (2) allows 4, not step 22's 3: every refinement's second stage changes, so more seed
ends move and BFS routes shift more; 4 stays under the house rule's 5. Set now, before any run.

**Validity checks.** Each refinement's generator is keyed as in `J5D` (voxel, kind, call count);
the second stage continues it. The first stage is checked identical to `J5D`'s on `mid`. Every
evaluation is counted, the final one and any check included. `PH(500)` and `DEF` are inert under
`CB`: the first stage stops after three failed blocks (at most about 96 evaluations with no gain)
and the second stage has at most 300 proposals, so no refinement has 500 evaluations without a
gain; the runs must show zero checks and zero pauses, else stop. A unit test checks the ball (no
angle over r, median angle 0.794 r, axes even over the sphere). The arm reads no truth. Steps whose
taken end moves over 0.05 degree from `J5D`'s are listed with both ends. Guards on `held2`, `held7`
and `fine7`.

**Machine.** Runs go one at a time (owner rule of 16 Nov). Replay on eight samples (about 24
million). Gate: `J5DC` live on `mid` and `held`, each then the screen pass and the wide repair
(about 18 million). Then `J5SDC` live on `held6`, `grad`, `big` (about 25 million). Then `held7`
built (about 15 minutes); `S8` (about 68 million); `J5D` and its passes (about 23 million);
`J5DC` and its passes (about 17 million). In all about 175 million. Step 22's runs did 9 to 20
million an hour, so 9 to 19 process-hours, and the same in wall time, one run at a time. Caps:
480 minutes for `S8`, 300 the rest. `S8` and `J5SD` on `held7` stay as clean baselines.

**Stop rule:** stop and report if a stored count does not reproduce (both tables above), if
prediction 1 fails, if the gate fails, if any check or pause fires under `CB`, if `J5DC` runs on
any sample before `ar5/ball/frozen.json` is committed, if `held7`'s logs or maps are read before
its runs and passes end, if anything of `held2` is opened or listed, if `fine7` is read, or if a
run hits its cap.

**Phase A (script layer).** 1. Link step 22's cache read only; check the counts of both tables;
write only into this step's cache. 2. Tests first, `tests/test_ar5_ball.py`: the ball's draws; r
halves after exactly 10 failures in a row and never otherwise; the stage ends at the floor or at
300; strict-less acceptance moves the centre; the first stage untouched; the generator continued;
every evaluation counted; zero checks and pauses under `CB`; guards on `held2`, `held7`, `fine7`.
3. Commit `frozen.json` (rule, constants, gate, decision rules, hashes). 4. g0 on `mid`. 5. Replay
on the eight samples (g1). 6. Gate (g2). 7. `J5SDC` on `held6`, `grad`, `big`. 8. Build `held7`
(pilot first); `S8`, `J5D` with passes, `J5DC` with passes. Plot first (run evaluations against
wrong voxels: `S8`, `J5SD`, `J5SDC`, six panels); tables; predictions. "READY FOR KEY".

**Phase B (key).** One opt-in key `SeedIndexMCBall` (0 or 1, default 0, off), valid only with
`SeedByIndexing` and `SeedIndexStopHitRatio` above 0: the seed step's MC refinement gets the
converging second stage; the grid search's own refinement and the BFS fits are untouched. Parsed in
`config_file.py`, checked in `check_bfs_keys`. Unit tests; key-off suite with no golden change; key
on, a library run on `mid` equal to the script arm `J5DCl` (step 22's `J5Dl` with `CB`: the
library's draws, no `G2b`), the same map and evaluation count. Page `findings_ar5_ball.md` (plot
first), history, pareto rows, README.

**Files:** `scripts/autoresearch/ar5/ball_*.py`, `tests/test_ar5_ball.py`,
`benchmarks/autoresearch/ar5/ball/` (under 3 MB tracked; the cache untracked).

**Library changes:** none in phase A (the script layer replaces the stage through a subclass, as
steps 11 to 16 did); phase B the key, if the step passes. The C++ code is not changed. With the key
off the library does exactly what it does now.

**What would make us stop the line.** The gate fails; `held7` saves under 10% against `J5SD` or
breaks more than 4 net; or lost 0.95 ends turn into rescues and grid searches that eat the saving.

**Not in this step.** The stock MC fit in BFS (it does not run in `J5SD`; it matters for `S8` and
C++ only); the restart-width fix; the first stage's schedule; the grid search's own refinement; a
pure stop at the check; a smaller `PH` k; CMA settings; learned proposals; the boundary-pass skip;
realistic images; other crystals; a C++ change; another random seed; `held2`; `fine7`.

**Choices made.**
- This idea: replace the second stage, which holds most of the tail, with a stage that ends once
  its step is under 0.01 degree. A general early stop (stop when the search has converged, not
  when a budget that grows with the cost variance runs out), no model, no threshold on the hit
  ratio, symmetry-agnostic. Ceiling 20% to 33% of `J5SD` on the `held` types, 16% to 20% on `mid`.
  Weak side: the stock stage's wide restarts (0.3 to 1.6 degrees out) find the right basin late
  in some refinements, and the ball gives up that reach; lost 0.95 ends move to the next
  candidate, the rescue or a grid search, which only the replay and live runs can size. Small on
  `grad` (5% to 8%) and `big` (3% to 4%). Chosen looking at eight disclosed samples, so `held7` is
  the only clean test. About 9 to 19 hours, one run at a time.
- Runners-up: a ball for the whole refinement, from the candidate (saves up to about 6% more, the
  first stage's 208 evaluations a refinement, but drops the stage where most 0.95 ends are found);
  a fixed budget of 200 for the second stage with no growth (an exact prefix, replayable, but a
  budget, not a schedule: its cube draws still rarely land near the start); the restart-width fix
  alone (no gain in 160 cases at the 5-degree box; at the 0.33-degree box it keeps the variance
  high and makes the stage longer); a pure stop at the check (23% to 31% of `J5S`, about 400
  changed seed steps a sample; a stop, not a schedule); a smaller `PH` k (tuning a rule on the
  disclosed samples); the stock MC fit in BFS (0 evaluations in `J5SD`); learned grid proposals
  (grid 8% to 24% on the `held` types; step 10's cut failed).

### AR5 step 24: the converging ball first, the stock second stage only when no ball end hits 0.95

Task `ar5x`. Approved: D-20261119-N3. Starts from `feature/eval-reduction` at 92005ed2 or later.
Reuses step 23's replay rows (`ar5/ball/cache/replay_on_*`, main checkout) and step 22's cache
(worktree `agent-aba6e1bdb95dff72a`: the `J5D` runs, their `e1_log`, `mc_traces`, `run.json`, and
the `J5SD` passes), read only. Phase A ends at "READY FOR KEY" or at a gate; phase B only after a
recorded decision.

**Why.** The owner asks for fewer evaluations and less wall-clock time at equal accuracy, and for
general early stops. Step 23's ball (rule `CB`) cut the seed step's index stage by 45% to 64%, but
8% to 10% of steps on the `held` types lost their 0.95 end. Every lost end was found by the stock
second stage after its first block, mostly over 0.33 degree away. So the stock stage is needed in
some steps, but not in most. The ball is cheap and finds the same 0.95 end in most steps. This step
runs the ball first on every candidate, and runs the stock second stage only when no ball end
reaches 0.95. It is step 22's idea (do not pay for a slow search while cheaper options remain)
with an earlier pause point: the ball's convergence instead of 500 evaluations with no gain.

**What was read and run.** PLAN steps 22 and 23; step 23's `info.json`, `info_tables.md`,
`frozen.json`, `replay_partial.json`, its report and review
(`.claude/reports/20261119-0119-ar5w.md`, `20261119-1010-review-ar5w.md`); the replay rows of
all eight disclosed samples (22,087 steps); `J5D`'s `e1_log`, `mc_traces` and `run.json`;
`variance_minimizing_optimize`,
`_zero_temp_with_variance`, `refine_from_candidates` and `_index_seed` in the library;
`ball_replay.py`, `ball_rules.py`, `defer_pick.py`. Arithmetic only, in scratch (`size24*.py`); no
run. Nothing of `held2`; nothing of `fine7`.

**What the record says.** Rule `BF` (below) costs exactly what the stored rows say, because a step
either ends in the ball pass (step 23's replay measured it) or falls back to `J5D`'s own stage 2
(step 22's traces measured it). The accounting identity (costing + CMA fits + MC refinements =
`J5D`'s step) holds in all 4,342 fallback steps. Reading line: per sample, `J5D`'s index stage; step
23's `CB` index stage; steps with no ball end at 0.95 (they fall back); `BF`'s index stage; `BF`
against `J5D` on the index stage; the saving as a share of `J5SD`'s total; lost 0.95 ends under
`BF`; seed steps under the 0.80 seed test for `J5D`, `CB` and `BF`.

| sample | J5D index | CB index | fall back | BF index | BF fewer | of J5SD | lost | seed test J5D/CB/BF |
|---|---|---|---|---|---|---|---|---|
| held5 | 8,690,640 | 3,398,238 | 773 | 7,376,267 | 15.1% | 6.55% | 0 | 104 / 540 / 96 |
| held4 | 10,882,299 | 3,884,631 | 899 | 9,076,725 | 16.6% | 8.73% | 0 | 70 / 558 / 68 |
| held3 | 8,840,378 | 3,287,015 | 813 | 7,856,813 | 11.1% | 5.36% | 0 | 90 / 528 / 85 |
| held | 9,129,537 | 3,433,723 | 777 | 7,666,627 | 16.0% | 8.55% | 0 | 66 / 480 / 64 |
| mid | 2,644,672 | 1,099,560 | 173 | 2,099,764 | 20.6% | 7.38% | 0 | 12 / 67 / 12 |
| grad | 1,255,626 | 690,503 | 97 | 1,023,858 | 18.5% | 3.78% | 0 | 11 / 35 / 10 |
| big | 261,653 | 143,064 | 6 | 192,919 | 26.3% | 2.80% | 0 | 0 / 1 / 0 |
| held6 | 10,054,329 | 3,684,745 | 804 | 8,368,681 | 16.8% | 7.43% | 0 | 63 / 504 / 62 |

Two facts from this table. First, step 23's 60% was not the whole cost of `CB`: its seed steps fell
under the 0.80 seed test 480 to 558 times on the `held` types, against 63 to 104 for `J5D` (66 on
`held` matches `run.json`'s `index_fallback_seed_test`). Each of those goes to the rescue and
perhaps a grid search (about 345,000 evaluations). `BF` keeps `J5D`'s count or lowers it. Second,
`BF` pays for the ball in the steps that fall back: in those steps the ball pass costs 2.3% to 3.1%
of `J5SD` on the `held` types (the ball's stage 2, and in lost steps the CMA fits and first stages
of later candidates that `J5D` never reached). The table's saving is net of that.

Where the saving sits. Reading line: per sample, seed steps; first-seed steps; `CB`'s saving on
first seeds alone, as a share of `J5SD`; lost ends on first seeds under `CB`; steps where the ball
finds an end that `J5D` did not (gained); steps whose taken end moves over 0.05 degree under `BF`.

| sample | steps | first seeds | first-seed saving | first lost | gained | moved over 0.05 |
|---|---|---|---|---|---|---|
| held5 | 3,556 | 282 | 1.13% | 7 | 83 | 90 |
| held4 | 4,109 | 338 | 1.37% | 19 | 107 | 124 |
| held3 | 3,391 | 293 | 1.38% | 10 | 72 | 79 |
| held | 3,685 | 297 | 1.52% | 11 | 80 | 92 |
| mid | 1,626 | 162 | 0.74% | 2 | 27 | 51 |
| grad | 1,281 | 348 | 0.52% | 0 | 20 | 19 |
| big | 330 | 39 | 0.07% | 0 | 0 | 5 |
| held6 | 4,109 | 339 | 1.19% | 19 | 113 | 105 |

Reseed steps hold about 95% of the index stage on the `held` types, so a rule for first seeds only
saves 0.7% to 1.5% of `J5SD`. Most moved ends are gained steps: a 0.95 end in place of a lower one.

The fallback steps are the largest pot left, and they cannot be cut by a threshold. Reading line:
per sample, fallback steps that `J5D` resolves (lost under `CB`) and that it does not (never); the
stock stage 2 inside the never steps, and the whole never steps, as shares of `J5SD`; the median
own hit ratio of the best ball end in lost and in never steps; the median hit ratio of `J5D`'s end
in never steps; lost steps whose 0.95 end came from the candidate with the best ball end.

| sample | lost | never | stock stage 2 in never | never steps | ball hit lost / never | J5D never | best cand. |
|---|---|---|---|---|---|---|---|
| held5 | 309 | 464 | 15.7% | 20.3% | 0.761 / 0.736 | 0.886 | 92 |
| held4 | 406 | 493 | 19.1% | 23.8% | 0.771 / 0.759 | 0.896 | 94 |
| held3 | 284 | 529 | 20.6% | 26.3% | 0.765 / 0.744 | 0.893 | 70 |
| held | 352 | 425 | 18.8% | 23.7% | 0.762 / 0.760 | 0.899 | 112 |
| mid | 70 | 103 | 11.1% | 13.9% | 0.847 / 0.817 | 0.888 | 18 |
| grad | 19 | 78 | 5.0% | 7.5% | 0.829 / 0.858 | 0.886 | 15 |
| big | 2 | 4 | 1.9% | 2.3% | 0.849 / 0.894 | 0.945 | 1 |
| held6 | 330 | 474 | 16.4% | 20.5% | 0.765 / 0.759 | 0.901 | 98 |

In never steps the stock stage does real work: it lifts the end from a median 0.74 to 0.76 to a
median 0.89 to 0.90, over the 0.80 seed test. The ball hit ratios of lost and never steps are
nearly the same, so a hit-ratio threshold cannot tell them apart. A fallback on the best ball
candidate only would keep 23% to 32% of the lost ends. Lost ends come from all five candidate
positions (the first holds 27% to 32%); on the `held` types the MC end lies a median 0.62 to 0.67
degree from the same candidate's CMA end (66% to 68% over 0.5 degree), and the refinement's last
gain comes a median 467 to 516 evaluations in (quartiles about 320 and 900).

Candidate (a), stop the stock stage's polishing between restarts, has almost nothing to cut. In
`variance_minimizing_optimize` every failed block (11 evaluations) restarts at once about the first
stage's best; polishing happens only after a gain. A bound of 12 evaluations per stage-2 gain (one
failed polishing block and the repeated start evaluation) is 1.5% to 1.9% of `J5SD` on the `held`
types, 1.1% on `mid`, 0.2% on `grad` and `big`.

**Terms.** **Hit ratio**: the share of simulated peaks that land on lit pixels; always the voxel's
**own** hit ratio. **Seed step**: index the voxel, cost 20 candidates, then for the best 5 in turn
a CMA fit and, under 0.95, an MC refinement from the same candidate; it stops at the first end at
0.95 or more, else takes the lowest-cost end. **Seed test**: the taken end's hit ratio must reach
0.80, else the rescue (`G1`, `G2b`) and then a grid search run. **First seed** / **reseed**: a seed
step on a voxel no fit has reached / on the centre voxel of a cluster that BFS left unfitted.
**MC refinement**, **first stage**, **second stage** (the stock variance stage), **box**, **gain**,
**ball of radius r**: as in step 23. **Ball stage**: step 23's rule `CB` stage (start at the first
stage's best, radius the box, halve after 10 proposals in a row with no gain, end under 0.01
degree or at 300 proposals). **Ball pass**: the seed step with the ball stage as each refinement's
second stage. **Kept state**: a refinement's best orientation, cost, evaluation count and generator
state at the end of its first stage. **Fallback**: `J5D`'s seed step run again from the kept
states. **Lost** / **gained**: `J5D` has an end at 0.95 and the arm has none / the other way round.
**Never step**: neither has one. **Check**, **`PH(500)`**, **`DEF`**: as in step 22. **Right** /
**wrong**: at most / over 1 degree from the truth (symmetry-reduced). **Fixed** / **broken**: by
voxel. **Arm `J5DB`**: `J5D` with rule `BF`; **`J5SDB`**: `J5DB` with step 21's screen pass and wide
repair. **Replay**: a stored seed step of `J5D` run again alone, from the same voxel and image.

**The rule, frozen now.** Rule `BF` (ball first). The seed step runs in two passes.
1. Ball pass. For the best 5 candidates in turn: the CMA fit as now; if its end is under 0.95, the
   MC refinement's first stage as now, then its kept state is stored, then the ball stage runs on
   a copy of the refinement's generator. The pass stops at the first end at 0.95 or more, and the
   step takes it.
2. Fallback, only if no end of the ball pass reaches 0.95. `J5D`'s seed step runs as it always
   does (`PH(500)`, `DEF`, candidate order, pauses and resumes), with two changes: each CMA fit is
   taken from the ball pass, not run again, and each MC refinement starts its stock second stage
   from its kept state, not from the candidate. Its result is `J5D`'s result for that step,
   evaluation by evaluation. If no end reaches 0.95, the lowest-cost end of the fallback is taken,
   as in `J5D`.
The ball stage's constants (10, 0.01 degree, 300, the box as start radius) are step 23's, not
changed and not tuned. No new threshold, no hit ratio inside the stages, no neighbour, no truth, no
crystal symmetry. Each voxel is fitted on its own, so grains with orientation gradients need
nothing special. Unchanged: costing, CMA fits, 0.80, 0.95, `G1`, `G2b`, the grid search with its own
final refinement, BFS, the screen pass and the wide repair. Every evaluation of both passes counts,
including the ball stage in steps that fall back. The rule was chosen looking at the eight
disclosed samples.

**Samples, and what has touched them.** `held5`, `held4`, `held3`, `held`, `mid`, `grad`, `big`
and `held6`: disclosed; the rule was chosen looking at their records. `held2`: not used, not
opened. `fine7`: not read. Clean test: a new sample `held7`, built by `held_make_sample.py` exactly
as `held` (250 Voronoi grains, 24,570 voxels, gradients 0 to 4 degrees) with generator seed 7, as
earlier steps built `held3` to `held6`: a new sample, not an extra random seed. Step 23 stopped
before building it, so nothing of it exists. Nobody reads its logs or maps before its runs and
passes end. Clean images, Q-max 8.

**Measured.** Seed-step evaluations, split into costing, CMA fits, first stage, ball stage and
fallback stage 2; steps that fall back, lost and gained steps (Wilson 95%); taken ends that move
over 0.05 degree; seed-test failures, rescues (`G1`, `G2b`) and grid searches; evaluations by stage;
wrong voxels (Wilson 95%); fixed and broken against `J5SD` and `S8`; median error over all 24,570
voxels; wall seconds of the `mid` runs under the timing preflight (reported, not gated). For a later
study, every fallback step also writes one feature row at no extra evaluation: the five ball ends'
and CMA ends' own hit ratios and costs, their spread in angle, and whether the fallback found 0.95.

**Predictions, written before any run.**
1. Replay on `mid`: the 1,453 steps that end in the ball pass equal step 23's stored rows (stop,
   taken end, evaluations); the 173 fallback steps equal `J5D`'s stored steps (fits, stop, taken
   end); the index stage is 2,099,764 evaluations exactly; 0 lost ends; 0 checks and 0 pauses in
   the ball pass.
2. The tracked sizing script reproduces the first table exactly on all eight samples.
3. Live `J5SDB` against `J5SD`: 4.5% to 10% fewer evaluations on the `held` types, 5% to 9.5% on
   `mid`, 2% to 5% on `grad`, 1.5% to 4% on `big`; broken minus fixed -3 to +3; grid searches within
   2 of `J5SD`'s; seed-test failures within 10 of `J5SD`'s; median error within 0.0005 degree.
4. `held7`: 4.5% to 10% fewer evaluations than `J5SD`; broken minus fixed -3 to +3; median error
   within 0.0005 degree; the house rule against `S8` passes.

**Decision rules, fixed now.** Gate, on the disclosed samples, before `held7` is built: (g0)
prediction 1, exactly; (g1) prediction 2, exactly, with lost 0.95 ends at most 3% of the steps on
every sample (by construction 0; any lost end is a bug and stops the step; the Wilson upper bound
is reported, for example 0 of 3,685 on `held`: 0% to 0.10%) and the index stage at least 10% fewer
evaluations than `J5D`'s on the `held` types and `mid`; (g2) live `J5SDB` on `mid` and `held`: at
least 4% fewer evaluations than `J5SD`, broken minus fixed at most 3, median error at most 0.001
degree higher, grid searches at most 2 more than `J5SD`'s, seed-test failures at most 10 more. Gate
fails: stop, report, no tuning (no change of the ball's constants or of the fallback), no `held7`.
Then `J5SDB` passes on a sample if, against `J5SD`: (1) at least 4% fewer evaluations on the `held`
types and `mid`, 1.5% on `grad`, 1% on `big`; (2) broken minus fixed at most 3; (3) median error at
most 0.001 degree higher; (4) the house rule against `S8`: at least 20% fewer evaluations, broken
minus fixed at most 5, median error at most 0.005 degree higher. The step passes if `J5SDB` passes
on `held7` and on the four disclosed samples run live (`mid`, `held`, `grad`, `big`). Fails on
`held7`: no key value; study merge. Rule (2) allows 3, as step 22: only steps that end in the ball
pass can change, and 79 to 124 taken ends a `held` sample move over 0.05 degree, fewer than step
22's 172 to 251 C steps. The 4% bar sits under the smallest stored saving on a `held` type (5.36%
on `held3`). Set now, before any run.

**Validity checks.** The first stage of every refinement equals `J5D`'s (`mid`). The ball stage
draws from a copy of the generator state at the end of the first stage, so its draws equal step
23's and the fallback's draws equal `J5D`'s; both are checked on `mid`. CMA fits and first stages
are never run twice in a step: the fallback reads them from the ball pass, and the script counts
them once. Every other evaluation is counted, checks and final evaluations included. `PH(500)` and
`DEF` cannot fire in the ball pass (the first stage stops after three failed blocks; the ball has
at most 300 proposals); the runs must show zero checks and pauses there, else stop. The arm reads
no truth. Steps whose taken end moves over 0.05 degree from `J5D`'s are listed with both ends.
Guards on `held2`, `held7` and `fine7`.

**Machine.** Runs go one at a time (owner rule of 16 Nov). Replay on `mid` (about 2.1 million;
step 23's replays ran 2.2 to 3.7 million an hour, so under 1 hour). Gate: `J5DB` live on `mid` and
`held`, each then the screen pass and the wide repair (about 22 million). Then `J5SDB` live on
`grad` and `big` (about 8 million). Then `held7` built (about 15 minutes); `S8` (about 68
million); `J5D` and its passes (about 23 million); `J5DB` and its passes (about 21 million). In all
about 144 million live. Step 22's live runs did 9 to 20 million an hour, so 7 to 16 hours, plus
about 1.25 hours for the replay and the build: about 8 to 17 hours of wall time, one run at a time.
Caps: 480 minutes for `S8`, 300 the rest. `S8` and `J5SD` on `held7` stay as clean baselines for
later steps.

**Stop rule:** stop and report if a stored count does not reproduce (the tables above), if
prediction 1 or 2 fails, if any lost end appears, if the gate fails, if any check or pause fires in
the ball pass, if `J5DB` runs on any sample before `ar5/bf/frozen.json` is committed, if `held7`'s
logs or maps are read before its runs and passes end, if anything of `held2` is opened or listed,
if `fine7` is read, or if a run hits its cap.

**Phase A (script layer).** 1. Link step 23's replay rows and step 22's cache read only; check the
counts of the three tables with `bf_size.py`; write only into this step's cache. 2. Tests first,
`tests/test_ar5_bf.py`: the ball pass stops at the first end at 0.95; the fallback runs only when
the ball pass has none; the fallback equals `defer_rules.run_seed_step` with cached CMA fits and
kept states (no CMA fit or first stage run twice); the generator copy (ball draws equal step 23's,
fallback draws equal `J5D`'s); every evaluation counted, the ball stage in fallback steps included;
zero checks and pauses in the ball pass; rule off gives `J5D`; guards on `held2`, `held7`,
`fine7`. 3. Commit `frozen.json` (rule, constants, gate, decision rules, hashes). 4. g0: replay on
`mid`. 5. g1: `bf_size.py` on the eight samples. 6. g2: `J5SDB` live on `mid` and `held`. 7.
`J5SDB` on `grad` and `big`. 8. Build `held7` (pilot first); `S8`, `J5D` with passes, `J5DB` with
passes. Plot first (run evaluations against wrong voxels: `S8`, `J5SD`, `J5SDB`, five panels);
tables; predictions. "READY FOR KEY".

**Phase B (key).** One opt-in key `SeedIndexMCBallFirst` (0 or 1, default 0, off), valid only with
`SeedByIndexing`, `SeedIndexStopHitRatio` above 0 and `SeedIndexMCDefer` on: the seed step runs the
ball pass, and falls back to the deferred stock refinement from the kept states. The grid search's
own refinement and the BFS fits are untouched. Parsed in `config_file.py`, checked in
`check_bfs_keys`. Unit tests; key-off suite with no golden change; key on, a library run on `mid`
equal to the script arm `J5DBl` (step 22's `J5Dl` with `BF`: the library's draws, no `G2b`), the
same map and evaluation count. Page `findings_ar5_bf.md` (plot first), history, pareto rows,
README.

**Files:** `scripts/autoresearch/ar5/bf_*.py`, `tests/test_ar5_bf.py`,
`benchmarks/autoresearch/ar5/bf/` (under 3 MB tracked; the cache untracked).

**Library changes:** none in phase A (the script layer builds the two passes from step 23's ball
stage and step 22's deferred step, as earlier steps did); phase B the key, if the step passes: the
MC refinement can return its kept state after the first stage and run the ball stage on a copy of
its generator, and the deferred seed step can start a refinement's stock stage from a kept state
(step 22's pause and resume already keep state and generator). The C++ code is not changed. With
the key off the library does exactly what it does now.

**What would make us stop the line.** The gate fails; `held7` saves under 4% against `J5SD` or
breaks more than 3 net; or the moved seed ends shift BFS so that seed-test failures, rescues and
grid searches eat the saving.

**Not in this step.** A stop for the never steps (20% to 26% of `J5SD` on the `held` types, the
largest pot; but the stock stage lifts their ends over the 0.80 seed test, and `CB` alone sent 480
to 558 seeds a `held` sample under it, so a skip needs a predictor and the rescue's cost; this
step's feature rows are for that study); a wider CMA fit or more candidates for reseeds; candidate
(a); first seeds alone; the grid search; `held6` live; the stock MC fit in BFS; the first stage's
schedule; the ball's constants; learned proposals; realistic images; other crystals; a C++ change;
another random seed; `held2`; `fine7`.

**Choices made.**
- This idea: keep step 23's ball where it works and pay for the wide stock stage only where it is
  needed. A general early stop (defer a slow search while a cheap one has not failed), no model, no
  new threshold, no tuned constant, symmetry-agnostic. Sized exactly from stored rows: 5.4% to 8.7%
  of `J5SD` on the `held` types, 7.4% on `mid`, 3.8% on `grad`, 2.8% on `big`; 0 lost ends by
  construction; seed-test failures equal to or below `J5D`'s. Weak side: the saving is modest; the
  ball costs 2.3% to 3.1% of `J5SD` in the steps that fall back; the never steps, the largest pot,
  are left as they are; live BFS routes change through 79 to 124 moved seed ends a `held` sample,
  which only live runs can judge. Chosen looking at eight disclosed samples, so `held7` is the only
  clean test. About 8 to 17 hours, one run at a time.
- Runners-up: (b) the ball on first seeds only (0.7% to 1.5% of `J5SD`; reseeds hold about 95% of
  the index stage); (a) stop the stock stage's polishing between restarts (at most 1.1% to 1.9%; a
  failed block already restarts at once); (c) more candidates or a wider CMA fit for reseeds (lost
  ends come from all five positions and lie a median 0.62 to 0.67 degree from the same candidate's
  CMA end, so more candidates do not target them; a wider CMA fit might, but no stored trace sizes
  it: it needs a replay of about 800 fallback steps a `held` sample); a fallback gated by the ball's
  hit ratio (lost and never steps have nearly the same ball hit ratios); a fallback on the best
  ball candidate only (keeps 23% to 32% of lost ends); a cap on the fallback's stage 2 (tuning a cap
  on the disclosed samples against the 3% limit); (d) the grid search (8% to 24%; step 10's cut
  failed and there is no new idea); a full replay of `held` (about 7.7 million at replay speed,
  2 to 3.5 hours; prediction 1 on `mid` already tests both halves exactly, and `held` runs live).

### AR5 step 25: predict at the stall point whether a seed step will still reach 0.95

Task `ar5y`. Starts from `feature/eval-reduction` after step 24. Reuses step 22's cache (`J5D_<s>`:
`e1_log`, `mc_traces`, `run.json`; `stored_J5SD_<s>`) through the links in `ar5/ball/cache` (main
checkout), read only. Phase 0 (stored rows) ends at gate s1; phase A only if s1 passes.

**Why.** The owner: "One of the biggest challenges we have is our lack of ability to predict if
further optimization will result in a more optimal voxel orientation. If you have a method to help
with this, then we can stop earlier. General early stopping algorithms are desired." Step 24 left
the "never" steps as the largest pot (20% to 26% of `J5SD` on the `held` types). First question, on
stored rows: at a fixed point in a seed step, do features known there tell a step that will still
reach 0.95 from one that never will? A live run follows only if so.

**What was read and run.** PLAN steps 11, 15, 22 to 24; step 24's report and `bf_size.py`; `J5D`'s
`e1_log`, `mc_traces`, `run.json` and `J5SD`'s `repair.json` on the eight disclosed samples.
Arithmetic and small fits in scratch (`size25.py`, `pred25.py`, `gb25.py`, `imp25b.py`); no run, no
truth, nothing of `held2` or `fine7`. Step 24's ball feature rows exist only for `mid` and `held`,
so they cannot train on `held3` to `held6`; not used.

**Terms.** **Stall point D**: the first moment in a `J5D` seed step when all 5 candidates have had
their CMA fit and every MC refinement has reached its first `PH(500)` check (500 evaluations with
no gain) or ended, with no end at 0.95. **End at D**: the lowest-cost end among the CMA ends and
the running bests at D. **Rescued** / **never step**: reaches D, and `J5D` does / does not find a
0.95 end after it. **AUC**: the chance that a random rescued step scores above a random never step
(0.5 is a coin). **LOSO**: train on three of `held3` to `held6`, test on the fourth. Never steps
here are all `J5D` steps with no 0.95 end (step 24's never plus gained steps: 505 = 425 + 80 on
`held`). Others as in steps 22 to 24.

**What the record says.** About one seed step in six reaches D. Most are never steps, and their
work after D is a large pot. Reading line: per sample, seed steps; steps reaching D; rescued and
never steps; evaluations after D in never and in rescued steps, as shares of `J5SD`; seed-test
failures if every step stopped at D (`J5D`'s own count in brackets).

| sample | steps | reach D | rescued | never | after D, never | after D, rescued | seed test if all stop |
|---|---|---|---|---|---|---|---|
| held3 | 3,391 | 686 | 85 | 601 | 14.48% | 1.65% | 100 (90) |
| held4 | 4,109 | 729 | 129 | 600 | 15.92% | 3.06% | 101 (70) |
| held5 | 3,556 | 638 | 91 | 547 | 11.49% | 1.66% | 86 (104) |
| held6 | 4,109 | 688 | 101 | 587 | 13.52% | 1.86% | 91 (63) |
| held | 3,685 | 625 | 120 | 505 | 15.05% | 2.71% | 77 (66) |
| mid | 1,626 | 149 | 19 | 130 | 9.55% | 1.22% | 8 (12) |
| grad | 1,281 | 108 | 10 | 98 | 2.61% | 0.20% | 6 (11) |
| big | 330 | 5 | 1 | 4 | 1.52% | 0.30% | 0 (0) |

In 65% to 68% of never steps on the `held` types the taken end does not change after D. But a plain
stop at D about doubles the seed-test failures (each cost 24,000 to 93,000 evaluations in rescue
and grid search here), so a stop must spare ends under 0.80 and pick out the never steps among the
rest.

The features at D do not separate rescued from never steps. Features (own only; no truth, neighbour
or symmetry): top two hit ratios at D, their gap and cost; best and mean CMA hit ratio; MC lift;
gains and last-gain evaluation; ends over 0.80; reseed; evaluations to D. Reading line: per sample,
AUC of a logistic model, a depth-3 tree, a boosted tree (LOSO; trained on all four for the rest);
then, for the logistic model with the 0.80 guard, its threshold set to stop at most 10% of the
training rescued steps: stored saving, rescued steps stopped.

| sample | AUC logistic | AUC tree | AUC boosted | saving (logistic) | rescued stopped |
|---|---|---|---|---|---|
| held3 | 0.588 | 0.570 | 0.599 | 2.34% | 6 of 85 |
| held4 | 0.574 | 0.575 | 0.564 | 1.59% | 11 of 129 |
| held5 | 0.589 | 0.559 | 0.604 | 1.40% | 7 of 91 |
| held6 | 0.560 | 0.566 | 0.561 | 2.37% | 11 of 101 |
| held | 0.547 | 0.563 | 0.582 | 2.41% | 13 of 120 |
| mid | 0.574 | 0.528 | 0.569 | 2.10% | 4 of 19 |
| grad | 0.349 | 0.594 | 0.540 | 1.92% | 7 of 10 |

No single feature beats 0.56 (the hit ratio at D: 0.50); the median hit ratio at D is 0.86 to 0.90
in both groups. This repeats steps 11 and 15 at a new point: a stall best does not tell whether the
search will still find 0.95. Another label is predictable: "the end improves by over 0.05 after D"
has AUC 0.75 to 0.82, but stopping on it drops 26 to 41 rescued steps a `held` sample.

A stop that predicts nothing has a pot, with step 15's kind of loss (lost ends broke voxels). `R0`:
stop at D whenever the end at D passes 0.80. Reading line: per sample, the stored saving (share of
`J5SD`, arithmetic), lost 0.95 ends (Wilson 95%, share of seed steps), changed ends.

| sample | R0 saving | lost 0.95 ends | changed ends |
|---|---|---|---|
| held3 | 12.42% | 62 (1.43% to 2.34%) | 173 |
| held4 | 15.79% | 99 (1.98% to 2.92%) | 220 |
| held5 | 10.28% | 71 (1.59% to 2.51%) | 164 |
| held6 | 12.98% | 77 (1.50% to 2.34%) | 190 |
| held | 15.19% | 93 (2.06% to 3.08%) | 207 |
| mid | 9.86% | 19 (0.75% to 1.82%) | 48 |

**The rule, frozen now** (live only if s1 passes). Rule `ES`. In `J5D`'s seed step the 5th
refinement also pauses at its first check. At D the step computes the features and the frozen
logistic model's probability p of a rescue (coefficients, scaling and threshold t fitted on `held3`
to `held6` only; t stops at most 10% of their rescued steps). If p < t and the end at D passes
0.80, the step takes the end at D. Else it goes on exactly as `J5D` (the 5th refinement first, then
the resumes). Each voxel on its own; no truth, neighbour or crystal symmetry. All else as `J5SD`.
Arms `J5DE`, `J5SDE`.

**Samples.** `held3` to `held6`: training. `held`, `mid`, `grad`, `big`: disclosed. Clean test:
`held7`, built by `held_make_sample.py` exactly as `held` with generator seed 7 (a new sample, not
an extra random seed; nothing of it exists). Not read before its runs end. Clean images, Q-max 8.

**Decision rules, fixed now.** (s0) the tracked script reproduces the three tables exactly. (s1)
stored rows, before any run: AUC at least 0.70 in every LOSO fold and on `held` and `mid`, and a
stored saving of at least 4% of `J5SD` in every fold and on `held` and `mid`. On the drafting
numbers s1 fails (AUC 0.55 to 0.59, saving 1.4% to 2.4%); if reproduced, the step stops there with
a negative result, and no other feature set, label or threshold is tried in it. Only if s1 passes:
(g0) `mid` replay: rule off equals `J5D`; rule on, each stopped step's end and count equal its
stored state at D. (g2) live `J5SDE` on `mid` and `held` against `J5SD`, pooled: at least 4% fewer
evaluations in the sum; neither sample above its `J5SD` total; on each, broken minus fixed at most
3 and median error at most 0.001 degree higher; in the sum, grid searches at most 2 more and
seed-test failures at most 10 more. Pooling is the route allowance, set now: one grid search (about
330,000) is 4.5% of `mid` but 1.3% of the pooled 24.5 million. The step passes if the five live
samples (`mid`, `held`, `grad`, `big`, `held7`) pooled save at least 4% with grid searches at most
3 more; `held7` alone saves at least 2% with broken minus fixed at most 3 and median error at most
0.001 degree higher; and the house rule against `S8` holds on `held7` (20% fewer, broken minus
fixed at most 5, median at most 0.005 degree higher). Gate fails: stop, no tuning, no `held7`.

**Validity checks.** No leakage: no model sees the sample it is scored on; the others use the one
model fitted on `held3` to `held6`, frozen in `ar5/es/frozen.json` (coefficients, scaling, t,
hashes) before any live run. Live features equal stored ones on `mid` (g0). No truth at run time.
Every evaluation counts, checks included. Guards on `held2`, `held7`, `fine7`.

**Machine.** s0 and s1: under 1 hour, no reconstruction. If s1 passes: `mid` replay (2.6 million);
`J5DE` with passes on `mid`, `held` (22 million), `grad`, `big` (8 million); `held7` built, `S8`
(68 million), `J5D` and `J5DE` with passes (44 million). About 145 million: 8 to 17 hours, one run
at a time. Caps: 480 minutes for `S8`, 300 the rest.

**Stop rule:** s0, s1, g0 or g2 fails; `J5DE` runs before `frozen.json` is committed; `held7` read
early; anything of `held2` opened or listed; `fine7` read; a run hits its cap.

**Known limits.** One route per sample. D is late (500 evaluations with no gain in each
refinement). A rescue may come from the stock stage's random wide restarts, which nothing at D can
foresee. The drafting counts may be off by one per refinement at a check; s0 settles it.

**Files:** `scripts/autoresearch/ar5/es_*.py`, `tests/test_ar5_es.py`,
`benchmarks/autoresearch/ar5/es/` (under 3 MB). No library change.

**Not in this step.** `R0` (a separate plan if the owner wants it); the ball rows (replays of
`held3` to `held6`, about 8 million each); other labels; deep models; `BF` again; realistic images;
other crystals; C++; another random seed; `held2`; `fine7`.

**Choices made.**
- This idea: a learned stop at one fixed point, on features already logged, trained on four samples
  and checked on the others, with a cheap stored-row gate before machine time. General,
  symmetry-agnostic, no new search. Weak side, measured while drafting: the features barely beat a
  coin (AUC 0.55 to 0.60, a boosted tree too), so the expected outcome is a stop at s1 within an
  hour. Worth that hour only as a record that closes this path; a live run is not worth it on these
  numbers (at most 2.4% stored).
- Pooling is the route allowance (one grid search moved step 24's verdict on `mid`); fixed before
  any run, it does not reopen `BF`.
- Runners-up: `R0`; a stop at each refinement's first check (a dropped refinement's rescue is not
  replayable); a sequence model on gain traces (about 2,700 training steps, 406 rescued).

Approved: D-20261119-N5, phase 0 only (s0 and s1); a live phase needs a new decision.

### AR5 step 26: a faster ray indexer with the same candidates (wall time; evaluations unchanged)

Task `ar5z`. Starts from `feature/eval-reduction` after step 25's merge. Reads step 22's cache
through the links in `ar5/ball/cache` (main checkout), read only. Phase 0 (a quiet profile) ends at
gate z0; phase A only if it passes; phase B only after a recorded decision.

**Why.** The owner counts wall-clock time as well as operations. Steps 23 to 25 tried to cut
evaluations in the seed step's second stage and earned no key. The stored runs show that most
seconds are not spent in evaluations. With `FusedCost` on, one evaluation takes 0.07 to 0.11 ms.
Ray indexing counts zero evaluations, yet takes 0.5 to 1.6 s a seed step: 41% to 77% of `J5SD`'s
wall time. In AR1 (`FusedCost` off, 64 seeds) it was under 1%; cheaper evaluations and many reseeds
moved it to the top. A faster indexer that returns the same candidates changes no search, no map
and no evaluation count, so no accuracy trade decides it.

**What was read and run.** PLAN steps 22 to 25; `J5D`'s `run.json`, `e1_log`, `polish.json`,
`records` and `J5SD`'s `repair.json` (eight samples); step 23's replay rows; `ray_indexing.py`;
shmem on AR1, AR3 step 9, steps 17, 18. Scratch `size26.py`; no run; nothing of `held2`, `fine7`.

**Terms.** **Evaluation**: one strict cost of one orientation at one voxel. **Seed step**: index
the voxel, cost 20 candidates, then for the best 5 a CMA fit and, under 0.95, an MC refinement.
**Index stage**: the seed step's evaluations. **Indexing**: the work before them (`index_voxel`:
attach rays, vote, cluster, refine); seconds in `index_s`. **Rescue**: `G1`, `G2b` after a failed
seed test. **Neighbour stage**: BFS neighbour local optimizations and revisits. **Candidates**: the
indexer's output (rotations, explained rays, votes, in order). **Contended**: timed while other
jobs ran (preflight refused); indicative. **Quiet**: the timing preflight (`--require`) passed.

**Where the evaluations go.** Reading line: per sample, `J5SD`'s evaluations (`J5D` stages plus
screen plus repair, equal to step 24's totals), then shares of the index stage, rescue, grid search
(count in brackets), neighbour stage, and screen pass plus repair.

| sample | J5SD evals | index stage | rescue | grid (n) | neighbour | screen + repair |
|---|---|---|---|---|---|---|
| held3 | 18,363,331 | 48.1% | 1.6% | 16.9% (9) | 24.1% | 9.2% |
| held4 | 20,685,850 | 52.6% | 1.3% | 13.3% (8) | 24.0% | 8.8% |
| held5 | 20,076,940 | 43.3% | 2.0% | 24.0% (14) | 22.8% | 7.9% |
| held6 | 22,690,096 | 44.3% | 1.3% | 24.4% (16) | 22.0% | 8.0% |
| held | 17,114,589 | 53.3% | 1.3% | 8.1% (4) | 27.4% | 9.9% |
| mid | 7,383,368 | 35.8% | 0.6% | 9.0% (2) | 43.0% | 11.5% |
| grad | 6,134,231 | 20.5% | 0.8% | 11.0% (2) | 53.1% | 14.6% |
| big | 2,452,596 | 10.7% | 0.0% | 0.0% (0) | 81.2% | 8.1% |

In the index stage, costing is 0.3% to 0.4% of `J5SD` and CMA fits 3.7% to 9.5%; MC is the rest
(about 15% first stage, 85% second, on refinements matched in step 23's replay).

**Where the seconds go (contended).** Reading line: per sample, wall seconds (BFS run, screen pass
and repair; one worker), the share in indexing and its seconds per seed step, the share in all
evaluations and in the neighbour stage, then ms per evaluation: index stage, neighbour stage,
rescue plus grid, repair, screen pass.

| sample | wall s | indexing | s/step | all evals | neighbour | ms/eval idx, nbr, grid, repair, screen |
|---|---|---|---|---|---|---|
| held3 | 5,926 | 71.7% | 1.25 | 26.5% | 8.1% | 0.077, 0.108, 0.068, 0.088, 0.310 |
| held4 | 8,138 | 76.7% | 1.52 | 21.6% | 6.6% | 0.075, 0.107, 0.068, 0.091, 0.300 |
| held5 | 7,536 | 76.4% | 1.62 | 22.2% | 6.5% | 0.076, 0.107, 0.067, 0.089, 0.303 |
| held6 | 7,416 | 72.6% | 1.31 | 25.5% | 7.3% | 0.075, 0.108, 0.068, 0.088, 0.301 |
| held | 5,601 | 72.0% | 1.09 | 25.8% | 8.8% | 0.073, 0.105, 0.066, 0.086, 0.299 |
| mid | 1,665 | 58.1% | 0.59 | 40.1% | 19.9% | 0.073, 0.104, 0.065, 0.086, 0.330 |
| grad | 1,258 | 50.1% | 0.49 | 48.2% | 27.7% | 0.081, 0.107, 0.069, 0.086, 0.327 |
| big | 438 | 40.8% | 0.54 | 58.3% | 47.8% | 0.073, 0.105, n/a, 0.087, 0.351 |

BFS bookkeeping is 0.8% to 2.2%; the one-time ray pairing (17 to 39 s) sits in "all evals".
Indexing is spread, not a tail: time outside evaluations per seed step has a median of 0.50 s on
`mid` and 0.80 s on `held` (10th to 90th percentile 0.23 to 1.10 s, 0.37 to 2.04 s). No voxel is
indexed twice, so a cache cannot help. AR1 split one indexing into voting 12%, clustering 48%,
refining 39%. Clustering counts the neighbours of every vote in a k-d tree, then loops in Python
over up to 150 seeds with 48 ball queries each; refining loops over up to 150 clusters, though 20
are kept.

**The rule, frozen now.** Rule `IX`: a new indexer returns, for every voxel, exactly the old
`index_voxel` output (rotations, explained-ray and vote counts, attached rays, in the same order).
Allowed: vectorizing, reordering, reusing arrays, skipping work that provably cannot change the
output (for example refining a cluster that cannot enter the top 20). Not allowed: any change of
`IndexParams` (tolerances, `n_seeds` 150, `max_pairs` 600,000, the pair subsample generator), of
`top_k` 20, or of anything after indexing. No truth, no threshold, nothing crystal-specific; any
symmetry group works.

**Samples.** Profiling and design: 200 seed voxels of `mid` only. Exactness: all stored seed voxels
of `mid` and `held` (5,311), and 300 per sample (generator 26) from `held3` to `held6`, `grad`,
`big`. Timing: `mid` and `held`. `held7` (generator seed 7 of the held family; nothing of it
exists) stays unused under tier 1: an exact change makes no accuracy claim. `held2` not opened;
`fine7` not read.

**Gates, fixed now.**
- z0 (quiet): the old indexer on the 200 `mid` voxels, scaled to 1,626 steps, is at least 30% of a
  quiet `J5SD` estimate on `mid` (stored evaluations at the session's quiet ms per evaluation).
- e1 (exactness, index only): old and new on every exactness voxel. Tier 1: bit-identical outputs.
  Tier 2, only if tier 1 fails on at most 0.1% of voxels: same order, counts and attached rays,
  rotation entries within 1e-12. Anything else: stop.
- e2 (exactness, live): the `J5SD` script arm with the new indexer on `mid` and `held`. Tier 1:
  maps and counts equal stored `J5SD` exactly (7,383,368 and 17,114,589). Tier 2: evaluations
  within 0.5% on each; at most 0.1% of voxels move over 0.001 degree; against stored `J5SD` judged
  on `mid` and `held` pooled, broken minus fixed at most 3 per sample and median error at most
  0.001 degree higher.
- w1 (seconds, quiet, index only): summed indexing seconds over the 5,311 voxels at least 50% lower
  than the old indexer's, same session, old and new in alternating blocks of 100 voxels.
- w2 (seconds, quiet, whole run): `J5SD` old and new, each run alone, on `mid` and `held`: pooled
  wall seconds (sum of both) at least 25% lower, and at least 15% lower on each sample.
- Tier 2 only, house rule on `held7` (built as `held`): new `J5SD` against `S8`, at least 20% fewer
  wall seconds and evaluations, broken minus fixed at most 5, median at most 0.005 degree higher.

**Predictions, written before any run.** e1 tier 1 holds. w1: indexing 2 to 5 times faster. w2:
`held` 35% to 55% less wall, `mid` 25% to 45%. Evaluation counts unchanged everywhere.

**Costs.** Phase 0: 0.2 hour. e1 and w1: 7,111 voxels, old and new, no evaluations, 2 to 3 hours.
e2 and w2 (the new runs serve both): four live runs, about 49 million evaluations, about 3 hours.
In all 5 to 6 hours, one run at a time. Tier 2 adds `held7`: about 91 million, 4 to 6 hours.

**Stop rules.** z0, e1 or e2 fails; an `IndexParams` value differs between arms; a timed run
without a passing preflight (rerun it; no contended number enters a gate); the first table does not
reproduce; `held2` opened or listed, `fine7` read, or `held7` built under tier 1.

**Validity checks.** The old indexer stays unchanged as the reference. Both arms read the same
`SmoothCost` and rays. One thread (`OMP`, `MKL`, torch), same power state, preflight JSON kept per
timed run; ms per evaluation by stage reported, to show evaluations did not slow. No truth is read.

**Known limits.** Contended stored seconds (z0 checks first); one machine; `S8` gains as well.

**Phase A.** 1. z0; commit `ix/frozen.json` (rule, gates, voxel lists, hashes). 2. Tests first,
`tests/test_ar5_ix.py`: new equals old on fixtures (under 2 rays, no votes, tied votes and
explained rays, more than 600,000 pairs, cubic and hexagonal symmetry). 3. e1, w1. 4. e2, w2. Plot
first (seconds per sample, old and new, indexing and evaluations), then tables. "READY FOR KEY".

**Phase B (key).** Opt-in `RayIndexFast` (default 0, off), valid with `SeedByIndexing`; parsed,
checked in `check_bfs_keys`; key-off suite with no golden change; key on, a library run on `mid`
equal to `J5Dl`'s map and counts. Page `findings_ar5_ix.md`, history, README.

**Files:** `scripts/autoresearch/ar5/ix_*.py`, `tests/test_ar5_ix.py`, `ar5/ix/` (under 3 MB).
**Library:** phase A a new function beside `index_voxel`, unused unless called; phase B the key.

**Not in this step.** Other indexer settings (an accuracy study); a learned proposer; batching;
the neighbour stage; the seed step's second stage; the grid search; C++; realistic images; another
random seed; `held2`; `fine7`.

**Choices made.**
- This idea: the largest wall-time pot (41% to 77% of `J5SD`, contended), exact, 5 to 6 hours.
  Weak side: contended shares; bit-identical k-d tree counts may be hard to keep (a borderline
  distance can flip), hence tier 2; the evaluation count does not move.
- Runners-up, sized from the same rows: (d) a learned proposer in place of ray indexing (same wall
  pot, but candidate quality drives seed-test failures, each 24,000 to 93,000 evaluations; 22,087
  stored seed steps to train on; many hours); (a) batching several voxels' evaluations in one
  tensor call (evaluations are 22% to 27% of wall on the `held` types, 40% `mid`, 58% `big`; the
  neighbour stage runs at 0.105 ms against 0.067 for batched grid calls, so at most about a third
  of its 6.5% to 48%; BFS order changes, so not exact); (b) a cheaper accept in the neighbour stage
  (22% to 27% of evaluations on the `held` types, 43% `mid`, 81% `big`; AR3 step 9's start test
  raised the median error 0.013 to 0.024 degree; steps 17 and 18's skips failed); (c) a smooth
  orientation field per grain, verified with one evaluation per voxel (neighbour stage and repair,
  30% to 89% of evaluations; new machinery, gradients of several degrees, tens of hours).

Approved: D-20261119-N6, phases 0 and A; phase B needs a new decision.

### AR5 step 27: a smaller ray-pair sample in the indexer, checked by the final cost (not exact)

Task `ar5aa`. Starts from `feature/eval-reduction` after step 26's merge. Reads step 22's cache
through the links in `ar5/ix/cache` (main checkout), read only. Phase A only after a recorded
decision; phase B (a key) needs a second one.

**Why.** Ray indexing counts no evaluations, yet it is 50% to 65% of `J5D`'s wall time even with
step 26's faster indexer (contended: `mid` 627 of 1,247 s, `held` 2,668 of 4,097 s). Clustering is
82% to 84% of the indexer, and it is bound by the number of votes. Votes grow with the ray pairs:
a median of 100,046 pairs on `mid` and 229,206 on `held`, under the cap of 600,000. On 400 stored
seed voxels, a cap of 50,000 pairs cut indexing seconds by 54% (`mid`) and 78% (`held`), and the
orientation the run finally accepted stayed among the top 20 candidates in every voxel where it
was before. The seed step's own strict costs and fits then check every candidate as now.

**What was read and run.** STATE; PLAN steps 22 to 26; step 26's report and review; shmem on step
26, reseeds, AR3 steps 6, 8, 9; `ray_indexing.py` (`index_voxel`, `index_voxel_fast`);
`reconstructor.py` (`_reseed_clusters`, `_bfs_expand`). Stored rows: `J5D` logs, `e1_log`,
`recon.mic`, `run.json` of `mid` and `held` (step 26's live runs, equal to step 22's). Scratch
`size27a.py` (reuse, stored rows only) and `size27b.py` (indexer settings on 200 seed voxels of
`mid` and 200 of `held`, a fixed stride, order of settings alternating per voxel, one thread,
contended: load 2.9 to 4.5, `cpptools` busy). No evaluation, no run; nothing of `held2`, `fine7`.

**Terms.** **Seed step**: index the voxel, cost 20 candidates, refine the best 5 (a CMA fit, then
MC under 0.95). **Reseed**: a seed step on the voxel nearest the centre of a cluster of rejected
voxels. **Ray pairs**: pairs of attached rays that vote; above `max_pairs` a fixed random subsample
(generator 0) is kept. **Votes**: rotations cast by the pairs. **Rank**: position of the first
candidate within 1 degree of the orientation the stored run finally accepted at that voxel.
**Evaluation**, **contended**, **quiet**: as step 26. **Paired timing**: old and new run in one
session, interleaved, so machine load falls on both alike.

**(B) Cheaper settings, measured.** Reading line: per setting, the share of 200 voxels with rank 0,
under 5, under 20, and the mean indexing seconds per voxel (contended), `mid` then `held`. All
settings use step 26's clustering, so `base` is today's output.

| setting | top 1 mid / held | top 5 | top 20 | s per voxel mid / held |
|---|---|---|---|---|
| base (600,000 pairs, 150 seeds) | 77.5% / 82.5% | 99.0% / 98.0% | 100% / 99.5% | 0.380 / 0.688 |
| 150,000 pairs | 77.5 / 82.0 | 98.5 / 98.0 | 100 / 99.0 | 0.368 / 0.379 |
| 50,000 pairs | 76.5 / 82.0 | 99.0 / 99.5 | 100 / 100 | 0.174 / 0.150 |
| 20,000 pairs | 77.5 / 82.5 | 99.5 / 98.0 | 99.5 / 99.0 | 0.092 / 0.099 |
| 25 seeds | 77.5 / 82.5 | 99.5 / 98.5 | 99.5 / 99.0 | 0.325 / 0.614 |
| refine top 40 clusters | 77.5 / 82.5 | 99.0 / 98.5 | 99.5 / 99.0 | 0.352 / 0.648 |
| 50,000 pairs and 25 seeds | 77.0 / 82.5 | 99.0 / 100 | 99.0 / 100 | 0.130 / 0.095 |

Paired by voxel, 50,000 pairs lost no voxel from the top 20 (gained 1) and moved 1 from rank 1 to
rank 5. The nearest candidate's angle to the accepted orientation is unchanged (median 0.368 and
0.389 degree, against 0.376 and 0.395). Fewer seeds or refined clusters save little: the rest of
the indexer is 8% to 10% of `base`. A coarse-to-fine vote was not built (pairs already do it).

**(A) Known orientations first, sized.** Reading line: per sample, the share of seed steps whose
accepted orientation lies within 0.5, 1, 2, 5 degrees of a voxel accepted EARLIER in the run
within 5 side lengths; in brackets, the share that also ended at hit ratio 0.95 or more.

| sample | steps | 0.5 | 1 | 2 | 5 |
|---|---|---|---|---|---|
| mid | 1,626 | 84.6% (78.5%) | 90.2% (83.3%) | 92.0% (85.0%) | 92.3% (85.3%) |
| held | 3,685 | 76.6% (67.7%) | 89.7% (78.7%) | 92.2% (80.8%) | 93.1% (81.4%) |

The 8 nearest accepted voxels hold one within 1 degree in 89.4% and 88.7% of steps. Why the BFS
failed there: in 1,339 of 1,354 (`mid`) and 2,873 of 2,901 (`held`) such steps, a voxel within 2
side lengths and 1 degree of the answer was accepted before the voxel's last rejected fit; the
voxel was then fitted a median 6 and 7 more times, best hit ratio a median 0.838 and 0.791, while
the seed step reached 0.990. The neighbour fit (CMA, 60 evaluations, first step 0.05 degree) stops
short of an answer 0.3 to 1 degree away (AR3 step 8's class A, untested there). Reuse would need the
seed step's longer fit, not the same start. Its share is an upper bound; the CMA result is not
measured without evaluations. Saved per success: indexing plus up to 781 (`mid`) and 1,274
(`held`) evaluations; cost per failure: about 209.

**(C), (D).** (C) Exact headroom: the vote count is about 90% of clustering (25 seeds cut it only
5% to 7%), so an exact gain needs a compiled neighbour count (an edge check adds a query); likely
30% to 50% of indexing. (D) A learned proposer: 22,087 stored seed steps to train on; it must
match top 20 in 99.5% to 100%. After (B), indexing is 0.10 to 0.17 s a step, about 30% of a wall
that is itself 27% (`mid`) to 51% (`held`) shorter (projected from the table, contended).

**The rule, frozen now.** Rule `PX`: the indexer runs with `max_pairs` 50,000 in place of 600,000.
Nothing else changes: the pair generator (seed 0), tolerances, `n_seeds` 150, `top_k` 20, the
seed step, rescue, BFS, screen and repair. Both arms use the same indexer function (`index_voxel`,
or `index_voxel_fast` in both if step 26's key is approved first). No truth, no threshold read
from data after this freeze, nothing crystal-specific. Chosen on 400 voxels of `mid` and `held`:
in-sample, so `held7` is the clean test.

**Samples.** Disclosed: `held3` to `held6`, `held`, `mid`, `grad`, `big`. Candidate quality: 300
stored seed voxels per sample (generator 27) of all eight. Gates: `mid` and `held`. Clean test:
`held7` (generator seed 7 of the held family, built as `held`; nothing of it exists). `held2` not
opened or listed; `fine7` not read.

**Gates, fixed now.** Counts first; they do not depend on machine load.
- q1 (candidates, no evaluations): over the 2,400 voxels, the final orientation leaves the top 20
  in at most 0.5% of voxels where it was in it under `base` (gained ones reported apart).
- c1 (indexer work): summed ray pairs and votes of the new runs at least 60% lower than those of
  stored `J5SD`'s seed voxels (indexed again, no evaluations), `mid` and `held` pooled.
- c2 (evaluations): new `J5SD` at most 2% more than stored (7,383,368 + 17,114,589), summed.
- c3 (accuracy, against stored `J5SD`, 23 and 45 wrong; `mid` and `held` pooled): broken minus
  fixed at most 3 per sample; median error at most 0.001 degree higher; seed steps that reach
  the grid search at most 2 more per sample.
- t1 (paired, index only): old and new settings on the 5,311 stored seed voxels of `mid` and
  `held`, alternating blocks of 100, one session: summed indexing seconds at least 50% lower.
- t2 (paired, whole run): `J5SD` old and new in the order old `mid`, new `mid`, new `held`, old
  `held`, one session, preflight before each, load sampled every 30 s: pooled wall at least 25%
  lower. t1 and t2 are judged with any load, and reported with the preflight verdicts; a quiet
  repeat of t2 is added if a quiet window comes, but is not required.
- h1, house rule on `held7` (a non-exact change): new `J5SD` against `S8`, at least 20% fewer
  evaluations, broken minus fixed at most 5, median error at most 0.005 degree higher.

**Predictions, before any run.** q1: 0 to 6 lost. c1: pairs 70% to 80% lower. c2: within 1%.
c3: -3 to +3. t1: 55% to 80% lower. t2: `mid` 25% to 40%, `held` 40% to 60% less. h1 passes.

**Costs.** q1: 2,400 voxels, two settings, no evaluations, about 0.7 hour. t1: about 1.5 hours.
c1 to c3 and t2: four live runs (the two old runs also check that stored `J5SD` reproduces), about
4 hours. h1: build `held7`, run `S8` and new `J5SD` on it, about 4 to 5 hours. In all 10 to 11
hours, one run at a time.

**Stop rules.** q1 fails (no live run); an old run does not reproduce stored `J5SD`; any other
`IndexParams` value differs; c2 or c3 fails (no `held7`); `held2` opened or listed, `fine7` read.
A gate is never loosened after a timed run; every refused preflight is kept.

**Validity checks.** `n_pairs`, `n_votes` logged per seed step (c1 read from the run). One thread.
Same `SmoothCost` and rays in both arms; ms per evaluation by stage reported. No truth read before
scoring (as step 24).

**Known limits.** In-sample choice (two samples); route noise; a voxel under 50,000 pairs is
unchanged; contended t1, t2 give a ratio, not quiet seconds; copper only (the rule is generic).

**Phase A.** 1. Commit `px/frozen.json` (rule, gates, voxel lists, hashes). 2. Tests first,
`tests/test_ar5_px.py`: `max_pairs` reaches `cast_votes`; under 50,000 pairs the output equals
`base`; hexagonal fixture runs. 3. q1, t1. 4. c1 to c3, t2. 5. h1. Plot first (wall per sample,
old and new, indexing and evaluations stacked; beside it broken and fixed), then the tables.

**Phase B (key).** Opt-in `RayIndexMaxPairs N` (default 600,000), checked in `check_bfs_keys`;
key-off suite, no golden change. Page `findings_ar5_px.md`, history, README.
**Files:** `scripts/autoresearch/ar5/px_*.py`, `tests/test_ar5_px.py`, `ar5/px/` (under 3 MB).
Library: phase A none (the setting passes through `IndexParams`); phase B the key.

**Not in this step.** Reuse (A); an exact count (C); a learned proposer (D); `n_seeds`, `top_k`;
the seed step's fits; the neighbour stage; batching; C++; realistic images; extra seeds.

**Choices made.**
- This idea: the largest wall saving per machine hour. Indexing seconds 54% to 78% lower on stored
  voxels with no top-20 loss, one parameter, no new search code; about 10 hours with `held7`.
  Weak side: not exact (route noise), in-sample choice, contended micro-benchmark.
- 50,000, not 20,000 pairs: 20,000 is 1.5 to 1.9 times faster again but lost 1 of 200 top-20
  voxels on each sample; one lost candidate can send a seed step to the grid search (24,000 to
  93,000 evaluations), which c2 would feel. 25 seeds on top saves little more on `held`.
- Runner-up (A), reuse known orientations (8 nearest accepted voxels, one 200-evaluation CMA fit,
  index only under 0.95): at most 79% to 83% of seed steps skip indexing, and evaluations may
  fall; but its success rate needs evaluations to measure, and new BFS code. After (B) it can
  still remove most of the remaining 30%. The natural next step.
- (C) 30% to 50% of indexing at best, compiled code. (D) the same pot as (A), hours of training.

Approved: D-20261120-N1, phase A, with `index_voxel` in both arms; phase B needs a new decision.

### AR5 step 28: a small pair sample first, the full sample again when the seed test fails

Task `ar5ab`. Starts from `feature/eval-reduction` after step 27's merge (its scripts are reused).
Reads step 27's run cache in the main checkout (`ar5/px/cache`), read only. Phase A only after a
recorded decision; a key (phase B) needs a second one.

**Why.** Step 27's rule `PX` (`max_pairs` 50,000) cut indexing seconds 75.8% and pooled wall 44.9%,
but failed two count gates: evaluations +2.19% (limit 2%) and, on `held`, 7 grid searches against 4
(limit 2 more). Each grid search costs about 346,000 evaluations and 23 s. The verdict stands; `PX`
alone is not tried again. The stored rows show the extra grid searches are NOT lost candidates:
the two settings give the fragile fits different starts, and they fail on different voxels. A
second try with the full sample would have kept out of the grid search every voxel `PX` sent there.

**What was read and run.** PLAN step 27, its report, `px_*.py`; `reconstructor.py` (`_index_seed`);
`defer_run.py`, `high500_run.py`, `polish_rules.make_wrapper`, `order_draws.py`. Stored rows of
step 27's old and new runs (`mid`, `held`): `e1_log`, `polish_log`, `run.json`, `px_index_run_*`,
`after.npy`, logs. Scratch `s28_*.py`: stored rows only, and 11 voxels indexed with both settings
(no evaluation, one thread, contended). No run; nothing of `held2`, `fine7`.

**Terms.** **Seed step**: index the voxel, cost 20 candidates, refine the best 5 (CMA, then MC
under 0.95). **Seed test** fails: no candidate, or the step's hit ratio is under 0.80
(`_index_seed`: `hit < min_acceleration_threshold`). Then come the **rescue** (`G1`: CMA from the
step's kept ends, accept at 0.95), the **polish** (`G2b`: local grid and CMA about each rescue end,
accept at 0.95) and the **grid search** (full search). **Old**, **new**: step 27's runs with
600,000 and 50,000 pairs (old = stored `J5SD`). **Shared voxel**: seeded once in both runs.

**What the runs did.** Reading line: per run, seed steps, steps failing the seed test, and where
those ended (rescued by `G1`, by `G2b`, or sent to the grid search).

| run | seed steps | failed seed test | `G1` | `G2b` | grid search |
|---|---|---|---|---|---|
| mid old (= stored) | 1,626 | 12 | 8 | 2 | 2 |
| mid new | 1,664 | 9 | 8 | 1 | 0 |
| held old (= stored) | 3,685 | 66 | 46 | 16 | 4 |
| held new | 3,679 | 63 | 38 | 18 | 7 |

The failure rate barely moves. But which voxels fail changes a lot: of 2,691 shared `held` voxels,
62 fail old and 60 fail new, yet only 20 fail in both (mid: 983 shared, 12 and 9, 4 in both).
Every grid voxel is a reseed. 6 of the 13 grid searches ended under 0.80 (rejected anyway).

**Lost candidates? No.** Reading line: per voxel that reached a grid search in either run, its
outcome old then new, and the rank among the 20 candidates of the orientation finally accepted
(old map / new map) at 600,000 and at 50,000 pairs; -1 would mean not in the top 20.

| voxel | old to new | rank, 600,000 | rank, 50,000 | voxel | old to new | rank, 600,000 | rank, 50,000 |
|---|---|---|---|---|---|---|---|
| mid 7167 | grid, pass | 0 / 0 | 0 / 0 | held 20758 | grid, pass | 1 / 0 | 1 / 0 |
| mid 16979 | grid, pass | 0 / 0 | 0 / 0 | held 15599 | pass, grid | 2 / 2 | 2 / 2 |
| held 7225 | grid, grid | 0 / 0 | 0 / 0 | held 19825 | pass, grid | 0 / 0 | 0 / 0 |
| held 3207 | grid, pass | 0 / 0 | 0 / 0 | held 18884 | `G1`, grid | 0 / 4 | 0 / 4 |
| held 13831 | grid, grid | 0 / 0 | 0 / 0 | held 13284 | `G1`, grid | 0 / 0 | 0 / 0 |
| | | | | held 11428 | pass, grid | 4 / 0 | 0 / 0 |

All 11 accepted orientations stay in the top 5 at both settings, mostly first. The nearest
candidate's angle to them changes by 0 to 0.57 degree, and fits from such starts end differently
(the cost is sharp within 0.2 to 0.5 degree). 5 voxels newly reach the grid search; 4 no longer do.
So step 28 does not restore lost candidates; it is a second try from other starts.

**Replay, sized.** Seed-step fits use a fresh `default_rng(0)`; rescue, polish and grid use draws
keyed by voxel. So a voxel's chain is fixed by its candidates (69 voxels under 50,000 pairs gave
equal evaluations in both runs), and `PXF` at a shared voxel is the new step, then, if it fails,
the stored old chain. Reading line: per sample, over shared voxels, evaluations of seed steps,
rescue, polish and grid; steps sent to the rescue; grid searches; ray pairs fewer than old.

| sample | shared | old | `PX` | `PXF` | `PXF` vs old | to rescue o/PX/PXF | grid o/PX/PXF | pairs fewer |
|---|---|---|---|---|---|---|---|---|
| mid | 983 | 2,483,991 | 1,919,642 | 1,955,036 | -21.3% | 12 / 9 / 4 | 2 / 0 / 0 | 52.4% |
| held | 2,691 | 8,719,207 | 10,010,675 | 8,375,442 | -3.9% | 62 / 60 / 20 | 4 / 7 / 2 | 76.8% |

Fallbacks on shared voxels: 9 (`mid`) and 60 (`held`); their steps cost 47,406 and 272,939
evaluations, their indexings 6.9 and 110.9 s (contended). Without grid searches `PXF` costs 7.5%
and 4.8% more than old: a 50,000-pair step costs 3% to 7% more on average, plus the fallbacks.
The gain rests on 2 to 3 fewer grid searches per run.

**The rule, frozen now.** Rule `PXF`. (1) A seed step (seed or reseed) indexes with `max_pairs`
50,000, as `PX`. (2) If it fails the seed test (hit ratio under 0.80, or no candidate) and that
indexing kept a subsample (the voxel had more than 50,000 ray pairs), the same voxel is indexed
again with 600,000 pairs and the whole seed step runs again on the new candidates, same pick, same
generators, no candidate skipped. (3) If that step passes, it is the result. Else `G1`, `G2b` and
the grid search run as now, from the second step's kept ends. (4) Both steps count as the seed's
evaluations; both indexings count in c1. Nothing else changes (pair generator 0, `n_seeds` 150,
`top_k` 20, tolerances, BFS, screen, repair). No truth, no threshold read from data, nothing
crystal-specific. Script layer only: a wrapper of `BFSReconstruction._index_seed`.

**Samples.** Gates: `mid`, `held` (in-sample); old arm there: step 27's old runs (= stored
`J5SD`), so only the new arm runs. Clean test: `held7` (generator seed 7 of the held family, built
as `held`; nothing of it exists). `held2` not opened or listed; `fine7` not read.

**Gates, fixed now.** Counts first; they do not depend on machine load. No candidate gate (q1).
- c1 (indexer work): ray pairs and votes of the new runs, fallback indexings included, at least
  50% lower than stored `J5SD`'s seed voxels indexed again (step 27's c1 rows), pooled.
- c2 (evaluations): new `J5SD` at most 2% more than stored (7,383,368 + 17,114,589), summed.
- c3 (accuracy, against stored `J5SD`, 23 and 45 wrong; per sample): broken minus fixed at most
  3; median error at most 0.001 degree higher; seed steps reaching the grid search at most 2 more.
- t2 (whole run): pooled wall of the new `J5SD` runs at least 25% lower than step 27's old runs
  (8,219 s), judged with any load; preflight verdicts and sampled load reported.
- h1, house rule on `held7`: new `J5SD` against `S8`, at least 20% fewer evaluations, broken minus
  fixed at most 5, median error at most 0.005 degree higher; wall seconds of both reported.

**Predictions, before any run.** c1: 68% to 76% fewer pairs. c2: -5% to +1% (centre -2.5%).
Fallbacks: `mid` 7 to 14, `held` 50 to 70; steps to the rescue `mid` 2 to 7, `held` 15 to 30.
c3: broken minus fixed -6 to +3; grid searches `mid` 0 to 1, `held` 1 to 4. t2: 35% to 50% lower
(step 27's `PX`: 44.9%). h1 passes.

**Costs.** Tests and wrapper 1 hour; new `J5SD` on `mid` and `held` 1.5 hours; h1 (build `held7`,
`S8`, new `J5SD`) 4 to 5 hours. In all 6.5 to 7.5 hours, one run at a time.

**Stop rules.** A test fails. A fallback at a voxel seeded once in step 27's old run does not equal
the stored old step (evaluations, fits, hit ratio): the replay premise is broken, stop. Any other
`IndexParams` value differs. c1, c2 or c3 fails: no `held7`. `held2` opened or listed, `fine7`
read. A gate is never loosened after a run; every refused preflight is kept.

**Validity checks.** Logged per seed step: `n_pairs`, `n_votes`, seconds of each indexing, a
fallback flag, both steps' evaluations and hit ratios. The stop rule's reproduction check on every
such voxel (about 5 on `mid`, 45 on `held`). One thread; same `SmoothCost` and rays as step 27.

**Known limits.** In-sample: the rule was designed after seeing step 27's `mid` and `held` runs,
and the replay uses them; `held7` is the clean test. The replay covers shared voxels only (60% and
73% of steps) and ignores route changes. Each grid search is 1.4% of pooled evaluations, so c2's
margin is 2 to 3 grid searches. t2 is not paired in one session. Copper only (rule generic).

**Phase A.** 1. Commit `pxf/frozen.json` (rule, gates, predictions, hashes). 2. Tests first,
`tests/test_ar5_pxf.py`: no fallback at or under 50,000 pairs; the fallback step equals a direct
600,000-pair step at the voxel; `G1` gets the second step's kept ends; both steps' evaluations
counted; hexagonal fixture runs. 3. New `J5SD` on `mid`, then `held`; score c1 to c3, t2. 4. h1.
Plot first (wall per sample, old, `PX` and `PXF`, indexing and evaluations stacked; beside it
grid searches and broken and fixed), then the tables.

**Phase B (key).** Opt-in `RayIndexMaxPairs N` and `RayIndexFallbackFull 1` in `check_bfs_keys`;
key-off suite, no golden change; `findings_ar5_pxf.md`, history, README. **Files:** `pxf_*.py`
(reusing `px_launch`, `px_score`, `px_build_held7`), `tests/test_ar5_pxf.py`, `ar5/pxf/` (< 3 MB).

**Not in this step.** `PX` alone; other pair counts; a second 50,000-pair try from another pair
generator (no stored rows to size it); skipping grid searches that end rejected (6 of 13); reuse
(A); `n_seeds`, `top_k`; the seed step's fits; C++; realistic images; extra seeds.

**Choices made.**
- This idea: it keeps `PX`'s indexing saving (pairs 73% fewer counting fallbacks) and, by the
  replay, takes steps to the rescue from 74 to 24 and grid searches from 6 to 2 on shared voxels.
  Weak side: the gain over old is 2 to 3 grid searches per run, a small count; in-sample design.
- No skip of candidates near ones already tried: the second step helps because its starts lie
  0 to 0.6 degree from the first; a skip within 1 degree would drop most and break reproduction.
- Full sample as the fallback, not another 50,000-pair sample: only it reproduces stored steps, so
  it is sized now and checked exactly in the run. Only after a subsampled first indexing (else the
  step would be identical).
- Runner-up (A), reuse known orientations with a longer fit: larger possible saving (79% to 83% of
  seed steps skip indexing), but its success rate needs evaluations and new BFS code.
- Old `J5SD` on `held7` too (about 1 hour) would give a clean count comparison; left as an option.

Approved: D-20261120-N2, phase A, plus plain `J5SD` on `held7` (reported, not a gate).
Phase B needs a new decision.

## Machine-time runs from the review page

Run when the machine is otherwise quiet: item 14 (both keys on E0 and the pilot, with the timing
A/B of item 3), item 5 (the full-sample baseline), item 18 (does C++ fuse the clip), and a measured
C++ time for the pilot region.

### Items 14 and 3: both speed keys together, with a timing A/B

Task `item14`. A measurement step: new scripts and one config only. Review page items 14 and 3,
delegated by D-20261109-N3 (item 14) and D-20261109-N2 (item 3); plan approved by
D-20261111-N25. The C++ time stays a separate step.

**Why.** `FusedCost` (Track S) makes one strict evaluation about 15 times cheaper. `SmoothCostRank`
(Track E) ranks the coarse-level candidates with a smooth cost that `FusedCost` does not cover. The
two keys have never run together. Every stored pilot wall time is contended. Item 5 (the full-sample
run) waits for this step, and N2 says a full-scale smooth arm runs only if item 14 shows a gain
larger than the pilot's 3 voxels.

**Terms.** **E0**: 200 hard test voxels, each solved by a **global search** (a search with no
starting guess). **Pilot**: the 2,000-voxel BFS region of the 500-grain sample, clean images, Q-max
8, BFS seed 0. **Full search**: a global search in a BFS run (seeds and reseeds). **Wrong**: over
1 degree from the truth. **Arms**: `off` (neither key), `F` (`FusedCost 1` only), `S`
(`SmoothCostRank E3+E4,tau=0.25` only), `FS` (both). **Quiet**: the timing preflight passed and the
load stayed low during the run (rule below); anything else is **contended**.

**Arms and what is stored.** Read each row as: what exists, and whether this step reads or runs it.
E0: the Task E harness (`scripts/eval_reduction/smooth_rank_probe.py`), files gitignored in the main
checkout's `benchmarks/eval_reduction/E/cache/rank/`. Pilot: `G/g2/bfs/runs/keyoff_rerun`,
`G/g3/bfs/runs/smooth_only` (`G` = `benchmarks/eval_reduction/G`), `ar2/key/runs/key_off_fused`.

| set | arm | stored files | wrong | evaluations | stored wall | this step |
|---|---|---|---|---|---|---|
| E0 | off | `e0_old_base` | 41 | 48,700 per voxel | 823 s, contended (2 workers) | read |
| E0 | S | `e0_old_e34t25` | 7 | 49,215 + 252 smooth per voxel | none usable | read |
| E0 | F | none in this harness | | | | run |
| E0 | FS | none | | | | run |
| pilot | off | `keyoff_rerun` | 37 | 25,547,935 | 8,997 s, contended | read |
| pilot | S | `smooth_only` | 34 | 24,490,331 + 980,594 smooth | 11,136 s, contended | read |
| pilot | F | `key_off_fused` | 24 | 25,602,383 | 1,633 s, contended | rerun, A |
| pilot | FS | none | | | | run, B |

Stored pairs: E0 off to S is 34 fixed, 0 broken (`E/e2e/tables.md`). Pilot off to S is 14 fixed, 11
broken; off to F is 20 fixed, 7 broken (`G/g3/tables.md`, `ar2/key/tables.md`). This last pair is
the scale of **route noise**: one pixel-edge flip changes the BFS path. S2's E0 (35 of 200 wrong
with and without `FusedCost`, `benchmarks/s2_results/e0.json`) used another harness (`e0_run.py`)
with another |sin eta| setting, so it is not an arm here. E0 `S` has no usable wall (resumed).

**Timing A/B (item 3).** Arm A = pilot `F`, arm B = pilot `FS` (N2: `FusedCost` on in both, one
worker, quiet machine). Order: an untimed warm-up (A config, `--max-voxels 5`, to fill the page
cache), then A1, B1, A2, B2, with a 120 s pause and a preflight before each. Two repeats per arm.
The expected effect (B about twice A) is far larger than a quiet machine's run-to-run spread, so a
second repeat gives the spread and the identity check; a third pair would add about 75 minutes
and decide nothing. Repeats are not random seeds: A2 must equal A1 and B2 must equal B1
(`recon.mic` byte for byte; `records.json` and the `run.json` counts equal once every `wall*` field
is dropped). A1 must equal the stored `key_off_fused` the same way. Only the clock may differ.

**Where.** No file under `icenine_py/icenine/` changes. Reused: `scripts/phase_d/pilot_bfs.py`,
`scripts/autoresearch/ar2/key_off_fused.config` (arm A), `smooth_rank_probe.py` (imported, not
edited), `scripts/speed_parity/s2/s2_search.py` (its section timers, imported),
`scripts/autoresearch/ar2/key_summary.py` (`load_run`), `scripts/common/stats.py` (`wilson`,
`paired_discordant`), `doc_tables.py`, `scripts/dev/timing_preflight.py`, `job_status.sh`.
New, under `scripts/autoresearch/item14/`:
- `both_keys.config`: `key_off_fused.config` plus the line `SmoothCostRank E3+E4,tau=0.25`.
- `e0_arms.py`: the Task E harness with `fused_cost` set and reset beside `smooth_rank` in each
  worker (its own pool initializer, since spawned workers do not see a parent's patch), re-ranking
  off; `run --spec off|E3+E4,tau=0.25 --fused`, and `check --n 20` (rerun stored arms, compare).
- `same_run.py A B`: the identity check above; exit 1 on the first difference, naming it.
- `split.py`: one full search on full-sample voxels 15901 and 17034, arms F and FS, timed by
  section (s2_search's timers plus one on `SmoothRanker.__call__`): screen, smooth, strict, rest.
- `run_item14.sh`: the one command (below). `summary.py`: summary.json, tables.md, figure.
- `tests/test_item14_scripts.py`, written and committed before the scripts.
Every run sets `PYTHONWARNINGS="error:FusedCost is not supported:RuntimeWarning"`, so a silent
fall-back to the torch path stops the run instead of giving a wrong F arm.

**Library changes:** none. Nothing under `icenine_py/icenine/` changes; no golden file moves.

**What is measured.** (1) Wrong voxels per arm with a Wilson interval, and paired fixed/broken by
voxel: E0 F against off, FS against S, FS against F; pilot FS against S, FS against F, F against
off. No p-values. (2) Evaluations: strict and smooth, per voxel on E0, per run and per full search
on the pilot. (3) Wall time of A and B: `wall_total_s` from `run.json` (the BFS; image reads happen
inside it, the stack is memory-mapped) and the whole process under `/usr/bin/time -p`. Median and
range over the two repeats, split into full searches (`wall_seed_s` + `wall_reseed_s`) and local
optimizations. (4) The share each key does not cover: from `split.py`, the part of a full search
spent in the batched screen and in the smooth evaluations (neither runs in the fused C call); from
the A/B, (B - A) / B, the smooth key's share of the pilot. E0 walls are multi-worker and contended.

**Predictions, written before any run.**
1. E0 `F`: 41 wrong (39 to 43), at most 2 broken against `off`; evaluations within 0.1% of 48,700.
2. E0 `FS`: 7 wrong (5 to 9), at most 2 broken against `S`.
3. A1 equals the stored `key_off_fused`: 24 wrong, 25,602,383 evaluations, 69 full searches.
4. Pilot `FS`: 17 to 31 wrong (24 plus or minus 7, the broken count of the route-noise pair); about
   980,000 smooth evaluations (65 full searches at about 15,100 each).
5. A quiet: 1,350 to 1,700 s. B quiet: 2,600 to 3,700 s. B / A between 1.7 and 2.4, because one
   smooth evaluation costs about 1.7 ms, near 100 fused strict evaluations (G1: +25 s per full
   search single-worker; S2: 15.9 us per fused call).
6. Repeats of one arm differ by under 5% in wall time; every repeat is identical in results.
7. Split, both keys: smooth evaluations 40 to 65% of a full search, the screen 20 to 40%, the fused
   strict evaluation under 15%.

**Decision rules, fixed now.**
- R1, the keys work together on E0: `FS` within 3 voxels of `S` (7) and at most 3 broken against it.
  Otherwise the page reports an interaction and does not recommend the pair.
- R2, identity: a repeat that differs from its first run, or A1 from the stored run, stops the step
  (stop rule); no time is reported.
- R3, the full-scale smooth arm (N2 item 5): queued only if wrong(A1) - wrong(B1) > 3. Pilot wrong
  counts are read against route noise (27 flips between `off` and `F`); no other claim on them.
- R4, quiet: each timed run starts only after `timing_preflight.py --require --json` passes. During
  it a sampler logs the 1-minute load every 30 s and the top CPU processes every 5 minutes. A run
  is quiet if no load sample exceeds 2.5 and no other process stays above 20% CPU. If a preflight
  refuses or a run is not quiet, the script stops; no time is reported, the step reports why.
  Phase A checks the limit of 2.5 first: it measures the load the run adds by itself (one
  worker and its threads, on the 30-voxel check). If the idle load plus that is over 2.5, phase
  A sets the limit to the preflight's own limit plus the run's own load, writes it down before
  any timed run, and the page says so. The limit is never changed after a timed run.

**Machine.** At most 6 heavy processes. E0 arms use 4 workers and run first. The warm-up, the timing
arms and `split.py` then run alone: no other study, no build, no test suite (AR5 step 6 waits).
Estimate: E0 checks and arms about 15 minutes; A/B about 9,200 s; split about 5 minutes; total about
2.9 h. All under `caffeinate -i`; logs and markers in `benchmarks/autoresearch/item14/cache`.

**Stop rule:** stop and report if the run passes 6 hours, if a stored arm does not reproduce
(`e0_arms.py check`, or A1 against `key_off_fused`), if a repeat is not identical, if a run falls
back from the fused path, or if the machine is not quiet (R4). Nothing is tuned to match a count.

**Phase A, prepare (implementer, own worktree).** `read_on_fetch.sh --save` first. Symlink the
gitignored data from the main checkout (`Examples/Example2.ManyGrains/ScatteringData_PhaseD`,
`icenine_py/scripts/findoptimal_robustness/cache`); never commit them. Build `_rasterize` with the
fused entry. Tests first, then the scripts. Short checks: `e0_arms.py check --n 20` on `off` and
`S` (identical to stored); E0 `F` and `FS` on 5 voxels; `pilot_bfs.py --max-voxels 30` with
`both_keys.config` (no fall-back, smooth evaluations above 0); `split.py` on one voxel; a dry run of
`run_item14.sh`. Then stop with `READY TO RUN` and the one command:
`nohup caffeinate -i bash scripts/autoresearch/item14/run_item14.sh > /dev/null 2>&1 &`
The script runs, in order: E0 `F` and `FS` (4 workers); the warm-up; A1, B1, A2, B2, each gated
and followed by `same_run.py`; `split.py`; a closing preflight. The main session starts it.

**Phase B, results.** `summary.py` writes `benchmarks/autoresearch/item14/summary.json`,
`tables.md` and `item14_time.png`; A1 and B1 are copied to `runs/` (repeats stay in the cache). Page
`docs/autoresearch/findings_item14.md`: first the plot (wall time per arm, the two repeats as dots,
stacked into full searches and local optimizations; beside it the wrong counts of the four arms on
E0 and on the pilot), then the tables, each with a reading line; predictions against results;
every wall time labelled quiet or contended. Numbers: `sync_doc_tables.py`, `audit_numbers.py`.
A `MIGRATION_HISTORY.md` section; a mailbox line to Ilium (S2's key ran with the smooth key).

**What this step does not test.** Accuracy on realistic images, the full sample or another layer.
`BatchedQuickMC`. The smooth key's cost with `FusedCost` off (N2: not run). The cause of any pilot
difference beyond route noise. Other worker counts or machines. C++.

**C++ time: out, a separate step.** The binary exists (Release, 2026-03-23, no `Src/` commit
since), but its BFS (`LazyBFS`) runs only in the MPI mode `p` (`RunLazyBFS` in
`Src/XDMServer.cpp`); the serial mode `r` has no BFS. That needs `mpirun` with one server and one
client (not run on this machine in the record) and a C++ config without the Python-only keys (the
parser refuses unknown lines, `Src/ConfigFile.cpp`), pointed at the pilot's `region_grid.mic`. Its
BFS has no reseed, so its work differs and needs its own wrong count. Its wall would include
reading the 360 ASCII image files. Plan it next, in the same quiet window if one remains.

**Choices made.**
- E0 in the Task E harness (stored `off` and `S`); not S2's harness, which has no smooth arm.
- Stored pilot `S` not rerun on today's code (about 3 h with `FusedCost` off); A1 checks the code.
- Two repeats per timing arm, not three. One warm-up run against a cold page cache.
- The uncovered share from a separate timed split, not from timers inside the timed runs.
- C++ time left out (above): a new MPI run, not a small addition.

### C++ time on the pilot: one client, beside the Python timing A/B

Task `cpptime`. A measurement step after "Items 14 and 3": new scripts and one generated config.
No file under `Src/` or `icenine_py/icenine/` changes; the stored binary is used as it is.

**Library changes:** none. Nothing under `icenine_py/icenine/` or `Src/` changes; no rebuild.

**Why.** The owner wants an honest Python-against-C++ wall time and work count for one region.
Item 14 times Python on the pilot (arm A, `FusedCost` on, about 1,500 s, one worker). The only
C++ timing on record is ThreeVoxels (3 voxels, March 2026: 2.25 s, 150,251 evaluations, 15 us
each; memory `project_reconstruction_benchmark`). C++ never ran on the 500-grain sample
(`MIGRATION_HISTORY.md` line 4600).

**Terms.** **Pilot**: the 2,000-voxel BFS region of the 500-grain sample
(`benchmarks/phase_d_pilot/region_grid.mic`, 60 grains), clean images, Q-max 8. **Global search**:
a search with no starting guess (C++ `ReconstructVoxel`, Python `reconstruct_voxel`). **Local
optimization**: the MC fit from a neighbour's orientation. **Server** and **client**: the two MPI
processes of C++ mode `p`; the server hands out seed voxels, the client fits. **Wrong**: over
1 degree from the truth (`rand_500grains_1mm_neworient_s0.mic`, `region_index.npy`).
**Start-up**: process start to the first BFS work (config, image reads, setup). **Quiet** and
**contended**: as rule R4 of item 14.

**What is known (checked 2026-11-12).**
- Binary `IceNine` at the root of the main checkout (gitignored, not in worktrees): arm64, built
  2026-03-23 19:56, `CMakeCache.txt` says `Release`, `-O3 -DNDEBUG`. Three sources are newer:
  exactly the files of commit `b9b566b` (the evaluation counter); the binary already holds its
  symbol `g_cost_eval_count` and its `VOXEL_SUMMARY` strings. Single-threaded: `#pragma omp` in
  `ReconstructionSetup.cpp:78` is inert (no OpenMP flag in `CMakeLists.txt`).
- MPI: `/opt/homebrew/bin/mpirun` is MPICH Hydra 5.0.0; the binary links MPICH (`otool -L`).
- Mode `p`: `Driver.cpp:94-96` to `ParallelDriver.h:65-83`: rank 0 is the server, every other rank
  a client. `-n 1` hangs: with no client the server waits forever (`BreadthFirstReconstructor
  .tmpl.cpp:374-408`). So `mpirun -n 2` = one server (fits nothing) + one fitting client: the
  match for Python's single worker. Several clients change the search (private grids) and make
  the order depend on timing.
- Config: the parser exits on an unknown line (`ConfigFile.cpp:704-708`). Of
  `scripts/autoresearch/ar2/key_off_fused.config` (arm A) only three lines are Python-only:
  `BFSRevisitRefit 1`, `BFSReseed 1`, `FusedCost 1` (checked against the keyword list,
  `ConfigFile.cpp:74-180`). `SampleFilename` must point at `region_grid.mic` (`pilot_bfs.py:76`).
  No C++ config for this sample exists (`Examples/Example2.ManyGrains/ConfigFiles/`: simulation).
- Every C++ failure exits with code 0 (`RUNTIME_ASSERT`, `XDM++/libXDM/Error.h:25-31`; config
  errors, `ConfigFile.cpp:707`), and one rank exiting can leave the other blocked in MPI.
- Output: `<OutStructureBasename>.mic.LBFS`, reduced to the fundamental zone
  (`XDMServer.cpp:219-222`, `:132-133`), voxels in grid order, unvisited voxels left out
  (`ReconstructionStrategies.tmpl.cpp:233-244`). Logs `ParallelReconLog_<rank>.log` in the cwd.
- Timers (whole seconds, `time_t`): server `Total Time elapsed`, from after the image reads (the
  constructor reads them, `XDMServer.h:116`; timer `XDMServer.cpp:405-425`); client `Total
  Processing Time`, BFS only (`XDMClient.cpp:64-103`); client `First Voxel Required: N Sec` per
  global search (`BreadthFirstReconstructor.tmpl.cpp:127-132`); one `QueueSize` line per local
  optimization (`:174`). Stderr prints `Before reading experimental data` and `After reading
  experimental data` (`ReconstructionSetup.cpp:124-126`): timestamped per line, they give the
  image-read time without a code change. Only the server reads; the client gets a broadcast.
- Evaluations: `g_cost_eval_count` counts each `CalculateDiffractionOverlap` call
  (`OverlapInfo.tmpl.cpp:160`) but only serial mode `r` prints it. Mode `p` reports no count.
- Determinism: `std::random_shuffle` and default-seeded `std::mt19937` (`ReconstructionStrategies
  .tmpl.cpp:67`, `XDM++/libXDM/Sampling.h:43-49`): one client should repeat exactly.

**The search is not the same work.** Read each row as: what C++ one-client does, then Python arm A.

| part | C++ LazyBFS, one client | Python arm A (`key_off_fused.config`) |
|---|---|---|
| seed order | `random_shuffle`, server grid (`ReconstructionStrategies.tmpl.cpp:67`) | `rng.shuffle`, seed 0 (`reconstructor.py:1840`) |
| seed test | reject if hit < 0.8 and not converged (`BreadthFirstReconstructor.tmpl.cpp:142`) | reject if hit < 0.8 (`reconstructor.py:2224-2225`) |
| neighbour test | ratio to best hit > 0.9 (`:163`) | the same (`reconstructor.py:2596`) |
| REFIT voxel later | never fitted again: client grid keeps marks (`ReconstructionStrategies.tmpl.cpp:322-354`) | revisited (`BFSRevisitRefit 1`) |
| reseed | none; server seeds only NOT_VISITED (`:151-193`) | `BFSReseed 1`: 32 reseeds in the stored run |
| cost | C++ overlap count | fused C call (`FusedCost 1`) |

Stored arm A: 37 seeds + 32 reseeds = 69 global searches, 25,602,383 evaluations, 24 wrong. The
closest Python run to C++ rules (`phase_d_pilot/runs/mc_clean`: revisit on, no reseed, no fused
cost): 34 seeds, 740 wrong, 12.6 million evaluations, 8,163 s contended.

**Stored, run, where.** Python: item 14's A1 and A2 (read, not rerun). C++: nothing stored; run
C0 (untimed warm-up and evaluation count), C1 and C2 (timed). New, under
`icenine_py/scripts/autoresearch/cpptime/`: `make_config.py` (writes
`cpp_pilot.config` from `key_off_fused.config`: drops the three lines, makes every path absolute,
sets `SampleFilename` to `region_grid.mic`, `OutStructureBasename cpp_pilot`, `SecondsBetweenSave
100000` so no partial file is written; the search lines are untouched); `run_cpp_time.sh` (the one
command); `parse_logs.py` (global searches, local optimizations, seconds, read time); `score_cpp.py`
(matches C++ voxels to the region by centroid and up/down within 1e-3 mm, then wrong per
`score_bfs.pair_misorientation_deg`; missing voxels counted apart); `same_cpp_run.py A B`;
`summary.py`; `tests/test_cpptime_scripts.py`, written first. Reused: `timing_preflight.py`,
`job_status.sh`, `common/stats.py` (`wilson`), `doc_tables.py`. Results:
`benchmarks/autoresearch/cpptime/`.

**What is measured.** (1) Wall: whole process under `/usr/bin/time -p` (also user+sys of both
ranks); start-up and image reads from the timestamped stderr; BFS time from the client timer;
split into global searches (sum of `First Voxel Required`) and the rest. Python: `wall_total_s`
(BFS only; memory-mapped `stacks/clean.npy`, `setup_s` about 0.1 s) and its process time.
(2) Work: global searches, local optimizations, voxels fitted, evaluations (C0 only). (3) Unit
cost: seconds per global search, per local optimization, microseconds per evaluation. (4) Wrong
voxels with a Wilson interval; voxels wrong in both, only in C++, only in Python. No p-values.
Not comparable: total evaluations as "the same job" (different rules), and C++ timers below 1 s.

**Rules.** Quiet: R4 of item 14 (preflight `--require --json` before each timed run; load sampled
every 30 s, top processes every 5 minutes; limit 2.5, set in phase A from the run's own load and
never changed after a timed run). The run's own processes (`mpirun`, `hydra_pmi_proxy`, both
`IceNine` ranks) are not "other" processes. Repeats: two timed runs; C1, C2 and C0 must be
identical (`.mic.LBFS` byte for byte, client log equal once time fields are dropped), or the
difference is reported. No extra random seeds. No p-values. Nothing is tuned to match a count.

**Predictions, written before any run.**
1. The generated config parses; the server log says 2,000 voxels in the queue.
2. C++ wrong: 740 to 1,300 (no revisit and no reseed; Python with revisit, no reseed: 740).
3. C++ global searches: 30 to 45 (Python without reseed: 34); local optimizations 1,900 to 2,000.
4. C++ evaluations: 10 to 16 million.
5. C++ BFS time, quiet: 100 to 600 s; 10 to 30 us per evaluation (ThreeVoxels: 15 us).
6. Start-up 10 to 120 s, most of it the 360 ASCII files (135 MB).
7. C1 and C2 identical in results; walls within 5%. The server uses about one core while waiting.
8. Python A (item 14) over C++ BFS time: 2.5 to 15; per evaluation: 2 to 6.

**Machine.** At most 2 busy processes plus the sampler; nothing else (no build, test suite or
study). After item 14's closing preflight, or a later quiet window. Estimate: C0 about 10
minutes, C1 and C2 2 to 12 minutes each with 120 s pauses; total under 45 minutes; under
`caffeinate -i`; logs and markers in `benchmarks/autoresearch/cpptime/cache`.

**Stop rule:** stop and report if the whole run passes 2 hours, one run passes 1 hour, a log does
not grow for 10 minutes (a blocked rank), stderr holds `RUNTIME_ASSERT` or `[ConfigFile] Error`,
the `.mic.LBFS` is missing, or the machine is not quiet (R4). A differing repeat is reported.

**Phase A, prepare (implementer, own worktree).** `read_on_fetch.sh --save` first. Symlink the
gitignored data from the main checkout (`Examples/Example2.ManyGrains/ScatteringData_PhaseD`); the
binary is called by its main-checkout path. Tests first, then scripts. Short checks:
`mpirun -n 2 hostname`; `make_config.py` and a C++ run on a 30-voxel sub-grid of the region
(parses, reads 360 files, writes `.mic.LBFS`, `score_cpp.py` matches every voxel); the run's own
load (sets the R4 limit); `lldb` reads `g_cost_eval_count` at a breakpoint on `MPI_Finalize` in the
client (`mpirun -n 1 IceNine p cfg : -n 1 lldb -b ...`); a few frames of `stacks/clean.npy` equal
the ASCII files; a dry run. Then stop with `READY TO RUN` and the one command, which the main
session starts with the Bash tool in the background (no `&`, no `|| true`, no `/dev/null`):
`caffeinate -i bash scripts/autoresearch/cpptime/run_cpp_time.sh`
The script logs to the cache, turns every check above into its own non-zero exit (IceNine always
exits 0), and runs C0, C1, C2 (each gated, then `same_cpp_run.py`) and a closing preflight.

**Phase B, results.** `summary.py` writes `summary.json`, `tables.md` and `cpp_time.png`. Page
`docs/autoresearch/findings_cpp_time.md`: first the plot (wall per run, Python A1/A2 and C++ C1/C2,
stacked into start-up, global searches and local optimizations; beside it wrong voxels and global
searches per side), then tables, each with a reading line: walls, work, unit costs, wrong voxels.
Predictions against results; every wall labelled quiet or contended. Optional orientation maps
(skill `mic-maps`). `sync_doc_tables.py`, `audit_numbers.py`; a `MIGRATION_HISTORY.md` section.

**What this step does not test.** C++ with several clients (its usual use). Python with C++ rules.
Realistic images, the full sample, another Q-max, layer or machine. Why unit costs differ.

**Choices made.**
- One client (`-n 2`), not several: repeatable and the match for one Python worker.
- Evaluations from the untimed C0 under `lldb`, not a print added at `XDMClient.cpp:103` and a
  rebuild (a `Src/` change). If `lldb` cannot read the counter, evaluations are reported as
  missing; the print is then an open choice for the owner.
- The stored binary, not a rebuild: Release, and the newer sources are only the counter commit.
- Python not rerun: item 14's A is the Python number. A Python arm with C++ rules (no revisit, no
  reseed) is not run; `pilot_bfs.py:77` requires `BFSRevisitRefit 1`, so it needs a new wrapper.
  Unit costs carry the like-for-like comparison instead.
- Image reads timed from stderr timestamps, not a new C++ timer; each side reads its usual images.
