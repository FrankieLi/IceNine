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
