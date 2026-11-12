# AR3 step 11: a held-out sample, and the S8 baseline on it (finished on the second run)

Step 11 of "AR3" in `docs/autoresearch/PLAN.md`, task `ar3k`. No file under `icenine/` changed.
Scripts are `scripts/autoresearch/ar3/held_*.py`. Results are in
`benchmarks/autoresearch/ar3/held/`. Every table below is filled from `tables.md` there.

**Short answer.** The sample `held` was built and its images were rendered. The base check on `big`
reproduced stored run `S8` exactly. The first run of `S8` on `held` reached its 120 minute cap and
did not finish, so the step missed its 120 minute budget. By decision D-20261111-N17 the same run
was made again with a cap of 480 minutes and nothing else changed. That second run finished. With
its boundary pass it needed 59,923,523 evaluations, 75 grid fallbacks, and left 95 wrong voxels
before the pass and 70 after it. Prediction 2 is met; predictions 1 and 3 are not. All wall times
are contended (other jobs were on the machine) and are not timing results. This step made no other
run on `held` and compared no methods. Runs of a later step on `held` are separate: the three runs
of AR5 step 4 (decision D-20261111-N18) started at 16:53 on 2026-11-11, during the second run of
this step. They belong to that step, and none of their results is read or used here.

## Terms

- **Evaluation**: one call of the strict cost function (a raw count). Indexing is not counted.
- **Seed step**: one search for a starting orientation at a voxel (a first seed or a reseed).
- **Grid fallback**: a seed step whose indexed result failed the seed test, so a full grid search
  ran. It is the most expensive event of the run.
- **Boundary pass**: the last step of the pipeline, which retries the orientation of voxels next to
  a different orientation. A voxel counts as a boundary voxel when the pass finds a neighbour in
  the map whose orientation differs from its own by more than 0.5 degree (symmetry-reduced); this
  is read from the map, not from the truth.
- **Wrong voxel**: a voxel more than 1 degree from the true orientation (symmetry-reduced).
- **Held out**: no choice of any method has been made by reading this sample.

## Why this sample, and why a new generator seed

Three samples exist (`big`, the first gradient sample, `mid`) and each has been read for a choice.
A rule picked on all the samples can fit them. The next whole-layer tests need one sample that no
choice has seen. `held` is that sample. It is a different sample for a carry-over test: does a
result found on the other samples hold on a new one? It is not a repeat with a new random seed to
make a comparison look firmer. For that reason it has its own generator seed (1, where `mid` used
0) for the centres, the orientations and the gradient levels.

## Written down before the run

Predictions of the plan, for `S8` on `held` with the pass:

1. Between 20 and 40 million evaluations.
2. Between 30 and 120 wrong voxels after the pass (0.12% to 0.49%).
3. Between 5 and 40 grid fallbacks.

## The sample

`held` is built like `mid` (`held_make_sample.py`; its pure parts are tested in
`tests/test_ar3_held_sample.py`): 250 Voronoi cells on the same 24,570 voxel centres, six gradient
levels (0, 0.2, 0.5, 1, 2 and 4 degrees), and a redraw of any orientation closer than 5 degrees to
a neighbour grain's (none was needed). Images are clean, rendered with one worker and packed
(`held_stack_note.py`); the stack matches the render record. The render script also always writes
a realistic variant: its files sit untracked in `ScatteringData_PhaseD/held/realistic_nodrop/`
and were not used, so no later step should read them by mistake. The render record's field
`sample_seed` (0) is the default of the render script's `--seed` option, which is unused when
`--mic` is given; it is not the generator seed of `held` (1).

Reading of the table: one row, the number of grains, the grain sizes in voxels (smallest, median,
largest), the lit pixels and spots of the images, and the largest gap between the measured spread
of a grain and its level in degrees.

<!-- table:ar3_held_sample -->
| grains | smin | smed | smax | lit | spots | err |
|---|---|---|---|---|---|---|
| 250 | 9 | 89.5 | 271 | 10,897,893 | 179,718 | 0.00005 |
<!-- /table:ar3_held_sample -->

Grains per level (0 to 4 degrees): 39, 40, 44, 37, 46 and 44.

## The base check on `big`

The script layer of AR5 step 2, with both of its changes off, was run on `big` first
(`base_check.json`). It reproduces stored `S8`: 7,201,863 evaluations, the same evaluations in
every stage, the same seed-step log (329 steps), the same map (largest difference 0.0) and 67 wrong
voxels before the pass. The check passed, so the run on `held` went ahead.

## The two runs on `held`

The first run was stopped by its 120 minute cap (7,200 seconds). Its counters at that moment are
kept in `cache/S8_held_cap120/` (untracked) and in `summary.json` under `first_run_capped`:
48,312,402 evaluations before the pass, 61 grid fallbacks, 305 first seeds and 2,932 reseeds. It
left no map. The second run, with the cap at 480 minutes, finished and needed 54,733,887
evaluations before the pass. It ran for 8,001.6 seconds (133.4 minutes) and the boundary pass for
583.1 seconds (9.7 minutes), both contended: cpptools used 101% CPU at the start of the run and of
the pass, and the three AR5 step 4 runs on `held` were on the machine from 16:53. These wall times
are not timing results and are not compared with the first run's 120 minutes. The
second run did 3,567 seed steps (305 first seeds and 3,262 reseeds); 2,801 of them stopped early
at a hit ratio of at least 0.95.

## Result: evaluations

Reading of the table: evaluations by stage, with the boundary pass in its own column, and the
number of grid fallbacks. The first three rows are the stored finished runs; the `held` row is the
second, finished run.

<!-- table:ar3_held_evals -->
| sample | grains | seed | grid | neighbour | revisit | pass_ | total | fallbacks |
|---|---|---|---|---|---|---|---|---|
| big | 32 | 2,856,712 | 0 | 3,491,437 | 853,714 | 1,138,565 | 8,340,428 | 0 |
| first gradient | 497 | 3,008,757 | 5,449,501 | 2,798,001 | 3,060,730 | 3,683,241 | 18,000,230 | 16 |
| mid | 120 | 12,007,352 | 5,002,750 | 3,859,661 | 2,517,272 | 3,056,039 | 26,443,074 | 15 |
| held | 250 | 21,257,244 | 25,915,118 | 3,004,901 | 4,556,624 | 5,189,636 | 59,923,523 | 75 |
<!-- /table:ar3_held_evals -->

## Result: wrong voxels

Reading of the table: wrong voxels before and after the pass, as a count with its Wilson 95%
interval over the 24,570 voxels; then the median error over all voxels in degrees; then the
boundary voxels the pass found and their share of the layer. The intervals count voxels as
independent, but wrong voxels cluster in grains, so they are too narrow for statements across
grains. On `held` the wrong voxels sit in 53 grains before the pass and in 34 after it.

<!-- table:ar3_held_wrong -->
| sample | before | after | median_b | median_a | boundary |
|---|---|---|---|---|---|
| big | 67 (0.21 to 0.35%) | 24 (0.07 to 0.15%) | 0.0320 | 0.0316 | 11054 (45.0%) |
| first gradient | 210 (0.75 to 0.98%) | 63 (0.20 to 0.33%) | 0.0229 | 0.0226 | 23126 (94.1%) |
| mid | 84 (0.28 to 0.42%) | 31 (0.09 to 0.18%) | 0.0302 | 0.0297 | 17859 (72.7%) |
| held | 95 (0.32 to 0.47%) | 70 (0.23 to 0.36%) | 0.0288 | 0.0286 | 21131 (86.0%) |
<!-- /table:ar3_held_wrong -->

On `held` the pass fixed 30 voxels, broke 5 and left 65 wrong both before and after. Reading of the
next table: the same counts by the gradient level of the voxel's true grain. `fixed` counts voxels
wrong before the pass and right after it; `broken` counts voxels right before and wrong after.

<!-- table:ar3_held_levels -->
| level | grains | voxels | wrong_before | wrong_after | fixed | broken |
|---|---|---|---|---|---|---|
| 0.0 | 39 | 4150 | 7 | 0 | 7 | 0 |
| 0.2 | 40 | 3562 | 7 | 0 | 7 | 0 |
| 0.5 | 44 | 4359 | 2 | 0 | 2 | 0 |
| 1.0 | 37 | 3680 | 9 | 2 | 7 | 0 |
| 2.0 | 46 | 4473 | 21 | 15 | 7 | 1 |
| 4.0 | 44 | 4346 | 49 | 53 | 0 | 4 |
<!-- /table:ar3_held_levels -->

At level 4 degrees the pass fixed no wrong voxel and broke 4 (49 wrong before, 53 after). At the
five lower levels together there were 46 wrong before and 17 after. Why is not tested here.

## Predictions against the result

1. 20 to 40 million evaluations: **not met.** With the pass, 59,923,523 (54,733,887 before it).
2. 30 to 120 wrong voxels after the pass: **met.** 70 (0.23 to 0.36%).
3. 5 to 40 grid fallbacks: **not met.** 75.

Why the guess was low is not tested here. It is consistent with the grid fallbacks: the first
gradient sample and `mid` had 16 and 15 of them, `held` had 75, and the grid searches cost
25,915,118 of its evaluations. Without them the total would be 34,008,405, inside prediction 1;
this is arithmetic, not a run. Whether the number of grains (250), the small grains (28 grains
under 40 voxels, the smallest 9 voxels) or something else drives the fallbacks needs a test that
was not made.

## What this decides

- The S8 baseline on `held` now exists, with the map before and after the pass, the seed-step log
  and the per-voxel records in the cache, so later steps can replay seed steps and pair voxels.
- Nothing is claimed about whether any method carries over: no method was compared on `held`.
- The sample stays unread for choices. Its image stack and sample files are untracked.

## How to reproduce

From `icenine_py/`, after `uv sync --extra dev`. The stored `S8` runs that the base check reads are
in two sibling worktrees (the folders named in `scripts/autoresearch/ar5/fit_common.py`); the
stacks of `big` and its truth map were linked in from them.

```
AR3_SAMPLE=held uv run python scripts/autoresearch/ar3/held_make_sample.py
uv run python scripts/phase_d/render_full.py --tag held --workers 1 --variants clean --mic <held truth .mic>
AR3_SAMPLE=held uv run python scripts/autoresearch/ar3/held_stack_note.py --build
scripts/autoresearch/ar3/held_bg.sh big base_S8x_big ar5/fit_run.py --name S8x
scripts/autoresearch/ar3/held_bg.sh big base_check ar3/held_base_check.py
scripts/autoresearch/ar3/held_chain.sh        # first run, cap 120 minutes (chain_held.*)
(cd benchmarks/autoresearch/ar3/held/cache && mv S8_held S8_held_cap120)
scripts/autoresearch/ar3/held_chain.sh 480 chain_held_480   # second run, cap 480 minutes
scripts/autoresearch/ar3/held_bg.sh held stats ar3/held_stats.py
uv run python scripts/autoresearch/ar3/held_pareto_rows.py
```

The cap argument of `held_chain.sh` was added after the runs. The second run itself was started by
the main session with the same commands the script now runs: from `icenine_py/`, with
`RESEED_SAMPLE=held RESEED_CACHE=<cache>` and under `caffeinate -i`,

```
uv run python scripts/autoresearch/ar5/fit_run.py --name S8 --cap-min 480
uv run python scripts/autoresearch/ar5/fit_pass.py --name S8   # only if no S8_held/capped.json
```

after the first run's folder was moved to `S8_held_cap120`. Its commit (e45bb0b) and start time
(2026-11-11 16:27:44) are in `cache/logs/chain_held_480.commit`.

## Deviations from the plan

- The 120 minute budget was missed. The first run reached the cap; the second run finished under a
  cap of 480 minutes, by decision D-20261111-N17.
- The render used one worker (the plan allows up to two), to keep the machine at two heavy
  processes while the base check ran.
- The plan allowed at most two heavy processes for this step. During the second run the three AR5
  step 4 runs on `held` (D-20261111-N18) started at 16:53 on the same machine, so four heavy
  processes shared it. Load changes the wall times, which are labelled contended. The counts and
  maps are expected not to depend on it (the base check on `big` reproduced a stored run's counts
  and map exactly), but this was not checked on `held`.
- The image render record `preflight_render_held_nodrop.json` is not committed (it holds an
  absolute home path).

## Limits and not claimed

One sample, one finished run, one voxel order. Wall times are contended and are not timing
results. No cause is claimed for the high number of grid fallbacks or for the pass at level 4
degrees fixing no wrong voxel and breaking 4.
