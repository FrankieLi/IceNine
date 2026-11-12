# AR3 step 11: a held-out sample, and the S8 baseline on it did not finish

Step 11 of "AR3" in `docs/autoresearch/PLAN.md`, task `ar3k`. No file under `icenine/` changed.
Scripts are `scripts/autoresearch/ar3/held_*.py`. Results are in
`benchmarks/autoresearch/ar3/held/`. Every table below is filled from `tables.md` there.

**Short answer.** The sample `held` was built and its images were rendered. The base check on `big`
reproduced stored run `S8` exactly. The run of `S8` on `held` then reached its 120 minute cap
before it finished. So there is no map, no count of wrong voxels and no boundary pass for `held`
yet. The counters at the cap already show that `held` costs far more than the plan predicted. The
plan's rule is that a run at its cap is reported as unfinished and is not extended. It was not
extended, and no other run was made on `held`.

## Terms

- **Evaluation**: one call of the strict cost function (a raw count). Indexing is not counted.
- **Seed step**: one search for a starting orientation at a voxel (a first seed or a reseed).
- **Grid fallback**: a seed step whose indexed result failed the seed test, so a full grid search
  ran. It is the most expensive event of the run.
- **Boundary pass**: the last step of the pipeline, which retries the orientation of voxels next to
  a different orientation.
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
(`held_stack_note.py`); the stack matches the render record.

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

## The run on `held`

`S8` on `held` ran for 7,200 seconds (120 minutes) and was stopped by the cap. The counters at that
moment are in `summary.json` (`capped`). Reading of the table: evaluations before the pass, by
stage. The three finished samples are complete runs; the `held` row is a count at the cap, so each
number can only grow.

<!-- table:ar3_held_capped -->
| sample | state | seed | grid | neighbour | revisit | before | fallbacks |
|---|---|---|---|---|---|---|---|
| big | finished | 2,856,712 | 0 | 3,491,437 | 853,714 | 7,201,863 | 0 |
| first gradient | finished | 3,008,757 | 5,449,501 | 2,798,001 | 3,060,730 | 14,316,989 | 16 |
| mid | finished | 12,007,352 | 5,002,750 | 3,859,661 | 2,517,272 | 23,387,035 | 15 |
| held | stopped at 120 minutes | 19,679,538 | 21,082,802 | 3,004,901 | 4,545,161 | 48,312,402 | 61 |
<!-- /table:ar3_held_capped -->

What the counters say about `held` at the cap:

- 305 first seeds and 2,932 reseeds, 3,237 seed steps in all.
- 61 grid fallbacks. They cost 21,082,802 evaluations (43.6% of the 48,312,402 before the pass),
  about 345,620 evaluations each.
- The run held 24,570 records, the number of voxels. The count of finished revisits at the cap is
  not known.
- Time in the reseeds was 6,236 of the 7,200 seconds; the grid searches took 1,428 seconds by
  their own timers (whether that lies inside the 6,236 was not checked). These wall times are contended (other jobs were on the machine) and are not
  timing results.

## Predictions against the result

1. 20 to 40 million evaluations: **not met.** The count before the pass is already 48,312,402 at the
   cap, and the pass is still to come.
2. 30 to 120 wrong voxels after the pass: **not evaluated.** There is no finished map.
3. 5 to 40 grid fallbacks: **not met.** There were 61 at the cap, and the run was not finished.

Why the guess was low is not tested here. It is consistent with the grid fallbacks: the first
gradient sample and `mid` had 16 and 15 of them, `held` had 61. Whether the number of grains
(250), the many small grains (smallest 9 voxels) or something else drives the fallbacks needs a
test that was not made.

## What this decides

- There is no S8 baseline on `held` yet. Later steps that need one must have it first.
- A longer cap would finish the run, but the plan forbids extending it. Whether to rerun with a
  longer cap (the run took more than 120 minutes) is the owner's decision.
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
scripts/autoresearch/ar3/held_chain.sh
scripts/autoresearch/ar3/held_bg.sh held stats ar3/held_stats.py
uv run python scripts/autoresearch/ar3/held_pareto_rows.py
```

## Deviations from the plan

- The run did not finish, so there is no pass, no wrong count, no map and no row pair in
  `pareto_rows.json`; the file lists the run under `excluded` instead.
- The render used one worker (the plan allows up to two), to keep the machine at two heavy
  processes while the base check ran.
- The image render record `preflight_render_held_nodrop.json` is not committed (it holds an
  absolute home path).

## Limits and not claimed

One sample, one run, one voxel order, and the run is unfinished. The counts at the cap are lower
bounds. No timing result. No cause is claimed for the high number of grid fallbacks.
