# RUN_LOG: Task G1, Task E's other 24 full-sample voxels and the SmoothCost cache (2026-10-28)

Branch `feature/eval-reduction-G1` off `feature/eval-reduction` (cd6a76f). Plan: `eval-reduction.md`, "Task G1" (steps 2b and
2c of Task G). Approval D-20261028-N2. Commands from `icenine_py/`.

## 1. What changed in the repository, and why

* `icenine/smooth_cost.py`: `SmoothCostCache`. A `SmoothCost` builds its lit-pixel trees and observed rays on first use (16.7 s on
  the full sample, section 4); before this change every `reconstruct_voxel` call with `SmoothCostRank` on built a new one. The
  cache keeps the most recent `SmoothCost` with the objects it was built from and returns it while a call names the same objects
  (identity of simulator, detector list, range map, `exp_data`, sample, structures), the same eta limit and `min_sin_eta` and an
  unchanged sample-to-lab transform (its bytes, so an in-place change is seen); otherwise it builds a new one. `make_ranker` takes
  an optional `cache`.
* `icenine/reconstructor.py` (the hook): `AdaptiveVoxelReconstructor.smooth_cost_cache = SmoothCostCache()` in `__init__`, passed
  to `make_ranker` in `_smooth_ranker`; `None` restores the per-call build. With `SmoothCostRank` off nothing is built and the
  code path is the one before (the import line, one attribute and one keyword argument are the whole diff).
* `tests/test_smooth_cost.py`: three tests (cache hits and rebuilds; a warm-cache key-on search equals a fresh build bit for bit,
  with all four arms so the trees, rays and attached rays are cached; key off builds no `SmoothCost`).
* `scripts/eval_reduction/g1_full24.py` (new; step 2b driver and summary), `g1_cache_timing.py` (new; step 2c timing and its
  tables), `smooth_rank_probe.py` (one line: the worker flag `no_rerank` skips measurement 2's re-ranking of the baseline lists,
  the part that made Task E's full-sample baseline take 6-14 minutes per voxel). `.gitignore`: `benchmarks/eval_reduction/G/cache/`.

## 2. Decisions taken before the data

* The 24 voxels are `VOXELS_30[6:]` of `benchmarks/phase_d_seed_diag/ray_overlap/voxels.json` (its first six are Task E's six
  interior voxels). Same probe, seed (`default_rng(1)`), clean full-sample images, Q-max 8 configuration and `min_sin_eta` 0 as
  Task E's full-sample runs.
* A fresh strict-search baseline on the 24 for both rules (Task E has none for them), without measurement 2's re-ranking.
* Wrong = error above 1 deg. Broken and fixed voxels traced with Task E's rule (`smooth_rank_summary.trace` / `fixed_trace`,
  D-20261014-N13).
* Decision default (plan): the smooth ranking key carried into G2 is the arm with more fixed minus broken voxels across both rules
  on the 24; if neither arm fixes more than it breaks, G2 runs with `SmoothCostRank` off. Tie-break, added after the first partial
  summary (one arm then led by one voxel, the other had not run that voxel): the cheaper arm per smooth evaluation (Task E's
  measurement 4), as Task E broke its own tie. It was not needed.

## 3. Runs

| What | Command | Wall | Label |
|---|---|---|---|
| pilot, voxel 3156, lenient, strict and E3+E4 | `g1_full24.py run --rule old --spec off --tag base --limit 1` and `--spec "E3+E4,tau=0.25" --tag e34t25 --limit 1`, under `/usr/bin/time -l` | 117 s and 171 s per voxel; peak RSS 1.1 GB and 2.7 GB | contended (2 processes plus a test run) |
| step 2b, 6 runs x 24 voxels | `benchmarks/eval_reduction/G/g1/run_full24.sh` (one process per rule x arm), then `--reverse` helpers for the two E1+E2+E3+E4 runs and the lenient E3+E4 run as slots freed (at most 6 processes) | 23:51 to 01:09 (78 min) | contended |
| summary | `g1_full24.py summary` | seconds | writes `summary.json`, `tables.md` |
| step 2c, whole-voxel timing | `g1_cache_timing.py --retries 30` (preflight passed on the 4th attempt, no busy process; `timing/preflight.json`) | 31 min | single-worker |
| timing tables | `g1_cache_timing.py --summary` | seconds | `timing/tables.md` |
| full suite | `uv run pytest tests/ --junitxml=<scratch>/junit.xml` | see section 6 | |

The pilot's estimate (6 runs x 24 voxels, about 1.2-1.5 h on 6 processes) was inside the plan's budget of about 3 h.

## 4. Results

Step 2b (`summary.json`, `tables.md`): the strict search is wrong on 1 of 24 voxels under the lenient rule (voxel 1729, 57.3 deg) and
on none under the C++ rule. E3+E4 (tau 0.25) fixes voxel 1729 and breaks none under either rule; E1+E2+E3+E4 fixes none and breaks
none (1729 stays wrong). Ranks below are 0-based, from `summary.json` `basin_ranks_of_wrong_voxels`. The fix is traced to the arm:
at level 1 of the arm's run, E3+E4 puts a truth-basin candidate at rank 27 of 3,694, where the strict cost's order of the same list
puts its first basin candidate at rank 2,964, below the cut at 923. In the strict search itself the basin candidate kept at level 0
(rank 43 of 10,332) gives no basin candidate after the level-1 quick MC, so the strict run's loss happens in the level-1 search from
that candidate, not at a cut; the two runs differ from level 1 on (their level-0 cuts keep different quarters). E1+E2+E3+E4 on 1729
keeps a basin candidate through level 1 (rank 11) and loses it after the level-2 quick MC.

Strict evaluations per voxel rise with the key on: +1.3% (E3+E4) and +5.0% (E1+E2+E3+E4) lenient, +3.6% and +7.7% C++, plus 15,119
and 15,762 (lenient) or 11,404 and 11,919 (C++) smooth evaluations per voxel.

Decision (default rule): **E3+E4, tau 0.25 deg is the smooth ranking key carried into G2** (net fixed minus broken across both rules
+1 against 0). The margin is one voxel of 48 voxel-rule pairs; the evidence for the key on full-sample images is thin, and the six
interior voxels of Task E gave 0 fixed and 0 broken.

How to read the next table: one row per rule and arm on the 24 voxels; "wrong" counts carry a 95% Wilson interval; fixed = wrong
in the strict search and right with the arm, broken = the reverse; evaluations are means per voxel.

<!-- table:evalG1_full24 -->
| rule | arm (SmoothCostRank) | n | wrong, strict | wrong, arm | fixed | broken | both wrong | broken: traced to arm / draw / untraced | fixed: traced to arm / draw | median error, right in both (deg), strict / arm | strict evaluations per voxel, strict search | strict evaluations per voxel, arm | smooth evaluations per voxel |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lenient | E3+E4,tau=0.25 | 24 | 1 (0.7-20.2%) | 0 (0.0-13.8%) | 1 | 0 | 0 | 0 / 0 / 0 | 1 / 0 | 0.0062 / 0.0112 | 330584 | 335047 | 15119 |
| lenient | E1+E2+E3+E4 | 24 | 1 (0.7-20.2%) | 1 (0.7-20.2%) | 0 | 0 | 1 | 0 / 0 / 0 | 0 / 0 | 0.0062 / 0.0077 | 330584 | 347207 | 15762 |
| C++ | E3+E4,tau=0.25 | 24 | 0 (0.0-13.8%) | 0 (0.0-13.8%) | 0 | 0 | 0 | 0 / 0 / 0 | 0 / 0 | 0.0096 / 0.0123 | 232124 | 240463 | 11404 |
| C++ | E1+E2+E3+E4 | 24 | 0 (0.0-13.8%) | 0 (0.0-13.8%) | 0 | 0 | 0 | 0 / 0 / 0 | 0 / 0 | 0.0096 / 0.0082 | 232124 | 249960 | 11919 |
<!-- /table:evalG1_full24 -->

How to read the next table: the final error (deg) of every voxel, strict search and each arm, per rule.

<!-- table:evalG1_full24_voxels -->
| voxel | strict, lenient (deg) | e34t25, lenient (deg) | e1234a, lenient (deg) | strict, C++ (deg) | e34t25, C++ (deg) | e1234a, C++ (deg) |
|---|---|---|---|---|---|---|
| 3156 | 0.006 | 0.001 | 0.004 | 0.004 | 0.005 | 0.001 |
| 14489 | 0.007 | 0.011 | 0.011 | 0.009 | 0.010 | 0.007 |
| 22800 | 0.006 | 0.020 | 0.008 | 0.019 | 0.014 | 0.007 |
| 18535 | 0.014 | 0.021 | 0.025 | 0.012 | 0.010 | 0.027 |
| 12260 | 0.002 | 0.010 | 0.002 | 0.010 | 0.009 | 0.002 |
| 11924 | 0.008 | 0.009 | 0.007 | 0.007 | 0.015 | 0.008 |
| 3285 | 0.010 | 0.035 | 0.025 | 0.040 | 0.012 | 0.013 |
| 13455 | 0.003 | 0.001 | 0.004 | 0.002 | 0.003 | 0.002 |
| 3188 | 0.013 | 0.014 | 0.007 | 0.013 | 0.013 | 0.016 |
| 3633 | 0.018 | 0.015 | 0.016 | 0.029 | 0.022 | 0.019 |
| 704 | 0.005 | 0.013 | 0.004 | 0.004 | 0.013 | 0.003 |
| 13333 | 0.006 | 0.002 | 0.007 | 0.003 | 0.125 | 0.049 |
| 17490 | 0.018 | 0.022 | 0.006 | 0.009 | 0.011 | 0.016 |
| 19572 | 0.024 | 0.020 | 0.025 | 0.055 | 0.040 | 0.029 |
| 14771 | 0.022 | 0.028 | 0.026 | 0.051 | 0.067 | 0.053 |
| 1729 | 57.346 | 0.005 | 57.347 | 0.002 | 0.004 | 0.001 |
| 9861 | 0.008 | 0.007 | 0.007 | 0.004 | 0.006 | 0.008 |
| 23300 | 0.004 | 0.007 | 0.009 | 0.013 | 0.003 | 0.012 |
| 19577 | 0.002 | 0.011 | 0.017 | 0.015 | 0.034 | 0.008 |
| 9864 | 0.006 | 0.008 | 0.003 | 0.034 | 0.004 | 0.001 |
| 11927 | 0.005 | 0.008 | 0.008 | 0.003 | 0.004 | 0.018 |
| 14492 | 0.003 | 0.015 | 0.012 | 0.004 | 0.013 | 0.005 |
| 12263 | 0.018 | 0.011 | 0.016 | 0.009 | 0.014 | 0.009 |
| 3634 | 0.005 | 0.001 | 0.003 | 0.019 | 0.022 | 0.002 |
<!-- /table:evalG1_full24_voxels -->

Step 2c (`timing/tables.md`, single-worker): the build takes 16.7 s (4.6 s lit-pixel trees, 12.1 s observed rays). The cache saves
17.1-17.8 s per voxel. Voxel 17034 with E3+E4: 161.6 s with a fresh build, 144.4 s from the warm cache, against 119.2 s with the key
off (+35.5% before, +21.1% with the cache). With E1+E2+E3+E4: 263.9 s and 246.8 s (+121.3% and +106.9%). The warm-cache searches
equal the fresh-build searches bit for bit (final orientation, cost, strict and smooth evaluation counts) on both voxels and both
arms. Voxel 15901's E3+E4 error (0.009509408798505606 deg) equals Task E's end-to-end result.

How to read the next two tables: whole `reconstruct_voxel` wall time in one process, one torch thread; "cold" builds the
`SmoothCost`, "warm" is served by the cache built on the other voxel or an earlier call, "fresh build" is the behaviour before the
cache; then the build on its own.

<!-- table:evalG1_cache_timing -->
| arm (SmoothCostRank) | timing | key off, 15901 / 17034 (s) | 15901: cold cache / warm cache (s) | 17034: fresh build / warm cache (s) | saved per voxel by the cache (s), A / B | extra time over key off, as before (fresh build), B | extra time over key off, warm cache, A / B | warm equals fresh (B), cold equals warm (A) |
|---|---|---|---|---|---|---|---|---|
| E3+E4,tau=0.25 | single-worker | 112.1 / 119.2 | 157.3 / 139.7 | 161.6 / 144.4 | 17.6 / 17.2 | +35.5% | +24.6% / +21.1% | yes / yes |
| E1+E2+E3+E4 | single-worker | 112.1 / 119.2 | 257.1 / 239.3 | 263.9 / 246.8 | 17.8 / 17.1 | +121.3% | +113.5% / +106.9% | yes / yes |
<!-- /table:evalG1_cache_timing -->

<!-- table:evalG1_cache_build -->
| timing | SmoothCost constructor (s) | lit-pixel trees (s) | observed rays (s) | build total (s) |
|---|---|---|---|---|
| single-worker | 0.002 | 4.6 | 12.1 | 16.7 |
<!-- /table:evalG1_cache_build -->

## 5. What went wrong, and deviations

* The worktree has no copy of the ignored full-sample images; `Examples/Example2.ManyGrains/ScatteringData_PhaseD` was linked
  from the main checkout (a symlink, untracked, not committed).
* The worktree isolation hook refused four Bash commands (a compound version-control command, two heredocs, a `uv run python -c`
  inside a compound command). Each refusal named the form of the command, not the action; the file edits were redone with the
  Edit and Write tools and the commands split into plain ones, as the hook's message asked.
* The first two pilot runs were launched before the cache existed; a run's import time decides which library it used, which does
  not change results (the cache is tested equal to a fresh build) but is noted.
* The plan expected the 24 voxels to include voxels the strict search gets wrong; on these images and settings it gets one wrong
  (lenient) and none (C++), so the test has little power.

## 6. Tests

* `uv run pytest tests/test_smooth_cost.py`: 15 passed (12 before plus 3 new), including the key-off base-commit fixture test.
* Full suite (`uv run pytest tests/ --junitxml=<scratch>/junit.xml`, 9 min 44 s, run alone after the timing): 1388 passed,
  44 skipped, 1 deselected.

## 7. Could be done better

* The 24 voxels carry almost no signal (one wrong voxel in 48 voxel-rule pairs). A set chosen for voxels the strict search gets
  wrong on full-sample images (for example from a strict BFS run) would test the key where it can act.
* The cache holds one image set. A caller that alternates between image sets (the E0 harness attaches new images per voxel) gets a
  build each time, as before; a small LRU would cover it if needed.
