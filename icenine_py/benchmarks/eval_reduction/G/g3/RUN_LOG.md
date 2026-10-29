# RUN_LOG: Task G3, the total and the consolidated findings (2026-10-29)

Branch `feature/eval-reduction-G3` off `feature/eval-reduction` (8a4e01f). Plan: `eval-reduction.md`, "Task G3" (steps 4
and 5 of Task G). Approval D-20261028-N2; the extra BFS arm is D-20261029-N3. Commands from `icenine_py/`. Clean images
only; `BatchedQuickMC`, `SearchDevice`, `CandidateDedup`, `FusedCost` and `OutwardSearch` off in the new arm. No library
change.

## 1. What changed in the repository, and why

* `benchmarks/eval_reduction/G/g3/smooth_only_mc_clean.config`: the stored pilot's configuration
  (`scripts/phase_d/configs/ReconstructPhaseD_mc_clean.config`) plus `BFSReseed 1` and `SmoothCostRank E3+E4,tau=0.25`.
  Checked: no other non-comment line differs from the stored config; against G2's integrated config only the
  `OutwardSearch 1` line is gone. So against G2's `keyoff_rerun` (the stored config with `--reseed`) the one change is the
  smooth key, and against G2's `integrated` the one change is `OutwardSearch`.
* `scripts/eval_reduction/g3_summary.py` (new): the BFS runs paired by voxel, the trace of every flip and every wrong voxel,
  the N13 trace from probe re-runs, the totals by stage (raw and strict-equivalent), the full-sample estimate, the
  single-lever table read from each task's own tables, the figure and the maps spec.
* `scripts/eval_reduction/g3_trace_probe.py` (new): Task E's probe (`smooth_rank_probe.task`) on the full-sample voxels of
  the full searches behind the arm's flips, key off and on.
* `docs/findings_eval_reduction_2026-10.md` filled in place; `G/maps/spec.json` gains the fifth panel; README, MIGRATION_HISTORY,
  ROADMAP row G.

## 2. Runs (from `icenine_py/`)

`G` is `benchmarks/eval_reduction/G`. Every run was started as its own plain `nohup caffeinate -i uv run ...` command (the worktree
hook refuses runner scripts and compound commands).

| What | Command | Wall | Label |
|---|---|---|---|
| preflight | `uv run python ../scripts/dev/timing_preflight.py --json G/g3/bfs/preflight.json` | seconds | quiet at the start (no busy process) |
| pilot (20-voxel limit; stopped after 4 seeds: the limit does not cap the BFS expansion) | `pilot_bfs.py --arm mc --variant clean --out G/cache/g3/pilot --name smooth_only_20 --config G/g3/smooth_only_mc_clean.config --max-voxels 20` | about 12 min | contended (a busy editor process; its preflight is kept in `G/cache/` because it records an absolute path) |
| the arm, 2,000 voxels | `pilot_bfs.py --arm mc --variant clean --out G/cache/g3/bfs --name smooth_only --config G/g3/smooth_only_mc_clean.config`; `recon.mic`, `records.json`, `run.json` copied to `G/g3/bfs/runs/smooth_only/` | 11,136 s | contended: one BFS run; other jobs not controlled; the probe runs did not overlap it |
| summary (first pass) | `uv run python scripts/eval_reduction/g3_summary.py` | about 1 min | |
| probe list | `uv run python scripts/eval_reduction/g3_trace_probe.py list` -> `G/g3/trace_probe_voxels.json` (2 voxels) | seconds | |
| probe, key off and on | `g3_trace_probe.py run --spec off --tag base --workers 2`; `... --spec "E3+E4,tau=0.25" --tag e34t25 --workers 2` | about 2 and 3 min | contended (4 workers plus the cache check) |
| cache check (scratch script, not committed; output `G/g3/cache_check.json`) | two full searches (region voxels 4, 19) on a `BFSReconstruction`'s reconstructor with the arm's config; records the cache's builds and hits | about 5 min | contended |
| summary (final), maps | `g3_summary.py`; `uv run python scripts/mic_maps_report.py G/maps/spec.json --out-dir G/maps` (`stats.json` not committed: absolute paths) | about 2 min | |
| doc tables | `scripts/dev/sync_doc_tables.py --doc docs/findings_eval_reduction_2026-10.md --tables G/g3/tables.md G/g2/tables.md G/g1/tables.md G/g1/timing/tables.md`, then `--check` (exit 0) | seconds | |

Compute: the pilot gave about 161 s per seed full search (G2's key-off rerun 116 s on the same voxel, +39%, like G2's integrated
seeds), so the estimate was about 3.2 h against the 3 h budget and the 8 h stop line.

## 3. Results (all numbers from `summary.json` / `tables.md`)

<!-- table:evalG3_total_sets -->
| set | n | wrong_before | wrong_after | fixed_broken | strict_before | strict_after | smooth | equiv_change |
|---|---|---|---|---|---|---|---|---|
| six full-sample voxels (sum) | 6 | 0 (0.0-39.0%) | 0 (0.0-39.0%) | 0 / 0 | 1,992,278 | 2,006,917 | 90,183 | +20.2% |
| Task E's other 24 full-sample voxels (per voxel) | 24 | 1 (0.7-20.2%) | 0 (0.0-13.8%) | 1 / 0 | 330,584 | 335,047 | 15,119 | +21.0% |
| E0 200 voxels (per voxel) | 200 | 41 (15.5-26.6%) | 7 (1.7-7.0%) | 34 / 0 | 48,700 | 49,215 | 252 | +1.7% |
| BFS pilot region, 2,000 voxels (whole run) | 2000 | 37 (1.3-2.5%) | 34 (1.2-2.4%) | 14 / 11 | 25,547,935 | 24,490,331 | 980,594 | +12.4% |
<!-- /table:evalG3_total_sets -->

<!-- table:evalG3_bfs_paired -->
| a | b | wrong_a | wrong_b | fixed | broken | both_wrong |
|---|---|---|---|---|---|---|
| mc_clean | keyoff_rerun | 37 | 37 | 14 | 14 | 23 |
| mc_clean | smooth_only | 37 | 34 | 13 | 10 | 24 |
| mc_clean | outward_only | 37 | 36 | 9 | 8 | 28 |
| mc_clean | integrated | 37 | 42 | 11 | 16 | 26 |
| keyoff_rerun | smooth_only | 37 | 34 | 14 | 11 | 23 |
| keyoff_rerun | outward_only | 37 | 36 | 18 | 17 | 19 |
| keyoff_rerun | integrated | 37 | 42 | 15 | 20 | 22 |
| smooth_only | integrated | 34 | 42 | 5 | 13 | 29 |
<!-- /table:evalG3_bfs_paired -->

<!-- table:evalG3_wrong_by_origin -->
| run | n | another grain | own full search | carried from a wrong full search | boundary crossing | untraced | max_distance_um | full_searches |
|---|---|---|---|---|---|---|---|---|
| stored pilot (2026-10-10) | 37 | 37 | 1 | 5 | 31 | 0 | 9.4 | 67 |
| stored configuration, today's code | 37 | 37 | 2 | 3 | 32 | 0 | 9.4 | 63 |
| SmoothCostRank only | 34 | 34 | 0 | 0 | 34 | 0 | 0.0 | 62 |
| OutwardSearch only | 36 | 36 | 5 | 4 | 27 | 0 | 9.4 | 68 |
| both (Task G2's integrated configuration) | 42 | 42 | 4 | 9 | 29 | 0 | 0.0 | 65 |
<!-- /table:evalG3_wrong_by_origin -->

<!-- table:evalG3_trace_n13 -->
| flip | n | traced_to_arm | traced_to_draw | untraced |
|---|---|---|---|---|
| broken | 11 | 0 | 11 | 0 |
| fixed | 14 | 0 | 14 | 0 |
<!-- /table:evalG3_trace_n13 -->

<!-- table:evalG3_full_search_time -->
| run | full_searches | wall_per_full_search_s | strict_per_full_search | smooth_per_full_search | implied_ms_per_smooth | timing |
|---|---|---|---|---|---|---|
| stored configuration, today's code | 69 | 114 | 331,636 | 0 | n/a | contended: 3 concurrent 2000-voxel BFS runs (Task G2) |
| SmoothCostRank only | 65 | 154 | 335,546 | 15,086 | 2.58 | contended: one 2000-voxel BFS run (Task G3); preflight quiet at the start, other jobs on the machine not controlled |
| both (Task G2's integrated configuration) | 70 | 147 | 317,739 | 14,233 | 2.66 | contended: 2 concurrent 2000-voxel BFS runs plus the E0 run (2 workers) and a 150-voxel pilot at the start (Task G2) |
<!-- /table:evalG3_full_search_time -->

<!-- table:evalG3_full_sample_estimate -->
| run | full_searches | full_search_wall_s | local_fit_wall_s | grain_scale | voxel_scale | estimate_h | timing |
|---|---|---|---|---|---|---|---|
| stored configuration, today's code | 69 | 7,848 | 1,149 | 8.28 | 12.29 | 22.0 | contended: 3 concurrent 2000-voxel BFS runs (Task G2) |
| SmoothCostRank only | 65 | 10,011 | 1,125 | 8.28 | 12.29 | 26.9 | contended: one 2000-voxel BFS run (Task G3); preflight quiet at the start, other jobs on the machine not controlled |
<!-- /table:evalG3_full_sample_estimate -->

Decision points: none in G3 beyond D-20261029-N3 (the arm). The proposals are in the findings doc, "Decisions for the owner".

## 4. What went wrong, and what could be done better

* `--max-voxels` of `pilot_bfs.py` does not limit the BFS expansion, so the pilot ran until it was stopped (4 seeds, enough to
  time a seed).
* The BFS wall times of G2 and G3 are contended at different levels and cannot be compared. A fair time comparison of the key in a
  BFS needs the arms run at the same time, or single-worker.
* The BFS full searches imply about 2.6 ms per smooth evaluation against Task E's 1.639 ms; the cause is open. `pilot_bfs.py`
  could record the cache's builds and hits and the smooth time per full search.
* The levers table reads other tasks' markdown tables by name and row; a renamed column there stops it with an assertion.
