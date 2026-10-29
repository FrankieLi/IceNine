# RUN_LOG: Task G2, the integrated configuration, no-start search and BFS pilot (2026-10-29)

Branch `feature/eval-reduction-G2` off `feature/eval-reduction` (4735f77). Plan: `eval-reduction.md`, "Task G2" (steps 1-3 of
Task G). Approval D-20261028-N2. Commands from `icenine_py/`. Clean images only; `BatchedQuickMC`, `SearchDevice`,
`CandidateDedup` and `FusedCost` off in every arm.

## 1. What changed in the repository, and why

* `icenine/config_file.py`, `icenine/orientation_search.py`: the opt-in key `OutwardSearch` (0/1, default 0), the field
  `outward_search` of `ConfigFile` and `SearchParameters`.
* `icenine/reconstructor.py` (the hook):
  * `AdaptiveVoxelReconstructor.outward_search(vertices, phase, start)`: the Task F fallback policy of D-20261028-N4. F4
    (`RayGaussNewton`, 1-frame window, tau 0.25 deg, acceptance strict cost <= 0.8); where F4 rejects, F5 (`ShellSearch`
    from the start, spacing 0.3 deg, outer shell 6 deg, best 32 by the screen re-ranked by E4, best 2 costed) and then F4
    from F5's point. The settings are Task F's (`outward_bench.py`, `outward_compare.py`), as module constants
    `OUTWARD_*`. The strict cost (pixel radius 0), the screen cost (pixel radius 3, Q-max 5, as Task F's
    `make_cf(.., "screen", "old")`) and the `SmoothCost` are built once per image set (`_outward_costs`; the SmoothCost
    comes from G1's `smooth_cost_cache`, so `SmoothCostRank` and the policy share one build).
  * `BFSReconstruction._local_fit` (BFS neighbours, revisits, restart pass): with the key on the policy runs from the
    inherited start; if it accepts, the local fit (the existing search) starts from the policy's end point, otherwise
    from the inherited start as before. The cma retry (wider sigma0) does not re-run the policy.
  * Reseeds (`_reseed_clusters` -> `_fit_from_seed(start=...)` -> `_outward_reseed`): the policy runs from the reseed
    voxel's current (failed) fit; if it accepts, one local fit from its end point and the usual seed test; if both pass,
    this replaces the full search. If the policy or the seed test rejects, the full search runs as before, so a voxel is
    still final only after a full search (the finality rule's premise).
  * Counting: the policy's strict and screen evaluations go into the voxel's `n_evals` and the phase bucket (as a seed's
    level screens do); the smooth (E4) evaluations, Jacobians and every other count go into `stats["outward_<phase>_*"]`
    and `BFSReconstruction.outward_log` (one entry per policy run). With the key off none of these keys exist.
* `tests/test_outward_bfs.py` (new, 12 tests): config key; key off never calls the policy and the scripted grid still
  matches its golden; accepted / rejected neighbours; the retry; reseeds (policy rejected -> full search; accepted ->
  no full search; seed test rejected -> full search); the real policy on the ThreeVoxels toy equals the hand-built
  composition of `RayGaussNewton` and `ShellSearch` and its counts equal the cost functions' counters.
* `scripts/phase_d/pilot_bfs.py`: `--config` (a config file in place of the arm's), `--name` (run folder), `--sample`
  (another grid); writes `outward_log.json` when the key is on. Default runs are unchanged.
* `scripts/eval_reduction/eval_ledger.py`: records `smooth_evals` (the ranking key's evaluations) when `SmoothCostRank` is on.
* `scripts/eval_reduction/g2_bfs.py` (new): `subregion` (a 150-voxel patch of the pilot region for short pilots),
  `compare` (key-off bit-identity of two pilot runs).
* `benchmarks/eval_reduction/G/integrated_mc_clean.config` (step 1), `G/g2/outward_only_mc_clean.config`.

## 2. Step 1: the configuration, and the decisions it rests on

| Lever | Choice | Source |
|---|---|---|
| overlap rule | lenient (library default) | D-20261013-I6 (`CppOverlapRule` stays off); D-20261012-I3 (the ray rule may become the default, flip not made); D-20261012-N1 (ray rule only at levels 0-2 if used) |
| candidate dedup | off | D-20261014-N4 |
| smooth ranking key, levels 0-2 | `SmoothCostRank E3+E4,tau=0.25`, with G1's cache | D-20261029-N1 |
| finisher | classic MC + VarianceMinimizing (no `LocalOptimizer` key) | finisher study "Decisions on record" (CMA-ES opt-in, off) |
| BFS neighbours and reseeds | `OutwardSearch 1` (F4, then F5 then F4 where F4 rejects; 1-frame window) | D-20261028-N4 |
| BFS | `BFSRevisitRefit 1`, `BFSReseed 1` | as the stored pilot (`phase_d_pilot_reseed_v2/runs/mc_clean`) |
| broken voxels | traced per D-20261014-N13 | D-20261014-N13 |

The overlap rule is the one choice no decision names for Task G. The plan's words "overlap rule as decided for the
coarse levels" were read as: no decision changed the coarse-level rule from the lenient default, and G1 chose the
smooth key on lenient and C++ runs. The ray rule R2 (approved as a possible default, D-20261012-I3) was not tried in G2.

Two implementation choices the plan leaves open, taken before the data:

* **What "accepted" does.** Task F found that F4 stops 0.02-0.05 deg short and "a finisher is still needed after it", so
  an accepted policy end point is the start of the existing local fit (the finisher), not the final orientation.
* **Reseeds.** A reseed's start is the voxel's current fit (a rejected neighbour fit, kept by the Push rule). If the
  policy accepts, one local fit from its end point is judged by the usual seed test; only if that passes does it
  replace the full search. A rejected policy or seed test falls through to the full search, so the per-voxel finality
  rule ("a failed full search is final") keeps its meaning.

## 3. Runs (from `icenine_py/`)

The plan's runner-script style (`nohup bash run.sh`) was refused by the worktree hook ("a construct too complex to
verify"), so every run was started as its own plain command (the commands below).

| What | Command | Wall | Label |
|---|---|---|---|
| key-off check and short pilot, 150-voxel sub-region (`g2_bfs.py subregion --n 150`) | `pilot_bfs.py --arm mc --variant clean --out G/cache/g2/pilot --sample G/cache/g2/subregion_150.mic` with `--name keyoff_base --reseed` (library of 4735f77 by `git archive` on `PYTHONPATH`), `--name keyoff_branch --reseed`, `--name integrated --config G/integrated_mc_clean.config`, `--name outward_only --config G/g2/outward_only_mc_clean.config` | 1013-1512 s | contended (4 runs, then 6 processes with the ledger) |
| key-off comparison | `g2_bfs.py compare G/cache/g2/pilot/keyoff_base G/cache/g2/pilot/keyoff_branch --json G/g2/keyoff_pilot_compare.json` | seconds | |
| step 2, six-voxel ledger | `eval_ledger.py run --voxel V --tag g2_base --dir G/g2/ledger` and `--tag g2_integrated --opt smooth_rank=E3+E4,tau=0.25`, V in the six | 118-180 s per voxel | contended (6 processes) |
| step 2, E0 200 | `smooth_rank_probe.py run --kind e0 --rule old --spec "E3+E4,tau=0.25" --tag g2 --workers 2` | about 15 min | contended |
| step 3, 2,000-voxel pilot region | `pilot_bfs.py --arm mc --variant clean --out G/cache/g2/bfs --name integrated --config G/integrated_mc_clean.config`, the same with `--name outward_only --config G/g2/outward_only_mc_clean.config`, and `--name keyoff_rerun --reseed` (the stored pilot's configuration on this branch, a reproduction check), each under `nohup caffeinate -i` | see section 4 | contended (3 BFS runs at once, plus the E0 run for its first 15 min) |
| summary | `g2_summary.py` | seconds | writes `summary.json`, `tables.md` |

Paths: `G` is `benchmarks/eval_reduction/G`. The 150-voxel `integrated` pilot ran before `pilot_bfs.py` recorded the smooth
evaluations of the full searches, so its `smooth_full_searches` is 0 (not measured), not zero work.

The first two 2,000-voxel runs were started as background tool commands, which stop after 30 minutes; they were killed
after 1-5 minutes and restarted under `nohup` (the integrated run once more, to add the smooth count). No result of
the killed starts is used.
