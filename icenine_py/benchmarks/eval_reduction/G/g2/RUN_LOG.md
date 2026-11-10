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

Compute: the pilot's estimate (about 0.6 h of policy time added to a run of about 3.8 h, plus about 35% on every full
search for the smooth key) was inside twice the 4 h budget, so the stop rule did not apply. The three 2,000-voxel runs
took 2.5 h (key off), 3.1 h (OutwardSearch only) and 3.8 h (integrated), in parallel, contended. Full suite:
1448 passed, 41 skipped, 1 deselected (`uv run pytest tests/ --junitxml=<scratch>/junit.xml`, 619 s).

## 4. Results (all numbers from `summary.json` / `tables.md`)

**Key off is bit-identical.** The 150-voxel sub-region run with the base commit's library and with this branch: same
`recon.mic`, records and stats (`keyoff_pilot_compare.json`: 4 seeds, 4 reseeds, 146 neighbour fits, 123 revisits). The
six-voxel ledger with the key off equals Task A's base ledger in every voxel (evaluations and error).

**The stored pilot is not today's code.** The stored pilot's configuration run on this branch (`keyoff_rerun`) is wrong
on 37 voxels like the stored pilot, but by voxel 14 fixed and 14 broken, with 40 seeds instead of 34. The stored pilot
dates from 2026-10-10; the centroid screen became the default on 2026-10-12. This is consistent with that change (not
tested which change it is). So the like-for-like baseline of step 3 is `keyoff_rerun`; the stored pilot is also shown.

Step 2, no-start search, six full-sample voxels (Task A's ledger; lenient rule):

<!-- table:evalG2_ledger_stages -->
| stage | task_a | integrated | change |
|---|---|---|---|
| screen | 465057 | 466695 | +0.4% |
| recost | 402506 | 404295 | +0.4% |
| quick | 1098888 | 1110036 | +1.0% |
| find | 17312 | 15792 | -8.8% |
| variance | 8509 | 10093 | +18.6% |
| final | 6 | 6 | +0.0% |
| all | 1992278 | 2006917 | +0.7% |
| smooth (E3+E4, outside the strict count) | 0 | 90183 | n/a |
<!-- /table:evalG2_ledger_stages -->

Step 2, E0 200 (paired against Task E's stored baseline; broken voxels traced per N13):

<!-- table:evalG2_e0 -->
| n | wrong_base | wrong_integrated | fixed | broken | both_wrong | broken_traced | evals_base | evals_integrated | smooth | identical_to_task_e |
|---|---|---|---|---|---|---|---|---|---|---|
| 200 | 41 (15.5-26.6%) | 7 (1.7-7.0%) | 34 | 0 | 7 | 0 / 0 / 0 | 48700 | 49215 | 252 | 200 of 200 |
<!-- /table:evalG2_e0 -->

Step 3, the 2,000-voxel pilot region (clean images, mc, revisit and reseed on, seed order 0):

<!-- table:evalG2_bfs_accuracy -->
| run | n | wrong | unresolved | grains_lost | median_err_right | seeds | reseeds | neighbor_fits | accepted_first | revisits |
|---|---|---|---|---|---|---|---|---|---|---|
| mc_clean | 2000 | 37 (1.3-2.5%) | 1 | 0 | 0.0118 | 34 | 35 (32 accepted) | 1966 | 770 | 2126 (1165 accepted) |
| keyoff_rerun | 2000 | 37 (1.3-2.5%) | 2 | 0 | 0.0101 | 40 | 29 (28 accepted) | 1960 | 769 | 2089 (1167 accepted) |
| integrated | 2000 | 42 (1.6-2.8%) | 0 | 1 | 0.0077 | 38 | 32 (32 accepted) | 1962 | 762 | 2140 (1173 accepted) |
| outward_only | 2000 | 36 (1.3-2.5%) | 1 | 1 | 0.0113 | 34 | 36 (35 accepted) | 1966 | 764 | 2166 (1168 accepted) |
<!-- /table:evalG2_bfs_accuracy -->

<!-- table:evalG2_bfs_paired -->
| a | b | wrong_a | wrong_b | fixed | broken | both_wrong |
|---|---|---|---|---|---|---|
| mc_clean | keyoff_rerun | 37 | 37 | 14 | 14 | 23 |
| mc_clean | integrated | 37 | 42 | 11 | 16 | 26 |
| mc_clean | outward_only | 37 | 36 | 9 | 8 | 28 |
| keyoff_rerun | integrated | 37 | 42 | 15 | 20 | 22 |
| keyoff_rerun | outward_only | 37 | 36 | 18 | 17 | 19 |
<!-- /table:evalG2_bfs_paired -->

<!-- table:evalG2_bfs_evals -->
| run | seed | reseed | neighbor | revisit | total | smooth_full_searches |
|---|---|---|---|---|---|---|
| mc_clean | 11274936 | 11598290 | 1290496 | 1398836 | 25562558 | 0 |
| keyoff_rerun | 13264873 | 9617992 | 1293861 | 1371209 | 25547935 | 0 |
| integrated | 12747993 | 9493718 | 19736111 | 18843288 | 60821110 | 996303 |
| outward_only | 11265144 | 10415877 | 19987665 | 19950797 | 61619483 | 0 |
<!-- /table:evalG2_bfs_evals -->

<!-- table:evalG2_bfs_policy -->
| run | phase | runs | accepted | fallback | fallback_accepted | strict | screen | smooth | jac | strict_equiv | wall_s | median_wall_s | max_wall_s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| integrated | neighbor | 1962 | 1523 | 455 | 16 | 23028 | 18424908 | 325119 | 13461 | 2391401 | 1218 | 0.11 | 3.0 |
| integrated | reseed | 32 | 30 | 2 | 0 | 150 | 83314 | 1836 | 163 | 14561 | 7 | 0.09 | 2.6 |
| integrated | revisit | 2140 | 1722 | 425 | 7 | 22283 | 17412731 | 312038 | 14757 | 2325773 | 1175 | 0.11 | 3.1 |
| outward_only | neighbor | 1966 | 1521 | 460 | 15 | 23290 | 18674835 | 325889 | 12560 | 2385335 | 1239 | 0.10 | 17.8 |
| outward_only | reseed | 36 | 33 | 3 | 0 | 201 | 124971 | 2599 | 209 | 20436 | 10 | 0.11 | 2.5 |
| outward_only | revisit | 2166 | 1722 | 449 | 5 | 23444 | 18495783 | 327026 | 14220 | 2420169 | 1237 | 0.11 | 3.0 |
<!-- /table:evalG2_bfs_policy -->

<!-- table:evalG2_bfs_wall -->
| run | seed | reseed | neighbor | revisit | total |
|---|---|---|---|---|---|
| mc_clean | 6025 | 6213 | 723 | 783 | 13744 |
| keyoff_rerun | 4553 | 3295 | 559 | 590 | 8997 |
| integrated | 5951 | 4325 | 1770 | 1775 | 13822 |
| outward_only | 3852 | 3489 | 1791 | 1848 | 10980 |
<!-- /table:evalG2_bfs_wall -->

Broken voxels traced (`summary.json` `bfs.paired`): against `keyoff_rerun`, 8 of the integrated run's 20 breaks and
8 of the OutwardSearch-only run's 17 end on the orientation of a reseed that took the outward-search path (traced to the
arm); 12 and 9 are untraced (they count against the arm). Every reseed that took that path ended wrong: 4 of 4
(integrated) and 5 of 5 (OutwardSearch only), 33 to 55 deg out. One of them is the reseed of grain 331, a one-voxel
grain in the region, which is lost in both arms.

## 5. What went wrong, and what could be done better

* The Task F acceptance test (strict cost at most 0.8) was chosen on starts within 5 deg of the truth. A BFS start is
  often a neighbouring grain's orientation, which on full-sample images can already score below 0.8 (a hit ratio of
  0.88 was seen). F4 then "accepts" a wrong grain. For neighbours this does no visible harm (the local fit and the BFS
  relative test follow), but a reseed starts from a failed fit by construction, and the absolute seed test (0.8) let
  5 of 36 such starts through, all wrong.
* F5 rarely helps a BFS neighbour: of 909 fall-backs in the OutwardSearch-only run (neighbours and revisits), 20 were
  accepted; each fall-back that finds nothing screens the whole 6 deg ball (about 2.4 s each, contended).
* The plan's runner-script pattern is refused by the worktree hook; background tool commands stop after 30 minutes, so
  long runs must be started with `nohup` directly.
* The integrated 150-voxel pilot ran before the smooth count was added; its smooth evaluations are not measured.
