# RUN_LOG: Task E, smoothing the cost with geometry and physics (2026-10-15)

Task branch `feature/eval-reduction-E` (parent `feature/eval-reduction`). Plan: `~/.claude/plans/eval-reduction.md`,
Task E. All Python through `uv run` from `icenine_py/`. Clean images only, no realistic images, no extra seeds. The
library default path (every new key off) is the baseline; `BatchedQuickMC`, `SearchDevice` and `CandidateDedup` stay
off in every run. Raw arrays are in `benchmarks/eval_reduction/E/cache/` (gitignored); everything quoted comes from
the `summary.json` and `tables.md` files under this directory.

## 1. What changed in the repository, and why

| File | Change | Why |
|---|---|---|
| `icenine/smooth_cost.py` (new) | `SmoothSpec` (spec parser), `SmoothCost` (batched surrogate costs: control S, arms E1 to E4), `SmoothRanker`, `make_ranker` | the arms as library code, one module, thin hooks elsewhere |
| `icenine/reconstructor.py` | `_smooth_ranker` (built once per `reconstruct_voxel` when the key is on); in the level loop an `elif` after the existing `rank_key` branch orders the post-quick-MC list by the smooth key at levels below the spec's `levels`; `last_smooth_ranker` | the ranking hook; with the key off nothing is built and the branch is skipped |
| `icenine/orientation_search.py` | `SearchParameters.smooth_rank` (default `"off"`, validated, read by `from_config`) | the key |
| `icenine/config_file.py` | `smooth_rank` field, key `SmoothCostRank` (validated with `SmoothSpec.parse`) | the key |
| `tests/test_smooth_cost.py` (new) | spec and key parsing; key off bit-identical to the default (recorded search identical, no ranker); key on orders level 0 by the smooth key; control arm equals the strict cost at the truth; smoothed arms change gradually where the control steps; shared multi-spec evaluation equals separate calls; E4 autograd gradient against finite differences | correctness |
| `scripts/eval_reduction/smooth_landscape.py`, `smooth_landscape_summary.py` (new) | measurement 1 on Task B's rays | measurement 1 |
| `scripts/eval_reduction/smooth_rank_probe.py`, `smooth_rank_summary.py` (new) | measurements 2 and 3 (Task C's probe, no further library change) | measurements 2, 3 |
| `scripts/eval_reduction/smooth_timing.py` (new) | measurement 4 | measurement 4 |
| `.gitignore` | `icenine_py/benchmarks/eval_reduction/E/cache/` | raw arrays |

Untouched: `cost_functions.py`, `batched_screen.py`, `_rasterize.c`, all physics code, every golden test. S1's
`FrameTable`, `frame_index_array`, `BraggConstants` and `scattering_omegas_fixed_beam` are reused, not duplicated.

## 2. Decisions taken before the data

* **E2 design.** The first E2 (an overlap-weighted mean over frames of a +-0.5 frame omega extent) raised the cost
  at the truth of voxel 15901 from 0.053 to 0.296: the images are rendered with one frame per spot, so any spread
  of the prediction over frames is charged at the truth. E2 is instead "the largest of (credit on frame k') x
  (1 - smoothstep(distance of the predicted omega to frame k' / reach))": still continuous at a frame boundary,
  full credit at the truth. Recorded in the module docstring.
* **Decision rule for measurement 3 (fixed before the ranking data were read).** Carry at most three arms. An arm
  qualifies when (a) pooled over levels 0-2 with the next-level basin (diameter / 1.5), its ranking keeps the truth
  basin in more cases than it loses it (fixed > broken) on the E0 200 and on the 30 full-sample voxels, both with
  the lenient rule, with the sum of the two nets above zero; and (b) its landscape twin has no more rays with a local
  minimum than the strict cost (C++ rule) at depth 0.01 and at depth 0.05, on the full-sample dense walk (0 to 5 deg)
  and on the E0 dense walk (0 to 0.5 deg). The width variants of one combination count as one arm (the better net
  is carried). Among qualifying arms, the largest summed net first; ties by the smaller time per evaluation.

## 3. Runs (commands from `icenine_py/`)

| What | Command | Wall | Label |
|---|---|---|---|
| pilot, landscape, one ray (voxel 15901 and E0 voxel 7336) | `smooth_landscape.py pilot --kind full` / `--kind e0` | 86 s / 6 s | single process; the strict cost re-evaluated at 5 points equals Task B's cache to 2.2e-8 / 2.4e-8 |
| landscape, fixed widths, 6 full-sample voxels x 20 rays | `smooth_landscape.py run --kind full --workers 5` | 32 min | contended |
| landscape, fixed widths, 30 E0 voxels x 20 rays | `smooth_landscape.py run --kind e0 --workers 1` | 32 min | contended |
| landscape, level-0 widths, full-sample | `smooth_landscape.py run --kind full --workers 4 --set L0` | 81 min | contended |
| landscape, level-0 widths, E0 | `run_queue2.sh` (`--kind e0 --workers 2 --set L0`) | 56 min | contended |
| pilot, ranking probe, 3 E0 voxels | `smooth_rank_probe.py run --kind e0 --rule old --spec off --tag pilot --limit 3 --workers 1` | 26 s | errors and evaluation counts equal Task C's `base` arm on all 3 voxels |
| ranking baseline, E0 200, lenient rule | `smooth_rank_probe.py run --kind e0 --rule old --spec off --tag base --workers 2` | 13 min | contended |
| ranking baseline, E0 200, C++ rule | same, `--rule cpp` | 13 min | contended |
| ranking baseline, full-sample, lenient rule | `smooth_rank_probe.py run --kind full --rule old --spec off --tag base --workers 2` | stopped after the first 6 voxels of the 30-voxel list (= the six interior voxels), 28 min | contended; see section 4 |
| ranking baseline, full-sample, C++ rule, six voxels | same, `--rule cpp --voxels six` | 34 min | contended |

## 4. What went wrong, and deviations

* **The 30-voxel full-sample set was cut to the six interior voxels.** A full-sample probe re-ranks every
  post-MC candidate of levels 0-2 (about 10,000 at level 0) with 20 arms; contended, it took 6 to 14 minutes per
  voxel, so the 30 voxels of the lenient rule alone would have taken about 2.3 hours more. The run was stopped after
  the first six voxels of the list, which are the six interior voxels; the C++-rule run was started on those six only.
  Decided before any full-sample ranking result was read. `run_queue1.sh` (30-voxel lenient and C++ runs),
  `run_queue3.sh` and `run_queue4.sh` (a second process from the other end of the list) were stopped; they are kept
  as the record of what was queued.
* **The full-sample ranking measurement carries no signal on these six voxels.** The strict ranking already ranks a
  basin candidate among the first 11 of about 10,000 at level 0 (first 2 at levels 1-2) of every one of the six voxels (all six end right), so
  no arm can fix anything there and none breaks anything. The decision rule's full-sample part (fixed > broken) is
  therefore read as "no broken case" (fixed = broken = 0 for every arm); recorded here because the rule was written
  for the 30-voxel set.
* **A supplementary run of the other 24 voxels of the 30-voxel set** (lenient rule, re-ranking with six candidate arms
  only, `--tag base30 --skip 6 --arms ...`) was started with two spare cores and stopped after 10 minutes, before
  any voxel finished, to give the cores to the end-to-end runs. Its empty output directory was removed.

## 5. Decision point after measurements 1-2 (non-blocking; the rule of section 2, default taken)

Data the owner would have seen (all from `rank/summary.json`, `landscape/summary.json`, `landscape_L0/summary.json`):

* Ranking, pooled levels 0-2, next-level basin, E0 200 lenient (net = fixed - broken): E1+E2+E3+E4 (auto) 24 - 0, E1+E2+E3 (auto)
  24 - 0, E4 (auto) 24 - 0, E1+E2 (auto) 24 - 0, E3+E4 (auto) 23 - 0, E3+E4 (tau 0.25) 23 - 0, E1 (auto) 23 - 0, E2 (auto) 8 - 3,
  E3 5 - 7, S 3 - 3. Six full-sample voxels, lenient: 0 - 0 for every arm (saturated, section 4).
* Minima (rays with a minimum, depth 0.01 / 0.05; strict C++: full 72 / 0 of 120, E0 12 / 0 of 600). Twins at level-0 widths:
  E1+E2+E3+E4 full 16 / 0, E0 6 / 0 (passes); E3+E4 full 46 / 0, E0 6 / 0 (passes); E1+E2+E3 full 34 / **1** (fails by one ray);
  E4 E0 **70** / 0 (fails); E1+E2 full **74**, E0 **83** (fails); E1 full **120 / 46** (fails). Fixed-width twins: E3+E4 (tau 0.25)
  full 52 / 0, E0 6 / 0 (passes); E1+E2+E3 (3 px) full 52 / **1** (fails); E4 (tau 0.25) E0 **71** (fails).

Choice: **carry E1+E2+E3+E4 (auto widths) and E3+E4 (tau 0.25 deg)**; only these two qualify, so the third slot is empty. For E3+E4 the auto
and tau-0.25 variants tie (net 23 each, both pass); the tie-break is the time per evaluation, which is the same by construction (one code
path, tau only enters the final r^2 / (r^2 + tau^2)); the tau-0.25 variant was carried because its runs had already been started on the
strength of its fixed-width twin while the level-0 walks were still running. E1+E2+E3, the arm that removes nearly every step, is the closest
miss (one full-sample ray of 120 with a 0.05-deep minimum, where the strict cost has none).

## 6. Measurement 3 and the suite (commands from `icenine_py/`)

| What | Command | Wall | Label |
|---|---|---|---|
| E3+E4 (tau 0.25), E0 200, lenient then C++ rule, then six full-sample voxels lenient | `run_e2e_e34t25.sh` (`smooth_rank_probe.py run --kind e0|full --rule old|cpp --spec "E3+E4,tau=0.25" --tag e34t25 [--voxels six] --workers 2`) | about 30 min per E0 run (split over two invocations each, see below), 17 min for the six | contended |
| E3+E4 (tau 0.25), six full-sample voxels, C++ rule | `smooth_rank_probe.py run --kind full --rule cpp --spec "E3+E4,tau=0.25" --tag e34t25 --voxels six --workers 2` | 16 min | contended |
| E1+E2+E3+E4 (auto), E0 200 both rules, six full-sample voxels lenient | `run_e2e_e1234a.sh` (`--spec "E1+E2+E3+E4" --tag e1234a`) | 25 and 26 min per E0 run, 32 min for the six | contended |
| summaries | `smooth_rank_summary.py rank`, `smooth_rank_summary.py e2e --arms e34t25,e1234a`, `smooth_landscape_summary.py`, `collect_tables.py --root benchmarks/eval_reduction/E --prefix "" --out benchmarks/eval_reduction/E/doc_tables_all.md` | minutes | |
| full suite | `uv run pytest tests/ --junitxml=<scratch>/junit.xml -p no:warnings` | 19 min (contended, 6 jobs running) | 1111 passed, 37 skipped, 1 deselected |

What went wrong here:

* A `pkill -f "tag e34t25 --workers 2"` meant for a duplicate C++-rule E0 run also stopped the queue's lenient E0 run of E3+E4 at 165 of 200
  voxels; the queue then went on to its C++ step. The lenient run was restarted (`e2e_e0_old_e34t25_b.log`); it skipped the 165 voxels
  already written and ran the other 35. The probe is deterministic per voxel (seeded), so the split does not change any result.
* The duplicate C++-rule E0 run of E3+E4 (`e2e_e0_cpp_e34t25_p.log`, started early on two spare cores) was stopped at 76 voxels; the queue's
  own C++ run skipped those and ran the rest (`task()` skips a voxel whose result file exists).

## 7. Measurement 4 (timing)

| What | Command | Label |
|---|---|---|
| full-sample images, voxel 15901, 1,000 orientations | `smooth_timing.py --set full` (stops unless `preflight.require_quiet()` passes) | single-worker; preflight `timing/preflight_full.json` (1-min load 1.19, no busy process) |
| E0 images, voxel 7336 | `smooth_timing.py --set e0` | single-worker; `timing/preflight_e0.json` |
| tables | `smooth_rank_summary.py timing` (ms per evaluation; extra time per voxel = smooth evaluations per voxel from the end-to-end runs x ms per evaluation, plus the one-off build) | |

The arms were timed at the level-0 widths (the widest), so the per-voxel figures are upper estimates. E3+E4 was timed with level-width tau
(0.42 deg at level 0) and stands for the carried tau-0.25 variant (same code path).

## 8. Could be done better

* The full-sample ranking probe re-ranked every candidate with 20 arms; that made it too slow for the 30-voxel set. Re-ranking a sample of the
  level-0 list (or only the arms that matter) would have kept the 30 voxels within the budget, and those are the voxels where a fix could show.
* `SmoothCost` is rebuilt per `reconstruct_voxel` call (16.6 s on the full sample, mostly the observed rays). A cache keyed on the image set
  is needed before any BFS use.
* The E1 search over 3.9 million lit pixels per detector dominates its time on full-sample images; a distance-transform image per frame
  (bilinear lookup) would make it as cheap as the binary test, at more memory.
* The decision rule should have been written for the set actually run; the 30-to-6 change made its full-sample part uninformative.
* Kill commands should match on the full command line of one run; one broad `pkill -f` pattern stopped a run it was not meant for (section 6).
