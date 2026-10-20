# RUN_LOG: Task F, searching outward from a start (2026-10-20)

Task branch `feature/eval-reduction-F` (parent `feature/eval-reduction`). Plan: `~/.claude/plans/eval-reduction.md`,
Task F. All Python through `uv run` from `icenine_py/`. Clean images only, no extra seeds. `BatchedQuickMC`,
`SearchDevice` and `CandidateDedup` stay off (no run here touches the reconstructor). Raw arrays are in
`benchmarks/eval_reduction/F/cache/` (gitignored); everything quoted comes from the `summary_*.json`, `tables_*.md`,
`pilot/`, `timing/` and `diag/` files under this directory.

## 1. What changed in the repository, and why

| File | Change | Why |
|---|---|---|
| `icenine/outward_search.py` (new) | `RayGaussNewton` (F4), `ShellSearch` (F5), `shell_rotations`, `shell_size`, `shell_directions`, `screen_fn_from`, `smooth_fn_from`, `angle_deg`, `rotate_left` | the two new methods, library only |
| `tests/test_outward_search.py` (new) | shell geometry and coverage, `ShellSearch` on a synthetic cost (counts, early stop, not found, no smooth stage), `RayGaussNewton` on the ThreeVoxels toy (Jacobian against finite differences, recovery from 1 and 3 deg, strict cost evaluated twice, callback, frame window, e4_start regression) | correctness |
| `scripts/eval_reduction/outward_bench.py`, `outward_summary.py`, `outward_pilot.py`, `outward_diag.py` (new) | benchmark, summary (tables, figures, choice rule, policy), settings pilots, failure diagnostic | measurements |
| `.gitignore` | `icenine_py/benchmarks/eval_reduction/F/cache/` | raw arrays |
| `README.md`, `MIGRATION_HISTORY.md`, `docs/ROADMAP.md` (row F) | rows, section, status | docs |

Untouched: `cost_functions.py`, `orientation_search.py`, `batched_screen.py`, `reconstructor.py`, `_rasterize.c`, all
physics code, every golden test. `SmoothCost` (Task E) and `ScreenEvaluator` are used, not edited.

## 2. Decisions taken before the data

* **Acceptance test.** A run is accepted when the strict cost of its end point is 0.8 or less. Chosen on the pilot
  voxels (section 3): the strict cost at the truth is 0.06 to 0.1, at 0.11 deg off it is half-way to 1 (E0), and a
  shell point within the spacing of the truth costs 0.5 to 0.9, so a tighter test (0.5) leaves F5 without an accepted
  shell in most runs at 0.3 deg spacing.
* **Capture thresholds.** Within 1, 0.1 and 0.05 deg of the truth (the plan asked for 0.05 and 1). A tolerance of
  0.001 deg is applied (`CAPTURE_TOL`): R_true is float32, so a start exactly 1 deg out measures 0.9999 deg and would
  count as captured at 1 deg with no work.
* **Choice rule (decision point).** The share within 1 deg per 1,000 strict-equivalent evaluations to stop, at 2 and
  5 deg; ranks at the two offsets added; the two lowest sums are carried; ties by the larger mean score. The
  candidates are the plan's methods (F1, F2 in three variants, F4, F5, F5+F4); F3 cannot run on the full sample.
  Written before the full-sample runs were read.
* **Counting.** Raw (every evaluation once) and strict-equivalent (weights from the single-worker timing run).
* **Voxel sets.** E0: all 8 near-axis voxels (r_perp below 120 um) of the 200 with clean images, and 16 far-from-axis
  (r_perp 250 um or more) spread evenly; 4 more (2 near, 2 far) as pilot voxels, not used in the benchmark. Full
  sample: the six interior voxels. 20 directions (seed 7000 + voxel, the same directions at every offset).
* **F6 not run.** The plan allows it only as a reference at the smooth costs "if the per-run cost is small"; Task E's
  smooth costs are not wired into the Adam optimizer, and Phase C has it behind CMA-ES from 1.5 deg.

## 3. Pilots (settings, on the four pilot voxels; `pilot/`)

| What | Command | Result |
|---|---|---|
| F4 robust-scale schedules | `outward_pilot.py schedule` | one scale of 0.25 deg is best or equal at every offset and the cheapest (`pilot/f4_schedule.txt`) |
| F5 spacing and acceptance | `outward_pilot.py spacing` | 0.3 deg with acceptance 0.8 found a passing shell in all 16 runs at 2 and 5 deg; 0.45 deg left 12% (`pilot/f5_spacing.txt`) |
| F4 frame window | `outward_pilot.py window` | from 3 frames upward 100% within 1 deg from 5 deg; 5 frames chosen (`pilot/f4_window.txt`). Run after the main run showed the failure (section 6). |

The strict cost at the truth, the screen's half-way rise (0.78 deg on E0, from Task B's cache, the median over the
walks) and the strict cost's half-way rise (0.11 deg on E0, 0.33 deg full sample) were read from Task B's ray walks
(`benchmarks/eval_reduction/B/cache`).

## 4. Runs (commands from `icenine_py/`)

| What | Command | Wall | Label |
|---|---|---|---|
| pilots | `outward_pilot.py schedule|spacing|window` | 2 to 6 min each | one process, others running (contended) |
| E0, 24 voxels, 20 dirs, 6 offsets, all methods | `benchmarks/eval_reduction/F/run_all.sh` (`outward_bench.py run --kind e0 --workers 5 --dirs 20 --n-near 8 --n-far 16 --tag e0`) | 891 s | contended (5 workers) |
| full sample, six voxels, offsets 1, 2, 5 deg | `outward_bench.py run --kind full --workers 5 --dirs 20 --offsets 2,3,5 --tag full --methods f1,f2_05,f2_1,f2_2,f4,f5,f54` | 191 s | contended |
| F4 with the 5-frame window, E0 and full | `outward_bench.py run ... --tag e0_f4w --methods f4w` (and `full_f4w`) | 2 min | contended |
| timing, E0 (2 voxels, 5 dirs, offsets 1, 2, 5 deg) | `outward_bench.py timing --kind e0 --dirs 5` | 1 min | **single-worker**, preflight passed (`timing/preflight_e0.json`); the first attempts stopped on a busy preflight (load) and were retried |
| timing, full (voxel 15901) | `outward_bench.py timing --kind full --dirs 5` | 25 s | **single-worker**, `timing/preflight_full.json` |
| diagnostic of F4's failures | `outward_diag.py e0 --n-voxels 12`, `outward_diag.py full` | 5 min | contended |
| summaries | `outward_summary.py e0`, `... full`, `... merge` | 1 min | |
| reproducibility | the E0 run was done twice (once with a tolerance of 1e-6, then 1e-3): all 13,536 non-timing arrays are identical; a third check of one voxel after the library fix (`_weights`) is identical too | | |

## 5. What went wrong

* **The 1 deg column of "within 1 deg".** The first summary counted a start exactly 1 deg out (measured 0.9999 deg
  because R_true is float32) as captured at 1 deg. Fixed with a 0.001 deg tolerance, and both runs repeated (15 + 3
  minutes). After the fix the 1 deg column of Table 1 counts only runs that moved in.
* **First `e4_value` call of a fresh solver.** `RayGaussNewton._weights` read the detector size before the observed
  rays (which set it) were built, so the very first E4 value of a new solver was 1.0 and `e4_start` was wrong. No
  table is affected (the cost used in the loop is the second call; `e4_start` is not tabulated). Fixed, with a
  regression assertion; one E0 voxel re-run with the fixed library gave identical arrays.
* **The E0 failure mode was understood after the run.** F4 with one frame of matching window loses 5% of the E0 runs
  from 5 deg. The window variant (f4w) was added afterwards, chosen on pilot voxels, and reported as an addition, not as
  one of the plan's methods; the carry rule was taken without it.
* **The findings doc.** The write of `docs/findings_outward_search_2026-10.md` was refused by the tool ("Subagents should
  return findings as text, not write report files"). The draft is returned in the report; see the report.

## 6. Why F4 fails from far starts (`diag/`)

`diag/e0_w1.json`: at 5 deg, 10 of 240 runs on 12 E0 voxels end more than 1 deg out; they end with an E4 cost near
1 (almost no peak has a ray within one frame), the share of peaks that start on the truth's ray is 0.02 against 0.25
in the captured runs. With a 5-frame window all 240 are captured (`diag/e0_w5.json`). On the full sample the same
holds partly: 39 of 120 not captured with one frame, 31 of 120 with five; a voxel has a median of 217.5 attached
rays for 58 predicted peaks (E0: 58 for 58), consistent with a wide window picking up other grains' rays (untested).

## 7. Could be done better

* The Jacobian (autograd, one batched backward pass) costs about 26 strict evaluations. Three forward differences
  would cost about 3 E4 values (about 5 strict evaluations) and are accurate enough for Gauss-Newton; forward-mode
  autograd was tried and was slower (13 ms against 9 ms on one voxel). Not done, because the plan asks for an
  analytic Jacobian.
* F5 spends most of its screen evaluations on shells the truth is not in. A coarse-to-fine spacing would cut that.
* The timing runs cover two E0 voxels and one full-sample voxel; more would tighten the strict-equivalent weights.
* F3 on the full sample needs windows cut from the image stack (not built).
* The policy "F4, then F5 and F4 on rejection" is in-sample; it needs a run on voxels it was not read from (Task G's
  other 24 full-sample voxels).
