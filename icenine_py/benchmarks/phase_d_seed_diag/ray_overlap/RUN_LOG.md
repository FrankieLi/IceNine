# RUN_LOG: ray-consistency overlap rule, steps 1-2, 2026-10-11

Task branch `feature/finisher-mc-study-ray-overlap` (parent `feature/finisher-mc-study`), plan section "Phase D task:
ray-consistency overlap rule". Steps 1 (peak level) and 2 (coarse stage) are done; step 3 (E0 end to end, 6
full-sample voxels) is NOT run, as instructed. All timings are "contended" (several workers on the shared laptop; the
`preflight.json` here was taken while the step-1 job was already running, so it shows the load of that job).

## 1. What changed

| File | Change |
|---|---|
| `icenine/ray_overlap.py` | new: `RayRule`, `ray_flags` (observed ray cast through the bright detector-0 pixels, tolerance window arms R and A), `apply_cases_python` (twin of the C block), `ray_rule_from_config` |
| `icenine/_rasterize.c` | `stage_d_overlap` takes optional `obs_in`, `obs_bright`, `ray_f1`; cases a-f applied after the per-detector loop; returns 8 case counters as an extra tuple only when the arrays are given; `HAS_RAY_RULE` |
| `icenine/cost_functions.py` | `OverlapInfo.ray_cases`; `calculate_diffraction_overlap_batched(ray_rule=, voxel_size=)` (implies `cpp_rule`; fallback Python loop too); `VoxelCostFunction(ray_rule=)` |
| `icenine/reconstructor.py` | level loop: `_ray_rule_for_level(level)`; levels 0-2 get a ray global cost fn, a ray local cost fn (re-cost) and their own `MCOptimizer` (same rng stream) for the quick MC; FindOptimal, VarianceMinimizing and the final cost keep the strict local cost fn; local evaluation count includes the coarse local fns |
| `icenine/orientation_search.py`, `config_file.py` | `SearchParameters.ray_rule`; config keys `RayOverlapRule off|R<n>|A` |
| `tests/test_ray_overlap.py` | 8 tests |
| `scripts/phase_d/ray_overlap_study.py`, `ray_overlap_summary.py` | steps 1-2 and the summary |

Default off; library behaviour unchanged. Update 2026-10-12: case f1 and the `RayOverlapF1` key were dropped by owner decision (f1 kept 0.966 vs f0 0.994 of evaluations at the truth, arm A, clean screen, step 1; no benefit in steps 1-2; may revisit); case f is always forgiven. The result files and tables keep the f1 rows as the historical record. `ray_flags` now has a C version (`_rasterize.ray_flags`), see the speed section below. `peak/` (17 MB of per-evaluation rows, regenerated in about 1 min by `run --step 1`) and `cache/` are not committed. The extension must be rebuilt (`HAS_RAY_RULE`); it is gitignored.

## 2. Decisions the plan left open (flag if wrong)

- Ray origin: voxel centroid at the peak's omega. Detector-0 "bright pixels the spot overlaps": bright pixels of the
  +-3 window (screen) or of the bounding box of the projected triangle (local cost, radius 0), at most 64.
- n_det = 2 here; for n_det > 2 the same straight line is intersected with each later detector (not re-aimed through
  bright pixels at each stage). Untested.
- Arm A: every term is a half-width; d_omega = range-map step (2 deg); r = distance of the voxel from the rotation axis.
- Predicted out of bounds with the observed ray bright: left as the C++ rule (counted as `f_bright`).
- "Levels 0-2": max_local_resolution is 3, so there are levels 0-3; level 3 uses the strict cost.
- When `CppOverlapRule` is 0 the strict part of the search is the OLD rule; the experiments here set it on for the
  `cpp` and ray arms and off for the `old` arm.
- Owner review 2026-10-12: the half-width reading and levels 0-2 confirmed; no objection to the rest. Arm A terms are half-widths (intended reading); windows stay as
  set, untuned ("WINDOW DECISION" note in ray_overlap.py and _rasterize.c); R2 and A go to step 3; level 3 stays strict.

## 3. Runs

| Run | Command (from `icenine_py/`) | Workers | Wall |
|---|---|---|---|
| pilot (2 voxels, 3 arms) | `ray_overlap_study.py pilot` | 1 | 607 s step 2 for one voxel x 3 arms |
| step 1, 30 voxels x 2 variants | `run --step 1 --workers 10` (`cache/run_step1.sh`) | 10 | about 1 min |
| step 2, 12 voxels x 2 variants x 12 arms | `run --step 2 --workers 10 --limit 12` (`cache/run_step2.sh`) | 10 | about 2.5 h, contended |

Tests: new file 8 passed; full suite recorded in section 5.

## 4. Results (tables in `tables.md`, numbers in `summary.json`; reading in MIGRATION_HISTORY "Ray-consistency overlap rule")

Step 1 (30 voxels, clean images, screen, pooled peaks): valid share of detector-0 overlaps at random orientations 0.209
(C++ rule), 0.244 (R1) to 0.435 (A, f0); at the truth 0.994 (C++), 0.994 (A, f0), 0.966 (A, f1); at 1 deg 0.656 (C++),
0.902 (R1), 0.927 (A, f0). The ray rule filters LESS than the C++ rule and more than the old rule.
Step 2 (12 voxels): level-0 candidates 10240 old, 7832 C++, 8207 to 9336 ray arms (clean); the truth basin survives
the cut in 12 of 12 voxels in all arms on clean images, 11 of 12 in two realistic cells (3 deg threshold effect, the
closest candidate was inside the kept set); 12 voxels cannot separate the arms (interval 0.76 to 1.00).

## 5. Verification

- C stage equals the Python twin on random input, both f variants (test).
- Cases a-f on hand-built C-stage input (test); geometry on synthetic parallel detectors (test).
- Full suite: `uv run pytest tests/`: 864 passed, 34 skipped, 1 deselected.
- Number audit: `audit_numbers.py --doc ../../../../MIGRATION_HISTORY.md --section "Ray-consistency overlap rule" --sources summary.json tables.md`.

## 6. Could be done better

1. Step 2 uses 12 voxels (compute was about 100 times the plan's estimate); 30 would narrow the survival interval.
2. The truth-basin criterion (3 deg after the quick MC) is noisy; a rank-of-closest-candidate measure on the
   discrete-stage candidates (before the MC moves them) would be cleaner, and a 1-2 deg threshold would test more.
3. Step 1 pools peaks over voxels; a voxel-level interval (bootstrap over voxels) is not computed.
4. (done 2026-10-12) `ray_flags` in C; see the speed section below.
5. No separate "rule only at the screen" versus "rule at the re-cost and quick MC" split: the plan applies it to all of
   levels 0-2, so which of the three uses carries the effect is untested.
6. Arm A's terms are untuned guesses; the measured case-f rate at the truth (0.029 of events) suggests the window may
   still be too small for f1 or that f1 is simply too strict.
7. Case shares use the screen cost function; the local-cost case table is in `summary.json` only.
8. n_det > 2 chaining untested.

## 7. Speed work, 2026-10-12 (C ray check; profile)

Single-worker, quiet machine (`preflight_speed.json`; load about 1.2 = this job), clean images, voxels 22867, 15901,
17034, levels 0-2, seconds per voxel (`speed_timings.json`):

| arm | before (Python ray check) | after (C ray check) |
|---|---|---|
| old | 162 / 164 / 168 (unchanged code) | same |
| cpp | 117 / 116 / 118 (unchanged code) | same |
| f0_R2 | 293 / 294 / 289 | 132 / 131 / 132 |
| f0_A | 340 / 337 / 335 | 149 / 149 / 148 |

"Before" is the Python reference `ray_flags` (switch `cost_functions.USE_C_RAY_FLAGS = False`) in the same build;
the per-level results (candidate counts, distances, ranks, evaluation counts) are identical before and after for both
ray arms on all three voxels. In a separate check 5784 evaluations (3 voxels, screen and local, arms R1, R2, R8, A)
gave identical flags, costs and case counts with the C and the Python check.

Profile (cProfile, cpp arm, voxel 22867, 133 s inflated by the profiler; 233,546 cost evaluations of one orientation each:
96,338 in the discrete screen, 137,208 in the quick MC): `calculate_diffraction_overlap_batched` 86 s cumulative of
which its own Python 46 s, `get_scattering_omegas_torch` 15 s, torch small-tensor ops (sum, reshape, einsum, where) about
15 s, the C `stage_d_overlap` 7 s (5%). So the baseline coarse stage is bound by Python / torch overhead of about
0.5 ms per single-orientation evaluation, not by the overlap test. The rule's cost after the C check is 12% (R2) and
27% (A) over the cpp arm. Not done: batching many orientations per call (the screen evaluates orientations one by
one inside `run_discrete_search`; the MC is sequential by construction). A batched screen would change float32
matmul shapes and is not guaranteed bit-identical, so it was left out.

## 8. 30-voxel coarse rerun and step 3, 2026-10-12

Both runs follow the owner decisions of 2026-10-12 (f1 dropped, windows untuned, ray rule at levels 0-2, arms old, cpp,
R2, A). All wall times are contended (10 workers) unless marked single-process; results in `coarse30/`, `summary30.json`,
`tables30.md` (reading in MIGRATION_HISTORY "Ray-consistency overlap rule", last subsection), `e2e_summary.json`,
`e2e_tables.md`, `preflight_step3.json`.

| Run | Command (from `icenine_py/`) | Workers | Wall |
|---|---|---|---|
| step 2 again, 30 voxels x 2 variants x 4 arms (all 240 re-run, to get the new fields) | `RAY_COARSE_DIR=coarse30 ray_overlap_study.py run --step 2 --workers 10 --arms old,cpp,f0_R2,f0_A` (`cache/run_step2_30.sh`) | 10 | about 75 min, contended |
| step 3a: six full-sample voxels, arms R2 and A (old and cpp reused from `runs/base_full_*`, `cpprule_full_*`) | `seed_diag.py run --voxel V --images full --tag rayR2|rayA --opt ray_rule=R2|A` (`scripts/phase_d/cache/ray_step3.sh`) | 1 each, quiet preflight | 12 runs, about 29 min, 117 to 158 s each (single-process) |
| step 3b: E0, 200 voxels x 2 variants, R2 and A (old and cpp reused from `cache/off_ms0`, `on_ms0`) | `e0_cpprule.py run --arm rayR2|rayA --msin 0 --workers 10` | 10 | 1124 s and 1149 s, contended |

Projection before step 3: about 17 min per E0 arm for old and cpp (stored), ray arms 12% to 27% more, about 1.2 h in all,
under the 3 h limit. Code changes: `ray_overlap_study.py` (pre-MC rank, survival at 1 and 2 deg, `RAY_COARSE_DIR`),
`e0_cpprule.py` (arms rayR2, rayA), `seed_diag.py` (`--opt ray_rule=`), two summary scripts.

Could be done better (additions): (a) old and cpp arms of step 3 are the stored runs, not re-run in this build (the code
path of those arms is unchanged, but not re-verified); (b) E0 is one seed; (c) the survival criterion at 1 deg at level
0 mostly measures the grid (5 of 30 voxels have a level-0 candidate within 1 deg); (d) no per-case or per-window
decomposition of why the arms differ (R2 has the smaller level-0 rank, A the fewer wrong on realistic E0); (e) the six
full-sample voxels are all right in every arm, so they measure cost only; (f) still untested: rule only at screen vs
also at re-cost and quick MC; n_det > 2.
