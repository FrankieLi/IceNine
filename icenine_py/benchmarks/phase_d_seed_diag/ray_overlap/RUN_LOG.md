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
| `icenine/orientation_search.py`, `config_file.py` | `SearchParameters.ray_rule`; config keys `RayOverlapRule off|R<n>|A`, `RayOverlapF1` |
| `tests/test_ray_overlap.py` | 8 tests |
| `scripts/phase_d/ray_overlap_study.py`, `ray_overlap_summary.py` | steps 1-2 and the summary |

Default off; library behaviour unchanged. `peak/` (17 MB of per-evaluation rows, regenerated in about 1 min by `run --step 1`) and `cache/` are not committed. The extension must be rebuilt (`HAS_RAY_RULE`); it is gitignored.

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
4. `ray_flags` is a per-peak Python loop (about 0.5 ms extra per evaluation near the truth); a vectorised or C
   version would cut the 2 to 3 times wall cost of the ray arms.
5. No separate "rule only at the screen" versus "rule at the re-cost and quick MC" split: the plan applies it to all of
   levels 0-2, so which of the three uses carries the effect is untested.
6. Arm A's terms are untuned guesses; the measured case-f rate at the truth (0.029 of events) suggests the window may
   still be too small for f1 or that f1 is simply too strict.
7. Case shares use the screen cost function; the local-cost case table is in `summary.json` only.
8. n_det > 2 chaining untested.
