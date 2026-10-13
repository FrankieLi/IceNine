# RUN_LOG: centroid start point everywhere, step 3 (clean images only), 2026-10-12 to 2026-10-13

Task branch `feature/finisher-mc-study-centroid` (parent `feature/finisher-mc-study`). Step 1 (commit d96c413) made the
voxel centroid the start point of the pixel-radius screen for every overlap rule; vertex 0 is only the opt-out
`cpp_center_screen=False`. Step 3 re-runs the clean-image experiments of the ray-overlap task with that build.
Owner decisions: clean images only (the realistic images are skipped for now), one seed (no extra seeds).
All wall times are "contended" (10 workers sharing the laptop; the six full-sample runs were six processes at once).
`preflight.json` was taken before the run started (load 1.2, no busy processes); the machine was not guarded
during the run, so nothing here is a single-worker timing claim.

## Before the re-run

The first attempt (stopped half-way) was deleted: `ray_overlap/coarse30_centroid/` (partial), and its `peak_centroid/`
and `centroid/` folders were removed too and made again, because the run is cheap and a clean restart was safer than
auditing a half-finished one. The two script edits it left (`e0_cpprule.py`: `--tag`, `--variant`;
`ray_overlap_study.py`: `RAY_VARIANTS`, `RAY_PEAK_DIR`, `RAY_VERTEX0`) were reviewed and kept; one line was added so the
coarse stage also follows `RAY_VERTEX0`.
Pilot (1 voxel, v22867, clean, one process, 8 torch threads): old arm 171 s, C++ arm 123 s for levels 0-2.

## Runs (one driver, `run_all.sh`, from `icenine_py/`)

| Part | What | Command | Wall |
|---|---|---|---|
| a | peak level, 30 voxels, clean, arms old, cpp, f0_R1, R2, R4, R8, A | `RAY_PEAK_DIR=peak_centroid ray_overlap_study.py run --step 1` | about 10 s |
| b | coarse stage levels 0-2, 30 voxels x 4 arms (old, cpp, f0_R2, f0_A), clean | `RAY_COARSE_DIR=coarse30_centroid ray_overlap_study.py run --step 2 --arms old,cpp,f0_R2,f0_A` | about 36 min, contended |
| c | E0, 200 voxels, clean, min_sin_eta 0, arms off (= old), on (= cpp), rayR2, rayA | `e0_cpprule.py run --arm X --msin 0 --tag centroid --variant clean` | about 9 min each (36 min in all), contended |
| d | six full-sample voxels, old rule | `seed_diag.py run --voxel V --images full --tag base_centroid` (six processes at once) | about 3 min, contended |

Total about 1 h 15 min (23:19 to 00:35, clock times of 2026-10-12 and 10-13). Logs in `cache/logs/`.

Why all four arms were re-run in parts b and c, not only the old arm: the plan allowed re-using the cpp, R2 and A
runs because their path already used the centroid. Re-running them cost little and gives a check: their coarse-stage
and E0 results come out the same as the stored runs (identical coarse-stage tables and medians, identical E0 mean evaluations and wrong counts, 0 voxels fixed or broken; individual final orientations were not compared), so
nothing else moved in the build between the stored runs and this one. Part d re-ran only the old rule (stored
comparison: `runs/lenientnow_full_*`, the lenient rule in the build of 2026-10-12).

## Where the results are

`summary.json` and `tables.md` come from `scripts/phase_d/centroid_step3_summary.py`. "Before" = stored runs with the
old rule's screen at vertex 0: `ray_overlap/peak`, `ray_overlap/coarse30`, `scripts/phase_d/cache/off_ms0` (and
`on_ms0`, `rayR2_ms0`, `rayA_ms0`), `runs/lenientnow_full_*`. "After" = `ray_overlap/peak_centroid`,
`ray_overlap/coarse30_centroid`, `scripts/phase_d/cache/*_ms0_centroid`, `runs/base_centroid_full_*`.
Reading: `icenine_py/docs/findings_ray_overlap_rule_2026-10.md` section 10 and MIGRATION_HISTORY "Centroid start point
everywhere".

## Could be done better

- One seed per voxel. The old rule's 30 changed voxels in E0 (14 fixed, 16 broken) show how much a small change of the
  candidate list moves a single run; they cannot say which of the two start points is better.
- The realistic-image numbers for the old rule in the earlier write-up are still the vertex-0 ones (not re-run, owner).
- Wall times are contended and for orientation only.
