# Q8 run log (2026-10-14, branch `feature/python-speed-parity-q8-nodrop`)

Commands from `icenine_py/`. Every wall time is contended (the machine never met `require_quiet`:
the desktop apps keep the 1-minute load at 2 to 3.5, and the Step 4 blocks follow one another,
each a 10-worker pool). Not timing claims.

| Step | Command | Result |
|---|---|---|
| 2 render | `uv run python scripts/phase_d/render_full.py --workers 6 --noise nodrop --variants realistic` | `Examples/Example2.ManyGrains/ScatteringData_PhaseD/full/realistic_nodrop/` (gitignored), `scripts/phase_d/results/render_full_nodrop.json` (49 s) |
| 2 sanity | `uv run python scripts/phase_d/sanity_checks.py --variant realistic_nodrop` | `scripts/phase_d/results/sanity_full_realistic_nodrop.json` (14 s). The clean and realistic variants were re-run with `--skip-containment` to check for code drift (cost numbers identical, only `load_s` differs) and the tracked files restored with `git checkout` |
| 2 truth, Phase D | `uv run python scripts/phase_d/q8_sanity.py run --workers 6`, then `summary` | `sanity_truth.json`, `sanity_summary.json`, `sanity_tables.md` (builds the uint8 stacks of the three variants under `ScatteringData_PhaseD/stacks/`, gitignored) |
| 4 images | `uv run python scripts/findoptimal_robustness/e0_run.py images --variants all realistic_nodrop --n-voxels 200 --workers 10` | 400 case images in `scripts/findoptimal_robustness/cache/images/` (all buildable; about 1 minute) |
| 4 pilot | `uv run python scripts/phase_d/e0_cpprule.py pilot --arm off --workers 10 --limit 10 --variant realistic_nodrop --tag q8 --meta-dir benchmarks/q8_nodrop/pilot` | 10 voxels in 31 s wall, 9 right (`pilot/`) |
| 4 runs | `nohup bash scripts/phase_d/run_q8_step4.sh` (log `scripts/phase_d/cache/logs/q8_step4.log`, marker `q8_step4.done`) | 18 blocks (6 arms x 3 variants), 8.1 to 9.4 minutes each, 157 minutes in all; `runs/run_meta_*.json`, `runs/preflight_*.json`; per-voxel results in `scripts/phase_d/cache/<arm>[_ms0]_q8/` (gitignored) and, compacted, in `e0_per_voxel.npz` |
| 4 truth, E0 | `uv run python scripts/phase_d/q8_e0_truth.py run --workers 10`, then `summary` | `e0_truth_summary.json`, `e0_truth_tables.md` (18 s) |
| 4 summary | `uv run python scripts/phase_d/q8_e0_summary.py` | `e0_summary.json`, `e0_per_voxel.npz`, `e0_tables.md` |
| tables | `cat sanity_tables.md e0_truth_tables.md e0_tables.md > tables.md`, then `scripts/dev/sync_doc_tables.py` (against a scratch file holding the doc's own blocks plus these, because `MIGRATION_HISTORY.md` holds many other table blocks) | `tables.md` |

Not run: Step 3 (Q-max 16 render; the two E0 rows do not use the Phase D images), the Phase D
pilot and reseed on no-drop images (waits with X1).

The preflight files in `runs/` between arms show the load of the previous 10-worker pool (the
state before the next arm), not an idle machine; `preflight_start.json` is the state at the first
arm (load about 1.8 to 2.3, no busy process).

## Steps 6-8 (2026-10-14, after the owner decisions D-20261014-I5 to I8)

Same conventions: commands from `icenine_py/`, every wall time contended (10-worker pools on a desktop that never met
`require_quiet`), none a timing claim.

| Step | Command | Result |
|---|---|---|
| 6 tests | `uv run pytest tests/ --junitxml=<scratch>/junit.xml` after the default flip | 981 passed, 42 skipped, 1 deselected (251 s); the 42 skips are the C++ ground-truth and benchmark skips of before (977 + 4 new tests) |
| 7 pilot | `uv run python scripts/phase_d/q8_step7.py pilot --workers 2` (voxel 7336, then 3 voxels for the width) | level 0: 0.8 s per task for the lenient and C++ rule (batched screen), 12 to 13 s for the ray arm (per-orientation screen); width 0.6 to 0.9 s per (voxel, images); estimated 100 CPU-minutes, about 10 to 15 minutes on 10 workers |
| 7 width | `uv run python scripts/phase_d/q8_step7.py width --workers 10` | 400 tasks (200 voxels x 2 image sets; costs at the truth and at 7 radii x 50 axes, 3 arms, level-0 and final cost and the level-0 screen), 71 s wall |
| 7 level 0 | `uv run python scripts/phase_d/q8_step7.py level0 --workers 10` | 1,200 tasks (200 voxels x 2 image sets x 3 arms), 1,075 s wall; `scripts/phase_d/cache/q8_step7/` (gitignored) |
| 7 summary | `uv run python scripts/phase_d/q8_step7.py summary` | `benchmarks/q8_nodrop/step7/{summary.json,tables.md,per_voxel.npz}` (12 s, bootstrap included) |

Step 7 cross-check: the level-0 records of the Step 4 runs (`n_l0`, `n_within`, `min_dist_l0`) equal this script's in all 1,200
(voxel, images, arm) tasks, although Step 7 used the library's default batched exact screen for the lenient and C++ arms (the
Step 4 runs forced the per-orientation screen); `summary.json` key `step4_check`.
