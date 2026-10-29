# RUN_LOG: Task F carry-option comparison (D-20261021-N1, 2026-10-21)

Task branch `feature/eval-reduction-Fcompare` (parent `feature/eval-reduction`). Request: `docs/DECISIONS.md` D-20261021-N1.
All Python through `uv run` from `icenine_py/`. Clean images only, no extra seeds. `BatchedQuickMC`, `SearchDevice` and
`CandidateDedup` stay off (no reconstructor run). Raw arrays are in `benchmarks/eval_reduction/Fcompare/cache/` (gitignored);
everything quoted comes from `summary_*.json`, `tables_*.md` and `timing/`.

## 1. What changed

| File | Change |
|---|---|
| `scripts/eval_reduction/outward_compare.py` (new) | `sets` (voxel lists), `run` (F4, F4 5-frame, F5, F5+F4, F5+F4 5-frame, CMA-ES 1 and 2 deg, from the same starts), `timing` (single worker, preflight-gated) |
| `scripts/eval_reduction/outward_compare_summary.py` (new) | per-set tables, paired counts, cheapest-per-level table, contention table, the figure, `merge` |
| `scripts/eval_reduction/outward_bench.py` | `run_f54` takes a `window` argument (default 1: no change to Task F's runs) |
| `.gitignore` | `icenine_py/benchmarks/eval_reduction/Fcompare/cache/` |
| docs | findings doc section, MIGRATION_HISTORY subsection, README script row, ROADMAP row F |

Nothing under `icenine_py/icenine/` was touched.

## 2. Decisions before the data

* **Near-axis group.** Only 10 of the 200 E0 voxels have r_perp below 120 um (`voxel_sets.json`, `notes`) and Task F used all 10
  (8 benchmark, 2 pilot). The 8 "near" voxels are the unused ones of smallest r_perp (122 to 150 um). The 16 far voxels
  (r_perp >= 250 um) are spread evenly (`np.linspace`) over the 138 unused far voxels. Full sample: the 24 voxels of
  `phase_d_seed_diag/ray_overlap/voxels.json` that are not `SIX_VOXELS`.
* **Policy 5.** The 5-frame policy uses the 5-frame window in both F4 parts (F5 then F4 with the window is `f54w`, the
  fall-back of that policy). Task F's `policy_f4w` used the 1-frame F4 after F5.
* **Budget.** Pilot (`--take 3`, 1 voxel per set, 4 offsets, all methods): E0 voxel 2910 17 s, full voxel 3156 36 s, so the whole run
  was far under the 3 h cut-off and all 24 full-sample voxels were kept.
* **Weights.** Task F's current strict-equivalent weights (`F/timing/{e0,full}.json`, `{e0,full}_f5.json`), not new ones.

## 3. Runs (from `icenine_py/`)

| What | Command | Wall | Label |
|---|---|---|---|
| sets | `outward_compare.py sets` | seconds | |
| pilot | `outward_compare.py run --kind e0 --workers 1 --take 3 --voxels 2910 --tag pilot_e0`, same for full voxel 3156 | 19 s, 36 s | two processes at once |
| main, both sets | `benchmarks/eval_reduction/Fcompare/run_all.sh` (E0 `--workers 2`, full `--workers 4`) | E0 1277 s, full 898 s | **contended** (6 processes) |
| timing | `benchmarks/eval_reduction/Fcompare/run_timing.sh` (`outward_compare.py timing --kind e0|full`) | about 1 min each | **single-worker**, preflight-gated |
| summaries | `outward_compare_summary.py e0`, `full`, `merge`, `figure` | seconds | |

Timing preflight (`timing/preflight_e0.json`, `preflight_full.json`): a first manual attempt, right after the main runs, was
refused for both sets (1-minute load 1.53 against the limit 1.5, one busy process). Then `run_timing.sh` retried: E0 passed on
its first attempt (load 1.35), the full-sample pass was refused twice (one refusal seen: load 1.83, left over from the E0 pass)
and passed on its third (load 1.48, no busy processes, 1 torch thread). Both passes ran alone; nothing else was started while
they ran.

## 4. Checks

* The timing pass reproduces the main pass: for every option and set, the 40 runs (2 voxels x 5 directions x 4 offsets) have
  the same end-point error (1e-6 deg) and the same strict-equivalent count as the same runs in the main pass (`summary_*.json`,
  `timing_reproduces_main`).
* The strict-equivalent counts of F4, F4 5-frame, F5 and F5+F4 equal Task F's `outward_summary.strict_equiv` (asserted in
  the summary script).
* Doc numbers audited with `scripts/dev/audit_numbers.py` against the tables and summaries; unmatched: the decision id,
  the CMA-ES budget (4,000), the 200 and 250 definitions of the voxel groups, the seed 7000. The contended-to-single-worker
  ratio range in the findings doc was first read off the wrong column (medians of 480 runs against 10) and corrected to the
  same-run table (1.07 to 1.29).

## 5. Not done

* A second single-worker timing pass or one on more voxels: the seconds rest on 10 runs per option, offset and set (two voxels per set).
* Realistic images (out of scope).
