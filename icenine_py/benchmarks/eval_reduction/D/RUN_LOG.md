# RUN_LOG: Task D, spot pairing across detectors and peak-to-voxel attribution (2026-10-14)

Task branch `feature/eval-reduction-D` (parent `feature/eval-reduction`). Plan: `~/.claude/plans/eval-reduction.md`, Task D.
All Python through `uv run` from `icenine_py/`. Clean full-sample images only (sample seed 0, Q-max 8), no realistic images, no
extra seeds, `BatchedQuickMC`, `SearchDevice` and `CandidateDedup` at their defaults. Raw arrays are in `cache/` (gitignored);
every quoted number comes from `step1/summary.json`, `step2/summary.json`, `step3/attach_summary.json`,
`step3/separate_summary.json` and the `tables*.md` beside them.

## 1. What changed

| File | Change | Why |
|---|---|---|
| `icenine/forward_simulation.py` | opt-in `ForwardSimulation.peak_hook` (default `None`; called once per rasterized peak with voxel, reflection, image index, omega) | names the source of each triangle; no physics change, no output change when off |
| `icenine/spot_pairing.py` (new) | detector geometry, `line_origin`, `SampleHull`, `PairingParams`, `pair_frame`, `assign_pairs`, `attribute_peaks` | the pairing |
| `tests/test_spot_pairing.py` (new) | synthetic two-detector patterns (round trip, origin recovery, gates, assignment, attribution) and the hook (off by default; hook calls match rasterized triangles) | correctness |
| `scripts/phase_d/render_provenance.py` (new) | sibling of `render_full.py`; recording image; writes provenance.npz | step 1; `render_full.py` untouched |
| `scripts/eval_reduction/pairing_*.py` (new) | common loaders, run, summary, attribution, separation summary, figures | steps 2 and 3 |
| `.gitignore` | `icenine_py/benchmarks/eval_reduction/D/cache/` | raw arrays |

Not touched: `cost_functions.py`, `orientation_search.py`, `batched_screen.py`, `reconstructor.py`, `_rasterize.c`, any physics code, any golden test.

## 2. Runs (commands from `icenine_py/`)

| What | Command | Wall | Label |
|---|---|---|---|
| pilot, 500 centre voxels | `render_provenance.py --workers 10 --n-voxels 500 --tag pilot --check` | 4.5 s | 10 workers, machine quiet at start |
| full render with provenance | `render_provenance.py --workers 10 --tag full --check` | 35 s | contended (10 workers; load 2.2 at start, preflight warned) |
| pairing, 10 presets + 2 binary-centroid runs | `pairing_run.py --presets all`; `--presets ring_tight,shape_tight --centroid binary` | about 2 s each | single process, not preflight-gated; no timing claim |
| full suite | `uv run pytest tests/ --junitxml=...` | 331 s | 1027 passed, 37 skipped, 1 deselected (the 10 new tests are `tests/test_spot_pairing.py`) |
| full suite after the first review's fixes (a8c9b5c) | `uv run pytest tests/ --junitxml=...` | 338 s | 1030 passed, 37 skipped, 1 deselected (13 tests in `tests/test_spot_pairing.py`, 3 added by the fixes) |
| full suite after the second review's fixes (d967595) | `uv run pytest tests/ --junitxml=...` | 338 s | 1030 passed, 37 skipped, 1 deselected (unchanged) |
| attribution | `pairing_attribution.py attach` (330 voxels), `separate` (30 voxels, 91 orientations each) | 1 to 2 min, 10 s | single process |

The check compares every recorded image with the existing `clean/` file: all 360 of 360 images have the same lit pixels, with intensities equal to the stored precision; 57,389 components equal `spots_total`
of `render_full.json`. 2,957,243 triangle records.

## 3. What went wrong, and decisions taken (non-blocking, plan default)

* Single-source spots do not exist in this sample (10 and 20 of 57,389), so the plan's decision number is undefined. The rule was applied to the nearest
  class (one grain, several voxels): recall 0.987 (0.990 before the review fixes), so **Step 4 not run**.
* The first size gate (log ratio with a +1 pixel offset) rejected 7% of true pairs by one-pixel quantisation of thin spots; changed to a +4 pixel offset
  before the numbers were collected (the earlier numbers are not quoted).
* The line through two centroids is a weak constraint alone (76% of spots have several candidates after the review fixes, 80% before); the ring angle of the crystal was added as a gate. This
  uses the structure (eight Cu rings at Q-max 8), not the orientation.
* The attribution by the origin cannot work for union spots, so a footprint rule B was added next to the plan's radius rule A.
* The 30 voxels of Step 3 are the 30-voxel ray-study set (it contains the six); the plan's "six and 30 others" is 30, not 36.
* Direction tolerance of the attributed cost: default 0.05° set from the pixel tolerance, then 0.02, 0.1 and 0.2° were scanned and are reported.

* The findings document `docs/findings_spot_pairing_2026-10.md` the plan lists could not be written by the implementer (the write tool refused every file with
  that name); its content is the MIGRATION_HISTORY subsection "Task D" and the owner report. The main session needs to create the file.

## 3b. Review fixes (D-20261014-N9, 2026-10-14)

All 16 required fixes of `.claude/reports/20261014-2208-review-D.md` applied. Two change numbers:

* Fix 12: `DetectorGeometry.distance_x()` returned the lab x of the pixel-(0, 0) corner (3.263 / 5.290 mm), not the detector distance
  (3.362 / 5.387 mm); L1/L0 1.621 -> 1.603, voxel term 5.82 -> 5.65 um (`plan` presets only).
* Fix 13: the detector-1 centroid tolerance was 1 px x L1/L0, counting the detector-0 error twice (the line fit already propagates it); now 1 px
  of its own grid. This tightens every gate that uses the tolerance (hull, ring).

Rerun (commands from `icenine_py/scripts/eval_reduction/`, single process, logs in `cache/logs/fix_*.log`): `pairing_run.py --presets all`,
`pairing_run.py --presets ring_tight,shape_tight --centroid binary`, `pairing_summary.py`, `pairing_attribution.py attach`,
`pairing_attribution.py separate` (30 voxels, 12 s), `pairing_separation_summary.py`, `pairing_figures.py`, then `collect_tables.py`.
`ring_tight`: precision 0.854 -> 0.887, link recall 0.889 -> 0.875, one-grain recall 0.990 -> 0.987, several-grain recall 0.349 -> 0.275;
origin-error median 2.667 -> 2.659 um, p90 15.66 -> 15.57 um (23,252 -> 23,205 pairs). Step 1 is unchanged. The step-3 truth, sharpness and
separation tables are unchanged at their printed precision (completeness at 0.05 degrees is identical per record; purity and the 0.1 and 0.2 degree
columns move). `A r=15.7 um` keeps its radius (the old p90). Fix 16 rebuilt `cache/geometry_s0.npz` with the frame check (arrays identical).

## 3c. Second review fixes (D-20261014-N10, 2026-10-14)

All 10 required fixes of `.claude/reports/20261014-2243-review-D-2.md` applied; they are doc fixes. One run: `pairing_summary.py`
(single process, about 5 s, log `cache/logs/fix2_summary.log`). It now also writes `step1/summary.json` `partner_size_ratio` and, in
`step2/summary.json` `derived_geometry`, the measured radial ray-direction error of the correct one-grain `ring_tight` pairs and the
both-detector 1-pixel range. Every other key and all tables are unchanged; `doc_tables_all.md` regenerates byte-identical.

## 4. Could be done better

* Peel unions of several grains (split by the detector-to-detector translation) before pairing; the 4,386 several-grain pairs are 15.7% of the true pairs.
* A per-pair attribution radius from the predicted origin tolerance.
* Trace the 3,890 detector-0 spots with no detector-1 partner.
* Time the pairing under the preflight gate if it ever matters (about 2 s for 180 frames).
