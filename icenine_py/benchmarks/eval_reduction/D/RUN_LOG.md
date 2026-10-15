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
| pairing, 12 presets | `pairing_run.py --presets all` | about 2 s each | single process, not preflight-gated; no timing claim |
| full suite | `uv run pytest tests/ --junitxml=...` | 331 s | 1027 passed, 37 skipped, 1 deselected (the 10 new tests are `tests/test_spot_pairing.py`) |
| attribution | `pairing_attribution.py attach` (330 voxels), `separate` (30 voxels, 91 orientations each) | 1 to 2 min, 10 s | single process |

The check compares every recorded image with the existing `clean/` file: 360 of 360 frames equal; 57,389 components equal `spots_total`
of `render_full.json`. 2,957,243 triangle records.

## 3. What went wrong, and decisions taken (non-blocking, plan default)

* Single-source spots do not exist in this sample (10 and 20 of 57,389), so the plan's decision number is undefined. The rule was applied to the nearest
  class (one grain, several voxels): recall 0.990, so **Step 4 not run**.
* The first size gate (log ratio with a +1 pixel offset) rejected 7% of true pairs by one-pixel quantisation of thin spots; changed to a +4 pixel offset
  before the numbers were collected (the earlier numbers are not quoted).
* The line through two centroids is a weak constraint alone (80% of spots have several candidates); the ring angle of the crystal was added as a gate. This
  uses the structure (eight Cu rings at Q-max 8), not the orientation.
* The attribution by the origin cannot work for union spots, so a footprint rule B was added next to the plan's radius rule A.
* The 30 voxels of Step 3 are the 30-voxel ray-study set (it contains the six); the plan's "six and 30 others" is 30, not 36.
* Direction tolerance of the attributed cost: default 0.05° set from the pixel tolerance, then 0.02, 0.1 and 0.2° were scanned and are reported.

* The findings document `docs/findings_spot_pairing_2026-10.md` the plan lists could not be written by the implementer (the write tool refused every file with
  that name); its content is the MIGRATION_HISTORY subsection "Task D" and the owner report. The main session needs to create the file.

## 4. Could be done better

* Peel unions of several grains (split by the detector-to-detector translation) before pairing; the 4,386 several-grain pairs are 15.7% of the pairs.
* A per-pair attribution radius from the predicted origin tolerance.
* Trace the 3,890 detector-0 spots with no detector-1 partner.
* Time the pairing under the preflight gate if it ever matters (about 2 s for 180 frames).
