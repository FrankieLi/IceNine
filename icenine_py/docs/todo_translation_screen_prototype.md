---
title: "TODO: prototype a table-driven level-0 screen (ROADMAP S4)"
subtitle: "Plan only: is a per-sample table of affine spot coefficients faster than BatchedScreen?"
date: "2026-10-13"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Planned, not started** (owner, 2026-10-13: plan the prototype and add it to the roadmap; execute
later when the roadmap reaches it). Track S, Ilium Agent. Work on a task branch
`feature/python-speed-parity-s4-translation-screen` (`scripts/dev/start_task.sh`); agent review
when done (`/review-task`). No library change in Steps 0-1.

# What

Build, once per sample, a table that holds for every level-0 orientation and every one of its
peaks the frame and the affine map from the voxel centroid to the projected pixel on each detector,
$\mathbf{u}_{pd} = A_{pd}\mathbf{x} + \mathbf{c}_{pd}$ (`docs/spot_translation_derivation.md`).
Screen each voxel from the table: one $2\times3$ matrix-vector product per peak and detector, the
ray-direction sign, then the existing stage-D pixel-radius overlap. Measure the time per voxel
against `BatchedScreen` and check that `peak_overlap` and the pass sets agree.

# Why it might help

- The level-0 orientation set does not depend on the voxel: `setup.fz_orientations`
  (4,886 FZ points) times `generate_local_grid_multi_level(local_grid_radius, 0, 1)` (9 points),
  43,974 orientations (`icenine/reconstructor.py:498,513`). Today every no-start search recomputes
  the Bragg solve, eta filter, frame lookup and ray geometry for all of them
  (`ScreenEvaluator.overlaps_exact`, `icenine/batched_screen.py`).
- For a fixed orientation the spot moves rigidly with the voxel (checked to $10^{-3}$ px,
  `scripts/speed_parity/translation_check.py`), so all of that is reusable across voxels.

**Expected ceiling (read this first).** With `BatchedScreen` (default since 2026-10-13) the level-0
screen costs about 24 µs per orientation (MIGRATION_HISTORY "Task A", quiet single-worker table:
24.1-24.2 µs), so about 1.1 s of a 76-96 s no-start search: **about 1%**. In the full-sample BFS only
seed voxels run a no-start search (350 to 504 seeds per run, MIGRATION_HISTORY seed-pass estimate),
so the whole-run saving is minutes. S4 is worth building only when S2 (fused C evaluation) and S3
(batched quick MC) have made the rest of the search much cheaper, which raises the screen's share.

# Framing

Split the screen into an orientation-only part (per sample) and a position-only part (per voxel).
Keep every operation up to the frame lookup identical to `ScreenEvaluator._counts`, so the peak sets
are the same by construction; only the ray geometry changes form (affine coefficients, computed in
float64), which can move a pixel coordinate across an integer boundary in rare cases.

# Experiments (sketch)

**Step 0: measure before building (go / no-go).** On Example2.ThreeVoxels (here) and one Phase D
full-sample voxel (Newt Agent's machine), for `BatchedScreen` at level 0:

- wall time per voxel, split into the torch geometry and the per-orientation Python-to-C stage-D
  loop (`batched_screen.py:262-277` makes one `_c_stage_d_overlap` call per orientation, about 44k
  calls per voxel);
- its share of a no-start search, today and projected with the S2/S3 estimates.

Go if the projected share after S2/S3 is large enough to matter (owner to set the bar, for example
10%); otherwise park S4 with the numbers.

**Step 1: prototype (script only, `scripts/speed_parity/translation_screen.py`).**

- `build_table`: the torch operations of `ScreenEvaluator._counts` up to the frame lookup, then per
  peak and detector, in float64: the parallel flag, the ray parameter $t(\mathbf{x}) =
  \mathbf{t}_{\text{lin}}\cdot\mathbf{x} + t_0$ (hit when $t > 0$), $A$ and $\mathbf{c}$ for the
  column and row coordinates; peaks sorted by orientation; the image list built once.
- `screen_voxel`: $\mathbf{u} = A\mathbf{x} + \mathbf{c}$, the hit test, the same `int32`
  truncation as today, then the same stage-D C call per orientation.
- Report: table build time and memory, time per voxel (geometry and stage-D loop separately) against
  `BatchedScreen` run as the reconstructor runs it (32 FZ points per call).

**Step 2 (only if Step 1 shows the per-orientation stage-D loop dominates).** One C call for all
orientations (a segmented `stage_d_overlap`, returning `peak_overlap` per orientation). That is a new
function in `icenine/_rasterize.c` (no change to existing paths) or a scratch extension for the
prototype; decide at that point.

**Equality checks (every step).** Per orientation `peak_overlap` and pass (`peak_overlap > 0`) for
the three ThreeVoxels voxels and random in-layer positions (several hundred µm), and for a few
full-sample voxels; list every orientation whose pass status differs, with the pixel distance of the
coordinate that flipped.

**Success criteria.** Pass sets identical apart from listed float-boundary cases; level-0 time per
voxel at least 5 times below `BatchedScreen` (excluding the build); build time no more than a few
voxels' worth of `BatchedScreen` time; table memory measured and acceptable for the full sample
(rough size: orientations x peaks per orientation x detectors x 12 numbers, of order 100-250 MB,
to be measured).

# Risks

- **Small payoff** (see the ceiling above) until S2 and S3 land.
- **Float differences:** float64 coefficients against today's float32 path flip a pixel at an integer
  boundary in rare cases; Task A saw the same class of difference between C++ and Python.
- **Memory** of the table for the full sample; float32 storage halves it at the cost of more boundary
  flips.
- **Scope:** `CppEdgeRules` off only (Track E's task A-exact changes the edge rules; extend after it
  lands); `CppOmegaInterp` off (an interpolated-omega table is possible, also orientation-only, not in
  scope); no ray rule (`BatchedScreen` does not support it either).

# Relation to other work

- ROADMAP Track S: S2 and S3 decide whether S4 is worth it; S5 (BFS neighbour's first evaluation by
  translation) uses the same coefficients; S6 (full-sample timing).
- `docs/spot_translation_derivation.md` (the affine map and its check).
- MIGRATION_HISTORY "Task A" (`BatchedScreen` timing and equality), "Python speed parity" (the
  per-evaluation split).
- Track E: unaffected (S4 changes how level 0 is computed, not which candidates it returns).
