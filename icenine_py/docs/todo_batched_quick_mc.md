---
title: "TODO: batched quick MC with a GPU option (ROADMAP S3)"
subtitle: "Plan only: run every candidate's quick MC of a level in lockstep through one batched local cost"
date: "2026-10-13"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Planned, not started** (owner, 2026-10-13: a GPU path is wanted for the search cost, so S3 goes
before S2; write the S3 plan now). Track S, Ilium Agent. Order in Track S: S1, then S3. Before
starting: merge `feature/eval-reduction` (A-exact changes `cost_functions.py` and `_rasterize.c`
under `CppEdgeRules`), read Newt Agent's mailbox, and push a `heads-up` (S3 touches the shared files
`reconstructor.py` and `orientation_search.py`). Task branch `feature/python-speed-parity-s3`
(`scripts/dev/start_task.sh s3`). An `implementer` agent builds it; `/review-task` (Opus) when done,
required fixes to a fresh implementer (Track S rule). `/handoff`, `/clear`, `/resume` afterwards
(protocol rule 5).

**Owner decision needed at approval:** the random-draw assignment (A or B below).

# What

At each level of `AdaptiveVoxelReconstructor.reconstruct_voxel`, the quick MC runs one candidate at
a time (`icenine/reconstructor.py:622-633`), each through `MCOptimizer.optimize`
(`icenine/orientation_search.py:606`), each evaluation one `VoxelCostFunction.evaluate` call of
about 320 µs, of which about 10 µs is physics and pixel work (MIGRATION_HISTORY "Python speed
parity"). S3 runs all candidates of the level in lockstep: one batched local-cost call scores the
current trial of every candidate, so the fixed per-call overhead is paid once per step, not once per
candidate. A device option runs the torch stages on a GPU (MPS on Ilium Agent's machine, CUDA
elsewhere).

# Facts from the code (checked on 2c498b7)

- **One block per quick MC.** The quick MC calls `optimize(max_mc_steps=10, max_restarts=5)` with
  `angular_step = trial_box_width * MCRadiusScaleFactor`. All 12 configs set `MCRadiusScaleFactor
  0.4`, so `n_min_ergodic = int(2 (1/0.4)^3) = 31` (`orientation_search.py:664`) and the run is one
  block of 10 trial steps. The step, and so the trial radius `tan(step)/sqrt(12)`
  (`orientation_search.py:590`), is the same for every candidate of a level.
- **Evaluations per candidate: 12.** The start orientation as a matrix (in `optimize`), the start
  again after a quaternion round trip (first line of `_mc_block`), and 10 trials. The two start
  evaluations can differ in the last bit of the matrix, so the batched path keeps both.
- **Random draws per candidate: 30 or 33, not always 30.** Each trial draws 3 uniforms. When no trial
  improves (`tmp_cost >= global_min_cost`, `orientation_search.py:678`), the failure branch draws 3
  more for a restart that is never used (the step budget is spent). The ROADMAP S3 row said "always
  draws 30", so pre-drawing 30 per candidate keeps today's trial rotations only for the candidates
  before the first failed one. See "Random draws" below.
- **The local cost.** `_make_local_cost_fn` (`reconstructor.py:759`): `mode="hard"`, `pixel_radius
  0`, plus `CppOverlapRule`, `CppEdgeRules` and the ray rule when set. The cost is `1 - quality`; the
  quality is a running mean over peaks accumulated in the C stage D (`_rasterize.c`
  `stage_d_overlap`, `cost_functions.py:771-800`). The torch stages A-C already run in float32
  (`evaluate` casts the orientation with `.float()`).
- **What exists.** `ScreenEvaluator` (`icenine/batched_screen.py`) batches the pixel-radius screen
  (about 24 µs per orientation, MIGRATION_HISTORY "Task A") but returns only the screen pass, not
  the local cost. No batched local cost exists.

# Random draws (owner decision)

- **A. Exact (today's draws).** Run the lockstep assuming no earlier candidate fails; when candidate
  j's block fails, every later candidate used draws 3 positions too early, so re-run those
  candidates in another lockstep round with the corrected offset. Trial rotations equal today's for
  every candidate. Cost: one extra round per failed block, each re-running all candidates after it;
  if many blocks fail (likely at the finer levels), it falls back towards the sequential cost.
- **B. Block assignment (recommended).** Draw the 30 trial uniforms of all N candidates first, in
  candidate order (`rng.uniform(-r, r, size=(N, 10, 3))`), then 3 per failed candidate in candidate
  order. The generator consumes the same number of draws as today, so its state after the level, and
  every later level and FindOptimal draw, is unchanged. The candidates before the first failure keep
  today's trial rotations; the others get different but equally distributed ones, so per-voxel
  results can differ, and the check is by wrong counts. Tested exactly by injecting today's
  per-candidate draws (Step 2).

Under either option the batched path is opt-in, so goldens and default results do not change.

# Steps

**Step 0: go/no-go measurement (script only, `scripts/speed_parity/`, no library change).**
On Example2.ThreeVoxels voxel 0 and 10 voxels of the E0 clean set:

- quick-MC candidates per level (the batch sizes) and the fraction of failed blocks per level;
- a throwaway batched forward of stages A-C for B = 16, 64, 256, 1,024 orientations on CPU and on
  MPS: µs per orientation, and the number of host syncs per call (masking with data-dependent shapes
  forces one per op on a GPU).

Go if the batched cost at the measured batch sizes is at least 3 times cheaper per orientation than
`evaluate` on CPU. Report the MPS figures either way; a GPU that loses at these batch sizes is a
finding, not a failure. Estimate: if the quick MC (about 58% of evaluations with the C++ rule) falls
from about 320 µs to 30-40 µs per evaluation, a no-start search gets up to about 2 times faster;
Step 0 replaces this estimate with measured shares.

**Step 1: `BatchedCostEvaluator` (new module `icenine/batched_cost.py`).**
Built from one `VoxelCostFunction`; `evaluate_many(R (B,3,3), voxel_vertices, phase_index) ->
List[OverlapInfo]`, adds B to `cost_fn.eval_count`. Stages A-C batched over orientations with a fixed
shape (all 2K omega rows per orientation, validity as a mask, so the GPU path has no data-dependent
shapes until the hand-off to C). Stage D stays in C on CPU: one `stage_d_overlap` call per
orientation at first (about 8 µs each, quality accumulated exactly as today); a single batched C
entry only if Step 4 shows stage D matters. Covers the default cost, `CppOverlapRule` and
`CppEdgeRules`; the ray rule falls back to `evaluate`. Reuses S1's lookup and image tables if S1 has
landed. `device` argument: `cpu` (default), `mps`, `cuda`.

Tests (`tests/test_batched_cost.py`): on a fixed seeded set of orientations near and far from the
truth, every `OverlapInfo` field equal to `evaluate`'s on CPU. If some differ by a float32 last bit
in another batch shape (as for `BatchedScreen`), count them, name the stage, and stop for the owner
rather than loosening the test. On MPS (skipped when unavailable): count the orientations whose
counts differ.

**Step 2: lockstep quick MC (new function in a new module `icenine/batched_quick_mc.py`).**
Inputs: the level's candidates, `trial_box_width`, `trial_step`, the cost evaluator, the generator.
Per step: build the N trial rotations from the pre-drawn uniforms (`get_near_identity_point`,
`_quat_multiply`, `quaternion_to_matrix` as today; vectorise only where a test shows bit-identity,
else loop per candidate, which is cheap), one `evaluate_many`, strict `<` acceptance per candidate.
Fills `orientation`, `cost`, `overlap_info` and a per-candidate `last_run` as `optimize` does.
Applies only when the run is a single block (`n_min_ergodic >= max_mc_steps`) and `level_mc` is a
plain `MCOptimizer` (`type(level_mc) is MCOptimizer`; study subclasses such as `LoggedMC` in
`scripts/finisher_diagnosis/diagnose.py` keep the sequential path); otherwise sequential, with a
warning only when the key was set explicitly (the `BatchedScreen` fallback rule).

Tests (`tests/test_batched_quick_mc.py`):

- With an injected draw sequence (a generator stub that serves each candidate today's draws), the
  lockstep result equals sequential `optimize` for every candidate (orientation, cost, counts,
  `last_run`), first with `evaluate` as the cost (logic only), then with `BatchedCostEvaluator` on CPU.
- Pre-drawn `rng.uniform(-r, r, size=(N, 10, 3))` equals N x 30 scalar draws.
- Option B: generator state after the level equals the sequential run's.

**Step 3: wiring.** Config keys `BatchedQuickMC` (default 0) and `SearchDevice` (default `cpu`) in
`config_file.py` and `SearchParameters`; the level loop in `reconstructor.py` calls the lockstep
function when the key is on. Pin `BatchedQuickMC` off where evaluations are counted by patching
`evaluate` (`scripts/eval_reduction/eval_ledger.py` and the other scripts that pin `batched_screen`,
as eb1a98d did), or count `evaluate_many` there; evaluation counts must stay the same with the key on.

**Step 4: measurements.** Quiet machine, single worker, key off vs on, CPU vs MPS:

- ThreeVoxels voxel 0 no-start search: wall time, quick-MC time per level, µs per evaluation
  (`scripts/speed_parity/`);
- E0 clean set, 200 voxels, seed 0, key off vs on (CPU): wrong counts and evaluations side by side,
  per-level candidate counts, and the number of voxels whose answer changed. Three seeds only if seed
  0 differs by more than a few voxels.

**Step 5: docs.** MIGRATION_HISTORY section (results, the draw option, fallbacks), README (new
modules, keys), ROADMAP S3 row.

# Rules (ROADMAP "How to use this file")

New keys opt-in and off; no default flips and no golden regeneration (list them for the owner);
evaluations and wrong counts side by side; no p-values; clean images only (Q8 open). Merge into
`feature/python-speed-parity` with `--no-ff` after review; push and post mailbox entries.

# Out of scope

The batched screen on a GPU (a natural follow-up: level 0 has 43,974 orientations per voxel, a
better fit for a GPU than the quick MC's batches); FindOptimal and VarianceMinimizing (sequential by
nature); multi-block quick MC (no current config uses it); S2.
