---
title: "TODO: batched quick MC with a GPU option (ROADMAP S3)"
subtitle: "Plan only: run every candidate's quick MC of a level in lockstep through one batched local cost"
date: "2026-10-13"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Approved 2026-10-13, draw option A (exact); the next Track S task** (owner: a GPU path is wanted,
so S3 goes before S2; S3 goes before S1 because it does not wait for Newt Agent). Track S, Ilium
Agent. Task branch `feature/python-speed-parity-s3` (`scripts/dev/start_task.sh s3`). An
`implementer` agent builds it; `/review-task` (Opus) when done, required fixes to a fresh
implementer (Track S rule). `/handoff`, `/clear`, `/resume` afterwards (protocol rule 5).

**Why S3 can start before A-exact lands.** A-exact changes `cost_functions.py` and `_rasterize.c`
only under the opt-in `CppEdgeRules` (default path bit-identical). S3 adds new modules and does not
edit `cost_functions.py`; its batched cost calls the existing C `stage_d_overlap`, so A-exact's C
change reaches it when merged. Until then the batched path covers the default cost and
`CppOverlapRule`, and falls back to `evaluate` for `CppEdgeRules` (and the ray rule).

**Before starting:** `git fetch`, merge `feature/eval-reduction` if it moved, read Newt Agent's
mailbox, and push the `heads-up` posted with this plan (Step 3 edits the shared files
`reconstructor.py` and `orientation_search.py`). Run Step 0 and Step 4 timings on a quiet machine
(not during a Phase D render or a Q8 run).

**After A-exact merges into `feature/python-speed-parity`:** add `CppEdgeRules` to the batched cost,
extend the Step 1 equality test to it, and drop that fallback (a small follow-up commit, reviewed).

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

# Random draws: option A, exact (owner 2026-10-13)

Today candidate i's 30 trial uniforms start at stream position 30 i + 3 F_i, where F_i is the number
of failed blocks among candidates 0 .. i-1, so a candidate's draws depend on the outcomes of all
candidates before it. Option B (pre-draw 30 per candidate in a block) was not chosen. Option A keeps
today's draws for every candidate by speculation:

- **Buffer.** Save the generator state, draw a buffer of 33 N standard uniforms, map them to the trial
  radius `tan(step)/sqrt(12)` (and the unused restart radius) as `Generator.uniform` does, then
  restore the state and consume exactly 30 N + 3 F draws at the end, so the generator leaves the level
  in today's state (later levels and FindOptimal unchanged).
- **Rounds.** Each unresolved candidate runs its 10 trials in lockstep at an assumed offset: the
  failures known so far, plus 0 .. k-1 extra (speculation depth k, so k chains per candidate per
  round). After the round, walk the candidates in order: a candidate is resolved when its true
  offset, now known from the resolved candidates before it, is one of the offsets it ran; its chain
  at that offset is the result. Unresolved candidates go to the next round. Rounds end when all are
  resolved; with k = 1 that is one round per failed block, plus one.
- **Start evaluations** (2 per candidate) do not depend on the draws: computed once, in round 1.
- **Counting.** `eval_count` counts only the evaluations of each candidate's resolved chain (12 per
  candidate, as today), so evaluation ledgers do not change; speculative evaluations are counted
  separately (`computed_evals`) and reported in Step 4.
- **k** is chosen from Step 0's failure rates (fixed per run, a config value with default 1 if Step 0
  shows few failures; otherwise the smallest k that keeps the round count small). Larger k means
  larger batches, which suits the GPU.

The batched path is opt-in, so goldens and default results do not change; with the key on, the
quick MC's trial rotations equal today's for every candidate.

# Steps

**Step 0: go/no-go measurement (script only, `scripts/speed_parity/`, no library change).**
On Example2.ThreeVoxels voxel 0 and 10 voxels of the E0 clean set:

- quick-MC candidates per level (the batch sizes), the fraction of failed blocks per level, and the
  round count option A would need for k = 1, 2, 4 (computed from the sequential outcomes);
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
entry only if Step 4 shows stage D matters. Covers the default cost and `CppOverlapRule`;
`CppEdgeRules` after A-exact merges (see Status), the ray rule falls back to `evaluate`. Reuses S1's lookup and image tables if S1 has
landed. `device` argument: `cpu` (default), `mps`, `cuda`.

Tests (`tests/test_batched_cost.py`): on a fixed seeded set of orientations near and far from the
truth, every `OverlapInfo` field equal to `evaluate`'s on CPU. If some differ by a float32 last bit
in another batch shape (as for `BatchedScreen`), count them, name the stage, and stop for the owner
rather than loosening the test. On MPS (skipped when unavailable): count the orientations whose
counts differ.

**Step 2: lockstep quick MC (new function in a new module `icenine/batched_quick_mc.py`).**
Inputs: the level's candidates, `trial_box_width`, `trial_step`, the cost evaluator, the generator,
k. Implements the rounds of option A above. Per step: build the trial rotations of every running
chain from the buffered uniforms (`get_near_identity_point`,
`_quat_multiply`, `quaternion_to_matrix` as today; vectorise only where a test shows bit-identity,
else loop per candidate, which is cheap), one `evaluate_many`, strict `<` acceptance per candidate.
Fills `orientation`, `cost`, `overlap_info` and a per-candidate `last_run` as `optimize` does.
Applies only when the run is a single block (`n_min_ergodic >= max_mc_steps`) and `level_mc` is a
plain `MCOptimizer` (`type(level_mc) is MCOptimizer`; study subclasses such as `LoggedMC` in
`scripts/finisher_diagnosis/diagnose.py` keep the sequential path); otherwise sequential, with a
warning only when the key was set explicitly (the `BatchedScreen` fallback rule).

Tests (`tests/test_batched_quick_mc.py`):

- With the same seeded generator, the lockstep result equals sequential `optimize` for every
  candidate (orientation, cost, `OverlapInfo`, `last_run`) for k = 1, 2 and N, first with `evaluate`
  as the cost (logic only), then with `BatchedCostEvaluator` on CPU. Include levels with no failed
  block, one, and many (a cost stub that forces failures covers the many-failure case).
- The generator state after the level equals the sequential run's.
- The buffered, mapped uniforms equal today's scalar `rng.uniform(-r, r)` draws bit for bit.
- `eval_count` equals the sequential run's; `computed_evals` is at least that.

**Step 3: wiring.** Config keys `BatchedQuickMC` (default 0) and `SearchDevice` (default `cpu`) in
`config_file.py` and `SearchParameters`; the level loop in `reconstructor.py` calls the lockstep
function when the key is on. Pin `BatchedQuickMC` off where evaluations are counted by patching
`evaluate` (`scripts/eval_reduction/eval_ledger.py` and the other scripts that pin `batched_screen`,
as eb1a98d did), or count `evaluate_many` there; evaluation counts must stay the same with the key on.

**Step 4: measurements.** Quiet machine, single worker, key off vs on, CPU vs MPS:

- ThreeVoxels voxel 0 no-start search: wall time, quick-MC time per level, µs per evaluation
  (`scripts/speed_parity/`);
- E0 clean set (per-voxel images, `scripts/findoptimal_robustness/e0_run.py`), 200 voxels, seed 0, key
  off vs on (CPU): wrong counts, counted and computed evaluations side by side,
  per-level candidate counts, and the number of voxels whose answer changed (expected 0 with option
  A unless Step 1 found batched costs that differ from `evaluate`; trace any change).

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
