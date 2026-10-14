---
title: "TODO: batched quick MC with a GPU option (ROADMAP S3)"
subtitle: "Plan: run every candidate's quick MC of a level in lockstep through one batched local cost"
date: "2026-10-13, revised 2026-10-14 (draw option B)"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Approved 2026-10-13; revised 2026-10-14 to draw option B (owner, after Step 0); Steps 1-5 approved
to run. Steps 1-5 built and measured 2026-10-14 (MIGRATION_HISTORY "S3: batched quick MC"); awaiting review.** Step 0 (done, 4686f4d, report `.claude/reports/20261014-0010-s3_step0.md`) found option A
(exact draws by speculation) a no-go: 30-70% of quick-MC blocks fail, so the speculative rounds cost
21-67 times the counted evaluations (5.5-13 times even with reuse between offsets) and the quick MC
gets slower, not faster. The owner chose option B: each candidate's draws are fixed in advance and do
not depend on other candidates' outcomes. With the key on, results differ from today's, so Step 4
compares wrong counts (an A/B), not equality. Track S, Ilium Agent. Task branch
`feature/python-speed-parity-s3` (exists; Step 0 committed). An `implementer` agent builds Steps 1-5;
`/review-task` (Opus) when done, required fixes to a fresh implementer (Track S rule). `/handoff`,
`/clear`, `/resume` afterwards (protocol rule 5).

**Why S3 can start before A-exact lands.** A-exact changes `cost_functions.py` and `_rasterize.c`
only under the opt-in `CppEdgeRules` (default path bit-identical). S3 adds new modules and does not
edit `cost_functions.py`; its batched cost calls the existing C `stage_d_overlap`, so A-exact's C
change reaches it when merged. Until then the batched path covers the default cost and
`CppOverlapRule`, and falls back to `evaluate` for `CppEdgeRules` (and the ray rule).

**Before starting:** done (fetched; `feature/eval-reduction` merged at 1729aea; the heads-up on the
shared files `reconstructor.py` and `orientation_search.py` is pushed). Timings: single worker, not
during the A/B runs, a Phase D render or a Q8 run. The desktop apps keep the 1-minute load at 2-3.5
on Ilium Agent's machine, so `require_quiet` refuses: run with `--allow-busy`, label the timings
contended, and repeat them twice (as Step 0 did).

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

# Step 0 results (2026-10-14; contended, one torch thread)

- **Batch sizes.** Quick-MC candidates per level, median: E0 clean (10 voxels x 3 seeds, Q-max 8)
  168 / 44 / 11 / 2; ThreeVoxels Q8 397 / 105 / 27 / 6.
- **Failed blocks** (no trial below the start cost): E0 32% / 56% / 61% / 37% by level; ThreeVoxels
  Q8 47-68%. Every quick MC drew 30 uniforms, or 33 when its block failed.
- **Quick-MC share** at Q-max 8: 50-60% of wall time, 69-74% of `evaluate` calls.
- **Batched cost against `evaluate`** (326-349 µs at Q-max 8, K = 112 reflections):
  - CPU compact variant (data-dependent shapes): 55-72 µs per orientation at B >= 16 (4.8-5.9
    times), 246 µs at B = 1.
  - Fixed-shape variant on CPU: 95-109 µs (3.0-3.6 times).
  - MPS fixed shape: 29-34 µs at B >= 256 (10-11 times), 56-62 µs at B = 64, about 2.3 ms per call
    at B = 1 (slower than the CPU variant below a batch of about 64).
  - All reflections (K = 868): CPU compact 1.7-2.0 times, MPS 3.8 times.
- **Equality.** Both CPU variants equal `evaluate` in every `OverlapInfo` field on 1024/1024
  orientations, provided g is computed with one 3x3 gemm per orientation against the `g_hkl.T` view
  and kept component-major (3, B K). A batched matmul gives 29/1024 differing orientations, a
  contiguous copy of `g_hkl.T` 4/1024. `scripts/speed_parity/s3_step0.py` holds the working forward.
  MPS differs from `evaluate` in a count field on 112-118/1024 orientations (391/1024 at K = 868)
  and in the cost on 5-31/1024 (max |difference| 0.011), so MPS results are compared by wrong
  counts only.
- **Model, not a measurement.** Without speculation the quick MC gets 5.2 times faster on E0 (CPU and
  MPS) and 4.6 / 6.6 times on ThreeVoxels Q8 (CPU / MPS); the whole no-start search 1.7-2.0 times.

# Random draws: option B, a fixed slice per candidate (owner 2026-10-14)

Today candidate i's 30 trial uniforms start at stream position 30 i + 3 F_i, F_i the number of
failed blocks among candidates 0 .. i-1, so each candidate's draws depend on the outcomes before it.
Option A (exact, by speculation) was dropped after Step 0. Option B:

- **Buffer.** At the start of a level's quick MC, draw one buffer of 33 N standard uniforms
  (`rng.random(33 * N)`), N the level's candidate count. Candidate i owns positions 33 i .. 33 i + 32:
  the first 30 are its 10 trials (3 each, in today's order within a candidate), the last 3 the restart
  that the failure branch draws and never uses. The generator advances by exactly 33 N per level
  whatever the outcomes (today 30 N + 3 F), so later levels and FindOptimal also see a stream that
  differs from today's.
- **Mapping.** A standard uniform u maps to the trial range exactly as `Generator.uniform(-r, r)` does
  (`low + (high - low) * u`), r = `tan(step)/sqrt(12)`; a test checks bit identity with scalar
  `rng.uniform` on the same stream.
- **As built (changed in Step 2).** The buffer is `rng.uniform(-r, r, size=33 * N)`, not
  `rng.random(33 * N)` mapped by `low + (high - low) * u`: the two numpy roundings of that expression
  differ from the scalar `rng.uniform(-r, r)` in 10,714 of 20,000 draws (numpy's C code rounds once;
  `benchmarks/s3_steps1_5/draw_check.json`). The stream consumption is the same 33 N doubles, and the
  slices equal the scalar draws bit for bit. The two bullets above describe the plan, not the code.
- **No speculation.** Each step is one batched call over the N candidates; computed evaluations equal
  counted ones, 12 per candidate as today, so evaluation ledgers keep their meaning.
- **Sequential reference** (tests only): today's `MCOptimizer.optimize` run per candidate with its
  generator replaced by a replay of that candidate's slice (a small test helper serving the `uniform`
  calls `optimize` makes, with the same mapping). The lockstep result must equal it exactly. The
  library's sequential path does not change.

Option B is on only with `BatchedQuickMC`; goldens and default results do not change.

# Steps

**Step 0: done** (see "Step 0 results"; script `scripts/speed_parity/s3_step0.py`, data
`benchmarks/s3_step0/`).

**Step 1: `BatchedCostEvaluator` (new module `icenine/batched_cost.py`).**
Built from one `VoxelCostFunction`; `evaluate_many(R (B,3,3), voxel_vertices, phase_index) ->
List[OverlapInfo]`, adds B to `cost_fn.eval_count`. Two forward variants, both taken from Step 0's
script:

- **CPU:** the compact variant, with the gemm layout Step 0 found bit-identical.
- **GPU** (`mps`, `cuda`): the fixed-shape variant (all 2K omega rows per orientation, validity as a
  mask, one packed transfer to the host), used only for batches of at least `GPU_MIN_BATCH`
  orientations (module constant, 64 from Step 0's crossover; Step 4 checks it). Smaller batches use
  the CPU variant.

Stage D stays in C on CPU: one `stage_d_overlap` call per orientation, quality accumulated exactly as
today; a single batched C entry only if Step 4 shows stage D dominates. Covers the default cost and
`CppOverlapRule`; `CppEdgeRules` after A-exact merges (see Status); the ray rule falls back to
`evaluate`. `device` argument: `cpu` (default), `mps`, `cuda`. Reuses S1's tables if S1 has landed.

Tests (`tests/test_batched_cost.py`): on a fixed seeded set of orientations near and far from the
truth, at B = 1, 7 and 64, every `OverlapInfo` field equal to `evaluate`'s for the CPU variant and for
the fixed-shape variant run on CPU. If any differ, count them, name the stage, and stop for the owner
rather than loosen the test. On MPS (skipped when unavailable): count and report the orientations
whose fields differ; no equality assertion.

**Step 2: lockstep quick MC (new module `icenine/batched_quick_mc.py`).**
Inputs: the level's candidates, `trial_box_width`, `trial_step`, the cost evaluator, the generator.
Draws the level's buffer (option B); evaluates both start evaluations of every candidate (as today:
the matrix, then the quaternion round trip); then 10 steps, each: build every candidate's trial
rotation from its slice (`get_near_identity_point`, `_quat_multiply`, `quaternion_to_matrix` as
today; vectorise only where a test shows bit identity, else loop per candidate), one `evaluate_many`
over N, strict `<` acceptance per candidate. Fills `orientation`, `cost`, `overlap_info` and a
per-candidate `last_run` as `optimize` does. Applies only when the run is a single block
(`n_min_ergodic >= max_mc_steps`) and `type(level_mc) is MCOptimizer` (study subclasses such as
`LoggedMC` in `scripts/finisher_diagnosis/diagnose.py` keep the sequential path); otherwise
sequential, with a warning only when the key was set explicitly (the `BatchedScreen` fallback rule).

Tests (`tests/test_batched_quick_mc.py`):

- The lockstep result equals the replayed sequential reference for every candidate (orientation,
  cost, `OverlapInfo`, `last_run`), first with `evaluate` as the cost (logic only), then with
  `BatchedCostEvaluator` on CPU. Include levels with no failed block, one, and many (a cost stub that
  forces failures covers the many-failure case).
- The generator state after the level equals the state after `rng.random(33 * N)` from the level's
  start state.
- The mapped uniforms equal scalar `rng.uniform(-r, r)` draws bit for bit.
- `eval_count` equals the sequential reference's (12 per candidate).

**Step 3: wiring.** Config keys `BatchedQuickMC` (default 0) and `SearchDevice` (default `cpu`) in
`config_file.py` and `SearchParameters`; the level loop in `reconstructor.py` calls the lockstep
function when the key is on. Pin `BatchedQuickMC` off where evaluations are counted by patching
`evaluate` (`scripts/eval_reduction/eval_ledger.py` and the other scripts that pin `batched_screen`,
as eb1a98d did), or count `evaluate_many` there. The full test suite passes with the key off and the
goldens do not change.

**Step 4: measurements.**

- **Timing first** (single worker, none of our jobs running; contended label, two repeats):
  ThreeVoxels voxel 0 no-start search and 10 E0 clean voxels, key off vs on, CPU vs MPS: wall time,
  quick-MC time per level, µs per evaluation, and the MPS crossover batch size (check
  `GPU_MIN_BATCH`).
- **A/B on wrong counts** (after the timing): E0 clean set, 200 voxels x seeds 0-2, key off vs key
  on (CPU), with `scripts/findoptimal_robustness/e0_run.py` (add a key-on flag and a separate cache
  tag; reuse its images). At most 10 workers, nohup, resumable. Ilium Agent's machine has no E0
  cache, so run both arms; the key-off arm should reproduce the published E0 clean result (204/600
  wrong, MIGRATION_HISTORY "Diagnosis"); say whether it does. Report wrong (> 1 deg) per seed and
  pooled with Wilson 95% intervals, counted evaluations per run, per-level candidate counts, wall time
  per run, and the voxels whose answer flips (right to wrong, wrong to right). No p-values.
- **MPS arm:** seed 0 on the first 50 of the 200 voxels, key on with `SearchDevice mps`, 1-2 workers:
  wrong counts against the CPU key-on arm on the same runs, and the time per voxel.
- **Budget:** pilot 5 voxels first and estimate. If an arm of 600 runs would take over 4 h on 10
  workers, run seed 0 only and say so.

**Step 5: docs.** MIGRATION_HISTORY section (Step 0 and Step 4 results, the draw option and why A
was dropped, fallbacks), README (new modules, keys), ROADMAP S3 row.

# Rules (ROADMAP "How to use this file")

New keys opt-in and off; no default flips and no golden regeneration (list them for the owner);
evaluations and wrong counts side by side; no p-values; clean images only (Q8 open). No large traces
in new commits (Step 0's committed data stays as it is until the owner decides). Merge into
`feature/python-speed-parity` with `--no-ff` after review; push and post mailbox entries.

# Out of scope

The batched screen on a GPU (a natural follow-up: level 0 has 43,974 orientations per voxel, a
better fit for a GPU than the quick MC's batches); FindOptimal and VarianceMinimizing (sequential by
nature); multi-block quick MC (no current config uses it); S2.
