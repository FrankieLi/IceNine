---
title: "TODO: S2: fused C evaluation behind an opt-in key"
subtitle: "One C call for the Bragg solve, eta filter, frame lookup, ray hits and stage D"
date: "2026-10-21"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Approved 2026-10-21 (D-20261021-I4). Step 0 done 2026-10-21 (speed go: 14.4 us per call against
234.1; parity rule revised); Steps 1-6 run under owner delegation (D-20261027-I1), in four tasks
(below), each reviewed and merged before the next that depends on it.** Opus at high effort
implements every task. S2 picked by the owner 2026-10-21 (D-20261021-I2). S3, which
had to come first (D-20261013-I5), is merged. Written by Ilium Agent from a code investigation;
the file and line references are on `feature/python-speed-parity` at ac2b49f. Built under the Track
S rule: implementer on a task branch, `/review-task`, required fixes to a fresh implementer. Task
branch `feature/python-speed-parity-s2-fused`. Shared files are touched, so a heads-up goes in the
Ilium mailbox before the build starts (protocol rule 3(a)).

# What

An opt-in config key `FusedCost` (default 0). With it on, one C call replaces the torch work of
`VoxelCostFunction.evaluate` for the modes it covers, and a batched C entry replaces the torch
forward of the batched quick-MC cost.

**Into C** (from `evaluate`, `icenine/cost_functions.py:1132-1260`):

- g = R g_hkl (1171-1172) and the Bragg solve (1177-1182; `diffraction_core.py:288-310`);
- observable rows, plane normals and the two-solution expansion (1185-1202);
- Rz, reflection, eta and the `min_sin_eta` filter (1204-1236);
- the frame lookup (`FrameTable`, 493-544 and 644-659);
- stages B and C, including `_pixel_coordinate` (662-747, 547-555);
- stage D: the per-peak loop body of `stage_d_overlap` (`_rasterize.c:759-909`) moves into one
  static function shared by `stage_d_overlap` and the fused entry, so stage D exists once.

**Batched entry:** a second C function scores B orientations of one voxel in one call;
`BatchedCostEvaluator._forward_compact` (`batched_cost.py:388`) uses it when the key is on.

**Stays in Python or torch:** key dispatch, `eval_count`, building `OverlapInfo`; the lazy
`FrameTable` fill (C returns "needs frames F", Python calls `ensure` and repeats the call); every
mode the C path does not cover; the GPU fixed-shape variant, `ScreenEvaluator`, `SmoothCost`
(`smooth_cost.py`) and `differentiable_cost.py`.

**GPU coexistence:** S2 is a CPU path. With `SearchDevice` set to `mps` or `cuda` and a batch at
or above `GPU_MIN_BATCH`, the torch fixed-shape variant keeps priority (the owner wants a GPU path,
D-20261013-I5).

**Built once per cost function** (a context object, rebuilt on the checks `frame_table()` uses,
1110-1130, and `BraggConstants.is_for`, `diffraction_core.py:249-261`, plus a geometry check):
per phase float32 g_hkl (K x 3), sin theta, |g| + eps, chi_Laue sin and cos; beam direction,
sample-to-lab rotation and translation; eta limit, `min_sin_eta`, rule bits, pixel radius,
centroid flag; frame index array; detector geometry (n_det x 20, the layout at
`_rasterize.c:933`), frame size and a cached image-pointer array. **Per call:** R (3 x 3), the
voxel vertices (3 x 3), the phase.

**Covered in v1:** pixel radius 0 and > 0 (centroid and vertex 0), the default rule,
`CppOverlapRule`, `CppEdgeRules`, `min_sin_eta`, hard mode. **Falls back to torch**, with a
one-time `RuntimeWarning` naming the reason as `batched_cost.py:429-462` does: `RayOverlapRule`,
soft mode, a missing or stale extension. A test fails when `VoxelCostFunction.__init__`
(1010-1027) gains an argument that the fused path's "unsupported" check does not classify, so a
new Track E mode cannot slip into the fused path untested. Unaffected: `CandidateDedup`,
`BatchedQuickMC` draws, `SmoothCostRank`, `CppOmegaInterp` and `BatchedScreen` (screen only).

# Why

- After S1, a no-start search on ThreeVoxels voxel 0 (default keys) takes 5.5 s, of which 84.9%
  is 8,311 `evaluate` calls at 564 us each (MIGRATION_HISTORY "S1", table `s1_search`). Called in
  isolation the three arms take 215-252 us per call (`s1_per_call`).
- The physics plus stage D is about 10 us of that (stage D alone 8.3-8.5 us); C++ takes 42-55 us
  per evaluation (MIGRATION_HISTORY "Python speed parity: per-evaluation cost"). The rest is
  Python and torch dispatch (about 700 ATen operations per call, `s1_per_call`).
- At 25 us per call the same search would take about 1.0 s. **Modelled, not measured**; Step 0
  measures it.
- The batched quick MC (opt-in, S3) costs 77 us per evaluation on CPU: 50-54 us of torch plus
  10-11 us of stage D at B >= 64 (`benchmarks/s3_step0/tables.md`).
- The batched screen is already 17.5 us per orientation (Newt mailbox 2026-10-20), so it is not in
  v1.

# Steps

## Step 0: go/no-go prototype (scripts only, no library change; its own implementer run)

**Library changes:** none in Step 0. All code goes under `scripts/speed_parity/s2/`, including the
scratch C extension and its build; nothing under `icenine_py/icenine/`, `setup.py` unchanged.

- 0a. Time a no-start search split by caller (screen, re-cost, quick MC, FindOptimal,
  VarianceMinimizing, final) on ThreeVoxels voxel 0 and 10 E0 voxels, with default keys and with
  `BatchedQuickMC 1`, wrapping the callers as `scripts/speed_parity/s4/step0_timing.py` and S1's
  search split do. This gives the share S2 can reach, which is not measured after S1.
- 0b. A scratch extension under `scripts/speed_parity/s2/` covering the Q-max 8 triangle cost and
  the radius-3 screen arm, default rule and `CppOverlapRule`.
- 0c. Single worker, after `scripts/dev/timing_preflight.py`, 5 interleaved rounds against torch
  on S1's 400 near-truth orientations; time per call includes building `OverlapInfo`.
- 0d. Differences per intermediate (omega, frame, pixel coordinate, eta) in ulp, and of every
  `OverlapInfo` field, on S1's 1,000 orientations of voxel 0.
- Go: median <= 30 us per call on the triangle arm and every difference traced (see Parity).
  30-50 us: report to the owner before Step 1. Above 50 us: the stop rule.

## Task A check with the prototype (scripts only; runs in parallel with Steps 1-2)

Decided D-20261027-I4. The Step 0 prototype (`scripts/speed_parity/s2/`) against C++ on Task A's
2,230 orientations, as Parity item 2 describes, for the screen and the complete-cost arms, beside
torch with `CppEdgeRules`. First check that `parity_inputs.py` regenerates the stored inputs. Also
trace the one 1.50 px float64-reference outlier of Step 0 (voxel 18558, edge mode), which torch and
the prototype share. A failure here (the prototype agrees with C++ on fewer orientations than torch)
stops S2: Steps 1-2 are not merged.

**Library changes:** none. Scripts and results under `scripts/speed_parity/s2/` and
`benchmarks/s2_taska/`, and a subsection of MIGRATION_HISTORY "S2".

**Process limit exception:** at most 3 heavy processes at once (the Task A check and Steps 1-2
run at the same time and share the machine's 6).

## Steps 1-2: the fused kernel and its wiring (library)

1. The C kernel and the shared stage-D function in `_rasterize.c`. With the key off,
   `scripts/speed_parity/s1_parity.py compare` is exact against the base commit. Try the plain
   operation order for g_lab (Step 0: matches torch on 0.929 of elements against 0.804 for the fmaf
   chain) and keep it only if it gives fewer flips under the parity rule (D-20261027-I5). Portable
   C only, no Accelerate (D-20261027-I3).
2. Wiring: the key in `config_file.py` (as the keys at 472-480), `SearchParameters`
   (`orientation_search.py:156-170`), `_cost_flags` (`reconstructor.py:897-911`); fallbacks with a
   one-time warning; the frame-refill protocol; the classification test of the "unsupported"
   check. Tests for key off (bit-identical), key on (parity rule on a small set) and fallbacks.

**Library changes:** new `icenine/_fused_eval.c`, compiled into the `_rasterize` module (added to
the sources in `setup.py`); in `_rasterize.c` the shared stage-D function, the method table entry
and a `HAS_FUSED` module constant; `cost_functions.py` (fused dispatch in `evaluate`),
`config_file.py`, `orientation_search.py`, `reconstructor.py` (key plumbing only) and the README
config table. Not `batched_cost.py` (Step 3). No new dependency. Do not edit MIGRATION_HISTORY
(the parallel Task A check writes there; Steps 4-6 write the S2 results).

**Process limit exception:** at most 3 heavy processes at once (the Task A check and Steps 1-2
run at the same time and share the machine's 6).

## Step 3: the batched entry (library)

The many-orientation C entry, used by `BatchedCostEvaluator._forward_compact`
(`batched_cost.py:388`) when the key is on; on `mps` or `cuda` at or above `GPU_MIN_BATCH` the torch
fixed-shape variant keeps priority. Tests: key off unchanged, key on equal to the single-orientation
fused path for every orientation of a batch.

**Library changes:** `batched_cost.py`, and the batched entry in `icenine/_fused_eval.c` and the
`_rasterize` method table. No new dependency.

## Steps 4-6: parity, timings, E0 and docs

4. Parity, as below, on the library path. Trace every flipped coordinate, not only the
   worst-by-ulp one per (row, detector) as `s2_parity.trace` did (Task A check review, fix 5).
5. Step 0c repeated on the library path (single-orientation and batched), plus quick-MC us per
   evaluation and the whole search, single-worker and preflight-gated.
6. E0 A/B (key off against on, with and without `BatchedQuickMC`), paired and traced; then
   MIGRATION_HISTORY "S2", the ROADMAP S2 row and the README.

**Library changes:** none expected; a parity or timing finding that needs a library fix is
reported, not fixed in this task.

**Stop rule:** stop and report if the Step 0 prototype takes more than 50 us per call on the Q-max
8 triangle arm; if any fused-versus-torch difference is not a pixel-edge flip within the parity
rule's bound (or is a frame, eta, observability or t > 0 flip); if with the key off anything differs
from the base commit in `s1_parity.py`; or if the fused path agrees with C++ on fewer Task A
orientations than torch does.

# Parity

`evaluate` works in float32: the g vectors and the orientation are cast with `.float()`
(`cost_functions.py:1102-1107, 1171`). Bit-identity with torch is unlikely: torch's sin, cos,
atan2 and asin are vectorised ATen kernels, the 3 x K product goes through BLAS, and batch shape
already flips last bits (MIGRATION_HISTORY, S3 draw notes). The C code uses float32 in torch's
operation order with `FP_CONTRACT OFF` (as `_rasterize.c:935-936`). If 0d finds every bit equal,
S1's exact rule applies; otherwise:

1. **Fused against torch** (rule revised D-20261027-I2, after Step 0 found float32 rounding flips
   up to 312 ulp and 0.0034 px, which the owner judged immaterial). `s1_parity.py` dumps every
   `OverlapInfo` field on S1's set (ThreeVoxels voxel 0 and Phase D voxels 15901, 17034, 18558;
   1,000 orientations each) for the covered modes, which the test lists. Every differing
   orientation must be traced by a per-peak diff to pixel-edge flips: a peak whose torch or fused
   pixel coordinate lies within 0.01 px of an integer (a bug guard, not a precision target). A
   frame, eta, observability or t > 0 flip fails. Report, for the flipped peaks, which side is
   nearer the float64 reference (Step 0: torch 75, fused 65).
2. **Against C++.** Task A's 2,230 orientations (`scripts/eval_reduction/parity_inputs.py`; voxels
   15901, 22867, 5242) against the stored C++ aggregates
   (`benchmarks/eval_reduction/A/step0/results.json.gz`). The fused path must match at least as
   often as torch with `CppEdgeRules` (2,230 of 2,230 complete cost, 2,229 of 2,230 screen;
   MIGRATION_HISTORY "A-exact"). Unverified: that `parity_inputs.py` regenerates the stored inputs
   exactly; the implementer checks this first.
3. **E0.** 200 voxels, seed 0, key off against key on; wrong flips paired and traced (D-20261014-N13
   rule: a break counts against the key only when traced to it). Seed-to-seed scale: 14-23 flips
   each way.

# Rules

- Key off is bit-identical to the base commit (`s1_parity.py`, E0 per-voxel records).
- No golden regenerated and no default flipped (rules for every task, D-20261012-I1).
- Timings single-worker, run alone, preflight-gated, labelled; never during E0 or other runs. If
  the preflight never passes, the timings stay labelled "contended" (D-20261014-I16).
- Report evaluations and wrong counts side by side; no p-values.

# Risks

- **Numerical flips at thresholds.** Task A's one remaining screen difference is a peak 3 ulp
  from a frame boundary; the fused path can add more. Mitigation: the traced-tolerance rule and
  the C++ comparison, which must not get worse.
- **Two implementations to maintain.** Every new Track E cost mode is torch-only until ported; the
  classification test forces the fallback.
- **Geometry changed under a live context.** `Detector.set_position` and `set_orientation`
  (`detector.py:268-294`) can move a detector; the context's geometry check covers it, and the
  implementer greps for callers on search paths (unverified whether any exist).
- **A floor outside `evaluate`.** About 4.2 us per evaluation is spent outside `evaluate` (S0),
  which caps the gain of the batched path.
- **S4.** S4 merged as a study with its library path parked (D-20261021-I3); if that path is
  revived, it takes over the level-0 screen and re-cost and shrinks S2's share.
- **Stale builds.** A stale `_rasterize` build falls back with a warning; README says to rebuild.

# Decisions

- Plan approved, Opus at high effort (D-20261021-I4).
- After Step 0 (D-20261027-I1 to I5): the owner delegates the rest of S2 to Ilium Agent and reviews
  it in full at the end (I1); the parity rule above (I2); portable C, no Accelerate (I3); the Task A
  check with the prototype runs in parallel with Steps 1-2 and gates their merge (I4); the plain
  order for g_lab is tried in Step 1 and kept only if it gives fewer flips (I5).

# Relation to other work

ROADMAP Track S rows S1, S3, S4 and "Coordination between Tracks E and S"; decisions
D-20261013-I5, D-20261021-I2, D-20261021-I3; Task A's harness
(`scripts/eval_reduction/parity_*.py`);
`docs/todo_s1_evaluate_loops.md` and `docs/todo_batched_quick_mc.md` (the S1 and S3 plans).
