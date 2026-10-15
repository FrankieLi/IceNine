---
title: "TODO: remove the per-call Python loops in evaluate (ROADMAP S1)"
subtitle: "Plan: stage-A lookup table, one frame table per cost function, hoisted Bragg constants; bit-identical"
date: "2026-10-14"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Built 2026-10-14 (Steps 0 to 6, Step 4 included), reviewed (REQUEST CHANGES) and the review fixes applied the
same day; not merged.** Quiet-machine repeat not done (preflight not quiet), owner decision D-20261014-I16: the
timings stay labelled single-worker, contended. Results: MIGRATION_HISTORY
"S1: per-call loops out of evaluate" and `benchmarks/s1_loops/`. Approved 2026-10-14, Step 4 included
(D-20261014-I15). The owner picked S1 as the next
Track S item on 2026-10-14 (D-20261014-I14). A-exact and S3 are both merged into `feature/python-speed-parity`, so nothing
blocks it. Branch: task branch `feature/python-speed-parity-s1-loops` off
`feature/python-speed-parity`. An `implementer` agent builds it, `/review-task` reviews it, and a
fresh implementer applies the fixes (Track S rule). `icenine/cost_functions.py` is a shared file, so
the heads-up goes in the Ilium mailbox and is pushed before the task starts.

# What

Make `VoxelCostFunction.evaluate` cheaper with three changes that leave its results bit-identical:

1. **Stage-A lookup table.** `calculate_diffraction_overlap_batched` maps each omega to a frame with
   a Python loop over the peaks (`cost_functions.py:561-564`). Replace the loop with one vectorised
   lookup in a frame-index array built once per cost function. `BatchedCostEvaluator._wedge`
   (`batched_cost.py:227-231`) already does this. Both should use one shared helper.
2. **One frame table per cost function.** Every call now builds a list of binary frames: it runs
   `np.unique`, a Python re-indexing loop, and one `exp_data.get_image(w, d).get_binary_numpy()` per
   frame and detector (174 calls on the S0 arm; `cost_functions.py:729-740`). Instead, keep one
   table indexed by `[frame * n_det + det]`, filled lazily per frame. Pass it to `stage_d_overlap`
   and `ray_flags_c` with the global frame indices. `BatchedCostEvaluator` already does this
   (`_table`, `_filled`, `_ensure_frames`, `batched_cost.py:131-243`). Move that code into one
   helper used by both. The C extension needs no change, because it only indexes the list.
3. **Hoisted Bragg constants.** On every call, `get_scattering_omegas_torch` rebuilds tensors that
   do not change between calls: the beam energy and chi as tensors, their expansion, the wavenumber,
   `sin_theta = |g| / 2k`, `|g| + eps`, and `sin` and `cos` of chi. Compute these once per phase and
   pass them to a new function next to the old one in `diffraction_core.py`. The new function also
   skips zeroing the omegas of unobservable rows, because `evaluate` drops those rows anyway.
   `get_scattering_omegas_torch` stays unchanged, since `differentiable_cost.py`,
   `batched_screen.py` and `batched_cost.py` also use it.

There is no new config key. The default path changes code but not results, and Step 5 proves that on
the default path and on every opt-in mode.

# Why

Of the ~318 µs of one re-cost, quick-MC or final `evaluate` call, three parts can go without
touching the physics (MIGRATION_HISTORY "Python speed parity: per-evaluation cost (2026-10-13)",
stage split):

| part | µs now | after (estimate) |
|---|---|---|
| stage A Python loop over peaks | 18-20 | ~2 (measured for a lookup table) |
| frame list loop (174 `get_image` calls) | 29 | a few (one identity check, one fill check; the C call's pointer prefetch, see Risks) |
| Bragg solve | 41 | ~25-30 (about 10 of its ~30 tensor ops are constant) |

That is about 55-60 µs of ~318 per call, so the ROADMAP's figure of **about 1.2x per call** holds.
Under the default keys `evaluate` takes about 89% of a no-start search (`BatchedScreen` on), so the
whole search should get roughly 1.15x faster. The default quick MC, re-cost, FindOptimal and final
cost all call `evaluate`. The default screen does not, so S1 does not speed it up. The rest of the
per-call cost (about 200 µs of tensor dispatch in the eta filter and stages B-C) is S2's, not S1's.

# Steps

**Step 0: baseline on the base commit (before any library change).** First write
`scripts/speed_parity/s1_parity.py`, the parity harness. It records every `OverlapInfo` field
(`pixel_overlap`, `pixel_on_detector`, `peak_overlap`, `peak_on_detector`, `quality`,
`n_quality_points`, `detectors_overlap`, `ray_cases`) of `evaluate` for a fixed orientation set and
saves it as `.npz` under `scripts/speed_parity/cache/` (gitignored). It runs single-threaded
(`torch.set_num_threads(1)`).

- Orientations: 500 within 3 degrees of the truth and 500 uniform over the fundamental zone, with a
  fixed seed. Data: Example2.ThreeVoxels (voxel 0) and three interior full-sample voxels on the clean
  Phase D images (Q-max 8). 1,000 orientations of ~220 omega rows should reach the edge bins of
  the scan, `[low - width, low)` and above the range, many times. The harness counts those rows so
  the coverage is shown, not assumed.
- Modes, each a separate `VoxelCostFunction`:
  - the local cost (pixel radius 0, the default rule)
  - `CppOverlapRule`, and `CppEdgeRules` with and without it
  - the screen (pixel radius 3, Q-max 5) with the centroid screen and with vertex 0
  - the ray rule (C flags, and `USE_C_RAY_FLAGS = False`)
  - `min_sin_eta > 0`
  - soft mode
  - the Python fallback path (`_HAS_C_RASTERIZE` patched off; 50 orientations, because it is slow)
- Record the timings on a quiet machine (`scripts/dev/timing_preflight.py`): `prof_eval.py` (µs per
  call for the three arms), `stage_split.py`, and `search_split.py` (whole no-start search, evaluation
  count, result).
- E0 baseline: `e0_run.py run --variants clean --n-seeds 1 --tag _s1base` (200 voxels, seed 0, clean,
  default keys; about 9 minutes with 10 workers). Run it after the timings, not during them.
- Commit the script and a small fixture, `tests/data/s1_evaluate_ref.npz` (64 orientations × the
  modes above, from ThreeVoxels, under 200 KB). It is new, not a regenerated golden.

**Step 1: stage-A lookup table.** Write one helper, `frame_of_omegas(omegas, low, width, index,
cpp_edge_rules)`, built on `_omega_bins`. It returns -1 when an omega has no frame: a bin outside
`[0, n_bins)` or an `index_list` entry of `None`.
- `VoxelCostFunction.__init__` builds the `int64` index array once and passes it to
  `calculate_diffraction_overlap_batched` through a new optional keyword. A caller that passes no
  array gets one built from `range_map` (vectorised, without the per-peak loop).
- `BatchedCostEvaluator._wedge` calls the same helper.
- The Python fallback path uses the helper too.

**Step 2: frame table.** Write a `FrameTable` helper in `cost_functions.py` holding the
`[frame * n_det + det]` list of binary frames, a filled-per-frame flag array, and the detector-0
frame size. Today's code reads the frame size from detector 0's first frame; the fill asserts that
every detector-0 frame has that same size. The fill takes `get_image(w, d).get_binary_numpy()` per
missing frame, as `_ensure_frames` does.
- **Identity:** the table belongs to one `exp_data` object. `evaluate` checks
  `cost_fn.exp_data is table.source` and rebuilds the table when they differ.
  `scripts/optimizer_sweep.py:442` swaps `exp_data` on a live cost function, and the swap must keep
  working.
- **Contract:** the images do not change after the first `evaluate` of a cost function. The
  implementer greps `icenine/`, `scripts/` and `tests/` for image mutators (`set_pixel`,
  `add_triangle`, `_binary_cache` writes, ...) and `exp_data` swaps applied to a live cost function,
  and lists every hit with how the table handles it.
- **Fast path:** pass the table and the global frame indices (`int32`) to `_c_stage_d_overlap` and
  `ray_flags_c`. The `np.unique`, local re-index and `get_image` loop go.
- **`BatchedCostEvaluator`:** uses the cost function's table instead of its own `_table`, `_filled`
  and `_ensure_frames`. Same arrays, so its results are unchanged.
- The slow Python path stays as it is apart from Step 1. It needs `ImageData` objects for soft mode.

**Step 3: hoisted Bragg constants.** Add `BraggConstants` and `scattering_omegas_fixed_beam(g,
consts)` to `diffraction_core.py`.
- The constants are built per phase in `VoxelCostFunction.__init__` from `g_mag`, `beam_energy` and
  `beam_deflection_chi`, with the same float32 expressions in the same order as
  `get_scattering_omegas_torch`.
- `evaluate` rebuilds them only if `simulator.beam_energy` or `beam_deflection_chi` changed (two
  float comparisons per call).
- The function returns `observable`, `omega1` and `omega2`. On observable rows it equals
  `get_scattering_omegas_torch` bit for bit. Unobservable rows are unspecified, and the docstring
  says so.

**Step 4 (owner: in or out): `CppEdgeRules` in the batched local cost.** This is S3's open
follow-up. Today `make_batched_cost_evaluator` refuses `CppEdgeRules`, and the quick MC falls back to
the sequential loop with a warning.
- Change: `_stage_bc` applies the "coordinate below 0 is pixel -1" rule (`_pixel_coordinate`), the
  frame lookup already floors (Step 1 helper), and the vertex-lit rule is already in `rule_bits`.
  `unsupported_reason` then needs only a `_C_HAS_VERTEX_LIT` build.
- Check: the compact CPU variant, and the fixed variant on the CPU, equal `evaluate` in every field
  with `CppEdgeRules`, alone and with `CppOverlapRule`, on the Step 0 orientation set.
- Effect: runs with both opt-in keys (`BatchedQuickMC` and `CppEdgeRules`) then use the lockstep
  quick MC and its draws (option B), so their results change. Neither key is on by default.
- Leave it out if the owner wants S1 to be strictly bit-identical everywhere.

**Step 5: checks.**
- `s1_parity.py compare`: every field of every mode equal to the Step 0 dump, exactly (`==`, also for
  `quality`), on all four voxels.
- The new tests pass:
  - the lookup table against the old loop (kept in the test as the reference) on synthetic ranges
    with `None` entries, edge omegas and both edge rules
  - `FrameTable` lazy fill and rebuild after an `exp_data` swap
  - `scattering_omegas_fixed_beam` against `get_scattering_omegas_torch` on observable rows (random
    g, chi 0 and not 0, two energies)
  - `evaluate` against the fixture, exactly
  - with Step 4, the batched cost against `evaluate` with `CppEdgeRules`
- The full suite passes (`uv run pytest tests/`).
- The same timings as Step 0, on a quiet machine, reported side by side: µs per call per arm, stage
  split, whole search. The search's evaluation count and result must be identical.
- E0 `--tag _s1` against `_s1base`: identical per-voxel records (orientation, cost, evaluations per
  level, wrong flag). Report wrong counts and evaluations side by side (both unchanged).

**Step 6: docs.** Add a MIGRATION_HISTORY section "S1: per-call loops out of evaluate" with the
timing table before and after, the parity result and the mutation list from Step 2. Update ROADMAP
row S1 (and S3's open list if Step 4 runs), and the README speed-parity paragraph (`s1_parity.py`,
the fixture, the two helpers). After the merge, post a mailbox heads-up.

# Risks

- **Stage D pointer prefetch.** `stage_d_overlap` type-checks and prefetches every entry of the
  image list on each call (`_rasterize.c:717-740`): n_frames × n_det entries instead of only the
  frames in use. That is 360 entries for a 180-frame, two-detector scan, probably 1-2 µs. Step 5
  measures it. If a long scan makes it a real cost, report it. S1 does not touch `_rasterize.c`; a
  cached pointer array would be a later item.
- **Images changed under a live cost function.** Today's code reads the binary caches on every call,
  so it would see a mutation; the table would not. Step 2's grep lists every case. The swap of the
  whole `exp_data` object is handled.
- **Float order.** Bit-identity depends on computing the hoisted constants with the same float32
  operations in the same order. Steps 3 and 5 test this directly rather than assume it. Thread count
  can matter for the 3x3 @ 3xK product, so parity runs use a fixed thread count on both sides.
- **The gain is small.** About 1.2x per call. The large part of the gap (tensor dispatch) is S2's,
  and S3 already batches the default-off quick MC. If Step 0 finds the three parts much smaller than
  the S0 split, report that before going on.

# Rules

- No new switch and no default flip: the change is bit-identical, and Step 5 proves it.
- No golden regeneration; the fixture is a new file.
- `_rasterize.c` is unchanged, so no rebuild is needed.
- Track E's shared files: only `cost_functions.py` changes. `batched_screen.py`'s per-call frame
  list (Newt Agent's `ScreenEvaluator`) is a candidate for the same table, but is left for Track E
  and listed in the mailbox.
- Timings only on a quiet machine, with nothing else running. E0 runs never overlap them.

# Relation to other work

- ROADMAP Track S rows S1 to S3.
- MIGRATION_HISTORY "Python speed parity: per-evaluation cost (2026-10-13)" (the stage split) and
  "S3: batched quick MC with a GPU option" (`BatchedCostEvaluator`, whose table and lookup this
  shares).
- Task A-exact (`CppEdgeRules`, `HAS_VERTEX_LIT`) for Step 4.
