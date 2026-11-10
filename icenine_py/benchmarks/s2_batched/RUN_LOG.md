# S2 Step 3 (the batched fused entry): run log

Branch `feature/python-speed-parity-s2-batched`; plan `docs/todo_s2_fused_c_eval.md`, section
"Step 3". The S2 results are written up in MIGRATION_HISTORY "S2" by Steps 4-6.

S1's parity set throughout: ThreeVoxels voxel 0 (`three0`) and Phase D voxels 15901, 17034,
18558, 1,000 orientations each (500 near the truth, 500 from the fundamental zone;
`s1_parity.orientation_set`), single-threaded.

## Single-orientation fused path unchanged by the refactor

`_fused_eval.c`'s per-orientation body moved into shared static functions. The single-orientation
fused `evaluate` (key on) dumped before (base `82a54da`) and after, every covered mode:

    uv run python scripts/speed_parity/s2/s2_batched_check.py dump-fused --out .../fused_base.npz
    uv run python scripts/speed_parity/s2/s2_batched_check.py dump-fused --out .../fused_new.npz
    uv run python scripts/speed_parity/s1_parity.py compare .../fused_base.npz .../fused_new.npz

Result `fused_single_compare.txt`: bit-identical, 36 groups (4 datasets x 9 modes).

## Batched fused against single fused (key on)

    uv run python scripts/speed_parity/s2/s2_batched_check.py batched --ref .../fused_new.npz \
        --json benchmarks/s2_batched/batched_vs_single.json

`BatchedCostEvaluator(fused_cost=True).evaluate_many` in chunks of B = 1, 7, 64 and 1,000 (a
fresh cost function per row, so the first call of each row fills the FrameTable), the 5 modes the
batched cost supports. Result: every field equal in 80 of 80 rows; every orientation went through
the fused entry and was counted once.

## Key off against the base commit

    uv run python scripts/speed_parity/s1_parity.py dump --out .../base.npz   # at 82a54da
    uv run python scripts/speed_parity/s1_parity.py dump --out .../new.npz
    uv run python scripts/speed_parity/s1_parity.py compare .../base.npz .../new.npz
    uv run python scripts/speed_parity/s1_parity.py batched --dump .../base.npz

`keyoff_compare.txt`: bit-identical, every field of every mode (60 groups: 4 datasets x 15
modes). The commit message of `98faa56` and the implementer report say 61 groups in error; the
file has 60 rows, and the review's re-run also gives 60. `keyoff_batched.txt`: the compact and
fixed variants equal the base `evaluate` in every field (40 rows).

(`...` = `scripts/speed_parity/cache/s2_batched`, gitignored.)

## Deviations and deferrals

Design choices that differ from the plan's wording (accepted in the review and by D-20261028-I1):

- The fused path is dispatched in `BatchedCostEvaluator.evaluate_many` (a new `_evaluate_fused`),
  not inside `_forward_compact`: that function returns the stage-D inputs, while the fused entry
  includes stage D.
- `variant="fused"` forces the fused path (even with the key off), as "compact" and "fixed" force
  theirs; it falls back to the compact variant where the fused evaluate does not cover the cost.
- `HAS_FUSED` stays 1; `batched_cost.py` imports `fused_evaluate_many` in a guarded import, and a
  Steps 1-2 build without it falls back with its own warning reason.
- One declaration line for `fused_evaluate_many` was added to `_stage_d.h` (where the other fused
  entry points are declared), so that `_rasterize.c`'s method table sees the prototype.

Deferred to Step 6: the fallback warning wording "using the torch evaluate"
(`cost_functions.py`, outside Step 3's scope; on the batched path the fallback is the torch compact
variant), README rows `_fused_eval.c` (add `fused_evaluate_many`) and `batched_cost.py` (the fused
variant; "the CPU compact variant equals `evaluate`" holds only with the key off), and the README
`FusedCost` wording "Falls back to the torch evaluate". Deferred to Step 5: all timings, including
whether the GPU fixed variant should keep priority at B >= `GPU_MIN_BATCH` (64, measured against
the torch compact variant) when `FusedCost` is on.
