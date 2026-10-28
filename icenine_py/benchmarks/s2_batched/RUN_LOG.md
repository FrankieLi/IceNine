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

`keyoff_compare.txt`: bit-identical, every field of every mode (61 groups). `keyoff_batched.txt`:
the compact and fixed variants equal the base `evaluate` in every field (40 rows).

(`...` = `scripts/speed_parity/cache/s2_batched`, gitignored.)
