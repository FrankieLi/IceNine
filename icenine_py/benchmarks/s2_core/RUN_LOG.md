# S2 Steps 1-2 (fused kernel and wiring): run log

Branch `feature/python-speed-parity-s2-core`; plan `docs/todo_s2_fused_c_eval.md`, section
"Steps 1-2". The S2 results are written up in MIGRATION_HISTORY "S2" by Steps 4-6.

## Key off: `s1_parity.py` against the base commit

- Base dump: commit 832db36 (the parent `feature/python-speed-parity` before this task),
  `uv run python scripts/speed_parity/s1_parity.py dump --out scripts/speed_parity/cache/s2core/base.npz`.
- New dump: commit 9e99d63 (FusedCost wiring, default off), same command with `new_final.npz`.
- `uv run python scripts/speed_parity/s1_parity.py compare base.npz new_final.npz`:
  `keyoff_compare.txt` (every field of every mode bit-identical, 4 datasets x 15 modes). The same
  comparison after the stage-D refactor alone (commit e8699fd's C code) was also bit-identical.

## Key on: the library fused path against torch, per g_lab order (D-20261027-I5)

`uv run python scripts/speed_parity/s2/s2_core_parity.py --out benchmarks/s2_core/glab_order.json`
at commit 9e99d63 (clean tree), torch threads 1: S1's orientation set (500 near the truth + 500 FZ
points) on `three0`, `full15901`, `full17034`, `full18558`, the 9 covered modes of
`s2_parity.MODES`, orders `fma` (the Step 0 prototype's fmaf chain) and `plain` (the library
default since this step). Every differing orientation traced with `s2_parity.trace`; the parity
rule (D-20261027-I2) is checked with a 0.01 px bound. `totals` gives the counts; the `fma` order is
also compared bit for bit with the Step 0 prototype (`prototype_differs`).
