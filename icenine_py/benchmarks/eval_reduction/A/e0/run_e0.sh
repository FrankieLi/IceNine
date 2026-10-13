#!/bin/bash
# Task A step 3: E0 clean (200 voxels, per-voxel images, seed 0, mc arm, min_sin_eta 0), 10 workers per arm,
# arms edge, interp, lenint (see scripts/phase_d/e0_cpprule.py). After the ledger runs; wall times are
# contended multi-worker and not a timing claim. Run from icenine_py/.
set -u
D=benchmarks/eval_reduction/A
while [ ! -f $D/ledger/logs/run_switches.done ]; do sleep 20; done
uv run python ../scripts/dev/timing_preflight.py --json $D/e0/preflight_e0.json || true
for arm in interp edge lenint; do
  uv run python scripts/phase_d/e0_cpprule.py run --arm $arm --workers 10 --msin 0 --tag a --variant clean
done
touch $D/e0/run_e0.done
