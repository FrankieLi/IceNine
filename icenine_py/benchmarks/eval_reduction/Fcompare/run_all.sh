#!/bin/bash
# Task F comparison (D-20261021-N1): main runs, contended (6 processes at once), then the
# single-worker timing pass (run by hand after this finishes, alone on the machine).
# From icenine_py/. Voxel lists: Fcompare/voxel_sets.json (outward_compare.py sets).
set -u
D=benchmarks/eval_reduction/Fcompare
mkdir -p $D/cache/logs
nohup uv run python scripts/eval_reduction/outward_compare.py run --kind e0 --workers 2 --tag e0 \
    > $D/cache/logs/e0.log 2>&1 &
nohup uv run python scripts/eval_reduction/outward_compare.py run --kind full --workers 4 --tag full \
    > $D/cache/logs/full.log 2>&1 &
wait
