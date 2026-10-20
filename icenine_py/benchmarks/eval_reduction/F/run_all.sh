#!/bin/bash
# Task F: the E0 run (8 near-axis + 16 far-from-axis voxels, 20 directions, all methods), then the
# full-sample run (six interior voxels, 1, 2 and 5 deg, every method except F3). Commands from
# icenine_py/. Logs in benchmarks/eval_reduction/F/cache/logs/.
set -e
L=benchmarks/eval_reduction/F/cache/logs
uv run python scripts/eval_reduction/outward_bench.py run --kind e0 --workers 5 --dirs 20 \
  --n-near 8 --n-far 16 --tag e0 > $L/e0.log 2>&1
uv run python scripts/eval_reduction/outward_bench.py run --kind full --workers 5 --dirs 20 \
  --offsets 2,3,5 --tag full --methods f1,f2_05,f2_1,f2_2,f4,f5,f54 > $L/full.log 2>&1
# the extra method, run after seeing why F4 fails from far starts (RUN_LOG section 6)
uv run python scripts/eval_reduction/outward_bench.py run --kind e0 --workers 5 --dirs 20 \
  --n-near 8 --n-far 16 --tag e0_f4w --methods f4w > $L/e0_f4w.log 2>&1
uv run python scripts/eval_reduction/outward_bench.py run --kind full --workers 5 --dirs 20 \
  --offsets 2,3,5 --tag full_f4w --methods f4w > $L/full_f4w.log 2>&1
