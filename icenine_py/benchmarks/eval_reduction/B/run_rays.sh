#!/bin/bash
# Task B steps 1 and 4: rays from the truth (full-sample voxels, then E0 voxels). From icenine_py/.
set -e
L=benchmarks/eval_reduction/B/cache/logs
uv run python scripts/eval_reduction/cost_edges_run.py run --kind full --workers 8 > $L/rays_full.log 2>&1
uv run python scripts/eval_reduction/cost_edges_run.py run --kind e0 --workers 8 > $L/rays_e0.log 2>&1
touch $L/rays.done
