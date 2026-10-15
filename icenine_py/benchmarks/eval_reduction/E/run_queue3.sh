#!/bin/bash
# Task E: C++-rule full-sample baseline on the six voxels (replaces the 30-voxel C++ run of
# run_queue1.sh, which was stopped before it started: a full-sample probe with 20 re-ranking arms
# takes about 6 minutes per voxel). Starts when the level-0-width landscape is done. From icenine_py/.
L=benchmarks/eval_reduction/E/cache/logs
while [ ! -f $L/landscape_full_L0.done ]; do sleep 30; done
uv run python scripts/eval_reduction/smooth_rank_probe.py run --kind full --rule cpp --spec off \
    --tag base --voxels six --workers 2 > $L/rank_full_cpp_base.log 2>&1
