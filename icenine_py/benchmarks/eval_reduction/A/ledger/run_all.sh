#!/bin/bash
# Task A step 1: the evaluation ledger, 6 voxels x 4 arms, one process at a time (single-worker timing).
# Run from icenine_py/:  nohup bash benchmarks/eval_reduction/A/ledger/run_all.sh > benchmarks/eval_reduction/A/ledger/logs/run_all.log 2>&1 &
set -u
D=benchmarks/eval_reduction/A/ledger
mkdir -p $D/runs $D/logs
# batched_screen=0: the recorded arms ran the per-orientation screen (BatchedScreen is default on since 2026-10-13)
for arm in "base:--opt batched_screen=0" "cpp:--opt cpp_overlap_rule=1 --opt batched_screen=0" "r2:--opt ray_rule=R2 --opt batched_screen=0" "cpp_edge:--opt cpp_overlap_rule=1 --opt cpp_edge_rules=1 --opt batched_screen=0"; do
  tag=${arm%%:*}; opts=${arm#*:}
  for v in 15901 17034 18558 19889 22867 5242; do
    [ -f $D/runs/${tag}_v$v.json ] && continue
    uv run python scripts/eval_reduction/eval_ledger.py run --voxel $v --tag $tag $opts --dir $D
  done
done
touch $D/logs/run_all.done
