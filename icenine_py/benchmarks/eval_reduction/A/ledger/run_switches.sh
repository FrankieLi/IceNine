#!/bin/bash
# Task A step 3: the ledger with the screen switches on (6 voxels x 4 arms), after run_all.sh, one process at a time.
# Run from icenine_py/ after run_all.done exists.
set -u
D=benchmarks/eval_reduction/A/ledger
while [ ! -f $D/logs/run_all.done ]; do sleep 20; done
for arm in "cpp_edge_interp_batched:--opt cpp_overlap_rule=1 --opt cpp_edge_rules=1 --opt cpp_omega_interp=1 --opt batched_screen=1" \
           "cpp_edge_batched:--opt cpp_overlap_rule=1 --opt cpp_edge_rules=1 --opt batched_screen=1" \
           "cpp_edge_interp:--opt cpp_overlap_rule=1 --opt cpp_edge_rules=1 --opt cpp_omega_interp=1" \
           "lenient_interp_batched:--opt cpp_edge_rules=1 --opt cpp_omega_interp=1 --opt batched_screen=1"; do
  tag=${arm%%:*}; opts=${arm#*:}
  for v in 15901 17034 18558 19889 22867 5242; do
    [ -f $D/runs/${tag}_v$v.json ] && continue
    uv run python scripts/eval_reduction/eval_ledger.py run --voxel $v --tag $tag $opts --dir $D
  done
done
touch $D/logs/run_switches.done
