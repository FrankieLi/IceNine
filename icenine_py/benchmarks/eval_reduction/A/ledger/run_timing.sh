#!/bin/bash
# Task A step 3: quiet single-worker repeat of the two ends of the speed table (C++ rule + edge rules, switches off and on),
# on voxels 15901 and 5242, each run started only when the 1-minute load is below 1.4 so the preflight passes.
# Run from icenine_py/ with nothing else running.
set -u
D=benchmarks/eval_reduction/A/ledger
wait_quiet() { until uv run python -c "import os,sys; sys.exit(0 if os.getloadavg()[0] < 1.4 else 1)"; do sleep 10; done; }
for v in 15901 5242; do
  for arm in "timing_cpp_edge:--opt cpp_overlap_rule=1 --opt cpp_edge_rules=1 --opt batched_screen=0" \
             "timing_cpp_edge_interp_batched:--opt cpp_overlap_rule=1 --opt cpp_edge_rules=1 --opt cpp_omega_interp=1 --opt batched_screen=1"; do
    tag=${arm%%:*}; opts=${arm#*:}
    [ -f $D/runs/${tag}_v$v.json ] && continue
    wait_quiet
    uv run python scripts/eval_reduction/eval_ledger.py run --voxel $v --tag $tag $opts --dir $D
  done
done
touch $D/logs/run_timing.done
