#!/bin/bash
# Task C step 1: candidate overlap on the six full-sample voxels (6 parallel single-thread processes) and
# the first 30 E0 voxels (10 workers), lenient and C++ rule. Counts only, so contended running is fine.
# Run from icenine_py/:  nohup bash benchmarks/eval_reduction/C/run_step1.sh > benchmarks/eval_reduction/C/cache/logs/step1.log 2>&1 &
set -u
D=benchmarks/eval_reduction/C
for rule in lenient cpp; do
  for v in 15901 17034 18558 19889 22867 5242; do echo "$v $rule"; done
done | xargs -P 6 -L 1 bash -c 'f=benchmarks/eval_reduction/C/cache/step1_raw/full_$1_v$0.npz; [ -f $f ] || uv run python scripts/eval_reduction/candidate_overlap.py run-full --voxel $0 --rule $1'
for rule in lenient cpp; do
  uv run python scripts/eval_reduction/candidate_overlap.py run-e0 --rule $rule --workers 10
done
touch $D/cache/logs/step1.done
