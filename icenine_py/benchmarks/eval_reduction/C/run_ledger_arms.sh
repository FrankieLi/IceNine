#!/bin/bash
# Task C step 3, ledger part: the six full-sample voxels (clean images, seed 1) per arm, one process at a time
# (single-worker). Between runs wait (up to 16 s) for the 1-minute load to fall under 1.3 so the preflight inside
# eval_ledger.py has a chance to pass; runs that still fail it are labelled "contended" by the tool.
# Run from icenine_py/:  nohup bash benchmarks/eval_reduction/C/run_ledger_arms.sh > benchmarks/eval_reduction/C/cache/logs/ledger.log 2>&1 &
set -u
D=benchmarks/eval_reduction/C
mkdir -p $D/ledger/runs $D/ledger/logs
uv run python ../scripts/dev/timing_preflight.py --json $D/ledger/preflight_start.json || true
for arm in base c0a_k1 c0a_k2 c0b_k1 c0b_k2 c1_t025 c1_t05 c2_t025 c2_t05 c3 c1_t05_qf c3_qf c0bL1_k1 c0bL1_k2; do
  spec=$(grep -v '^#' $D/arms.txt | awk -v a=$arm '$1==a {print $2}')
  for v in 15901 17034 18558 19889 22867 5242; do
    [ -f $D/ledger/runs/${arm}_v$v.json ] && continue
    for i in 1 2; do
      l=$(sysctl -n vm.loadavg | awk '{print $2}')
      awk -v l=$l 'BEGIN{exit !(l<1.3)}' && break
      sleep 8
    done
    uv run python scripts/eval_reduction/eval_ledger.py run --voxel $v --tag $arm --opt "candidate_dedup=$spec" --dir $D/ledger
  done
done
touch $D/ledger/logs/ledger.done
