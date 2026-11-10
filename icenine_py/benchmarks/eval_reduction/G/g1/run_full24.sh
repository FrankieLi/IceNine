#!/bin/bash
# Task G1 step 2b: the other 24 voxels of the 30-voxel full-sample set, strict search and the two
# carried SmoothCostRank arms, lenient and C++ rule; one process per run (6 heavy processes).
# From icenine_py/.
L=benchmarks/eval_reduction/G/cache/logs
P=scripts/eval_reduction/g1_full24.py
mkdir -p $L
for rule in old cpp; do
  nohup caffeinate -i uv run python $P run --rule $rule --spec off --tag base \
      > $L/g1_${rule}_base.log 2>&1 &
  nohup caffeinate -i uv run python $P run --rule $rule --spec "E3+E4,tau=0.25" --tag e34t25 \
      > $L/g1_${rule}_e34t25.log 2>&1 &
  nohup caffeinate -i uv run python $P run --rule $rule --spec "E1+E2+E3+E4" --tag e1234a \
      > $L/g1_${rule}_e1234a.log 2>&1 &
done
wait
touch $L/g1_full24.done
