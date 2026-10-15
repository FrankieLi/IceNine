#!/bin/bash
# Task E measurement 3, arm E1+E2+E3+E4 (level widths "auto"): E0 200 lenient and C++ rule after the
# E3+E4 E0 lenient run, and the six full-sample voxels (lenient) after the E0 level-0 landscape.
# From icenine_py/.
L=benchmarks/eval_reduction/E/cache/logs
P=scripts/eval_reduction/smooth_rank_probe.py
S="E1+E2+E3+E4"
(
  while [ ! -f $L/landscape_e0_L0.done ]; do sleep 30; done
  uv run python $P run --kind full --rule old --spec "$S" --tag e1234a --voxels six --workers 2 \
      > $L/e2e_full_old_e1234a.log 2>&1
) &
while [ ! -f $L/rank_e0_old_e34t25.done ]; do sleep 30; done
uv run python $P run --kind e0 --rule old --spec "$S" --tag e1234a --workers 2 > $L/e2e_e0_old_e1234a.log 2>&1
uv run python $P run --kind e0 --rule cpp --spec "$S" --tag e1234a --workers 2 > $L/e2e_e0_cpp_e1234a.log 2>&1
wait
