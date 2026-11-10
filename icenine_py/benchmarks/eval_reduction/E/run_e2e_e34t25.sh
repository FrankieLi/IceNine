#!/bin/bash
# Task E measurement 3, arm E3+E4 (tau 0.25 deg): E0 200 lenient and C++ rule, then the six
# full-sample voxels, lenient rule. From icenine_py/.
L=benchmarks/eval_reduction/E/cache/logs
P=scripts/eval_reduction/smooth_rank_probe.py
S="E3+E4,tau=0.25"
uv run python $P run --kind e0 --rule old --spec "$S" --tag e34t25 --workers 2 > $L/e2e_e0_old_e34t25.log 2>&1
uv run python $P run --kind e0 --rule cpp --spec "$S" --tag e34t25 --workers 2 > $L/e2e_e0_cpp_e34t25.log 2>&1
uv run python $P run --kind full --rule old --spec "$S" --tag e34t25 --voxels six --workers 2 > $L/e2e_full_old_e34t25.log 2>&1
