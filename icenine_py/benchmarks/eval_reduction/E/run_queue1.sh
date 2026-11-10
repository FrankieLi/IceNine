#!/bin/bash
# Task E: baseline probe runs for measurement 2, chained behind the jobs already running
# (at most about 6 processes at once). Run from icenine_py/.
L=benchmarks/eval_reduction/E/cache/logs
P=scripts/eval_reduction/smooth_rank_probe.py
while [ ! -f $L/rank_e0_cpp_e0base.flag ] && [ ! -f $L/rank_e0_cpp_base.done ]; do sleep 30; done
uv run python $P run --kind full --rule old --spec off --tag base --workers 2 > $L/rank_full_old_base.log 2>&1
uv run python $P run --kind full --rule cpp --spec off --tag base --workers 2 > $L/rank_full_cpp_base.log 2>&1
