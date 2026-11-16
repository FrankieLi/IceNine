#!/bin/zsh
# GO4c (PLAN.md): larger ladders for 2-5 deg (I37); the gated global search (I36) is blocked on cost
run() { uv run python scripts/grain_opt/run_pipeline.py --name "$1" --pipeline "$2" --sample "$3" --start centre --group true > benchmarks/grain_opt/runs/$1.log 2>&1; echo "$1 done $?"; }
for smp in corr3 rand; do
  run go4c_${smp}_L244x3 L/244/3 $smp
  run go4c_${smp}_L366x3 L/366/3 $smp
  run go4c_${smp}_L183_N30 L/183,N/30/0.3 $smp
done
