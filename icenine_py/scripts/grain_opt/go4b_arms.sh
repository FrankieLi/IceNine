#!/bin/zsh
# GO4b (PLAN.md): far starts on rand and corr3, grain-centre start, true grain ids
run() { uv run python scripts/grain_opt/run_pipeline.py --name "$1" --pipeline "$2" --sample "$3" --start centre --group true > benchmarks/grain_opt/runs/$1.log 2>&1; echo "$1 done $?"; }
for smp in corr3 rand; do
  run go4b_${smp}_L61 L/61 $smp
  run go4b_${smp}_L183 L/183 $smp
  run go4b_${smp}_N30 N/30/0.3 $smp
  run go4b_${smp}_L61_N30 L/61,N/30/0.3 $smp
  run go4b_${smp}_Als_As_L61_N30 Als/150,Astrict/100,L/61,N/30/0.3 $smp
done
