#!/bin/zsh
# GO4 (PLAN.md): robustness on rand and corr3, grain-centre start, true grain ids
run() { uv run python scripts/grain_opt/run_pipeline.py --name "$1" --pipeline "$2" --sample "$3" --start centre --group true > benchmarks/grain_opt/runs/$1.log 2>&1; echo "$1 done $?"; }
for smp in rand corr3; do
  run go4_${smp}_As100 Astrict/100 $smp
  run go4_${smp}_As100_V61 Astrict/100,V/61 $smp
  run go4_${smp}_As100_V61s03 Astrict/100,V/61/0.3 $smp
  run go4_${smp}_V61s03 V/61/0.3 $smp
  run go4_${smp}_Als150_As100_V61s03 Als/150,Astrict/100,V/61/0.3 $smp
done
