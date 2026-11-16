#!/bin/zsh
# GO2 arms 5 and 6 (PLAN.md GO2): Sinkhorn / LS as basin finder before the strict cost
set -e
run() { uv run python scripts/grain_opt/run_pipeline.py --name "$1" --pipeline "$2" > benchmarks/grain_opt/runs/$1.log 2>&1; echo "$1 done"; }
run go2_Ask150_As100_V20 Asink/150,Astrict/100,V/20
run go2_Aski150_As100_V20 Asinki/150,Astrict/100,V/20
run go2_Als150_As100_V20 Als/150,Astrict/100,V/20
run go2_Ask150_Bsk40_As50 Asink/150,Bsink/40,Astrict/50
run go2_Aski150 Asinki/150
