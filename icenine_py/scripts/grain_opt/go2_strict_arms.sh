#!/bin/zsh
# GO2 arms on the strict cost (PLAN.md GO2, arms 1 to 4); run from icenine_py/
set -e
run() { uv run python scripts/grain_opt/run_pipeline.py --name "$1" --pipeline "$2" > benchmarks/grain_opt/runs/$1.log 2>&1; echo "$1 done"; }
mkdir -p benchmarks/grain_opt/runs
run go2_As50 Astrict/50
run go2_As100 Astrict/100
run go2_As200 Astrict/200
run go2_As400 Astrict/400
run go2_As200o2 Astrict/200/o2
run go2_As100_V10 Astrict/100,V/10
run go2_As100_V30 Astrict/100,V/30
run go2_As100_V61 Astrict/100,V/61
run go2_V61 V/61
run go2_V20 V/20
