#!/bin/zsh
# GO2 arms 5-6, reduced (2 h step cap): after go2_Ask150_As100_V20 (20 grains) ends, the
# intensity-Sinkhorn and LS basin arms on the 10-grain subset (2 per level, inside G20), plus the
# strict reference on the same subset.
while pgrep -f "run_pipeline.py --name go2_Ask150_As100_V20" >/dev/null; do sleep 20; done
run() { uv run python scripts/grain_opt/run_pipeline.py --name "$1" --pipeline "$2" --grains 10 > benchmarks/grain_opt/runs/$1.log 2>&1; echo "$1 done"; }
run go2s_Aski150_As100_V20 Asinki/150,Astrict/100,V/20
run go2s_Als150_As100_V20 Als/150,Astrict/100,V/20
run go2s_As100_V20 Astrict/100,V/20
