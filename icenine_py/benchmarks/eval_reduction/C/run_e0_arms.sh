#!/bin/bash
# Task C step 3: E0 clean (200 voxels, per-voxel images, seed 0, mc arm, min_sin_eta 0), 10 workers per arm, one arm
# after the other. Arms and specs from arms.txt (or the arms named on the command line). Wall times are contended
# multi-worker and not a timing claim. Run from icenine_py/:
#   nohup bash benchmarks/eval_reduction/C/run_e0_arms.sh [arm ...] > benchmarks/eval_reduction/C/cache/logs/e0_arms.log 2>&1 &
set -u
D=benchmarks/eval_reduction/C
mkdir -p $D/cache/logs $D/step3
uv run python ../scripts/dev/timing_preflight.py --json $D/step3/preflight_e0.json || true
want=" $* "
grep -v '^#' $D/arms.txt | while read arm spec; do
  [ -z "$arm" ] && continue
  if [ "$want" != "  " ] && [[ "$want" != *" $arm "* ]]; then continue; fi
  uv run python scripts/eval_reduction/candidate_arms.py run --arm $arm --spec "$spec" --workers 10
done
touch $D/cache/logs/e0_arms.done
