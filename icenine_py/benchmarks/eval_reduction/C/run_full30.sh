#!/bin/bash
# Task C step 3, extension: the 24 other voxels of the 30-voxel ray-study set (benchmarks/phase_d_seed_diag/ray_overlap/
# voxels.json; its first six are the ledger voxels, whose single-worker runs are reused) for the baseline and the arms with
# a saving on the six, 8 processes at a time (contended; counts and errors only, no timing claim).
# Run from icenine_py/:  nohup bash benchmarks/eval_reduction/C/run_full30.sh > benchmarks/eval_reduction/C/cache/logs/full30.log 2>&1 &
set -u
D=benchmarks/eval_reduction/C
mkdir -p $D/full30/runs
uv run python ../scripts/dev/timing_preflight.py --json $D/full30/preflight_start.json || true
VOX=$(uv run python -c "import json;print(' '.join(str(v) for v in json.load(open('benchmarks/phase_d_seed_diag/ray_overlap/voxels.json'))[6:]))")
for arm in base c0a_k1 c0a_k2 c0b_k2 c1_t05 c3; do
  spec=$(grep -v '^#' $D/arms.txt | awk -v a=$arm '$1==a {print $2}')
  for v in $VOX; do echo "$arm $spec $v"; done
done | xargs -P 8 -L 1 bash -c '[ -f benchmarks/eval_reduction/C/full30/runs/$0_v$2.json ] || uv run python scripts/eval_reduction/eval_ledger.py run --voxel $2 --tag $0 --opt "candidate_dedup=$1" --dir benchmarks/eval_reduction/C/full30'
touch $D/cache/logs/full30.done
