#!/bin/bash
# S3 Step 4 timings, run in sequence (single worker, nothing else of ours running): ThreeVoxels
# voxel 0, the first 10 E0 clean voxels, and the MPS crossover. Contended label when the machine is
# not quiet (--allow-busy). Log: scripts/speed_parity/cache/logs/s3_timing.log, marker s3_timing.done.
# Usage (from icenine_py/): nohup bash scripts/speed_parity/run_s3_timing.sh &
cd "$(dirname "$0")/../.." || exit 1
OUT=benchmarks/s3_steps1_5
LOGD=scripts/speed_parity/cache/logs
mkdir -p "$OUT" "$LOGD"
rm -f "$LOGD/s3_timing.done"
{
  uv run python scripts/speed_parity/s3_timing.py search --set three --seeds 3 --repeats 2 \
    --out "$OUT" --allow-busy
  uv run python scripts/speed_parity/s3_timing.py search --set e0 --n-voxels 10 --seeds 1 \
    --repeats 2 --out "$OUT" --allow-busy
  uv run python scripts/speed_parity/s3_timing.py crossover --out "$OUT" --allow-busy
} > "$LOGD/s3_timing.log" 2>&1
touch "$LOGD/s3_timing.done"
