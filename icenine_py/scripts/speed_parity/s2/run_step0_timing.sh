#!/bin/bash
# S2 Step 0 timings, run alone and one after the other: 0c (s2_timing.py, 5 interleaved rounds) and
# 0a (s2_search.py: ThreeVoxels voxel 0, 3 repeats; 10 E0 voxels, 2 repeats). Before each step the
# machine must be quiet (scripts/dev/timing_preflight.py --require); the runner waits up to
# S2_WAIT_MIN minutes (default 30) for it; every script runs with --allow-busy and its own preflight
# labels its JSON "single-worker" or "contended" (D-20261014-I16), saved beside its results.
#
# Usage (from icenine_py/):
#   nohup caffeinate -i bash scripts/speed_parity/s2/run_step0_timing.sh > /dev/null 2>&1 &
# Log: scripts/speed_parity/cache/logs/s2_timing.log, marker s2_timing.done.
cd "$(dirname "$0")/../../.." || exit 1
OUT=${S2_OUT:-benchmarks/s2_step0}
LOGD=scripts/speed_parity/cache/logs
WAIT_MIN=${S2_WAIT_MIN:-30}
mkdir -p "$OUT" "$LOGD"
rm -f "$LOGD/s2_timing.done"

wait_quiet() {
  local i=0
  while [ "$i" -lt "$WAIT_MIN" ]; do
    if uv run python ../scripts/dev/timing_preflight.py --require > /dev/null 2>&1; then
      echo "quiet"
      return 0
    fi
    i=$((i + 1))
    sleep 60
  done
  echo "not-quiet"
}

{
  echo "start $(date) $(uptime | sed 's/.*load/load/')"
  flag=$(wait_quiet)
  echo "0c timing (preflight wait: $flag) $(date)"
  uv run python scripts/speed_parity/s2/s2_timing.py --out "$OUT" --allow-busy
  flag=$(wait_quiet)
  echo "0a search three (preflight wait: $flag) $(date)"
  uv run python scripts/speed_parity/s2/s2_search.py --set three --repeats 3 --out "$OUT" \
    --allow-busy
  flag=$(wait_quiet)
  echo "0a search e0 (preflight wait: $flag) $(date)"
  uv run python scripts/speed_parity/s2/s2_search.py --set e0 --n-voxels 10 --repeats 2 \
    --out "$OUT" --allow-busy
  echo "end $(date) $(uptime | sed 's/.*load/load/')"
} > "$LOGD/s2_timing.log" 2>&1
touch "$LOGD/s2_timing.done"
