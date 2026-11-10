#!/bin/bash
# S2 Step 5 timings on the library path, run alone and one after the other: single-orientation and
# batched us per call (s2_lib_timing.py, 5 interleaved rounds), then the no-start search split
# (s2_search.py --lib: ThreeVoxels voxel 0, 3 repeats; 10 E0 voxels, 2 repeats). Before each step
# the machine must be quiet (scripts/dev/timing_preflight.py --require); the runner waits up to
# S2_WAIT_MIN minutes (default 30) for it; every script runs with --allow-busy and its own preflight
# labels its JSON "single-worker" or "contended" (D-20261014-I16), saved beside its results.
#
# Usage (from icenine_py/):
#   nohup caffeinate -i bash scripts/speed_parity/s2/run_results_timing.sh > /dev/null 2>&1 &
# Log: scripts/speed_parity/cache/logs/s2_results_timing.log, marker s2_results_timing.done.
cd "$(dirname "$0")/../../.." || exit 1
OUT=${S2_OUT:-benchmarks/s2_results}
LOGD=scripts/speed_parity/cache/logs
WAIT_MIN=${S2_WAIT_MIN:-30}
mkdir -p "$OUT" "$LOGD"
rm -f "$LOGD/s2_results_timing.done"

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
  echo "single (preflight wait: $flag) $(date)"
  uv run python scripts/speed_parity/s2/s2_lib_timing.py single --out "$OUT" --allow-busy
  flag=$(wait_quiet)
  echo "batched (preflight wait: $flag) $(date)"
  uv run python scripts/speed_parity/s2/s2_lib_timing.py batched --out "$OUT" --allow-busy
  flag=$(wait_quiet)
  echo "search three (preflight wait: $flag) $(date)"
  uv run python scripts/speed_parity/s2/s2_search.py --set three --repeats 3 --lib --out "$OUT" \
    --allow-busy
  flag=$(wait_quiet)
  echo "search e0 (preflight wait: $flag) $(date)"
  uv run python scripts/speed_parity/s2/s2_search.py --set e0 --n-voxels 10 --repeats 2 --lib \
    --out "$OUT" --allow-busy
  echo "end $(date) $(uptime | sed 's/.*load/load/')"
} > "$LOGD/s2_results_timing.log" 2>&1
touch "$LOGD/s2_results_timing.done"
