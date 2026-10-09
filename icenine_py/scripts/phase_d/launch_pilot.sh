#!/bin/bash
# Launch the five Phase D pilot arms concurrently (one process each).
# usage (from anywhere): scripts/phase_d/launch_pilot.sh <outdir> <tag> [extra pilot_bfs.py args]
# Logs: scripts/phase_d/cache/pilot/logs/<tag>_<arm>_<variant>.{log,done}.
# NOTE: the pilot of 2026-10-08 was launched with an earlier version of this script that did NOT
# set the thread variables below (torch default of 8 threads); pass --timing-label with the real
# number of concurrent jobs.
cd "$(dirname "$0")/../.." || exit 1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
OUT=$1; TAG=$2; shift 2
mkdir -p scripts/phase_d/cache/pilot/logs
while read -r A V; do
  (nohup uv run python scripts/phase_d/pilot_bfs.py --arm "$A" --variant "$V" --out "$OUT" "$@" \
     > "scripts/phase_d/cache/pilot/logs/${TAG}_${A}_${V}.log" 2>&1
   touch "scripts/phase_d/cache/pilot/logs/${TAG}_${A}_${V}.done") &
done <<LIST
mc clean
cma clean
cma_noretry clean
mc realistic_q16
cma realistic_q16
LIST
