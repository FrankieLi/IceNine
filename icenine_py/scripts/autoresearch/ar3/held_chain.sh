#!/bin/bash
# AR3 step 11: the one S8 baseline on held, then its boundary pass (skipped if the run is capped).
#   scripts/autoresearch/ar3/held_chain.sh [cap_min] [log_name]     (run from icenine_py/)
# cap_min: the run's cap in minutes (default 120, the first run); log_name: default chain_held.
# The second, finished run (D-20261111-N17) is `held_chain.sh 480 chain_held_480`, started after
# the first run's folder cache/S8_held was moved to cache/S8_held_cap120.
# Log: benchmarks/autoresearch/ar3/held/cache/logs/<log_name>.log; marker <log_name>.done;
# <log_name>.commit holds the commit and the start time.
set -u
cap="${1:-120}"
name="${2:-chain_held}"
case "$cap" in '' | *[!0-9]*) echo "cap_min must be a whole number of minutes" >&2; exit 2 ;; esac
cache="$PWD/benchmarks/autoresearch/ar3/held/cache"
mkdir -p "$cache/logs"
{ git rev-parse --short HEAD; date "+%Y-%m-%d %H:%M:%S"; } > "$cache/logs/$name.commit"
(
  export RESEED_SAMPLE=held RESEED_CACHE="$cache"
  nohup caffeinate -i bash -c '
    uv run python scripts/autoresearch/ar5/fit_run.py --name S8 --cap-min '"$cap"' &&
    test ! -e "'"$cache"'/S8_held/capped.json" &&
    uv run python scripts/autoresearch/ar5/fit_pass.py --name S8' \
    > "$cache/logs/$name.log" 2>&1
  echo "exit $?" > "$cache/logs/$name.done"
) > /dev/null 2>&1 &
echo "started $name (cap $cap minutes)"
