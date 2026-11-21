#!/bin/bash
# AR5 step 19, phase A step 5: the repair on the five stored maps and on J5L (gate g1, g2), two
# chains side by side (at most 2 heavy processes). Run from icenine_py/:
#   nohup scripts/autoresearch/ar5/repair_gate_runs.sh > /dev/null 2>&1 &
# Logs and markers: <cache>/logs/gate_<sample>_<arm>.log and .done ("exit=<code>").
set -u
cd "$(dirname "$0")/../../.." || exit 2
cache="${REPAIR_CACHE:-$PWD/benchmarks/autoresearch/ar5/repair/cache}"
export REPAIR_CACHE="$cache"
unset VIRTUAL_ENV
mkdir -p "$cache/logs"
one() { # one <sample> <arm>
  RESEED_SAMPLE="$1" uv run python scripts/autoresearch/ar5/repair_run.py --arm "$2" \
    > "$cache/logs/gate_$1_$2.log" 2>&1
  echo "exit=$?" > "$cache/logs/gate_$1_$2.done"
}
( one held3 J5; one grad J5; one big J5 ) &
( one held J5; one mid J5; one mid J5L ) &
wait
echo "exit=0" > "$cache/logs/gate_runs.done"
