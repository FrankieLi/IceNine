#!/bin/bash
# AR5 step 21, phase A step 4: the gate runs on the six disclosed samples, three chains side by side
# (at most 3 heavy processes). Per sample, in order: the full-budget re-run of J5's pass tries with
# the gain trace (screen_pass.py --mode full), the full-tries replay of J5C (screen_run.py --arm
# J5C --source full) and arm J5S (screen_run.py --arm J5S). Run from icenine_py/ after
# ``screen_link_cache.sh`` and the committed freeze:
#   nohup caffeinate -i scripts/autoresearch/ar5/screen_gate_runs.sh > /dev/null 2>&1 &
# Logs and markers: <cache>/logs/gate_<sample>_<step>.log and gate_<sample>.done ("exit=<code>"),
# gate_runs.done at the end. Then: uv run python scripts/autoresearch/ar5/screen_gate.py
set -u
cd "$(dirname "$0")/../../.." || exit 2
cache="${SCREEN_CACHE:-$PWD/benchmarks/autoresearch/ar5/screen/cache}"
export SCREEN_CACHE="$cache"
unset VIRTUAL_ENV
mkdir -p "$cache/logs"
one() { # one <sample>
  s="$1"
  rc=0
  RESEED_SAMPLE="$s" uv run python scripts/autoresearch/ar5/screen_pass.py --mode full \
    > "$cache/logs/gate_${s}_full.log" 2>&1 || rc=$?
  if [ "$rc" = "0" ]; then
    RESEED_SAMPLE="$s" uv run python scripts/autoresearch/ar5/screen_run.py --arm J5C --source full \
      > "$cache/logs/gate_${s}_J5C.log" 2>&1 || rc=$?
  fi
  if [ "$rc" = "0" ]; then
    RESEED_SAMPLE="$s" uv run python scripts/autoresearch/ar5/screen_run.py --arm J5S \
      > "$cache/logs/gate_${s}_J5S.log" 2>&1 || rc=$?
  fi
  echo "exit=$rc" > "$cache/logs/gate_${s}.done"
}
( one held4; one big ) &
( one held3; one mid ) &
( one held; one grad ) &
wait
echo "exit=0" > "$cache/logs/gate_runs.done"
