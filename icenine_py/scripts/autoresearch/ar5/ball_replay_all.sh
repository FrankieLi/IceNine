#!/bin/bash
# AR5 step 23, phase A steps 4 and 5: the replay runs, ONE AT A TIME (owner rule of 16 Nov 2026).
#
#   cd icenine_py && nohup caffeinate -i scripts/autoresearch/ar5/ball_replay_all.sh --set g0 \
#       > /dev/null 2>&1 &
#   --set g0     mid with rule CB off, then mid with CB on (gate g0)
#   --set rest   CB on, the other seven samples (held5 held4 held3 held grad big held6)
#   --dry-run    the commands only
#
# Each job writes <cache>/logs/replay_<mode>_<s>.log and .rc; the set writes
# <cache>/logs/replay_<set>.done ("exit=<code>"). A job whose summary.json and replay.done exist is
# kept (resume skips it); an interrupted job resumes from its rows.jsonl. Before each job: no other
# heavy run (IceNine, pilot_bfs, guard.py, mpirun, a ball_ or defer_ run) is active, else it waits
# (it never stops anything).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${BALL_CACHE:-$here/benchmarks/autoresearch/ar5/ball/cache}"
set_name=""
dry=0
while [ $# -gt 0 ]; do
  case "$1" in
    --set) shift; set_name="$1" ;;
    --dry-run) dry=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
case "$set_name" in
  g0) jobs="off:mid on:mid" ;;
  rest) jobs="on:held5 on:held4 on:held3 on:held on:grad on:big on:held6" ;;
  *) echo "--set g0|rest is required" >&2; exit 2 ;;
esac
export SCREEN_CACHE="$cache" DEFER_CACHE="$cache" BALL_CACHE="$cache"
unset VIRTUAL_ENV
logs="$cache/logs"
mkdir -p "$logs"
busy() { # a heavy run of anyone else (the C++ binary, pilot_bfs, guard.py, mpirun)
  ps -axo comm,args | grep -E '^(/[^ ]*/)?IceNine |pilot_bfs|guard\.py|mpirun' | grep -v grep > /dev/null
}
rc_all=0
for j in $jobs; do
  mode="${j%%:*}"; s="${j#*:}"
  cmd="RESEED_SAMPLE=$s uv run python scripts/autoresearch/ar5/ball_replay.py --mode $mode --resume"
  if [ "$dry" = "1" ]; then echo "$cmd > $logs/replay_${mode}_$s.log"; continue; fi
  if [ -e "$cache/replay_${mode}_$s/replay.done" ]; then echo "kept $mode $s"; continue; fi
  while busy; do sleep 60; done
  eval "$cmd" >> "$logs/replay_${mode}_$s.log" 2>&1
  rc=$?
  echo "$rc" > "$logs/replay_${mode}_$s.rc"
  if [ "$rc" != "0" ]; then rc_all=1; break; fi
done
if [ "$dry" = "0" ]; then
  { echo "exit=$rc_all"; date +%Y-%m-%dT%H:%M:%S%z; } > "$logs/replay_$set_name.done"
fi
exit $rc_all
