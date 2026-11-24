#!/bin/bash
# AR5 step 23, information only: after the replay set "rest" ends, run ball_where.py on each of the
# eight samples (the lost steps again, one sample at a time), then ball_info.py. Waits while any
# other heavy run is active (never stops one). Marker logs/info.done ("exit=<code>").
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${BALL_CACHE:-$here/benchmarks/autoresearch/ar5/ball/cache}"
logs="$cache/logs"
export SCREEN_CACHE="$cache" DEFER_CACHE="$cache" BALL_CACHE="$cache"
unset VIRTUAL_ENV
exec >> "$logs/info_chain.log" 2>&1
busy() { ps -axo comm,args | grep -E '^(/[^ ]*/)?IceNine |pilot_bfs|guard\.py|mpirun' | grep -v grep > /dev/null; }
while [ ! -e "$logs/replay_rest.done" ]; do sleep 60; done
[ "$(head -1 "$logs/replay_rest.done")" = "exit=0" ] || { echo "exit=1" > "$logs/info.done"; exit 1; }
rc=0
for s in mid held5 held4 held3 held grad big held6; do
  [ -e "$cache/where_$s.json" ] && continue
  while busy; do sleep 60; done
  RESEED_SAMPLE=$s uv run python scripts/autoresearch/ar5/ball_where.py || { rc=1; break; }
done
[ "$rc" = "0" ] && { uv run python scripts/autoresearch/ar5/ball_info.py || rc=1; }
echo "exit=$rc" > "$logs/info.done"
exit $rc
