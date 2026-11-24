#!/bin/bash
# AR5 step 23, phase A steps 4 to 6 as one detached chain, stopping at the first failed gate (the
# plan's stop rule: no tuning, no further run). ONE JOB AT A TIME (each script below runs one job at
# a time and waits for any heavy run of anyone else).
#
#   cd icenine_py && nohup caffeinate -i scripts/autoresearch/ar5/ball_chain.sh > /dev/null 2>&1 &
#
# 1. wait for logs/replay_g0.done (mid with CB off, then on; started by ball_replay_all.sh --set g0)
# 2. g0 judged (ball_replay_collect.py g0)       -> ar5/ball/g0.json
# 3. the replay of the other seven samples (ball_replay_all.sh --set rest)
# 4. g1 judged (ball_replay_collect.py g1)       -> ar5/ball/replay.json, replay_tables.md
# 5. live J5DC with the boundary pass and the screen arm on mid and held (ball_run_all.sh --stage gate)
# 6. g2 judged (ball_gate.py)                    -> ar5/ball/gate.json
# It writes logs/chain.log and, at its end, logs/chain.done ("exit=<code>", "reason=<text>"). It
# NEVER writes logs/gate_passed: the implementer does that after reading gate.json.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${BALL_CACHE:-$here/benchmarks/autoresearch/ar5/ball/cache}"
logs="$cache/logs"
mkdir -p "$logs"
export SCREEN_CACHE="$cache" DEFER_CACHE="$cache" BALL_CACHE="$cache"
unset VIRTUAL_ENV
exec >> "$logs/chain.log" 2>&1
fin() { echo "exit=$1" > "$logs/chain.done"; echo "reason=$2" >> "$logs/chain.done"; date +%Y-%m-%dT%H:%M:%S%z >> "$logs/chain.done"; echo "chain end: exit=$1 $2"; exit "$1"; }
echo "chain start $(date +%Y-%m-%dT%H:%M:%S%z)"
while [ ! -e "$logs/replay_g0.done" ]; do sleep 30; done
[ "$(head -1 "$logs/replay_g0.done")" = "exit=0" ] || fin 1 "the g0 replay jobs failed (see replay_*_mid.log)"
uv run python scripts/autoresearch/ar5/ball_replay_collect.py g0 || fin 1 "g0 FAILED (stop rule: no tuning, no held7)"
scripts/autoresearch/ar5/ball_replay_all.sh --set rest
[ "$(head -1 "$logs/replay_rest.done")" = "exit=0" ] || fin 1 "a replay job of the rest set failed (see replay_on_*.log)"
uv run python scripts/autoresearch/ar5/ball_replay_collect.py g1 || fin 1 "g1 FAILED or a check or first-stage mismatch (stop rule: no tuning, no held7)"
scripts/autoresearch/ar5/ball_run_all.sh --stage gate --resume
[ "$(head -1 "$logs/ball_gate.done")" = "exit=0" ] || fin 1 "the gate runs failed (see ball_gate.done and the logs)"
uv run python scripts/autoresearch/ar5/ball_gate.py || fin 1 "g2 FAILED (stop rule: no tuning, no held7)"
fin 0 "g0, g1 and g2 passed; the gate's verdict is in ar5/ball/gate.json (gate_passed not written)"
