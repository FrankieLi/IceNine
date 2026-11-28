#!/bin/bash
# C++ time on the pilot: one client, beside the Python timing A/B (PLAN.md, "C++ time on the pilot").
#
#   cd icenine_py && caffeinate -i bash scripts/autoresearch/cpptime/run_cpp_time.sh
#   bash scripts/autoresearch/cpptime/run_cpp_time.sh --dry-run
#       the dry run is real but small: a 30-voxel sub-grid, no pause, the preflight and the quiet
#       verdict are reported and do not stop the script
#
# Order: stack check; C0 (untimed warm-up, also the first identity reference); C1 and C2 (each: 120 s
# pause, preflight, load sampler, /usr/bin/time -p, quiet verdict); identity C1 vs C0 and C2 vs C1;
# closing preflight. Every run is `mpirun -n 2 IceNine p cpp_pilot.config` (rank 0 server, rank 1 the
# one fitting client) under guard.py, which stops it on the plan's stop rules. IceNine exits 0 on
# every failure, so each run is also checked by its files.
#
# Runs in the foreground (no self-detach); the caller backgrounds it. Logs: <cache>/logs/; runs:
# <cache>/runs/<C0|C1|C2>/; the marker <cache>/cpptime.done is written when the script ends, with
# its exit code and reason; <cache>/exit_codes.txt has one line per step.
#   exit 0 all done; 2 a step failed (stop rule of a run: guard code 124 passed 1 hour, 125 a log did
#   not grow for 10 minutes, 126 RUNTIME_ASSERT / [ConfigFile] Error; the .mic.LBFS or a timer line is
#   missing; scoring failed); 3 a preflight refused; 4 passed 2 hours; 5 a timed run was not quiet
#   (R4: neither strict nor relaxed-ok); 6 a repeat differs from its reference (identity, stops at
#   once) or, at the very end, a labelled run's C1/C2 differ by more than 2% in real time (R4b).
#   Verdict code 3 (strict check failed, relaxed ok: a macOS daemon, median load <= 2.5) does not
#   stop: the label is logged and kept in quiet/<run>.json and pair_C.json. C0 is untimed (no pair).
# In every non-zero case the script stops and no wall time is to be reported for the step.
set -u
cd "$(dirname "$0")/../../.." || exit 99   # icenine_py/

DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1
CT=scripts/autoresearch/cpptime
QUIET=scripts/autoresearch/item14/quiet.py   # the one shared quiet.py (R4, R4b)
C=benchmarks/autoresearch/cpptime/cache
[ "$DRY" = 1 ] && C=$C/dry
LIMIT=$(cat $CT/r4_limit.txt)      # R4 load limit, fixed in phase A before any timed run
PAUSE=120
RUN_LIMIT_S=3600
STALL_S=600
MAXS=7200
[ "$DRY" = 1 ] && PAUSE=0
# the binary sits at the root of the main checkout (not in worktrees)
REPO=$(cd .. && pwd)
case "$REPO" in
  */.claude/worktrees/*) REPO=${REPO%%/.claude/worktrees/*} ;;
esac
BIN=$REPO/IceNine
DATA=../Examples/Example2.ManyGrains/ScatteringData_PhaseD
mkdir -p "$C/logs" "$C/runs" "$C/quiet"
MAIN=$C/logs/run_cpp_time.log
CODES=$C/exit_codes.txt
MARKER=$C/cpptime.done
: > "$CODES"
rm -f "$MARKER"
T0=$(date +%s)

say() { echo "$(date '+%F %T') $*" >> "$MAIN"; }

finish() { # rc reason
  say "STOP/END: $2 (exit $1)"
  echo "exit=$1 reason=$2 elapsed_s=$(( $(date +%s) - T0 )) codes=$CODES log=$MAIN" > "$MARKER"
  exit "$1"
}

check_time() { [ $(( $(date +%s) - T0 )) -gt $MAXS ] && finish 4 "passed 2 hours"; return 0; }

step() { # name cmd...   (untimed non-C++ step; a non-zero exit stops the script)
  local n=$1; shift
  say "START $n: $*"
  "$@" > "$C/logs/$n.log" 2>&1
  local rc=$?
  echo "$n rc=$rc" >> "$CODES"
  say "END $n rc=$rc"
  if [ $rc -ne 0 ]; then
    tail -n 8 "$C/logs/$n.log" >> "$MAIN"
    finish 2 "step $n failed (rc $rc); see $C/logs/$n.log"
  fi
  check_time
}

gate() { # name  (preflight; refusal stops the script; report-only in the dry run)
  say "preflight before $1 (load limit of the preflight 1.5)"
  if [ "$DRY" = 1 ]; then
    uv run python ../scripts/dev/timing_preflight.py --json "$C/preflight_$1.json" \
      > "$C/logs/preflight_$1.log" 2>&1
    say "DRY: preflight (report only) rc=$?"
    return 0
  fi
  uv run python ../scripts/dev/timing_preflight.py --require --json "$C/preflight_$1.json" \
    > "$C/logs/preflight_$1.log" 2>&1
  local rc=$?
  echo "preflight_$1 rc=$rc" >> "$CODES"
  if [ $rc -ne 0 ]; then
    tail -n 6 "$C/logs/preflight_$1.log" >> "$MAIN"
    finish 3 "preflight refused before $1; no time is reported"
  fi
}

check_files() { # name  (IceNine exits 0 on every failure: look at the files)
  local R=$C/runs/$1
  [ -s "$R/cpp_pilot.mic.LBFS" ] || return 11
  grep -q "Total Processing Time" "$R/ParallelReconLog_1.log" 2>/dev/null || return 12
  grep -q "Total Time elapsed" "$R/ParallelReconLog_0.log" 2>/dev/null || return 13
  if grep -l -e "RUNTIME_ASSERT" -e "\[ConfigFile\] Error" "$R"/ParallelReconLog_*.log \
    "$C/logs/$1.log" "$R/stderr.log" > /dev/null 2>&1; then return 14; fi
  return 0
}

run_cpp() { # name timed(0|1)  (one C++ run; a stop-rule trip or a failed check stops the script)
  local n=$1 timed=$2
  local R=$C/runs/$n
  rm -rf "$R"; mkdir -p "$R"
  local RA CFG pid spid rc q
  RA="$(pwd)/$R"; CFG="$(pwd)/$C/cpp_pilot.config"
  if [ "$timed" = 1 ]; then
    [ "$PAUSE" -gt 0 ] && { say "pause ${PAUSE}s before $n"; sleep "$PAUSE"; }
    gate "$n"
  fi
  if [ "$timed" = 1 ]; then say "START timed $n (single-worker, one client)"; else say "START untimed $n"; fi
  uv run python $CT/guard.py --limit-s $RUN_LIMIT_S --stall-s $STALL_S --cwd "$RA" \
    --stderr "$RA/stderr.log" --watch 'ParallelReconLog_*.log' -- \
    /usr/bin/time -p -o "$RA/time.txt" mpirun -n 2 "$BIN" p "$CFG" \
    > "$C/logs/$n.log" 2>&1 &
  pid=$!
  uv run python $QUIET sample --out "$C/quiet/$n.jsonl" --watch-pid $pid --period \
    "$([ "$DRY" = 1 ] && echo 5 || echo 30)" > "$C/logs/$n.sampler.log" 2>&1 &
  spid=$!
  wait $pid
  rc=$?
  wait $spid
  echo "$n guard_rc=$rc" >> "$CODES"
  say "END $n guard_rc=$rc, time: $(tr '\n' ' ' < "$R/time.txt" 2>/dev/null)"
  if [ $rc -ne 0 ]; then
    tail -n 5 "$R/stderr.log" >> "$MAIN" 2>/dev/null
    finish 2 "run $n stopped by the stop rules or failed (guard rc $rc: 124 passed 1 hour, 125 no log grew for 10 min, 126 RUNTIME_ASSERT/[ConfigFile] Error); see $R"
  fi
  check_files "$n"
  rc=$?
  echo "$n check_files=$rc" >> "$CODES"
  [ $rc -ne 0 ] && finish 2 "run $n failed its file checks (code $rc: 11 .mic.LBFS missing, 12 no client timer, 13 no server timer, 14 RUNTIME_ASSERT or [ConfigFile] Error in a log); see $R"
  uv run python $CT/parse_logs.py "$R" --time-file "$R/time.txt" --out "$R/parse.json" \
    > "$C/logs/$n.parse.log" 2>&1
  echo "$n parse rc=$?" >> "$CODES"
  local rows=()
  [ "$DRY" = 1 ] && rows=(--rows "$C/sub_rows.npy")
  uv run python $CT/score_cpp.py "$R" ${rows[@]+"${rows[@]}"} --out "$R/score.json" \
    > "$C/logs/$n.score.log" 2>&1
  rc=$?
  echo "$n score rc=$rc" >> "$CODES"
  [ $rc -ne 0 ] && finish 2 "scoring of $n failed (rc $rc); see $C/logs/$n.score.log"
  if [ "$timed" = 1 ]; then
    uv run python $QUIET verdict --samples "$C/quiet/$n.jsonl" --limit "$LIMIT" \
      --out "$C/quiet/$n.json" >> "$MAIN" 2>&1
    q=$?
    echo "quiet_$n rc=$q" >> "$CODES"
    cp "$R/time.txt" "$C/logs/$n.time"   # quiet.py pair reads <cache>/logs/<run>.time
    if [ $q -eq 3 ]; then  # R4b (D-20261114-N1): strict check failed, relaxed check passed: go on, label it
      say "LABEL $n: $(grep -o '"label": "[^"]*"' "$C/quiet/$n.json" | cut -d'"' -f4)"
      q=0
    fi
    if [ $q -ne 0 ] && [ "$DRY" = 1 ]; then
      say "DRY: run $n was not quiet; a real run would stop here (exit 5)"; q=0
    fi
    [ $q -ne 0 ] && finish 5 "run $n was not quiet (limit $LIMIT, see $C/quiet/$n.json); no time is reported"
  fi
  check_time
}

same() { # name dirA dirB   (identity check; a difference stops the script)
  say "START same_$1: $2 vs $3"
  uv run python $CT/same_cpp_run.py "$2" "$3" > "$C/logs/same_$1.log" 2>&1
  local rc=$?
  echo "same_$1 rc=$rc" >> "$CODES"
  say "END same_$1 rc=$rc: $(tail -n 1 "$C/logs/same_$1.log")"
  [ $rc -ne 0 ] && finish 6 "repeat differs (same_$1): $(tail -n 1 "$C/logs/same_$1.log"); reported, no timing claim"
  return 0
}

PAIR_FAIL=""
pair() { # arm run1 run2   (R4b: a labelled run needs its repeat within 2% in real time)
  uv run python $QUIET pair --arm "$1" --runs "$2" "$3" --cache "$C" > "$C/logs/pair_$1.log" 2>&1
  local rc=$?
  echo "pair_$1 rc=$rc" >> "$CODES"
  say "PAIR $1: $(tail -n 1 "$C/logs/pair_$1.log")"
  [ $rc -ne 0 ] && PAIR_FAIL="$PAIR_FAIL $1"
  return 0
}

say "cpptime start (dry=$DRY), R4 load limit $LIMIT, host load $(uptime | sed 's/.*averages: //')"
[ -x "$BIN" ] || finish 99 "binary not found at $BIN"
[ -e "$DATA/full/clean" ] || finish 99 "data not found at $DATA (symlink ScatteringData_PhaseD from the main checkout)"

# --- the config (and in the dry run a 30-voxel sub-grid)
if [ "$DRY" = 1 ]; then
  step subgrid uv run python $CT/make_subgrid.py --n 30 --out "$C/sub.mic" --rows-out "$C/sub_rows.npy"
  step config uv run python $CT/make_config.py --out "$C/cpp_pilot.config" --sample "$C/sub.mic"
else
  step config uv run python $CT/make_config.py --out "$C/cpp_pilot.config"
fi
step stacks_check uv run python $CT/stacks_check.py --out "$C/stacks_check.json"

# --- C0 untimed (warm-up, evaluation count would need lldb: not available here), then C1 and C2
run_cpp C0 0
run_cpp C1 1
same C1_C0 "$C/runs/C1" "$C/runs/C0"
run_cpp C2 1
same C2_C1 "$C/runs/C2" "$C/runs/C1"
pair C C1 C2

# --- closing preflight (report only)
uv run python ../scripts/dev/timing_preflight.py --json "$C/preflight_final.json" \
  > "$C/logs/preflight_final.log" 2>&1
say "closing preflight written"
if [ -n "$PAIR_FAIL" ]; then
  finish 6 "no time reported for C1/C2 (a labelled run, repeat differs by more than 2%; see $C/pair_C.json)"
fi
finish 0 "all steps done"
