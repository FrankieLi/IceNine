#!/bin/bash
# fine1: the slice runs and the whole-layer runs, one arm after another (never two at once).
#
#   cd icenine_py
#   nohup caffeinate -i bash scripts/autoresearch/fine1/fine1_launch.sh --stage slice \
#       --cpp-bin <C++ binary for CS> > /dev/null 2>&1 &
#   nohup caffeinate -i bash scripts/autoresearch/fine1/fine1_launch.sh --stage whole \
#       --cpp-bin <binary> --cap-h CS=<h>,J5S=<h>,PS=<h> --resume > /dev/null 2>&1 &
#   bash scripts/autoresearch/fine1/fine1_launch.sh --stage slice --cpp-bin <binary> --dry-run
#
# Options
#   --stage slice|whole   slice: the 2,000 voxels nearest the centre (fine7_slice2000.mic);
#                         whole: the 98,280-voxel layer
#   --arms CS,J5S,PS      order of the runs (default CS,J5S,PS, then CF if --cf-bin is given)
#   --cpp-bin PATH        the C++ binary of arm CS (required for CS)
#   --cf-bin PATH         the fixed C++ binary of arm CF (an extra labelled run; CF is added)
#   --ranks N             MPI ranks of the parallel arms CSP and CFP (default 13 = one server and 12
#                         fitting clients; CS and CF always use 2); CSP and CFP are run only when
#                         named in --arms (CSP: the CS binary and config, CFP: the CF ones)
#   --cap-h A=h,B=h       wall-clock cap per arm in hours (the plan: 1.5 x the arm's high
#                         projection); required for --stage whole; slice default 1 h each
#   --pause S             seconds of rest before each timed run (default 120)
#   --allow-contended     start a run even if the timing preflight refuses (labelled contended);
#                         default: the preflight refusal stops the launcher (exit 3)
#   --continue-contended  after a not-quiet arm (marker exit 5) go on to the next arm; the launcher
#                         then exits 5 at the end (default: stop at once)
#   --resume              skip arms whose <arm>.done marker exists; an arm that was started and has
#                         no marker is moved to interrupted_<time>/ and run again from the start
#                         (neither program can restart half-way: a run is one piece)
#   --dry-run             write the configs under <cache>/dry, run the checks and print every
#                         command; starts no arm, writes no marker
#
# Per arm: freeze check, 120 s rest, timing preflight (--require), then the run under guard.py
# (stop rules: the cap, no growth of the watched log for 30 minutes, RUNTIME_ASSERT or
# "[ConfigFile] Error" in stderr), /usr/bin/time -p and the load sampler of item14/quiet.py; then
# its file checks, the quiet verdict (load limit of cpptime/r4_limit.txt), and the marker
# <cache>/runs/<stage>/<arm>.done ("exit=.. reason=.. elapsed_s=.. label=.."). The next arm starts
# only if the previous one ended with exit 0 and a marker. A not-quiet run is relabelled
# "contended", its marker says so and the launcher stops (exit 5; no rerun).
# The launcher reads no map and no result log of fine7. In the slice stage only (counts and seconds
# are what the slice is for) cpptime/parse_logs.py reads the C++ run's log; in the whole stage
# nothing is parsed or scored until all end markers exist (phase D).
#
# Exit: 0 all arms done; 2 a run failed or hit a stop rule; 3 preflight refused; 5 not quiet;
#       64 bad arguments; 99 a precondition is missing (freeze, binary, data).
set -u
cd "$(dirname "$0")/../../.." || exit 99   # icenine_py/

STAGE=""; ARMS=""; CPP_BIN=""; CF_BIN=""; CAPS=""; RANKS=13; PAUSE=120; ALLOW=0; RESUME=0; DRY=0; CONT=0; NOTQUIET=0
while [ $# -gt 0 ]; do
  case "$1" in
    --stage) STAGE=$2; shift 2;;
    --arms) ARMS=$2; shift 2;;
    --cpp-bin) CPP_BIN=$2; shift 2;;
    --cf-bin) CF_BIN=$2; shift 2;;
    --cap-h) CAPS=$2; shift 2;;
    --pause) PAUSE=$2; shift 2;;
    --ranks) RANKS=$2; shift 2;;
    --allow-contended) ALLOW=1; shift;;
    --continue-contended) CONT=1; shift;;
    --resume) RESUME=1; shift;;
    --dry-run) DRY=1; shift;;
    *) echo "unknown argument $1" >&2; exit 64;;
  esac
done
case "$STAGE" in
  slice) GRID=slice; N_EXPECT=2000;;
  whole) GRID=full; N_EXPECT=98280;;
  *) echo "--stage slice|whole is required" >&2; exit 64;;
esac
if [ -z "$ARMS" ]; then
  ARMS=CS,J5S,PS
  [ -n "$CF_BIN" ] && ARMS=$ARMS,CF
fi
if [ "$STAGE" = whole ] && [ -z "$CAPS" ]; then
  echo "--cap-h is required for --stage whole (1.5 x each arm's high projection)" >&2; exit 64
fi

F1=scripts/autoresearch/fine1
CT=scripts/autoresearch/cpptime
QUIET=scripts/autoresearch/item14/quiet.py
C=benchmarks/autoresearch/fine1/cache
LIMIT=$(cat $CT/r4_limit.txt)
RUNS=$C/runs/$STAGE
[ "$DRY" = 1 ] && RUNS=$C/dry/$STAGE
LOGS=$C/logs
mkdir -p "$RUNS" "$LOGS"
MAIN=$LOGS/launch_${STAGE}$([ "$DRY" = 1 ] && echo _dry).log
T0=$(date +%s)
STALL_S=1800

say() { echo "$(date '+%F %T') $*" >> "$MAIN"; echo "$*"; }

cap_s() { # arm -> seconds
  local arm=$1 h=""
  if [ -n "$CAPS" ]; then
    h=$(echo "$CAPS" | tr ',' '\n' | grep "^$arm=" | cut -d= -f2)
  fi
  [ -z "$h" ] && [ "$STAGE" = slice ] && h=1
  [ -z "$h" ] && { say "no --cap-h for arm $arm"; exit 64; }
  awk -v h="$h" 'BEGIN { printf "%d", h * 3600 }'
}

ARM_T0=$T0
marker() { # arm rc reason label
  local arm=$1
  echo "exit=$2 reason=$3 elapsed_s=$(( $(date +%s) - ARM_T0 )) label=$4 run=$RUNS/$arm log=$LOGS/${STAGE}_${arm}.log" > "$RUNS/$arm.done"
}

# ---- preconditions
say "fine1 launch: stage $STAGE grid $GRID arms $ARMS dry=$DRY resume=$RESUME continue-contended $([ "$CONT" = 1 ] && echo on || echo off) host load $(uptime | sed 's/.*averages: //')"
uv run python $F1/fine1_freeze.py --check --full $([ "$DRY" = 1 ] || echo --outputs) \
  > "$LOGS/freeze_check_$STAGE.log" 2>&1
RC=$?
if [ $RC -ne 0 ]; then
  say "freeze check failed (rc $RC): $(tail -n 1 "$LOGS/freeze_check_$STAGE.log")"
  [ "$DRY" = 1 ] || exit 99
fi

for arm in $(echo "$ARMS" | tr ',' ' '); do
  case "$arm" in
    CS|CSP) BIN=$CPP_BIN;;
    CF|CFP) BIN=$CF_BIN;;
    *) BIN="";;
  esac
  case "$arm" in
    CS|CSP|CF|CFP) [ -x "$BIN" ] || { say "arm $arm needs an executable binary ($BIN)"; exit 99; };;
  esac
  case "$arm" in CF|CFP) [ -z "$CF_BIN" ] && { say "arm $arm needs --cf-bin"; exit 64; };; esac
done

is_cpp() { case "$1" in CS|CSP|CF|CFP) return 0;; *) return 1;; esac; }
ranks_of() { case "$1" in CSP|CFP) echo "$RANKS";; *) echo 2;; esac; }  # one server + clients

run_arm() { # arm
  local arm=$1 R=$RUNS/$1 BIN="" cap pid spid rc q
  local -a CMD
  case "$arm" in CS|CSP) BIN=$CPP_BIN;; CF|CFP) BIN=$CF_BIN;; esac
  cap=$(cap_s "$arm") || exit 64
  if [ "$RESUME" = 1 ] && [ -f "$RUNS/$arm.done" ]; then
    say "arm $arm: marker exists, skipped ($(cat "$RUNS/$arm.done"))"; return 0
  fi
  if [ -e "$R" ] && [ "$DRY" = 0 ]; then
    if [ "$RESUME" = 1 ]; then
      mkdir -p "$RUNS/interrupted_$(date +%Y%m%d_%H%M%S)"
      mv "$R" "$RUNS/interrupted_$(date +%Y%m%d_%H%M%S)/$arm" && say "arm $arm: partial run moved aside"
      rm -f "$RUNS/$arm.done"
    else
      say "arm $arm: $R exists; use --resume (moves it aside) or remove it"; exit 99
    fi
  fi
  mkdir -p "$R"
  local LOG=$LOGS/${STAGE}_${arm}.log
  local CFG=$(pwd)/$R/arm.config
  case "$arm" in
    CS|CSP|CF|CFP)
      uv run python $F1/fine1_config.py --arm "$arm" --grid $GRID --out "$CFG" > "$LOG" 2>&1 \
        || { say "arm $arm: config not written"; exit 99; };;
  esac
  if is_cpp "$arm"; then
    CMD=(uv run python $CT/guard.py --limit-s "$cap" --stall-s $STALL_S --cwd "$(pwd)/$R"
      --stderr "$(pwd)/$R/stderr.log" --watch 'ParallelReconLog_*.log' --
      /usr/bin/time -p -o "$(pwd)/$R/time.txt" mpirun -n "$(ranks_of "$arm")" "$BIN" p "$CFG")
  else
    CMD=(uv run python $CT/guard.py --limit-s "$cap" --stall-s $STALL_S
      --stderr "$(pwd)/$R/stderr.log" --watch "$(pwd)/$R/progress.jsonl" --
      /usr/bin/time -p -o "$(pwd)/$R/time.txt"
      uv run python $F1/fine1_run_py.py --arm "$arm" --grid $GRID --out "$(pwd)/$R")
  fi
  say "arm $arm: cap ${cap} s, stall ${STALL_S} s; command: ${CMD[*]}"
  if [ "$DRY" = 1 ]; then
    if is_cpp "$arm"; then
      say "DRY: config written to $CFG; binary $BIN sha256 $(shasum -a 256 "$BIN" | cut -d' ' -f1)"
    else
      uv run python $F1/fine1_run_py.py --arm "$arm" --grid $GRID --out "$R" --dry-run >> "$LOG" 2>&1
      say "DRY: $(tail -n 2 "$LOG" | tr '\n' ' ')"
    fi
    uv run python ../scripts/dev/timing_preflight.py --json "$R/preflight.json" > "$R/preflight.log" 2>&1
    say "DRY: preflight (report only) rc=$?"
    return 0
  fi
  uv run python $F1/fine1_freeze.py --check || { marker "$arm" 99 "freeze check failed" none; exit 99; }
  [ "$PAUSE" -gt 0 ] && { say "arm $arm: rest ${PAUSE} s"; sleep "$PAUSE"; }
  local LABEL="single-worker, run alone"
  uv run python ../scripts/dev/timing_preflight.py --require --json "$R/preflight.json" \
    > "$R/preflight.log" 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then
    if [ "$ALLOW" = 1 ]; then
      LABEL="contended: timing preflight refused"; say "arm $arm: $LABEL (allowed)"
    else
      tail -n 6 "$R/preflight.log" >> "$MAIN"
      rm -rf "$R"; say "arm $arm: preflight refused; nothing started"; exit 3
    fi
  fi
  say "arm $arm: START ($LABEL)"
  local T1=$(date +%s)
  ARM_T0=$T1
  "${CMD[@]}" > "$LOG" 2>&1 &
  pid=$!
  uv run python $QUIET sample --out "$R/quiet.jsonl" --watch-pid $pid --period 30 \
    > "$R/sampler.log" 2>&1 &
  spid=$!
  wait $pid; rc=$?
  wait $spid
  say "arm $arm: END guard_rc=$rc after $(( $(date +%s) - T1 )) s; $(tr '\n' ' ' < "$R/time.txt" 2>/dev/null)"
  if [ $rc -ne 0 ]; then
    tail -n 5 "$R/stderr.log" >> "$MAIN" 2>/dev/null
    marker "$arm" 2 "guard rc $rc (124 passed the cap, 125 no log growth for 30 min, 126 RUNTIME_ASSERT or [ConfigFile] Error)" "$LABEL"
    exit 2
  fi
  # file checks
  if is_cpp "$arm"; then
    if [ ! -s "$R/fine1_cpp.mic.LBFS" ] \
       || ! grep -q "Total Processing Time" "$R/ParallelReconLog_1.log" 2>/dev/null \
       || ! grep -q "Total Time elapsed" "$R/ParallelReconLog_0.log" 2>/dev/null \
       || grep -l -e "RUNTIME_ASSERT" -e "\[ConfigFile\] Error" "$R"/ParallelReconLog_*.log \
            "$R/stderr.log" > /dev/null 2>&1; then
      marker "$arm" 2 "C++ file checks failed (.mic.LBFS, client or server timer, RUNTIME_ASSERT)" "$LABEL"
      exit 2
    fi
    [ "$STAGE" = slice ] && uv run python $CT/parse_logs.py "$R" --time-file "$R/time.txt" \
      --out "$R/parse.json" > "$LOGS/${STAGE}_${arm}.parse.log" 2>&1
  else
    if [ ! -f "$R/run.done" ] || [ ! -s "$R/run.json" ]; then
      marker "$arm" 2 "run.done or run.json missing" "$LABEL"; exit 2
    fi
  fi
  uv run python $QUIET verdict --samples "$R/quiet.jsonl" --limit "$LIMIT" --out "$R/quiet.json" \
    >> "$MAIN" 2>&1
  q=$?
  cp "$R/time.txt" "$LOGS/${STAGE}_${arm}.time"
  if [ $q -eq 0 ]; then
    :
  elif [ $q -eq 3 ]; then   # strict check failed, relaxed check passed (R4b): go on, labelled
    LABEL="$(grep -o '"label": "[^"]*"' "$R/quiet.json" | cut -d'"' -f4)"
    say "arm $arm: LABEL $LABEL"
  else
    LABEL="contended (the run was not quiet: load limit $LIMIT, see $R/quiet.json)"
    marker "$arm" 5 "not quiet; relabelled contended; stop, no rerun" "$LABEL"
    say "arm $arm: $LABEL"
    if [ "$CONT" = 1 ]; then NOTQUIET=1; return 0; fi
    exit 5
  fi
  marker "$arm" 0 "done" "$LABEL"
  say "arm $arm: DONE ($(cat "$RUNS/$arm.done"))"
}

for arm in $(echo "$ARMS" | tr ',' ' '); do
  run_arm "$arm"
done
say "all arms of stage $STAGE finished"
[ "$NOTQUIET" = 1 ] && { say "at least one arm was not quiet (labelled contended)"; exit 5; }
[ "$DRY" = 1 ] && say "dry run: nothing was started"
exit 0
