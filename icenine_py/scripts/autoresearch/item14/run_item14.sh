#!/bin/bash
# Item 14 and 3: both speed keys together on E0 and the pilot, with a timing A/B (PLAN.md).
#
#   nohup caffeinate -i bash scripts/autoresearch/item14/run_item14.sh > /dev/null 2>&1 &
#   bash scripts/autoresearch/item14/run_item14.sh --dry-run     # plumbing only, no heavy work
#
# Order: E0 check (stored arms off and S, 20 voxels); E0 F and FS (4 workers); untimed warm-up;
# A1, B1, A2, B2 (120 s pause + preflight before each, load sampler during, quiet verdict after,
# identity check after); split.py; closing preflight.
#
# Runs in the foreground (no self-detach); the caller backgrounds it. Everything is logged under
# benchmarks/autoresearch/item14/cache/logs/; the marker cache/item14.done is written when the script
# ends, with its exit code and reason; cache/exit_codes.txt has one line per step.
#   exit 0 all done; 2 a step failed (E0 check, a repeat not identical, a run fell back from the
#   fused path, ...); 3 a preflight refused; 4 passed 6 hours; 5 a run was neither strictly quiet
#   nor relaxed-ok (R4b); 6 (at the very end, after all arms ran) a labelled run's repeat differs
#   by more than 2% in real time: no time is reported for that arm (cache/pair_<arm>.json).
# R4b (D-20261114-N1): quiet.py verdict exit 3 (strict check failed, median load within the limit,
# every busy process a macOS daemon) does not stop the script; the label is logged and kept in
# cache/quiet/<run>.json and cache/pair_<arm>.json. Exits 2-5 stop at once and no time is reported.
set -u
cd "$(dirname "$0")/../../.." || exit 99   # icenine_py/

DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1
C=benchmarks/autoresearch/item14/cache
[ "$DRY" = 1 ] && C=$C/dry
IT=scripts/autoresearch/item14
STORED_A=benchmarks/autoresearch/ar2/key/runs/key_off_fused
LIMIT=$(cat $IT/r4_limit.txt)      # R4 load limit, fixed in phase A before any timed run
PAUSE=120
MAXS=21600
[ "$DRY" = 1 ] && PAUSE=0
mkdir -p "$C/logs" "$C/runs" "$C/quiet"
MAIN=$C/logs/run_item14.log
CODES=$C/exit_codes.txt
MARKER=$C/item14.done
[ "$DRY" = 1 ] && MARKER=$C/item14_dry.done
: > "$CODES"
rm -f "$MARKER"
T0=$(date +%s)
export PYTHONWARNINGS="error:FusedCost is not supported:RuntimeWarning"

say() { echo "$(date '+%F %T') $*" >> "$MAIN"; }

finish() { # rc reason
  say "STOP/END: $2 (exit $1)"
  echo "exit=$1 reason=$2 elapsed_s=$(( $(date +%s) - T0 )) codes=$CODES log=$MAIN" > "$MARKER"
  exit "$1"
}

check_time() { [ $(( $(date +%s) - T0 )) -gt $MAXS ] && finish 4 "passed 6 hours"; return 0; }

step() { # name cmd...   (untimed step; a non-zero exit stops the script)
  local n=$1; shift
  say "START $n: $*"
  if [ "$DRY" = 1 ]; then say "DRY: $n not run"; echo "$n rc=dry" >> "$CODES"; return 0; fi
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

gate() { # name  (preflight; refusal stops the script)
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

timed() { # name cmd...   (gated, sampled, time -p; not quiet or failing stops the script)
  local n=$1; shift
  [ "$PAUSE" -gt 0 ] && { say "pause ${PAUSE}s before $n"; sleep "$PAUSE"; }
  gate "$n"
  say "START timed $n (single-worker): $*"
  local period=30
  if [ "$DRY" = 1 ]; then set -- sleep 3; period=1; fi
  /usr/bin/time -p -o "$C/logs/$n.time" "$@" > "$C/logs/$n.log" 2>&1 &
  local pid=$!
  uv run python $IT/quiet.py sample --out "$C/quiet/$n.jsonl" --watch-pid $pid --period $period \
    > "$C/logs/$n.sampler.log" 2>&1 &
  local spid=$!
  wait $pid
  local rc=$?
  wait $spid
  echo "$n rc=$rc" >> "$CODES"
  say "END timed $n rc=$rc, time: $(tr '\n' ' ' < "$C/logs/$n.time")"
  if [ $rc -ne 0 ]; then
    tail -n 8 "$C/logs/$n.log" >> "$MAIN"
    finish 2 "timed run $n failed (rc $rc), e.g. a fall-back from the fused path; see $C/logs/$n.log"
  fi
  uv run python $IT/quiet.py verdict --samples "$C/quiet/$n.jsonl" --limit "$LIMIT" \
    --out "$C/quiet/$n.json" >> "$MAIN" 2>&1
  local q=$?
  echo "quiet_$n rc=$q" >> "$CODES"
  if [ $q -eq 3 ] && [ "$n" = split ]; then  # a single run has no repeat for R4b (iii)
    if [ "$DRY" != 1 ]; then finish 5 "run split was not strictly quiet and has no repeat (R4b needs one); no time is reported"; fi
  fi
  if [ $q -eq 3 ]; then  # R4b: strict check failed, relaxed check passed; go on, with the label
    say "LABEL $n: $(grep -o '"label": "[^"]*"' "$C/quiet/$n.json" | cut -d'"' -f4)"
    q=0
  fi
  if [ $q -ne 0 ] && [ "$DRY" = 1 ]; then say "DRY: run $n was not quiet; a real run would stop here (exit 5)"; q=0; fi
  if [ $q -ne 0 ]; then finish 5 "run $n was not quiet (limit $LIMIT, see $C/quiet/$n.json); no time is reported"; fi
  check_time
}

same() { # name dirA dirB   (identity check, also in the dry run; a difference stops the script)
  say "START same_$1: $2 vs $3"
  uv run python $IT/same_run.py "$2" "$3" > "$C/logs/same_$1.log" 2>&1
  local rc=$?
  echo "same_$1 rc=$rc" >> "$CODES"
  say "END same_$1 rc=$rc: $(tail -n 1 "$C/logs/same_$1.log")"
  [ $rc -ne 0 ] && finish 2 "identity check same_$1 failed (stop rule): $(tail -n 1 "$C/logs/same_$1.log")"
  return 0
}

PAIR_FAIL=""
pair() { # arm run1 run2   (R4b: a labelled run needs its repeat within 2% in real time)
  uv run python $IT/quiet.py pair --arm "$1" --runs "$2" "$3" --cache "$C" > "$C/logs/pair_$1.log" 2>&1
  local rc=$?
  echo "pair_$1 rc=$rc" >> "$CODES"
  say "PAIR $1: $(tail -n 1 "$C/logs/pair_$1.log")"
  [ $rc -ne 0 ] && PAIR_FAIL="$PAIR_FAIL $1"
  return 0
}

PC=()
pilot() { # name config  -> PC (the command words of one pilot BFS run)
  PC=(uv run python scripts/phase_d/pilot_bfs.py --arm mc --variant clean --out "$C/runs"
    --name "$1" --config "$2" --timing-label "single-worker; quiet verdict in quiet/$1.json")
}

say "item14 start (dry=$DRY), R4 load limit $LIMIT, host load $(uptime | sed 's/.*averages: //')"
[ "$DRY" = 1 ] && say "DRY RUN: heavy steps are echoed, not run"

# --- E0: the stored arms reproduce; arms F and FS (4 workers, contended, no timing claim)
step e0_check uv run python $IT/e0_arms.py check --n 20 --workers 4
step e0_F uv run python $IT/e0_arms.py run --spec off --fused --workers 4
step e0_FS uv run python $IT/e0_arms.py run --spec "E3+E4,tau=0.25" --fused --workers 4

# --- warm-up (untimed): arm A config, 5 voxels, fills the page cache
CFG_A=scripts/autoresearch/ar2/key_off_fused.config
CFG_B=$IT/both_keys.config
pilot warmup $CFG_A
step warmup "${PC[@]}" --max-voxels 5

# --- timing A/B: A1, B1, A2, B2
pilot A1 $CFG_A; timed A1 "${PC[@]}"
if [ "$DRY" = 1 ]; then same A1_stored "$STORED_A" "$STORED_A"; else same A1_stored "$C/runs/A1" "$STORED_A"; fi
pilot B1 $CFG_B; timed B1 "${PC[@]}"
pilot A2 $CFG_A; timed A2 "${PC[@]}"
if [ "$DRY" = 1 ]; then same A2_A1 "$STORED_A" "$STORED_A"; else same A2_A1 "$C/runs/A2" "$C/runs/A1"; fi
pair A A1 A2
pilot B2 $CFG_B; timed B2 "${PC[@]}"
if [ "$DRY" = 1 ]; then same B2_B1 "$STORED_A" "$STORED_A"; else same B2_B1 "$C/runs/B2" "$C/runs/B1"; fi
pair B B1 B2

# --- split of a full search by section (alone, gated)
timed split uv run python $IT/split.py --voxels 15901,17034 --out "$C/split" --label "single-worker"

# --- closing preflight (report only)
uv run python ../scripts/dev/timing_preflight.py --json "$C/preflight_final.json" \
  > "$C/logs/preflight_final.log" 2>&1
say "closing preflight written"
if [ -n "$PAIR_FAIL" ]; then
  finish 6 "no time reported for arm(s):$PAIR_FAIL (labelled run, repeat differs by more than 2%; see $C/pair_*.json)"
fi
finish 0 "all steps done"
