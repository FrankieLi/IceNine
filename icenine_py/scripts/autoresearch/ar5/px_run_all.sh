#!/bin/bash
# AR5 step 27 (task ar5aa): the live runs of gates c1 to c3 and t2, and the clean test on held7.
# A FOREGROUND script, ONE JOB AT A TIME (owner rule of 16 Nov 2026); start it from icenine_py/
# detached:
#
#   nohup caffeinate -i scripts/autoresearch/ar5/px_run_all.sh --stage live > /dev/null 2>&1 &
#
#   --stage pilot  old and new J5D on mid with --smoke 40 into <cache>_pilot (checks the wrapper)
#   --stage live   the four J5SD runs in the t2 order: old mid, new mid, new held, old held.
#                  Each: the J5D run (px_launch.py defer_run.py), J5's boundary pass on its map
#                  (high500_pass.py --name J5D --aligned 1), then step 21's screen pass and wide
#                  repair (screen_run.py --arm J5S) on a sub-cache of links, as step 22 did.
#   --stage build  build held7 (px_build_held7.py): only after logs/gate_passed exists. Nothing of
#                  held7 is read; its hashes are committed afterwards (held7_frozen.json).
#   --stage clean  S8 then the new J5SD on held7: only after logs/gate_passed AND the committed
#                  ar5/px/held7_frozen.json. Writes logs/bf_runs.done after both end (the marker
#                  that opens held7's logs and maps to the scoring scripts).
#   --dry-run      guards and commands only (nothing written)
#   --cache DIR    another cache folder
#
# Before each timed job: timing_preflight --require is RUN, its JSON KEPT in ar5/px/ and its verdict
# recorded in logs/<job>.label (single-worker if it passed, contended if it refused); a refusal is
# recorded and the job runs anyway (t1 and t2 are judged with any load). The load average is
# sampled every 30 s into logs/<job>.load while the job runs. Caps (minutes): S8 480, all else 300;
# a run at its cap writes capped.json and exits 3 (stop rule). One thread (OMP/MKL = 1). Writes
# logs/<job>.log and logs/<job>.done ("run=<rc> pass=<rc> arm=<rc>"), logs/<job>.t0 and .t1 (epoch
# seconds of the whole job), logs/px_<stage>.done at the end. Refuses if an output exists.
# Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
case "$here" in
  */.claude/worktrees/*) main="${here%%/.claude/worktrees/*}" ;;
  *) main="$here" ;;
esac
cache="${PX_CACHE:-$main/icenine_py/benchmarks/autoresearch/ar5/px/cache}"
[ -z "${PX_CACHE:-}" ] && [ "$main" = "$here" ] && cache="$here/benchmarks/autoresearch/ar5/px/cache"
dry=0
stage=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    --stage) shift; stage="$1" ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
case "$stage" in pilot|live|build|clean) ;; *) echo "--stage pilot|live|build|clean is required" >&2; exit 2 ;; esac
case "$cache" in /*) ;; *) cache="$here/$cache" ;; esac
case "$cache" in *held2*|*fine7*) echo "held2 / fine7 are not used by this step" >&2; exit 2 ;; esac
logs="$cache/logs"
pxd="$here/benchmarks/autoresearch/ar5/px"
export PX_CACHE="$cache" SCREEN_CACHE="$cache" DEFER_CACHE="$cache" BF_CACHE="$cache"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
unset VIRTUAL_ENV
L="uv run python scripts/autoresearch/ar5/px_launch.py"
pf="uv run python ../scripts/dev/timing_preflight.py"
cap_s8=480; cap_rest=300

refuse() {
  echo "refusing to start: $1" >&2
  if [ "$dry" = "0" ]; then mkdir -p "$logs"; echo "exit=2 reason=refused: $1" > "$logs/px_run_all.refused.$(date +%s)"; fi
  exit 2
}
if [ "$stage" = "build" ] || [ "$stage" = "clean" ]; then
  [ -e "$logs/gate_passed" ] || refuse "logs/gate_passed does not exist (q1, c1, c2, c3 not judged and passed)"
fi
if [ "$stage" = "clean" ]; then
  [ -e "$here/benchmarks/autoresearch/ar5/px/held7_frozen.json" ] || refuse "held7_frozen.json does not exist"
  if [ -z "$(git log -1 --format=%H -- benchmarks/autoresearch/ar5/px/held7_frozen.json)" ]; then refuse "held7_frozen.json is not committed"; fi
  if [ -n "$(git status --porcelain -- benchmarks/autoresearch/ar5/px/held7_frozen.json)" ]; then refuse "held7_frozen.json differs from its commit"; fi
fi
ferr="$(mktemp)"
if ! uv run python scripts/autoresearch/ar5/px_freeze.py --check > /dev/null 2>"$ferr"; then
  msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
  refuse "${msg:-px freeze check failed}"
fi
rm -f "$ferr"
if pgrep -f "autoresearch/ar5/px_(launch|q1|t1|build_held7)\.py" > /dev/null 2>&1; then
  refuse "a px job is already running (one at a time)"
fi

# ---- jobs
# name | arm | sample | cache dir | smoke
job_name() { echo "$1_$2"; }

sampler() { # sampler <file>: the load average every 30 s, until killed
  while true; do echo "$(date +%s) $(sysctl -n vm.loadavg)" >> "$1"; sleep 30; done
}

run_stages() { # run_stages <arm old|new|none> <sample> <cache dir> <job> [smoke]
  a="$1"; s="$2"; c="$3"; job="$4"; smoke="${5:-0}"
  sm=""; [ "$smoke" != "0" ] && sm="--smoke $smoke"
  unset prc arc
  env="RESEED_SAMPLE=$s DEFER_CACHE=$c SCREEN_CACHE=$c PX_ARM=$a"
  eval "$env PX_INDEX_LOG=$c/px_index_run_$s.json.gz $L defer_run.py --cap-min $cap_rest $sm" > "$logs/$job.log" 2>&1; rc=$?
  echo "run=$rc" > "$logs/$job.rc"
  if [ "$smoke" != "0" ]; then echo "run=$rc pass=skipped arm=skipped" > "$logs/$job.done"; return; fi
  if [ "$rc" = "0" ] && [ -e "$c/J5D_$s/run.json" ] && [ ! -e "$c/J5D_$s/capped.json" ]; then
    eval "$env PX_INDEX_LOG=$c/px_index_pass_$s.json.gz $L high500_pass.py --name J5D --aligned 1" > "$logs/$job.pass.log" 2>&1; prc=$?
    echo "pass=$prc" >> "$logs/$job.rc"
    if [ "$prc" = "0" ]; then
      mkdir -p "$c/j5sd/passes_$s/D4_J5"
      ln -sfn "$c/passes_$s/before_J5D.npy" "$c/j5sd/passes_$s/before_J5.npy"
      for f in tries.npz pass.json; do ln -sfn "$c/passes_$s/D4_J5D/$f" "$c/j5sd/passes_$s/D4_J5/$f"; done
      eval "SCREEN_PRODUCER=1 RESEED_SAMPLE=$s DEFER_CACHE=$c SCREEN_CACHE=$c/j5sd PX_ARM=$a PX_INDEX_LOG=$c/px_index_arm_$s.json.gz $L screen_run.py --arm J5S --cap-min $cap_rest" > "$logs/$job.arm.log" 2>&1; arc=$?
      echo "arm=$arc" >> "$logs/$job.rc"
    else arc=skipped; fi
  else prc=skipped; arc=skipped; fi
  echo "run=$rc pass=${prc:-skipped} arm=${arc:-skipped}" > "$logs/$job.done"
}

timed_job() { # timed_job <arm> <sample> <cache dir> <job> [smoke]
  a="$1"; s="$2"; c="$3"; job="$4"; smoke="${5:-0}"
  [ -e "$logs/$job.done" ] && { echo "refusing: $logs/$job.done exists" >&2; return 2; }
  for p in "$c/J5D_$s" "$c/passes_$s" "$c/j5sd/J5S_$s"; do
    [ -e "$p" ] && { echo "refusing: $p exists" >&2; return 2; }
  done
  mkdir -p "$c" "$c/j5sd"
  label=contended
  if $pf --require --json "$pxd/preflight_$job.json" > /dev/null 2>&1; then label=single-worker; fi
  echo "$label" > "$logs/$job.label"
  date +%s > "$logs/$job.t0"
  sampler "$logs/$job.load" & spid=$!
  run_stages "$a" "$s" "$c" "$job" "$smoke"
  kill "$spid" 2>/dev/null; wait "$spid" 2>/dev/null
  date +%s > "$logs/$job.t1"
}

s8_job() { # held7 only
  s="held7"; c="$cache/live_new"; job="S8_$s"
  [ -e "$logs/$job.done" ] && { echo "refusing: $logs/$job.done exists" >&2; return 2; }
  [ -e "$c/S8_$s" ] && { echo "refusing: $c/S8_$s exists" >&2; return 2; }
  mkdir -p "$c"
  label=contended
  if $pf --require --json "$pxd/preflight_$job.json" > /dev/null 2>&1; then label=single-worker; fi
  echo "$label" > "$logs/$job.label"
  date +%s > "$logs/$job.t0"
  sampler "$logs/$job.load" & spid=$!
  env="RESEED_SAMPLE=$s DEFER_CACHE=$c SCREEN_CACHE=$c PX_ARM=none"
  eval "$env $L fit_run.py --name S8 --cap-min $cap_s8" > "$logs/$job.log" 2>&1; rc=$?
  echo "run=$rc" > "$logs/$job.rc"
  if [ "$rc" = "0" ] && [ -e "$c/S8_$s/run.json" ] && [ ! -e "$c/S8_$s/capped.json" ]; then
    eval "$env $L fit_pass.py --name S8" > "$logs/$job.pass.log" 2>&1; prc=$?
  else prc=skipped; fi
  echo "run=$rc pass=$prc arm=skipped" > "$logs/$job.done"
  kill "$spid" 2>/dev/null; wait "$spid" 2>/dev/null
  date +%s > "$logs/$job.t1"
}

if [ "$dry" = "1" ]; then
  echo "dry run: guards passed (stage $stage, cache $cache, one job at a time)"
  case "$stage" in
    pilot) for a in old new; do echo "  pilot $a mid: PX_ARM=$a $L defer_run.py --smoke 40 (cache ${cache}_pilot/$a)"; done ;;
    live) for x in old:mid new:mid new:held old:held; do
        echo "  $x: preflight (kept) ; $L defer_run.py --cap-min $cap_rest ; high500_pass.py --name J5D --aligned 1 ; screen_run.py --arm J5S (sub-cache)"
        echo "     cache $cache/live_${x%%:*}  load sampled every 30 s"
      done ;;
    build) echo "  uv run python scripts/autoresearch/ar5/px_build_held7.py --resume" ;;
    clean) echo "  S8 held7: $L fit_run.py --name S8 --cap-min $cap_s8 ; fit_pass.py --name S8   (cache $cache/live_new)"
           echo "  new J5SD held7: PX_ARM=new defer_run.py ; high500_pass.py ; screen_run.py   (cache $cache/live_new)"
           echo "  then logs/bf_runs.done" ;;
  esac
  exit 0
fi

mkdir -p "$logs" "$pxd"
trap '' HUP
exec < /dev/null
exec >> "$logs/px_run_all.out" 2>&1
echo "start $(date +%Y-%m-%dT%H:%M:%S%z) stage=$stage"
fail=0
case "$stage" in
  pilot)
    pc="${cache}_pilot"; mkdir -p "$pc/logs"
    for a in old new; do
      logs="$pc/logs"
      timed_job "$a" mid "$pc/$a" "pilot_$a" 40 || fail=1
    done ;;
  live)
    for x in old:mid new:mid new:held old:held; do
      a="${x%%:*}"; s="${x#*:}"
      timed_job "$a" "$s" "$cache/live_$a" "live_${a}_$s" || { echo "job live $a $s ended with $?"; fail=1; }
    done ;;
  build)
    uv run python scripts/autoresearch/ar5/px_build_held7.py --resume > "$logs/held7_driver.log" 2>&1
    brc=$?; echo "$brc" > "$logs/held7_build.rc"
    [ "$brc" != "0" ] && fail=1 ;;
  clean)
    s8_job || fail=1
    timed_job new held7 "$cache/live_new" "live_new_held7" || fail=1
    for j in S8_held7 live_new_held7; do
      d="$(cat "$logs/$j.done" 2>/dev/null)"
      case "$d" in "run=0 pass=0"*) ;; *) fail=1 ;; esac
    done
    {
      echo "exit=$fail"
      date +%Y-%m-%dT%H:%M:%S%z
      for j in S8_held7 live_new_held7; do echo "$j $(cat "$logs/$j.done" 2>/dev/null || echo none)"; done
    } > "$logs/bf_runs.done.tmp" && mv "$logs/bf_runs.done.tmp" "$logs/bf_runs.done" ;;
esac
echo "exit=$fail" > "$logs/px_$stage.done"
echo "end $(date +%Y-%m-%dT%H:%M:%S%z) fail=$fail"
exit "$fail"
