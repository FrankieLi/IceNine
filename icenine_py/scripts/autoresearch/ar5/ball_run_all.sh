#!/bin/bash
# AR5 step 23, phase A steps 6 to 8: the live runs of arm J5DC and arm J5SDC, and the clean test on
# held7. A FOREGROUND script, ONE JOB AT A TIME (owner rule of 16 Nov 2026): the main session starts
# it from icenine_py/ detached:
#
#   cd icenine_py && nohup caffeinate -i scripts/autoresearch/ar5/ball_run_all.sh --stage gate \
#       --resume > /dev/null 2>&1 &
#   ... --stage rest    held6, grad, big: only after the gate was judged and logs/gate_passed exists
#   ... --stage clean   the held7 build, pilot, S8, J5D (+ J5SD) and J5DC (+ J5SDC): only after
#                       logs/gate_passed (and the live runs of stage rest) exist
#   ... --dry-run       guards and commands only (nothing written)
#   ... --cache DIR     another cache folder
#   ... --resume        after an interruption or a FAILED earlier attempt: a job whose .done says it
#                       ended 0 with its outputs is kept; the leftovers of every other job and an
#                       old marker are MOVED to <cache>/interrupted_<time>/ and those jobs run from
#                       scratch. Safe on a first start too.
#
# A job is: <arm>'s run (J5DC: ball_launch.py ball_run.py; J5D: defer_run.py), J5's boundary pass
# on its map (high500_pass.py --name <arm> --aligned 1), then step 21's screen pass and wide repair
# on that map (screen_run.py --arm J5S, run from a sub-cache <cache>/<j5sd|j5sdc> whose
# before_J5.npy and D4_J5/ are links to the arm's pass files, so that step 21's script is imported
# unchanged): J5SDC = <cache>/j5sdc/J5S_<s>, J5SD = <cache>/j5sd/J5S_<s>. Markers
# logs/ball_<stage>.done ("exit=<code>", "reason=<text>", end time, one line a job). Stage clean
# writes logs/ball_runs.done (which opens held7's logs and maps to this step's scripts: nobody reads
# them before it exists).
# Caps (minutes): S8 480, everything else 300. A run at its cap writes capped.json and exits 3: its
# pass is skipped and the marker says so (stop rule).
# Guards (exit 2, nothing started; a note logs/ball_run_all.refused.<time> holds the reason): an
# output of a job or the stage's marker exists (not with --resume); ball_freeze.py --check fails;
# logs/gate_passed missing (rest, clean); a heavy run of anyone (the C++ binary, pilot_bfs,
# guard.py, mpirun) or of this study's own scripts is running (the script WAITS for the former, it
# never stops a process); an input is missing. Nothing of held2 or fine7 is used.
# Plain bash 3.2 (macOS default).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${BALL_CACHE:-$here/benchmarks/autoresearch/ar5/ball/cache}"
dry=0
resume=0
stage=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    --stage) shift; stage="$1" ;;
    --resume) resume=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
case "$stage" in gate|rest|clean) ;; *) echo "--stage gate|rest|clean is required (got '$stage')" >&2; exit 2 ;; esac
case "$cache" in /*) ;; *) cache="$here/$cache" ;; esac
logs="$cache/logs"
marker="$logs/ball_$stage.done"
pilot_cache="${cache}_pilot"
export SCREEN_CACHE="$cache" DEFER_CACHE="$cache" BALL_CACHE="$cache"
unset VIRTUAL_ENV
cap_s8=480; cap_rest=300
case "$stage" in
  gate) arms="J5DC:mid J5DC:held" ;;
  rest) arms="J5DC:held6 J5DC:grad J5DC:big" ;;
  clean) arms="S8:held7 J5D:held7 J5DC:held7" ;;
esac

end() { # end <exit> <reason>
  code="$1"; reason="$2"
  if [ "$dry" = "0" ]; then
    mkdir -p "$logs"
    {
      echo "exit=$code"
      echo "reason=$reason"
      date +%Y-%m-%dT%H:%M:%S%z
      for a in $arms; do echo "$a $(cat "$logs/${a%%:*}_${a#*:}.done" 2>/dev/null || echo none)"; done
    } > "$marker.tmp" && mv "$marker.tmp" "$marker"
    if [ "$stage" = "clean" ]; then cp "$marker" "$logs/ball_runs.done"; fi
  else
    echo "dry run: would end with exit=$code reason=$reason"
  fi
  exit "$code"
}
refuse() {
  echo "refusing to start: $1" >&2
  if [ "$dry" = "0" ]; then mkdir -p "$logs"; echo "exit=2 reason=refused: $1" > "$logs/ball_run_all.refused.$(date +%s)"; fi
  exit 2
}
if [ "$stage" != "gate" ] && [ ! -e "$logs/gate_passed" ]; then
  refuse "logs/gate_passed does not exist (the gate has not been judged and passed): $logs/gate_passed"
fi
if [ "$stage" = "clean" ] && [ ! -e "$logs/ball_rest.done" ]; then
  refuse "logs/ball_rest.done does not exist (the live runs of stage rest have not ended)"
fi
if [ "$stage" = "gate" ]; then
  [ -e "$here/benchmarks/autoresearch/ar5/ball/replay.json" ] || refuse "ar5/ball/replay.json does not exist (g1 not run)"
  [ -e "$here/benchmarks/autoresearch/ar5/ball/g0.json" ] || refuse "ar5/ball/g0.json does not exist (g0 not run)"
fi

files() { # files <arm> <sample>: the outputs of a job, relative to the cache
  a="$1"; s="$2"
  case "$a" in
    S8) echo "S8_$s passes_$s/D4_S8 passes_$s/before_S8.npy passes_$s/after_S8.npy logs/S8_$s.log logs/S8_$s.run_rc logs/pass_S8_$s.log logs/pass_S8_$s.rc logs/S8_$s.done" ;;
    *) sub="$(subname "$a")"
       echo "${a}_$s passes_$s/D4_$a passes_$s/before_$a.npy passes_$s/after_$a.npy $sub/passes_$s $sub/J5S_$s logs/${a}_$s.log logs/${a}_$s.run_rc logs/pass_${a}_$s.log logs/pass_${a}_$s.rc logs/arm_${a}_$s.log logs/arm_${a}_$s.rc logs/${a}_$s.done" ;;
  esac
}
subname() { if [ "$1" = "J5DC" ]; then echo j5sdc; else echo j5sd; fi; }
finished() { # finished <arm> <sample>
  a="$1"; s="$2"
  d="$(cat "$logs/${a}_$s.done" 2>/dev/null)"
  [ -e "$cache/${a}_$s/run.json" ] && [ ! -e "$cache/${a}_$s/capped.json" ] \
    && [ -e "$cache/passes_$s/after_$a.npy" ] && [ -e "$cache/passes_$s/D4_$a/pass.json" ] || return 1
  if [ "$a" = "S8" ]; then [ "$d" = "run=0 pass=0 arm=skipped" ]; return; fi
  sub="$(subname "$a")"
  [ "$d" = "run=0 pass=0 arm=0" ] && [ -e "$cache/$sub/J5S_$s/repair.json" ] && [ -e "$cache/$sub/J5S_$s/after.npy" ]
}

# ---- guards
kept=""; redo=""
if [ -e "$marker" ] && [ "$resume" = "0" ]; then refuse "the marker $marker already exists (use --resume)"; fi
for a in $arms; do
  arm="${a%%:*}"; s="${a#*:}"
  if [ "$resume" = "1" ]; then
    if finished "$arm" "$s"; then kept="$kept ${arm}_$s"; else redo="$redo ${arm}_$s"; fi
  else
    for p in $(files "$arm" "$s"); do
      if [ -e "$cache/$p" ]; then refuse "$cache/$p already exists"; fi
    done
  fi
done
ferr="$(mktemp)"
if ! uv run python scripts/autoresearch/ar5/ball_freeze.py --check > /dev/null 2>"$ferr"; then
  msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
  refuse "${msg:-ball freeze check failed}"
fi
rm -f "$ferr"
if pgrep -f "autoresearch/ar5/(ball|defer|screen|high500|fit)_(run|pass|launch|replay|build_held7)\.py" > /dev/null 2>&1; then
  refuse "a ball / defer / screen / high500 / fit run, pass, launch, replay or build process is already running"
fi
st="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
need=""
for a in $arms; do
  s="${a#*:}"
  [ "$s" = "held7" ] && continue   # held7's files do not exist before its build
  if [ "$s" = "grad" ]; then g="$here/benchmarks/autoresearch/ar3/gradient/sample/gradient_grid.mic"
  else g="$here/benchmarks/autoresearch/ar3/$s/gradient/sample/gradient_grid.mic"; fi
  need="$need $g $st/${s}_clean.npy"
done
for p in $need; do [ -e "$p" ] || refuse "input missing: $p"; done

# ---- the commands
ev() { echo "RESEED_SAMPLE=$1 DEFER_CACHE=$cache BALL_CACHE=$cache SCREEN_CACHE=${2:-$cache}"; }
L="uv run python scripts/autoresearch/ar5/ball_launch.py"
run_cmd() { # run_cmd <arm> <sample> [cap] [cache] [extra]
  a="$1"; s="$2"; cap="${3:-$cap_rest}"; c="${4:-$cache}"
  case "$a" in
    S8) echo "RESEED_SAMPLE=$s DEFER_CACHE=$cache BALL_CACHE=$cache SCREEN_CACHE=$c $L fit_run.py --name S8 --cap-min $cap ${5:-}" ;;
    J5D) echo "RESEED_SAMPLE=$s DEFER_CACHE=$cache BALL_CACHE=$cache SCREEN_CACHE=$c $L defer_run.py --cap-min $cap ${5:-}" ;;
    J5DC) echo "RESEED_SAMPLE=$s DEFER_CACHE=$cache BALL_CACHE=$cache SCREEN_CACHE=$c $L ball_run.py --cap-min $cap ${5:-}" ;;
  esac
}
pass_cmd() {
  if [ "$1" = "S8" ]; then echo "$(ev "$2") $L fit_pass.py --name S8"
  else echo "$(ev "$2") $L high500_pass.py --name $1 --aligned 1"; fi
}
arm_cmd() { echo "SCREEN_PRODUCER=1 $(ev "$2" "$cache/$(subname "$1")") $L screen_run.py --arm J5S --cap-min $cap_rest"; }

if [ "$dry" = "1" ]; then
  echo "dry run: guards passed (stage $stage, cache $cache, one job at a time)"
  [ "$resume" = "1" ] && echo "  resume: KEPT:${kept:- none} RERUN (leftovers moved to $cache/interrupted_<time>/):${redo:- none}"
  if [ "$stage" = "clean" ]; then
    echo "  build: uv run python scripts/autoresearch/ar5/ball_build_held7.py --resume"
    echo "  pilot J5DC: $(run_cmd J5DC held7 5 "$pilot_cache" "--smoke 30")"
  fi
  for a in $arms; do
    arm="${a%%:*}"; s="${a#*:}"; cap=$cap_rest; [ "$arm" = "S8" ] && cap=$cap_s8
    echo "  $arm $s: $(run_cmd "$arm" "$s" "$cap")"
    echo "    then: $(pass_cmd "$arm" "$s")"
    [ "$arm" != "S8" ] && echo "    then (sub-cache links to the pass files): $(arm_cmd "$arm" "$s")"
  done
  echo "  then: marker $marker"
  exit 0
fi

mkdir -p "$logs"
trap '' HUP
exec < /dev/null
exec >> "$logs/ball_run_all.out" 2>&1
if [ "$resume" = "1" ]; then
  idir="$cache/interrupted_$(date +%Y%m%d-%H%M%S)"
  if [ -e "$marker" ]; then mkdir -p "$idir/logs"; mv "$marker" "$idir/logs/$(basename "$marker")"; fi
  for n in $redo; do
    for p in $(files "${n%%_*}" "${n#*_}"); do
      if [ -e "$cache/$p" ] || [ -L "$cache/$p" ]; then mkdir -p "$idir/$(dirname "$p")"; mv "$cache/$p" "$idir/$p"; fi
    done
  done
  echo "resume $(date +%Y-%m-%dT%H:%M:%S%z) kept:${kept:- none} rerun:${redo:- none}" >> "$logs/ball_run_all.resume"
else
  date +%Y-%m-%dT%H:%M:%S%z > "$logs/ball_run_all.start"
fi
echo "start $(date +%Y-%m-%dT%H:%M:%S%z) stage=$stage resume=$resume"

busy() { ps -axo comm,args | grep -E '^(/[^ ]*/)?IceNine |pilot_bfs|guard\.py|mpirun' | grep -v grep > /dev/null; }
wait_quiet() { while busy; do echo "waiting: another heavy run is active ($(date +%H:%M))"; sleep 120; done; }

mk_sub() { # mk_sub <arm> <sample>
  a="$1"; s="$2"; sub="$cache/$(subname "$a")"
  mkdir -p "$sub/passes_$s/D4_J5"
  ln -sfn "$cache/passes_$s/before_$a.npy" "$sub/passes_$s/before_J5.npy"
  for f in tries.npz pass.json; do ln -sfn "$cache/passes_$s/D4_$a/$f" "$sub/passes_$s/D4_J5/$f"; done
}

job() { # job <arm> <sample>
  a="$1"; s="$2"; cap=$cap_rest; [ "$a" = "S8" ] && cap=$cap_s8
  case " $kept " in *" ${a}_$s "*) return ;; esac
  wait_quiet
  eval "$(run_cmd "$a" "$s" "$cap")" > "$logs/${a}_$s.log" 2>&1
  rc=$?
  echo "$rc" > "$logs/${a}_$s.run_rc"
  d="$cache/${a}_$s"
  prc=skipped; arc=skipped
  if [ "$rc" = "0" ] && [ -e "$d/run.json" ] && [ ! -e "$d/capped.json" ]; then
    wait_quiet
    eval "$(pass_cmd "$a" "$s")" > "$logs/pass_${a}_$s.log" 2>&1
    prc=$?
    echo "$prc" > "$logs/pass_${a}_$s.rc"
    if [ "$prc" = "0" ] && [ "$a" != "S8" ]; then
      mk_sub "$a" "$s"
      wait_quiet
      eval "$(arm_cmd "$a" "$s")" > "$logs/arm_${a}_$s.log" 2>&1
      arc=$?
      echo "$arc" > "$logs/arm_${a}_$s.rc"
    fi
  else
    echo skipped > "$logs/pass_${a}_$s.rc"
  fi
  echo "run=$rc pass=$prc arm=$arc" > "$logs/${a}_$s.done"
}

fail=0; capped=""
if [ "$stage" = "clean" ]; then
  rm -f "$logs/held7_pilot.ok" "$logs/held7_pilot.failed"
  wait_quiet
  uv run python scripts/autoresearch/ar5/ball_build_held7.py --resume > "$logs/held7_driver.log" 2>&1
  brc=$?
  echo "$brc" > "$logs/held7_build.rc"
  if [ "$brc" != "0" ]; then echo "build failed ($brc)" > "$logs/held7_pilot.failed"; end 1 "the held7 build failed (see held7_driver.log; held7 not read)"; fi
  mkdir -p "$pilot_cache/logs"
  wait_quiet
  eval "SCREEN_CACHE=$pilot_cache $(run_cmd J5DC held7 5 "$pilot_cache" "--smoke 30")" > "$logs/pilot_J5DC.log" 2>&1
  prc=$?
  echo "$prc" > "$logs/pilot_J5DC.rc"
  if { [ "$prc" != "0" ] && [ "$prc" != "3" ]; } || [ ! -e "$pilot_cache/J5DC_held7_smoke" ]; then
    echo "pilot failed (exit $prc or no folder; not read for results)" > "$logs/held7_pilot.failed"
    end 1 "the held7 pilot failed (see pilot_J5DC.rc; held7 not read)"
  fi
  echo ok > "$logs/held7_pilot.ok"
fi
for a in $arms; do job "${a%%:*}" "${a#*:}"; done

for a in $arms; do
  arm="${a%%:*}"; s="${a#*:}"
  case " $kept " in *" ${arm}_$s "*) continue ;; esac
  rc="$(cat "$logs/${arm}_$s.run_rc" 2>/dev/null || echo 99)"
  prc="$(cat "$logs/pass_${arm}_$s.rc" 2>/dev/null || echo 99)"
  arc="$(cat "$logs/arm_${arm}_$s.rc" 2>/dev/null || echo skipped)"
  if [ "$rc" = "3" ]; then capped="${arm}_$s"; fail=1; elif [ "$rc" != "0" ]; then fail=1; fi
  [ "$prc" != "0" ] && fail=1
  if [ "$arm" != "S8" ]; then [ "$arc" != "0" ] && fail=1; [ "$arc" = "3" ] && capped="arm_${arm}_$s"; fi
done
date +%Y-%m-%dT%H:%M:%S%z > "$logs/ball_run_all.end"
if [ -n "$capped" ]; then end 1 "STOP RULE: a run reached its cap: $capped"; fi
if [ "$fail" = "1" ]; then end 1 "a run, a pass or an arm failed (see the rc files; not read for results)"; fi
end 0 "all runs, passes and arms finished"
