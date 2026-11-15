#!/bin/bash
# AR5 step 20, phase A step 5: the held4 build, arms S8 and J5 on held4 with their boundary passes,
# then the CMA repair on held4's J5 map. A FOREGROUND script: the main session starts it from
# icenine_py/ detached:
#
#   cd icenine_py && nohup caffeinate -i scripts/autoresearch/ar5/cmarepair_run_all.sh --resume \
#       > /dev/null 2>&1 &
#   ... --dry-run      guards and commands only (nothing written)
#   ... --cache DIR    another cache folder (tests)        ... --jobs N   at most N heavy (<= 4)
#   ... --smoke        test form: short jobs (--smoke 30, a cap of --smoke-cap M minutes, default
#                      2; no pass, no repair): S8 and J5 on held4, in <cache_smoke>; the held4 build
#                      must exist. Exit 3 (the cap) is these jobs' normal end.
#   ... --resume       after an interruption or a FAILED earlier attempt: the held4 build skips its
#                      finished stages and moves an unfinished stage's output aside; a job whose
#                      .done says its run and pass ended 0 (the repair: exit=0) with its outputs is
#                      kept; the leftovers of every other job and an old marker (failed or not) are
#                      MOVED to <cache>/interrupted_<time>/ and those jobs run from scratch. It is
#                      safe to use on a first start too (nothing to keep).
#
# Order (at most 4 heavy processes at once; the real run uses at most 2):
#   held4 build (cmarepair_build_held4.py: sample, render with 1 worker, stack; about 15 minutes);
#   then, side by side: S8 held4 (the long pole: fit_run.py --name S8, then fit_pass.py --name S8)
#   and J5 held4 (high500_run.py, then high500_pass.py --name J5 --aligned 1, then the repair:
#   cmarepair_run.py, arm J5 -> J5C). Stats are left for later.
# Caps (--cap-min): S8 480, J5 300. A run at its cap writes capped.json and exits 3; its pass is
# skipped; the marker says so (stop rule). The repair is a few minutes; no cap.
#
# Guard 0 (exit 2, nothing written): the real path of an output lies outside the cache
# (stallhit_cache_guard.py). Guards (exit 2, nothing started; a note
# logs/cmarepair_run_all.refused.<time> holds the reason): an output of a job or the final marker
# exists (not with --resume); ar5/cmarepair/gate.json is missing (the gate passed on the disclosed
# samples); cmarepair_freeze.py --check fails (frozen.json committed, clean, hashes agree); a
# cmarepair / high500 / fit run, pass, launch or build process is running; an input is missing.
# Nothing of held2 is used (cmarepair_common refuses RESEED_SAMPLE=held2).
#
# Files, in <cache>/logs (cache = benchmarks/autoresearch/ar5/cmarepair/cache, untracked):
#   <arm>_held4.log, .run_rc, pass_<arm>_held4.log, pass_<arm>_held4.rc, <arm>_held4.done
#   ("run=<code> pass=<code|skipped>"), repair_held4.log / .rc / .done ("exit=<code>"),
#   held4_driver.log / .rc and held4_<stage>.log / .rc (the build), cmarepair_run_all.start / .end,
#   and the MARKER cmarepair_runs.done: first line "exit=<code>", then "reason=<text>", the end
#   time and one line per job. The script's exit code is the marker's. Run folders <cache>/<arm>_
#   held4/, passes in <cache>/passes_held4/, the repair in <cache>/J5C_held4/. The readers of this
#   step's layer (cmarepair_guard) open nothing of held4 before the marker exists; the repair job
#   alone (CMAREPAIR_PRODUCER=1, set here for it) opens its own inputs by name. NOBODY reads held4's
#   logs or maps before the marker says all runs and the repair ended (plan stop rule). Wall times
#   are contended and are not timing results.
# Plain bash 3.2 (macOS default).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${CMAREPAIR_CACHE:-$here/benchmarks/autoresearch/ar5/cmarepair/cache}"
dry=0
resume=0
jobs=4
smoke=0
smoke_cap=2
cache_given=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1"; cache_given=1 ;;
    --jobs) shift; jobs="$1" ;;
    --smoke) smoke=1 ;;
    --smoke-cap) shift; smoke_cap="$1" ;;
    --resume) resume=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
if [ "$jobs" -gt 4 ] || [ "$jobs" -lt 1 ]; then echo "--jobs must be 1 to 4" >&2; exit 2; fi
if [ "$smoke" = "1" ] && [ "$cache_given" = "0" ]; then
  cache="$here/benchmarks/autoresearch/ar5/cmarepair/cache_smoke"
fi
logs="$cache/logs"
marker="$logs/cmarepair_runs.done"
export CMAREPAIR_CACHE="$cache"
unset VIRTUAL_ENV

gate="$here/benchmarks/autoresearch/ar5/cmarepair/gate.json"
st="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
grid="$here/benchmarks/autoresearch/ar3/held4/gradient/sample/gradient_grid.mic"
if [ "$smoke" = "1" ]; then
  joblist="S8:$smoke_cap J5:$smoke_cap"
  smoke_arg=" --smoke 30"
else
  joblist="S8:480 J5:300"
  smoke_arg=""
fi

end() { # end <exit> <reason>: write the marker and exit
  code="$1"; reason="$2"
  if [ "$dry" = "0" ]; then
    mkdir -p "$logs"
    {
      echo "exit=$code"
      echo "reason=$reason"
      date +%Y-%m-%dT%H:%M:%S%z
      echo "held4_build driver=$(cat "$logs/held4_driver.rc" 2>/dev/null || echo none)"
      for j in $joblist; do
        a="${j%%:*}"
        echo "${a}_held4 run=$(cat "$logs/${a}_held4.run_rc" 2>/dev/null || echo none) pass=$(cat "$logs/pass_${a}_held4.rc" 2>/dev/null || echo none)"
      done
      if [ "$smoke" = "0" ]; then
        echo "repair_held4 exit=$(cat "$logs/repair_held4.rc" 2>/dev/null || echo none)"
      fi
    } > "$marker.tmp" && mv "$marker.tmp" "$marker"
  else
    echo "dry run: would end with exit=$code reason=$reason"
  fi
  exit "$code"
}

refuse() { # a refusal never writes the final marker (that would block the real start)
  echo "refusing to start: $1" >&2
  if [ "$dry" = "0" ]; then
    mkdir -p "$logs"
    echo "exit=2 reason=refused: $1" > "$logs/cmarepair_run_all.refused.$(date +%s)"
  fi
  exit 2
}

rundir() { if [ "$smoke" = "1" ]; then echo "$1_held4_smoke"; else echo "$1_held4"; fi; }
jfiles() { # jfiles <arm>: the outputs of a job (paths relative to the cache)
  if [ "$smoke" = "1" ]; then
    echo "$(rundir "$1") logs/$1_held4.log logs/$1_held4.done logs/$1_held4.run_rc logs/pass_$1_held4.log logs/pass_$1_held4.rc"
  else
    echo "$1_held4 passes_held4/D4_$1 passes_held4/before_$1.npy passes_held4/after_$1.npy logs/$1_held4.log logs/$1_held4.done logs/$1_held4.run_rc logs/pass_$1_held4.log logs/pass_$1_held4.rc"
  fi
}
repair_files() { echo "J5C_held4 logs/repair_held4.log logs/repair_held4.rc logs/repair_held4.done"; }
job_finished() { # job_finished <arm>
  d="$(cat "$logs/$1_held4.done" 2>/dev/null)"
  if [ "$smoke" = "1" ]; then
    { [ "$d" = "run=0 pass=skipped" ] && [ -e "$cache/$(rundir "$1")/run.json" ]; } \
      || { [ "$d" = "run=3 pass=skipped" ] && [ -e "$cache/$(rundir "$1")/capped.json" ]; }
    return
  fi
  [ "$d" = "run=0 pass=0" ] && [ -e "$cache/$1_held4/run.json" ] && [ ! -e "$cache/$1_held4/capped.json" ] \
    && [ -e "$cache/passes_held4/after_$1.npy" ] && [ -e "$cache/passes_held4/D4_$1/pass.json" ] \
    && [ -e "$cache/passes_held4/before_$1.npy" ]
}
repair_finished() {
  [ "$(cat "$logs/repair_held4.done" 2>/dev/null)" = "exit=0" ] && [ -e "$cache/J5C_held4/repair.json" ] \
    && [ -e "$cache/J5C_held4/after.npy" ]
}
build_finished() {
  [ "$(cat "$logs/held4_driver.rc" 2>/dev/null)" = "0" ] && [ -e "$grid" ] && [ -e "$st/held4_clean.npy" ]
}

# ---- guard 0: every output must really lie in this cache, not behind a link into another cache
outs="logs"
for j in $joblist; do
  a="${j%%:*}"; outs="$outs $(jfiles "$a")"
done
[ "$smoke" = "0" ] && outs="$outs passes_held4 $(repair_files)"
if ! uv run python scripts/autoresearch/ar5/stallhit_cache_guard.py --cache "$cache" $outs; then
  echo "refusing to start: an output's real path lies outside the cache $cache" >&2
  exit 2
fi

# ---- guard 1: no output exists (--resume: finished jobs are kept, leftovers and an old marker are
# moved aside later)
kept=""; redo=""
if [ -e "$marker" ] && [ "$resume" = "0" ]; then refuse "the marker $marker already exists (use --resume)"; fi
for j in $joblist; do
  a="${j%%:*}"
  if [ "$resume" = "1" ]; then
    if job_finished "$a"; then kept="$kept ${a}_held4"; else redo="$redo ${a}_held4"; fi
    continue
  fi
  for p in $(jfiles "$a"); do
    if [ -e "$cache/$p" ]; then refuse "$cache/$p already exists"; fi
  done
done
if [ "$smoke" = "0" ]; then
  if [ "$resume" = "1" ]; then
    if repair_finished; then kept="$kept repair_held4"; else redo="$redo repair_held4"; fi
  else
    for p in $(repair_files); do
      if [ -e "$cache/$p" ]; then refuse "$cache/$p already exists"; fi
    done
  fi
fi

# ---- guard 2: the gate and the freeze
[ -e "$gate" ] || refuse "ar5/cmarepair/gate.json does not exist (the gate has not passed)"
ferr="$(mktemp)"
if ! uv run python scripts/autoresearch/ar5/cmarepair_freeze.py --check > /dev/null 2>"$ferr"; then
  msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
  refuse "${msg:-cmarepair freeze check failed}"
fi
rm -f "$ferr"

# ---- guard 3: nothing of this step (or an arm it builds on) is running
if pgrep -f "autoresearch/ar5/(cmarepair|high500|fit)_(run|pass|launch|build_held4|make_held3)\.py" > /dev/null 2>&1; then
  refuse "a cmarepair / high500 / fit run, pass, launch or build process is already running"
fi

# ---- guard 4: inputs
need="$here/benchmarks/autoresearch/ar3/reseed/remedies_pick.json"
if [ "$smoke" = "1" ]; then need="$need $grid $st/held4_clean.npy"; fi
for p in $need; do
  if [ ! -e "$p" ]; then refuse "input missing: $p"; fi
done

run_cmd() { # run_cmd <arm> <cap>
  case "$1" in
    S8) echo "RESEED_SAMPLE=held4 uv run python scripts/autoresearch/ar5/cmarepair_launch.py fit_run.py --name S8 --cap-min $2$smoke_arg" ;;
    J5) echo "RESEED_SAMPLE=held4 uv run python scripts/autoresearch/ar5/cmarepair_launch.py high500_run.py --cap-min $2$smoke_arg" ;;
  esac
}
pass_cmd() { # pass_cmd <arm>
  case "$1" in
    S8) echo "RESEED_SAMPLE=held4 uv run python scripts/autoresearch/ar5/cmarepair_launch.py fit_pass.py --name S8" ;;
    J5) echo "RESEED_SAMPLE=held4 uv run python scripts/autoresearch/ar5/cmarepair_launch.py high500_pass.py --name J5 --aligned 1" ;;
  esac
}
repair_cmd() {
  echo "CMAREPAIR_PRODUCER=1 RESEED_SAMPLE=held4 uv run python scripts/autoresearch/ar5/cmarepair_run.py --arm J5"
}
build_cmd() { echo "uv run python scripts/autoresearch/ar5/cmarepair_build_held4.py --resume"; }

if [ "$dry" = "1" ]; then
  echo "dry run: guards passed (cache $cache, at most $jobs heavy processes)"
  if [ "$resume" = "1" ]; then
    echo "  resume: KEPT (finished, not rerun):${kept:- none}"
    echo "  resume: RERUN from scratch (leftovers moved to $cache/interrupted_<time>/):${redo:- none}"
    [ -e "$marker" ] && echo "    would move $marker (old marker: $(head -1 "$marker"))"
  fi
  if [ "$smoke" = "0" ]; then echo "  $(build_cmd) > $logs/held4_driver.log 2>&1"; fi
  for j in $joblist; do
    a="${j%%:*}"; cap="${j#*:}"
    echo "  $(run_cmd "$a" "$cap") > $logs/${a}_held4.log 2>&1"
    if [ "$smoke" = "0" ]; then echo "    then, if it finished: $(pass_cmd "$a") > $logs/pass_${a}_held4.log 2>&1"; fi
  done
  if [ "$smoke" = "0" ]; then echo "    after J5's pass: $(repair_cmd) > $logs/repair_held4.log 2>&1"; fi
  echo "  then: marker $marker (exit=<code>, reason=<text>)"
  exit 0
fi

mkdir -p "$logs"
# safe to start detached: ignore a hangup, read no terminal
trap '' HUP
exec < /dev/null
if [ "$resume" = "1" ]; then
  stamp="$(date +%Y%m%d-%H%M%S)"
  idir="$cache/interrupted_$stamp"
  if [ -e "$marker" ]; then   # an old marker, failed or not, never blocks a resume
    mkdir -p "$idir/logs"; mv "$marker" "$idir/logs/cmarepair_runs.done"
  fi
  for n in $redo; do
    case "$n" in
      repair_held4) fl="$(repair_files)" ;;
      *) fl="$(jfiles "${n%_held4}")" ;;
    esac
    for p in $fl; do
      if [ -e "$cache/$p" ]; then
        mkdir -p "$idir/$(dirname "$p")"
        mv "$cache/$p" "$idir/$p"
      fi
    done
  done
  echo "resume $(date +%Y-%m-%dT%H:%M:%S%z) kept:${kept:- none} rerun:${redo:- none}" >> "$logs/cmarepair_run_all.resume"
else
  date +%Y-%m-%dT%H:%M:%S%z > "$logs/cmarepair_run_all.start"
fi

# ---- the held4 build (a finished stage is skipped by the builder itself)
if [ "$smoke" = "0" ]; then
  eval "$(build_cmd)" > "$logs/held4_driver.log" 2>&1
  brc=$?
  [ "$brc" = "0" ] && { build_finished || brc=4; }
  echo "$brc" > "$logs/held4_driver.rc"
  if [ "$brc" != "0" ]; then
    date +%Y-%m-%dT%H:%M:%S%z > "$logs/cmarepair_run_all.end"
    end 1 "the held4 build failed (see held4_driver.log; not read for results)"
  fi
fi

run_job() { # run_job <arm> <cap>   (a subshell: run, then pass, then, for J5, the repair)
  a="$1"; cap="$2"
  case " $kept " in *" ${a}_held4 "*) ;; *)
    eval "$(run_cmd "$a" "$cap")" > "$logs/${a}_held4.log" 2>&1
    rc=$?
    echo "$rc" > "$logs/${a}_held4.run_rc"
    d="$cache/$(rundir "$a")"
    if [ "$rc" = "0" ] && [ -e "$d/run.json" ] && [ ! -e "$d/capped.json" ] && [ "$smoke" = "0" ]; then
      eval "$(pass_cmd "$a")" > "$logs/pass_${a}_held4.log" 2>&1
      prc=$?
      echo "$prc" > "$logs/pass_${a}_held4.rc"
    else
      prc=skipped
      echo skipped > "$logs/pass_${a}_held4.rc"
    fi
    echo "run=$rc pass=$prc" > "$logs/${a}_held4.done"
  ;; esac
  if [ "$a" = "J5" ] && [ "$smoke" = "0" ]; then
    case " $kept " in *" repair_held4 "*) ;; *)
      if job_finished J5; then
        eval "$(repair_cmd)" > "$logs/repair_held4.log" 2>&1
        rrc=$?
      else
        rrc=skipped
      fi
      echo "$rrc" > "$logs/repair_held4.rc"
      echo "exit=$rrc" > "$logs/repair_held4.done"
    ;; esac
  fi
}

pids=""
for j in $joblist; do
  a="${j%%:*}"; cap="${j#*:}"
  run_job "$a" "$cap" &
  pids="$pids $!"
done
for p in $pids; do wait "$p"; done

capped=""
fail=0
for j in $joblist; do
  a="${j%%:*}"
  rc="$(cat "$logs/${a}_held4.run_rc" 2>/dev/null || echo 99)"
  prc="$(cat "$logs/pass_${a}_held4.rc" 2>/dev/null || echo 99)"
  if [ "$rc" = "3" ]; then
    capped="${a}_held4"
    [ "$smoke" = "1" ] || fail=1
  elif [ "$rc" != "0" ]; then fail=1
  fi
  if [ "$smoke" = "0" ] && [ "$prc" != "0" ]; then fail=1; fi
done
if [ "$smoke" = "0" ]; then
  [ "$(cat "$logs/repair_held4.rc" 2>/dev/null || echo 99)" != "0" ] && fail=1
fi
reason="all runs, passes and the repair finished"
[ "$smoke" = "1" ] && reason="smoke finished (no pass, no repair)"
date +%Y-%m-%dT%H:%M:%S%z > "$logs/cmarepair_run_all.end"
if [ "$smoke" = "0" ] && [ -n "$capped" ]; then end 1 "STOP RULE: a run reached its cap: $capped"; fi
if [ "$fail" = "1" ]; then
  end 1 "a run, a pass or the repair failed (see the run_rc, pass rc and repair rc files; not read for results)"
fi
end 0 "$reason"
