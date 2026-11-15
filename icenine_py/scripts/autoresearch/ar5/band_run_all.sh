#!/bin/bash
# UNTESTED end to end: the gate failed, no form was chosen, so the guards stop it before any job
# (D-20261114-N7: no held3, no whole layer). Only its guard and its parts on big were tried.
# AR5 step 15, phase A steps 6 and 7: build held3, then the whole-layer runs of arms S8, J and JB on
# held3 and of arm JB (the form the gate chose) on held, mid, grad and big, each followed by its
# boundary pass in the same slot.
# A FOREGROUND script: the main session starts it exactly once, from icenine_py/ (it does not
# detach; the caller keeps it alive, e.g. with nohup setsid).
#
#   cd icenine_py && scripts/autoresearch/ar5/band_run_all.sh             # the real thing
#   cd icenine_py && scripts/autoresearch/ar5/band_run_all.sh --dry-run   # guards and commands only
#   ... --cache DIR     another cache folder (tests)      ... --jobs N   at most N heavy (<= 4)
#   ... --smoke-cap M   test form: four short jobs (--smoke 30, a cap of M minutes each, no pass,
#                       no held3): JB on held and big, J and S8 on big through band_launch.py;
#                       the cache defaults to benchmarks/autoresearch/ar5/band/cache_smoke
#   ... --resume        after an interruption (the script was killed): the held3 build skips its
#                       finished stages; a job whose .done says "run=0 pass=0" with its run.json and
#                       pass outputs is kept; the leftovers of every other job are MOVED to
#                       <cache>/interrupted_<time>/ and the job is run from scratch;
#                       band_run_all.start is kept, band_run_all.resume gets a line; the marker
#                       covers all jobs. With --dry-run it lists kept / moved and changes nothing.
#                       Start it detached: nohup setsid scripts/autoresearch/ar5/band_run_all.sh --resume
#
# Order (at most 4 heavy processes at once; a job is a run, then its pass in the same slot):
#   held3 build (band_make_held3.py: sample, render with 1 worker, stack), in a slot, with
#   JB held, JB mid, JB grad beside it; then, when the build has ended and a slot is free,
#   S8 held3 (the long pole), J held3, JB held3; JB big last.
# Arms (each followed by its boundary pass):
#   S8  band_launch.py fit_run.py --name S8 (step 8's S8: stock seed step), then fit_pass.py --name S8
#   J   band_launch.py polish_run.py (arm J of step 14: arm H of step 12 with rule G2b), then
#       order_pass.py --name J --aligned 1
#   JB  band_run.py --form <chosen> (arm J with rule PB(k) in place of PH(3000)), then
#       band_pass.py --name JB --aligned 1
# <chosen> is read from ar5/band/chosen.json (committed by the gate step; band_chosen.py).
# Caps (--cap-min): S8 on held3 480, every other run 300. A run at its cap writes capped.json and
# exits 3; its pass is skipped; no new job is then started and the marker says so (stop rule).
#
# Guard 0 (exit 2, nothing written, reason on stderr): the real path of an output lies outside the
# cache (stallhit_cache_guard.py). Guards (exit 2, nothing started; a note logs/band_run_all.
# refused.<time> holds the reason): an output of a job or the final marker exists (not with
# --resume); band_freeze.py --check fails (frozen.json committed, clean, hashes agree) or
# chosen.json is not committed; a band, stallhit, rescue, order, fit, patience or polish run, pass,
# replay or build process is running; an input is missing (samples, image stacks, remedies_pick.json,
# the stored arm-J runs linked in the cache). held3's inputs are made by the build. Nothing of
# held2 is used: no job, no input, no path (band_common refuses RESEED_SAMPLE=held2).
#
# Files, in <cache>/logs (cache = benchmarks/autoresearch/ar5/band/cache, untracked):
#   <arm>_<sample>.log, .run_rc, pass_<arm>_<sample>.log, pass_<..>.rc ("skipped" or the code),
#   <arm>_<sample>.done ("run=<code> pass=<code|skipped>"), held3_driver.log / .rc and
#   held3_<stage>.log / .rc (the build), band_run_all.start / .end, and the MARKER band_runs.done:
#   first line "exit=<code>", then "reason=<text>", the end time and one line per job. The script's
#   exit code is the marker's. Run folders <cache>/<arm>_<sample>/, passes in <cache>/passes_<sample>/.
#   Wall times are contended and are not timing results. NOBODY reads held3's logs before the
#   marker says all runs ended (plan stop rule).
# Plain bash 3.2 (macOS default).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${BAND_CACHE:-$here/benchmarks/autoresearch/ar5/band/cache}"
dry=0
resume=0
jobs=4
smoke_cap=0
cache_given=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1"; cache_given=1 ;;
    --jobs) shift; jobs="$1" ;;
    --smoke-cap) shift; smoke_cap="$1" ;;
    --resume) resume=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
if [ "$jobs" -gt 4 ] || [ "$jobs" -lt 1 ]; then echo "--jobs must be 1 to 4" >&2; exit 2; fi
if [ "$smoke_cap" != "0" ] && [ "$cache_given" = "0" ]; then
  cache="$here/benchmarks/autoresearch/ar5/band/cache_smoke"
fi
logs="$cache/logs"
marker="$logs/band_runs.done"
export BAND_CACHE="$cache"
export POLISH_CACHE="$cache"
export STALLHIT_CACHE="$cache"
export RESCUE_CACHE="$cache"
export ORDER_CACHE="$cache"   # order_pass.py / order_common read this one
export RESEED_CACHE="$cache"  # fit_run.py / fit_pass.py read this one
unset VIRTUAL_ENV

# job list: sample:arm:cap, in start order ("held3:build:0" is the held3 build)
if [ "$smoke_cap" != "0" ]; then
  joblist="held:JB:$smoke_cap big:JB:$smoke_cap big:J:$smoke_cap big:S8:$smoke_cap"
  smoke_arg=" --smoke 30"
  have_build=0
else
  joblist="held3:build:0 held:JB:300 mid:JB:300 grad:JB:300 held3:S8:480 held3:J:300 held3:JB:300 big:JB:300"
  smoke_arg=""
  have_build=1
fi

jname() { # jname <sample> <arm>  ->  <arm>_<sample>   (the build: held3_driver)
  if [ "$2" = "build" ]; then echo "held3_driver"; else echo "$2_$1"; fi
}

end() { # end <exit> <reason>: write the marker and exit
  code="$1"; reason="$2"
  if [ "$dry" = "0" ]; then
    mkdir -p "$logs"
    {
      echo "exit=$code"
      echo "reason=$reason"
      date +%Y-%m-%dT%H:%M:%S%z
      for j in $joblist; do
        s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
        if [ "$a" = "build" ]; then
          echo "held3_build driver=$(cat "$logs/held3_driver.rc" 2>/dev/null || echo none) build=$(cat "$logs/held3_build.rc" 2>/dev/null || echo none) render=$(cat "$logs/held3_render.rc" 2>/dev/null || echo none) stack=$(cat "$logs/held3_stack.rc" 2>/dev/null || echo none)"
        else
          echo "$n run=$(cat "$logs/$n.run_rc" 2>/dev/null || echo none) pass=$(cat "$logs/pass_$n.rc" 2>/dev/null || echo none)"
        fi
      done
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
    echo "exit=2 reason=refused: $1" > "$logs/band_run_all.refused.$(date +%s)"
  fi
  exit 2
}

# the files of a job (paths relative to the cache)
jfiles() { # jfiles <sample> <arm> <name>
  if [ "$2" = "build" ]; then
    echo "logs/held3_driver.log logs/held3_driver.rc logs/held3_build.log logs/held3_build.rc logs/held3_render.log logs/held3_render.rc logs/held3_stack.log logs/held3_stack.rc"
  else
    echo "$3 passes_$1/D4_$2 passes_$1/before_$2.npy passes_$1/after_$2.npy logs/$3.log logs/$3.done logs/$3.run_rc logs/pass_$3.log logs/pass_$3.rc"
  fi
}
job_finished() { # job_finished <sample> <arm> <name>: the outputs say it ended well
  if [ "$2" = "build" ]; then
    [ "$(cat "$logs/held3_driver.rc" 2>/dev/null)" = "0" ] \
      && [ -e "$here/benchmarks/autoresearch/ar3/held3/gradient/sample/gradient_grid.mic" ] \
      && [ -e "$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/held3_clean.npy" ]
    return
  fi
  [ "$(cat "$logs/$3.done" 2>/dev/null)" = "run=0 pass=0" ] && [ -e "$cache/$3/run.json" ] \
    && [ ! -e "$cache/$3/capped.json" ] && [ -e "$cache/passes_$1/after_$2.npy" ] \
    && [ -e "$cache/passes_$1/D4_$2/pass.json" ] && [ -e "$cache/passes_$1/before_$2.npy" ]
}
kept=""   # names of the jobs kept by --resume
redo=""   # names of the jobs run (again) with leftovers moved

# ---- guard 0: every output (logs, run folders, pass folders and files) must really lie in this
# cache, not behind a link into another worktree's cache. Runs first, also in a dry run.
outs="logs"
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  if [ "$a" != "build" ]; then outs="$outs passes_$s"; fi
  outs="$outs $(jfiles "$s" "$a" "$n")"
done
if ! uv run python scripts/autoresearch/ar5/stallhit_cache_guard.py --cache "$cache" $outs; then
  echo "refusing to start: an output's real path lies outside the cache $cache" >&2
  exit 2
fi

# ---- guard 1: no output exists (--resume: finished jobs are kept, leftovers are moved later)
if [ -e "$marker" ]; then refuse "the marker $marker already exists"; fi
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  if [ "$resume" = "1" ]; then
    if job_finished "$s" "$a" "$n"; then kept="$kept $n"; else redo="$redo $n"; fi
    continue
  fi
  for p in $(jfiles "$s" "$a" "$n"); do
    if [ -e "$cache/$p" ]; then refuse "$cache/$p already exists"; fi
  done
done

# ---- guard 2: the freeze and the chosen form
ferr="$(mktemp)"
if ! uv run python scripts/autoresearch/ar5/band_freeze.py --check > /dev/null 2>"$ferr"; then
  msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
  refuse "${msg:-band freeze check failed}"
fi
form="$(uv run python scripts/autoresearch/ar5/band_chosen.py 2>"$ferr")"
if [ $? -ne 0 ] || [ -z "$form" ]; then
  msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
  refuse "${msg:-no chosen form}"
fi
rm -f "$ferr"
case "$form" in PB1000|PB500) ;; *) refuse "unknown chosen form: $form" ;; esac

# ---- guard 3: nothing of this step (or an arm it builds on) is running
if pgrep -f "autoresearch/ar5/(band|order|rescue|patience|stallhit|fit|polish)_(run|pass|replay|record|paths|launch|make_held3|make_held2|probe)\.py" > /dev/null 2>&1; then
  refuse "a band / stallhit / rescue / order / fit / polish run, replay, recording, pass or build process is already running"
fi

# ---- guard 4: inputs (link the untracked ones first; links only)
if [ "$dry" = "1" ]; then
  scripts/autoresearch/ar5/order_link_inputs.sh --dry-run > /dev/null 2>&1 || refuse "a source of the input links is missing (order_link_inputs.sh --dry-run)"
  scripts/autoresearch/ar5/band_link_cache.sh --dry-run --cache "$cache" big grad mid held > /dev/null 2>&1 || refuse "a source of the cache links is missing (band_link_cache.sh --dry-run)"
else
  scripts/autoresearch/ar5/order_link_inputs.sh > /dev/null 2>&1 || refuse "a source of the input links is missing (order_link_inputs.sh)"
  scripts/autoresearch/ar5/band_link_cache.sh --cache "$cache" big grad mid held > /dev/null 2>&1 || refuse "a source of the cache links is missing (band_link_cache.sh)"
fi
st="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
need="$here/benchmarks/autoresearch/ar3/reseed/remedies_pick.json"
for j in $joblist; do
  s="${j%%:*}"
  case "$s" in
    big) need="$need $here/benchmarks/autoresearch/ar3/big/gradient/sample/gradient_grid.mic $st/big_clean.npy" ;;
    mid) need="$need $here/benchmarks/autoresearch/ar3/mid/gradient/sample/gradient_grid.mic $st/mid_clean.npy" ;;
    grad) need="$need $here/benchmarks/autoresearch/ar3/gradient/sample/gradient_grid.mic $st/grad_clean.npy" ;;
    held) need="$need $here/benchmarks/autoresearch/ar3/held/gradient/sample/gradient_grid.mic $st/held_clean.npy" ;;
    held3) ;;   # made by the build
  esac
done
if [ "$smoke_cap" = "0" ]; then   # the linked stored arm-J runs (existence only)
  for s in big grad mid held; do
    need="$need $cache/J_$s/run.json $cache/passes_$s/after_J.npy $cache/passes_$s/D4_J/pass.json"
  done
fi
for p in $need; do
  if [ ! -e "$p" ]; then refuse "input missing: $p"; fi
done

run_cmd() { # run_cmd <sample> <arm> <cap>: the run command (text, for the dry run and the log)
  case "$2" in
    S8) echo "RESEED_SAMPLE=$1 BAND_CACHE=$cache uv run python scripts/autoresearch/ar5/band_launch.py fit_run.py --name S8 --cap-min $3$smoke_arg" ;;
    J) echo "RESEED_SAMPLE=$1 BAND_CACHE=$cache uv run python scripts/autoresearch/ar5/band_launch.py polish_run.py --cap-min $3$smoke_arg" ;;
    JB) echo "RESEED_SAMPLE=$1 BAND_CACHE=$cache uv run python scripts/autoresearch/ar5/band_run.py --form $form --cap-min $3$smoke_arg" ;;
  esac
}
pass_cmd() { # pass_cmd <sample> <arm>
  case "$2" in
    S8) echo "RESEED_SAMPLE=$1 BAND_CACHE=$cache uv run python scripts/autoresearch/ar5/band_launch.py fit_pass.py --name S8" ;;
    J) echo "RESEED_SAMPLE=$1 BAND_CACHE=$cache uv run python scripts/autoresearch/ar5/band_launch.py order_pass.py --name J --aligned 1" ;;
    JB) echo "RESEED_SAMPLE=$1 BAND_CACHE=$cache uv run python scripts/autoresearch/ar5/band_pass.py --name JB --aligned 1" ;;
  esac
}
build_cmd() { # build_cmd: the held3 build (text)
  local r=""; [ "$resume" = "1" ] && r=" --resume"
  echo "BAND_CACHE=$cache uv run python scripts/autoresearch/ar5/band_make_held3.py$r"
}

if [ "$dry" = "1" ]; then
  echo "dry run: guards passed (cache $cache, chosen form $form, at most $jobs heavy processes)"
  if [ "$resume" = "1" ]; then
    echo "  resume: KEPT (finished, not rerun):${kept:- none}"
    echo "  resume: RERUN from scratch (leftovers moved to $cache/interrupted_<time>/):${redo:- none}"
    for n in $redo; do
      for j in $joblist; do
        s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"
        [ "$(jname "$s" "$a")" = "$n" ] || continue
        for p in $(jfiles "$s" "$a" "$n"); do
          [ -e "$cache/$p" ] && echo "    would move $p"
        done
      done
    done
  fi
  for j in $joblist; do
    s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; cap="${rest#*:}"; n="$(jname "$s" "$a")"
    case " $kept " in *" $n "*) continue ;; esac
    if [ "$a" = "build" ]; then
      echo "  $(build_cmd) > $logs/held3_driver.log 2>&1"
      continue
    fi
    echo "  $(run_cmd "$s" "$a" "$cap") > $logs/$n.log 2>&1"
    [ "$smoke_cap" = "0" ] && echo "    then, if it finished: $(pass_cmd "$s" "$a") > $logs/pass_$n.log 2>&1"
  done
  echo "  then: marker $marker (exit=<code>, reason=<text>)"
  exit 0
fi

mkdir -p "$logs"
# safe to start detached (setsid + nohup): ignore a hangup, read no terminal, and let caffeinate
# wait on this script's pid only (it ends with it)
trap '' HUP
exec < /dev/null
caffeinate -i -w $$ > /dev/null 2>&1 &
if [ "$resume" = "1" ]; then
  stamp="$(date +%Y%m%d-%H%M%S)"
  for n in $redo; do
    for j in $joblist; do
      s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"
      [ "$(jname "$s" "$a")" = "$n" ] || continue
      for p in $(jfiles "$s" "$a" "$n"); do
        if [ -e "$cache/$p" ]; then
          mkdir -p "$cache/interrupted_$stamp/$(dirname "$p")"
          mv "$cache/$p" "$cache/interrupted_$stamp/$p"
        fi
      done
    done
  done
  echo "resume $(date +%Y-%m-%dT%H:%M:%S%z) kept:${kept:- none} rerun:${redo:- none}" >> "$logs/band_run_all.resume"
else
  date +%Y-%m-%dT%H:%M:%S%z > "$logs/band_run_all.start"
fi

build_job() { # a subshell: the held3 build, then its rc ("0" only if the inputs exist)
  eval "$(build_cmd)" > "$logs/held3_driver.log" 2>&1
  brc=$?
  if [ "$brc" = "0" ]; then
    for p in "$here/benchmarks/autoresearch/ar3/held3/gradient/sample/gradient_grid.mic" "$st/held3_clean.npy"; do
      [ -e "$p" ] || brc=4
    done
  fi
  echo "$brc" > "$logs/held3_driver.rc"
}

run_job() { # run_job <sample> <arm> <cap>   (a subshell: run, then pass)
  s="$1"; a="$2"; cap="$3"; n="$(jname "$s" "$a")"
  eval "$(run_cmd "$s" "$a" "$cap")" > "$logs/$n.log" 2>&1
  rc=$?
  echo "$rc" > "$logs/$n.run_rc"
  d="$cache/$n${smoke_arg:+_smoke}"
  if [ "$rc" = "0" ] && [ -e "$d/run.json" ] && [ ! -e "$d/capped.json" ] && [ -z "$smoke_arg" ]; then
    eval "$(pass_cmd "$s" "$a")" > "$logs/pass_$n.log" 2>&1
    prc=$?
    echo "$prc" > "$logs/pass_$n.rc"
  else
    prc=skipped
    echo skipped > "$logs/pass_$n.rc"
  fi
  echo "run=$rc pass=$prc" > "$logs/$n.done"
}

any_capped() { # prints the name of a run that reached its cap (exit 3), if any
  for k in $joblist; do
    ks="${k%%:*}"; kr="${k#*:}"; ka="${kr%%:*}"
    [ "$ka" = "build" ] && continue
    if [ "$(cat "$logs/$(jname "$ks" "$ka").run_rc" 2>/dev/null)" = "3" ]; then
      echo "$(jname "$ks" "$ka")"; return
    fi
  done
}

capped=""
build_failed=""
pids=""
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; cap="${rest#*:}"; n="$(jname "$s" "$a")"
  case " $kept " in *" $n "*) continue ;; esac
  # a held3 run waits for the build to end (and to have ended well)
  if [ "$s" = "held3" ] && [ "$a" != "build" ] && [ "$have_build" = "1" ]; then
    while [ ! -e "$logs/held3_driver.rc" ]; do sleep 20; done
    if [ "$(cat "$logs/held3_driver.rc")" != "0" ]; then build_failed=1; break; fi
  fi
  # wait for a free slot
  while :; do
    alive=""; cnt=0
    for p in $pids; do
      if kill -0 "$p" 2>/dev/null; then alive="$alive $p"; cnt=$((cnt + 1)); fi
    done
    pids="$alive"
    [ "$cnt" -lt "$jobs" ] && break
    sleep 20
  done
  # stop rule: after a run hit its cap, no new job is started
  capped="$(any_capped)"
  [ -n "$capped" ] && break
  if [ "$a" = "build" ]; then
    build_job &
  else
    run_job "$s" "$a" "$cap" &
  fi
  pids="$pids $!"
done
for p in $pids; do wait "$p"; done   # the jobs only (not caffeinate)
[ -z "$capped" ] && capped="$(any_capped)"

fail=0
reason="all runs and passes finished"
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  if [ "$a" = "build" ]; then
    [ "$(cat "$logs/held3_driver.rc" 2>/dev/null || echo 99)" != "0" ] && fail=1
    continue
  fi
  rc=$(cat "$logs/$n.run_rc" 2>/dev/null || echo 99)
  prc=$(cat "$logs/pass_$n.rc" 2>/dev/null || echo 99)
  if [ "$rc" != "0" ]; then fail=1; fi
  if [ "$prc" != "0" ] && [ -z "$smoke_arg" ]; then fail=1; fi
done
[ -n "$smoke_arg" ] && reason="smoke finished (no pass)"
[ -n "$build_failed" ] && reason="the held3 build failed (see held3_driver.log; not read for results)"
[ -n "$capped" ] && reason="STOP RULE: a run reached its cap: $capped"
date +%Y-%m-%dT%H:%M:%S%z > "$logs/band_run_all.end"
if [ "$fail" = "1" ]; then
  if [ "$reason" = "all runs and passes finished" ]; then reason="a run or a pass failed (see the run_rc and pass rc files)"; fi
  end 1 "$reason"
fi
end 0 "$reason"
