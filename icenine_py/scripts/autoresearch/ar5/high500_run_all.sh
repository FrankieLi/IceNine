#!/bin/bash
# AR5 step 16, phase A steps 5 and 6: the held3 pilot and build, then the whole-layer runs of arms
# S8, J and J5 on held3 and of arm J5 on held, mid, grad and big, each followed by its boundary
# pass in the same slot.
# A FOREGROUND script: the main session starts it from icenine_py/ (it does not detach; start it
# with nohup ... &, see below).
#
#   cd icenine_py && nohup scripts/autoresearch/ar5/high500_run_all.sh > /dev/null 2>&1 &
#   (macOS has no setsid; nohup and the script's own trap on HUP are enough)
#   ... --dry-run      guards and commands only (nothing written)
#   ... --cache DIR    another cache folder (tests)        ... --jobs N   at most N heavy (<= 4)
#   ... --smoke        test form: the held3 build (REAL: the sample is the one the runs use; a
#                      finished stage is skipped later), then short jobs (--smoke 30, a cap of
#                      --smoke-cap M minutes each, default 20, no pass): J5 on held and on big, J and
#                      S8 on big through high500_launch.py, and the held3 pilot; cache defaults to
#                      benchmarks/autoresearch/ar5/high500/cache_smoke
#   ... --resume       after an interruption or a FAILED earlier attempt: the held3 build skips its
#                      finished stages and moves an unfinished stage's output aside; a job whose
#                      .done says "run=0 pass=0" (a smoke job: "run=0 pass=skipped") with its run.json
#                      is kept; the leftovers of every other job, an old marker (failed or not) and
#                      an old held3 marker are MOVED to <cache>/interrupted_<time>/ and those jobs
#                      run from scratch; the start file is kept, high500_run_all.resume gets a line.
#                      With --dry-run it lists kept / moved and changes nothing.
#
# Order (at most 4 heavy processes at once; a job is a run, then its pass in the same slot):
#   held3 build (high500_make_held3.py: sample, render with 1 worker, stack), in a slot, with J5
#   held, J5 mid, J5 grad beside it; then, when the build has ended: the held3 pilot (J, --smoke
#   30, cap 4 minutes, no pass; its logs and results are not read; a BFS floods the whole layer
#   whatever --smoke says, so the pilot ends at its cap: exit 3 is the pilot's success, a crash is
#   not); S8 held3 (the long pole), J held3 and J5 held3 start only after the pilot ended with exit
#   0 or 3; J5 big last.
# Arms (each followed by its boundary pass):
#   S8  high500_launch.py fit_run.py --name S8 (step 8's S8: stock seed step), then fit_pass.py
#   J   high500_launch.py polish_run.py (arm J of step 14), then order_pass.py --name J --aligned 1
#   J5  high500_run.py (arm J with rule PH(500) in place of PH(3000)), then
#       high500_pass.py --name J5 --aligned 1
# Caps (--cap-min): S8 on held3 480, every other run 300. A run at its cap writes capped.json and
# exits 3; its pass is skipped; no new job is then started and the marker says so (stop rule).
# Not for the pilot, and not for the smoke jobs on held and held3 (their cap is the design: 4 and
# 5 minutes): exit 3 is their normal end.
#
# Guard 0 (exit 2, nothing written, reason on stderr): the real path of an output lies outside the
# cache (stallhit_cache_guard.py). Guards (exit 2, nothing started; a note logs/high500_run_all.
# refused.<time> holds the reason): an output of a job or the final marker exists (not with
# --resume: there an old marker is moved aside); high500_freeze.py --check fails (frozen.json
# committed, clean, hashes agree); a high500, stallhit, rescue, order, fit, patience, band or polish
# run, pass, replay or build process is running; an input is missing (samples, image stacks,
# remedies_pick.json, the stored arm-J runs linked in the cache). held3's inputs are made by the
# build. Nothing of held2 is used: no job, no input, no path (high500_common refuses
# RESEED_SAMPLE=held2). Older steps' freeze checks that cannot agree any more are replaced inside
# high500_common only (see there); step 14's check stays active.
#
# Files, in <cache>/logs (cache = benchmarks/autoresearch/ar5/high500/cache, untracked):
#   <arm>_<sample>.log, .run_rc, pass_<arm>_<sample>.log, pass_<..>.rc ("skipped" or the code),
#   <arm>_<sample>.done ("run=<code> pass=<code|skipped>"), held3_driver.log / .rc and
#   held3_<stage>.log / .rc (the build), pilot_held3.* (the pilot), high500_run_all.start / .end,
#   and the MARKER high500_runs.done: first line "exit=<code>", then "reason=<text>", the end time
#   and one line per job. The script's exit code is the marker's. When every held3 job has ended
#   (any outcome) <cache>/held3_runs_done.json is written (the readers of this step open held3's
#   logs only after it exists; high500_guard). Run folders <cache>/<arm>_<sample>/, passes in
#   <cache>/passes_<sample>/. Wall times are contended and are not timing results. NOBODY reads
#   held3's logs before the marker says all runs ended (plan stop rule).
# Plain bash 3.2 (macOS default).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${HIGH500_CACHE:-$here/benchmarks/autoresearch/ar5/high500/cache}"
dry=0
resume=0
jobs=4
smoke=0
smoke_cap=20
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
  cache="$here/benchmarks/autoresearch/ar5/high500/cache_smoke"
fi
logs="$cache/logs"
marker="$logs/high500_runs.done"
h3marker="$cache/held3_runs_done.json"
export HIGH500_CACHE="$cache"
export POLISH_CACHE="$cache"
export STALLHIT_CACHE="$cache"
export RESCUE_CACHE="$cache"
export ORDER_CACHE="$cache"   # order_pass.py / order_common read this one
export RESEED_CACHE="$cache"  # fit_run.py / fit_pass.py read this one
unset VIRTUAL_ENV

# job list: sample:arm:cap, in start order ("held3:build:0" is the held3 build, "held3:pilot:C" the
# held3 pilot)
if [ "$smoke" = "1" ]; then
  joblist="held3:build:0 held:J5:5 big:J5:$smoke_cap big:J:$smoke_cap big:S8:$smoke_cap held3:pilot:4"
  smoke_arg=" --smoke 30"
else
  joblist="held3:build:0 held:J5:300 mid:J5:300 grad:J5:300 held3:pilot:4 held3:S8:480 held3:J:300 held3:J5:300 big:J5:300"
  smoke_arg=""
fi

jname() { # jname <sample> <arm>  ->  <arm>_<sample>   (the build: held3_driver)
  if [ "$2" = "build" ]; then echo "held3_driver"; else echo "$2_$1"; fi
}
is_short() { # is_short <arm>: the pilot and every smoke job: no pass, a *_smoke run folder
  [ "$smoke" = "1" ] || [ "$1" = "pilot" ]
}
cap_ok() { # cap_ok <sample> <arm>: exit 3 (the cap) is this job's normal end
  [ "$2" = "pilot" ] && return 0
  [ "$smoke" = "1" ] && { [ "$1" = "held" ] || [ "$1" = "held3" ]; } && return 0
  return 1
}
rundir() { # rundir <sample> <arm> <name>: the run folder (relative to the cache)
  if [ "$2" = "pilot" ]; then echo "J_$1_smoke"
  elif is_short "$2"; then echo "$3_smoke"
  else echo "$3"; fi
}

end() { # end <exit> <reason>: write the markers and exit
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
    if [ "$smoke" = "0" ]; then   # every held3 job has ended: its logs may be read
      echo "{\"exit\": $code, \"reason\": \"$reason\"}" > "$h3marker"
    fi
  else
    echo "dry run: would end with exit=$code reason=$reason"
  fi
  exit "$code"
}

refuse() { # a refusal never writes the final marker (that would block the real start)
  echo "refusing to start: $1" >&2
  if [ "$dry" = "0" ]; then
    mkdir -p "$logs"
    echo "exit=2 reason=refused: $1" > "$logs/high500_run_all.refused.$(date +%s)"
  fi
  exit 2
}

# the files of a job (paths relative to the cache)
jfiles() { # jfiles <sample> <arm> <name>
  if [ "$2" = "build" ]; then
    echo "logs/held3_driver.log logs/held3_driver.rc logs/held3_build.log logs/held3_build.rc logs/held3_render.log logs/held3_render.rc logs/held3_stack.log logs/held3_stack.rc"
  elif is_short "$2"; then
    echo "$(rundir "$1" "$2" "$3") logs/$3.log logs/$3.done logs/$3.run_rc logs/pass_$3.log logs/pass_$3.rc"
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
  rd="$(rundir "$1" "$2" "$3")"
  if is_short "$2"; then
    d="$(cat "$logs/$3.done" 2>/dev/null)"
    if [ "$d" = "run=0 pass=skipped" ] && [ -e "$cache/$rd/run.json" ] && [ ! -e "$cache/$rd/capped.json" ]; then
      return 0
    fi
    cap_ok "$1" "$2" && [ "$d" = "run=3 pass=skipped" ] && [ -e "$cache/$rd/capped.json" ]
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
  if [ "$a" != "build" ] && ! is_short "$a"; then outs="$outs passes_$s"; fi
  outs="$outs $(jfiles "$s" "$a" "$n")"
done
if ! uv run python scripts/autoresearch/ar5/stallhit_cache_guard.py --cache "$cache" $outs; then
  echo "refusing to start: an output's real path lies outside the cache $cache" >&2
  exit 2
fi

# ---- guard 1: no output exists (--resume: finished jobs are kept, leftovers and an old marker
# are moved aside later)
if [ -e "$marker" ] && [ "$resume" = "0" ]; then refuse "the marker $marker already exists (use --resume)"; fi
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  if [ "$resume" = "1" ]; then
    if job_finished "$s" "$a" "$n"; then kept="$kept $n"; else redo="$redo $n"; fi
    continue
  fi
  [ "$a" = "build" ] && continue   # the build skips a finished stage and refuses nothing here
  for p in $(jfiles "$s" "$a" "$n"); do
    if [ -e "$cache/$p" ]; then refuse "$cache/$p already exists"; fi
  done
done

# ---- guard 2: the freeze
ferr="$(mktemp)"
if ! uv run python scripts/autoresearch/ar5/high500_freeze.py --check > /dev/null 2>"$ferr"; then
  msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
  refuse "${msg:-high500 freeze check failed}"
fi
rm -f "$ferr"

# ---- guard 3: nothing of this step (or an arm it builds on) is running
if pgrep -f "autoresearch/ar5/(high500|band|order|rescue|patience|stallhit|fit|polish)_(run|pass|replay|record|paths|launch|make_held3|make_held2|probe)\.py" > /dev/null 2>&1; then
  refuse "a high500 / band / stallhit / rescue / order / fit / polish run, replay, recording, pass or build process is already running"
fi

# ---- guard 4: inputs (link the untracked ones first; links only)
if [ "$dry" = "1" ]; then
  scripts/autoresearch/ar5/high500_link_cache.sh --dry-run --cache "$cache" big grad mid held > /dev/null 2>&1 || refuse "a source of the cache links is missing (high500_link_cache.sh --dry-run)"
else
  scripts/autoresearch/ar5/high500_link_cache.sh --cache "$cache" big grad mid held > /dev/null 2>&1 || refuse "a source of the cache links is missing (high500_link_cache.sh)"
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
if [ "$smoke" = "0" ]; then   # the linked stored arm-J runs (existence only)
  for s in big grad mid held; do
    need="$need $cache/J_$s/run.json $cache/passes_$s/after_J.npy $cache/passes_$s/D4_J/pass.json"
  done
fi
for p in $need; do
  if [ ! -e "$p" ]; then refuse "input missing: $p"; fi
done

run_cmd() { # run_cmd <sample> <arm> <cap>: the run command (text, for the dry run and the log)
  case "$2" in
    S8) echo "RESEED_SAMPLE=$1 HIGH500_CACHE=$cache uv run python scripts/autoresearch/ar5/high500_launch.py fit_run.py --name S8 --cap-min $3$smoke_arg" ;;
    J) echo "RESEED_SAMPLE=$1 HIGH500_CACHE=$cache uv run python scripts/autoresearch/ar5/high500_launch.py polish_run.py --cap-min $3$smoke_arg" ;;
    pilot) echo "RESEED_SAMPLE=$1 HIGH500_CACHE=$cache uv run python scripts/autoresearch/ar5/high500_launch.py polish_run.py --cap-min $3 --smoke 30" ;;
    J5) echo "RESEED_SAMPLE=$1 HIGH500_CACHE=$cache uv run python scripts/autoresearch/ar5/high500_run.py --cap-min $3$smoke_arg" ;;
  esac
}
pass_cmd() { # pass_cmd <sample> <arm>
  case "$2" in
    S8) echo "RESEED_SAMPLE=$1 HIGH500_CACHE=$cache uv run python scripts/autoresearch/ar5/high500_launch.py fit_pass.py --name S8" ;;
    J) echo "RESEED_SAMPLE=$1 HIGH500_CACHE=$cache uv run python scripts/autoresearch/ar5/high500_launch.py order_pass.py --name J --aligned 1" ;;
    J5) echo "RESEED_SAMPLE=$1 HIGH500_CACHE=$cache uv run python scripts/autoresearch/ar5/high500_pass.py --name J5 --aligned 1" ;;
  esac
}
build_cmd() { # build_cmd: the held3 build (text); always --resume: a finished stage is skipped
  echo "HIGH500_CACHE=$cache uv run python scripts/autoresearch/ar5/high500_make_held3.py --resume"
}

if [ "$dry" = "1" ]; then
  echo "dry run: guards passed (cache $cache, at most $jobs heavy processes)"
  if [ "$resume" = "1" ]; then
    echo "  resume: KEPT (finished, not rerun):${kept:- none}"
    echo "  resume: RERUN from scratch (leftovers moved to $cache/interrupted_<time>/):${redo:- none}"
    [ -e "$marker" ] && echo "    would move $marker (old marker: $(head -1 "$marker"))"
    [ -e "$h3marker" ] && echo "    would move $h3marker"
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
    if ! is_short "$a"; then echo "    then, if it finished: $(pass_cmd "$s" "$a") > $logs/pass_$n.log 2>&1"; fi
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
  idir="$cache/interrupted_$stamp"
  for f in "$marker" "$h3marker"; do   # an old marker, failed or not, never blocks a resume
    if [ -e "$f" ]; then
      mkdir -p "$idir/$(basename "$(dirname "$f")")"
      mv "$f" "$idir/$(basename "$(dirname "$f")")/$(basename "$f")"
    fi
  done
  for n in $redo; do
    for j in $joblist; do
      s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"
      [ "$(jname "$s" "$a")" = "$n" ] || continue
      for p in $(jfiles "$s" "$a" "$n"); do
        if [ -e "$cache/$p" ]; then
          mkdir -p "$idir/$(dirname "$p")"
          mv "$cache/$p" "$idir/$p"
        fi
      done
    done
  done
  echo "resume $(date +%Y-%m-%dT%H:%M:%S%z) kept:${kept:- none} rerun:${redo:- none}" >> "$logs/high500_run_all.resume"
else
  date +%Y-%m-%dT%H:%M:%S%z > "$logs/high500_run_all.start"
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
  d="$cache/$(rundir "$s" "$a" "$n")"
  if [ "$rc" = "0" ] && [ -e "$d/run.json" ] && [ ! -e "$d/capped.json" ] && ! is_short "$a"; then
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
    cap_ok "$ks" "$ka" && continue
    if [ "$(cat "$logs/$(jname "$ks" "$ka").run_rc" 2>/dev/null)" = "3" ]; then
      echo "$(jname "$ks" "$ka")"; return
    fi
  done
}

capped=""
build_failed=""
pilot_failed=""
pids=""
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; cap="${rest#*:}"; n="$(jname "$s" "$a")"
  case " $kept " in *" $n "*) continue ;; esac
  # a held3 run waits for the build to end (and to have ended well)
  if [ "$s" = "held3" ] && [ "$a" != "build" ]; then
    while [ ! -e "$logs/held3_driver.rc" ]; do sleep 20; done
    if [ "$(cat "$logs/held3_driver.rc")" != "0" ]; then build_failed=1; break; fi
  fi
  # S8, J and J5 on held3 wait for the pilot to have ended with exit 0 (a kept pilot counts)
  if [ "$s" = "held3" ] && [ "$a" != "build" ] && [ "$a" != "pilot" ] && [ "$smoke" = "0" ]; then
    while [ ! -e "$logs/pilot_held3.done" ]; do sleep 20; done
    pd="$(cat "$logs/pilot_held3.done")"
    if [ "$pd" != "run=0 pass=skipped" ] && [ "$pd" != "run=3 pass=skipped" ]; then pilot_failed=1; break; fi
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
  if [ "$rc" != "0" ]; then
    if [ "$rc" = "3" ] && cap_ok "$s" "$a"; then :; else fail=1; fi
  fi
  if [ "$prc" != "0" ] && ! is_short "$a"; then fail=1; fi
done
[ "$smoke" = "1" ] && reason="smoke finished (no pass)"
[ -n "$build_failed" ] && reason="the held3 build failed (see held3_driver.log; not read for results)"
[ -n "$pilot_failed" ] && reason="the held3 pilot failed (no held3 run was started)"
[ -n "$capped" ] && reason="STOP RULE: a run reached its cap: $capped"
date +%Y-%m-%dT%H:%M:%S%z > "$logs/high500_run_all.end"
if [ "$fail" = "1" ]; then
  if [ "$reason" = "all runs and passes finished" ]; then reason="a run or a pass failed (see the run_rc and pass rc files)"; fi
  end 1 "$reason"
fi
end 0 "$reason"
