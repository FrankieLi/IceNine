#!/bin/bash
# AR5 step 14, phase A step 6: the whole-layer runs of arm J (arm H with rule G2, the frozen form) on
# held2, held, mid, grad and big, each followed by its boundary pass in the same slot.
# A FOREGROUND script: the main session starts it exactly once, from icenine_py/ (it does not
# detach; the caller keeps it alive, e.g. with nohup setsid).
#
#   cd icenine_py && scripts/autoresearch/ar5/polish_run_all.sh             # the real thing
#   cd icenine_py && scripts/autoresearch/ar5/polish_run_all.sh --dry-run   # guards and commands only
#   ... --cache DIR     another cache folder (tests)      ... --jobs N   at most N heavy (<= 4)
#   ... --smoke-cap M   test form: one job, big, with --smoke 30 and a cap of M minutes
#   ... --resume        after an interruption (the script was killed): a job whose .done says
#                       "run=0 pass=0" with its run.json and pass outputs is kept; the leftovers of
#                       every other job are MOVED to <cache>/interrupted_<time>/ and the job is run
#                       from scratch; polish_run_all.start is kept, polish_run_all.resume gets a
#                       line; the marker covers all jobs. With --dry-run it lists kept / moved and
#                       changes nothing. Start it detached:
#                       nohup setsid scripts/autoresearch/ar5/polish_run_all.sh --resume
#
# Order (at most 4 heavy processes at once; a job is a run, then its pass in the same slot):
#   held2, held, mid, grad, big (longest first; the fifth starts as a slot frees)
#   J   polish_run.py (arm H with rule G2), then polish_pass.py --name J --aligned 1
# Caps (--cap-min): 300 minutes a run. A run at its cap writes capped.json and exits 3; its pass is
# skipped; no new job is then started and the marker says so (stop rule of the plan).
#
# Guard 0 (exit 2, nothing written, reason on stderr): the real path of an output lies outside the
# cache (stallhit_cache_guard.py). Guards (exit 2, nothing started; a note logs/polish_run_all.
# refused.<time> holds the reason): an output of a job or the final marker exists (not with
# --resume); polish_freeze.py --check fails (frozen.json committed, clean, hashes agree); a
# stallhit, rescue, order, fit, patience or polish run, pass or replay process is running; an input
# is missing (samples and image stacks, remedies_pick.json). held2 is read only after the freeze
# commit (polish_common). Nothing under step 13's passes_held2 is opened: this script's own cache
# holds its own passes_held2/D4_J outputs only, which it writes and does not read.
#
# Files, in <cache>/logs (cache = benchmarks/autoresearch/ar5/polish/cache, untracked):
#   J_<sample>.log, .run_rc, pass_J_<sample>.log, pass_J_<sample>.rc ("skipped" or the code),
#   J_<sample>.done ("run=<code> pass=<code|skipped>"), polish_run_all.start / .end, and the MARKER
#   polish_runs.done: first line "exit=<code>", then "reason=<text>", the end time and one line per
#   job. The script's exit code is the marker's. Run folders <cache>/J_<sample>/, passes in
#   <cache>/passes_<sample>/. Wall times are contended and are not timing results.
# Plain bash 3.2 (macOS default).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${POLISH_CACHE:-$here/benchmarks/autoresearch/ar5/polish/cache}"
dry=0
resume=0
jobs=4
smoke_cap=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    --jobs) shift; jobs="$1" ;;
    --smoke-cap) shift; smoke_cap="$1" ;;
    --resume) resume=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
if [ "$jobs" -gt 4 ] || [ "$jobs" -lt 1 ]; then echo "--jobs must be 1 to 4" >&2; exit 2; fi
logs="$cache/logs"
marker="$logs/polish_runs.done"
export POLISH_CACHE="$cache"
export STALLHIT_CACHE="$cache"
export RESCUE_CACHE="$cache"
export ORDER_CACHE="$cache"   # order_pass.py reads this one
export RESEED_CACHE="$cache"
unset VIRTUAL_ENV

if [ "$smoke_cap" != "0" ]; then
  joblist="big:J:$smoke_cap"
  smoke_arg=" --smoke 30"
else
  joblist="held2:J:300 held:J:300 mid:J:300 grad:J:300 big:J:300"
  smoke_arg=""
fi

jname() { echo "$2_$1"; }   # jname <sample> <arm>  ->  <arm>_<sample>

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
        echo "$n run=$(cat "$logs/$n.run_rc" 2>/dev/null || echo none) pass=$(cat "$logs/pass_$n.rc" 2>/dev/null || echo none)"
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
    echo "exit=2 reason=refused: $1" > "$logs/polish_run_all.refused.$(date +%s)"
  fi
  exit 2
}

jfiles() { # jfiles <sample> <arm> <name>: the files of a job (paths relative to the cache)
  echo "$3 passes_$1/D4_$2 passes_$1/before_$2.npy passes_$1/after_$2.npy logs/$3.log logs/$3.done logs/$3.run_rc logs/pass_$3.log logs/pass_$3.rc"
}
job_finished() { # job_finished <sample> <arm> <name>: .done says run=0 pass=0 and the outputs exist
  [ "$(cat "$logs/$3.done" 2>/dev/null)" = "run=0 pass=0" ] && [ -e "$cache/$3/run.json" ] \
    && [ ! -e "$cache/$3/capped.json" ] && [ -e "$cache/passes_$1/after_$2.npy" ] \
    && [ -e "$cache/passes_$1/D4_$2/pass.json" ] && [ -e "$cache/passes_$1/before_$2.npy" ]
}
kept=""
redo=""

# ---- guard 0: every output must really lie in this cache (not behind a link into another cache)
outs="logs"
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  outs="$outs passes_$s $(jfiles "$s" "$a" "$n")"
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

# ---- guard 2: the freeze (rule G2b and the scripts that define it)
ferr="$(mktemp)"
if ! uv run python scripts/autoresearch/ar5/polish_freeze.py --check > /dev/null 2>"$ferr"; then
  msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
  refuse "${msg:-polish freeze check failed}"
fi
rm -f "$ferr"

# ---- guard 3: nothing of this step (or an arm it builds on) is running
if pgrep -f "autoresearch/ar5/(order|rescue|patience|stallhit|fit|polish)_(run|pass|replay|probe|make_held2)\.py" > /dev/null 2>&1; then
  refuse "a polish / stallhit / rescue / order / fit run, replay, pass or probe process is already running"
fi

# ---- guard 4: inputs (existence only)
st="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
need="$here/benchmarks/autoresearch/ar3/reseed/remedies_pick.json"
for j in $joblist; do
  s="${j%%:*}"
  case "$s" in
    big|mid|held|held2) need="$need $here/benchmarks/autoresearch/ar3/$s/gradient/sample/gradient_grid.mic $here/benchmarks/autoresearch/ar3/$s/gradient/sample/gradient_truth.mic $st/${s}_clean.npy" ;;
    grad) need="$need $here/benchmarks/autoresearch/ar3/gradient/sample/gradient_grid.mic $st/grad_clean.npy" ;;
  esac
done
for p in $need; do
  if [ ! -e "$p" ]; then refuse "input missing: $p"; fi
done

run_cmd() { # run_cmd <sample> <arm> <cap>: the run command (text, for the dry run and the log)
  echo "RESEED_SAMPLE=$1 POLISH_CACHE=$cache uv run python scripts/autoresearch/ar5/polish_run.py --cap-min $3$smoke_arg"
}
pass_cmd() { # pass_cmd <sample> <arm>
  echo "RESEED_SAMPLE=$1 POLISH_CACHE=$cache uv run python scripts/autoresearch/ar5/polish_pass.py --name $2 --aligned 1"
}

if [ "$dry" = "1" ]; then
  echo "dry run: guards passed (cache $cache, at most $jobs heavy processes)"
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
    echo "  $(run_cmd "$s" "$a" "$cap") > $logs/$n.log 2>&1"
    echo "    then, if it finished: $(pass_cmd "$s" "$a") > $logs/pass_$n.log 2>&1"
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
  echo "resume $(date +%Y-%m-%dT%H:%M:%S%z) kept:${kept:- none} rerun:${redo:- none}" >> "$logs/polish_run_all.resume"
else
  date +%Y-%m-%dT%H:%M:%S%z > "$logs/polish_run_all.start"
fi

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

capped=""
pids=""
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; cap="${rest#*:}"
  case " $kept " in *" $(jname "$s" "$a") "*) continue ;; esac
  while :; do   # wait for a free slot
    alive=""; cnt=0
    for p in $pids; do
      if kill -0 "$p" 2>/dev/null; then alive="$alive $p"; cnt=$((cnt + 1)); fi
    done
    pids="$alive"
    [ "$cnt" -lt "$jobs" ] && break
    sleep 20
  done
  # stop rule: after a run hit its cap, no new job is started
  for k in $joblist; do
    ks="${k%%:*}"; kr="${k#*:}"; ka="${kr%%:*}"
    if [ "$(cat "$logs/$(jname "$ks" "$ka").run_rc" 2>/dev/null)" = "3" ]; then capped="$(jname "$ks" "$ka")"; fi
  done
  [ -n "$capped" ] && break
  run_job "$s" "$a" "$cap" &
  pids="$pids $!"
done
for p in $pids; do wait "$p"; done   # the jobs only (not caffeinate)
for k in $joblist; do
  ks="${k%%:*}"; kr="${k#*:}"; ka="${kr%%:*}"
  if [ "$(cat "$logs/$(jname "$ks" "$ka").run_rc" 2>/dev/null)" = "3" ]; then capped="$(jname "$ks" "$ka")"; fi
done

fail=0
reason="all runs and passes finished"
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  rc=$(cat "$logs/$n.run_rc" 2>/dev/null || echo 99)
  prc=$(cat "$logs/pass_$n.rc" 2>/dev/null || echo 99)
  if [ "$rc" != "0" ]; then fail=1; fi
  if [ "$prc" != "0" ] && [ -z "$smoke_arg" ]; then fail=1; fi
done
[ -n "$smoke_arg" ] && reason="smoke finished (no pass)"
[ -n "$capped" ] && reason="STOP RULE: a run reached its cap: $capped"
date +%Y-%m-%dT%H:%M:%S%z > "$logs/polish_run_all.end"
if [ "$fail" = "1" ]; then
  if [ "$reason" = "all runs and passes finished" ]; then reason="a run or a pass failed (see the run_rc and pass rc files)"; fi
  end 1 "$reason"
fi
end 0 "$reason"
