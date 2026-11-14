#!/bin/bash
# AR5 step 12, phase A step 5 and 6: build held2, then the whole-layer runs of arms S8, G and H on
# held2, H on held, and H on mid, grad, big, each followed by its boundary pass.
# A FOREGROUND script: the main session starts it exactly once, from icenine_py/ (it does not
# detach; the caller keeps it alive, e.g. with nohup setsid).
#
#   cd icenine_py && scripts/autoresearch/ar5/stallhit_run_all.sh             # the real thing
#   cd icenine_py && scripts/autoresearch/ar5/stallhit_run_all.sh --dry-run   # guards and commands only
#   ... --cache DIR     another cache folder (tests)      ... --jobs N   at most N heavy (<= 4)
#   ... --resume        after an interruption (the script was killed): the held2 build skips its
#                       finished stages; a job whose .done says "run=0 pass=0" with its run.json and
#                       pass outputs is kept; the leftovers of every other job are MOVED to
#                       <cache>/interrupted_<time>/ and the job is run from scratch;
#                       stallhit_run_all.start is kept, stallhit_run_all.resume gets a line; the
#                       marker covers all jobs. With --dry-run it lists kept / moved, changes nothing.
#                       Start it detached: nohup setsid scripts/autoresearch/ar5/stallhit_run_all.sh --resume
#
# Order (at most 4 heavy processes at once; a job is a run, then its pass in the same slot):
#   1. held2 build (stallhit_make_held2.py: sample, render with 1 worker, stack), alone
#   2. S8 held2, G held2, H held2, H held, side by side (4 slots)
#   3. H mid, H grad, H big, each as a slot frees
# Arms (each followed by its boundary pass):
#   S8  fit_run.py --name S8 (step 8's S8: stock seed step), then fit_pass.py --name S8
#   G   rescue_run.py (arm G of step 9: interleaved order, aligned draws, rule G1), then
#       order_pass.py --name G --aligned 1
#   H   stallhit_run.py (arm G with rule PH(3000)), then order_pass.py --name H --aligned 1
# Caps (--cap-min): S8 on held2 480 (held's S8 needed over 120), every other job 300. A run at its
# cap writes capped.json and exits 3; its pass is skipped; after a held2 run hits its cap no new job
# is started and the marker says so (stop rule of the plan).
#
# Guard 0 (exit 2, nothing written, reason on stderr): the real path of an output (logs, run folder,
# pass folder or file) lies outside the cache (stallhit_cache_guard.py; review fix 1, second review).
# Guards (exit 2, nothing started; a note logs/stallhit_run_all.refused.<time> holds the reason): an
# output of a job (run folder, pass files, log, marker) or the final marker exists; the freeze
# check fails (stallhit_freeze.py --check: stallhit/frozen.json is committed, clean, and the hashes
# of the scripts and of the library agree); a stallhit, rescue, order, fit or patience process
# (run, pass, replay) is running; an input is missing (samples and image stacks of held, mid, grad,
# big; remedies_pick.json; the stored arm A and G runs and passes, linked in the cache). The
# held2 inputs are made by step 1 and are not required before it.
# The untracked inputs are linked first by order_link_inputs.sh and stallhit_link_cache.sh (links
# only, never overwrites).
#
# Files, in <cache>/logs (cache = benchmarks/autoresearch/ar5/stallhit/cache, untracked):
#   <arm>_<sample>.log, .run_rc (exit code of the run), pass_<arm>_<sample>.log, pass_<..>.rc
#   ("skipped" or the code), <arm>_<sample>.done ("run=<code> pass=<code|skipped>"),
#   held2_<stage>.log / .rc (the build), stallhit_run_all.start / .end, and the MARKER
#   stallhit_run_all.done: first line "exit=<code>", then "reason=<text>", the end time and one
#   line per job. The script's exit code is the marker's. Run folders <cache>/<arm>_<sample>/,
#   passes in <cache>/passes_<sample>/. Wall times are contended and are not timing results.
# Plain bash 3.2 (macOS default).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${STALLHIT_CACHE:-$here/benchmarks/autoresearch/ar5/stallhit/cache}"
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
marker="$logs/stallhit_run_all.done"
export STALLHIT_CACHE="$cache"
export RESCUE_CACHE="$cache"
export ORDER_CACHE="$cache"   # order_pass.py / order_common read this one
export RESEED_CACHE="$cache"  # fit_run.py / fit_pass.py read this one
unset VIRTUAL_ENV

# job list: sample:arm:cap, in start order (the build is step 1, before them)
if [ "$smoke_cap" != "0" ]; then
  joblist="big:H:$smoke_cap"
else
  joblist="held2:S8:480 held2:G:300 held2:H:300 held:H:300 mid:H:300 grad:H:300 big:H:300"
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
      echo "held2_build build=$(cat "$logs/held2_build.rc" 2>/dev/null || echo none) render=$(cat "$logs/held2_render.rc" 2>/dev/null || echo none) stack=$(cat "$logs/held2_stack.rc" 2>/dev/null || echo none)"
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
    echo "exit=2 reason=refused: $1" > "$logs/stallhit_run_all.refused.$(date +%s)"
  fi
  exit 2
}

# the files of a job (paths relative to the cache)
jfiles() { # jfiles <sample> <arm> <name>
  echo "$3 passes_$1/D4_$2 passes_$1/before_$2.npy passes_$1/after_$2.npy logs/$3.log logs/$3.done logs/$3.run_rc logs/pass_$3.log logs/pass_$3.rc"
}
# a finished job (resume): .done says "run=0 pass=0" and the run and pass outputs exist
job_finished() { # job_finished <sample> <arm> <name>
  [ "$(cat "$logs/$3.done" 2>/dev/null)" = "run=0 pass=0" ] && [ -e "$cache/$3/run.json" ] \
    && [ ! -e "$cache/$3/capped.json" ] && [ -e "$cache/passes_$1/after_$2.npy" ] \
    && [ -e "$cache/passes_$1/D4_$2/pass.json" ] && [ -e "$cache/passes_$1/before_$2.npy" ]
}
kept=""   # names of the jobs kept by --resume
redo=""   # names of the jobs run (again) with leftovers moved

# ---- guard 0 (review fix 1 of the second review): every output (logs, run folders, pass folders
# and files) must really lie in this cache, not behind a link into another worktree's cache. Runs
# first, also in a dry run; on a refusal it writes no note (the logs folder may be the link).
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

# ---- guard 2: the freeze (not in test mode)
if [ "$smoke_cap" = "0" ]; then
  ferr="$(mktemp)"
  if ! uv run python scripts/autoresearch/ar5/stallhit_freeze.py --check > /dev/null 2>"$ferr"; then
    msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
    refuse "${msg:-stallhit freeze check failed}"
  fi
  rm -f "$ferr"
fi

# ---- guard 3: nothing of this step (or an arm it builds on) is running
if pgrep -f "autoresearch/ar5/(order|rescue|patience|stallhit|fit)_(run|pass|replay|make_held2)\.py" > /dev/null 2>&1; then
  refuse "a stallhit / patience / rescue / order / fit run, replay, pass or build process is already running"
fi

# ---- guard 4: inputs (link the untracked ones first; links only)
if [ "$dry" = "1" ]; then
  scripts/autoresearch/ar5/order_link_inputs.sh --dry-run > /dev/null 2>&1 || refuse "a source of the input links is missing (order_link_inputs.sh --dry-run)"
  scripts/autoresearch/ar5/stallhit_link_cache.sh --dry-run --cache "$cache" > /dev/null 2>&1 || refuse "a source of the cache links is missing (stallhit_link_cache.sh --dry-run)"
else
  scripts/autoresearch/ar5/order_link_inputs.sh > /dev/null 2>&1 || refuse "a source of the input links is missing (order_link_inputs.sh)"
  scripts/autoresearch/ar5/stallhit_link_cache.sh --cache "$cache" > /dev/null 2>&1 || refuse "a source of the cache links is missing (stallhit_link_cache.sh)"
fi
wt="$(cd "$here/../.." && pwd)"
st="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
need=""
for j in $joblist; do
  s="${j%%:*}"
  case "$s" in
    big) need="$need $here/benchmarks/autoresearch/ar3/big/gradient/sample/gradient_grid.mic $st/big_clean.npy" ;;
    mid) need="$need $here/benchmarks/autoresearch/ar3/mid/gradient/sample/gradient_grid.mic $st/mid_clean.npy" ;;
    grad) need="$need $here/benchmarks/autoresearch/ar3/gradient/sample/gradient_grid.mic $st/grad_clean.npy" ;;
    held) need="$need $here/benchmarks/autoresearch/ar3/held/gradient/sample/gradient_grid.mic $st/held_clean.npy" ;;
    held2) ;;   # made by step 1
  esac
done
need="$need $here/benchmarks/autoresearch/ar3/reseed/remedies_pick.json"
mixc="$wt/agent-a5a3601b746680640/icenine_py/benchmarks/autoresearch/ar5/mix/cache"
for s in big grad mid; do
  need="$need $mixc/S11_$s/run.json $mixc/passes_$s/after_S11.npy"
done
for s in big grad mid held; do   # the linked stored arm A and G runs (existence only)
  need="$need $cache/A_$s/run.json $cache/passes_$s/after_A.npy $cache/passes_$s/D4_A/pass.json"
  need="$need $cache/G_$s/run.json $cache/passes_$s/after_G.npy $cache/passes_$s/D4_G/pass.json"
done
need="$need $wt/agent-a635b2f904279bdc9/icenine_py/benchmarks/autoresearch/ar3/reseed/cache/S8_big/e1_log.json.gz"
need="$need $wt/agent-a635b2f904279bdc9/icenine_py/benchmarks/autoresearch/ar3/reseed/cache/S8_grad/e1_log.json.gz"
need="$need $wt/agent-a9c88b610fb10e1b6/icenine_py/benchmarks/autoresearch/ar3/start/cache/S8_mid/e1_log.json.gz"
if [ "$smoke_cap" = "0" ]; then
  hc="$wt/agent-a8a1ac8e80f8913cf/icenine_py/benchmarks/autoresearch/ar5/held/cache"
  need="$need $hc/S11_held/run.json $hc/passes_held/after_S11.npy"
  need="$need $wt/agent-af8961874e131bda2/icenine_py/benchmarks/autoresearch/ar3/held/cache/S8_held/e1_log.json.gz"
fi
for p in $need; do
  if [ ! -e "$p" ]; then refuse "input missing: $p"; fi
done

run_cmd() { # run_cmd <sample> <arm> <cap>: the run command (text, for the dry run and the log)
  case "$2" in
    S8) echo "RESEED_SAMPLE=$1 RESEED_CACHE=$cache uv run python scripts/autoresearch/ar5/fit_run.py --name S8 --cap-min $3" ;;
    G) echo "RESEED_SAMPLE=$1 RESCUE_CACHE=$cache uv run python scripts/autoresearch/ar5/rescue_run.py --cap-min $3" ;;
    H) echo "RESEED_SAMPLE=$1 STALLHIT_CACHE=$cache uv run python scripts/autoresearch/ar5/stallhit_run.py --cap-min $3" ;;
  esac
}
pass_cmd() { # pass_cmd <sample> <arm>
  case "$2" in
    S8) echo "RESEED_SAMPLE=$1 RESEED_CACHE=$cache uv run python scripts/autoresearch/ar5/fit_pass.py --name S8" ;;
    *) echo "RESEED_SAMPLE=$1 ORDER_CACHE=$cache uv run python scripts/autoresearch/ar5/order_pass.py --name $2 --aligned 1" ;;
  esac
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
  if [ "$smoke_cap" = "0" ]; then
    resarg=""; [ "$resume" = "1" ] && resarg=" --resume"
    echo "  held2 build first: uv run python scripts/autoresearch/ar5/stallhit_make_held2.py$resarg > $logs/held2_build.log 2>&1"
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
  echo "resume $(date +%Y-%m-%dT%H:%M:%S%z) kept:${kept:- none} rerun:${redo:- none}" >> "$logs/stallhit_run_all.resume"
else
  date +%Y-%m-%dT%H:%M:%S%z > "$logs/stallhit_run_all.start"
fi

# ---- step 1: the held2 build (alone), unless this is a smoke test
if [ "$smoke_cap" = "0" ]; then
  barg=""; [ "$resume" = "1" ] && barg="--resume"
  uv run python scripts/autoresearch/ar5/stallhit_make_held2.py $barg > "$logs/held2_build.log" 2>&1
  brc=$?
  if [ "$brc" != "0" ]; then
    date +%Y-%m-%dT%H:%M:%S%z > "$logs/stallhit_run_all.end"
    end 1 "the held2 build failed (exit $brc; see held2_build.log)"
  fi
  # the build made the held2 inputs: they must exist now
  for p in "$here/benchmarks/autoresearch/ar3/held2/gradient/sample/gradient_grid.mic" "$st/held2_clean.npy"; do
    if [ ! -e "$p" ]; then
      date +%Y-%m-%dT%H:%M:%S%z > "$logs/stallhit_run_all.end"
      end 1 "held2 input missing after the build: $p"
    fi
  done
fi

run_job() { # run_job <sample> <arm> <cap>   (a subshell: run, then pass)
  s="$1"; a="$2"; cap="$3"; n="$(jname "$s" "$a")"
  eval "$(run_cmd "$s" "$a" "$cap")" > "$logs/$n.log" 2>&1
  rc=$?
  echo "$rc" > "$logs/$n.run_rc"
  d="$cache/$n"
  if [ "$rc" = "0" ] && [ -e "$d/run.json" ] && [ ! -e "$d/capped.json" ]; then
    eval "$(pass_cmd "$s" "$a")" > "$logs/pass_$n.log" 2>&1
    prc=$?
    echo "$prc" > "$logs/pass_$n.rc"
  else
    prc=skipped
    echo skipped > "$logs/pass_$n.rc"
  fi
  echo "run=$rc pass=$prc" > "$logs/$n.done"
}

capped_held2=""
pids=""
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; cap="${rest#*:}"
  case " $kept " in *" $(jname "$s" "$a") "*) continue ;; esac
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
  # stop rule: after a held2 run hit its cap, no new job is started
  for h in S8_held2 G_held2 H_held2; do
    if [ "$(cat "$logs/$h.run_rc" 2>/dev/null)" = "3" ]; then capped_held2="$h"; fi
  done
  [ -n "$capped_held2" ] && break
  run_job "$s" "$a" "$cap" &
  pids="$pids $!"
done
for p in $pids; do wait "$p"; done   # the jobs only (not caffeinate)
for h in S8_held2 G_held2 H_held2; do
  if [ "$(cat "$logs/$h.run_rc" 2>/dev/null)" = "3" ]; then capped_held2="$h"; fi
done

fail=0
reason="all runs and passes finished"
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  rc=$(cat "$logs/$n.run_rc" 2>/dev/null || echo 99)
  prc=$(cat "$logs/pass_$n.rc" 2>/dev/null || echo 99)
  if [ "$rc" = "3" ]; then reason="a run reached its cap (unfinished): $n"; fi
  if [ "$rc" != "0" ] || [ "$prc" != "0" ]; then fail=1; fi
done
[ -n "$capped_held2" ] && reason="STOP RULE: a held2 run reached its cap: $capped_held2"
date +%Y-%m-%dT%H:%M:%S%z > "$logs/stallhit_run_all.end"
if [ "$fail" = "1" ]; then
  if [ "$reason" = "all runs and passes finished" ]; then reason="a run or a pass failed (see the run_rc and pass rc files)"; fi
  end 1 "$reason"
fi
end 0 "$reason"
