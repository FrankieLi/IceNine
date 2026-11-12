#!/bin/bash
# AR5 step 8, phase A step 5: the nine whole-layer runs and their boundary passes.
# A FOREGROUND script: the main session starts it exactly once, from icenine_py/ (it does not
# detach; the caller keeps it alive, e.g. with nohup or a background task).
#
#   cd icenine_py && scripts/autoresearch/ar5/order_run_all.sh             # the real thing
#   cd icenine_py && scripts/autoresearch/ar5/order_run_all.sh --dry-run   # guards and commands only
#   ... --cache DIR     another cache folder (tests)      ... --jobs N   at most N heavy (<= 4)
#   ... --smoke-cap M   test mode: only the three runs on big, each capped at M minutes (the cap
#                       path: they exit 3, passes are skipped); no freeze guard; held untouched
#
# Jobs, started in this order, at most 4 heavy processes at once (a job is a run, then its pass
# in the same slot):  held A, held Y, mid A, mid Y, grad A, grad Y, big A, big Y, big Yl
#   arms: A interleaved order, Y sweep order (both aligned draws), Yl sweep order with the
#   library's own draws (its pass too). Caps (--cap-min): held 300, the rest 240. A run at its cap
#   writes capped.json and exits 3; its pass is skipped.
#
# Guards (exit 2, nothing started; a note logs/order_run_all.refused.<time> holds the reason): an
# output of a job (run folder, pass files, log, marker) or the final marker exists; the freeze
# check fails (order_freeze.py --check: frozen.json committed, clean, hashes of the scripts and of
# the library agree); another order_run.py process is running; an input is missing (samples,
# image stacks, remedies_pick.json, the stored S8 / S11 runs and passes for the comparison).
# The untracked inputs are linked first by order_link_inputs.sh (links only, never overwrites).
#
# Files, in <cache>/logs (cache = benchmarks/autoresearch/ar5/order/cache, untracked):
#   <arm>_<sample>.log, .run_rc (exit code of the run), pass_<arm>_<sample>.log, pass_<..>.rc
#   ("skipped" or the code), <arm>_<sample>.done ("run=<code> pass=<code|skipped>"),
#   order_run_all.start / .end, and the MARKER order_run_all.done: first line "exit=<code>", then
#   "reason=<text>", the end time and one line per job. The script's exit code is the marker's.
#   Run folders <cache>/<arm>_<sample>/, passes in <cache>/passes_<sample>/.
# Wall times are contended and are not timing results. Plain bash 3.2 (macOS default).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${ORDER_CACHE:-$here/benchmarks/autoresearch/ar5/order/cache}"
dry=0
jobs=4
smoke_cap=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    --jobs) shift; jobs="$1" ;;
    --smoke-cap) shift; smoke_cap="$1" ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
if [ "$jobs" -gt 4 ] || [ "$jobs" -lt 1 ]; then echo "--jobs must be 1 to 4" >&2; exit 2; fi
logs="$cache/logs"
marker="$logs/order_run_all.done"
export ORDER_CACHE="$cache"
unset VIRTUAL_ENV

if [ "$smoke_cap" != "0" ]; then
  joblist="big:A:$smoke_cap big:Y:$smoke_cap big:Yl:$smoke_cap"
else
  joblist="held:A:300 held:Y:300 mid:A:240 mid:Y:240 grad:A:240 grad:Y:240 big:A:240 big:Y:240 big:Yl:240"
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
    echo "exit=2 reason=refused: $1" > "$logs/order_run_all.refused.$(date +%s)"
  fi
  exit 2
}

# ---- guard 1: no output exists
if [ -e "$marker" ]; then refuse "the marker $marker already exists"; fi
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  for p in "$cache/$n" "$cache/passes_$s/D4_$a" "$cache/passes_$s/before_$a.npy" \
           "$cache/passes_$s/after_$a.npy" "$logs/$n.log" "$logs/$n.done" "$logs/$n.run_rc" \
           "$logs/pass_$n.log" "$logs/pass_$n.rc"; do
    if [ -e "$p" ]; then refuse "$p already exists"; fi
  done
done

# ---- guard 2: the freeze (not in test mode)
if [ "$smoke_cap" = "0" ]; then
  ferr="$(mktemp)"
  if ! uv run python scripts/autoresearch/ar5/order_freeze.py --check > /dev/null 2>"$ferr"; then
    msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
    refuse "${msg:-freeze check failed}"
  fi
  rm -f "$ferr"
fi

# ---- guard 3: nothing of this step is running
if pgrep -f "autoresearch/ar5/order_(run|pass)\.py" > /dev/null 2>&1; then
  refuse "an order_run.py / order_pass.py process is already running"
fi

# ---- guard 4: inputs (link the untracked ones first; links only)
if [ "$dry" = "1" ]; then
  scripts/autoresearch/ar5/order_link_inputs.sh --dry-run > /dev/null 2>&1 || refuse "a source of the input links is missing (order_link_inputs.sh --dry-run)"
else
  scripts/autoresearch/ar5/order_link_inputs.sh > /dev/null 2>&1 || refuse "a source of the input links is missing (order_link_inputs.sh)"
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
  esac
done
need="$need $here/benchmarks/autoresearch/ar3/reseed/remedies_pick.json"
mixc="$wt/agent-a5a3601b746680640/icenine_py/benchmarks/autoresearch/ar5/mix/cache"
for s in big grad mid; do
  need="$need $mixc/S11_$s/run.json $mixc/passes_$s/after_S11.npy"
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

if [ "$dry" = "1" ]; then
  echo "dry run: guards passed (cache $cache, at most $jobs heavy processes)"
  for j in $joblist; do
    s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; cap="${rest#*:}"; n="$(jname "$s" "$a")"
    al=1; [ "$a" = "Yl" ] && al=0
    echo "  RESEED_SAMPLE=$s ORDER_CACHE=$cache uv run python scripts/autoresearch/ar5/order_run.py --arm $a --cap-min $cap > $logs/$n.log 2>&1"
    echo "    then, if it finished: ... order_pass.py --name $a --aligned $al > $logs/pass_$n.log 2>&1"
  done
  echo "  then: marker $marker (exit=<code>, reason=<text>)"
  exit 0
fi

mkdir -p "$logs"
caffeinate -i -w $$ &
date +%Y-%m-%dT%H:%M:%S%z > "$logs/order_run_all.start"

run_job() { # run_job <sample> <arm> <cap>   (a subshell: run, then pass)
  s="$1"; a="$2"; cap="$3"; n="$(jname "$s" "$a")"
  al=1; [ "$a" = "Yl" ] && al=0
  RESEED_SAMPLE="$s" uv run python scripts/autoresearch/ar5/order_run.py --arm "$a" --cap-min "$cap" \
    > "$logs/$n.log" 2>&1
  rc=$?
  echo "$rc" > "$logs/$n.run_rc"
  d="$cache/$n"
  if [ "$rc" = "0" ] && [ -e "$d/run.json" ] && [ ! -e "$d/capped.json" ]; then
    RESEED_SAMPLE="$s" uv run python scripts/autoresearch/ar5/order_pass.py --name "$a" --aligned "$al" \
      > "$logs/pass_$n.log" 2>&1
    prc=$?
    echo "$prc" > "$logs/pass_$n.rc"
  else
    prc=skipped
    echo skipped > "$logs/pass_$n.rc"
  fi
  echo "run=$rc pass=$prc" > "$logs/$n.done"
}

pids=""
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; cap="${rest#*:}"
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
  run_job "$s" "$a" "$cap" &
  pids="$pids $!"
done
for p in $pids; do wait "$p"; done   # the jobs only (not caffeinate)

fail=0
reason="all runs and passes finished"
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  rc=$(cat "$logs/$n.run_rc" 2>/dev/null || echo 99)
  prc=$(cat "$logs/pass_$n.rc" 2>/dev/null || echo 99)
  if [ "$rc" = "3" ]; then reason="a run reached its cap (unfinished): $n"; fi
  if [ "$rc" != "0" ] || [ "$prc" != "0" ]; then fail=1; fi
done
date +%Y-%m-%dT%H:%M:%S%z > "$logs/order_run_all.end"
if [ "$fail" = "1" ]; then
  if [ "$reason" = "all runs and passes finished" ]; then reason="a run or a pass failed (see the run_rc and pass rc files)"; fi
  end 1 "$reason"
fi
end 0 "$reason"
