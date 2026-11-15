#!/bin/bash
# AR5 step 17, phase A step 4: the recording run J5r on the five disclosed samples (held3, held, mid,
# grad, big): arm J5 of step 16 (the same code, aligned draws) with the hook of fitskip_hook.py
# recording every BFS local optimization, each followed by its boundary pass, by prediction 1
# (fitskip_p1.py: J5r equals the stored J5) and by the final-answer flags (fitskip_final.py), all in
# the same slot. Nothing of held2 or held4 is used.
# DETACHED-SAFE: start it exactly once, from icenine_py/, as
#   cd icenine_py && nohup setsid scripts/autoresearch/ar5/fitskip_record_all.sh > /dev/null 2>&1 &
# (it ignores a hangup, reads no terminal and keeps the machine awake with caffeinate).
#
#   scripts/autoresearch/ar5/fitskip_record_all.sh --dry-run    # guards and commands only
#   ... --cache DIR     another cache folder (tests)      ... --jobs N   at most N heavy (<= 4)
#   ... --smoke-cap M   test form: three short jobs (--smoke 30, a cap of M minutes each, no pass):
#                       J5r on big and on held, and the unhooked arm J5 on big (high500_run.py),
#                       then prediction 1 on big (the two short runs must agree); the cache
#                       defaults to benchmarks/autoresearch/ar5/fitskip/cache_smoke
#   ... --resume        after an interruption: a job whose .done says "run=0 pass=0 p1=0 final=0" with
#                       its outputs is kept; the leftovers of every other job are MOVED to
#                       <cache>/interrupted_<time>/ and the job is run from scratch. An old FAILED
#                       marker (exit not 0) is handled by the script itself: it is moved to
#                       logs/fitskip_record.done.failed.<time> and --resume is switched on. An old
#                       marker "exit=0" is refused (the recording is done).
#
# Order (at most 4 heavy processes at once): held3 (the long pole), held, mid, grad, then big when a
# slot frees. A job is: run (fitskip_run.py --name J5r, cap 300 minutes), then its pass
# (fitskip_pass.py --name J5r --aligned 1), then fitskip_p1.py, then fitskip_final.py. A run at its
# cap writes capped.json and exits 3; its pass is skipped; no new job is then started and the marker
# says so (stop rule). A failed prediction 1 (exit 1) is the stop rule too: the marker says so.
#
# Guards (exit 2, nothing started; a note logs/fitskip_record.refused.<time> holds the reason): the
# real path of an output lies outside the cache (stallhit_cache_guard.py); an output of a job
# exists (not with --resume); fitskip_freeze.py --check fails (frozen.json committed, clean, hashes
# agree); a fitskip run, pass, p1 or final process is already running; an input is missing (samples,
# image stacks, the stored arm-J5 files linked in the cache). held3's inputs and the stored files are
# linked read only by fitskip_link_cache.sh; the others by order_link_inputs.sh.
#
# Files, in <cache>/logs (cache = benchmarks/autoresearch/ar5/fitskip/cache, untracked):
#   J5r_<sample>.log, .run_rc, pass_J5r_<sample>.log / .rc, p1_J5r_<sample>.log / .rc,
#   final_J5r_<sample>.log / .rc, J5r_<sample>.done ("run=<c> pass=<c|skipped> p1=<c|skipped>
#   final=<c|skipped>"), fitskip_record.start / .end, and the MARKER fitskip_record.done: first line
#   "exit=<code>", then "reason=<text>", the end time and one line per job. The script's exit code is
#   the marker's. Rows: <cache>/record/J5r_<sample>.jsonl.gz and final_J5r_<sample>.json (NOBODY reads
#   them before ar5/fitskip/frozen.json is committed; the freeze guard enforces it for held, held3).
#   Wall times are contended and are not timing results.
# Plain bash 3.2 (macOS default).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${FITSKIP_CACHE:-$here/benchmarks/autoresearch/ar5/fitskip/cache}"
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
  cache="$here/benchmarks/autoresearch/ar5/fitskip/cache_smoke"
fi
logs="$cache/logs"
marker="$logs/fitskip_record.done"
export FITSKIP_CACHE="$cache"
export HIGH500_CACHE="$cache"
export POLISH_CACHE="$cache"
export STALLHIT_CACHE="$cache"
export RESCUE_CACHE="$cache"
export ORDER_CACHE="$cache"   # order_pass.py / order_common read this one
export RESEED_CACHE="$cache"  # fit_run.py / fit_pass.py read this one
export BAND_CACHE="$cache"
unset VIRTUAL_ENV

# job list: sample:arm:cap, in start order (arm J5r = the recording run; J5b = the unhooked J5, smoke only)
if [ "$smoke_cap" != "0" ]; then
  joblist="big:J5r:$smoke_cap held:J5r:$smoke_cap big:J5b:$smoke_cap"
  smoke_arg=" --smoke 30"
else
  joblist="held3:J5r:300 held:J5r:300 mid:J5r:300 grad:J5r:300 big:J5r:300"
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
        echo "$n run=$(cat "$logs/$n.run_rc" 2>/dev/null || echo none) pass=$(cat "$logs/pass_$n.rc" 2>/dev/null || echo none) p1=$(cat "$logs/p1_$n.rc" 2>/dev/null || echo none) final=$(cat "$logs/final_$n.rc" 2>/dev/null || echo none)"
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
    echo "exit=2 reason=refused: $1" > "$logs/fitskip_record.refused.$(date +%s)"
  fi
  exit 2
}

# an old marker: exit=0 means the recording is done; anything else is a failed or interrupted
# attempt, which this script handles itself (move it aside, resume)
if [ -e "$marker" ]; then
  first="$(head -1 "$marker" 2>/dev/null)"
  if [ "$first" = "exit=0" ]; then refuse "the marker $marker says exit=0 (the recording is done)"; fi
  if [ "$dry" = "0" ]; then
    mv "$marker" "$marker.failed.$(date +%Y%m%d-%H%M%S)"
    echo "old failed marker ($first) moved aside; --resume switched on"
  else
    echo "dry run: old failed marker ($first) would be moved aside; --resume would be switched on"
  fi
  resume=1
fi

# the files of a job (paths relative to the cache)
jfiles() { # jfiles <sample> <arm> <name>
  if [ "$2" = "J5r" ]; then
    echo "$3 passes_$1/D4_$2 passes_$1/before_$2.npy passes_$1/after_$2.npy record/$3.jsonl record/$3.jsonl.gz record/final_$3.json p1_$3.json logs/$3.log logs/$3.done logs/$3.run_rc logs/pass_$3.log logs/pass_$3.rc logs/p1_$3.log logs/p1_$3.rc logs/final_$3.log logs/final_$3.rc"
  else
    echo "J5_$1_smoke logs/$3.log logs/$3.done logs/$3.run_rc"
  fi
}
job_finished() { # job_finished <sample> <arm> <name>: the outputs say it ended well
  [ "$2" = "J5r" ] || return 1
  [ "$(cat "$logs/$3.done" 2>/dev/null)" = "run=0 pass=0 p1=0 final=0" ] && [ -e "$cache/$3/run.json" ] \
    && [ ! -e "$cache/$3/capped.json" ] && [ -e "$cache/passes_$1/after_$2.npy" ] \
    && [ -e "$cache/record/$3.jsonl.gz" ] && [ -e "$cache/record/final_$3.json" ]
}
kept=""   # names of the jobs kept by --resume
redo=""   # names of the jobs run (again) with leftovers moved

# ---- guard 0: every output must really lie in this cache, not behind a link into another cache.
outs="logs"
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  if [ "$a" = "J5r" ]; then outs="$outs passes_$s"; fi
  outs="$outs $(jfiles "$s" "$a" "$n")"
done
if ! uv run python scripts/autoresearch/ar5/stallhit_cache_guard.py --cache "$cache" $outs; then
  echo "refusing to start: an output's real path lies outside the cache $cache" >&2
  exit 2
fi

# ---- guard 1: no output exists (--resume: finished jobs are kept, leftovers are moved later)
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

# ---- guard 2: the freeze
ferr="$(mktemp)"
if ! uv run python scripts/autoresearch/ar5/fitskip_freeze.py --check > /dev/null 2>"$ferr"; then
  msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
  refuse "${msg:-fitskip freeze check failed}"
fi
rm -f "$ferr"

# ---- guard 3: nothing of this step is running
if pgrep -f "autoresearch/ar5/fitskip_(run|pass|p1|final)\.py" > /dev/null 2>&1; then
  refuse "a fitskip run, pass, p1 or final process is already running"
fi

# ---- guard 4: inputs (link the untracked ones first; links only)
lk=""
for j in $joblist; do lk="$lk ${j%%:*}"; done
lk="$(echo $lk | tr ' ' '\n' | sort -u | tr '\n' ' ')"
if [ "$dry" = "1" ]; then
  scripts/autoresearch/ar5/order_link_inputs.sh --dry-run > /dev/null 2>&1 || refuse "a source of the input links is missing (order_link_inputs.sh --dry-run)"
  scripts/autoresearch/ar5/fitskip_link_cache.sh --dry-run --cache "$cache" $lk > /dev/null 2>&1 || refuse "a source of the cache links is missing (fitskip_link_cache.sh --dry-run)"
else
  scripts/autoresearch/ar5/order_link_inputs.sh > /dev/null 2>&1 || refuse "a source of the input links is missing (order_link_inputs.sh)"
  scripts/autoresearch/ar5/fitskip_link_cache.sh --cache "$cache" $lk > /dev/null 2>&1 || refuse "a source of the cache links is missing (fitskip_link_cache.sh)"
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
    held3) need="$need $here/benchmarks/autoresearch/ar3/held3/gradient/sample/gradient_grid.mic $st/held3_clean.npy" ;;
  esac
done
if [ "$smoke_cap" = "0" ]; then   # the linked stored arm-J5 files (existence only)
  for s in held3 held mid grad big; do
    need="$need $cache/J5_$s/run.json $cache/J5_$s/recon.mic.gz $cache/passes_$s/after_J5.npy $cache/passes_$s/D4_J5/pass.json"
  done
fi
for p in $need; do
  if [ ! -e "$p" ]; then refuse "input missing: $p"; fi
done

run_cmd() { # run_cmd <sample> <arm> <cap>: the run command (text, for the dry run and the log)
  case "$2" in
    J5r) echo "RESEED_SAMPLE=$1 FITSKIP_CACHE=$cache uv run python scripts/autoresearch/ar5/fitskip_run.py --name J5r --cap-min $3$smoke_arg" ;;
    J5b) echo "RESEED_SAMPLE=$1 HIGH500_CACHE=$cache uv run python scripts/autoresearch/ar5/high500_run.py --cap-min $3$smoke_arg" ;;
  esac
}
pass_cmd() { echo "RESEED_SAMPLE=$1 FITSKIP_CACHE=$cache uv run python scripts/autoresearch/ar5/fitskip_pass.py --name J5r --aligned 1"; }
p1_cmd() { echo "RESEED_SAMPLE=$1 FITSKIP_CACHE=$cache uv run python scripts/autoresearch/ar5/fitskip_p1.py --name J5r"; }
final_cmd() { echo "RESEED_SAMPLE=$1 FITSKIP_CACHE=$cache uv run python scripts/autoresearch/ar5/fitskip_final.py --name J5r"; }
smoke_p1_cmd() { echo "RESEED_SAMPLE=big FITSKIP_CACHE=$cache uv run python scripts/autoresearch/ar5/fitskip_p1.py --name J5r --stored J5 --suffix _smoke"; }

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
    if [ "$smoke_cap" = "0" ]; then
      echo "    then, if it finished: $(pass_cmd "$s") > $logs/pass_$n.log 2>&1"
      echo "    then: $(p1_cmd "$s") > $logs/p1_$n.log 2>&1"
      echo "    then: $(final_cmd "$s") > $logs/final_$n.log 2>&1"
    fi
  done
  [ "$smoke_cap" != "0" ] && echo "  then: $(smoke_p1_cmd)"
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
  echo "resume $(date +%Y-%m-%dT%H:%M:%S%z) kept:${kept:- none} rerun:${redo:- none}" >> "$logs/fitskip_record.resume"
else
  date +%Y-%m-%dT%H:%M:%S%z > "$logs/fitskip_record.start"
fi

run_job() { # run_job <sample> <arm> <cap>   (a subshell: run, then pass, p1, final)
  s="$1"; a="$2"; cap="$3"; n="$(jname "$s" "$a")"
  eval "$(run_cmd "$s" "$a" "$cap")" > "$logs/$n.log" 2>&1
  rc=$?
  echo "$rc" > "$logs/$n.run_rc"
  prc=skipped; p1rc=skipped; frc=skipped
  d="$cache/$n${smoke_arg:+_smoke}"
  if [ "$a" = "J5r" ] && [ "$rc" = "0" ] && [ -e "$d/run.json" ] && [ ! -e "$d/capped.json" ] && [ -z "$smoke_arg" ]; then
    eval "$(pass_cmd "$s")" > "$logs/pass_$n.log" 2>&1
    prc=$?
    echo "$prc" > "$logs/pass_$n.rc"
    if [ "$prc" = "0" ]; then
      eval "$(p1_cmd "$s")" > "$logs/p1_$n.log" 2>&1
      p1rc=$?
      echo "$p1rc" > "$logs/p1_$n.rc"
      eval "$(final_cmd "$s")" > "$logs/final_$n.log" 2>&1
      frc=$?
      echo "$frc" > "$logs/final_$n.rc"
    fi
  fi
  echo "run=$rc pass=$prc p1=$p1rc final=$frc" > "$logs/$n.done"
}

any_stop() { # prints "<name> capped" or "<name> p1" for the first job that hit the stop rule
  for k in $joblist; do
    ks="${k%%:*}"; kr="${k#*:}"; ka="${kr%%:*}"; kn="$(jname "$ks" "$ka")"
    if [ "$(cat "$logs/$kn.run_rc" 2>/dev/null)" = "3" ]; then echo "$kn capped"; return; fi
    if [ "$(cat "$logs/p1_$kn.rc" 2>/dev/null || echo 0)" = "1" ]; then echo "$kn p1"; return; fi
  done
}

stopped=""
pids=""
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; cap="${rest#*:}"; n="$(jname "$s" "$a")"
  case " $kept " in *" $n "*) continue ;; esac
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
  # stop rule: after a run hit its cap or prediction 1 failed, no new job is started
  stopped="$(any_stop)"
  [ -n "$stopped" ] && break
  run_job "$s" "$a" "$cap" &
  pids="$pids $!"
done
for p in $pids; do wait "$p"; done   # the jobs only (not caffeinate)
[ -z "$stopped" ] && stopped="$(any_stop)"

# the smoke form: the two short runs of big must agree (prediction 1 on the short run)
smoke_rc=0
if [ -n "$smoke_arg" ]; then
  eval "$(smoke_p1_cmd)" > "$logs/p1_smoke_big.log" 2>&1
  smoke_rc=$?
  echo "$smoke_rc" > "$logs/p1_smoke_big.rc"
fi

fail=0
reason="all runs, passes, prediction-1 checks and final flags finished"
for j in $joblist; do
  s="${j%%:*}"; rest="${j#*:}"; a="${rest%%:*}"; n="$(jname "$s" "$a")"
  rc=$(cat "$logs/$n.run_rc" 2>/dev/null || echo 99)
  if [ "$rc" != "0" ]; then fail=1; fi
  if [ "$a" = "J5r" ] && [ -z "$smoke_arg" ]; then
    for k in pass p1 final; do
      v=$(cat "$logs/${k}_$n.rc" 2>/dev/null || echo 99)
      if [ "$v" != "0" ]; then fail=1; fi
    done
  fi
done
[ -n "$smoke_arg" ] && reason="smoke finished (no pass); prediction 1 on the short runs: rc=$smoke_rc"
[ "$smoke_rc" != "0" ] && fail=1
[ -n "$stopped" ] && reason="STOP RULE: $stopped"
date +%Y-%m-%dT%H:%M:%S%z > "$logs/fitskip_record.end"
if [ "$fail" = "1" ]; then
  if [ "$reason" = "all runs, passes, prediction-1 checks and final flags finished" ]; then reason="a run, pass, p1 or final step failed (see the rc files)"; fi
  end 1 "$reason"
fi
end 0 "$reason"
