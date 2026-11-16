#!/bin/bash
# AR5 step 21, phase A step 5: the held5 build and pilot, arms S8 and J5 on held5 with their
# boundary passes, then arms J5C and J5S on J5's stored map. A FOREGROUND script: the main session
# starts it from icenine_py/ detached:
#
#   cd icenine_py && nohup caffeinate -i scripts/autoresearch/ar5/screen_run_all.sh --resume \
#       > /dev/null 2>&1 &
#   ... --dry-run      guards and commands only (nothing written)
#   ... --cache DIR    another cache folder (tests)
#   ... --smoke        test form (needs --sample held4 or another disclosed sample whose inputs are
#                      linked in the cache): short jobs, no pass: S8 and J5 with --smoke 30 and a cap
#                      of --smoke-cap M minutes (default 3), then J5C --limit 5 and J5S --limit 20
#                      on the linked stored maps. The cache defaults to <screen>/cache_smoke. Exit 3
#                      (the cap) is the S8 / J5 smoke jobs' normal end. No build, no pilot.
#   ... --sample S     the sample of a smoke (default held5; any other only with --smoke)
#   ... --resume       after an interruption or a FAILED earlier attempt: the held5 build skips its
#                      finished stages and moves an unfinished stage's output aside; a job whose
#                      .done says it ended 0 with its outputs is kept; the leftovers of every other
#                      job and an old marker (failed or not) are MOVED to <cache>/interrupted_<time>/
#                      and those jobs run from scratch. Safe on a first start too.
#
# Order (at most 3 heavy processes at once):
#   1. the held5 build (screen_build_held5.py: sample, render with 1 worker, stack; about 15
#      minutes; one process);
#   2. a pilot: S8 and J5 on held5 with --smoke 30 and a cap of 5 minutes each, in <cache>_pilot, two
#      processes. Only their exit codes (0 or 3 = the cap, normal for a smoke) and the existence of
#      their folders are looked at, never their content; the elapsed time goes to logs/pilot.time;
#   3. side by side: S8 (fit_run.py --name S8, then fit_pass.py --name S8: the long pole) and J5
#      (high500_run.py, then high500_pass.py --name J5 --aligned 1, then, together, arm J5C
#      (screen_run.py --arm J5C --source stored) and arm J5S (screen_run.py --arm J5S: the screen
#      pass and the repair on the wide flag)): 3 processes at the most.
# Caps (--cap-min): S8 480, J5 300, J5C 300, J5S 300. A run at its cap writes capped.json and exits
# 3 (J5C / J5S: screen_run.py --cap-min); a run's pass is skipped then; the marker says so (stop
# rule). Stats and the summary are left for later.
#
# Guard 0 (exit 2, nothing written): the real path of an output lies outside the cache
# (stallhit_cache_guard.py). Guards (exit 2, nothing started; a note logs/screen_run_all.refused.
# <time> holds the reason): an output of a job or the final marker exists (not with --resume);
# ar5/screen/gate.json is missing (the gate passed on the disclosed samples); screen_freeze.py
# --check fails (frozen.json committed, clean, hashes agree); a screen / cmarepair / high500 / fit
# run, pass, launch or build process is running; an input is missing. Nothing of held2 is used
# (screen_common refuses RESEED_SAMPLE=held2).
#
# Files, in <cache>/logs (cache = benchmarks/autoresearch/ar5/screen/cache, untracked):
#   screen_run_all.out (this script's own output), <arm>_held5.log, .run_rc, pass_<arm>_held5.log,
#   pass_<arm>_held5.rc, <arm>_held5.done ("run=<code> pass=<code|skipped>"), J5C_held5.log / .rc /
#   .done and J5S_held5.log / .rc / .done ("exit=<code>"), held5_driver.log / .rc and
#   held5_<stage>.log / .rc (the build), pilot_<arm>.rc, pilot.time, screen_run_all.start / .end,
#   and the MARKER screen_runs.done: first line "exit=<code>", then "reason=<text>", the end time and
#   one line per job. The script's exit code is the marker's. Run folders <cache>/<arm>_held5/,
#   passes in <cache>/passes_held5/, the arms in <cache>/J5C_held5/ and <cache>/J5S_held5/. The
#   readers of this step's layer (screen_guard) open nothing of held5 before the marker exists;
#   the J5C and J5S jobs alone (SCREEN_PRODUCER=1, set here for them) open their own inputs by
#   name. NOBODY reads held5's logs or maps before the marker says all runs and repairs ended (plan
#   stop rule). Wall times are contended and are not timing results.
# Plain bash 3.2 (macOS default).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${SCREEN_CACHE:-$here/benchmarks/autoresearch/ar5/screen/cache}"
dry=0
resume=0
smoke=0
smoke_cap=3
cache_given=0
smp=held5
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1"; cache_given=1 ;;
    --smoke) smoke=1 ;;
    --smoke-cap) shift; smoke_cap="$1" ;;
    --sample) shift; smp="$1" ;;
    --resume) resume=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
case "$smp" in held2*|*held2*) echo "held2 is not used by this step" >&2; exit 2 ;; esac
if [ "$smp" != "held5" ] && [ "$smoke" = "0" ]; then echo "--sample other than held5 only with --smoke" >&2; exit 2; fi
case "$smp" in held5|held4|held3|held|mid|grad|big) ;; *) echo "unknown sample: $smp" >&2; exit 2 ;; esac
if [ "$smoke" = "1" ] && [ "$cache_given" = "0" ]; then
  cache="$here/benchmarks/autoresearch/ar5/screen/cache_smoke"
fi
logs="$cache/logs"
marker="$logs/screen_runs.done"
pilot_cache="${cache}_pilot"
export SCREEN_CACHE="$cache"
unset VIRTUAL_ENV

gate="$here/benchmarks/autoresearch/ar5/screen/gate.json"
st="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
grid="$here/benchmarks/autoresearch/ar3/$smp/gradient/sample/gradient_grid.mic"
if [ "$smoke" = "1" ]; then
  joblist="S8:$smoke_cap J5:$smoke_cap"
  smoke_arg=" --smoke 30"
else
  joblist="S8:480 J5:300"
  smoke_arg=""
fi
arms="J5C J5S"

end() { # end <exit> <reason>: write the marker and exit
  code="$1"; reason="$2"
  if [ "$dry" = "0" ]; then
    mkdir -p "$logs"
    {
      echo "exit=$code"
      echo "reason=$reason"
      date +%Y-%m-%dT%H:%M:%S%z
      echo "build driver=$(cat "$logs/held5_driver.rc" 2>/dev/null || echo none) pilot S8=$(cat "$logs/pilot_S8.rc" 2>/dev/null || echo none) J5=$(cat "$logs/pilot_J5.rc" 2>/dev/null || echo none)"
      for j in $joblist; do
        a="${j%%:*}"
        echo "${a}_${smp} run=$(cat "$logs/${a}_${smp}.run_rc" 2>/dev/null || echo none) pass=$(cat "$logs/pass_${a}_${smp}.rc" 2>/dev/null || echo none)"
      done
      for a in $arms; do
        echo "${a}_${smp} exit=$(cat "$logs/${a}_${smp}.rc" 2>/dev/null || echo none)"
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
    echo "exit=2 reason=refused: $1" > "$logs/screen_run_all.refused.$(date +%s)"
  fi
  exit 2
}

rundir() { if [ "$smoke" = "1" ]; then echo "$1_${smp}_smoke"; else echo "$1_${smp}"; fi; }
jfiles() { # jfiles <arm>: the outputs of an S8 / J5 job (paths relative to the cache)
  if [ "$smoke" = "1" ]; then
    echo "$(rundir "$1") logs/$1_${smp}.log logs/$1_${smp}.done logs/$1_${smp}.run_rc logs/pass_$1_${smp}.log logs/pass_$1_${smp}.rc"
  else
    echo "$1_${smp} passes_${smp}/D4_$1 passes_${smp}/before_$1.npy passes_${smp}/after_$1.npy logs/$1_${smp}.log logs/$1_${smp}.done logs/$1_${smp}.run_rc logs/pass_$1_${smp}.log logs/pass_$1_${smp}.rc"
  fi
}
afiles() { # afiles <arm>: the outputs of a J5C / J5S job
  if [ "$smoke" = "1" ]; then
    echo "$1_${smp}_smoke logs/$1_${smp}.log logs/$1_${smp}.rc logs/$1_${smp}.done"
  else
    echo "$1_${smp} logs/$1_${smp}.log logs/$1_${smp}.rc logs/$1_${smp}.done"
  fi
}
job_finished() { # job_finished <arm>   (S8, J5)
  d="$(cat "$logs/$1_${smp}.done" 2>/dev/null)"
  if [ "$smoke" = "1" ]; then
    { [ "$d" = "run=0 pass=skipped" ] && [ -e "$cache/$(rundir "$1")/run.json" ]; } \
      || { [ "$d" = "run=3 pass=skipped" ] && [ -e "$cache/$(rundir "$1")/capped.json" ]; }
    return
  fi
  [ "$d" = "run=0 pass=0" ] && [ -e "$cache/$1_${smp}/run.json" ] && [ ! -e "$cache/$1_${smp}/capped.json" ] \
    && [ -e "$cache/passes_${smp}/after_$1.npy" ] && [ -e "$cache/passes_${smp}/D4_$1/pass.json" ] \
    && [ -e "$cache/passes_${smp}/before_$1.npy" ]
}
arm_finished() { # arm_finished <arm>   (J5C, J5S)
  [ "$(cat "$logs/$1_${smp}.done" 2>/dev/null)" = "exit=0" ] \
    && [ -e "$cache/$(rundir "$1")/repair.json" ] && [ -e "$cache/$(rundir "$1")/after.npy" ]
}
build_finished() { # the build's outputs exist (the check that the first start of step 20 got wrong
  # by asking for a file that was written later: here only the sample's own files are asked)
  [ -e "$grid" ] && [ -e "$st/${smp}_clean.npy" ]
}

# ---- guard 0: every output must really lie in this cache, not behind a link into another cache
outs="logs"
for j in $joblist; do
  a="${j%%:*}"; outs="$outs $(jfiles "$a")"
done
for a in $arms; do outs="$outs $(afiles "$a")"; done
[ "$smoke" = "0" ] && outs="$outs passes_${smp}"
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
    if job_finished "$a"; then kept="$kept ${a}_${smp}"; else redo="$redo ${a}_${smp}"; fi
    continue
  fi
  for p in $(jfiles "$a"); do
    if [ -e "$cache/$p" ]; then refuse "$cache/$p already exists"; fi
  done
done
for a in $arms; do
  if [ "$resume" = "1" ]; then
    if arm_finished "$a"; then kept="$kept ${a}_${smp}"; else redo="$redo ${a}_${smp}"; fi
  else
    for p in $(afiles "$a"); do
      if [ -e "$cache/$p" ]; then refuse "$cache/$p already exists"; fi
    done
  fi
done

# ---- guard 2: the gate and the freeze
[ -e "$gate" ] || refuse "ar5/screen/gate.json does not exist (the gate has not passed)"
ferr="$(mktemp)"
if ! uv run python scripts/autoresearch/ar5/screen_freeze.py --check > /dev/null 2>"$ferr"; then
  msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
  refuse "${msg:-screen freeze check failed}"
fi
rm -f "$ferr"

# ---- guard 3: nothing of this step (or an arm it builds on) is running
if pgrep -f "autoresearch/ar5/(screen|cmarepair|high500|fit)_(run|pass|launch|build_held4|build_held5|make_held3)\.py" > /dev/null 2>&1; then
  refuse "a screen / cmarepair / high500 / fit run, pass, launch or build process is already running"
fi

# ---- guard 4: inputs
need="$here/benchmarks/autoresearch/ar3/reseed/remedies_pick.json"
if [ "$smoke" = "1" ]; then
  need="$need $grid $st/${smp}_clean.npy $cache/passes_${smp}/before_J5.npy $cache/passes_${smp}/D4_J5/tries.npz"
fi
for p in $need; do
  if [ ! -e "$p" ]; then refuse "input missing: $p"; fi
done

run_cmd() { # run_cmd <arm> <cap>
  case "$1" in
    S8) echo "RESEED_SAMPLE=$smp uv run python scripts/autoresearch/ar5/screen_launch.py fit_run.py --name S8 --cap-min $2$smoke_arg" ;;
    J5) echo "RESEED_SAMPLE=$smp uv run python scripts/autoresearch/ar5/screen_launch.py high500_run.py --cap-min $2$smoke_arg" ;;
  esac
}
pass_cmd() { # pass_cmd <arm>
  case "$1" in
    S8) echo "RESEED_SAMPLE=$smp uv run python scripts/autoresearch/ar5/screen_launch.py fit_pass.py --name S8" ;;
    J5) echo "RESEED_SAMPLE=$smp uv run python scripts/autoresearch/ar5/screen_launch.py high500_pass.py --name J5 --aligned 1" ;;
  esac
}
arm_cmd() { # arm_cmd <arm>
  lim=""
  if [ "$smoke" = "1" ]; then
    if [ "$1" = "J5C" ]; then lim=" --limit 5"; else lim=" --limit 20"; fi
  fi
  case "$1" in
    J5C) echo "SCREEN_PRODUCER=1 RESEED_SAMPLE=$smp uv run python scripts/autoresearch/ar5/screen_run.py --arm J5C --source stored --cap-min 300$lim" ;;
    J5S) echo "SCREEN_PRODUCER=1 RESEED_SAMPLE=$smp uv run python scripts/autoresearch/ar5/screen_run.py --arm J5S --cap-min 300$lim" ;;
  esac
}
pilot_cmd() { # pilot_cmd <arm>: a smoke of an S8 / J5 run on held5, in the pilot cache
  case "$1" in
    S8) echo "SCREEN_CACHE=$pilot_cache RESEED_SAMPLE=$smp uv run python scripts/autoresearch/ar5/screen_launch.py fit_run.py --name S8 --cap-min 5 --smoke 30" ;;
    J5) echo "SCREEN_CACHE=$pilot_cache RESEED_SAMPLE=$smp uv run python scripts/autoresearch/ar5/screen_launch.py high500_run.py --cap-min 5 --smoke 30" ;;
  esac
}
build_cmd() { echo "uv run python scripts/autoresearch/ar5/screen_build_held5.py --resume"; }

if [ "$dry" = "1" ]; then
  echo "dry run: guards passed (cache $cache, at most 3 heavy processes, sample $smp)"
  if [ "$resume" = "1" ]; then
    echo "  resume: KEPT (finished, not rerun):${kept:- none}"
    echo "  resume: RERUN from scratch (leftovers moved to $cache/interrupted_<time>/):${redo:- none}"
    [ -e "$marker" ] && echo "    would move $marker (old marker: $(head -1 "$marker"))"
  fi
  if [ "$smoke" = "0" ]; then
    echo "  $(build_cmd) > $logs/held5_driver.log 2>&1"
    echo "  pilot: $(pilot_cmd S8) > $logs/pilot_S8.log 2>&1"
    echo "  pilot: $(pilot_cmd J5) > $logs/pilot_J5.log 2>&1"
  fi
  for j in $joblist; do
    a="${j%%:*}"; cap="${j#*:}"
    echo "  $(run_cmd "$a" "$cap") > $logs/${a}_${smp}.log 2>&1"
    if [ "$smoke" = "0" ]; then echo "    then, if it finished: $(pass_cmd "$a") > $logs/pass_${a}_${smp}.log 2>&1"; fi
  done
  for a in $arms; do
    echo "    after J5's pass (or at once, in a smoke): $(arm_cmd "$a") > $logs/${a}_${smp}.log 2>&1"
  done
  echo "  then: marker $marker (exit=<code>, reason=<text>)"
  exit 0
fi

mkdir -p "$logs"
# safe to start detached: ignore a hangup, read no terminal, keep this script's own output in the log
trap '' HUP
exec < /dev/null
exec >> "$logs/screen_run_all.out" 2>&1
if [ "$resume" = "1" ]; then
  stamp="$(date +%Y%m%d-%H%M%S)"
  idir="$cache/interrupted_$stamp"
  if [ -e "$marker" ]; then   # an old marker, failed or not, never blocks a resume
    mkdir -p "$idir/logs"; mv "$marker" "$idir/logs/screen_runs.done"
  fi
  for n in $redo; do
    case "$n" in
      J5C_*|J5S_*) fl="$(afiles "${n%_${smp}}")" ;;
      *) fl="$(jfiles "${n%_${smp}}")" ;;
    esac
    for p in $fl; do
      if [ -e "$cache/$p" ]; then
        mkdir -p "$idir/$(dirname "$p")"
        mv "$cache/$p" "$idir/$p"
      fi
    done
  done
  echo "resume $(date +%Y-%m-%dT%H:%M:%S%z) kept:${kept:- none} rerun:${redo:- none}" >> "$logs/screen_run_all.resume"
else
  date +%Y-%m-%dT%H:%M:%S%z > "$logs/screen_run_all.start"
fi
echo "start $(date +%Y-%m-%dT%H:%M:%S%z) sample=$smp smoke=$smoke resume=$resume"

# ---- the held5 build (a finished stage is skipped by the builder itself) and the pilot
if [ "$smoke" = "0" ]; then
  eval "$(build_cmd)" > "$logs/held5_driver.log" 2>&1
  brc=$?
  [ "$brc" = "0" ] && { build_finished || brc=4; }
  echo "$brc" > "$logs/held5_driver.rc"
  if [ "$brc" != "0" ]; then
    date +%Y-%m-%dT%H:%M:%S%z > "$logs/screen_run_all.end"
    end 1 "the held5 build failed (see held5_driver.log; not read for results)"
  fi
  # pilot: two short smoke runs; only exit codes and time are looked at, never their content
  mkdir -p "$pilot_cache/logs"
  t0="$(date +%s)"
  pp=""
  for a in S8 J5; do
    ( eval "$(pilot_cmd "$a")" > "$logs/pilot_$a.log" 2>&1; echo "$?" > "$logs/pilot_$a.rc" ) &
    pp="$pp $!"
  done
  for p in $pp; do wait "$p"; done
  echo "pilot_seconds=$(( $(date +%s) - t0 ))" > "$logs/pilot.time"
  for a in S8 J5; do
    prc="$(cat "$logs/pilot_$a.rc" 2>/dev/null || echo 99)"
    if [ "$prc" != "0" ] && [ "$prc" != "3" ]; then
      date +%Y-%m-%dT%H:%M:%S%z > "$logs/screen_run_all.end"
      end 1 "the held5 pilot failed for $a (exit $prc; see pilot_$a.log; not read for results)"
    fi
    if [ ! -e "$pilot_cache/${a}_${smp}_smoke" ]; then
      date +%Y-%m-%dT%H:%M:%S%z > "$logs/screen_run_all.end"
      end 1 "the held5 pilot of $a left no folder (not read for results)"
    fi
  done
fi

run_arm() { # run_arm <arm>   (J5C or J5S: a subshell)
  a="$1"
  case " $kept " in *" ${a}_${smp} "*) return ;; esac
  eval "$(arm_cmd "$a")" > "$logs/${a}_${smp}.log" 2>&1
  arc=$?
  echo "$arc" > "$logs/${a}_${smp}.rc"
  echo "exit=$arc" > "$logs/${a}_${smp}.done"
}

run_job() { # run_job <arm> <cap>   (a subshell: run, then pass, then, for J5, the two arms)
  a="$1"; cap="$2"
  case " $kept " in *" ${a}_${smp} "*) ;; *)
    eval "$(run_cmd "$a" "$cap")" > "$logs/${a}_${smp}.log" 2>&1
    rc=$?
    echo "$rc" > "$logs/${a}_${smp}.run_rc"
    d="$cache/$(rundir "$a")"
    if [ "$rc" = "0" ] && [ -e "$d/run.json" ] && [ ! -e "$d/capped.json" ] && [ "$smoke" = "0" ]; then
      eval "$(pass_cmd "$a")" > "$logs/pass_${a}_${smp}.log" 2>&1
      prc=$?
      echo "$prc" > "$logs/pass_${a}_${smp}.rc"
    else
      prc=skipped
      echo skipped > "$logs/pass_${a}_${smp}.rc"
    fi
    echo "run=$rc pass=$prc" > "$logs/${a}_${smp}.done"
  ;; esac
  if [ "$a" = "J5" ] && [ "$smoke" = "0" ]; then
    if job_finished J5; then
      ( run_arm J5C ) &
      pa=$!
      ( run_arm J5S ) &
      pb=$!
      wait "$pa"; wait "$pb"
    else
      for x in $arms; do echo skipped > "$logs/${x}_${smp}.rc"; echo "exit=skipped" > "$logs/${x}_${smp}.done"; done
    fi
  fi
}

pids=""
for j in $joblist; do
  a="${j%%:*}"; cap="${j#*:}"
  run_job "$a" "$cap" &
  pids="$pids $!"
done
if [ "$smoke" = "1" ]; then   # the arms on the linked stored maps, at once (S8 and J5 smoke + 2)
  for a in $arms; do
    ( run_arm "$a" ) &
    pids="$pids $!"
  done
fi
for p in $pids; do wait "$p"; done

capped=""
fail=0
for j in $joblist; do
  a="${j%%:*}"
  rc="$(cat "$logs/${a}_${smp}.run_rc" 2>/dev/null || echo 99)"
  prc="$(cat "$logs/pass_${a}_${smp}.rc" 2>/dev/null || echo 99)"
  if [ "$rc" = "3" ]; then
    capped="${a}_${smp}"
    [ "$smoke" = "1" ] || fail=1
  elif [ "$rc" != "0" ]; then fail=1
  fi
  if [ "$smoke" = "0" ] && [ "$prc" != "0" ]; then fail=1; fi
done
for a in $arms; do
  arc="$(cat "$logs/${a}_${smp}.rc" 2>/dev/null || echo 99)"
  if [ "$arc" = "3" ]; then capped="${a}_${smp}"; fail=1
  elif [ "$arc" != "0" ]; then fail=1
  fi
done
reason="all runs, passes and arms finished"
[ "$smoke" = "1" ] && reason="smoke finished (no pass; J5C and J5S with a limit)"
date +%Y-%m-%dT%H:%M:%S%z > "$logs/screen_run_all.end"
if [ "$smoke" = "0" ] && [ -n "$capped" ]; then end 1 "STOP RULE: a run reached its cap: $capped"; fi
if [ "$fail" = "1" ]; then
  end 1 "a run, a pass or an arm failed or reached its cap (see the rc files; not read for results)"
fi
end 0 "$reason"
