#!/bin/bash
# AR5 step 22, phase A steps 5 to 7: the live runs of arm J5D and arm J5SD, and the clean test on
# held6. A FOREGROUND script: the main session starts it from icenine_py/ detached:
#
#   cd icenine_py && nohup caffeinate -i scripts/autoresearch/ar5/defer_run_all.sh --stage gate \
#       --resume > /dev/null 2>&1 &
#   ... --stage rest ...    only after the gate was judged and logs/gate_passed was written
#   ... --dry-run           guards and commands only (nothing written)
#   ... --cache DIR         another cache folder (tests)
#   ... --smoke             test form (needs --sample, a disclosed sample whose inputs are linked
#                           in the cache): short jobs, no pass: J5D / S8 / J5 with --smoke 30 and a
#                           cap of --smoke-cap M minutes (default 3), then the screen arm
#                           (--limit 20) on the linked STORED maps (stored_passes_<sample>). The
#                           cache defaults to <defer>/cache_smoke. Exit 3 (the cap) is the smoke
#                           jobs' normal end. In stage rest the sample stands in for held6 (the
#                           build command is "true"; the build check and the pilot's exit-code
#                           logic still run) and also is the only J5SD sample.
#   ... --sample S          the sample of a smoke (default mid)
#   ... --resume            after an interruption or a FAILED earlier attempt: a job whose .done
#                           says it ended 0 with its outputs is kept; the leftovers of every other
#                           job and an old marker (failed or not) are MOVED to
#                           <cache>/interrupted_<time>/ and those jobs run from scratch. Safe on a
#                           first start too.
#
# Stages (at most 3 heavy processes at once; a job is one process at a time):
#   gate  live J5D on mid and held (2 jobs). Each job: the J5D run (defer_launch.py defer_run.py),
#         then J5's boundary pass on its map (high500_pass.py --name J5D --aligned 1), then step
#         21's screen pass and wide repair on that map (screen_run.py --arm J5S, run from a
#         sub-cache <cache>/j5sd whose before_J5.npy and D4_J5/ are links to J5D's pass files, so
#         that step 21's script is imported unchanged): arm J5SD = <cache>/j5sd/J5S_<s>. Marker
#         logs/defer_gate.done ("exit=<code>", "reason=<text>", the end time, one line a job).
#   rest  refuses unless logs/gate_passed exists (written by the implementer after g2 and g3 were
#         judged from the gate stage). Then, in this queue order with a pool of 3: the held6 build
#         and pilot; J5SD live on held5 and held4; S8 on held6 (run, pass; cap 480 minutes); J5SD
#         on held3; J5 on held6 (run, pass, then arm J5S from its map); J5SD on held6; J5SD on
#         grad and big. The held6 jobs wait for the build and pilot (a failed build or pilot ends
#         them: nothing of held6 is run). The pilot is S8 and J5 on held6 with --smoke 30 and a cap
#         of 5 minutes each, one after the other; only exit codes (0 or 3 = the cap, normal for a
#         smoke) and the existence of their folders are looked at, never their content. Marker
#         logs/defer_runs.done.
# Caps (--cap-min): S8 480, everything else 300. A run at its cap writes capped.json and exits 3;
# its pass is skipped and the marker says so (stop rule).
#
# Guards (exit 2, nothing started; a note logs/defer_run_all.refused.<time> holds the reason): an
# output of a job or the stage's marker exists (not with --resume); ar5/defer/g0.json and
# replay.json are missing (gate stage) or logs/gate_passed is missing (rest stage);
# defer_freeze.py --check fails (frozen.json committed, clean, every script and library hash
# agrees); a defer / screen / cmarepair / high500 / fit run, pass, launch or build process is
# running; an input is missing; an output's real path lies outside the cache (stallhit_cache_
# guard.py). Nothing of held2 is used (defer_common refuses RESEED_SAMPLE=held2).
#
# Files, in <cache>/logs (cache = benchmarks/autoresearch/ar5/defer/cache, untracked):
#   defer_run_all.out (this script's own output), J5D_<s>.log / .run_rc, pass_J5D_<s>.log / .rc,
#   J5SD_<s>.log / .rc, J5D_<s>.done ("run=<code> pass=<code|skipped> arm=<code|skipped>"),
#   the same for S8_<s6> / J5_<s6> (S8: run, pass; J5: run, pass, arm = J5S), held6_<stage>.log /
#   .rc (the build), pilot_<arm>.rc, pilot.time, held6_pilot.ok, defer_run_all.start / .end and the
#   MARKER. The script's exit code is the marker's. NOBODY reads held6's logs or maps before the
#   marker defer_runs.done says all of its runs and repairs ended (plan stop rule; defer_guard).
#   Wall times are contended and are not timing results.
# Plain bash 3.2 (macOS default).
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${SCREEN_CACHE:-$here/benchmarks/autoresearch/ar5/defer/cache}"
dry=0
resume=0
smoke=0
smoke_cap=3
cache_given=0
stage=""
smp=mid
maxheavy=3
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1"; cache_given=1 ;;
    --smoke) smoke=1 ;;
    --smoke-cap) shift; smoke_cap="$1" ;;
    --sample) shift; smp="$1" ;;
    --stage) shift; stage="$1" ;;
    --resume) resume=1 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
case "$stage" in gate|rest) ;; *) echo "--stage gate|rest is required (got '$stage')" >&2; exit 2 ;; esac
case "$smp" in held2*|*held2*) echo "held2 is not used by this step" >&2; exit 2 ;; esac
case "$smp" in held5|held4|held3|held|mid|grad|big) ;; *) echo "unknown sample: $smp" >&2; exit 2 ;; esac
if [ "$smoke" = "1" ] && [ "$cache_given" = "0" ]; then
  cache="$here/benchmarks/autoresearch/ar5/defer/cache_smoke"
fi
case "$cache" in /*) ;; *) cache="$here/$cache" ;; esac   # the jobs read it from other places
logs="$cache/logs"
if [ "$stage" = "gate" ]; then marker="$logs/defer_gate.done"; else marker="$logs/defer_runs.done"; fi
pilot_cache="${cache}_pilot"
pilot_cap=5
if [ "$smoke" = "1" ]; then pilot_cap="$smoke_cap"; fi
export SCREEN_CACHE="$cache"
export DEFER_CACHE="$cache"
unset VIRTUAL_ENV
sub="$cache/j5sd"   # the sub-cache of the J5SD arm (links to J5D's pass files; its own outputs)

if [ "$smoke" = "1" ]; then
  smoke_arg=" --smoke 30"
  cap_s8="$smoke_cap"; cap_rest="$smoke_cap"
  lim=" --limit 20"
else
  smoke_arg=""
  cap_s8=480; cap_rest=300
  lim=""
fi
# the samples of each stage (a smoke: the one sample, standing in for held6 in stage rest)
if [ "$stage" = "gate" ]; then
  if [ "$smoke" = "1" ]; then dsamples="$smp"; else dsamples="mid held"; fi
  s6=""
else
  if [ "$smoke" = "1" ]; then dsamples="$smp"; s6="$smp"; else dsamples="held5 held4 held3 grad big held6"; s6="held6"; fi
fi

end() { # end <exit> <reason>: write the marker and exit
  code="$1"; reason="$2"
  if [ "$dry" = "0" ]; then
    mkdir -p "$logs"
    {
      echo "exit=$code"
      echo "reason=$reason"
      date +%Y-%m-%dT%H:%M:%S%z
      for s in $dsamples; do
        echo "J5D_${s} $(cat "$logs/J5D_${s}.done" 2>/dev/null || echo none)"
      done
      if [ -n "$s6" ]; then
        echo "build=$(cat "$logs/held6_build.rc" 2>/dev/null || echo none) pilot=$(cat "$logs/held6_pilot.ok" 2>/dev/null || echo none)"
        echo "S8_${s6} $(cat "$logs/S8_${s6}.done" 2>/dev/null || echo none)"
        echo "J5_${s6} $(cat "$logs/J5_${s6}.done" 2>/dev/null || echo none)"
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
    echo "exit=2 reason=refused: $1" > "$logs/defer_run_all.refused.$(date +%s)"
  fi
  exit 2
}

# ---- the gate_passed file comes first for stage rest (a plain file test, before any other work)
if [ "$stage" = "rest" ] && [ ! -e "$logs/gate_passed" ]; then
  refuse "logs/gate_passed does not exist (the gate has not been judged and passed): $logs/gate_passed"
fi

rundir() { if [ "$smoke" = "1" ]; then echo "$1_$2_smoke"; else echo "$1_$2"; fi; }
# ---- outputs of each job (paths relative to the cache)
dfiles() { # D job on sample $1
  echo "$(rundir J5D "$1") passes_$1/D4_J5D passes_$1/before_J5D.npy passes_$1/after_J5D.npy j5sd/passes_$1 $(rundir j5sd/J5S "$1") logs/J5D_$1.log logs/J5D_$1.run_rc logs/pass_J5D_$1.log logs/pass_J5D_$1.rc logs/J5SD_$1.log logs/J5SD_$1.rc logs/J5D_$1.done"
}
bfiles() { # B job: arm $1 (S8 or J5) on sample $2
  echo "$(rundir "$1" "$2") passes_$2/D4_$1 passes_$2/before_$1.npy passes_$2/after_$1.npy logs/$1_$2.log logs/$1_$2.run_rc logs/pass_$1_$2.log logs/pass_$1_$2.rc logs/$1_$2.done"
  if [ "$1" = "J5" ]; then echo "$(rundir J5S "$2") logs/J5S_$2.log logs/J5S_$2.rc"; fi
}
d_finished() { # d_finished <sample>
  d="$(cat "$logs/J5D_$1.done" 2>/dev/null)"
  if [ "$smoke" = "1" ]; then
    { [ "$d" = "run=0 pass=skipped arm=0" ] && [ -e "$cache/$(rundir J5D "$1")/run.json" ] && [ -e "$sub/$(rundir J5S "$1")/repair.json" ]; } \
      || { [ "$d" = "run=3 pass=skipped arm=0" ] && [ -e "$cache/$(rundir J5D "$1")/capped.json" ] && [ -e "$sub/$(rundir J5S "$1")/repair.json" ]; }
    return
  fi
  [ "$d" = "run=0 pass=0 arm=0" ] && [ -e "$cache/J5D_$1/run.json" ] && [ ! -e "$cache/J5D_$1/capped.json" ] \
    && [ -e "$cache/passes_$1/after_J5D.npy" ] && [ -e "$cache/passes_$1/D4_J5D/pass.json" ] \
    && [ -e "$cache/passes_$1/before_J5D.npy" ] && [ -e "$sub/J5S_$1/repair.json" ] && [ -e "$sub/J5S_$1/after.npy" ]
}
b_finished() { # b_finished <arm> <sample>
  d="$(cat "$logs/$1_$2.done" 2>/dev/null)"
  if [ "$smoke" = "1" ]; then
    { [ "$d" = "run=0 pass=skipped arm=skipped" ] && [ -e "$cache/$(rundir "$1" "$2")/run.json" ]; } \
      || { [ "$d" = "run=3 pass=skipped arm=skipped" ] && [ -e "$cache/$(rundir "$1" "$2")/capped.json" ]; }
    return
  fi
  if [ "$1" = "J5" ]; then want="run=0 pass=0 arm=0"; else want="run=0 pass=0 arm=skipped"; fi
  [ "$d" = "$want" ] && [ -e "$cache/$1_$2/run.json" ] && [ ! -e "$cache/$1_$2/capped.json" ] \
    && [ -e "$cache/passes_$2/after_$1.npy" ] && [ -e "$cache/passes_$2/D4_$1/pass.json" ] \
    && [ -e "$cache/passes_$2/before_$1.npy" ] \
    && { [ "$1" != "J5" ] || { [ -e "$cache/J5S_$2/repair.json" ] && [ -e "$cache/J5S_$2/after.npy" ]; }; }
}
build_finished() { # the held6 (or smoke sample's) own files exist; nothing is opened
  [ -e "$here/benchmarks/autoresearch/ar3/$s6/gradient/sample/gradient_grid.mic" ] \
    && [ -e "$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/${s6}_clean.npy" ]
}

# ---- guard 0: every output must really lie in this cache, not behind a link into another cache
outs="logs j5sd"
for s in $dsamples; do outs="$outs $(dfiles "$s")"; done
if [ -n "$s6" ]; then
  outs="$outs $(bfiles S8 "$s6") $(bfiles J5 "$s6")"
  [ "$smoke" = "0" ] && outs="$outs passes_$s6"
fi
if ! uv run python scripts/autoresearch/ar5/stallhit_cache_guard.py --cache "$cache" $outs; then
  echo "refusing to start: an output's real path lies outside the cache $cache" >&2
  exit 2
fi

# ---- guard 1: no output exists (--resume: finished jobs are kept, leftovers and an old marker are
# moved aside later)
kept=""; redo=""
if [ -e "$marker" ] && [ "$resume" = "0" ]; then refuse "the marker $marker already exists (use --resume)"; fi
for s in $dsamples; do
  if [ "$resume" = "1" ]; then
    if d_finished "$s"; then kept="$kept D_$s"; else redo="$redo D_$s"; fi
  else
    for p in $(dfiles "$s"); do
      if [ -e "$cache/$p" ]; then refuse "$cache/$p already exists"; fi
    done
  fi
done
if [ -n "$s6" ]; then
  for a in S8 J5; do
    if [ "$resume" = "1" ]; then
      if b_finished "$a" "$s6"; then kept="$kept B_${a}_$s6"; else redo="$redo B_${a}_$s6"; fi
    else
      for p in $(bfiles "$a" "$s6"); do
        if [ -e "$cache/$p" ]; then refuse "$cache/$p already exists"; fi
      done
    fi
  done
fi

# ---- guard 2: the evidence of the earlier gates, and the freeze (our own, every file pinned)
if [ "$stage" = "gate" ] && [ "$smoke" = "0" ]; then
  [ -e "$here/benchmarks/autoresearch/ar5/defer/replay.json" ] || refuse "ar5/defer/replay.json does not exist (g1 not run)"
  [ -e "$here/benchmarks/autoresearch/ar5/defer/g0.json" ] || refuse "ar5/defer/g0.json does not exist (g0 not run)"
fi
ferr="$(mktemp)"
if ! uv run python scripts/autoresearch/ar5/defer_freeze.py --check > /dev/null 2>"$ferr"; then
  msg="$(tail -1 "$ferr" 2>/dev/null)"; rm -f "$ferr"
  refuse "${msg:-defer freeze check failed}"
fi
rm -f "$ferr"

# ---- guard 3: nothing of this step (or an arm it builds on) is running
if pgrep -f "autoresearch/ar5/(defer|screen|cmarepair|high500|fit)_(run|pass|launch|build_held4|build_held5|build_held6|make_held3)\.py" > /dev/null 2>&1; then
  refuse "a defer / screen / cmarepair / high500 / fit run, pass, launch or build process is already running"
fi

# ---- guard 4: inputs (held6's own files do not exist before its build; nothing is opened)
need="$here/benchmarks/autoresearch/ar3/reseed/remedies_pick.json"
st="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
for s in $dsamples; do
  if [ "$s" = "$s6" ] && [ "$smoke" = "0" ]; then continue; fi   # held6's files do not exist yet
  if [ "$s" = "grad" ]; then g="$here/benchmarks/autoresearch/ar3/gradient/sample/gradient_grid.mic"
  else g="$here/benchmarks/autoresearch/ar3/$s/gradient/sample/gradient_grid.mic"; fi
  need="$need $g $st/${s}_clean.npy"
  if [ "$smoke" = "1" ]; then
    need="$need $cache/stored_passes_$s/before_J5.npy $cache/stored_passes_$s/D4_J5/tries.npz $cache/stored_passes_$s/D4_J5/pass.json"
  fi
done
for p in $need; do
  if [ ! -e "$p" ]; then refuse "input missing: $p"; fi
done

# ---- the commands
ev() { # ev <sample> [subcache]: the environment of a job
  echo "RESEED_SAMPLE=$1 DEFER_CACHE=$cache SCREEN_CACHE=${2:-$cache}"
}
L="uv run python scripts/autoresearch/ar5/defer_launch.py"
d_run_cmd() { echo "$(ev "$1") $L defer_run.py --cap-min $cap_rest$smoke_arg"; }
d_pass_cmd() { echo "$(ev "$1") $L high500_pass.py --name J5D --aligned 1"; }
d_arm_cmd() { echo "SCREEN_PRODUCER=1 $(ev "$1" "$sub") $L screen_run.py --arm J5S --cap-min $cap_rest$lim"; }
b_run_cmd() { # b_run_cmd <arm> <sample> <cap> [cache override]
  c="${4:-$cache}"
  if [ "$1" = "S8" ]; then echo "RESEED_SAMPLE=$2 DEFER_CACHE=$cache SCREEN_CACHE=$c $L fit_run.py --name S8 --cap-min $3$smoke_arg"
  else echo "RESEED_SAMPLE=$2 DEFER_CACHE=$cache SCREEN_CACHE=$c $L high500_run.py --cap-min $3$smoke_arg"; fi
}
b_pass_cmd() {
  if [ "$1" = "S8" ]; then echo "$(ev "$2") $L fit_pass.py --name S8"
  else echo "$(ev "$2") $L high500_pass.py --name J5 --aligned 1"; fi
}
b_arm_cmd() { echo "SCREEN_PRODUCER=1 $(ev "$1") $L screen_run.py --arm J5S --cap-min $cap_rest"; }
build_cmd() { if [ "$smoke" = "1" ]; then echo "true"; else echo "uv run python scripts/autoresearch/ar5/defer_build_held6.py --resume"; fi; }

if [ "$dry" = "1" ]; then
  echo "dry run: guards passed (stage $stage, cache $cache, at most $maxheavy heavy processes, smoke $smoke)"
  if [ "$resume" = "1" ]; then
    echo "  resume: KEPT (finished, not rerun):${kept:- none}"
    echo "  resume: RERUN from scratch (leftovers moved to $cache/interrupted_<time>/):${redo:- none}"
    [ -e "$marker" ] && echo "    would move $marker (old marker: $(head -1 "$marker"))"
  fi
  if [ -n "$s6" ]; then
    echo "  build: $(build_cmd) > $logs/held6_build.log 2>&1"
    echo "  pilot S8: SCREEN_CACHE=$pilot_cache $(b_run_cmd S8 "$s6" "$pilot_cap" "$pilot_cache" | sed 's/$/ --smoke 30/')"
    echo "  pilot J5: SCREEN_CACHE=$pilot_cache $(b_run_cmd J5 "$s6" "$pilot_cap" "$pilot_cache" | sed 's/$/ --smoke 30/')"
  fi
  for s in $dsamples; do
    echo "  D $s: $(d_run_cmd "$s") > $logs/J5D_$s.log 2>&1"
    echo "    then: $(d_pass_cmd "$s") > $logs/pass_J5D_$s.log 2>&1"
    echo "    then (sub-cache links to J5D's pass files): $(d_arm_cmd "$s") > $logs/J5SD_$s.log 2>&1"
  done
  if [ -n "$s6" ]; then
    echo "  B S8 $s6: $(b_run_cmd S8 "$s6" "$cap_s8") ; then $(b_pass_cmd S8 "$s6")"
    echo "  B J5 $s6: $(b_run_cmd J5 "$s6" "$cap_rest") ; then $(b_pass_cmd J5 "$s6") ; then $(b_arm_cmd "$s6")"
  fi
  echo "  then: marker $marker (exit=<code>, reason=<text>)"
  exit 0
fi

mkdir -p "$logs" "$sub"
# safe to start detached: ignore a hangup, read no terminal, keep this script's own output in the log
trap '' HUP
exec < /dev/null
exec >> "$logs/defer_run_all.out" 2>&1
if [ "$resume" = "1" ]; then
  stamp="$(date +%Y%m%d-%H%M%S)"
  idir="$cache/interrupted_$stamp"
  if [ -e "$marker" ]; then   # an old marker, failed or not, never blocks a resume
    mkdir -p "$idir/logs"; mv "$marker" "$idir/logs/$(basename "$marker")"
  fi
  for n in $redo; do
    case "$n" in
      D_*) fl="$(dfiles "${n#D_}")" ;;
      B_S8_*) fl="$(bfiles S8 "${n#B_S8_}")" ;;
      B_J5_*) fl="$(bfiles J5 "${n#B_J5_}")" ;;
    esac
    for p in $fl; do
      if [ -e "$cache/$p" ] || [ -L "$cache/$p" ]; then
        mkdir -p "$idir/$(dirname "$p")"
        mv "$cache/$p" "$idir/$p"
      fi
    done
  done
  echo "resume $(date +%Y-%m-%dT%H:%M:%S%z) kept:${kept:- none} rerun:${redo:- none}" >> "$logs/defer_run_all.resume"
else
  date +%Y-%m-%dT%H:%M:%S%z > "$logs/defer_run_all.start"
fi
echo "start $(date +%Y-%m-%dT%H:%M:%S%z) stage=$stage smoke=$smoke resume=$resume sample=$smp"

mk_sub() { # mk_sub <sample>: the sub-cache's inputs are links to J5D's pass files (a smoke: to the
  # stored J5 pass files), named as step 21's script reads J5's
  s="$1"
  if [ "$smoke" = "1" ]; then src="$cache/stored_passes_$s"; nm=J5; else src="$cache/passes_$s"; nm=J5D; fi
  mkdir -p "$sub/passes_$s/D4_J5"
  ln -sfn "$src/before_$nm.npy" "$sub/passes_$s/before_J5.npy"
  for f in tries.npz pass.json; do ln -sfn "$src/D4_$nm/$f" "$sub/passes_$s/D4_J5/$f"; done
}

d_job() { # d_job <sample>: the J5D run, J5's pass on its map, then the screen pass and repair
  s="$1"
  case " $kept " in *" D_$s "*) return ;; esac
  eval "$(d_run_cmd "$s")" > "$logs/J5D_$s.log" 2>&1
  rc=$?
  echo "$rc" > "$logs/J5D_$s.run_rc"
  d="$cache/$(rundir J5D "$s")"
  prc=skipped; arc=skipped
  if [ "$smoke" = "1" ]; then   # a smoke: the arm runs on the stored maps whatever the smoke run's end
    echo skipped > "$logs/pass_J5D_$s.rc"
    mk_sub "$s"
    eval "$(d_arm_cmd "$s")" > "$logs/J5SD_$s.log" 2>&1
    arc=$?
    echo "$arc" > "$logs/J5SD_$s.rc"
  elif [ "$rc" = "0" ] && [ -e "$d/run.json" ] && [ ! -e "$d/capped.json" ]; then
    eval "$(d_pass_cmd "$s")" > "$logs/pass_J5D_$s.log" 2>&1
    prc=$?
    echo "$prc" > "$logs/pass_J5D_$s.rc"
    if [ "$prc" = "0" ]; then
      mk_sub "$s"
      eval "$(d_arm_cmd "$s")" > "$logs/J5SD_$s.log" 2>&1
      arc=$?
      echo "$arc" > "$logs/J5SD_$s.rc"
    fi
  else
    echo skipped > "$logs/pass_J5D_$s.rc"
  fi
  echo "run=$rc pass=$prc arm=$arc" > "$logs/J5D_$s.done"
}

b_job() { # b_job <arm> <sample> <cap>: S8 or J5 on the new sample, its pass, (J5: arm J5S)
  a="$1"; s="$2"; cap="$3"
  case " $kept " in *" B_${a}_$s "*) return ;; esac
  eval "$(b_run_cmd "$a" "$s" "$cap")" > "$logs/${a}_$s.log" 2>&1
  rc=$?
  echo "$rc" > "$logs/${a}_$s.run_rc"
  d="$cache/$(rundir "$a" "$s")"
  prc=skipped; arc=skipped
  if [ "$rc" = "0" ] && [ -e "$d/run.json" ] && [ ! -e "$d/capped.json" ] && [ "$smoke" = "0" ]; then
    eval "$(b_pass_cmd "$a" "$s")" > "$logs/pass_${a}_$s.log" 2>&1
    prc=$?
    echo "$prc" > "$logs/pass_${a}_$s.rc"
    if [ "$a" = "J5" ] && [ "$prc" = "0" ]; then
      eval "$(b_arm_cmd "$s")" > "$logs/J5S_$s.log" 2>&1
      arc=$?
      echo "$arc" > "$logs/J5S_$s.rc"
    fi
  else
    echo skipped > "$logs/pass_${a}_$s.rc"
  fi
  echo "run=$rc pass=$prc arm=$arc" > "$logs/${a}_$s.done"
}

bp_job() { # the held6 build, then the pilot (S8 and J5 smokes, one after the other)
  rm -f "$logs/held6_pilot.ok" "$logs/held6_pilot.failed"
  eval "$(build_cmd)" > "$logs/held6_driver.log" 2>&1
  brc=$?
  [ "$brc" = "0" ] && { build_finished || brc=4; }
  echo "$brc" > "$logs/held6_build.rc"
  if [ "$brc" != "0" ]; then echo "build failed ($brc)" > "$logs/held6_pilot.failed"; return; fi
  mkdir -p "$pilot_cache/logs"
  t0="$(date +%s)"
  for a in S8 J5; do
    eval "SCREEN_CACHE=$pilot_cache $(b_run_cmd "$a" "$s6" "$pilot_cap" "$pilot_cache" | sed 's/$/ --smoke 30/')" > "$logs/pilot_$a.log" 2>&1
    prc=$?
    echo "$prc" > "$logs/pilot_$a.rc"
    if { [ "$prc" != "0" ] && [ "$prc" != "3" ]; } || [ ! -e "$pilot_cache/${a}_${s6}_smoke" ]; then
      echo "pilot $a failed (exit $prc or no folder; not read for results)" > "$logs/held6_pilot.failed"
      echo "pilot_seconds=$(( $(date +%s) - t0 ))" > "$logs/pilot.time"
      return
    fi
  done
  echo "pilot_seconds=$(( $(date +%s) - t0 ))" > "$logs/pilot.time"
  echo ok > "$logs/held6_pilot.ok"
}

nrun() { jobs -rp | wc -l | tr -d ' '; }
pool() { # pool <function> [args]: start when fewer than $maxheavy jobs run
  while [ "$(nrun)" -ge "$maxheavy" ]; do sleep "$([ "$smoke" = "1" ] && echo 2 || echo 30)"; done
  "$@" &
}

held6_ok=1
if [ "$stage" = "gate" ]; then
  for s in $dsamples; do pool d_job "$s"; done
else
  bp_pid=""
  if [ -n "$s6" ]; then pool bp_job; bp_pid=$!; fi
  start_held6() { # wait for the build and pilot (a short job), then start the held6 jobs' turn
    [ -n "$bp_pid" ] && wait "$bp_pid" 2>/dev/null
    bp_pid=""
    [ -e "$logs/held6_pilot.ok" ] || held6_ok=0
  }
  first=1
  # queue order: J5SD held5, held4, S8 held6, J5SD held3, J5 held6, J5SD held6, J5SD grad, big
  queue="D:held5 D:held4 B8:$s6 D:held3 B5:$s6 D6:$s6 D:grad D:big"
  if [ "$smoke" = "1" ]; then queue="B8:$s6 B5:$s6 D6:$s6"; fi
  for q in $queue; do
    k="${q%%:*}"; s="${q#*:}"
    case "$k" in
      D) pool d_job "$s" ;;
      B8|B5|D6)
        if [ "$first" = "1" ]; then start_held6; first=0; fi
        [ "$held6_ok" = "1" ] || continue
        case "$k" in
          B8) pool b_job S8 "$s" "$cap_s8" ;;
          B5) pool b_job J5 "$s" "$cap_rest" ;;
          D6) pool d_job "$s" ;;
        esac ;;
    esac
  done
  [ "$first" = "1" ] && start_held6
fi
wait

# ---- the verdict
capped=""
fail=0
for s in $dsamples; do
  rc="$(cat "$logs/J5D_$s.run_rc" 2>/dev/null || echo 99)"
  prc="$(cat "$logs/pass_J5D_$s.rc" 2>/dev/null || echo 99)"
  arc="$(cat "$logs/J5SD_$s.rc" 2>/dev/null || echo skipped)"
  case " $kept " in *" D_$s "*) continue ;; esac
  if [ "$rc" = "3" ]; then capped="J5D_$s"; [ "$smoke" = "1" ] || fail=1
  elif [ "$rc" != "0" ]; then fail=1; fi
  if [ "$smoke" = "0" ] && { [ "$prc" != "0" ] || [ "$arc" != "0" ]; }; then fail=1; fi
  if [ "$smoke" = "1" ] && [ "$arc" != "0" ]; then fail=1; fi
  if [ "$arc" = "3" ]; then capped="J5SD_$s"; fail=1; fi
done
if [ -n "$s6" ]; then
  [ "$held6_ok" = "1" ] || fail=1
  if [ "$held6_ok" = "1" ]; then
    for a in S8 J5; do
      case " $kept " in *" B_${a}_$s6 "*) continue ;; esac
      rc="$(cat "$logs/${a}_${s6}.run_rc" 2>/dev/null || echo 99)"
      prc="$(cat "$logs/pass_${a}_${s6}.rc" 2>/dev/null || echo 99)"
      if [ "$rc" = "3" ]; then capped="${a}_$s6"; [ "$smoke" = "1" ] || fail=1
      elif [ "$rc" != "0" ]; then fail=1; fi
      if [ "$smoke" = "0" ] && [ "$prc" != "0" ]; then fail=1; fi
    done
    if [ "$smoke" = "0" ]; then
      case " $kept " in *" B_J5_$s6 "*) ;; *)
        [ "$(cat "$logs/J5S_${s6}.rc" 2>/dev/null || echo 99)" = "0" ] || fail=1 ;;
      esac
    fi
  fi
fi
reason="all runs, passes and arms finished"
[ "$smoke" = "1" ] && reason="smoke finished (no pass; the screen arm with a limit on the stored maps)"
date +%Y-%m-%dT%H:%M:%S%z > "$logs/defer_run_all.end"
if [ "$smoke" = "0" ] && [ -n "$capped" ]; then end 1 "STOP RULE: a run reached its cap: $capped"; fi
if [ -n "$s6" ] && [ "$held6_ok" = "0" ]; then end 1 "the held6 build or pilot failed (see held6_driver.log, held6_pilot.failed; held6 not read)"; fi
if [ "$fail" = "1" ]; then
  end 1 "a run, a pass or an arm failed or reached its cap (see the rc files; not read for results)"
fi
end 0 "$reason"
