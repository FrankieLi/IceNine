#!/bin/bash
# AR5 step 26 (task ar5z): the exactness and timing jobs, one at a time, detached.
#
#   scripts/autoresearch/ar5/ix_run_all.sh --stage e1|live [--dry-run] [--cache DIR]
#       [--samples "mid held"] [--allow-contended] [--tag TAG]
#
#   --stage e1    ix_exact.py on each sample of --samples (default: mid held held3 held4 held5
#                 held6 grad big), one after the other: the old and the new indexer on the
#                 exactness voxels (gate e1) and, on mid and held, the alternating-block seconds
#                 (gate w1). With --require-quiet semantics unless --allow-contended: a refused
#                 timing_preflight ends the job without a run (exit 3).
#   --stage live  the four live J5SD runs, one at a time, in the order: old mid, new mid, old
#                 held, new held (gates e2 and w2). Each run: the J5D run (ix_launch.py
#                 defer_run.py), J5's boundary pass on its map (high500_pass.py --name J5D
#                 --aligned 1), then step 21's screen pass and repair (screen_run.py --arm J5S) on
#                 a sub-cache of links, as step 22's defer_run_all.sh did. Caches <cache>/live_old
#                 and <cache>/live_new. Before each run: timing_preflight --require, its JSON kept
#                 next to the run; refused and not --allow-contended: no run (exit 3).
#
# Writes <cache>/logs/<job>.log and <job>.done ("exit=<code>"); refuses if an output exists.
# One thread (OMP/MKL = 1). Never runs two jobs at once. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${IX_CACHE:-$here/benchmarks/autoresearch/ar5/ix/cache}"
dry=0
stage=""
samples=""
allow=0
tag=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    --stage) shift; stage="$1" ;;
    --samples) shift; samples="$1" ;;
    --allow-contended) allow=1 ;;
    --tag) shift; tag="$1" ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
case "$stage" in e1|live) ;; *) echo "--stage e1|live is required" >&2; exit 2 ;; esac
case "$cache" in /*) ;; *) cache="$here/$cache" ;; esac
case "$cache$samples$tag" in *held2*|*fine7*|*held7*) echo "held2 / fine7 / held7 are not used by this step" >&2; exit 2 ;; esac
if [ -z "$samples" ]; then
  if [ "$stage" = "e1" ]; then samples="mid held held3 held4 held5 held6 grad big"; else samples="mid held"; fi
fi
for s in $samples; do
  case "$s" in mid|held|held3|held4|held5|held6|grad|big) ;; *) echo "unknown sample: $s" >&2; exit 2 ;; esac
done
logs="$cache/logs"
ixd="$here/benchmarks/autoresearch/ar5/ix"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
unset VIRTUAL_ENV
L="uv run python scripts/autoresearch/ar5/ix_launch.py"
pf="uv run python ../scripts/dev/timing_preflight.py"

quiet_check() { # quiet_check <json>: 0 when the preflight passes (JSON kept either way)
  if $pf --require --json "$1" > /dev/null 2>&1; then return 0; fi
  return 1
}

e1_job() { # e1_job <sample>
  s="$1"; job="ix_exact_${s}${tag:+_$tag}"
  [ -e "$logs/$job.done" ] && { echo "refusing: $logs/$job.done exists" >&2; return 2; }
  if [ "$allow" = "0" ]; then
    if ! quiet_check "$ixd/preflight_$job.json"; then echo "exit=3 refused by the preflight" > "$logs/$job.done"; return 3; fi
  fi
  RESEED_SAMPLE="$s" uv run python scripts/autoresearch/ar5/ix_exact.py ${tag:+--tag "$tag"} > "$logs/$job.log" 2>&1
  echo "exit=$?" > "$logs/$job.done"
}

live_job() { # live_job <arm old|new> <sample>
  a="$1"; s="$2"; c="$cache/live_$a"; job="live_${a}_${s}"
  [ -e "$logs/$job.done" ] && { echo "refusing: $logs/$job.done exists" >&2; return 2; }
  for p in "$c/J5D_$s" "$c/passes_$s" "$c/j5sd/J5S_$s"; do
    [ -e "$p" ] && { echo "refusing: $p exists" >&2; return 2; }
  done
  mkdir -p "$c" "$c/j5sd/passes_$s/D4_J5"
  label=contended
  if quiet_check "$ixd/preflight_$job.json"; then label=single-worker
  elif [ "$allow" = "0" ]; then echo "exit=3 refused by the preflight" > "$logs/$job.done"; return 3
  fi
  echo "$label" > "$logs/$job.label"
  env="RESEED_SAMPLE=$s DEFER_CACHE=$c SCREEN_CACHE=$c IX_INDEXER=$a"
  eval "$env $L defer_run.py --cap-min 300" > "$logs/$job.log" 2>&1; rc=$?
  echo "run=$rc" > "$logs/$job.rc"
  if [ "$rc" = "0" ] && [ -e "$c/J5D_$s/run.json" ] && [ ! -e "$c/J5D_$s/capped.json" ]; then
    eval "$env $L high500_pass.py --name J5D --aligned 1" > "$logs/$job.pass.log" 2>&1; prc=$?
    echo "pass=$prc" >> "$logs/$job.rc"
    if [ "$prc" = "0" ]; then
      ln -sfn "$c/passes_$s/before_J5D.npy" "$c/j5sd/passes_$s/before_J5.npy"
      for f in tries.npz pass.json; do ln -sfn "$c/passes_$s/D4_J5D/$f" "$c/j5sd/passes_$s/D4_J5/$f"; done
      eval "SCREEN_PRODUCER=1 RESEED_SAMPLE=$s DEFER_CACHE=$c SCREEN_CACHE=$c/j5sd IX_INDEXER=$a $L screen_run.py --arm J5S --cap-min 300" > "$logs/$job.arm.log" 2>&1; arc=$?
      echo "arm=$arc" >> "$logs/$job.rc"
    fi
  fi
  echo "exit=$rc label=$label" > "$logs/$job.done"
}

if [ "$dry" = "1" ]; then
  echo "dry run: stage $stage, cache $cache, samples: $samples, allow-contended $allow"
  if [ "$stage" = "e1" ]; then
    for s in $samples; do echo "  e1 $s: RESEED_SAMPLE=$s uv run python scripts/autoresearch/ar5/ix_exact.py ${tag:+--tag $tag} > $logs/ix_exact_${s}${tag:+_$tag}.log"; done
  else
    for s in $samples; do for a in old new; do
      c="$cache/live_$a"
      echo "  live $a $s: $L defer_run.py --cap-min 300 ; high500_pass.py --name J5D --aligned 1 ; screen_run.py --arm J5S (sub-cache $c/j5sd)"
      echo "     outputs: $c/J5D_$s $c/passes_$s $c/j5sd/J5S_$s"
    done; done
  fi
  echo "  one job at a time; preflight --require before each timed job: $pf --require"
  exit 0
fi

mkdir -p "$logs"
trap '' HUP
exec < /dev/null
exec >> "$logs/ix_run_all_$stage.out" 2>&1
echo "start $(date +%Y-%m-%dT%H:%M:%S%z) stage=$stage samples=$samples allow=$allow"
if pgrep -f "autoresearch/ar5/(ix_exact|ix_launch|ix_profile)\.py" > /dev/null 2>&1; then
  echo "refusing: an ix job is already running" >&2; exit 2
fi
fail=0
if [ "$stage" = "e1" ]; then
  for s in $samples; do e1_job "$s" || { echo "job e1 $s ended with $?"; fail=1; }; done
else
  for s in $samples; do for a in old new; do live_job "$a" "$s" || { echo "job live $a $s ended with $?"; fail=1; }; done; done
fi
echo "exit=$fail" > "$logs/ix_run_all_$stage.done"
echo "end $(date +%Y-%m-%dT%H:%M:%S%z) fail=$fail"
exit "$fail"
