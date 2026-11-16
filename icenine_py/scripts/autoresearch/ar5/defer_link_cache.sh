#!/bin/bash
# AR5 step 22: link, read only and per file, what this step reads from earlier steps' worktrees into
# this worktree (nothing is copied, nothing is written to the sources). Each source is resolved with
# realpath and linked here directly (no whole-folder link: a write through one would land in the
# source):
#   * <cache>/J5_<s>/{e1_log.json.gz,mc_traces.json.gz,run.json,polish.json,high500.json}
#       J5's stored run: step 16 (held3 held mid grad big), step 20 (held4), step 21 (held5)
#   * <cache>/stored_passes_<s>/{D4_J5/{pass.json,tries.npz},before_J5.npy,after_J5.npy,after_S8.npy}
#       and <cache>/stored_J5S_<s>/{after.npy,repair.json}: step 21's, to compare with
#   * the image stack of each sample (stacks/<s>_clean.npy) and the sample files that are not
#     tracked (gradient_grid.mic, gradient_truth.mic; held3 to held5 also grainmap.npy, sample.json,
#     held_check.json)
#   scripts/autoresearch/ar5/defer_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...]
# SAMPLE in big grad mid held held3 held4 held5 ONLY. Anything else (held2, held6, ...) is refused:
# nothing of held2 is linked, listed or opened by this step; held6 is built by this step's own
# build script. Creates only where missing. No file is opened here. Exit 2 if a source is missing
# or an argument is refused. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${SCREEN_CACHE:-$here/benchmarks/autoresearch/ar5/defer/cache}"
dry=0
samples=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    big|grad|mid|held|held3|held4|held5) samples="$samples $1" ;;
    *) echo "refused argument (only big grad mid held held3 held4 held5): $1" >&2; exit 2 ;;
  esac
  shift
done
wt="$(cd "$here/../.." && pwd)"
w16="${DEFER_SRC_WT16:-$wt/agent-a8d3953c88c593c42}"   # step 16: J5 runs of held3 held mid grad big
w20="${DEFER_SRC_WT20:-$wt/agent-a8b2d0066621fd0f8}"   # step 20: held4's J5, the data files
w21="${DEFER_SRC_WT21:-$wt/agent-a8c1c16b70cf43e39}"   # step 21: held5, the stored passes and J5S
c16="$w16/icenine_py/benchmarks/autoresearch/ar5/high500/cache"
c20="$w20/icenine_py/benchmarks/autoresearch/ar5/cmarepair/cache"
c21="$w21/icenine_py/benchmarks/autoresearch/ar5/screen/cache"
rc=0
link() { # link <dest> <source>
  dest="$1"; from="$2"
  case "$from$dest" in *held2*|*held6*) echo "refusing held2 / held6: $from" >&2; rc=2; return ;; esac
  if [ -e "$dest" ] || [ -L "$dest" ]; then echo "present  $dest"; return; fi
  if [ ! -e "$from" ]; then echo "MISSING SOURCE $from" >&2; rc=2; return; fi
  real="$(realpath "$from")"
  case "$real" in *held2*|*held6*) echo "refusing held2 / held6 (resolved): $real" >&2; rc=2; return ;; esac
  if [ "$dry" = "1" ]; then echo "would link $dest -> $real"; return; fi
  mkdir -p "$(dirname "$dest")"
  ln -s "$real" "$dest" && echo "linked   $dest"
}
ex="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
ar3="$here/benchmarks/autoresearch/ar3"
for s in $samples; do
  case "$s" in
    held4) jsrc="$c20/J5_held4" ;;
    held5) jsrc="$c21/J5_held5" ;;
    *) jsrc="$c16/J5_$s" ;;
  esac
  for f in e1_log.json.gz mc_traces.json.gz run.json polish.json high500.json; do
    link "$cache/J5_$s/$f" "$jsrc/$f"
  done
  for f in D4_J5/pass.json D4_J5/tries.npz before_J5.npy after_J5.npy after_S8.npy; do
    link "$cache/stored_passes_$s/$f" "$c21/passes_$s/$f"
  done
  for f in after.npy repair.json; do
    link "$cache/stored_J5S_$s/$f" "$c21/J5S_$s/$f"
  done
  link "$ex/${s}_clean.npy" "$w21/Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/${s}_clean.npy"
  case "$s" in
    big|mid|held|held3|held4)
      for f in gradient_grid.mic gradient_truth.mic; do
        link "$ar3/$s/gradient/sample/$f" "$w20/icenine_py/benchmarks/autoresearch/ar3/$s/gradient/sample/$f"
      done ;;
    held5)
      for f in gradient_grid.mic gradient_truth.mic; do
        link "$ar3/$s/gradient/sample/$f" "$w21/icenine_py/benchmarks/autoresearch/ar3/$s/gradient/sample/$f"
      done ;;
  esac
  case "$s" in
    held3|held4)
      for f in grainmap.npy sample.json held_check.json; do
        link "$ar3/$s/gradient/sample/$f" "$w20/icenine_py/benchmarks/autoresearch/ar3/$s/gradient/sample/$f"
      done ;;
    held5)
      for f in grainmap.npy sample.json held_check.json; do
        link "$ar3/$s/gradient/sample/$f" "$w21/icenine_py/benchmarks/autoresearch/ar3/$s/gradient/sample/$f"
      done ;;
  esac
done
exit $rc
