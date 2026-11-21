#!/bin/bash
# AR5 step 19: link, read only and per file, what this step reads from other worktrees into this
# worktree (nothing is copied, nothing is written to the sources):
#   * step 16's stored arm-J5 pass files  -> <cache>/passes_<s>/ (before_J5, after_J5, D4_J5/{tries.npz,pass.json})
#   * the stored S8 end map after its pass -> <cache>/passes_<s>/after_S8.npy (held3: step 16's cache;
#     held, mid, grad, big: where step 16's stats read them, fit_common.STORED_PASS / held's AR3 cache)
#   * step 18's live J5L run of mid        -> <cache>/passes_mid/ (before_J5L, after_J5L, D4_J5L/...) with "j5l"
#   * the sample files and image stacks a live repair reads (untracked in git)
#   scripts/autoresearch/ar5/repair_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...] [j5l]
# SAMPLE in big grad mid held held3 ONLY. Anything else (held2, held4, ...) is refused: nothing of
# held2 is linked, listed or opened by this step; held4 is built by this step's own build script.
# Creates only where missing. No whole-folder link (a write through one would land in the source);
# no file is opened here. Exit 2 if a source is missing or an argument is refused. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${REPAIR_CACHE:-$here/benchmarks/autoresearch/ar5/repair/cache}"
dry=0
samples=""
j5l=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    j5l) j5l=1 ;;
    big|grad|mid|held|held3) samples="$samples $1" ;;
    *) echo "refused argument (only big grad mid held held3 j5l): $1" >&2; exit 2 ;;
  esac
  shift
done
wt="$(cd "$here/../.." && pwd)"
w16="${HIGH500_SRC_WT:-$wt/agent-a8d3953c88c593c42}"   # step 16: J5 passes, S8 on held3, held3's inputs
w18="${LEARNSKIP_SRC_WT:-$wt/agent-adedc789a656fa5ed}" # step 18: J5L (mid), the sample links
src16="$w16/icenine_py/benchmarks/autoresearch/ar5/high500/cache"
src18="$w18/icenine_py/benchmarks/autoresearch/ar5/learnskip/cache"
s8_reseed="${RESEED_SRC:-$wt/agent-a635b2f904279bdc9/icenine_py/benchmarks/autoresearch/ar3/reseed/cache}"
s8_start="${START_SRC:-$wt/agent-a9c88b610fb10e1b6/icenine_py/benchmarks/autoresearch/ar3/start/cache}"
s8_held="${HELD_SRC:-$wt/agent-af8961874e131bda2/icenine_py/benchmarks/autoresearch/ar3/held/cache}"
rc=0
link() { # link <dest> <source>
  dest="$1"; from="$2"
  case "$from$dest" in *held2*|*held4*) echo "refusing held2 / held4: $from" >&2; rc=2; return ;; esac
  if [ -e "$dest" ] || [ -L "$dest" ]; then echo "present  $dest"; return; fi
  if [ ! -e "$from" ]; then echo "MISSING SOURCE $from" >&2; rc=2; return; fi
  real="$(realpath "$from")"
  case "$real" in *held2*|*held4*) echo "refusing held2 / held4 (resolved): $real" >&2; rc=2; return ;; esac
  if [ "$dry" = "1" ]; then echo "would link $dest -> $real"; return; fi
  mkdir -p "$(dirname "$dest")"
  ln -s "$real" "$dest" && echo "linked   $dest"
}
ex="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
ar3="$here/benchmarks/autoresearch/ar3"
for s in $samples; do
  for f in D4_J5/pass.json D4_J5/tries.npz before_J5.npy after_J5.npy; do
    link "$cache/passes_$s/$f" "$src16/passes_$s/$f"
  done
  case "$s" in
    held3) s8="$src16/passes_held3" ;;
    held) s8="$s8_held/passes_held" ;;
    mid) s8="$s8_start/passes_mid" ;;
    grad|big) s8="$s8_reseed/passes_$s" ;;
  esac
  link "$cache/passes_$s/after_S8.npy" "$s8/after_S8.npy"
  # the image stack of the sample (a live repair reads it); grad's lives in the main checkout
  link "$ex/${s}_clean.npy" "$w18/Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/${s}_clean.npy"
  # sample files that git does not hold (the grid and truth maps of big, mid, held, held3)
  case "$s" in
    big|mid|held|held3)
      for f in gradient_grid.mic gradient_truth.mic; do
        link "$ar3/$s/gradient/sample/$f" "$w18/icenine_py/benchmarks/autoresearch/ar3/$s/gradient/sample/$f"
      done ;;
  esac
  if [ "$s" = "held3" ]; then
    for f in grainmap.npy sample.json held_check.json; do
      link "$ar3/held3/gradient/sample/$f" "$w18/icenine_py/benchmarks/autoresearch/ar3/held3/gradient/sample/$f"
    done
  fi
done
if [ "$j5l" = "1" ]; then
  for f in D4_J5L/pass.json D4_J5L/tries.npz before_J5L.npy after_J5L.npy; do
    link "$cache/passes_mid/$f" "$src18/passes_mid/$f"
  done
fi
exit $rc
