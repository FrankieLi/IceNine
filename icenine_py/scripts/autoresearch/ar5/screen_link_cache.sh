#!/bin/bash
# AR5 step 21: link, read only and per file, what this step reads from step 20's worktree into this
# worktree (nothing is copied, nothing is written to the sources). Step 20's own cache holds links
# (to step 16's J5 maps and pass tries, the S8 maps) and, for held4, its own files; each is
# resolved with realpath and linked here directly:
#   * <cache>/passes_<s>/{before_J5.npy,after_J5.npy,after_S8.npy,D4_J5/{pass.json,tries.npz}}
#   * <cache>/stored_J5C_<s>/{after.npy,repair.json}   (step 20's J5C, to compare with)
#   * the image stack of each sample (stacks/<s>_clean.npy) and the sample files that are not
#     tracked (gradient_grid.mic, gradient_truth.mic; held3 and held4 also grainmap.npy,
#     sample.json, held_check.json)
#   scripts/autoresearch/ar5/screen_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...]
# SAMPLE in big grad mid held held3 held4 ONLY. Anything else (held2, held5, ...) is refused:
# nothing of held2 is linked, listed or opened by this step; held5 is built by this step's own
# build script. Creates only where missing. No whole-folder link (a write through one would land in
# the source); no file is opened here. Exit 2 if a source is missing or an argument is refused.
# Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${SCREEN_CACHE:-$here/benchmarks/autoresearch/ar5/screen/cache}"
dry=0
samples=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    big|grad|mid|held|held3|held4) samples="$samples $1" ;;
    *) echo "refused argument (only big grad mid held held3 held4): $1" >&2; exit 2 ;;
  esac
  shift
done
wt="$(cd "$here/../.." && pwd)"
w20="${SCREEN_SRC_WT:-$wt/agent-a8b2d0066621fd0f8}"   # step 20: every stored input, as links
src20="$w20/icenine_py/benchmarks/autoresearch/ar5/cmarepair/cache"
rc=0
link() { # link <dest> <source>
  dest="$1"; from="$2"
  case "$from$dest" in *held2*|*held5*) echo "refusing held2 / held5: $from" >&2; rc=2; return ;; esac
  if [ -e "$dest" ] || [ -L "$dest" ]; then echo "present  $dest"; return; fi
  if [ ! -e "$from" ]; then echo "MISSING SOURCE $from" >&2; rc=2; return; fi
  real="$(realpath "$from")"
  case "$real" in *held2*|*held5*) echo "refusing held2 / held5 (resolved): $real" >&2; rc=2; return ;; esac
  if [ "$dry" = "1" ]; then echo "would link $dest -> $real"; return; fi
  mkdir -p "$(dirname "$dest")"
  ln -s "$real" "$dest" && echo "linked   $dest"
}
ex="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
ar3="$here/benchmarks/autoresearch/ar3"
for s in $samples; do
  for f in D4_J5/pass.json D4_J5/tries.npz before_J5.npy after_J5.npy after_S8.npy; do
    link "$cache/passes_$s/$f" "$src20/passes_$s/$f"
  done
  for f in after.npy repair.json; do
    link "$cache/stored_J5C_$s/$f" "$src20/J5C_$s/$f"
  done
  link "$ex/${s}_clean.npy" "$w20/Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/${s}_clean.npy"
  case "$s" in
    big|mid|held|held3|held4)
      for f in gradient_grid.mic gradient_truth.mic; do
        link "$ar3/$s/gradient/sample/$f" "$w20/icenine_py/benchmarks/autoresearch/ar3/$s/gradient/sample/$f"
      done ;;
  esac
  case "$s" in
    held3|held4)
      for f in grainmap.npy sample.json held_check.json; do
        link "$ar3/$s/gradient/sample/$f" "$w20/icenine_py/benchmarks/autoresearch/ar3/$s/gradient/sample/$f"
      done ;;
  esac
done
exit $rc
