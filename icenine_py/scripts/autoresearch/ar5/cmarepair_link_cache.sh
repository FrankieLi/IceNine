#!/bin/bash
# AR5 step 20: link, read only and per file, what this step reads from step 19's worktree into this
# worktree (nothing is copied, nothing is written to the sources). Step 19's own cache holds links
# (to step 16's J5 maps and pass tries, the S8 maps, step 18's J5L, the sample files and the image
# stacks); each is resolved with realpath and linked here directly:
#   * <cache>/passes_<s>/{before_J5.npy,after_J5.npy,after_S8.npy,D4_J5/{pass.json,tries.npz}}
#   * with "j5l": <cache>/passes_mid/{before_J5L.npy,after_J5L.npy,D4_J5L/{pass.json,tries.npz}}
#   * the image stack of each sample (stacks/<s>_clean.npy) and the sample files that git does not
#     hold (gradient_grid.mic, gradient_truth.mic; held3 also grainmap.npy, sample.json,
#     held_check.json)
#   scripts/autoresearch/ar5/cmarepair_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...] [j5l]
# SAMPLE in big grad mid held held3 ONLY. Anything else (held2, held4, ...) is refused: nothing of
# held2 is linked, listed or opened by this step; held4 is built by this step's own build script.
# Creates only where missing. No whole-folder link (a write through one would land in the source);
# no file is opened here. Exit 2 if a source is missing or an argument is refused. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${CMAREPAIR_CACHE:-$here/benchmarks/autoresearch/ar5/cmarepair/cache}"
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
w19="${CMAREPAIR_SRC_WT:-$wt/agent-af5e5060b23027220}"   # step 19: every stored input, as links
src19="$w19/icenine_py/benchmarks/autoresearch/ar5/repair/cache"
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
  for f in D4_J5/pass.json D4_J5/tries.npz before_J5.npy after_J5.npy after_S8.npy; do
    link "$cache/passes_$s/$f" "$src19/passes_$s/$f"
  done
  link "$ex/${s}_clean.npy" "$w19/Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/${s}_clean.npy"
  case "$s" in
    big|mid|held|held3)
      for f in gradient_grid.mic gradient_truth.mic; do
        link "$ar3/$s/gradient/sample/$f" "$w19/icenine_py/benchmarks/autoresearch/ar3/$s/gradient/sample/$f"
      done ;;
  esac
  if [ "$s" = "held3" ]; then
    for f in grainmap.npy sample.json held_check.json; do
      link "$ar3/held3/gradient/sample/$f" "$w19/icenine_py/benchmarks/autoresearch/ar3/held3/gradient/sample/$f"
    done
  fi
done
if [ "$j5l" = "1" ]; then
  for f in D4_J5L/pass.json D4_J5L/tries.npz before_J5L.npy after_J5L.npy; do
    link "$cache/passes_mid/$f" "$src19/passes_mid/$f"
  done
fi
exit $rc
