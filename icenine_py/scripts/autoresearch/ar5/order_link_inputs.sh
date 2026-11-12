#!/bin/bash
# AR5 step 8: link the untracked inputs of the sample runs from the kept worktrees (read only).
#   scripts/autoresearch/ar5/order_link_inputs.sh [--dry-run]
# Creates, only where missing, symlinks (never overwrites, never writes into another worktree):
#   benchmarks/autoresearch/ar3/<big|mid|held>/gradient/sample/gradient_{grid,truth}.mic
#   ../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/<big|grad|mid|held>_clean.npy
# The targets are the real files behind the links in the kept worktrees of steps 3 and 4. Prints
# one line per link; exit 2 if a source is missing. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
wt="$(cd "$here/../.." && pwd)"
dry=0
[ "${1:-}" = "--dry-run" ] && dry=1
mix="$wt/agent-a5a3601b746680640"
held="$wt/agent-a8a1ac8e80f8913cf"
rc=0

link() { # link <dest> <source-link-or-file>
  dest="$1"; src="$2"
  if [ -e "$dest" ]; then echo "present  $dest"; return; fi
  if [ ! -e "$src" ]; then echo "MISSING SOURCE $src" >&2; rc=2; return; fi
  real="$(realpath "$src")"
  if [ "$dry" = "1" ]; then echo "would link $dest -> $real"; return; fi
  mkdir -p "$(dirname "$dest")"
  ln -s "$real" "$dest" && echo "linked   $dest"
}

for s in big mid; do
  for f in gradient_grid.mic gradient_truth.mic; do
    link "$here/benchmarks/autoresearch/ar3/$s/gradient/sample/$f" \
         "$mix/icenine_py/benchmarks/autoresearch/ar3/$s/gradient/sample/$f"
  done
done
for f in gradient_grid.mic gradient_truth.mic; do
  link "$here/benchmarks/autoresearch/ar3/held/gradient/sample/$f" \
       "$held/icenine_py/benchmarks/autoresearch/ar3/held/gradient/sample/$f"
done
st="Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
for s in big grad mid; do
  link "$here/../$st/${s}_clean.npy" "$mix/$st/${s}_clean.npy"
done
link "$here/../$st/held_clean.npy" "$held/$st/held_clean.npy"
exit $rc
