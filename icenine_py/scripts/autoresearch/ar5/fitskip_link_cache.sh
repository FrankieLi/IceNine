#!/bin/bash
# AR5 step 17: link step 16's stored arm-J5 files into step 17's own cache, per entry, read only;
# and link the input files of sample held3 (built by step 16) that a fresh worktree lacks.
#   scripts/autoresearch/ar5/fitskip_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...]
# SAMPLE in big grad mid held held3 ONLY. Anything else (held2, held4, ...) is refused: nothing of
# held2 is linked, listed or opened by this step; held4 does not exist yet and is built by a later
# step. Creates, only where missing, <cache>/J5_<s>/<file> and <cache>/passes_<s>/<file>: a few
# files per sample (run.json, recon.mic.gz, polish.json, high500.json, the pass's pass.json and its
# before/after maps). No whole-folder link (a write through one would land in step 16's cache); no
# file is opened here; the logs, records and e1 logs of step 16 are not linked. Exit 2 if a source
# is missing or a sample is refused. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${FITSKIP_CACHE:-$here/benchmarks/autoresearch/ar5/fitskip/cache}"
dry=0
samples=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    big|grad|mid|held|held3) samples="$samples $1" ;;
    *) echo "refused argument (only big grad mid held held3): $1" >&2; exit 2 ;;
  esac
  shift
done
wt="$(cd "$here/../.." && pwd)"
sw="${HIGH500_SRC_WT:-$wt/agent-a8d3953c88c593c42}"
src="${HIGH500_SRC:-$sw/icenine_py/benchmarks/autoresearch/ar5/high500/cache}"
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
for s in $samples; do
  for f in run.json recon.mic.gz polish.json high500.json; do
    link "$cache/J5_$s/$f" "$src/J5_$s/$f"
  done
  for f in D4_J5/pass.json before_J5.npy after_J5.npy; do
    link "$cache/passes_$s/$f" "$src/passes_$s/$f"
  done
done
# held3's inputs (made by step 16, untracked): linked per file from step 16's worktree
exs="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
ar3="$here/benchmarks/autoresearch/ar3"
for s in $samples; do
  if [ "$s" = "held3" ]; then
    link "$exs/held3_clean.npy" "$sw/Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/held3_clean.npy"
    for f in gradient_grid.mic gradient_truth.mic grainmap.npy sample.json held_check.json; do
      link "$ar3/held3/gradient/sample/$f" "$sw/icenine_py/benchmarks/autoresearch/ar3/held3/gradient/sample/$f"
    done
  fi
done
exit $rc
