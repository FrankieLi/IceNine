#!/bin/bash
# AR5 step 27: link, read only and per file, what this step reads into this worktree and into this
# step's cache (nothing is copied, nothing is written to a source). Each source is resolved with
# realpath and linked directly (no whole-folder link):
#   * <cache>/J5D_<s>/{e1_log.json.gz,mc_traces.json.gz,run.json,polish.json,defer.json,
#       defer_steps.json.gz,high500.json}, <cache>/stored_J5SD_<s>/{after.npy,repair.json},
#       <cache>/stored_S8_<s>/after.npy: from the main checkout's step 26 cache (ar5/ix/cache),
#       which links step 22's stored runs
#   * the image stack of each sample (stacks/<s>_clean.npy) and the untracked sample files
#     (gradient_grid.mic, gradient_truth.mic, grainmap.npy, sample.json, check files), from step
#     22's worktree (PX_SRC_WT22) or the main checkout
#   scripts/autoresearch/ar5/px_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...]
# SAMPLE in big grad mid held held3 held4 held5 held6 ONLY. Anything else (held2, held7, fine7 ...)
# is refused. Creates only where missing. No file is opened here. Exit 2 if a source is missing or
# an argument is refused. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${PX_CACHE:-$here/benchmarks/autoresearch/ar5/px/cache}"
dry=0
samples=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    big|grad|mid|held|held3|held4|held5|held6) samples="$samples $1" ;;
    *) echo "refused argument (only big grad mid held held3 held4 held5 held6): $1" >&2; exit 2 ;;
  esac
  shift
done
wt="$(cd "$here/../.." && pwd)"
main="$(cd "$wt/../.." && pwd)"   # the main checkout (worktrees live in <main>/.claude/worktrees)
w22="${PX_SRC_WT22:-$wt/agent-aba6e1bdb95dff72a}"
cix="${PX_SRC_IX:-$main/icenine_py/benchmarks/autoresearch/ar5/ix/cache}"
rc=0
link() { # link <dest> <source>
  dest="$1"; from="$2"
  case "$from$dest" in *held2*|*held7*|*fine7*) echo "refusing held2 / held7 / fine7: $from" >&2; rc=2; return ;; esac
  if [ -e "$dest" ] || [ -L "$dest" ]; then echo "present  $dest"; return; fi
  if [ ! -e "$from" ]; then echo "MISSING SOURCE $from" >&2; rc=2; return; fi
  real="$(realpath "$from")"
  case "$real" in *held2*|*held7*|*fine7*) echo "refusing (resolved): $real" >&2; rc=2; return ;; esac
  if [ "$dry" = "1" ]; then echo "would link $dest -> $real"; return; fi
  mkdir -p "$(dirname "$dest")"
  ln -s "$real" "$dest" && echo "linked   $dest"
}
ex="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
ar3="$here/benchmarks/autoresearch/ar3"
src3="$w22/icenine_py/benchmarks/autoresearch/ar3"
for s in $samples; do
  for f in e1_log.json.gz mc_traces.json.gz run.json polish.json defer.json defer_steps.json.gz high500.json; do
    link "$cache/J5D_$s/$f" "$cix/J5D_$s/$f"
  done
  for f in after.npy repair.json; do
    link "$cache/stored_J5SD_$s/$f" "$cix/stored_J5SD_$s/$f"
  done
  link "$cache/stored_S8_$s/after.npy" "$cix/stored_S8_$s/after.npy"
  link "$ex/${s}_clean.npy" "$w22/Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/${s}_clean.npy"
  if [ "$s" = "grad" ]; then d="gradient/sample"; else d="$s/gradient/sample"; fi
  for f in gradient_grid.mic gradient_truth.mic grainmap.npy sample.json held_check.json mid_check.json big_check.json; do
    if [ -e "$src3/$d/$f" ]; then link "$ar3/$d/$f" "$src3/$d/$f"; fi
  done
done
exit $rc
