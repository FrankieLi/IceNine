#!/bin/bash
# AR5 step 24: link, read only and per file, what this step reads from step 22's worktree into this
# worktree and into this step's cache (nothing is copied, nothing is written to the sources). Each
# source is resolved with realpath and linked directly (no whole-folder link: a write through one
# would land in the source):
#   * <cache>/J5D_<s>/{e1_log.json.gz,mc_traces.json.gz,run.json,polish.json,defer.json,
#       defer_steps.json.gz,high500.json}: step 22's live J5D run (the replay's input)
#   * <cache>/stored_J5SD_<s>/{after.npy,repair.json}: step 22's J5SD arm (to compare with)
#   * <cache>/stored_S8_<s>/after.npy: the S8 map (step 21's stored pass; held6: step 22's)
#   * the image stack of each sample (stacks/<s>_clean.npy) and the sample files that are not
#     tracked (gradient_grid.mic, gradient_truth.mic; also grainmap.npy, sample.json and the check
#     file where the worktree lacks them)
#   scripts/autoresearch/ar5/bf_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...]
# SAMPLE in big grad mid held held3 held4 held5 held6 ONLY. Anything else (held2, held7, fine7 ...)
# is refused: nothing of held2 is linked, listed or opened by this step; held7 is built by this
# step's own build script. Creates only where missing. No file is opened here. Exit 2 if a source
# is missing or an argument is refused. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${BF_CACHE:-$here/benchmarks/autoresearch/ar5/bf/cache}"
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
w22="${BF_SRC_WT22:-$wt/agent-aba6e1bdb95dff72a}"   # step 22: J5D runs, J5SD arms, held6, data
c22="$w22/icenine_py/benchmarks/autoresearch/ar5/defer/cache"
main="$(cd "$wt/../.." && pwd)"   # the main checkout (worktrees live in <main>/.claude/worktrees)
cb23="${BF_SRC_BALL:-$main/icenine_py/benchmarks/autoresearch/ar5/ball/cache}"  # step 23's replay rows
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
    link "$cache/J5D_$s/$f" "$c22/J5D_$s/$f"
  done
  for f in rows.jsonl summary.json replay.done; do
    link "$cache/replay_on_$s/$f" "$cb23/replay_on_$s/$f"
  done
  for f in after.npy repair.json; do
    link "$cache/stored_J5SD_$s/$f" "$c22/j5sd/J5S_$s/$f"
  done
  if [ "$s" = "held6" ]; then
    link "$cache/stored_S8_$s/after.npy" "$c22/passes_$s/after_S8.npy"
  else
    link "$cache/stored_S8_$s/after.npy" "$c22/stored_passes_$s/after_S8.npy"
  fi
  link "$ex/${s}_clean.npy" "$w22/Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/${s}_clean.npy"
  if [ "$s" = "grad" ]; then d="gradient/sample"; else d="$s/gradient/sample"; fi
  for f in gradient_grid.mic gradient_truth.mic grainmap.npy sample.json held_check.json mid_check.json big_check.json; do
    if [ -e "$src3/$d/$f" ]; then link "$ar3/$d/$f" "$src3/$d/$f"; fi
  done
done
exit $rc
