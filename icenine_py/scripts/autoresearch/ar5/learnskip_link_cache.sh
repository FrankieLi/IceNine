#!/bin/bash
# AR5 step 18: link, read only and per file, what this step reads from other worktrees into this
# worktree (nothing is copied, nothing is written to the sources):
#   * step 16's stored arm-J5 run and pass files  -> <cache>/J5_<s>/, <cache>/passes_<s>/
#   * step 16's S8 on held3                       -> <cache>/S8_held3/, <cache>/passes_held3/D4_S8
#   * step 17's recording J5r (rows, final flags, the fit-count checks p1_*) -> <cache>/record/, <cache>/
#   * the sample files and image stacks a run reads (untracked in git; as step 17's worktree has them)
#   scripts/autoresearch/ar5/learnskip_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...]
# SAMPLE in big grad mid held held3 ONLY. Anything else (held2, held4, ...) is refused: nothing of
# held2 is linked, listed or opened by this step; held4 is built by this step's own build script.
# Creates only where missing. No whole-folder link (a write through one would land in the source);
# no file is opened here. Exit 2 if a source is missing or an argument is refused. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${LEARNSKIP_CACHE:-$here/benchmarks/autoresearch/ar5/learnskip/cache}"
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
w16="${HIGH500_SRC_WT:-$wt/agent-a8d3953c88c593c42}"   # step 16: J5, S8 on held3, held3's inputs
w17="${FITSKIP_SRC_WT:-$wt/agent-aeba15b7c3850b007}"   # step 17: the recording J5r, sample links
src16="$w16/icenine_py/benchmarks/autoresearch/ar5/high500/cache"
src17="$w17/icenine_py/benchmarks/autoresearch/ar5/fitskip/cache"
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
  for f in run.json recon.mic.gz polish.json high500.json; do
    link "$cache/J5_$s/$f" "$src16/J5_$s/$f"
  done
  for f in D4_J5/pass.json before_J5.npy after_J5.npy; do
    link "$cache/passes_$s/$f" "$src16/passes_$s/$f"
  done
  # the recording J5r: rows, final-answer flags and the fit-count check of step 17
  link "$cache/record/J5r_$s.jsonl.gz" "$src17/record/J5r_$s.jsonl.gz"
  link "$cache/record/final_J5r_$s.json" "$src17/record/final_J5r_$s.json"
  link "$cache/p1_J5r_$s.json" "$src17/p1_J5r_$s.json"
  # the image stack of the sample (a run reads it)
  link "$ex/${s}_clean.npy" "$w17/Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/${s}_clean.npy"
  # sample files that git does not hold (the grid and truth maps of big, mid, held, held3)
  case "$s" in
    big|mid|held|held3)
      for f in gradient_grid.mic gradient_truth.mic; do
        link "$ar3/$s/gradient/sample/$f" "$w17/icenine_py/benchmarks/autoresearch/ar3/$s/gradient/sample/$f"
      done ;;
  esac
  if [ "$s" = "held3" ]; then
    for f in grainmap.npy sample.json held_check.json; do
      link "$ar3/held3/gradient/sample/$f" "$w17/icenine_py/benchmarks/autoresearch/ar3/held3/gradient/sample/$f"
    done
  fi
done
case "$samples" in
  *held3*)
    for f in run.json recon.mic.gz e1_log.json.gz; do
      link "$cache/S8_held3/$f" "$src16/S8_held3/$f"
    done
    link "$cache/passes_held3/D4_S8/pass.json" "$src16/passes_held3/D4_S8/pass.json"
    link "$cache/passes_held3/before_S8.npy" "$src16/passes_held3/before_S8.npy"
    link "$cache/passes_held3/after_S8.npy" "$src16/passes_held3/after_S8.npy" ;;
esac
exit $rc
