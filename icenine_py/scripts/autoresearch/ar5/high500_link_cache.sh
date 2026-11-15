#!/bin/bash
# AR5 step 16: link step 14's stored arm-J files and step 15's recording rows into step 16's own
# cache, per entry, read only; and link the held, mid, grad and big input files (stacks, truth and
# grid maps) that a fresh worktree lacks.
#   scripts/autoresearch/ar5/high500_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...]
# SAMPLE in big grad mid held ONLY. Anything else (held2, held3, ...) is refused: nothing of held2
# is linked, listed or opened by this step; held3 is built by this step and has no stored run.
# Creates, only where missing, <cache>/J_<s>/<file>, <cache>/passes_<s>/<file> and
# <cache>/record/<row file>: no whole-folder link (a write through one would land in another step's
# cache); no file is opened here. Exit 2 if a source is missing or a sample is refused.
# Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${HIGH500_CACHE:-$here/benchmarks/autoresearch/ar5/high500/cache}"
dry=0
samples=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    big|grad|mid|held) samples="$samples $1" ;;
    *) echo "refused argument (only big grad mid held): $1" >&2; exit 2 ;;
  esac
  shift
done
wt="$(cd "$here/../.." && pwd)"
src="${POLISH_SRC:-$wt/agent-a1bd2e7a7dbfe1efd/icenine_py/benchmarks/autoresearch/ar5/polish/cache}"
rsrc="${BAND_SRC:-$wt/agent-a38701500191fde36/icenine_py/benchmarks/autoresearch/ar5/band/cache}"
rc=0
link() { # link <dest> <source>
  dest="$1"; from="$2"
  case "$from$dest" in *held2*) echo "refusing held2: $from" >&2; rc=2; return ;; esac
  if [ -e "$dest" ] || [ -L "$dest" ]; then echo "present  $dest"; return; fi
  if [ ! -e "$from" ]; then echo "MISSING SOURCE $from" >&2; rc=2; return; fi
  real="$(realpath "$from")"
  case "$real" in *held2*) echo "refusing held2 (resolved): $real" >&2; rc=2; return ;; esac
  if [ "$dry" = "1" ]; then echo "would link $dest -> $real"; return; fi
  mkdir -p "$(dirname "$dest")"
  ln -s "$real" "$dest" && echo "linked   $dest"
}
for s in $samples; do
  for f in e1_log.json.gz mc_traces.json.gz run.json records.json.gz recon.mic.gz \
           polish_log.json.gz polish.json stallhit.json order.json; do
    link "$cache/J_$s/$f" "$src/J_$s/$f"
  done
  for f in D4_J/pass.json before_J.npy after_J.npy; do
    link "$cache/passes_$s/$f" "$src/passes_$s/$f"
  done
  case "$s" in
    held) rows="record_held_s0of2.jsonl record_held_s1of2.jsonl" ;;
    *) rows="record_${s}_s0of1.jsonl" ;;
  esac
  for f in $rows; do
    link "$cache/record/$f" "$rsrc/record/$f"
  done
done
# input files a fresh worktree lacks (gitignored or untracked): linked from step 14's worktree
# (which points at the real files), per file
exs="$here/../Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
exsrc="$wt/agent-a1bd2e7a7dbfe1efd/Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks"
ar3="$here/benchmarks/autoresearch/ar3"
ar3src="$wt/agent-a1bd2e7a7dbfe1efd/icenine_py/benchmarks/autoresearch/ar3"
for s in $samples; do
  link "$exs/${s}_clean.npy" "$exsrc/${s}_clean.npy"
  case "$s" in
    held|mid|big)
      link "$ar3/$s/gradient/sample/gradient_truth.mic" "$ar3src/$s/gradient/sample/gradient_truth.mic"
      link "$ar3/$s/gradient/sample/gradient_grid.mic" "$ar3src/$s/gradient/sample/gradient_grid.mic"
      ;;
  esac
done
exit $rc
