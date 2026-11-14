#!/bin/bash
# AR5 step 11: link the stored arm A and arm G runs, passes and the per-level record (not needed here, linked for tools) into the
# patience cache (read only; never overwrites; never writes into another worktree).
#   scripts/autoresearch/ar5/patience_link_cache.sh [--dry-run] [--cache DIR]
# Cache: benchmarks/autoresearch/ar5/patience/cache (untracked). Creates, only where missing:
#   <cache>/G_<sample>            -> step 9's arm G run          (big grad mid held)
#   <cache>/A_<sample>            -> step 8's arm A run          (via step 9's cache link)
#   <cache>/levels, levels_<mid|grad>.json -> step 9's per-level record (not needed here, linked for tools)
#   <cache>/passes_<sample>/{D4_A,before_A.npy,after_A.npy,D4_G,before_G.npy,after_G.npy}
#                                 -> the stored passes (a real folder: the pass of C is written here)
# Only names are linked: no file of a sample (held included) is opened. Exit 2 if a source is
# missing. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${PATIENCE_CACHE:-$here/benchmarks/autoresearch/ar5/patience/cache}"
dry=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
wt="$(cd "$here/../.." && pwd)"
src="$wt/agent-a8d20b3821621fd41/icenine_py/benchmarks/autoresearch/ar5/rescue/cache"
rc=0

link() { # link <dest> <source>
  dest="$1"; from="$2"
  if [ -e "$dest" ] || [ -L "$dest" ]; then echo "present  $dest"; return; fi
  if [ ! -e "$from" ]; then echo "MISSING SOURCE $from" >&2; rc=2; return; fi
  real="$(realpath "$from")"
  if [ "$dry" = "1" ]; then echo "would link $dest -> $real"; return; fi
  mkdir -p "$(dirname "$dest")"
  ln -s "$real" "$dest" && echo "linked   $dest"
}

for s in big grad mid held; do
  link "$cache/G_$s" "$src/G_$s"
  link "$cache/A_$s" "$src/A_$s"
  for f in D4_A before_A.npy after_A.npy D4_G before_G.npy after_G.npy; do
    link "$cache/passes_$s/$f" "$src/passes_$s/$f"
  done
done
link "$cache/levels" "$src/levels"
link "$cache/levels_mid.json" "$src/levels_mid.json"
link "$cache/levels_grad.json" "$src/levels_grad.json"
exit $rc
