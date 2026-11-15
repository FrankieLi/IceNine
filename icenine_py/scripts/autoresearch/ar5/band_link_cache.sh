#!/bin/bash
# AR5 step 15: link step 14's stored arm-J files into step 15's own cache, per entry, read only.
#   scripts/autoresearch/ar5/band_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...]
# SAMPLE in big grad mid held ONLY. Anything else (held2, held3, ...) is refused: step 14's held2
# files and step 13's passes_held2 are never linked, listed or opened by this step; held3 is built
# by this step and has no stored run. Creates, only where missing, <cache>/J_<s>/<file> -> step 14's
# file and <cache>/passes_<s>/<file> likewise. No whole-folder link (a write through one would land
# in step 14's cache); no file is opened here. Exit 2 if a source is missing or a sample is refused.
# Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${BAND_CACHE:-$here/benchmarks/autoresearch/ar5/band/cache}"
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
rc=0
link() { # link <dest> <source>
  dest="$1"; from="$2"
  case "$from$dest" in *held2*) echo "refusing held2: $from" >&2; rc=2; return ;; esac
  if [ -e "$dest" ] || [ -L "$dest" ]; then echo "present  $dest"; return; fi
  if [ ! -e "$from" ]; then echo "MISSING SOURCE $from" >&2; rc=2; return; fi
  real="$(realpath "$from")"
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
done
exit $rc
