#!/bin/bash
# AR5 step 14: link step 12's stored arm-H files into step 14's own cache, per entry, read only.
#   scripts/autoresearch/ar5/polish_link_cache.sh [--dry-run] [--cache DIR] SAMPLE [SAMPLE ...]
# SAMPLE in big grad mid held (gate 1) and held2 (only after ar5/polish/frozen.json is committed:
# ``polish_freeze.py --check`` must pass). Creates, only where missing, <cache>/H_<s>/<file> ->
# step 12's file, for the files below. No whole-folder link (a write through one would land in step
# 12's cache), no file of any passes_held2 folder, no file is opened here. Exit 2 if a source is
# missing or the sample is held2 without a passing freeze check. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${POLISH_CACHE:-$here/benchmarks/autoresearch/ar5/polish/cache}"
dry=0
samples=""
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry=1 ;;
    --cache) shift; cache="$1" ;;
    big|grad|mid|held|held2) samples="$samples $1" ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
wt="$(cd "$here/../.." && pwd)"
src="${STALLHIT_SRC:-$wt/agent-a0ffe667a1a6ee068/icenine_py/benchmarks/autoresearch/ar5/stallhit/cache}"
rc=0
link() { # link <dest> <source>
  dest="$1"; from="$2"
  case "$from$dest" in *passes_held2*) echo "refusing passes_held2: $from" >&2; rc=2; return ;; esac
  if [ -e "$dest" ] || [ -L "$dest" ]; then echo "present  $dest"; return; fi
  if [ ! -e "$from" ]; then echo "MISSING SOURCE $from" >&2; rc=2; return; fi
  real="$(realpath "$from")"
  if [ "$dry" = "1" ]; then echo "would link $dest -> $real"; return; fi
  mkdir -p "$(dirname "$dest")"
  ln -s "$real" "$dest" && echo "linked   $dest"
}
for s in $samples; do
  if [ "$s" = "held2" ]; then
    if ! uv run python scripts/autoresearch/ar5/polish_freeze.py --check >/dev/null 2>&1; then
      echo "held2 is not linked before the freeze commit (polish_freeze.py --check fails)" >&2
      exit 2
    fi
  fi
  for f in rescue_log.json.gz records.json.gz recon.mic.gz run.json rescue.json; do
    link "$cache/H_$s/$f" "$src/H_$s/$f"
  done
done
exit $rc
