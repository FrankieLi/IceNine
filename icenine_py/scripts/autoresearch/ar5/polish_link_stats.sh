#!/bin/bash
# AR5 step 14, phase A step 6 (stats): link more of step 12's stored arm-H files, per entry, read only.
#   scripts/autoresearch/ar5/polish_link_stats.sh [--dry-run]
# For big grad mid held: <cache>/H_<s>/{fit.json,e1_log.json.gz,mc_traces.json.gz,order.json,
# stallhit.json} and <cache>/passes_<s>/{D4_H/pass.json,after_H.npy,before_H.npy}.
# For held2 (after the freeze commit only): <cache>/H_held2/{fit.json,e1_log.json.gz,order.json,
# stallhit.json}. NEVER anything of a passes_held2 folder (step 13's unopened test): refused here.
# Only names are linked; no file is opened. Exit 2 if a source is missing. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
here="$PWD"
cache="${POLISH_CACHE:-$here/benchmarks/autoresearch/ar5/polish/cache}"
dry=0
[ "${1:-}" = "--dry-run" ] && dry=1
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
for s in big grad mid held; do
  for f in fit.json e1_log.json.gz mc_traces.json.gz order.json stallhit.json; do
    link "$cache/H_$s/$f" "$src/H_$s/$f"
  done
  for f in D4_H/pass.json after_H.npy before_H.npy; do
    link "$cache/passes_$s/$f" "$src/passes_$s/$f"
  done
done
for f in fit.json e1_log.json.gz order.json stallhit.json; do
  link "$cache/H_held2/$f" "$src/H_held2/$f"
done
exit $rc
