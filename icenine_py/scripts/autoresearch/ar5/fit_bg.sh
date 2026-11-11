#!/bin/bash
# AR5 step 2: start one study script in the background, with a log and a .done marker.
#   scripts/autoresearch/ar5/fit_bg.sh <sample: big|grad|mid> <name> <script.py> [args...]
# Run from icenine_py/. The cache is benchmarks/autoresearch/ar5/fit/cache (untracked); the log is
# cache/logs/<name>.log and the marker cache/logs/<name>.done (written when the script ends).
# Needs RESEED_SAMPLE for the AR3 modules; the order rule for mid is enforced in fit_common.py.
set -u
sample="$1"; name="$2"; script="$3"; shift 3
cache="$PWD/benchmarks/autoresearch/ar5/fit/cache"
mkdir -p "$cache/logs"
(
  export RESEED_SAMPLE="$sample" RESEED_CACHE="$cache"
  nohup caffeinate -i uv run python "scripts/autoresearch/ar5/$script" "$@" \
    > "$cache/logs/$name.log" 2>&1
  echo "exit $?" > "$cache/logs/$name.done"
) > /dev/null 2>&1 &
echo "started $name"
