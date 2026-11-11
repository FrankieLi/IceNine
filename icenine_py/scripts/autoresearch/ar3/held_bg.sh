#!/bin/bash
# AR3 step 11: start one study script in the background, with a log and a .done marker.
#   scripts/autoresearch/ar3/held_bg.sh <sample: big|held> <name> <dir/script.py> [args...]
# Run from icenine_py/. The cache folder is benchmarks/autoresearch/ar3/held/cache (untracked); the
# log is cache/logs/<name>.log and the marker cache/logs/<name>.done (written when the script ends).
set -u
sample="$1"; name="$2"; script="$3"; shift 3
cache="$PWD/benchmarks/autoresearch/ar3/held/cache"
mkdir -p "$cache/logs"
(
  export RESEED_SAMPLE="$sample" RESEED_CACHE="$cache"
  nohup caffeinate -i uv run python "scripts/autoresearch/$script" "$@" \
    > "$cache/logs/$name.log" 2>&1
  echo "exit $?" > "$cache/logs/$name.done"
) > /dev/null 2>&1 &
echo "started $name"
