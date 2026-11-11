#!/bin/bash
# AR3 step 11: the one S8 baseline on held (cap 120 minutes), then its boundary pass.
#   scripts/autoresearch/ar3/held_chain.sh        (run from icenine_py/)
# Log: benchmarks/autoresearch/ar3/held/cache/logs/chain_held.log; marker chain_held.done.
set -u
cache="$PWD/benchmarks/autoresearch/ar3/held/cache"
mkdir -p "$cache/logs"
(
  export RESEED_SAMPLE=held RESEED_CACHE="$cache"
  nohup caffeinate -i bash -c '
    uv run python scripts/autoresearch/ar5/fit_run.py --name S8 --cap-min 120 &&
    test ! -e "'"$cache"'/S8_held/capped.json" &&
    uv run python scripts/autoresearch/ar5/fit_pass.py --name S8' \
    > "$cache/logs/chain_held.log" 2>&1
  echo "exit $?" > "$cache/logs/chain_held.done"
) > /dev/null 2>&1 &
echo "started chain_held"
