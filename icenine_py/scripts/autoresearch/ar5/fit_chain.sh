#!/bin/bash
# AR5 step 2, Part 2: one run, its boundary pass and its numbers, one after the other.
#   scripts/autoresearch/ar5/fit_chain.sh <sample: big|grad|mid> <name> [fit_run.py options...]
# Run from icenine_py/. Log: cache/logs/chain_<name>_<sample>.log, marker: chain_<name>_<sample>.done.
# If the run does not finish (the 90 minute cap), the pass is skipped and the stats script writes
# "unfinished".
set -u
sample="$1"; name="$2"; shift 2
cache="$PWD/benchmarks/autoresearch/ar5/fit/cache"
mkdir -p "$cache/logs"
(
  export RESEED_SAMPLE="$sample" RESEED_CACHE="$cache"
  run="$cache/${name}_${sample}"
  nohup caffeinate -i bash -c "
    uv run python scripts/autoresearch/ar5/fit_run.py --name $name $* &&
    uv run python scripts/autoresearch/ar5/fit_pass.py --name $name;
    uv run python scripts/autoresearch/ar5/fit_part2_stats.py --name $name
  " > "$cache/logs/chain_${name}_${sample}.log" 2>&1
  echo "exit $?" > "$cache/logs/chain_${name}_${sample}.done"
) > /dev/null 2>&1 &
echo "started chain $name $sample"
