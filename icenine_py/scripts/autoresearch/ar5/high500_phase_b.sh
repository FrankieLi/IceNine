#!/bin/bash
# AR5 step 12, phase B: the script arm J5l and the key-on library run LibJ5 on mid (one job each).
#   cd icenine_py && nohup perl -MPOSIX -e 'POSIX::setsid(); exec @ARGV' \
#       scripts/autoresearch/ar5/high500_phase_b.sh <J5l|LibJ5> > /dev/null 2>&1 < /dev/null &
# Jobs (one heavy process each, under caffeinate -i):
#   J5l    script arm J5l (high500_hl.py: step 9's arm Gl + rule PH(500)), then its pass
#         (order_pass.py --name J5l --aligned 0, the library's single generator)
#   LibJ5  library run, keys of S11L + SeedRescueCMA 0.25 1000 + SeedIndexMCStallHit 500
#         (high500_lib_run.py; its boundary pass runs in the run)
# Resumable: a job whose <job>_mid.done marker reads exit=0 is not rerun; a job whose run folder
# exists without that marker is refused (move the folder away to rerun it).
# Logs and markers in <cache>/logs: <job>_mid.log, <job>_mid.done ("exit=<code> ..."). Caps 240
# minutes per run. Wall times are contended, not timing results. Refuses (exit 2, nothing written)
# when an output's real path lies outside the cache (stallhit_cache_guard.py). Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
job="${1:-}"
cache="${HIGH500_CACHE:-$PWD/benchmarks/autoresearch/ar5/high500/cache}"
logs="$cache/logs"
export HIGH500_CACHE="$cache" STALLHIT_CACHE="$cache" RESCUE_CACHE="$cache" ORDER_CACHE="$cache" RESEED_CACHE="$cache" RESEED_SAMPLE=mid
unset VIRTUAL_ENV
A5=scripts/autoresearch/ar5
case "$job" in
  J5l) outs="J5l_mid passes_mid passes_mid/D4_J5l passes_mid/before_J5l.npy passes_mid/after_J5l.npy" ;;
  LibJ5) outs="LibJ5_mid" ;;
  *) echo "usage: $0 J5l|LibJ5" >&2; exit 2 ;;
esac
# review fix 1 (second review): refuse before any write when the real path of an output (run
# folder, pass folder and files, logs) lies outside this cache, e.g. through a whole-folder link
# into another worktree's cache. Nothing is written on a refusal; the reason goes to stderr.
uv run python $A5/stallhit_cache_guard.py --cache "$cache" logs $outs || exit 2
mkdir -p "$logs"
if grep -q "^exit=0" "$logs/${job}_mid.done" 2>/dev/null; then
  echo "kept: $job finished before" >> "$logs/phase_b.resume"; exit 0
fi
if [ -e "$cache/${job}_mid" ]; then
  echo "refusing: $cache/${job}_mid exists without a finished marker" > "$logs/${job}_mid.refused"
  exit 2
fi
date +%Y-%m-%dT%H:%M:%S%z > "$logs/${job}_mid.start"
if [ "$job" = "J5l" ]; then
  caffeinate -i uv run python $A5/high500_hl.py --cap-min 300 > "$logs/J5l_mid.log" 2>&1
  rc=$?
  prc=skipped
  if [ "$rc" = "0" ]; then
    caffeinate -i uv run python $A5/order_pass.py --name J5l --aligned 0 \
      > "$logs/pass_J5l_mid.log" 2>&1
    prc=$?
  fi
  echo "exit=$rc pass=$prc" > "$logs/J5l_mid.done"
else
  caffeinate -i uv run python $A5/high500_lib_run.py --cap-min 300 > "$logs/LibJ5_mid.log" 2>&1
  echo "exit=$?" > "$logs/LibJ5_mid.done"
fi
