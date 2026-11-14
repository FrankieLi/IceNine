#!/bin/bash
# AR5 step 12, phase B: the script arm Hl and the key-on library run LibH on mid (one job each).
#   cd icenine_py && nohup perl -MPOSIX -e 'POSIX::setsid(); exec @ARGV' \
#       scripts/autoresearch/ar5/stallhit_phase_b.sh <Hl|LibH> > /dev/null 2>&1 < /dev/null &
# Jobs (one heavy process each, under caffeinate -i):
#   Hl    script arm Hl (stallhit_hl.py: step 9's arm Gl + rule PH(3000)), then its pass
#         (order_pass.py --name Hl --aligned 0, the library's single generator)
#   LibH  library run, keys of S11L + SeedRescueCMA 0.25 1000 + SeedIndexMCStallHit 3000
#         (stallhit_lib_run.py; its boundary pass runs in the run)
# Resumable: a job whose <job>_mid.done marker reads exit=0 is not rerun; a job whose run folder
# exists without that marker is refused (move the folder away to rerun it).
# Logs and markers in <cache>/logs: <job>_mid.log, <job>_mid.done ("exit=<code> ..."). Caps 240
# minutes per run. Wall times are contended, not timing results. Refuses (exit 2, nothing written)
# when an output's real path lies outside the cache (stallhit_cache_guard.py). Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
job="${1:-}"
cache="${STALLHIT_CACHE:-$PWD/benchmarks/autoresearch/ar5/stallhit/cache}"
logs="$cache/logs"
export STALLHIT_CACHE="$cache" RESCUE_CACHE="$cache" ORDER_CACHE="$cache" RESEED_SAMPLE=mid
unset VIRTUAL_ENV
A5=scripts/autoresearch/ar5
case "$job" in
  Hl) outs="Hl_mid passes_mid passes_mid/D4_Hl passes_mid/before_Hl.npy passes_mid/after_Hl.npy" ;;
  LibH) outs="LibH_mid" ;;
  *) echo "usage: $0 Hl|LibH" >&2; exit 2 ;;
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
if [ "$job" = "Hl" ]; then
  caffeinate -i uv run python $A5/stallhit_hl.py --cap-min 240 > "$logs/Hl_mid.log" 2>&1
  rc=$?
  prc=skipped
  if [ "$rc" = "0" ]; then
    caffeinate -i uv run python $A5/order_pass.py --name Hl --aligned 0 \
      > "$logs/pass_Hl_mid.log" 2>&1
    prc=$?
  fi
  echo "exit=$rc pass=$prc" > "$logs/Hl_mid.done"
else
  caffeinate -i uv run python $A5/stallhit_lib_run.py --cap-min 240 > "$logs/LibH_mid.log" 2>&1
  echo "exit=$?" > "$logs/LibH_mid.done"
fi
