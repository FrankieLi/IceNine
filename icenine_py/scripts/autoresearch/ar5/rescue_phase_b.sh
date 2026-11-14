#!/bin/bash
# AR5 step 9, phase B: the key-on library run on mid and the script arm it must equal, side by side.
#   cd icenine_py && nohup perl -MPOSIX -e 'POSIX::setsid(); exec @ARGV' \
#       scripts/autoresearch/ar5/rescue_phase_b.sh > /dev/null 2>&1 < /dev/null &
# Jobs (3 heavy processes at once; each under caffeinate -i):
#   Gl    script arm Gl (rescue_gl.py: arm Al of step 8 + rule G1), then its pass
#         (order_pass.py --name Gl --aligned 0, the library's single generator)
#   LibG  library run, keys of S11L + SeedRescueCMA 0.25 1000 (rescue_lib_run.py --kind LibG)
#   LibS11 library run with the key off (control only: it must equal the stored script run S11 of
#         AR5 step 3 on mid; it says where a LibG/Gl difference comes from, if one appears)
# The names differ from Gl by more than case (macOS file names are case-insensitive).
# Logs and markers in <cache>/logs: <job>.log, <job>.done ("exit=<code>"); final marker
# phase_b.done with one line per job. Refuses (exit 2) if a run folder or the marker exists.
# Caps 240 minutes per run. Wall times are contended, not timing results. Plain bash 3.2.
set -u
cd "$(dirname "$0")/../../.." || exit 2
cache="${RESCUE_CACHE:-$PWD/benchmarks/autoresearch/ar5/rescue/cache}"
logs="$cache/logs"
mkdir -p "$logs"
export RESCUE_CACHE="$cache" ORDER_CACHE="$cache" RESEED_SAMPLE=mid
unset VIRTUAL_ENV
for f in "$cache/Gl_mid" "$cache/LibG_mid" "$cache/LibS11_mid" "$logs/phase_b.done"; do
  if [ -e "$f" ]; then echo "refusing: $f exists" > "$logs/phase_b.refused"; exit 2; fi
done
date +%Y-%m-%dT%H:%M:%S%z > "$logs/phase_b.start"
A5=scripts/autoresearch/ar5

( caffeinate -i uv run python $A5/rescue_gl.py --cap-min 240 > "$logs/Gl_mid.log" 2>&1
  rc=$?
  prc=skipped
  if [ "$rc" = "0" ]; then
    caffeinate -i uv run python $A5/order_pass.py --name Gl --aligned 0 \
      > "$logs/pass_Gl_mid.log" 2>&1
    prc=$?
  fi
  echo "exit=$rc pass=$prc" > "$logs/Gl_mid.done" ) &
( caffeinate -i uv run python $A5/rescue_lib_run.py --kind LibG --cap-min 240 \
    > "$logs/LibG_mid.log" 2>&1
  echo "exit=$?" > "$logs/LibG_mid.done" ) &
( caffeinate -i uv run python $A5/rescue_lib_run.py --kind LibS11 --cap-min 240 \
    > "$logs/LibS11_mid.log" 2>&1
  echo "exit=$?" > "$logs/LibS11_mid.done" ) &
wait
{
  date +%Y-%m-%dT%H:%M:%S%z
  for j in Gl_mid LibG_mid LibS11_mid; do echo "$j $(cat "$logs/$j.done" 2>/dev/null)"; done
} > "$logs/phase_b.done"
