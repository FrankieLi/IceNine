#!/bin/bash
# AR5 step 5: check 4 (both keys on the first gradient sample against the stored S11) and the
# passes of check 5 (the library pass from the stored maps of S11 on the gradient sample and on
# big), one after the other, in the background.
#
#   scripts/autoresearch/ar5/keys_run_grad.sh [--dry]
#
# Run from icenine_py/. Steps: keys_run.py --kind S11L (gradient sample), keys_check.py,
# keys_pass.py on the gradient sample, keys_pass.py on big. Log:
# benchmarks/autoresearch/ar5/keys/cache/logs/grad_chain.log; marker: grad_chain.done (one line per
# step with its exit code, then "all done"). The whole chain is capped at 90 minutes: a step still
# running at the deadline is killed and the chain stops ("capped"). Refuses to start if the output
# folder of the run (cache/S11L_grad) or a marker already exists. "--dry" prints what would run and
# exits before anything is started or created.
set -u
dry=0
[ "${1:-}" = "--dry" ] && dry=1
cache="$PWD/benchmarks/autoresearch/ar5/keys/cache"
log="$cache/logs/grad_chain.log"
done_file="$cache/logs/grad_chain.done"
cap=5400  # seconds, for the whole chain
steps=(
  "run_grad|uv run python scripts/autoresearch/ar5/keys_run.py --kind S11L"
  "check_grad|uv run python scripts/autoresearch/ar5/keys_check.py --kind S11L"
  "pass_grad|uv run python scripts/autoresearch/ar5/keys_pass.py --kind S11L"
  "pass_big|env AR3_SAMPLE=big uv run python scripts/autoresearch/ar5/keys_pass.py --kind S11L"
)
if [ ! -f scripts/autoresearch/ar5/keys_run.py ]; then
  echo "run from icenine_py/" >&2; exit 2
fi
if [ -e "$cache/S11L_grad" ] || [ -e "$done_file" ] || [ -e "$cache/pass_check_S11L_grad" ] \
   || [ -e "$cache/pass_check_S11L_big" ]; then
  echo "refusing to start: an output folder or marker already exists under $cache" >&2; exit 3
fi
if [ "$dry" = 1 ]; then
  echo "dry run: nothing started, nothing created. Chain (cap ${cap} s in total, log $log):"
  for s in "${steps[@]}"; do echo "  ${s%%|*}: ${s#*|}"; done
  exit 0
fi
mkdir -p "$cache/logs"
(
  start=$(date +%s)
  : > "$done_file.tmp"
  for s in "${steps[@]}"; do
    name="${s%%|*}"; cmd="${s#*|}"
    left=$(( cap - ($(date +%s) - start) ))
    if [ "$left" -le 0 ]; then echo "$name capped (no time left)" >> "$done_file.tmp"; break; fi
    echo "=== $name: $cmd" >> "$log"
    nohup caffeinate -i $cmd >> "$log" 2>&1 &
    pid=$!
    ( sleep "$left"; kill "$pid" 2>/dev/null; pkill -P "$pid" 2>/dev/null ) &
    killer=$!
    wait "$pid"; rc=$?
    kill "$killer" 2>/dev/null
    echo "$name exit $rc" >> "$done_file.tmp"
    [ "$rc" -ne 0 ] && { echo "stopped after $name" >> "$done_file.tmp"; break; }
  done
  echo "all done" >> "$done_file.tmp"
  mv "$done_file.tmp" "$done_file"
) > /dev/null 2>&1 &
echo "started the chain (log $log, marker $done_file)"
