#!/bin/bash
# Single-worker timing pass (Fcompare), preflight-gated: each command retries (every 60 s, up to 40
# times) until the machine is quiet. Run alone on the machine, from icenine_py/.
D=benchmarks/eval_reduction/Fcompare
for kind in e0 full; do
  for i in $(seq 1 40); do
    uv run python scripts/eval_reduction/outward_compare.py timing --kind $kind > $D/cache/logs/timing_$kind.log 2>&1 && break
    grep -q MachineBusyError $D/cache/logs/timing_$kind.log || break
    sleep 60
  done
  echo "$kind: attempts $i" >> $D/cache/logs/timing_attempts.txt
done
touch $D/cache/logs/timing.done
