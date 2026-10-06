#!/bin/bash
# after F1b: remaining re-run fixes, then the classifier end-to-end runs
cd "$(dirname "$0")/../.."
L=scripts/findoptimal_robustness/cache
while pgrep -f "fixes_run.py run --fix F1b" > /dev/null; do sleep 10; done
for F in F2a F3a F2b F3b; do
  uv run python scripts/findoptimal_robustness/fixes_run.py run --fix $F --seeds 0 --workers 10 > $L/fix_$F.log 2>&1
done
uv run python scripts/findoptimal_robustness/e2_endtoend.py final --model GBT --workers 10 > $L/e2_final.log 2>&1
uv run python scripts/findoptimal_robustness/e2_endtoend.py rerank --model GBT --workers 10 > $L/e2_rerank.log 2>&1
echo CHAIN2DONE > $L/chain2.done
