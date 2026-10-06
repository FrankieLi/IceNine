#!/bin/bash
# F1 post hoc, E2 dataset, then the re-run fixes (seed 0), all with 10 workers
cd "$(dirname "$0")/../.."
L=scripts/findoptimal_robustness/cache
uv run python scripts/findoptimal_robustness/f1_run.py run --workers 10 > $L/f1.log 2>&1
uv run python scripts/findoptimal_robustness/e2_dataset.py run --workers 10 > $L/e2_dataset.log 2>&1
for F in F1b F2a F2b F3a F3b F3c; do
  uv run python scripts/findoptimal_robustness/fixes_run.py run --fix $F --seeds 0 --workers 10 > $L/fix_$F.log 2>&1
done
echo CHAIN1DONE > $L/chain1.done
