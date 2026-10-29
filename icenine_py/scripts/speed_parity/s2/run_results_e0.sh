#!/bin/bash
# S2 Step 6: the E0 A/B of FusedCost (E0 clean, 200 voxels, seed 0), one arm after the other with
# S2_E0_WORKERS workers (default 6, the machine limit): key off (_s2off), FusedCost (_s2on),
# BatchedQuickMC (_s2bq), both (_s2bqon); then the key-off arm of the base commit from a git archive
# copy (BASE = its icenine_py directory, built in place; PYTHONPATH points the workers at it), into
# BASE's own cache (e0_s2base). Then the aggregates (s2_e0.py), the traces of the wrong flips
# (s2_e0_trace.py) and compare-e0 of the key-off arm against the base.
#
# Usage (from icenine_py/):
#   S2_BASE=/path/to/base/icenine_py nohup caffeinate -i bash scripts/speed_parity/s2/run_results_e0.sh \
#       > /dev/null 2>&1 &
# Log: scripts/speed_parity/cache/logs/s2_results_e0.log, marker s2_results_e0.done.
cd "$(dirname "$0")/../../.." || exit 1
ICE="$(pwd)"
BASE=${S2_BASE:?set S2_BASE to the base copy's icenine_py directory}
W=${S2_E0_WORKERS:-6}
OUT=benchmarks/s2_results
LOGD=scripts/speed_parity/cache/logs
mkdir -p "$OUT" "$LOGD"
rm -f "$LOGD/s2_results_e0.done"
RUN="uv run python scripts/findoptimal_robustness/e0_run.py run --variants clean --n-seeds 1"
RUN="$RUN --n-voxels 200 --workers $W"
{
  echo "start $(date) $(uptime | sed 's/.*load/load/')"
  $RUN --tag _s2off
  echo "off done $(date)"
  $RUN --fused-cost --tag _s2on
  echo "on done $(date)"
  $RUN --batched-quick-mc --tag _s2bq
  echo "bq done $(date)"
  $RUN --batched-quick-mc --fused-cost --tag _s2bqon
  echo "bqon done $(date)"
  PYTHONPATH="$BASE" uv run python "$BASE/scripts/findoptimal_robustness/e0_run.py" run \
    --variants clean --n-seeds 1 --n-voxels 200 --workers "$W" --tag _s2base
  echo "base done $(date)"
  uv run python scripts/speed_parity/s1_parity.py compare-e0 \
    scripts/findoptimal_robustness/cache/e0_s2off "$BASE/scripts/findoptimal_robustness/cache/e0_s2base" \
    > "$OUT/e0_keyoff_compare.txt" 2>&1
  uv run python scripts/speed_parity/s2/s2_e0.py --workers "$W" --out "$OUT/e0.json" \
    --base "$BASE/scripts/findoptimal_robustness/cache/e0_s2base"
  uv run python scripts/speed_parity/s2/s2_e0_trace.py --out "$OUT/e0_trace.json"
  echo "end $(date) $(uptime | sed 's/.*load/load/')"
} > "$LOGD/s2_results_e0.log" 2>&1
touch "$LOGD/s2_results_e0.done"
