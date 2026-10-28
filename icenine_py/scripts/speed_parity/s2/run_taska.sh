#!/usr/bin/env bash
# S2 Task A check: regenerate Task A's inputs, check them, run torch and the fused prototype on
# them, trace the Step 0 float64 outlier, and write the summary (benchmarks/s2_taska/).
# Run from icenine_py/. About 1 minute; one process at a time. No timing is taken.
set -euo pipefail
cd "$(dirname "$0")/../../.."
ICE="$(pwd)"
D="$ICE/scripts/speed_parity/cache/s2_taska/inputs"
L="$ICE/scripts/speed_parity/cache/logs"
O="$ICE/benchmarks/s2_taska"
mkdir -p "$D" "$L" "$O"
rm -f "$L/s2_taska.done"
uv run python scripts/eval_reduction/parity_inputs.py --out "$D" > "$L/s2_taska_parity_inputs.log" 2>&1
uv run python scripts/speed_parity/s2/taska_check.py inputs --dir "$D" --out "$O/inputs.json" \
    > "$L/s2_taska_inputs.log" 2>&1
uv run python scripts/speed_parity/s2/taska_check.py run --dir "$D" --out "$O/run.json.gz" \
    > "$L/s2_taska_run.log" 2>&1
uv run python scripts/speed_parity/s2/taska_check.py outlier --min-px 0.02 --out "$O/outlier.json" \
    > "$L/s2_taska_outlier.log" 2>&1
uv run python scripts/speed_parity/s2/taska_summary.py --dir "$O" > "$L/s2_taska_summary.log" 2>&1
date > "$L/s2_taska.done"
