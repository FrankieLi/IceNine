#!/bin/bash
# Task E: the level-0-width landscape on E0, after the full-sample one. Run from icenine_py/.
L=benchmarks/eval_reduction/E/cache/logs
while [ ! -f $L/landscape_full_L0.done ]; do sleep 30; done
uv run python scripts/eval_reduction/smooth_landscape.py run --kind e0 --workers 2 --set L0 > $L/landscape_e0_L0.log 2>&1
