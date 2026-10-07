#!/bin/bash
# Single-worker timing (U0, U1/U2), then the 10-worker contention runs. Run from icenine_py/ with
# nothing else heavy on the machine; one job at a time. Logs go to scripts/profiling/cache/logs.
set -e
mkdir -p scripts/profiling/cache/logs
L=scripts/profiling/cache/logs
uv run python scripts/profiling/prof_u0.py run --workers 1 --repeats 3 --tag w1 > $L/u0_w1.log 2>&1
uv run python scripts/profiling/prof_seeded.py run --workers 1 --repeats 3 --tag w1 > $L/seeded_w1.log 2>&1
uv run python scripts/profiling/prof_seeded.py run --workers 10 --repeats 1 --dirs 2 --tag w10 --warm light > $L/seeded_w10.log 2>&1
uv run python scripts/profiling/prof_u0.py run --workers 10 --repeats 1 --tag w10 --warm light --limit 20 > $L/u0_w10.log 2>&1
echo done > $L/run_all.done
