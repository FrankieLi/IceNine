#!/bin/bash
# Task F: the E0 run (8 near-axis + 16 far-from-axis voxels, 20 directions, all methods), then the
# full-sample run (six interior voxels, 1, 2 and 5 deg, every method except F3). Commands from
# icenine_py/. Logs in benchmarks/eval_reduction/F/cache/logs/.
set -e
L=benchmarks/eval_reduction/F/cache/logs
uv run python scripts/eval_reduction/outward_bench.py run --kind e0 --workers 5 --dirs 20 \
  --n-near 8 --n-far 16 --tag e0 > $L/e0.log 2>&1
uv run python scripts/eval_reduction/outward_bench.py run --kind full --workers 5 --dirs 20 \
  --offsets 2,3,5 --tag full --methods f1,f2_05,f2_1,f2_2,f4,f5,f54 > $L/full.log 2>&1
# the extra method, run after seeing why F4 fails from far starts (RUN_LOG section 6)
uv run python scripts/eval_reduction/outward_bench.py run --kind e0 --workers 5 --dirs 20 \
  --n-near 8 --n-far 16 --tag e0_f4w --methods f4w > $L/e0_f4w.log 2>&1
uv run python scripts/eval_reduction/outward_bench.py run --kind full --workers 5 --dirs 20 \
  --offsets 2,3,5 --tag full_f4w --methods f4w > $L/full_f4w.log 2>&1
# settings pilots on the four pilot voxels (RUN_LOG 3), the half-rise medians from Task B's walks
uv run python scripts/eval_reduction/outward_pilot.py schedule > $L/pilot_schedule.log 2>&1
uv run python scripts/eval_reduction/outward_pilot.py spacing > $L/pilot_spacing.log 2>&1
uv run python scripts/eval_reduction/outward_pilot.py window > $L/pilot_window.log 2>&1
uv run python scripts/eval_reduction/outward_pilot.py halfrise > $L/pilot_halfrise.log 2>&1
# single-worker timing, preflight-gated (stops on a busy machine; retry when quiet)
uv run python scripts/eval_reduction/outward_bench.py timing --kind e0 --dirs 5 > $L/timing_e0.log 2>&1
uv run python scripts/eval_reduction/outward_bench.py timing --kind full --dirs 5 > $L/timing_full.log 2>&1
uv run python scripts/eval_reduction/outward_bench.py timing-f5 --kind e0 --dirs 5 > $L/timing_f5_e0.log 2>&1
uv run python scripts/eval_reduction/outward_bench.py timing-f5 --kind full --dirs 5 > $L/timing_f5_full.log 2>&1
# why F4 fails from far starts (RUN_LOG 6)
uv run python scripts/eval_reduction/outward_diag.py e0 --n-voxels 12 > $L/diag_e0.log 2>&1
uv run python scripts/eval_reduction/outward_diag.py full > $L/diag_full.log 2>&1
# summaries, tables and figures
uv run python scripts/eval_reduction/outward_summary.py e0
uv run python scripts/eval_reduction/outward_summary.py full
uv run python scripts/eval_reduction/outward_summary.py merge
