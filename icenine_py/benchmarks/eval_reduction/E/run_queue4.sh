#!/bin/bash
# Task E: a second process on the 30-voxel lenient baseline, from the end of the list, once the
# E0 level-0-width landscape is done (keeps about 6 processes at once). From icenine_py/.
L=benchmarks/eval_reduction/E/cache/logs
while [ ! -f $L/landscape_e0_L0.done ]; do sleep 30; done
uv run python scripts/eval_reduction/smooth_rank_probe.py run --kind full --rule old --spec off \
    --tag base --workers 2 --reverse > $L/rank_full_old_base_rev.log 2>&1
