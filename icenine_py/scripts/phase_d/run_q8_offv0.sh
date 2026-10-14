#!/bin/bash
# Q8 step 8: why do the lenient wrong counts of the Step 4 re-runs differ from the stored ones
# (35 / 27 against 40 / 25 at min_sin_eta 0.3; 41 / 25 against 39 / 18 at 0)? Hypothesis: the
# stored lenient runs used the vertex-0 level screen (library default before 2026-10-12); the
# re-runs use the centroid screen. This runs the lenient rule with the vertex-0 screen (arm
# `offv0` of e0_cpprule.py) on the clean and old realistic ("all") E0 images, 200 voxels, seed 0,
# at both min_sin_eta settings. One pool of 10 workers at a time ("contended multi-worker").
# Resumable: a finished (msin, variant) block is skipped.
#
#   cd icenine_py && nohup bash scripts/phase_d/run_q8_offv0.sh > scripts/phase_d/cache/logs/q8_offv0.log 2>&1 &
#
# Results: scripts/phase_d/cache/offv0[_ms0]_q8/v<voxel>_<variant>.npz (gitignored);
# benchmarks/q8_nodrop/runs/run_meta_run_offv0*_q8_*.json (tracked).
set -u
cd "$(dirname "$0")/../.."
LOGS=scripts/phase_d/cache/logs
META=benchmarks/q8_nodrop/runs
mkdir -p "$LOGS" "$META"
for msin in 0.3 0; do
    suffix=""
    [ "$msin" = "0" ] && suffix="_ms0"
    for var in clean all; do
        if [ -f "$META/run_meta_run_offv0${suffix}_q8_${var}.json" ]; then
            echo "skip offv0$suffix $var (done)"
            continue
        fi
        echo "== $(date '+%F %T') arm offv0 msin $msin variant $var"
        uv run python scripts/phase_d/e0_cpprule.py run --arm offv0 --msin "$msin" \
            --variant "$var" --tag q8 --workers 10 --meta-dir "$META" 2>&1 |
            grep -v "Warning\|warn(\|_read_step_size_file\|sample.rotate\|Loaded crystal\|Max Q calc"
        echo "== $(date '+%F %T') finished offv0$suffix $var"
    done
done
echo done > "$LOGS/q8_offv0.done"
