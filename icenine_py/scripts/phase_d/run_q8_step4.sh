#!/bin/bash
# Q8 step 4: E0 (200 voxels, seed 0) rows of docs/todo_q8_noise_nodrop.md on three image variants
# (clean, all = old realistic, realistic_nodrop), arms off/on (min_sin_eta 0.3: Q3) and
# off/on/rayR2/rayA (min_sin_eta 0: ray-overlap study). One pool of 10 workers at a time, so every
# wall is "contended multi-worker". Resumable: a finished (arm, variant) block is skipped.
#
#   cd icenine_py && nohup bash scripts/phase_d/run_q8_step4.sh > scripts/phase_d/cache/logs/q8_step4.log 2>&1 &
#
# Results: scripts/phase_d/cache/<arm>[_ms0]_q8/v<voxel>_<variant>.npz (gitignored);
# benchmarks/q8_nodrop/runs/run_meta_*.json and preflight_<arm>.json (tracked).
set -u
cd "$(dirname "$0")/../.."
LOGS=scripts/phase_d/cache/logs
META=benchmarks/q8_nodrop/runs
mkdir -p "$LOGS" "$META"
# wait (up to 10 min) for a quiet machine before the first arm; record the state either way
for i in $(seq 1 20); do
    uv run python ../scripts/dev/timing_preflight.py --require --json "$META/preflight_start.json" \
        > /dev/null 2>&1 && break
    sleep 30
done
ARMS="off:0.3 on:0.3 off:0 on:0 rayR2:0 rayA:0"
for spec in $ARMS; do
    arm=${spec%%:*}
    msin=${spec##*:}
    suffix=""
    [ "$msin" = "0" ] && suffix="_ms0"
    uv run python ../scripts/dev/timing_preflight.py --json "$META/preflight_${arm}${suffix}.json" \
        > /dev/null 2>&1
    for var in clean all realistic_nodrop; do
        if [ -f "$META/run_meta_run_${arm}${suffix}_q8_${var}.json" ]; then
            echo "skip $arm$suffix $var (done)"
            continue
        fi
        echo "== $(date '+%F %T') arm $arm msin $msin variant $var"
        uv run python scripts/phase_d/e0_cpprule.py run --arm "$arm" --msin "$msin" \
            --variant "$var" --tag q8 --workers 10 --meta-dir "$META" 2>&1 |
            grep -v "Warning\|warn(\|_read_step_size_file\|sample.rotate\|Loaded crystal\|Max Q calc"
        echo "== $(date '+%F %T') finished $arm$suffix $var"
    done
done
echo done > "$LOGS/q8_step4.done"
