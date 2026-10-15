#!/bin/bash
# S1 Step 5 timings: the base commit (a git worktree) and S1 run interleaved, so both arms share the
# machine's conditions: prof_eval (us per evaluate call, three arms), stage_split (the S1 parts) and
# search_split (whole no-start search: wall, evaluation count, result), ROUNDS rounds, the arm order
# alternating per round. Single worker (torch threads 1 in every script). The preflight state is
# saved with the results; label the timings "contended" unless preflight_exit.txt holds 0.
#
# Usage (from icenine_py/), BASE_WT = a worktree of the base commit (needs the built _rasterize
# extension copied in and Examples/ linked or copied; the scripts are copied in from this tree):
#   git worktree add --detach "$BASE_WT" 42672d1
#   S1_BASE_WT="$BASE_WT" nohup bash scripts/speed_parity/run_s1_timing.sh &
# Log: scripts/speed_parity/cache/logs/s1_timing.log, marker s1_timing.done.
cd "$(dirname "$0")/../.." || exit 1
BASE_WT=${S1_BASE_WT:?set S1_BASE_WT to a worktree of the base commit}
ROUNDS=${S1_ROUNDS:-5}
HERE=$(pwd)
OUT=scripts/speed_parity/cache/s1/timing
LOGD=scripts/speed_parity/cache/logs
mkdir -p "$OUT" "$LOGD"
rm -f "$LOGD/s1_timing.done" "$OUT"/*.txt
cp scripts/speed_parity/prof_eval.py scripts/speed_parity/stage_split.py \
  scripts/speed_parity/search_split.py "$BASE_WT/icenine_py/scripts/speed_parity/"
{
  uv run python ../scripts/dev/timing_preflight.py --require --json "$OUT/preflight_start.json" \
    > /dev/null 2>&1
  echo $? > "$OUT/preflight_exit.txt"
  for r in $(seq 1 "$ROUNDS"); do
    for tool in prof_eval stage_split search_split; do
      if [ $((r % 2)) -eq 1 ]; then arms="base s1"; else arms="s1 base"; fi
      for arm in $arms; do
        if [ "$arm" = base ]; then dir="$BASE_WT/icenine_py"; else dir="$HERE"; fi
        echo "round $r $tool $arm: $(uptime | sed 's/.*load/load/')"
        (cd "$dir" && uv run --project "$HERE" python "scripts/speed_parity/$tool.py") \
          > "$OUT/${arm}_${tool}_r${r}.txt" 2>&1
      done
    done
  done
  uv run python ../scripts/dev/timing_preflight.py --json "$OUT/preflight_end.json" > /dev/null 2>&1
} > "$LOGD/s1_timing.log" 2>&1
touch "$LOGD/s1_timing.done"
