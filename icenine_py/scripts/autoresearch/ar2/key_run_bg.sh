#!/bin/bash
# AR2 step 2: one BFS pilot run in the background, detached, with a log and a .done marker.
#   scripts/autoresearch/ar2/key_run_bg.sh NAME CONFIG [SAMPLE.mic]
# Run from icenine_py/. Results: benchmarks/autoresearch/ar2/key/cache/runs/NAME (gitignored by
# not being added), log: benchmarks/autoresearch/ar2/key/cache/logs/NAME.log.
set -u
NAME="$1"
CONFIG="$2"
SAMPLE="${3:-}"
OUT=benchmarks/autoresearch/ar2/key/cache
mkdir -p "$OUT/logs" "$OUT/runs"
EXTRA=()
if [ -n "$SAMPLE" ]; then EXTRA=(--sample "$SAMPLE"); fi
LABEL="contended: timing preflight refused (an editor process at about 200% CPU), single run"
nohup caffeinate -i bash -c "uv run python scripts/phase_d/pilot_bfs.py --arm mc --variant clean --out $OUT/runs --name $NAME --config $CONFIG --timing-label '$LABEL' ${EXTRA[*]:-} > $OUT/logs/$NAME.log 2>&1; echo done > $OUT/logs/$NAME.done" > /dev/null 2>&1 &
echo "started $NAME"
