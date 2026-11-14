#!/bin/bash
# AR3 step 2: one run_layer.py run in the background, detached, with a log and a .done marker.
#   scripts/autoresearch/ar3/run_layer_bg.sh NAME KIND SAMPLE [extra run_layer.py args]
# Run from icenine_py/. Log: benchmarks/autoresearch/ar3/layer/cache/logs/NAME.log
set -u
NAME="$1"
KIND="$2"
SAMPLE="$3"
shift 3
L=benchmarks/autoresearch/ar3/layer/cache/logs
mkdir -p "$L"
rm -f "$L/$NAME.done"
nohup caffeinate -i bash -c "uv run python scripts/autoresearch/ar3/run_layer.py --name $NAME --kind $KIND --sample $SAMPLE $* > $L/$NAME.log 2>&1; echo done > $L/$NAME.done" > /dev/null 2>&1 &
echo "started $NAME"
