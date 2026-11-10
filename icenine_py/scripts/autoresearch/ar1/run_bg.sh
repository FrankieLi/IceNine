#!/bin/sh
# Run one AR1 script in the background-safe way: log under cache/logs/<name>.log, marker <name>.done.
# usage (from icenine_py/): nohup scripts/autoresearch/ar1/run_bg.sh <name> <script.py> [args...] &
name="$1"
shift
logs=scripts/autoresearch/ar1/cache/logs
mkdir -p "$logs"
uv run python "$@" > "$logs/$name.log" 2>&1
touch "$logs/$name.done"
