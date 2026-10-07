#!/bin/bash
# Usage: checkpoint.sh [note]  -- save job_status output to .claude/checkpoints/<stamp>.md
ROOT=$(git rev-parse --show-toplevel) || exit 1
cd "$ROOT" || exit 1
DIR=.claude/checkpoints
mkdir -p "$DIR"
OUT="$ROOT/$DIR/$(date +%Y%m%d-%H%M).md"
{
    echo "# Checkpoint $(date '+%Y-%m-%d %H:%M:%S')"
    echo
    [ -n "$1" ] && printf 'Note: %s\n\n' "$1"
    echo '```'
    "$ROOT/scripts/dev/job_status.sh"
    echo '```'
} >"$OUT"
echo "$OUT"
