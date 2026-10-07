#!/bin/bash
# Usage: start_task.sh <task-name> [--push]
# Create <current-feature-branch>-<task> off the current feature branch and record its parent.
set -euo pipefail
TASK=""
PUSH=0
for a in "$@"; do
    case "$a" in
        --push) PUSH=1 ;;
        -*) echo "unknown option: $a" >&2; exit 2 ;;
        *) TASK="$a" ;;
    esac
done
[ -n "$TASK" ] || { echo "usage: start_task.sh <task-name> [--push]" >&2; exit 2; }
CUR=$(git rev-parse --abbrev-ref HEAD)
case "$CUR" in
    feature/*) ;;
    *) echo "error: current branch '$CUR' is not a feature/* branch" >&2; exit 1 ;;
esac
NEW="$CUR-$TASK"
git checkout -b "$NEW"
git config "branch.$NEW.parent" "$CUR"
if [ $PUSH -eq 1 ]; then
    git push -u origin "$NEW"
fi
echo "task branch: $NEW (parent: $CUR)"
