#!/bin/bash
# Read-only snapshot of running jobs, recent logs, done markers and git state.
ROOT=$(GIT_OPTIONAL_LOCKS=0 git rev-parse --show-toplevel 2>/dev/null) || ROOT=$(pwd)
cd "$ROOT" || exit 0

echo "== running jobs (pid, elapsed, command) =="
ps -Ao pid,etime,command |
    grep -E "uv run python scripts/|run_all\.sh" | grep -v grep || echo "(none)"

echo
echo "== logs modified in the last 48 h =="
found=0
for log in $(find icenine_py/scripts -path '*/cache/logs/*.log' -mtime -2 2>/dev/null | sort); do
    found=1
    echo "-- $log ($(stat -f '%Sm' -t '%Y-%m-%d %H:%M' "$log" 2>/dev/null || stat -c '%y' "$log"))"
    tail -n 3 "$log" | sed 's/^/   /'
done
[ $found -eq 1 ] || echo "(none)"

echo
echo "== done markers =="
find icenine_py/scripts -name '*.done' 2>/dev/null | sort | head -50
echo

echo "== git =="
BRANCH=$(GIT_OPTIONAL_LOCKS=0 git rev-parse --abbrev-ref HEAD 2>/dev/null)
echo "branch: $BRANCH"
STATUS=$(GIT_OPTIONAL_LOCKS=0 git status --short 2>/dev/null)
echo "status: $(printf '%s' "$STATUS" | grep -c . ) changed/untracked entries"
printf '%s\n' "$STATUS" | head -20
PARENT=$(GIT_OPTIONAL_LOCKS=0 git config "branch.$BRANCH.parent" 2>/dev/null)
echo "-- commits not on ${PARENT:-develop}:"
GIT_OPTIONAL_LOCKS=0 git log --oneline "${PARENT:-develop}..HEAD" 2>/dev/null | head -20
echo "-- commits not on develop:"
GIT_OPTIONAL_LOCKS=0 git log --oneline develop..HEAD 2>/dev/null | head -20
echo "-- unpushed commits:"
GIT_OPTIONAL_LOCKS=0 git log --oneline '@{upstream}..HEAD' 2>/dev/null | head -20 || true
echo "-- worktrees:"
GIT_OPTIONAL_LOCKS=0 git worktree list 2>/dev/null
exit 0
