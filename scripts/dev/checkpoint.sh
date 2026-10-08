#!/bin/bash
# Usage: checkpoint.sh [--memory-file FILE] [note]
# Write .claude/checkpoints/<stamp>.md (gitignored): branch, parent, HEAD, commits ahead, status,
# worktrees, running jobs, active plans and latest reports, plus three sections the main session
# fills by hand. With --memory-file, one pointer line to the checkpoint is kept in that file.
# Env overrides (tests): CHECKPOINT_PLANS_DIR (default ~/.claude/plans).
ROOT=$(git rev-parse --show-toplevel) || exit 1
cd "$ROOT" || exit 1
MEM=""
NOTE=""
while [ $# -gt 0 ]; do
    case "$1" in
        --memory-file) shift; MEM="${1:-}" ;;
        *) NOTE="$1" ;;
    esac
    shift
done
PLANS="${CHECKPOINT_PLANS_DIR:-$HOME/.claude/plans}"
DIR=.claude/checkpoints
mkdir -p "$DIR"
OUT="$ROOT/$DIR/$(date +%Y%m%d-%H%M%S).md"
BRANCH=$(git rev-parse --abbrev-ref HEAD)
PARENT=$(git config "branch.$BRANCH.parent" 2>/dev/null)
[ -n "$PARENT" ] || PARENT="(none recorded)"
BASE="$PARENT"
git rev-parse --verify -q "$BASE" >/dev/null || BASE=develop
JOBS=$("$ROOT/scripts/dev/job_status.sh" 2>/dev/null |
    sed -n '/== running jobs/,/== logs modified/p' | grep -v '^==' | grep -v '^$' | head -10)
{
    echo "# Checkpoint $(date '+%Y-%m-%d %H:%M:%S')"
    [ -n "$NOTE" ] && printf '\nNote: %s\n' "$NOTE"
    echo
    echo "## State"
    echo "- branch: $BRANCH"
    echo "- parent: $PARENT"
    echo "- HEAD: $(git log -1 --format='%h %s' | cut -c1-100)"
    echo "- worktree: $ROOT"
    echo
    echo "## Commits ahead of $BASE (max 15)"
    echo '```'
    git log --oneline "$BASE..HEAD" 2>/dev/null | head -15
    echo '```'
    echo
    echo "## git status --short (max 20)"
    echo '```'
    git status --short | head -20
    echo '```'
    echo
    echo "## Worktrees"
    echo '```'
    git worktree list
    echo '```'
    echo
    echo "## Running jobs"
    echo '```'
    printf '%s\n' "$JOBS"
    echo '```'
    echo
    echo "## Active plans (newest 3)"
    if [ -d "$PLANS" ]; then
        ls -t "$PLANS"/*.md 2>/dev/null | head -3 | sed 's/^/- /'
    fi
    echo
    echo "## Latest reports (newest 5)"
    ls -t "$ROOT"/.claude/reports/*.md 2>/dev/null | head -5 | sed "s|^$ROOT/|- |"
    echo
    echo "## Done this session"
    echo "<fill by hand, 3-6 bullets>"
    echo
    echo "## Waiting on the owner"
    echo "<fill by hand: decisions, commits only the owner may make>"
    echo
    echo "## Next step"
    echo "<fill by hand: exact commands or agent briefs>"
} >"$OUT"
if [ -n "$MEM" ]; then
    touch "$MEM"
    { grep -v '^- Latest checkpoint:' "$MEM" || true; } >"$MEM.tmp"
    printf -- '- Latest checkpoint: %s (resume with /resume)\n' "$OUT" >>"$MEM.tmp"
    mv "$MEM.tmp" "$MEM"
fi
if [ -n "$JOBS" ] && [ "$JOBS" != "(none)" ]; then
    echo "WARNING: background job(s) still running; /clear does not stop them:" >&2
    printf '%s\n' "$JOBS" >&2
fi
echo "$OUT"
