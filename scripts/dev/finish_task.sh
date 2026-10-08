#!/bin/bash
# Usage: finish_task.sh [--no-tests] [--yes] [--push] [-m <message>] [--next <task>]
# Run tests/build if relevant files changed, merge the current task branch into its parent with
# --no-ff, delete the local task branch. With --push, push the parent and then delete the remote
# task branch. A dirty CLAUDE.md is ignored for the cleanliness check (never staged). Any other
# worktree holding the task branch is detached and removed; --next <task> then starts the next
# task off the parent. Review stays with the caller. The commands can be overridden for testing with
# FINISH_TASK_TEST_CMD and FINISH_TASK_BUILD_CMD.
set -euo pipefail
SELF_DIR=$(cd "$(dirname "$0")" && pwd)
NO_TESTS=0
YES=0
PUSH=0
MSG=""
NEXT=""
while [ $# -gt 0 ]; do
    case "$1" in
        --no-tests) NO_TESTS=1 ;;
        --yes) YES=1 ;;
        --push) PUSH=1 ;;
        -m) shift; MSG="${1:-}" ;;
        --next) shift; NEXT="${1:-}" ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done
ROOT=$(git rev-parse --show-toplevel)
cd "$ROOT"
TASK_BRANCH=$(git rev-parse --abbrev-ref HEAD)
case "$TASK_BRANCH" in
    feature/*) ;;
    *) echo "error: current branch '$TASK_BRANCH' is not a feature/* task branch" >&2; exit 1 ;;
esac
PARENT=$(git config "branch.$TASK_BRANCH.parent" || true)
if [ -z "$PARENT" ]; then
    if [ $YES -eq 1 ]; then
        echo "error: no recorded parent for '$TASK_BRANCH'; --yes does not guess" >&2
        exit 1
    fi
    PARENT="${TASK_BRANCH%-*}"
    if [ "$PARENT" = "$TASK_BRANCH" ] || ! git rev-parse --verify -q "$PARENT" >/dev/null; then
        echo "error: no recorded parent and guessed parent '$PARENT' does not exist" >&2
        exit 1
    fi
    printf "No recorded parent; guess '%s'. Proceed? [y/N] " "$PARENT"
    read -r ans
    [ "$ans" = "y" ] || exit 1
fi
case "$PARENT" in
    develop | master | main) echo "error: refusing to merge a task into '$PARENT'" >&2; exit 1 ;;
esac
if [ -n "$(git status --porcelain --untracked-files=no -- . ':(exclude)CLAUDE.md')" ]; then
    echo "error: uncommitted changes" >&2
    exit 1
fi
CHANGED=$(git diff --name-only "$PARENT...HEAD")
if [ $NO_TESTS -ne 1 ]; then
    if echo "$CHANGED" | grep -q '\.py$'; then
        bash -c "${FINISH_TASK_TEST_CMD:-cd icenine_py && uv run pytest tests/ -q}" || exit 1
    fi
    if echo "$CHANGED" | grep -qE '\.(cpp|h|hpp)$|(^|/)CMakeLists\.txt$'; then
        bash -c "${FINISH_TASK_BUILD_CMD:-cmake -DCMAKE_BUILD_TYPE=Release . && make -j8}" || exit 1
    fi
fi
TASK_NAME="${TASK_BRANCH#"$PARENT"-}"
[ -n "$MSG" ] || MSG="Merge task $TASK_NAME: $(git log -1 --format=%s "$TASK_BRANCH")"
# Where to merge: the worktree that already holds the parent, else this one.
wt_of() {  # print the worktree path holding refs/heads/$1, if any
    git worktree list --porcelain |
        awk -v b="refs/heads/$1" '/^worktree /{w=substr($0,10)} $1=="branch" && $2==b {print w}'
}
CUR_WT=$(pwd -P)
PARENT_WT=$(wt_of "$PARENT")
if [ -n "$PARENT_WT" ] && [ "$(cd "$PARENT_WT" && pwd -P)" != "$CUR_WT" ]; then
    MERGE_WT=$(cd "$PARENT_WT" && pwd -P)
    if [ -n "$(git -C "$MERGE_WT" status --porcelain --untracked-files=no -- . ':(exclude)CLAUDE.md')" ]; then
        echo "error: uncommitted changes in the parent worktree" >&2
        exit 1
    fi
else
    MERGE_WT="$CUR_WT"
    git checkout "$PARENT"
fi
git -C "$MERGE_WT" merge --no-ff "$TASK_BRANCH" -m "$MSG"
# free the task branch from any other worktree (detach, then remove if clean)
TASK_WT=$(wt_of "$TASK_BRANCH")
if [ -n "$TASK_WT" ] && [ "$(cd "$TASK_WT" && pwd -P)" != "$MERGE_WT" ]; then
    git -C "$TASK_WT" checkout --detach -q
    cd "$MERGE_WT"
    git worktree remove "$TASK_WT" || echo "warning: could not remove worktree $TASK_WT" >&2
fi
cd "$MERGE_WT"
git branch -d "$TASK_BRANCH"
git config --remove-section "branch.$TASK_BRANCH" 2>/dev/null || true
if [ $PUSH -eq 1 ]; then
    git push origin "$PARENT"
    if git ls-remote --exit-code --heads origin "$TASK_BRANCH" >/dev/null 2>&1; then
        git push origin --delete "$TASK_BRANCH"
    fi
fi
git log --oneline -5
if [ -n "$NEXT" ]; then
    "$SELF_DIR/start_task.sh" "$NEXT"
fi
