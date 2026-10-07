#!/bin/bash
# Usage: finish_task.sh [--no-tests] [--yes] [--push] [-m <message>]
# Run tests/build if relevant files changed, merge the current task branch into its parent with
# --no-ff, delete the local task branch. With --push, push the parent and then delete the remote
# task branch. Review stays with the caller. The commands can be overridden for testing with
# FINISH_TASK_TEST_CMD and FINISH_TASK_BUILD_CMD.
set -euo pipefail
NO_TESTS=0
YES=0
PUSH=0
MSG=""
while [ $# -gt 0 ]; do
    case "$1" in
        --no-tests) NO_TESTS=1 ;;
        --yes) YES=1 ;;
        --push) PUSH=1 ;;
        -m) shift; MSG="${1:-}" ;;
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
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
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
git checkout "$PARENT"
git merge --no-ff "$TASK_BRANCH" -m "$MSG"
git branch -d "$TASK_BRANCH"
git config --remove-section "branch.$TASK_BRANCH" 2>/dev/null || true
if [ $PUSH -eq 1 ]; then
    git push origin "$PARENT"
    if git ls-remote --exit-code --heads origin "$TASK_BRANCH" >/dev/null 2>&1; then
        git push origin --delete "$TASK_BRANCH"
    fi
fi
git log --oneline -5
