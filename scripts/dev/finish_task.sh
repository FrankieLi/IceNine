#!/bin/bash
# Usage: finish_task.sh [--no-tests] [--yes] [--push] [-m <message>]
# Run tests/build if relevant files changed, merge the current task branch into its parent with
# --no-ff, delete the task branch (local, and remote if present). Review stays with the caller.
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
PARENT=$(git config "branch.$TASK_BRANCH.parent" || true)
if [ -z "$PARENT" ]; then
    PARENT="${TASK_BRANCH%-*}"
    if [ "$PARENT" = "$TASK_BRANCH" ] || ! git rev-parse --verify -q "$PARENT" >/dev/null; then
        echo "error: no recorded parent and guessed parent '$PARENT' does not exist" >&2
        exit 1
    fi
    if [ $YES -ne 1 ]; then
        printf "No recorded parent; guess '%s'. Proceed? [y/N] " "$PARENT"
        read -r ans
        [ "$ans" = "y" ] || exit 1
    fi
fi
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "error: uncommitted changes" >&2; exit 1; }
CHANGED=$(git diff --name-only "$PARENT...HEAD")
if [ $NO_TESTS -ne 1 ]; then
    if echo "$CHANGED" | grep -q '\.py$'; then
        (cd icenine_py && uv run pytest tests/ -q)
    fi
    if echo "$CHANGED" | grep -qE '\.(cpp|h|hpp|txt)$|CMakeLists'; then
        cmake -DCMAKE_BUILD_TYPE=Release . && make -j8
    fi
fi
TASK_NAME="${TASK_BRANCH#"$PARENT"-}"
[ -n "$MSG" ] || MSG="Merge task $TASK_NAME: $(git log -1 --format=%s "$TASK_BRANCH")"
git checkout "$PARENT"
git merge --no-ff "$TASK_BRANCH" -m "$MSG"
git branch -d "$TASK_BRANCH"
git config --remove-section "branch.$TASK_BRANCH" 2>/dev/null || true
if git ls-remote --exit-code --heads origin "$TASK_BRANCH" >/dev/null 2>&1; then
    git push origin --delete "$TASK_BRANCH"
fi
if [ $PUSH -eq 1 ]; then
    git push origin "$PARENT"
fi
git log --oneline -5
