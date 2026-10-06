---
name: finish-task
description: Finish current task branch - run tests, review changes, merge back to parent feature branch.
user-invocable: true
---

Finish the current task branch and merge it back to the parent feature branch. Task branches are
named `feature/<parent>-<task>`; the parent is recorded in `git config branch.<task>.parent`.

## Step 1: Review changes (needs judgement, stays here)
- Find the parent: `git config branch.$(git rev-parse --abbrev-ref HEAD).parent`. If none is
  recorded, ask the user to confirm the guess (the branch name without its last `-<segment>`).
  If the branch is not a task branch, abort and suggest `/finish-feature`.
- Run `git diff <parent>...HEAD`, review it, summarise what changed and any concerns. Ask the
  user whether to proceed.

## Step 2: Document what was accomplished
- Append a brief summary to `icenine_py/MIGRATION_HISTORY.md` under the appropriate section.
- If a new module or feature was added, update `icenine_py/README.md`.
- If a plan file exists in `.claude/plans/`, delete it; the summary replaces it.
- Commit the documentation updates (stage files by explicit path) before merging.

## Step 3: Tests and merge (deterministic, after the user confirms)
Run `scripts/dev/finish_task.sh [--yes] [--push]`. It runs the Python suite
(`cd icenine_py && uv run pytest tests/ -q`) if `.py` files changed and the C++ Release build if
C++ files changed, stops on failure, merges into the parent with `--no-ff` and the message
"Merge task <task>: <subject of the last commit>" (override with `-m`), deletes the task branch
locally and remotely, and pushes the parent only with `--push`. Use `--no-tests` only when the
tests were just run.

Confirm the merge and show `git log --oneline -5`.
