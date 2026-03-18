---
name: finish-task
description: Finish current task branch - run tests, review changes, merge back to parent feature branch.
user-invocable: true
allowed-tools: Bash(git *), Bash(pytest *), Bash(cd *), Bash(make *)
---

Finish the current task branch and merge it back to the parent feature branch.

Follow these steps in order:

## Step 1: Identify branches
- Get the current branch name
- The parent branch is everything up to the last `/` segment (e.g., `feature/foo/fix-bar` → parent is `feature/foo`)
- If the current branch is not a task branch (doesn't have 3+ path segments), abort and suggest using `/finish-feature` instead

## Step 2: Run tests
- If there are Python files changed (`git diff --name-only <parent>...HEAD | grep '\.py$'`), run:
  ```
  cd icenine_py && pytest tests/ -v
  ```
- If there are C++ files changed, run:
  ```
  cmake -DCMAKE_BUILD_TYPE=Release . && make -j8
  ```
- If tests fail, STOP and report the failures. Do NOT proceed to merge.

## Step 3: Review changes
- Run `git diff <parent>...HEAD` and review the changes
- Provide a brief summary of what changed and any concerns
- Ask the user if they want to proceed with the merge

## Step 4: Merge (only after user confirms)
- Checkout the parent branch: `git checkout <parent>`
- Merge the task branch: `git merge --no-ff <task-branch> -m "Merge <task-branch>: <brief summary>"`
- Push the parent branch: `git push origin <parent>`
- Delete the task branch locally and remotely:
  ```
  git branch -d <task-branch>
  git push origin --delete <task-branch>
  ```
- Confirm the merge and show `git log --oneline -5`
