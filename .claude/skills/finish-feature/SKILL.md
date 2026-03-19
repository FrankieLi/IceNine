---
name: finish-feature
description: Finish a feature branch - run tests, review, create PR to develop, and merge via GitHub.
user-invocable: true
---

Finish the current feature branch and merge into develop via a GitHub PR.

**CRITICAL**: NEVER merge locally and push develop directly. Always create the PR first, then merge through GitHub so there is a permanent record.

Follow these steps in order:

## Step 1: Identify branches
- Get the current branch name
- Verify it's a `feature/*` branch (not a task sub-branch). If it's a task branch, abort and suggest `/finish-task`.

## Step 2: Run full test suite
- Run Python tests:
  ```
  cd icenine_py && uv run pytest tests/ -v
  ```
- If C++ files were changed (`git diff develop...HEAD --name-only | grep -E '\.(cpp|h|tmpl\.cpp)$'`), build and verify:
  ```
  cmake -DCMAKE_BUILD_TYPE=Release . && make -j8
  ```
- If any tests fail, STOP and report. Do NOT proceed.

## Step 3: Code review (separate context)
- Spawn the `code-reviewer` agent using the Agent tool with:
  - `subagent_type`: `"code-reviewer"`
  - `model`: `"opus"`
  - `prompt`: `"Review the feature branch <branch> against develop. Run: git diff develop...HEAD --stat, git log --oneline develop..HEAD, then read and review all changed files."`
- Wait for the review to complete
- Present the review results to the user
- If the verdict is REQUEST CHANGES, STOP and work with the user to address the issues before proceeding

## Step 4: Document what was accomplished
- Update `icenine_py/MIGRATION_HISTORY.md`:
  - Add or update the relevant phase section with a summary of all work done in this feature
  - Include key decisions, gotchas discovered, and modules added/changed
- If new modules or features were added, update `icenine_py/README.md` (project structure, dependencies, usage)
- If any plan files exist in `.claude/plans/`, delete them — MIGRATION_HISTORY.md is the permanent record
- Commit the documentation updates

## Step 5: Push and create PR
- Push latest changes: `git push -u origin <branch>`
- Create PR to develop using `gh pr create`:
  - Title: concise summary of the feature (under 70 chars)
  - Body: use this format:
    ```
    ## Summary
    <bullet points of key changes>

    ## Test results
    <pass/fail counts, integration test results>

    ## Notes
    <any gotchas, breaking changes, or follow-up items>
    ```
  - Base branch: `develop`
- Return the PR URL to the user

## Step 6: Merge PR via GitHub (after user confirms)
- Merge the PR using: `gh pr merge <number> --merge --delete-branch`
  - This performs a merge commit (matching gitflow's `--no-ff`) and deletes the remote branch
- Pull develop locally: `git checkout develop && git pull origin develop`
- Delete the local feature branch: `git branch -d <branch>`
- Confirm with `git log --oneline -5`
