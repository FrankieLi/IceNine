---
name: finish-feature
description: Finish a feature branch - run tests, review, and create a PR to develop.
user-invocable: true
allowed-tools: Bash(git *), Bash(gh *), Bash(pytest *), Bash(cd *), Bash(make *)
---

Finish the current feature branch and create a PR to merge into develop.

Follow these steps in order:

## Step 1: Identify branches
- Get the current branch name
- Verify it's a `feature/*` branch. If not, abort.

## Step 2: Run full test suite
- Run Python tests:
  ```
  cd icenine_py && pytest tests/ -v
  ```
- If C++ files were changed (`git diff develop...HEAD --name-only | grep -E '\.(cpp|h|tmpl\.cpp)$'`), build and verify:
  ```
  cmake -DCMAKE_BUILD_TYPE=Release . && make -j8
  ```
- If any tests fail, STOP and report. Do NOT proceed.

## Step 3: Review all changes
- Show the full diff against develop: `git diff develop...HEAD --stat`
- Show commit history: `git log --oneline develop..HEAD`
- Provide a summary of all changes, organized by topic
- Flag any concerns (large files, potential issues, missing tests)

## Step 4: Document what was accomplished
- Update `icenine_py/MIGRATION_HISTORY.md`:
  - Add or update the relevant phase section with a summary of all work done in this feature
  - Include key decisions, gotchas discovered, and modules added/changed
- If new modules or features were added, update `icenine_py/README.md` (project structure, dependencies, usage)
- If any plan files exist in `.claude/plans/`, delete them — MIGRATION_HISTORY.md is the permanent record
- Commit the documentation updates

## Step 5: Create PR (after user confirms)
- Push latest changes: `git push origin <branch>`
- Create PR to develop using `gh pr create`:
  - Title: concise summary of the feature
  - Body: summary of changes, test results, any notes
  - Base branch: `develop`
- Return the PR URL
