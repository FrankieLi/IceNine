---
name: start-task
description: Start a sub-task branch off the current feature branch. Use for breaking work into smaller pieces.
user-invocable: true
allowed-tools: Bash(git *)
argument-hint: "<task-name>"
---

Create a sub-task branch for: $ARGUMENTS

Follow these steps exactly:

1. Verify the current branch is a `feature/*` branch. If not, abort and tell the user to checkout a feature branch first.
2. Save the current branch name as the parent (e.g., `feature/experiment_setup_port`)
3. Create the task branch: `git checkout -b <parent-branch>/$ARGUMENTS`
   - Example: if on `feature/experiment_setup_port` and task is `fix-euler`, create `feature/experiment_setup_port/fix-euler`
4. Push and set upstream: `git push -u origin <new-branch>`
5. Confirm and show the branch hierarchy

The task branch name should use kebab-case.
