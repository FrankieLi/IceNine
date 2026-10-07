---
name: start-task
description: Start a sub-task branch off the current feature branch. Use for breaking work into smaller pieces.
user-invocable: true
allowed-tools: Bash(git *), Bash(scripts/dev/start_task.sh *)
argument-hint: "<task-name> [--push]"
---

Create a sub-task branch for: $ARGUMENTS

Run `scripts/dev/start_task.sh $ARGUMENTS` (kebab-case task name). The script:

1. checks the current branch is `feature/*` and aborts otherwise;
2. creates `<current-feature-branch>-<task>` (git cannot hold `feature/x/task` while `feature/x`
   exists) and records the parent with `git config branch.<new>.parent <current>`;
3. pushes with upstream only when `--push` is given.

Then confirm the new branch name and its parent to the user.
