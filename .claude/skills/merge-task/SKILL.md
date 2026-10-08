---
name: merge-task
description: Merge the current task branch into its parent (--no-ff), push the parent, delete the branch and worktree, optionally start the next task.
user-invocable: true
allowed-tools: Bash(git *), Bash(scripts/dev/finish_task.sh *)
argument-hint: "[next-task]"
---

Only run when the user has asked to merge. Review stays with `/review-task` and `/finish-task`.

1. Show `git status --short` and `git log --oneline <parent>..HEAD`. `scripts/dev/finish_task.sh`
   ignores a dirty `CLAUDE.md` (never stages it) and refuses any other uncommitted change; if it
   refuses, stop and report, do not stash or commit on your own.
2. Run `scripts/dev/finish_task.sh --yes --push ${ARGUMENTS:+--next $ARGUMENTS}`. Add
   `--no-tests` only when the suite was just run; otherwise changed `.py`/C++ files
   trigger tests/build. The script merges with `--no-ff` and the message "Merge task <task>:
   <last commit subject>", pushes the parent, deletes the task branch (local and remote), detaches
   and removes any other worktree holding the task branch, and with `--next` starts the next task
   off the parent. It requires a recorded parent and never merges into develop or master.
3. Show `git log --oneline -5` and the current branch.
4. If this finished a phase (no further task planned), suggest `/handoff`.
