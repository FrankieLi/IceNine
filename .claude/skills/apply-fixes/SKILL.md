---
name: apply-fixes
description: Brief a fresh implementer to apply a review report's required fixes (make_brief.sh --fixes) in its own worktree, then run verify_handback.sh with --fixes.
user-invocable: true
allowed-tools: Bash(scripts/dev/make_brief.sh *), Bash(scripts/dev/verify_handback.sh *), Bash(git *), Agent, SendMessage
argument-hint: "<review report> <plan> <task> [section] [model] [effort]"
---

Track S rule: the fixes go to a fresh implementer with no owner review in between.

1. Find the current branch and the task branch the review covered (the review report names it;
   ask the user if it does not, never guess). Run `scripts/dev/make_brief.sh --plan <plan> --task
   <task> --fixes <review report> [--section "<section>"] --trailer "<the Co-Authored-By line
   from the current system reminder>"`. If it refuses, show the error and stop.
2. Free the task branch, so the new agent can check it out in its own worktree (a subagent's Edit
   and Write tools work only inside its own worktree; never point it at another worktree's path).
   Find the worktree holding the branch (`git worktree list --porcelain`). If that worktree has
   uncommitted changes (`git -C <worktree> status --short`), stop and ask the user: those files
   would not reach the new agent. Uncommitted `.claude/**` files (left for the owner) are the
   exception: list them, and if a required fix touches them, apply that fix here instead and say
   so. Then run `git -C <worktree> checkout --detach`.
3. If the session has a keep-awake tool (in the desktop app, `request_keep_awake`), request it
   until the session is idle.
4. Launch a fresh `implementer` agent in the background with `isolation: "worktree"` and the brief
   exactly as printed. Model and effort as in `/implement-plan`. If it stops before handing back,
   resume it as in `/implement-plan` step 4.
5. When it returns, run `scripts/dev/verify_handback.sh <parent> <fix report> --branch <task
   branch> --fixes <review report>` (the parent is the task branch's recorded parent,
   `git config branch.<task branch>.parent`). Show its summary and the fix report path, and any
   refused tool or hook action the agent reports, first.
6. Suggest `/review-task` again when there is no FAIL.
