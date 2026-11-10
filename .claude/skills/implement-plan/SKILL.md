---
name: implement-plan
description: Build the brief for an approved plan with make_brief.sh, launch the implementer in its own worktree in the background, then run verify_handback.sh on the hand-back.
user-invocable: true
allowed-tools: Bash(scripts/dev/make_brief.sh *), Bash(scripts/dev/verify_handback.sh *), Bash(git *), Agent, SendMessage
argument-hint: "<plan> <task> [section] [model] [effort]"
---

Only run for a plan the owner has approved (`make_brief.sh` refuses otherwise and prints why).

1. Find the current branch (`git rev-parse --abbrev-ref HEAD`); it is the parent of the task
   branch, which is named `<current>-<task>`. Run `scripts/dev/make_brief.sh --plan <plan> --task
   <task> [--section "<section>"] --trailer "<the Co-Authored-By line from the current system
   reminder>"`. If it refuses (unapproved plan, `icenine.maxHeavyProcs` unset, no trailer), show
   the error and stop.
2. If the session has a keep-awake tool (in the desktop app, `request_keep_awake`), request it
   until the session is idle, before launching: sleep stops background agents mid-run.
3. Launch the `implementer` agent in the background with `isolation: "worktree"`, with the brief
   exactly as printed (no additions, no edits); the brief's branch line creates the task branch in
   that worktree. Model and effort: those the plan or its approval decision names, else the
   agent's own default; pass the user's `[model]` and `[effort]` if given.
4. If the agent stops before handing back (sleep, app restart, API error): check its worktree
   (`git -C <worktree> log --oneline <current>..HEAD`, `git -C <worktree> status --short`) and
   resume the same agent with SendMessage, saying what is committed and what is not. Never start a
   second agent on the same branch.
5. When it returns, take the report path from its final message and run
   `scripts/dev/verify_handback.sh <parent> <report> --branch <task branch>`. Show its summary and
   the report path, with the agent's status line and headline numbers only. If the agent reports
   that a tool or hook refused an action, show that first: it is the owner's to resolve.
6. Suggest `/review-task` when there is no FAIL; with a FAIL, say so and give the failing check.
