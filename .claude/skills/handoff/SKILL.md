---
name: handoff
description: Checkpoint the session to a file so the user can /clear and resume cheaply with /resume.
user-invocable: true
allowed-tools: Bash(scripts/dev/checkpoint.sh *), Bash(scripts/dev/job_status.sh), Bash(git *), Read, Edit
argument-hint: "[note]"
---

Write a checkpoint so a fresh session can continue with `/resume`. Keep it short.

1. Pick the memory file for the current work (usually the `project_*.md` file named in
   `~/.claude/projects/-Users-sfli-Research-IceNine/memory/MEMORY.md` for this branch).
2. Run `scripts/dev/checkpoint.sh --memory-file <that file> "$ARGUMENTS"`. It prints the checkpoint
   path (under `.claude/checkpoints/`, gitignored) and fills branch, parent, HEAD, commits ahead
   of the parent, `git status --short`, worktrees, running jobs, active plans in `~/.claude/plans/`
   and the latest `.claude/reports/` files. It also keeps one line `- Latest checkpoint: <path>`
   in the memory file; that pointer is the only memory edit.
3. Edit the three placeholder sections of the checkpoint by hand, briefly (the whole file stays
   under about 80 lines):
   - "Done this session": 3-6 bullets, with commit hashes and report paths, no result prose.
   - "Waiting on the owner": decisions, commits only the owner may make (`.claude/` files).
   - "Next step": exact commands or agent briefs to run (briefs may say "follow
     study-conventions" and name the plan and report paths).
4. Tell the user: the checkpoint path, and that it is safe to `/clear` and then run `/resume`.
   If the script printed a WARNING, or any background agent is still running, say so first:
   clearing does not stop them and their reports will land in `.claude/reports/`.
