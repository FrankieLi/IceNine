---
name: resume
description: First action in a fresh session. Read the newest checkpoint and confirm state cheaply.
user-invocable: true
allowed-tools: Bash(git *), Bash(ls *), Bash(scripts/dev/job_status.sh), Read, Grep
---

Resume from the latest checkpoint with minimal context use.

1. Read only these:
   - the newest file in `.claude/checkpoints/` (`ls -t .claude/checkpoints/*.md | head -1`; if the
     directory is empty, say so and stop);
   - the plan file(s) the checkpoint names, only the one for the next step;
   - `git status --short` and `git log --oneline -5`;
   - the report files the checkpoint lists, only if the next step needs them (read the headline
     and "Required fixes", not the whole file).
2. Do not read `MIGRATION_HISTORY.md`, the README or other large docs wholesale. When a section is
   needed, find it with `grep -n "<heading>" <doc>` and read that line range only.
3. Reply in at most 10 lines: branch and parent, whether git state matches the checkpoint, running
   jobs (`scripts/dev/job_status.sh` only if the checkpoint lists jobs), what is waiting on the
   owner, and the proposed next step.
4. Ask before starting any long work (agent launches, background jobs, merges).
