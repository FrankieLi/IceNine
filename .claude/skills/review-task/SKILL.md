---
name: review-task
description: Launch the code-reviewer (opus) on the current task branch against its parent; findings go to a report file.
user-invocable: true
allowed-tools: Bash(git *), Bash(ls *), Agent
argument-hint: "[extra focus]"
---

1. Find the branch and parent: `git rev-parse --abbrev-ref HEAD` and
   `git config branch.<branch>.parent` (if none is recorded, ask the user; never guess).
2. Launch the `code-reviewer` agent with model opus, in the foreground unless the user wants it
   in the background. Brief (keep it this short): "Follow study-conventions. Review <branch>
   against <parent> with the standard claims audit (audit_numbers, sync_doc_tables --check, then
   causal wording, costs, small counts, timing labels, leakage and pairing, C++ parity checked in
   the C++ source). Extra focus: $ARGUMENTS. Write the full review and a numbered 'Required fixes'
   list to .claude/reports/<YYYYMMDD-HHMM>-review-<task>.md; return at most 250 words."
3. When it returns, show the user the verdict, the report path and the numbered fixes only.
4. If there are required fixes, tell the main session to relay them to the implementer by path, not
   by pasting: "Apply the required fixes in <report path>, follow study-conventions, re-run the
   affected tests, and add the new commit hashes to your report."
