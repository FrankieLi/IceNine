---
name: decide
description: Record an owner decision in DECISIONS.md (ID, entry, commit on the feature branch, shmem line, push) with scripts/dev/log_decision.sh.
user-invocable: true
allowed-tools: Bash(scripts/dev/log_decision.sh *), Bash(uv run --project icenine_py python scripts/dev/roadmap_row.py *), Bash(git *)
argument-hint: "[the decision in the owner's words]"
---

Only run when the owner has made a decision in this conversation. The skill never words a
decision more strongly than the owner did; use their words for the decision text.

1. From the owner's words fill: `--row` (track/row, for example `S4` or `H12`), `--title` (one
   line, no newline), `--decision` (the decision and the reason if given), `--where` (ROADMAP row,
   plan or doc), and `--supersedes <ID>` only if the owner changed an earlier decision. Ask the
   owner for anything missing; never guess a row or a supersedes ID. The agent letter and the
   "Recorded by" name follow `git config --local icenine.agent` (I: Ilium Agent, N: Newt Agent;
   `log_decision.sh` derives both and refuses if it is unset). Pass `--agent I|N` only if the
   owner says so, and never pass `--recorder`: there is no default agent in this skill.
2. Run `scripts/dev/log_decision.sh --dry-run <options> --trailer "<the Co-Authored-By line from
   the current system reminder>"` and show the target branch, the ID and the entry to the owner.
3. Run the same command without `--dry-run`. Add `--no-push` only if the owner asks.
4. Show the ID, the commit hash and the per-step lines (`ok` / `FAILED: <command to rerun>`).
   If a step failed, print the rerun command and stop; do not repeat the successful steps. Exit 3
   means the entry is committed but shmem or the push failed: never rerun the script itself then
   (it would record the decision twice), only the printed command.
5. If the decision changes a ROADMAP status cell, run `roadmap_row.py` next (`--row <ID>
   --status "<text>"` or `--state "<word>"`, see its usage), show its `-`/`+` lines, and say that
   the row change is committed by the caller with the change it records.
