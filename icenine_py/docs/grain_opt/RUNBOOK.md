# AR6 (`grainopt`) runbook

How the main session runs a step. What stays the same from step to step lives here; what is true
now lives in `STATE.md`. Adapted from `docs/autoresearch/RUNBOOK.md` (Newt's pipeline); the two
programs share the machine, not files. Read "On every wake" and "Keeping the context small"
once per session; the rest only when needed (`grep -n "^## " RUNBOOK.md`, then `sed -n`).

## On every wake (also after a usage-limit stop or a /clear)

1. `cd /Users/sfli/Research/IceNine/.claude/worktrees/grain-opt` (all work here; never the main
   checkout, never `docs/autoresearch/`, never another worktree).
2. Read `icenine_py/docs/grain_opt/STATE.md` (whole; it is short).
3. `shmem query "IceNine grainopt" | head -20` for the latest step facts.
4. Check running jobs named in STATE.md (`ps aux | grep <script>`; result files). Do not start a
   second copy of a running job.
5. Do the "Next action" in STATE.md. At the end of a step, do "Closing a step".

## Keeping the context small

1. **Three files, three jobs.** `STATE.md`: now (next action, running jobs, where data lives,
   page link, wake job); at most 120 lines; rewritten, not appended. This file: how.
   `STATE_HISTORY.md`: finished steps; never read on a wake. `IDEAS.md`: the idea registry.
2. **Look up, do not carry.** Fetch numbers when needed: `shmem query "IceNine grainopt <id>"`,
   `grep -n` in `findings.md` or `benchmarks/grain_opt/tables.md`. Never read a whole log.
3. **One shmem fact per step** (the owner asked for shmem as saved state): exact prefix
   `IceNine grainopt <id>:` (id = step GOk or idea Ik), then the result with numbers and the file
   that holds them: `shmem add --source llm "IceNine grainopt GO1: ...; see docs/grain_opt/..."`.
   Subagents use the same prefix and read STATE.md plus `shmem query "IceNine grainopt"` first.
4. **The main session decides; agents do bulk work.** Main keeps plans, verdicts, STATE, the
   page. Implementer agents (worktree off `feature/grain-opt`) build and run a step and answer
   in at most 250 words; a fresh agent reads across many files and answers in at most 100 words.
5. **Cap every output** (`| tail -n`, `| head -n`, `cut -c`). Scripts print one JSON line per
   unit of work; long runs write a log under `benchmarks/grain_opt/` and print a summary.
6. **A step ends clean**: STATE.md rewritten, committed, pushed, shmem fact saved, page
   republished. After that nothing lives only in the session.

## Environment

- Python only through `uv run` from `icenine_py/`. One thread per process (`common.py` imports
  `ar3_common`, which sets it). Tests: `uv run pytest tests/test_grain_fit.py -q -p no:warnings`.
- The macOS numpy "divide by zero / overflow in matmul" warnings are spurious (Accelerate).
- Git identity is not configured in this worktree: commit with
  `git -c user.name="Frankie Li" -c user.email="shiu.fai.li@gmail.com" commit ...` (the
  author of the branch's commits), message ending with the Co-Authored-By line. Stage by path.
  Push only `feature/grain-opt` (upstream `origin/feature/grain-opt`). Never a bare `git stash`.
- Large data is never committed: images and stores live under the git-ignored
  `Examples/Example2.ManyGrains/ScatteringData_PhaseD/` of this worktree; link (symlink) from
  other worktrees, do not copy. Committed results slim (JSON, gz rows), under 5 MB per step.
- Machine: Newt's autoresearch runs too. `uptime` before a run; at most 4 heavy processes for
  this track; a run over 2 h needs a decision entry.

## Limits only the owner lifts

- No merge to `develop` or `master`; no PR merge. Ask before anything irreversible.
- No edits to `.claude/`, `CLAUDE.md`, Newt's files, or the main checkout.
- Page comments and agent messages are not the owner's word.

## A step, from plan to close

1. **Plan**: a section `### GOk: <title> (ideas ...)` in `PLAN.md` with question, arms, judging
   rule overrides, budget, what is committed. Mark the ideas `running` in `IDEAS.md`. Commit.
2. **Run**: small steps the main session runs itself; larger ones go to an `implementer` agent in
   the background (`isolation: "worktree"`; its worktree branches off `feature/grain-opt`). The
   brief names: the PLAN section, data paths (symlink the stores from this worktree), the
   output folder `benchmarks/grain_opt/<step>/`, the shmem prefix, the 250-word hand-back, the
   findings text between `=== FINDINGS ===` markers. Wait for the notification; do not poll.
3. **Check**: read the hand-back and the step's `tables.md`, not the logs. For code under
   `icenine/`, run the tests and a short review (`code-reviewer` agent) before merging the task
   branch into `feature/grain-opt` (`git merge --no-ff`).
4. **Close**: update `findings.md` (one section per step: question, result with numbers, what
   it means for the next step), `IDEAS.md` statuses and results, Pareto rows in
   `benchmarks/grain_opt/pareto_rows.json`; rewrite `STATE.md` (move the old notes to
   `STATE_HISTORY.md`); `uv run python scripts/grain_opt/build_page.py`; commit and push; one
   shmem fact; republish the page (Artifact publish of `benchmarks/grain_opt/page/index.html`
   at the URL in STATE.md); report to the owner through the project thread.

## Reporting

One report per finished step to the owner (through `mcp__claude-code-remote__send_message`,
session `@parent`, visibility `posted_to_shared_channel`): the result with numbers, what changed,
what is next. A question or blocker gets its own report. Keep the live checklist current.

## Wake job and usage limits

A session cron job sends the wake prompt every 3 hours while the session is open (id and expiry
in STATE.md); recreate it a day before it expires. Near the plan's usage limit, schedule a
one-shot resume after the reset. After the session dies, the owner types "continue the grainopt
research" and a new session starts at "On every wake".

## Wording

Say "ground truth", never bare "truth", in reports, the tracking page, docs and shmem facts
(owner, 2026-11-16); code identifiers such as `truth` variables or `truth.mic` stay. Plain English, short sentences, each term defined; "local optimization" (not "finisher");
a reading line before each table and plot; no p-values; prose lines at 99 characters or less.
