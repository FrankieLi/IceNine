---
title: "TODO: Ilium Agent workflow tooling: scripts and thin skills"
subtitle: "Turn the recurring manual steps of the Ilium Agent workflow into deterministic scripts"
date: "2026-10-16"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Approved 2026-10-18, in progress** (D-20261018-I2; owner answers to O1-O9 in D-20261018-I1 to
I5). Written 2026-10-16 by Ilium Agent after the S4 session, at the owner's request ("which parts
could become skills or deterministic scripts"). Owner answers: O1 keep the written ID rule (n from 1
per day and agent); O2 approved with the three-task grouping; O3 H12 row placed by Ilium Agent; O4
at most 6 heavy processes on this machine, exceptions allowed, kept per machine
(`icenine.maxHeavyProcs`); O5-O9 as proposed below.

# What

Five scripts and four thin skills, in priority order:

| # | Tool | Replaces |
|---|---|---|
| 1 | `scripts/dev/log_decision.sh` + `/decide` | hand-written DECISIONS.md entries, shmem line, commit, push |
| 2 | `scripts/dev/roadmap_row.py` | hand (or `sed`) edits of a Track S status cell |
| 3 | `scripts/dev/verify_handback.sh` | the main session's manual check of an implementer hand-back |
| 4 | `scripts/dev/read_on_fetch.sh` | protocol rule 2 by hand: fetch, log, mailbox, Newt's decisions |
| 5 | `scripts/dev/make_brief.sh` + `/implement-plan`, `/apply-fixes` | implementer briefs typed per task |

Lower priority: `sync_doc_tables.py --section`, a timed-run wrapper that writes the preflight label
into the result JSON, a shared study-helper module from the S4 scripts, and a mailbox writer.

The scripts do the deterministic work and fail loudly; the skills only collect the arguments from
the conversation, call the script, and show the result. No skill re-implements a script step.

# Why

In the S4 session the main session did these steps by hand:

- four DECISIONS.md entries (I19 to I22), each with an ID to work out, a template to copy, a shmem
  line, a commit and a push;
- ROADMAP row edits, one of which (a `sed` on the S4 status cell) matched nothing and failed
  silently;
- implementer briefs that were mostly boilerplate repeated from study-conventions and the plan;
- hand-back checks (commits since the parent, files touched, report sections, the tests it names).

Each of these is mechanical, repeats every task, and costs context. One of them already went wrong
silently. Two more problems found while writing this plan show the cost of doing it by hand:

- **Decision IDs.** DECISIONS.md says n counts "from 1 per day and agent". The entries up to
  D-20261014-I21 follow that (I1 to I21 are all on 2026-10-14), but the 2026-10-16 entry is
  `D-20261016-I22`, a running counter, not `I1`. The log is append-only, so I22 stays; the script
  needs one rule (owner decision O1).
- **Track S table.** The S3 row of the Track S table has no closing `|`. A row editor must accept
  both forms; a naive split on `|` would drop or misplace that cell.

# Existing tooling (not duplicated)

Read on `feature/python-speed-parity` at bbce6bb:

- `scripts/dev/`: `start_task.sh` (task branch, records `branch.<b>.parent`), `finish_task.sh`
  (tests, `--no-ff` merge in the parent's worktree, push, cleanup; has the `wt_of` helper that finds
  the worktree holding a branch), `checkpoint.sh`, `job_status.sh`, `new_todo.sh`,
  `timing_preflight.py`, `sync_doc_tables.py`, `audit_numbers.py`, `precommit_check.py`,
  `install_hooks.sh`.
- `.claude/skills/`: `start-task`, `finish-task`, `merge-task`, `review-task`, `handoff`, `resume`,
  `start-feature`, `finish-feature`, `mic-maps`.
- `.claude/hooks/check-bash-command.sh` (denies bare python/pytest/pip, `git add -A/./-u`,
  `git commit -a/--no-verify`, `ALLOW_*`), `run-tests-on-change.sh`.
- `.claude/agents/`: `implementer.md` (report format with eight fields), `code-reviewer.md`,
  `study-conventions.md` (uv, explicit staging, `.claude/**` left uncommitted for the owner,
  preflight and labels, at most 10 workers, report path and 250-word hand-back).
- `icenine_py/README.md` "Developer tooling" table; tests in `icenine_py/tests/test_dev_tooling.py`
  (a temporary-repository fixture `repo`, and env overrides such as `FINISH_TASK_TEST_CMD`).

What none of them does: write a decision, edit a ROADMAP row, check a hand-back, read the other
agent's branch, or assemble a brief. Each new script reuses what exists: `log_decision.sh` reuses
`finish_task.sh`'s worktree lookup and parent rules; `make_brief.sh` points at study-conventions
instead of repeating it; `verify_handback.sh` checks the report fields defined in `implementer.md`.

# Steps

## 1. `scripts/dev/log_decision.sh` and `/decide`

**Interface.**

```
log_decision.sh --row "<track/row>" --title "<one-line decision>" --decision "<text>"
                --where "<ROADMAP row, plan or doc>" [--supersedes <ID>] [--date YYYY-MM-DD]
                [--agent I|N] [--recorder "Ilium Agent"] [--no-push] [--dry-run]
```

**Behaviour, in order.**

1. Find the target branch: the current branch's recorded parent chain (`git config
   branch.<b>.parent`), walked up to the first `feature/*` branch with no recorded parent; if the
   current branch is that branch, it is the target. Refuse `develop`, `master`, `main`, a detached
   HEAD, and a chain with no recorded parent (no guessing, as `finish_task.sh --yes`).
2. Take a lock (`mkdir "$(git rev-parse --git-common-dir)/log_decision.lock"`) so two sessions on
   this clone cannot pick the same ID; release it on exit (`trap`).
3. `git fetch origin <target>`; refuse if `origin/<target>` has commits the local target lacks
   (print the merge command; the caller merges, then reruns).
4. Read DECISIONS.md **from the target branch tip** (`git show <target>:<path>`), not the
   working tree, and compute the ID by the rule the owner picks (O1). Default, the written rule:
   n = 1 + the largest n of `### D-<YYYYMMDD>-<agent><n>` for that date and agent, else 1.
5. Validate: no newline in `--title`; `--supersedes` names an existing heading; no existing entry
   on that date has the same title (refuse a duplicate run).
6. Build the entry exactly in the file's format, wrapped at 98 columns with a two-space
   hanging indent (the docs keep lines to 100):

   ```
   ### <ID> · <track/row> · <title>
   - <date>, owner. Recorded by <recorder>. Supersedes: <ID or none>.
   - <decision>. Where: <where>.
   ```

   `--dry-run` prints the target branch, the ID, the entry and the shmem line, and stops here.
7. Commit on the target branch without switching the caller's checkout:
   - if a worktree already holds the target (the `wt_of` lookup of `finish_task.sh`), append and
     commit there, refusing if that worktree has uncommitted changes to DECISIONS.md;
   - otherwise `git worktree add` a temporary worktree of the target under `$TMPDIR`, append and
     commit there, and remove it on exit.

   Stage the one file by path; message `DECISIONS <ID>: <title>` plus the Co-Authored-By trailer
   (taken from `--trailer` or the `LOG_DECISION_TRAILER` environment variable; the `/decide` skill
   passes the one in the current system reminder). The git pre-commit hook runs as usual.
   Rejected alternative: building the commit with `git commit-tree` and a temporary index. It needs
   no worktree but skips the pre-commit hook.
8. Run `shmem add --source user "Owner decision <ID> (<date>, IceNine <track/row>): <title>.
   <decision>"` (command overridable with `LOG_DECISION_SHMEM_CMD` for tests).
9. `git push origin <target>` unless `--no-push` (protocol rule 3(d)).
10. Print the ID, commit hash, and one line per step (`ok` / `FAILED: <command to rerun>`). Exit 0
    only if all steps passed; a failed shmem or push leaves the commit in place and prints the exact
    command to rerun.

The task branch gets the entry when it is merged into the target (union merge, no conflict). The
caller never commits DECISIONS.md on a task branch; `verify_handback.sh` flags it if one does.

**`/decide` skill** (`.claude/skills/decide/SKILL.md`, about 20 lines, `allowed-tools:
Bash(scripts/dev/log_decision.sh *), Bash(git *)`): from the owner's words, fill row, title,
decision, where and supersedes; run with `--dry-run` and show the entry; run for real; report the ID
and commit. If the decision changes a ROADMAP status, call `roadmap_row.py` next (step 2) and say
so. The skill never words a decision more strongly than the owner did.

**Tests** (`tests/test_dev_tooling.py`, on the `repo` fixture plus a bare repository as `origin`):
ID on a new date, on a date with entries, ignoring the other agent's entries; entry format and
wrapping; commit lands on the parent from a task branch, both with the parent in another worktree
and with none; the caller's branch and working tree unchanged; refusals (develop, no parent,
duplicate, unknown supersedes, origin ahead); shmem stub called with the exact line; `--no-push`;
lock released after a failure.

## 2. `scripts/dev/roadmap_row.py`

**Interface.**

```
roadmap_row.py --row S4 (--status "<text>" | --state "<bold word(s)>" | --show) [--check]
               [--doc icenine_py/docs/ROADMAP.md]
```

- Finds the `# Track S:` section (up to the next `# ` heading) and the table rows in it whose first
  cell is exactly the row ID. Exactly one must match; zero or several exit 2 with a message naming
  the section and the IDs it does have.
- Refuses a row outside Ilium Agent's tracks (protocol rule 4): the allowed set is a constant,
  `{"S"}` (O7 decides whether Q or H rows the owner gave to Ilium are added).
- `--status` replaces the whole last cell; `--state` replaces only the first `**...**` span in it
  (for example `**in progress**` to `**parked**`) and keeps the detail; it exits 2 if the cell has
  no bold span.
- Handles rows with and without a closing `|` (the S3 row has none) and writes the row back in the
  form it had. Refuses text that contains a newline or an unescaped `|`.
- Prints the old and new row (`-`/`+`). `--check` writes nothing and exits 1 if the file would
  change. It never commits; the caller commits the row with the change it records.

**Tests:** replace, `--state` keep-detail, missing row, duplicate row, wrong track, row without a
closing pipe, pipe in the text, idempotent second run, `--show`.

## 3. `scripts/dev/verify_handback.sh`

**Interface.** `verify_handback.sh <parent> <report> [--branch B] [--fixes REVIEW_REPORT]
[--no-tests]`. `--branch` defaults to the current branch.

**Checks**, each printed as `ok`, `FLAG` (the main session decides) or `FAIL` (exit 1):

- Commits: `git log --oneline <parent>..<branch>`; FAIL if none. Every commit hash on the branch is
  named in the report and every hash the report names is on the branch (FLAG either way).
- Working tree of the branch's worktree: uncommitted changes outside `CLAUDE.md` and `.claude/**`
  are a FLAG (the `.claude/**` files are expected to be left for the owner and are listed).
- Changed files (`git diff --name-only <parent>...<branch>`), each class listed separately:
  - the shared files of "Coordination between Tracks E and S" (`icenine/cost_functions.py`,
    `orientation_search.py`, `batched_screen.py`, `reconstructor.py`, `_rasterize.c`): FLAG, with
    the reminder that a heads-up is due (protocol rule 3(a));
  - any other `icenine_py/icenine/**`: FLAG (library change; the plan must allow it);
  - `Src/**` (C++): FLAG;
  - `pyproject.toml`, `uv.lock`: FLAG (new dependency);
  - `icenine_py/tests/data/*golden*`: FLAG (no golden regeneration inside a task);
  - `icenine_py/docs/DECISIONS.md`, `CLAUDE.md`, `.claude/**` committed: FAIL;
  - `icenine_py/docs/ROADMAP.md`: FLAG if a line outside Track S or "From Ilium Agent" changed.
- Report: the eight fields of `implementer.md` ("Branch / commits", "Files", "Tests", "Results",
  "Criteria", "Deviations from the plan", "Suspicious / unexpected", "Left undone") each present
  (FAIL if one is missing); length over about 80 lines is a FLAG.
- Timing labels: a Results line with a time unit (`µs`, `us`, `ms`, ` s`) and none of
  `single-worker`, `contended`, `quiet` is a FLAG (a heuristic; the reviewer still checks).
- With `--fixes`: every number of the review's "Required fixes" list appears in the hand-back
  report (FLAG per missing number).
- Tests: the `tests/...py[::name]` paths named on the report's "Tests" line, run in the branch's
  worktree as `cd icenine_py && uv run pytest <paths> -q` (output through `tail -20`); FAIL on a
  failure; FLAG if the report names none. `--no-tests` skips this.

The output is at most about 30 lines. It is a gate before `/review-task`, not a review.

**Tests:** a fixture branch with a shared-file change, a missing report field, a committed
`.claude` file, a passing and a failing named test (via a `VERIFY_HANDBACK_PYTEST_CMD` override),
hash cross-check both ways, `--fixes` numbering.

## 4. `scripts/dev/read_on_fetch.sh`

**Interface.** `read_on_fetch.sh [--other feature/eval-reduction] [--other-agent N] [--since REF]
[--mark] [--save]`.

1. `git fetch origin`.
2. Last-read point: `git config --local icenine.readOnFetch.<other>` (shared by all worktrees of the
   clone), else `--since`, else the merge base with the current branch. Print
   `git log --oneline <last>..origin/<other>` and its `--stat` summary (at most 40 lines).
3. Mailbox: extract the "From Newt Agent" subsection of ROADMAP.md at `<last>` and at
   `origin/<other>` (`git show`), print a unified diff of the two.
4. Decisions: list the `### D-*-N*` entries of `origin/<other>:icenine_py/docs/DECISIONS.md` whose
   ID does not occur in `shmem query "<ID>"`, and print for each a ready line
   `shmem add --source user "Owner decision <ID> (<date>, IceNine <track/row>): <title>. <text>"`.
   `--save` runs them instead of printing (rule 6 asks for it; off by default until O8).
5. Print, not run, the merge command for protocol rule 2 (`git merge origin/<other>`): a conflict
   outside own rows needs the owner, so the merge stays a visible step.
6. `--mark` stores `origin/<other>`'s hash as the new last-read point; without it the point does not
   move, so an interrupted session reads the same entries again.

**Tests:** a fixture with a bare origin and an "other" branch that adds a commit, a mailbox line and
an N decision; shmem stub that knows one of two IDs; `--mark` moves the point and a second run is
empty.

## 5. `scripts/dev/make_brief.sh`, `/implement-plan`, `/apply-fixes`

**`make_brief.sh`** prints the brief; the skills pass it to the `implementer` agent unchanged, so
the brief is deterministic and testable.

```
make_brief.sh --plan docs/todo_<slug>.md --task <task> [--section "<heading>"]
              [--fixes <review report>] [--trailer "<Co-Authored-By line>"]
```

Contents, in this order (about 20 lines):

- "Follow study-conventions." (uv only, explicit-path staging since the hook blocks `git add -A`,
  preflight and timing labels, report format are there and are not repeated).
- Plan path and the section to implement; with `--fixes`, "Apply the required fixes in <report>"
  (the sentence `/review-task` already uses) and "name the fix numbers in your commit subjects".
- Branch: `scripts/dev/start_task.sh <task>` from `<parent>` (or, for fixes, the existing task
  branch or a `-fixes` branch when the task was already merged, as Track E did for Task C).
- Scope, taken from the plan: library changes allowed or "no change under `icenine_py/icenine/`";
  "no new dependency (`pyproject.toml`, `uv.lock` unchanged)"; no default flip, no golden change.
- The plan's stop rule, quoted from its `**Stop rule:**` lines (`grep`); "stop and report" if the
  plan has none.
- Compute: "timings run alone; at most 6 heavy processes at once on this 16 GB machine" (see O4).
- The commit trailer from `--trailer` (the skill passes the one in the current system reminder).
- Report path, computed now: `.claude/reports/<YYYYMMDD-HHMM>-<task>[-fixes].md`.

It refuses a plan with no DECISIONS.md entry naming it as approved (`grep` for the plan file name
and "approved"), and prints the entry it found.

**`/implement-plan <plan> [section] [model] [effort]`** (`.claude/skills/implement-plan/SKILL.md`,
modelled on `review-task`): run `make_brief.sh`; launch `implementer` in the background with the
model and effort the plan or its approval decision names (default: the agent's own); when it
returns, run `verify_handback.sh <parent> <report>` and show its summary with the report path; then
suggest `/review-task`.

**`/apply-fixes <review report>`**: run `make_brief.sh --fixes <review report>` for the task
branch, launch a **fresh** implementer (Track S rule: no owner review in between), then
`verify_handback.sh <parent> <fix report> --fixes <review report>`.

**Tests:** `make_brief.sh` on fixture plans with and without a stop rule, with and without an
approval entry, with `--fixes`; golden text comparison of the brief.

## Lower priority

- **`sync_doc_tables.py --section "<heading>"`**: sync or `--check` only the markers inside that
  heading's range, so a task that owns one section does not need generated blocks for every other
  section of the doc (today a doc marker with no generated block is exit 2).
- **Timed-run wrapper** `scripts/dev/timed_run.py --out RESULT.json [--workers N] [--require] --
  <command>`: run the preflight (`preflight.py`), refuse or record a busy machine, run the command,
  and write `preflight` and the labels into `RESULT.json` (or a `RESULT.preflight.json` sidecar for
  non-JSON output). Replaces S4's `common.preflight_json` and `label_of`. The labels need one
  vocabulary first (O5): study-conventions says "single-worker" or "contended", S4 wrote "quiet" or
  "contended"; proposed: two fields, `workers` and `machine` (`quiet` / `contended`).
- **Shared study module** (for example `icenine_py/scripts/common/study_data.py`) from
  `scripts/speed_parity/s4/`: `DataSet` and `load_set`, `capture_level0`, `reference_voxels`,
  `timed_voxels`, `screen_reference` (all in `common.py`); the flip listing (`trunc_c`, `flip_dist`
  in `translation_correctness.py`, the flip explanation in `recost.py`); the linearity check (in
  `orderings.py cmd_time`). The S4 scripts would import from it so their results still reproduce.
  This waits for S4: its branch is unmerged (D-20261016-I22).
- **Mailbox writer**: `roadmap_row.py mailbox add "<tag>" "<text>"` and `mailbox rm <date>
  <prefix>`, editing only "From Ilium Agent" (dated line, tag checked against heads-up / request /
  answer).

# Rules for the implementer

- Task branches off the branch O6 picks: one task per group, reviewed by `/review-task` as any Track
  S task. Suggested grouping: (a) items 1 and 2; (b) items 3 and 4; (c) item 5. A Sonnet
  implementer (the agent default) suffices: shell and stdlib Python with tests.
- Shell in the style of `finish_task.sh` (`set -euo pipefail`, usage line, `error:` to stderr,
  override env variables for tests); Python stdlib only, Black at 100.
- Scripts stage by explicit path and never pass `--no-verify`.
- `.claude/skills/*` files are written but left uncommitted (study-conventions; the pre-commit check
  needs `ALLOW_CLAUDE_CONFIG=1`, which only the owner sets). The report lists them for the owner.
- Every tool gets a row in the README "Developer tooling" table and tests in
  `tests/test_dev_tooling.py`; no network (a bare repository stands in for `origin`), no real
  shmem writes in tests.
- No change to `icenine_py/icenine/`, `Src/`, `pyproject.toml` or `uv.lock`.

# Instruction changes and sharing with Newt Agent

The owner wants the tools and the instructions to use them shared with Newt Agent (2026-10-16). Each
machine has its own memory, shmem store and checkpoints, so nothing in them reaches the other agent.
Every instruction therefore goes into a tracked file, and Newt gets it by merging
`feature/python-speed-parity` (protocol rule 2). Per-machine memory holds only a pointer.

| Instruction | Where it goes | Who edits |
|---|---|---|
| Record owner decisions with `log_decision.sh` (or `/decide`), never by hand | `study-conventions.md`; ROADMAP rule 6, one sentence | owner (`.claude/`, shared protocol section) |
| Do rule 2 with `read_on_fetch.sh`, then merge | ROADMAP rule 2, one sentence | owner (shared protocol section) |
| Edit a track's status cell with `roadmap_row.py`, not `sed` | `study-conventions.md` | owner |
| Main session: `verify_handback.sh` before `/review-task`; briefs via `/implement-plan`, `/apply-fixes` | `study-conventions.md`, the skills | owner |
| Process limit and other machine limits | per machine (see below), not shared | each agent |
| Usage and options of each tool | README "Developer tooling" table, script usage lines | implementer |
| ID rule, if O1 changes it | a DECISIONS entry, then the DECISIONS.md header | owner |

- **Agent-neutral tools.** The agent letter (`I` or `N`), the other agent's branch and the mailbox
  names come from `git config --local icenine.agent` (set once per clone: `I` here, `N` on Newt's
  machine) or from `--agent` / `--other`. A script refuses to run when neither is set; it never
  defaults to Ilium. Tests run both settings.
- **Machine limits stay per machine.** "At most 6 heavy processes" is a fact about this 16 GB
  machine; Newt's machine may differ. Proposed: `git config --local icenine.maxHeavyProcs <n>`, read
  by `make_brief.sh`, with study-conventions saying "the per-machine limit `icenine.maxHeavyProcs`".
  This replaces the choice in O4.
- **Skills.** The skills reach Newt only once the owner commits them (`ALLOW_CLAUDE_CONFIG=1`).
- **Rollout.** A heads-up in the "From Ilium Agent" mailbox when this plan is posted (2026-10-16),
  and again when each tool group merges, naming the tools and the instruction lines added. Each
  agent then adds one pointer line to its own memory ("decisions, rows, hand-backs and read on fetch
  go through `scripts/dev/`; see study-conventions"). Newt writes its own.
- **Owner text needed (O9).** The study-conventions lines and the one-sentence additions to rules 2
  and 6 above are owner edits. The implementer of group (c) drafts them in its report; the owner
  commits them.

# Risks

- **Writing on a branch the caller has not checked out.** `log_decision.sh` commits on the feature
  branch from a task branch. Mitigations: lock, fetch-and-refuse when origin is ahead, temporary
  worktree removed on exit, refusal on a dirty DECISIONS.md in the parent's worktree, tests for both
  worktree cases.
- **shmem as an existence check** (item 4) is keyword search; a miss only causes a duplicate fact,
  which is harmless. An unreachable shmem must be reported, not read as "not saved".
- **Heuristics** (timing labels, report hashes) can flag correct reports; they are FLAG, never FAIL.
- **Skills drift from scripts.** Keep each skill to calling one script and showing its output.
- **The other agent.** A tool that silently assumes Ilium would write Newt's decisions under `I`.
  The scripts refuse to run without `icenine.agent` or `--agent` (see "Instruction changes and
  sharing with Newt Agent").

# Owner decisions needed

- **O1 Decision IDs:** keep the written rule (n from 1 per day and agent; next 2026-10-16 entry I23,
  next day I1), or change the rule to a running counter per agent (I22 then I23 on any date) with a
  new DECISIONS entry. The script follows whichever is picked.
- **O2 Approval of this plan** and of the grouping into three tasks.
- **O3 ROADMAP row.** Protocol rule 4 lets only the track's agent edit a track's rows, on the
  track's branch; Track H's branch (`feature/housekeeping`) is not started and Track H is "either
  agent". So this plan adds no row. Proposed row, for the owner to place:
  `| H12 | scripts/dev/, .claude/skills/ | workflow tooling: log_decision, roadmap_row,
  verify_handback, read_on_fetch, make_brief and the /decide, /implement-plan, /apply-fixes skills
  (docs/todo_ilium_workflow_tooling.md); plan waiting for approval |`.
- **O4 Process limit:** study-conventions says at most 10 workers; D-20261014-I20 and the S4 plan
  say at most 6 heavy processes on this 16 GB machine. Make 6 the standing rule in study-conventions
  (an owner edit of `.claude/`), so `make_brief.sh` can drop that line, or keep it per plan.
- **O5 Timing label vocabulary** for the timed-run wrapper (two fields, as proposed above, or one).
- **O6 Branch** for the tooling tasks: task branches off `feature/python-speed-parity` (usable by
  Ilium now, reaches Newt with Track S) or a new `feature/housekeeping` (Track H).
- **O7 `roadmap_row.py` tracks:** Track S only, or also Q/H rows the owner assigns to Ilium (as Q8
  was).
- **O8 `read_on_fetch.sh --save`:** may it save Newt's decisions to shmem by default (rule 6 asks
  for it), or only print the lines.
- **O9 Shared instructions:** approve where each instruction goes (table in "Instruction changes
  and sharing with Newt Agent"), the per-machine process limit in place of O4, and the owner edits
  to study-conventions and protocol rules 2 and 6.

# Relation to other work

- ROADMAP "Between agents" rules 2, 3, 4 and 6 (what items 1, 2 and 4 automate).
- `.claude/agents/implementer.md` (report fields checked by item 3) and `study-conventions.md` (what
  item 5's briefs point at instead of repeating).
- `.claude/skills/review-task` (the model for the new skills; its step 4 sentence is reused by
  `/apply-fixes`).
- `docs/todo_translation_screen_prototype.md` (S4: source of the stop rule, process limit and the
  study helpers in "Lower priority").
