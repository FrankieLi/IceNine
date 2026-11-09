# Autoresearch stage: saved state

Read this file first after any restart. It says where the work stands and what to do next. The plan
is `docs/autoresearch/PLAN.md`. The mandate is D-20261109-N1 in `docs/DECISIONS.md`.

## How to restart

1. In the repository root, check `git status` and the branch (`feature/eval-reduction`
   or one of its task branches).
2. Read "Next action" below. For more context: `shmem query "IceNine autoresearch"` (`shmem ask`
   was broken on 2026-11-09; try it, fall back to `query`).
3. Before a new task: `scripts/dev/read_on_fetch.sh --save`, merge
   `origin/feature/python-speed-parity`, rerun with `--mark`, read Ilium's mailbox in
   `docs/ROADMAP.md`.
4. Do one step. Then update this file, commit, push, and save one shmem fact
   (`shmem add --source llm "IceNine autoresearch <step>: <result>"`).

Rules that hold for every step: Python only through `uv run`; one fresh implementer agent per
task in its own worktree (brief from `scripts/dev/make_brief.sh --plan docs/autoresearch/PLAN.md
--task <task> --section "<heading>"`), `scripts/dev/verify_handback.sh`, then `/review-task`, fixes
applied, merge with `scripts/dev/finish_task.sh`; at most 6 heavy processes; no merge to `develop`
without the owner's go; no change to `.claude/` files; no p-values; no extra seeds; findings in
plain English with the plot before the table.

## The wake-up job

A session cron job (id 6f6deed3, created 2026-11-09, every 3 hours at :23) sends the wake-up
prompt while the session is open and idle. It expires on 2026-11-16: create it again on or after
2026-11-15 and write the new date here. It dies with the session. After a closed session the owner
types `/resume` or "continue the autoresearch"; this file is enough to go on.

## Status (2026-11-09)

- Review page: 24 items decided (D-20261109-N2 to N4). Table on
  `docs/owner_review_2026-10-30.md`, "Decisions of 2026-11-09".
- Plan approved under the mandate (D-20261109-N5).
- Three agents were started on 2026-11-09 about 14:40 (see "Running jobs").

## Next action

1. Item 6: read `.claude/reports/20261109-review-review-fixes.md`. If it lists required fixes, send
   them to a fresh Opus implementer by path. Then merge the branch (`scripts/dev/finish_task.sh`
   from its worktree `.claude/worktrees/agent-a26f49085a74cdfcd`) and push.
2. AR0 and AR1: when a report exists (`.claude/reports/20261109-1439-ar0.md`, `-ar1.md`, in the
   agent's worktree under `.claude/worktrees/`), run `scripts/dev/verify_handback.sh`, then
   `/review-task`, apply fixes, merge, and write the result under "Results so far", in
   MIGRATION_HISTORY and in the README (new module `icenine/ray_indexing.py`).
3. If an agent died with the session (a task branch exists but no report): look at
   `git log feature/eval-reduction..feature/eval-reduction-ar0` (or `-ar1`) in its worktree, and
   start a fresh implementer with the same brief (regenerate it with `make_brief.sh`) and the note
   "continue from the commits on the branch".
4. After AR1: decide the gate (recall at 20 of at least 90% on the pilot's full-search voxels),
   then start AR2 and AR3, or AR4 first if the voting misses too much.

## Queue

Reading the table: one row per step, in order, with its state.

| step | what | state |
|---|---|---|
| item 6 | review and merge the review fixes | review running (Opus reviewer) |
| AR0 | stop rules from the stored results | implementer running, branch `feature/eval-reduction-ar0` |
| AR1 | recall of the ray-indexing proposer | implementer running, branch `feature/eval-reduction-ar1` |
| item 7 | open the PR to `develop` (merge waits for the owner) | after item 6 |
| item 14 + 3 | both keys on E0 and the pilot, timing A/B | not started; needs a quiet machine |
| item 18 | does C++'s clang fuse the clip's multiply-adds | not started |
| C++ time | measured C++ time on the pilot region (owner's question on finding 6) | not started |
| AR2 | indexing as the BFS seed search | waits for AR1 |
| AR3 | index the layer and paint grains | waits for AR1 |
| AR4 | the learned proposer | waits for AR1 |
| AR5 | the stop rule | waits for AR0 and AR1 |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |

## Running jobs

Agents of the session of 2026-11-09 (they die with the session; `git worktree list` shows their
worktrees, and their branches keep the commits): the reviewer of the review-fix branch, the AR0
implementer and the AR1 implementer. No detached compute job is running.

## Results so far

None.

## For the owner

- Still yours: the merge to `develop` (item 7), items 12, 22 and 23, and Ilium's plan approvals.
- `shmem ask` fails here: the API key in `~/.shmem/.env` is refused and the local embedder is
  not running.
