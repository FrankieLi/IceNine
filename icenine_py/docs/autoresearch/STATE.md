# Autoresearch stage: saved state

Read this file first after any restart or wake. It holds only what is true now and stays under
120 lines: rewrite it, do not append. Then read "Keeping the main context small" in
`docs/autoresearch/RUNBOOK.md`, and after that only the runbook section for the step at hand.

| what | where |
|---|---|
| mandate | D-20261109-N1 in `docs/DECISIONS.md` |
| plan | `docs/autoresearch/PLAN.md` |
| how to run a step, limits, tools | `docs/autoresearch/RUNBOOK.md` |
| finished steps, old notes (archive, use `grep`) | `docs/autoresearch/STATE_HISTORY.md` |
| numbers of a finished step | its findings page, or `shmem query "IceNine autoresearch <step>"` |
| review page for the owner | https://claude.ai/artifact/A63khAHf7nEx6NBZioSnGW |
| open PR | #43 (AR3 steps 2 to 7), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below. Do every step that is ready, one after another, and independent
   jobs side by side (owner, 2026-11-10). After each step: rewrite this file, commit, push, one
   shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action

Step 8 of AR3 (why the seed step costs so much in grains with strong gradients, task `ar3h`;
plan in `PLAN.md` under "Step 8: ...", approved by D-20261111-N4). An implementer is running.

1. Wait for its notification. Then: `verify_handback.sh feature/eval-reduction <report>
   --branch feature/eval-reduction-ar3h`, `/review-task`, the decision on the fixes
   (D-20261111-N5 or the day's next ID), a fixes round, `finish_task.sh --yes --push` in the
   worktree that holds the branch, and the close (runbook, "A step, from plan to merge",
   points 5 to 9). Copy the fixes report out of the fixes worktree before the merge: after
   the merge of step 7 that worktree was removed, and its untracked report with it.
2. At the close, add the step's whole-layer runs (run `S8`, if Part 3 ran) to
   `docs/autoresearch/pareto/pareto_rows.json`, and check the plan's three predictions on the
   review page.
3. Then: opt-in keys for the seed step and the boundary pass (library change, tests first;
   both numbers as settings; the recommended seed step is the owner's choice). Then AR4, AR5.

## Running jobs

- Step 8 implementer (Sonnet, own worktree, branch `feature/eval-reduction-ar3h`), started
  2026-11-11. Its report: `.claude/reports/20261111-0423-ar3h.md` (maybe in its own worktree;
  the findings page may sit next to it if its Write was refused).
- The review page is at version 34 (step 7 done, step 8 running, decision N4). No page job
  is open.

## Queue

Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0, AR1, AR2, AR3 steps 1 to 7) are in `STATE_HISTORY.md` under "Queue".

| step | what | state |
|---|---|---|
| AR3 step 8 | where the seed step's evaluations go; remedies (task `ar3h`) | implementer running |
| AR3 later | opt-in keys for the seed step and the boundary pass | with settings; the owner picks |
| AR4 | the learned proposer | after step 8 |
| AR5 | the stop rule inside the grid search | after step 8 |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B | needs a quiet machine |
| C++ time | measured C++ time on the pilot region | not started |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py` misses numbers with thousands separators | small, not started |

## Last finished step

AR3 step 7 (larger grains, task `ar3g`): merged into `feature/eval-reduction` at b98a14b on
2026-11-11, pushed; suite 1517 passed, 47 skipped, 1 deselected. Page: `docs/autoresearch/findings_ar3_big.md`. Decisions
D-20261111-N2 (plan), N3 (review: 10 fixes, no rerun). Sample `big`: 32 grains on the
24,570-voxel grid, gradients 0 to 4 degrees. Run plus pass: `D` 65.87 million evaluations, 26
wrong; `D4` 11.94 million, 21; `S5` 13.00 million, 30; `P4R4` 19.96 million, 31. Two of the
plan's three predictions do not agree: the seed step stays the larger part (57.5% of `S5`),
and `S5` costs 8.8% more than `D4` here (it saved 47% on the first gradient sample). The
plan's rule does not pick a seed step. Review report:
`.claude/reports/20261111-0500-review-ar3g.md`. More in `STATE_HISTORY.md`, "AR3 step 7".

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-a38ba6f5a25aa033e` (step 7, detached at b762441): its caches and runs under
  `ar3/big/`, and the report `.claude/reports/20261111-0151-ar3g.md`.
- `agent-a6db3cabf8ab230e9` (step 6, detached at 6c4665c): `seed/cache/`, `seed/runs/`, a
  `boundary/cache/`, the image-stack link, and the report `.claude/reports/20261110-2308-ar3f.md`.
- `agent-a27b4a932b0dd99b6`: step 6's fixes report `.claude/reports/20261111-0128-ar3f-fixes.md`.
- `agent-a4fd6282107a76186` (detached at 3813d69): step 5's `boundary/cache/` (the maps after
  the pass, `after_<map>.npy`, and the run environments `env_grad`, `env_layer`), the link to
  the image stacks `Examples/Example2.ManyGrains/ScatteringData_PhaseD`.
- `agent-ad88d4ade1ca5a039` (step 5 fixes): a second copy of `boundary/cache/`, and the fixes
  report `.claude/reports/20261110-2233-ar3e-fixes.md`.
- Steps 3 and 4: `agent-ae68d32fb972a8ad4` and `agent-ac609f1c28c619ed0` (see
  `grep -n "worktrees/agent-" docs/autoresearch/STATE_HISTORY.md`).

## For the owner (open items only; details by `grep` in `STATE_HISTORY.md`, "For the owner")

- Your word on PR #43 ("merge PR 43" in the session). It now holds AR3 steps 2 to 7.
- A choice: which seed step the opt-in key recommends. Step 7 does not pick one: step 6's
  (`S5`) saved 47% on the first gradient sample and cost 8.8% more on larger grains. My plan:
  the key takes both numbers as settings, so both are reachable.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (the commands are in the
  history file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Questions that block nothing: the shape and size of real orientation gradients; your
  local-optimization settings on real data; how many voxels a grain spans in your maps.
- A note to read: C++ parallel mode keeps the higher confidence for a voxel that two clients
  fit; that rule did harm in the Python test of step 5. A reading of the source, not a test.
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 to N4 (all
  on the review page). The gradient levels of step 7's sample are my choice. Step 5's result
  files add about 50 MB to git; say so if you want them pruned.
- Not done unless you ask: a repo-wide rename of "finisher" to "local optimization".
- Ilium's plan approvals, when posted, are yours.

## Bookkeeping

- Next decision ID: D-20261111-N5.
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`.
