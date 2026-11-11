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
| open PR | #43 (AR3 steps 2 to 6), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below. Do every step that is ready, one after another, and independent
   jobs side by side (owner, 2026-11-10). After each step: rewrite this file, commit, push, one
   shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action

Step 7 of AR3 (the same pipelines on larger grains, task `ar3g`; plan in `PLAN.md` under
"Step 7: the same pipelines on larger grains", approved by D-20261111-N2).

1. Wait for the implementer's notification (see "Running jobs"). Do not poll its output file.
2. Hand-back: `scripts/dev/verify_handback.sh feature/eval-reduction <its report> --branch
   feature/eval-reduction-ar3g`. Then `/review-task` (the `code-reviewer` agent, model Opus),
   a decision on the fixes (next ID below), the fixes round, `finish_task.sh --yes --push`, and
   the close (runbook, "A step, from plan to merge", points 5 to 9).
3. At the close, add step 7's runs to `docs/autoresearch/pareto/pareto_rows.json` as a new
   sample (`big`, a new panel in `pareto_build.py`), and check the plan's three predictions on
   the review page.
4. Then plan the next step from step 7's result: opt-in keys for the seed step that step 7
   supports and for the boundary pass (library change, so tests first); AR4 (learned proposer)
   for the voxels every seed step leaves to the grid search; AR5 (stop rule inside the grid
   search).

If an agent died with the session (a task branch exists but no report): look at
`git log feature/eval-reduction..<task branch>`, then start a fresh implementer with the same
brief and the note "continue from the commits on the branch".

## Running jobs

- Step 7 implementer (Sonnet, own worktree `.claude/worktrees/agent-a38ba6f5a25aa033e`, branch
  `feature/eval-reduction-ar3g`), started 2026-11-11 at about 01:50. Its report:
  `.claude/reports/20261111-0151-ar3g.md` (maybe in its own worktree). Four whole-layer runs,
  caps 180 and 90 minutes, so the hand-back can take several hours.
- The review page is at version 32 (step 6 done, step 7 running; "grid search" defined). No page job is open.

## Queue

Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0, AR1, AR2, AR3 steps 1 to 6) are in `STATE_HISTORY.md` under "Queue".

| step | what | state |
|---|---|---|
| AR3 step 7 | the runs of steps 3 to 6 on larger grains (task `ar3g`) | implementer running |
| AR3 later | opt-in keys for the seed step and the boundary pass | after step 7 |
| AR4 | the learned proposer | after step 7 |
| AR5 | the stop rule inside the grid search | after step 7 |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B | needs a quiet machine |
| C++ time | measured C++ time on the pilot region | not started |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py` misses numbers with thousands separators | small, not started |

## Last finished step

AR3 step 6 (the seed cost, task `ar3f`): merged into `feature/eval-reduction` at 73992b8 on
2026-11-11, pushed; suite 1520 passed, 44 skipped, 1 deselected. Page:
`docs/autoresearch/findings_ar3_seed.md`. Decisions D-20261110-N8 (plan), D-20261111-N1
(review: 11 fixes, no rerun). On 175 of 189 grid-fallback voxels a near-truth candidate was
costed and not refined. Run plus pass on the gradient sample: `S5` (refine the best 5 of 20)
15 fallbacks, 20.6 million evaluations, 62 wrong; `S` (neighbour tries first) 11, 16.1
million, 100 wrong; `K2` 53, 30.3 million, 96; `D4` 76, 38.8 million, 74. Review report:
`.claude/reports/20261111-0330-review-ar3f.md`. More in `STATE_HISTORY.md`, "AR3 step 6".

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-a4fd6282107a76186` (detached at 3813d69): step 5's `boundary/cache/` (the maps after
  the pass, `after_<map>.npy`, and the run environments `env_grad`, `env_layer`), the link to
  the image stacks `Examples/Example2.ManyGrains/ScatteringData_PhaseD`.
- `agent-ad88d4ade1ca5a039` (step 5 fixes): a second copy of `boundary/cache/`, and the fixes
  report `.claude/reports/20261110-2233-ar3e-fixes.md`.
- `agent-a6db3cabf8ab230e9` (step 6, detached at 6c4665c): `seed/cache/`, `seed/runs/`, a
  `boundary/cache/`, the image-stack link, and the report `.claude/reports/20261110-2308-ar3f.md`.
- `agent-a27b4a932b0dd99b6`: step 6's fixes report `.claude/reports/20261111-0128-ar3f-fixes.md`.
- `agent-a38ba6f5a25aa033e` (step 7, in use): its caches and runs under `ar3/big/`.
- Steps 3 and 4: `agent-ae68d32fb972a8ad4` and `agent-ac609f1c28c619ed0` (see
  `grep -n "worktrees/agent-" docs/autoresearch/STATE_HISTORY.md`).

## For the owner (open items only; details by `grep` in `STATE_HISTORY.md`, "For the owner")

- Your word on PR #43 ("merge PR 43" in the session). It now holds AR3 steps 2 to 6.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (the commands are in the
  history file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Questions that block nothing: the shape and size of real orientation gradients; your
  local-optimization settings on real data; how many voxels a grain spans in your maps.
- A choice: which seed step goes behind an opt-in key. Default: run `S5`'s, after step 7.
- A note to read: C++ parallel mode keeps the higher confidence for a voxel that two clients
  fit; that rule did harm in the Python test of step 5. A reading of the source, not a test.
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 and N2 (all
  on the review page). The gradient levels of step 7's sample are my choice. Step 5's result
  files add about 50 MB to git; say so if you want them pruned.
- Not done unless you ask: a repo-wide rename of "finisher" to "local optimization".
- Ilium's plan approvals, when posted, are yours.

## Bookkeeping

- Next decision ID: D-20261111-N3.
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`.
