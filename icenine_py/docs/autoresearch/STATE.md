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
| open PR | #43 (AR3 steps 2 to 8), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below. Do every step that is ready, one after another, and independent
   jobs side by side (owner, 2026-11-10). After each step: rewrite this file, commit, push, one
   shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action

Step 9 of AR3 (a start test for the neighbour fits, the rejected fits that end right, a
third sample; task `ar3i`; plan in `PLAN.md` under "Step 9: ...", approved by D-20261111-N6).
Handed back (branch `feature/eval-reduction-ar3i`, head 86e3ffa, worktree
`agent-a9c88b610fb10e1b6`; check 0 FAIL). The review runs; its report goes to the main
`.claude/reports/20261111-0810-review-ar3i.md`.

1. Wait for the review. Then runbook, "A step, from plan to merge", points 7 to 9: the
   decision on the fixes (D-20261111-N7 or the day's next ID), a fixes round (Opus; free the
   branch first with `git -C <worktree> checkout --detach`), `finish_task.sh --yes --push`, close.
   Copy the reports to the main `.claude/reports/` first: a merged worktree can vanish.
2. At the close: add the runs on the third sample (`mid`) and run `S9`, if it ran, to
   `docs/autoresearch/pareto/pareto_rows.json` (a new panel for `mid`); say on the review page
   whether run `S8` is confirmed (prediction 3 of step 9) and check the other predictions.
3. Then: opt-in keys for the seed step and the boundary pass (library change, tests first,
   settings off by default; the owner chooses what a key recommends). Then AR5, then AR4.

## Running jobs

- Step 9 review (`code-reviewer`, Opus), started 2026-11-11. The implementer is done.
- The review page is at version 35 (run `S8` is on its plots). No page job is open.

## Queue

Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0, AR1, AR2, AR3 steps 1 to 8) are in `STATE_HISTORY.md` under "Queue".

| step | what | state |
|---|---|---|
| AR3 step 9 | start test, rejected fits that end right, third sample (task `ar3i`) | review running |
| AR3 later | opt-in keys for the seed step and the boundary pass | with settings; the owner picks |
| AR4 | the learned proposer | after AR5 |
| AR5 | the stop rule inside the grid search (the grid fallbacks) | after step 9 |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B | needs a quiet machine |
| C++ time | measured C++ time on the pilot region | not started |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py` misses numbers with thousands separators | small, not started |

## Last finished step

AR3 step 8 (where the seed step's evaluations go, task `ar3h`): merged into
`feature/eval-reduction` at ce9b106 on 2026-11-11, pushed; suite 1517 passed, 47 skipped, 1
deselected. Page: `docs/autoresearch/findings_ar3_reseed.md`. Decisions D-20261111-N4 (plan), N5
(review: 10 fixes, no rerun). Refinements are 99.9% of a seed step on `big`. The plan's early
stop (hit ratio 0.80) failed on stored cases; a variant (0.95), chosen after results, gives run
`S8`. `S8` plus pass: 8.34 million evaluations and 24 wrong on `big` (`D4` 11.94 and 21, `S5`
13.00 and 30); 18.00 million and 63 on the first gradient sample (`S5` 20.63 and 62). Prediction
3 is met by the variant only: not confirmed until step 9's third sample. 1,069 of 1,471 voxels
in reseed clusters had a fit that ended right and was rejected (class A); the cause is not
tested. More in `STATE_HISTORY.md`, "AR3 step 8".

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-a635b2f904279bdc9` (step 8, detached at 2aaa72a): `ar3/reseed/cache/` (the logged
  rerun and the `S8` maps), and links to the image stacks of both gradient samples.
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

- Your word on PR #43 ("merge PR 43" in the session). It now holds AR3 steps 2 to 8.
- A choice: which seed step the opt-in key recommends. Three candidates: the old one (`D4`),
  step 6's (`S5`), and `S5` with the early stop (`S8`). `S8` was cheapest on both samples but
  its stop test was chosen after results; I recommend nothing before step 9 has tested it.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (the commands are in the
  history file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Questions that block nothing: the shape and size of real orientation gradients; your
  local-optimization settings on real data; how many voxels a grain spans in your maps.
- A note to read: C++ parallel mode keeps the higher confidence for a voxel that two clients
  fit; that rule did harm in the Python test of step 5. A reading of the source, not a test.
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 to N6 (all
  on the review page). The gradient levels of step 7's sample are my choice. Step 5's result
  files add about 50 MB to git; say so if you want them pruned.
- Not done unless you ask: a repo-wide rename of "finisher" to "local optimization".
- Ilium's plan approvals, when posted, are yours.

## Bookkeeping

- Next decision ID: D-20261111-N7.
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`.
