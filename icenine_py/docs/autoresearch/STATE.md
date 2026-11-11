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
| open PR | #43 (AR3 steps 2 to 5), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below. Do every step that is ready, one after another, and independent
   jobs side by side (owner, 2026-11-10). After each step: rewrite this file, commit, push, one
   shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action

Step 6 of AR3 (the seed cost, task `ar3f`; plan in `PLAN.md` under "Step 6: the seed cost",
approved by D-20261110-N8).

1. Wait for the fixes implementer's notification (see "Running jobs"). Do not poll its output
   file. Then `verify_handback.sh feature/eval-reduction <its report> --branch
   feature/eval-reduction-ar3f --fixes .claude/reports/20261111-0330-review-ar3f.md --no-tests`
   and the key-by-key check of the summary files against 6c4665c.
2. Merge and close: runbook section "A step, from plan to merge", points 8 and 9. Add step 6's
   runs (S, S5, K2, each with the boundary pass) to the cost-against-accuracy panels.
3. Then plan the next step from step 6's result. Candidates, in the plan's order: the same
   check on a sample with larger grains, then opt-in keys for the remedy and for the boundary
   pass; AR4 (learned proposer) for the voxels every remedy leaves to the grid search; AR5
   (stop rule inside the grid search).

If an agent died with the session (a task branch exists but no report): look at
`git log feature/eval-reduction..<task branch>`, then start a fresh implementer with the same
brief and the note "continue from the commits on the branch".

## Running jobs

- Step 6 is handed back (2026-11-11): task branch `feature/eval-reduction-ar3f` at 6c4665c, in
  the worktree `.claude/worktrees/agent-a6db3cabf8ab230e9` (report in its
  `.claude/reports/20261110-2308-ar3f.md`). The findings page `findings_ar3_seed.md` was placed
  by the main session from the implementer's text (its Write was refused). `verify_handback.sh`:
  0 FAIL, 2 expected flags.
- Step 6 review (2026-11-11): request changes, 11 fixes, text, tables and plot label only, no
  rerun. Report `.claude/reports/20261111-0330-review-ar3f.md`. Decision D-20261111-N1: all 11
  fixes; additions 7a, 7b, S5 kept as marked; default for the next work is the seed step of
  run S5 (best 5 of 20, no neighbour tries), after a larger-grain check; the owner can overrule.
- Step 6 fixes implementer (Opus, own worktree), started 2026-11-11. Its report:
  `.claude/reports/20261111-0128-ar3f-fixes.md` (maybe in its own worktree). The old worktree
  is detached at 6c4665c, so the task branch is free for it.
- The review page is at version 30 (step 5 done, step 6 running). No page job is open. Version
  30 added the section "Cost against accuracy, all methods" (owner's request, 2026-11-11): one
  plot per sample with the frontier marked. When a step ends, add its rows to
  `docs/autoresearch/pareto/pareto_rows.json` and rebuild the section (runbook, "The review
  page").

## Queue

Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0, AR1, AR2, AR3 steps 1 to 5) are in `STATE_HISTORY.md` under "Queue".

| step | what | state |
|---|---|---|
| AR3 step 6 | the seed cost (task `ar3f`) | reviewed; fixes round running (D-20261111-N1) |
| AR3 later | larger-grain check; opt-in keys for the boundary pass and the remedy | after step 6 |
| AR4 | the learned proposer | after step 6 |
| AR5 | the stop rule inside the grid search | after step 6 |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B | needs a quiet machine |
| C++ time | measured C++ time on the pilot region | not started |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py` misses numbers with thousands separators | small, not started |

## Last finished step

AR3 step 5 (boundary pass, task `ar3e`): merged into `feature/eval-reduction` at b69a66f on
2026-11-10, pushed; suite 1520 passed, 44 skipped, 1 deselected. Page:
`docs/autoresearch/findings_ar3_boundary.md`. Decisions D-20261110-N6 (plan), N7 (review).
With the pass on both sides: 74 against 88 wrong at 38.8 against 122.1 million evaluations
(BFS, gradient sample); 82 against 78 at 54.6 against 94.7 million (paint, gradient sample);
0 against 14 at 16.6 against 36.1 million (step 2's sample). Run `K` did not finish in its cap;
step 6 repeats it as `K2`. Reports: `.claude/reports/20261110-1836-ar3e.md`,
`20261110-2220-review-ar3e.md`; the fixes report is in the fixes worktree (below).

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-a4fd6282107a76186` (detached at 3813d69): step 5's `boundary/cache/` (the maps after
  the pass, `after_<map>.npy`, and the run environments `env_grad`, `env_layer`), the link to
  the image stacks `Examples/Example2.ManyGrains/ScatteringData_PhaseD`.
- `agent-ad88d4ade1ca5a039` (step 5 fixes): a second copy of `boundary/cache/`, and the fixes
  report `.claude/reports/20261110-2233-ar3e-fixes.md`.
- `agent-a6db3cabf8ab230e9` (step 6, detached at 6c4665c): `seed/cache/`, `seed/runs/`, a
  `boundary/cache/`, the image-stack link, and the report `.claude/reports/20261110-2308-ar3f.md`.
- Steps 3 and 4: `agent-ae68d32fb972a8ad4` and `agent-ac609f1c28c619ed0` (see
  `grep -n "worktrees/agent-" docs/autoresearch/STATE_HISTORY.md`).

## For the owner (open items only; details by `grep` in `STATE_HISTORY.md`, "For the owner")

- Your word on PR #43 ("merge PR 43" in the session). It now holds AR3 steps 2 to 5.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (the commands are in the
  history file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Questions that block nothing: the shape and size of real orientation gradients; your
  local-optimization settings on real data; how many voxels a grain spans in your maps.
- A note to read: C++ parallel mode keeps the higher confidence for a voxel that two clients
  fit; that rule did harm in the Python test of step 5. A reading of the source, not a test.
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8 and D-20261111-N1 (the
  review page lists them up to N8). Step 5's result files add about 50 MB to git; say so if you want them pruned.
- Not done unless you ask: a repo-wide rename of "finisher" to "local optimization".
- Ilium's plan approvals, when posted, are yours.

## Bookkeeping

- Next decision ID: D-20261111-N2.
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`); mailbox in
  `docs/ROADMAP.md`.
