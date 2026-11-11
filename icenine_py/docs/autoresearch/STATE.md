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
| open PR | #43 (AR3 steps 2 to 9), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below. Do every step that is ready, one after another, and independent
   jobs side by side (owner, 2026-11-10). After each step: rewrite this file, commit, push, one
   shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action

Step 10 of AR3 (opt-in keys in the library: the indexed seed step, its early stop, the
CMA-first neighbour fit, the boundary pass; all off by default; task `ar3j`; plan in `PLAN.md`
under "Step 10: ...", approved by D-20261111-N8).
The implementer runs (Sonnet, own worktree, branch `feature/eval-reduction-ar3j` off 4ff3ce6;
report `.claude/reports/20261111-0834-ar3j.md`).

1. Wait for the implementer's hand-back. Then runbook, "A step, from plan to merge", points 5
   to 9: `verify_handback.sh`, review (`code-reviewer`, Opus), the decision on the fixes
   (D-20261111-N9 or the day's next ID), a fixes round (Opus; it may change the library files
   the step touched, D-20261111-N8), `finish_task.sh --yes --push`, close. Copy every report
   to the main `.claude/reports/` before a merge. In a brief, name the report fields as the
   hand-back check does (branch / commits; files; tests; criteria; suspicious / unexpected).
2. At the close: `D4L`, `S5L`, `S8L` must equal the stored runs, so they add no point to the
   plots; say so on the review page. Then AR5 (the stop rule inside the grid search), then AR4.

## Running jobs

- Step 10 implementer (`implementer`, Sonnet, own worktree), started 2026-11-11.
- The review page is at version 36 (`mid` and run `S9` are on its plots). No page job.

## Queue

Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0, AR1, AR2, AR3 steps 1 to 9) are in `STATE_HISTORY.md` under "Queue".

| step | what | state |
|---|---|---|
| AR3 step 10 | opt-in keys: seed step, early stop, CMA-first fit, pass (task `ar3j`) | implementer running |
| AR4 | the learned proposer | after AR5 |
| AR5 | the stop rule inside the grid search (the grid fallbacks) | after step 10 |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B | needs a quiet machine |
| C++ time | measured C++ time on the pilot region | not started |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py` misses numbers with thousands separators | small, not started |

## Last finished step

AR3 step 9 (a start test, the rejected fits that end right, a third sample; task `ar3i`): merged
into `feature/eval-reduction` at 0b081fb on 2026-11-11, pushed; suite 1517 passed, 47 skipped, 1
deselected. Page: `docs/autoresearch/findings_ar3_start.md`. Decisions D-20261111-N6 (plan), N7
(review: 13 fixes, no rerun). Third sample `mid`, with the pass: `D4` 40.56 million evaluations and
48 wrong, `S5` 39.09 and 34, `S8` 26.44 and 31; `S8` counts as confirmed on one new sample (N7).
Run `S9` (`S8` plus a start test at hit ratio 0.90) saves 14% to 32% but breaks more voxels than it
fixes on all three samples: not a candidate for a key. One more fit from a rejected end passes the
relative test in 807 of 1,069 and 451 of 589 voxels; the remedy built on it (`F2`) failed its test
(4.8% pooled, a third asked). More in `STATE_HISTORY.md`, "AR3 step 9".

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-a9c88b610fb10e1b6` (step 9, detached at 86e3ffa): `ar3/start/cache/`, the image
  stack and the samples of `mid`, and links.
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

- Your word on PR #43 ("merge PR 43" in the session). It now holds AR3 steps 2 to 9.
- A choice: which seed step the documentation of the keys recommends. Step 10 adds the keys
  and recommends none. My recommendation: the settings of `S8` (20 candidates costed, 5
  refined, stop at hit ratio 0.95): cheapest on all three samples, wrong counts not worse.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (the commands are in the
  history file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Questions that block nothing: the shape and size of real orientation gradients; your
  local-optimization settings on real data; how many voxels a grain spans in your maps.
- A note to read: C++ parallel mode keeps the higher confidence for a voxel that two clients
  fit; that rule did harm in the Python test of step 5. A reading of the source, not a test.
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 to N8 (all
  on the review page). The gradient levels of step 7's sample are my choice. Step 5's result
  files add about 50 MB to git; say so if you want them pruned.
- Not done unless you ask: a repo-wide rename of "finisher" to "local optimization".
- Ilium's plan approvals, when posted, are yours.

## Bookkeeping

- Next decision ID: D-20261111-N9.
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`.
