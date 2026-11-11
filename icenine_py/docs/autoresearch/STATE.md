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
| open PR | #43 (AR3 steps 2 to 10, AR5 step 1), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below. Do every step that is ready, one after another, and independent
   jobs side by side (owner, 2026-11-10). After each step: rewrite this file, commit, push, one
   shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action

Two implementers run (Sonnet, own worktrees). Each hands back by message; then runbook, "A
step, from plan to merge", points 5 to 9 for each. Both add a section to
`MIGRATION_HISTORY.md`: merge the parent into the second task branch before its fixes round.

1. AR5 step 2 (seed refinement by CMA, no MC fallback after CMA, run `S10`; task `ar5b`;
   scripts only; plan D-20261111-N12). Report `.claude/reports/20261111-1134-ar5b.md` in its
   worktree. Finished `S10` runs add rows to the plots of cost against accuracy
   (`docs/autoresearch/pareto/pareto_rows.json`). Its branch was cut before AR5 step 1 merged.
2. AR4 step 1 (a learned verifier against the hit ratio, offline; task `ar4a`; plan
   D-20261111-N13). Report `.claude/reports/20261111-1208-ar4a.md` in its worktree.
3. After both: plan the next step from their results (a library key for `S10`'s changes if it
   passes; the verifier inside the BFS if a model passes).

## Running jobs

- AR5 step 2 implementer (started 2026-11-11 about 11:40), at most 3 heavy processes.
- AR4 step 1 implementer (started 2026-11-11 about 12:40), at most 2 heavy processes.
- The review page is at version 38 (AR5 step 1 is on it; no new plot point). No page job.

## Queue

Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0 to AR2, AR3 steps 1 to 10, AR5 step 1) are in `STATE_HISTORY.md`, "Queue".

| step | what | state |
|---|---|---|
| AR5 step 2 | CMA seed refinement, no MC fallback, run `S10` (task `ar5b`) | running |
| AR4 step 1 | learned verifier against the hit ratio, offline (task `ar4a`) | running |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B | needs a quiet machine |
| C++ time | measured C++ time on the pilot region | not started |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py`: thousands separators; false misses in wide table rows | small |

## Last finished step

AR5 step 1 (stop rules for the boundary pass, the grid fallbacks and the finished map; task
`ar5a`; scripts only): merged into `feature/eval-reduction` at e51e883 on 2026-11-11, pushed;
suite 1596 passed, 47 skipped, 1 deselected. Page: `docs/autoresearch/findings_ar5_stop.md`.
Decisions D-20261111-N10 (plan), N11 (review: 12 fixes to the documents). A negative result: no
rule passed, run `T1` was not run. More in `STATE_HISTORY.md`, "AR5 step 1".

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-a68538c58665829a1` (AR5 step 1, detached at 29e9b56): `ar5/cache/` and links.
- `agent-a6c7c67b4a03d1843` (step 10, detached at 84ce8b5): `ar3/keys/cache/` (28 MB) and
  links; `agent-a3145f902bf0a61fb` (step 10 fixes): a copy of it plus the two reruns.
- `agent-a9c88b610fb10e1b6` (step 9, detached at 86e3ffa): `ar3/start/cache/`, the image
  stack and the samples of `mid`, and links.
- `agent-a635b2f904279bdc9` (step 8, detached at 2aaa72a): `ar3/reseed/cache/` (the logged
  rerun and the `S8` maps), and links to the image stacks of both gradient samples.
- `agent-a38ba6f5a25aa033e` (step 7, detached at b762441): its caches and runs under
  `ar3/big/`, and the report `.claude/reports/20261111-0151-ar3g.md`.
- `agent-a6db3cabf8ab230e9` (step 6, detached at 6c4665c): `seed/cache/`, `seed/runs/`, a
  `boundary/cache/`, the image-stack link, and its report; `agent-a27b4a932b0dd99b6` holds
  step 6's fixes report.
- `agent-a4fd6282107a76186` (detached at 3813d69): step 5's `boundary/cache/` (the maps after
  the pass, `after_<map>.npy`, and the run environments), the image-stack link.
- `agent-ad88d4ade1ca5a039` (step 5 fixes): a second `boundary/cache/` and the fixes report.
- Steps 3 and 4: `agent-ae68d32fb972a8ad4` and `agent-ac609f1c28c619ed0` (see
  `grep -n "worktrees/agent-" docs/autoresearch/STATE_HISTORY.md`).

## For the owner (open items only; details by `grep` in `STATE_HISTORY.md`, "For the owner")

- Your word on PR #43 ("merge PR 43" in the session). It now holds AR3 steps 2 to 10 and AR5 step 1.
- A choice: which seed step the documentation of the keys recommends. Step 10 added the
  keys and recommends none. My recommendation: the settings of `S8` (`SeedIndexCandidates 20
  5`, `SeedIndexStopHitRatio 0.95`): cheapest on all three samples, wrong counts not worse.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (the commands are in the
  history file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Three questions that block nothing (gradient shapes, your settings, voxels per grain).
- A note: C++ parallel mode keeps the higher confidence for a voxel that two clients fit;
  that rule did harm in the Python test of step 5. A reading of the source, not a test.
- A choice: try test `P2 m=2` as an opt-in key (no gain lost; cuts 17% to 39% of the pass).
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 to N13 (all
  on the review page). N13 recasts AR4 from a learned proposer to a learned verifier.
  The gradient levels of step 7's sample are my choice. Step 5's result files add about 50
  MB to git; say so if you want them pruned.
- Not done unless you ask: a repo-wide rename of "finisher". Ilium's plan approvals are yours.

## Bookkeeping

- Next decision ID: D-20261111-N14.
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`.
