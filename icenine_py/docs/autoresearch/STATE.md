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
| open PR | #43 (AR3 steps 2 to 11, AR5 steps 1 to 3, AR4 step 1), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below. Do every step that is ready, one after another, and independent
   jobs side by side (owner, 2026-11-10). After each step: rewrite this file, commit, push, one
   shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action

Two things are open. Runbook: "A step, from plan to merge", points 5 to 9. Merge the parent
into a task branch before its fixes round (shared files conflict).

1. AR5 step 4 (task `ar5d`; N16, N18, N20): on `held` with the pass `S11` needs 23.3% fewer
   evaluations than `S8` and passes; `S10` 3.6% (fails). Both reviewed. `S11v` was a defective
   arm. Fixes round (Opus), phase A done: hook repaired in scripts, checked on `mid`. The
   rerun `S11v2` goes now. When it ends: phase B (`SendMessage` to the fixes implementer, or
   a new Opus one): statistics, page, predictions 5 and 6. Place the page if its Write is
   refused. Then hand-back check, merge (`finish_task.sh`), close (page with the `held`
   points, PR #43).
2. AR5 step 5 (task `ar5e`; N21): keys `SeedIndexCMAFirst`, `BFSNeighborCMANoFallback`.
   Phase A done (tip fb5ecc3): 62 tests, suite 1718 passed, checks 1 to 3 on `big` exact.
   Check 4 (first gradient sample) runs now. When it ends: phase B (`SendMessage` to the
   implementer): results, page `findings_ar5_keys.md`. Then hand-back check, review, N22.

## Running jobs

- `S11v2` on `held` (cap 480 minutes): background command of the main session, started
  2026-11-11 19:54 in worktree `agent-a8a1ac8e80f8913cf` (branch `ar5d`, tip 1983a40; report
  `20261111-2130-ar5d-fixes.md`). Marker there: `ar5/held/cache/logs/S11v2_held.done`.
- AR5 step 5 chain (check 4 and the pass checks, cap 90 minutes), started 20:13 in worktree
  `agent-afc3f195cb4fffb30` (report `20261111-2030-ar5e.md`); marker `grad_chain.done` in its logs.
- Review page: version 43 (builder `run43.sh`, scratchpad).

## Queue

Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0 to AR2, AR3 1 to 11, AR5 1 to 3, AR4 1) are in `STATE_HISTORY.md`, "Queue".

| step | what | state |
|---|---|---|
| AR5 step 4 | `S10`, `S11`, `S11v` on the held-out sample (task `ar5d`) | fixes round, one rerun |
| AR5 step 5 | library keys for `S11` (task `ar5e`) | check 4 runs |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B | needs a quiet machine |
| C++ time | measured C++ time on the pilot region | not started |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py`: thousands separators; false misses in wide table rows | small |

## Last finished step

AR3 step 11 (sample `held`, 250 grains, with its `S8` baseline; task `ar3k`; scripts only):
merged into `feature/eval-reduction` at 504fadb on 2026-11-11, pushed; suite 1656 passed, 47
skipped, 1 deselected. `S8` on `held`: 59,923,523 evaluations with the pass, 75 grid
fallbacks, 70 wrong after the pass. Page: `docs/autoresearch/findings_ar3_held.md`.
Decisions D-20261111-N15 (plan), N17 (cap 480 minutes), N19 (review). Before it: AR5 step 3
(ab7e9fc). More: `STATE_HISTORY.md`.

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-acc78ac2a483ed86f` (AR5 step 2), `agent-ab349f0b6193c2ca1` (AR4 step 1): caches, links.
- `agent-a5a3601b746680640` (AR5 step 3): the `S11` caches and maps under `ar5/mix/`, links.
- `agent-af8961874e131bda2` (AR3 step 11): the sample `held`, its image stack, the baseline
  cache (`S8_held/`, `S8_held_cap120/`).
- AR5 step 4: `agent-a814d0121e137b01c` (three run caches on `held`), `agent-a8a1ac8e80f8913cf`.
- AR5 step 5: `agent-afc3f195cb4fffb30` (the library runs under `ar5/keys/cache/`).
- `agent-a68538c58665829a1` (AR5 step 1, detached at 29e9b56): `ar5/cache/` and links.
- Step 10: `agent-a6c7c67b4a03d1843` (`ar3/keys/cache/`), `agent-a3145f902bf0a61fb` (two reruns).
- `agent-a9c88b610fb10e1b6` (step 9): `ar3/start/cache/`, the image stack and samples of `mid`.
- `agent-a635b2f904279bdc9` (step 8): `ar3/reseed/cache/` (logged rerun, `S8` maps), stack links.
- `agent-a38ba6f5a25aa033e` (step 7): its caches and runs under `ar3/big/`, and its report.
- `agent-a6db3cabf8ab230e9` (step 6): `seed/cache/`, `seed/runs/`, a `boundary/cache/`, its
  report; `agent-a27b4a932b0dd99b6`: step 6's fixes report.
- Step 5: `agent-a4fd6282107a76186` (`boundary/cache/`, the maps after the pass) and
  `agent-ad88d4ade1ca5a039` (fixes). Steps 3 and 4: `agent-ae68d32fb972a8ad4`,
  `agent-ac609f1c28c619ed0`. More: `grep -n "worktrees/agent-" docs/autoresearch/STATE_HISTORY.md`.

## For the owner (open items only; details by `grep` in `STATE_HISTORY.md`, "For the owner")

- Your word on PR #43 ("merge PR 43" in the session). It holds AR3 2 to 11, AR5 1 to 3, AR4 1.
- A choice: which seed step the documentation of the keys recommends. Step 10 added the
  keys and recommends none. My recommendation: the settings of `S8` (`SeedIndexCandidates 20
  5`, `SeedIndexStopHitRatio 0.95`): cheapest on all three samples, wrong counts not worse.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (the commands are in the
  history file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Three questions that block nothing (gradient shapes, your settings, voxels per grain).
- A note: C++ parallel mode keeps the higher confidence for a voxel that two clients fit;
  that rule did harm in the Python test of step 5. A reading of the source, not a test.
- A choice: try test `P2 m=2` as an opt-in key (no gain lost; cuts 17% to 39% of the pass).
- A choice: an opt-in key for "no MC fallback" (`S10`) or the mixed seed refinement (`S11`)
  now; N14 and N16 say not before the test on `held`.
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 to N21 (all
  on the review page). N15 adds `held` (new seed); N17, N18: run caps of 480 minutes.
  The gradient levels of step 7's sample are my choice. Step 5's result files add about 50
  MB to git; say so if you want them pruned.

## Bookkeeping

- Next decision ID: D-20261111-N22.
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`. Its plan approvals are the owner's; so is a "finisher" rename.
