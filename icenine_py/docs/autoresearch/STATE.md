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
| open PR | #43 (AR3 steps 2 to 11, AR5 steps 1 to 3 and 5, AR4 step 1), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below. Do every step that is ready, one after another, and independent
   jobs side by side (owner, 2026-11-10). After each step: rewrite this file, commit, push, one
   shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action

One step is open. Runbook: "A step, from plan to merge", points 5 to 9. Merge the parent
into a task branch before its fixes round (shared files conflict).

1. AR5 step 4 (task `ar5d`; N16, N18, N20, N23): on `held` with the pass `S11` needs 23.3%
   fewer evaluations than `S8`; the rerun `S11v2` (valid arm) 37.4% fewer, 18.4% fewer than
   `S11`, 69 wrong, broken minus fixed 4 against `S11` (limit 5). Phase B re-review: approve,
   three fixes (N23). Fix 3 is done (fac79f9 on `ar5d`).
   BLOCKED on fix 1: the page text is in `.claude/reports/20261111-2330-ar5d-phaseB-page.md`;
   the implementer's Write and the main session's copy were refused (auto mode). No way round:
   the owner places it or allows it. Then: sync check, commit, fix 2 (MIGRATION_HISTORY names
   that commit), merge, close (page version 45, PR #43).
2. AR5 step 6 (task `ar5f`; plan in `PLAN.md`, 2ed6618; N24): the frozen `N1` rule (`S11n`)
   and a plain hit-ratio rule (`S11h`) on `mid` and `grad`. Its implementer (Sonnet,
   worktree) starts only after step 4 is merged.
3. Items 14 + 3 (task `item14`; plan in `PLAN.md`; N25): phase A is with the implementer.
   On "READY TO RUN": read `run_item14.sh`, start it as a background command, wait for its
   exit; then phase B by `SendMessage`. While the timing arms run (about 2.6 hours): no
   merge, no test suite, no other run. The C++ time is a separate step, not planned yet.

## Running jobs

- Items 14 + 3 phase A: `implementer` agent `a012b8f461d1630a1` (Sonnet, own worktree,
  started 22:35 on 2026-11-11); reports by hand-back, report `20261111-2300-item14.md`.
  Branch `ar5d` is held by worktree `agent-a8a1ac8e80f8913cf` (tip fac79f9).
- Review page: version 44 (`run43.sh`, then `patch_v44.py`, scratchpad).

## Queue

Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0 to AR2, AR3 1 to 11, AR5 1 to 3 and 5, AR4 1) are in `STATE_HISTORY.md`, "Queue".

| step | what | state |
|---|---|---|
| AR5 step 4 | `S10`, `S11`, `S11v2` on the held-out sample (task `ar5d`) | waits for the owner (page) |
| AR5 step 6 | `S11n`, `S11h` on `mid` and `grad` (task `ar5f`) | planned; after step 4 |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B (task `item14`) | phase A |
| C++ time | measured C++ time on the pilot region (MPI mode, own config) | not planned |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py`: thousands separators; false misses in wide table rows | small |

## Last finished step

AR5 step 5 (two opt-in library keys for run `S11`; task `ar5e`): merged at a494bcd on
2026-11-11; suite 1724 passed, 47 skipped, 1 deselected. Page:
`docs/autoresearch/findings_ar5_keys.md`. Decisions N21, N22. More: `STATE_HISTORY.md`.

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
- Step 6: `agent-a6db3cabf8ab230e9` (`seed/`, a `boundary/cache/`), `agent-a27b4a932b0dd99b6`.
- Step 5: `agent-a4fd6282107a76186` (`boundary/cache/`, the maps after the pass) and
  `agent-ad88d4ade1ca5a039` (fixes). Steps 3 and 4: `agent-ae68d32fb972a8ad4`,
  `agent-ac609f1c28c619ed0`. More: `grep -n "worktrees/agent-" docs/autoresearch/STATE_HISTORY.md`.

## For the owner (open items only; details by `grep` in `STATE_HISTORY.md`, "For the owner")

- Place the page of AR5 step 4 (or allow me to): `.claude/reports/20261111-2330-ar5d-phaseB-page.md`
  goes to `icenine_py/docs/autoresearch/findings_ar5_held.md` in worktree `agent-a8a1ac8e80f8913cf`.
- Your word on PR #43 ("merge PR 43" in the session): AR3 2 to 11, AR5 1 to 3 and 5, AR4 1.
- A choice: which seed step the documentation of the keys recommends. Step 10 added the
  keys and recommends none. My recommendation: the settings of `S8` (`SeedIndexCandidates 20
  5`, `SeedIndexStopHitRatio 0.95`): cheapest on all three samples, wrong counts not worse.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (the commands are in the
  history file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Three questions that block nothing (gradient shapes, your settings, voxels per grain).
- A note: C++ parallel mode keeps the higher confidence for a voxel that two clients fit;
  that rule did harm in the Python test of step 5. A reading of the source, not a test.
- A choice: try test `P2 m=2` as an opt-in key (no gain lost; cuts 17% to 39% of the pass).
- A choice: whether the documentation recommends the two `S11` keys. They were added off by
  default under N21, after `S11` passed on `held`; the documentation recommends none now.
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 to N25 (all
  on the review page). N15 adds `held` (new seed); N17, N18: run caps of 480 minutes.
  The gradient levels of step 7's sample are my choice. Step 5's result files add about 50
  MB to git; say so if you want them pruned.

## Bookkeeping

- Next decision ID: D-20261111-N26.
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`. Its plan approvals are the owner's; so is a "finisher" rename.
