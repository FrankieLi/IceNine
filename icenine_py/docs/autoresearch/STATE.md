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

One step is open. Runbook: "A step, from plan to merge", points 5 to 9. Merge the parent into
a task branch before its fixes round (shared files conflict).

1. AR5 step 4 (task `ar5d`; N16, N18; plan in `PLAN.md`, predictions at 17a27c1): `S10`,
   `S11`, `S11v` (= `S11` plus verifier `N1`) on `held`. The runs and the statistics are done
   (branch tip 72be082; the page was placed by the main session; report
   `20261111-1640-ar5d.md`). Reported, not yet reviewed: with the pass `S11` needs 23.3% fewer
   evaluations than `S8` and passes; `S10` 3.6% (fails); `S11v` 36.3% but 82 wrong against 70
   (fails). Predictions 3, 4 met; 1, 2, 5, 6 not. Now: read the review, decision N20, Opus
   fixes round, merge, close (page with the `held` points of the three runs, PR #43).
2. Then, if `S11` still passes after the review: plan the library step (one opt-in key pair:
   mixed seed refinement and "no MC fallback"), with predictions; decision; Sonnet implementer.

## Running jobs

- Opus review of AR5 step 4, started 2026-11-11 about 19:35; report
  `.claude/reports/20261111-2000-review-ar5d.md`. No machine-time run is going.
- The review page is at version 42; version 43 (AR3 step 11) is being built
  (`build_v43.py`, `run43.sh` in the scratchpad). PR #43's text: see its last section.

## Queue

Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0 to AR2, AR3 1 to 11, AR5 1 to 3, AR4 1) are in `STATE_HISTORY.md`, "Queue".

| step | what | state |
|---|---|---|
| AR5 step 4 | `S10`, `S11`, `S11v` on the held-out sample (task `ar5d`) | in review |
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
  cache (`S8_held/`, `S8_held_cap120/`). `agent-a9fcb8cf18fe801c4`: links only.
- `agent-a814d0121e137b01c` (AR5 step 4, branch `ar5d`): the three run caches on `held`.
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
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 to N19 (all
  on the review page). N15 adds `held` (new seed); N17, N18: run caps of 480 minutes.
  The gradient levels of step 7's sample are my choice. Step 5's result files add about 50
  MB to git; say so if you want them pruned.

## Bookkeeping

- Next decision ID: D-20261111-N20.
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`. Its plan approvals are the owner's; so is a "finisher" rename.
