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
| open PR | #43 (AR3 steps 2 to 10, AR5 steps 1 and 2, AR4 step 1), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below. Do every step that is ready, one after another, and independent
   jobs side by side (owner, 2026-11-10). After each step: rewrite this file, commit, push, one
   shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action

Two agents run. At each hand-back: runbook, "A step, from plan to merge", points 5
to 9. Merge the parent into a task branch before its fixes round (shared files conflict).

1. AR5 step 3 (task `ar5c`, N14): handed back 2026-11-11 15:00, in review (Opus; report
   `.claude/reports/20261111-1503-review-ar5c.md`). Branch `feature/eval-reduction-ar5c` in
   worktree `agent-a5a3601b746680640` (page placed by the main session, 6174af7). Unreviewed:
   pick `X200`; `S11` with the pass 51.9 / 31.6 / 39.2% fewer evaluations than `S8`, wrong 18
   / 62 / 36; passes on all three (`mid` at the limit). Then decision N16, fixes, merge, close.
2. AR3 step 11 (task `ar3k`, plan approved by N15; Sonnet; at most 2 heavy processes): the
   sample `held` (250 grains) and its `S8` baseline. Report `20261111-1351-ar3k.md`. No other
   run may read `held` before the plan of the step that uses it is committed.
3. After `S11` and `held`: plan the test on `held` (`S8` against `S11`, and `S11` plus the
   verifier `N1` as an added seed-accept rule), predictions first; a key only after it.

## Running jobs

- AR3 step 11 implementer (started 2026-11-11 about 13:55) and the AR5 step 3 reviewer
  (started 15:05); each hands back by message.
- No other job. The review page is at artifact version 41 (its small print says 40; the next
  build is 42). The text of PR #43 holds AR5 step 2 and AR4 step 1.

## Queue

Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0 to AR2, AR3 1 to 10, AR5 1 and 2, AR4 1) are in `STATE_HISTORY.md`, "Queue".

| step | what | state |
|---|---|---|
| AR5 step 3 | mixed seed refinement (CMA, MC under 0.95), run `S11` (task `ar5c`) | in review |
| AR3 step 11 | sample `held` (250 grains) with its `S8` baseline (task `ar3k`) | running |
| test on `held` | `S8`, `S11`, `S11` + verifier `N1` on the held-out sample | after both |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B | needs a quiet machine |
| C++ time | measured C++ time on the pilot region | not started |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py`: thousands separators; false misses in wide table rows | small |

## Last finished step

AR4 step 1 (a learned verifier against the hit ratio, offline; task `ar4a`; scripts only; no
model passed its rule): merged into `feature/eval-reduction` at 1ffee63 on 2026-11-11, pushed;
suite 1635 passed, 47 skipped, 1 deselected. Page: `docs/autoresearch/findings_ar4_verifier.md`.
Decisions D-20261111-N13 (plan), N15 (review). Before it: AR5 step 2 (c688191). More:
`STATE_HISTORY.md`.

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-acc78ac2a483ed86f` (AR5 step 2), `agent-ab349f0b6193c2ca1` (AR4 step 1): caches, links.
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

- Your word on PR #43 ("merge PR 43" in the session). It holds AR3 2 to 10, AR5 1 and 2, AR4 1.
- A choice: which seed step the documentation of the keys recommends. Step 10 added the
  keys and recommends none. My recommendation: the settings of `S8` (`SeedIndexCandidates 20
  5`, `SeedIndexStopHitRatio 0.95`): cheapest on all three samples, wrong counts not worse.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (the commands are in the
  history file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Three questions that block nothing (gradient shapes, your settings, voxels per grain).
- A note: C++ parallel mode keeps the higher confidence for a voxel that two clients fit;
  that rule did harm in the Python test of step 5. A reading of the source, not a test.
- A choice: try test `P2 m=2` as an opt-in key (no gain lost; cuts 17% to 39% of the pass).
- A choice: an opt-in key for "no MC fallback" alone (`S10`; N14 says not now).
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 to N15 (all
  on the review page). N13 makes AR4 a learned verifier; N15 adds the sample `held` (new seed).
  The gradient levels of step 7's sample are my choice. Step 5's result files add about 50
  MB to git; say so if you want them pruned.
- Not done unless you ask: a repo-wide rename of "finisher". Ilium's plan approvals are yours.

## Bookkeeping

- Next decision ID: D-20261111-N16.
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`.
