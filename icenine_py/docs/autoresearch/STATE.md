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
| open PR | #43 (AR3 steps 2 to 11, AR5 steps 1 to 3, 5, 7 to 12, AR4 step 1), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below: every ready step in turn, independent jobs side by side (owner,
   2026-11-10). After each step: rewrite this file, commit, push, one shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action

Runbook: "A step, from plan to merge", points 5 to 9. Merge the parent into a task branch
before its fixes round (shared files conflict).

1. AR5 step 4 (task `ar5d`; N16, N18, N20, N23): `S11v2` on `held` 37.4% fewer evaluations
   than `S8`, 69 wrong; approved with three fixes (N23); fix 3 done (fac79f9). BLOCKED on fix 1
   (the page; agent writes refused): the owner places it or allows it. Then: sync check,
   commit, fix 2 (MIGRATION_HISTORY names that commit), merge, close (PR #43).
2. AR5 step 6 (`ar5f`; N3): fixed, waits. Steps 9 to 12: merged (10, 11 failed studies, no key).
   Step 12 (`ar5l`; key `SeedIndexMCStallHit`, off): merged at fe3d277, 13 Nov; page waits.
   Step 13 (`ar5m`): BLOCKED FOR THE OWNER (Queue). Nothing else is ready for an agent.
3. Items 14 + 3 (task `item14`; N25): runs 1 and 2 (13 Nov) stopped with exit 5: B1 not quiet
   (macOS indexing daemon). Rule amended (D-20261114-N1, R4b); implementer `a012b8f461d1630a1`
   changes the scripts (resumed 01:00, 14 Nov). On its report: check the commit, start run 3
   detached (command below; marker `item14.done`). After exit 0: phase B, then `cpptime`.

## Running jobs
- Items 14 + 3: `caffeinate -i bash scripts/autoresearch/item14/run_item14.sh` in `icenine_py/`
  of worktree `agent-a012b8f461d1630a1` (run 3 not started yet); marker `item14.done`.
  Branch `ar5d` is held by worktree `agent-a8a1ac8e80f8913cf` (tip fac79f9).
- Page v49 live (17:38, 13 Nov; v50 after the next finished step). Step 13: blocked (Queue).

## Queue
Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0 to AR2, AR3 1 to 11, AR5 1 to 3, 5, 7 to 12, AR4 1): `STATE_HISTORY.md`.

| step | what | state |
|---|---|---|
| AR5 step 4 | `S10`, `S11`, `S11v2` on the held-out sample (task `ar5d`) | waits for the owner (page): `.claude/reports/20261111-2330-ar5d-phaseB-page.md` (main checkout) goes to `icenine_py/docs/autoresearch/findings_ar5_held.md` in worktree `agent-a8a1ac8e80f8913cf` |
| AR5 step 6 | `S11n`, `S11h` on `mid` and `grad` (task `ar5f`) | reviewed and fixed; waits for the owner (page) and for step 4's merge. Branch tip d643eff, held by the fixes worktree `agent-a9420a98460078e81`; the first worktree `agent-a89458addd73332d9` (detached, 4c18740) holds the run cache `ar5/ver/cache/`: keep both. Opus review (`.claude/reports/20261112-0510-review-ar5f.md`): every headline number recomputed, none differs; 11 wording fixes, no rerun (N3). Fixes round done (7f272ce, 18260d5, d643eff; report `.claude/reports/20261112-0517-ar5f-fixes.md` in the fixes worktree; hand-back check: only the expected format flag; no number changed, `predictions[4].met` now false; 14 tests pass). Numbers: on `mid` `S11n` 14,105,137 evaluations (12.3% fewer than `S11`, 6 seeds kept, 0 wrong), `S11h` 15,680,025 (2.5%, 1 kept); both pass; "plain does as well" not met on `mid` (20.4% of the saving; 62.2% on `grad`, in-sample); six of seven predictions met. PAGE: not in the repo (Write refused). Corrected text: `.claude/reports/20261112-ar5f-page-v2.md` in the fixes worktree; the owner places it at `icenine_py/docs/autoresearch/findings_ar5_ver.md` there (or allows it). Then: sync check, commit, merge after step 4 (`finish_task.sh` in the fixes worktree), close (pareto points, page, PR #43); then a decision on an opt-in key for the `N1` rule |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B (task `item14`) | E0 done; timing started detached ~17:45, 13 Nov (see "Next action") |
| C++ time | measured C++ time on the pilot region (MPI mode, own config) | D-20261112-N1, task `cpptime`: phase A done (agent and worktree `ad44bb90c441eca35`, tip 5d6cec2, report in that worktree's `.claude/reports/20261112-0012-cpptime.md`). Run, after item 14's timing and on a quiet machine: `caffeinate -i bash scripts/autoresearch/cpptime/run_cpp_time.sh` in its `icenine_py/` (background command); marker `benchmarks/autoresearch/cpptime/cache/cpptime.done`; then phase B by `SendMessage`. R4 limit 3.5 (own load 2.0: the server rank spins). `lldb` cannot run (developer mode off): evaluations missing unless the owner enables it or accepts a print in `Src/` |
| AR5 step 7 | opt-in key `BoundaryPassSkipHitRatio` (task `ar5g`) | MERGED at 50c3b9c (N4, N5, N6), see "Last finished step". Open: only its page, for the owner to place: text `.claude/reports/20261112-ar5g-page-v2.md` (main checkout; the merge removed the fixes worktree, the text was recovered exactly from the recorded write and a table sync) goes to `icenine_py/docs/autoresearch/findings_ar5_gate.md` on `feature/eval-reduction`, then a sync check and a commit. Keep worktree `agent-a331e7882919e74aa` (stored passes, `ar5/gate/`). Step 6's page text and reports are now also copied to the main `.claude/reports/` |
| AR5 step 8 | sweep order of the seed step (task `ar5h`) | MERGED at 949dc61 as a study; failed rule 5 (N9), no key; see "Last finished step". Open: only its page for the owner: `.claude/reports/20261112-ar5h-page.md` (main checkout) goes to `icenine_py/docs/autoresearch/findings_ar5_order.md`. Keep worktree `agent-abdf64302f4022c0a` (detached): the step 8 cache `ar5/order/cache` that step 9 links, and the MC traces. Reports also in the main `.claude/reports/`: `20261112-1419-ar5h.md`, `-corrected.md`, `20261112-2215-ar5h-fixes.md` |
| AR5 step 9 | rule `G1`, key `SeedRescueCMA` (task `ar5i`) | MERGED at 28b75fa. Open: its page for the owner: `.claude/reports/20261113-ar5i-page.md` goes to `icenine_py/docs/autoresearch/findings_ar5_rescue.md`. Keep worktrees `agent-a8d20b3821621fd41`, `agent-ab6ee3e8ef1058af1` (`ar5/rescue/cache`) |
| AR5 step 12 | rule PH(3000), key `SeedIndexMCStallHit` (task `ar5l`) | MERGED at fe3d277 (13-N7, N8, N10), see "Last finished step". Open: only its page for the owner: `.claude/reports/20261113-ar5l-page-v2.md` (main checkout) goes to `icenine_py/docs/autoresearch/findings_ar5_stallhit.md`. Keep worktrees `agent-a0ffe667a1a6ee068` (held2 sample, runs; `passes_held2/D4_H/tries.npz` stays unopened until step 13's freeze) and `agent-ab6f021fc3d6bed02` (LibH, Hl runs on `mid`). Reports in the main `.claude/reports/`: `20261113-1100-ar5l.md`, `-1600-review-ar5l.md`, `-1630-ar5l-B.md`, `-1730-review-ar5l-B.md`, `-1800-ar5l-fixes2.md` |
| AR5 step 13 | arm HK = H + pass skip key 0.95, best current recipe; clean test on held2 (task `ar5m`) | D-20261113-N9 (6d13c12), PLAN "AR5 step 13". BLOCKED FOR THE OWNER since 17:05, 13 Nov (agent `afef88a15cf24bade`, worktree `.claude/worktrees/agent-afef88a15cf24bade`, tip be8d2fc, report `.claude/reports/20261113-1700-ar5m.md` there; no freeze commit; held2 tries unopened). (1) The agent's run of `recipe_freeze.py` was refused by the permission classifier; the main session does NOT run it for the agent; the owner runs it or allows it. (2) Gate g2 cannot hold as written: the library pass uses one generator, H's stored passes used aligned draws per try (mid: same skipped set 13,782, same 15,675 kept tries and 956,175 evaluations, only 105 tries identical, wrong 37 live vs 35 replay). Recommended amendment (needs a decision before any held2 read): g2 compares a script-layer aligned-draw pass with the skip against the replay; the library pass is compared on skipped set, counts and own hit ratios. Phase B's precondition (step 12 merged) is now met. |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py`: thousands separators; false misses in wide table rows | small |

## Last finished step

AR5 step 12 (`ar5l`): merged at fe3d277, 13 Nov; suite 1976 passed. Rule PH(3000) passes rules
1 to 5 on 5 samples; opt-in key `SeedIndexMCStallHit` (off); 3.8% to 8.1% fewer evaluations.

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-acc78ac2a483ed86f` (AR5 step 2), `agent-ab349f0b6193c2ca1` (AR4 step 1): caches, links.
- `agent-a5a3601b746680640` (AR5 step 3): the `S11` caches and maps under `ar5/mix/`, links.
- `agent-af8961874e131bda2` (AR3 step 11): the sample `held`, its image stack, the baseline
  cache (`S8_held/`, `S8_held_cap120/`).
- AR5 step 4: `agent-a814d0121e137b01c` (caches on `held`), `agent-a8a1ac8e80f8913cf`. Step 5:
  `agent-afc3f195cb4fffb30`. Item 14: `agent-a012b8f461d1630a1`. C++: `agent-ad44bb90c441eca35`.
- Step 12: `agent-a0ffe667a1a6ee068` (`held2`, `ar5/stallhit/cache`),
  `agent-ab6f021fc3d6bed02` (`mid`).
- `agent-a68538c58665829a1` (AR5 step 1, detached at 29e9b56): `ar5/cache/` and links.
- Step 10: `agent-a9920a0aeaf5d1276` (`ar5/cut/`); 11: `agent-a69be86211b6c2797` (`ar5/patience/`).
- Step 10: `agent-a6c7c67b4a03d1843` (`ar3/keys/cache/`), `agent-a3145f902bf0a61fb` (two reruns).
- AR5 step 8: `agent-abdf64302f4022c0a` (`ar5/order/cache`, MC traces; step 9 links it).
- `agent-a9c88b610fb10e1b6` (step 9): `ar3/start/cache/`, the image stack and samples of `mid`.
- `agent-a635b2f904279bdc9` (step 8): `ar3/reseed/cache/` (logged rerun, `S8` maps), stack links.
- `agent-a38ba6f5a25aa033e` (step 7): its caches and runs under `ar3/big/`, and its report.
- Step 6: `agent-a6db3cabf8ab230e9` (`seed/`, a `boundary/cache/`), `agent-a27b4a932b0dd99b6`.
- Step 5: `agent-a4fd6282107a76186`, `agent-ad88d4ade1ca5a039` (fixes). Steps 3 and 4:
  `agent-ae68d32fb972a8ad4`, `agent-ac609f1c28c619ed0`. More: `grep -n "worktrees/agent-"
  docs/autoresearch/STATE_HISTORY.md`.

## For the owner (open items only; details by `grep` in `STATE_HISTORY.md`, "For the owner")

- Place six pages, or allow me to: AR5 steps 4, 6, 7, 8, 9 and 12. Paths: their Queue rows.
- Your word on PR #43 ("merge PR 43" in the session): AR3 2 to 11, AR5 1 to 3, 5, 7 to 12, AR4 1.
- A choice: which seed step the documentation of the keys recommends. Step 10 added the
  keys and recommends none. My recommendation: the settings of `S8` (`SeedIndexCandidates 20
  5`, `SeedIndexStopHitRatio 0.95`): cheapest on all three samples, wrong counts not worse.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (the commands are in the
  history file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Three questions that block nothing (gradient shapes, your settings, voxels per grain).
- A note: C++ parallel mode keeps the higher confidence for a voxel that two clients fit;
  that rule did harm in the Python test of step 5. A reading of the source, not a test.
- A choice: try test `P2 m=2` as an opt-in key (no gain lost; cuts 17% to 39% of the pass).
- A choice: whether the documentation recommends the two `S11` keys (N21) and the skip key
  of step 7 (N6; `held` is not a clean test for it). All are off; none is recommended now.
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 to N25 and
  D-20261112-N1 to N10, D-20261113-N1 to N10. N15 adds `held` (new seed); N17, N18: run caps
  of 480 minutes. The gradient levels of step 7's sample are my choice. Step 5's result files
  add about 50 MB to git; say so if you want them pruned.

## Bookkeeping

- Next decision ID: D-20261114-N2 (first of a day is N1).
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`. Its plan approvals are the owner's; so is a "finisher" rename.
