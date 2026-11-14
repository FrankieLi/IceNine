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
| open PR | #43 (AR3 steps 2 to 11, AR5 steps 1 to 3, 5, 7 to 11, AR4 step 1), https://github.com/FrankieLi/IceNine/pull/43 |

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
2. AR5 step 6 (`ar5f`; N3): fixed, waits. Steps 9 to 11: merged (10, 11 failed studies, no key).
   Step 12 (`ar5l`; 13-N7, N8): rule PH(3000) PASSES all rules on 5 samples; review APPROVE.
   Fixes + key phase (agent `ab6f021fc3d6bed02`, Opus) since 16:15, 13 Nov; then merge, v49.
3. Items 14 + 3 (task `item14`; N25): E0 arms done. The run stopped at 00:01 on 2026-11-12,
   exit 3: the preflight refused (`cpptools`, pid 18415). Once the owner frees the machine:
   rerun the command below (E0 resumes; warm-up 25 min). After exit 0: phase B by
   `SendMessage`. While timing arms run: no merge, no test suite, no other run.

## Running jobs
- Items 14 + 3: `caffeinate -i bash scripts/autoresearch/item14/run_item14.sh` in `icenine_py/`
  of worktree `agent-a012b8f461d1630a1` (stopped, exit 3, 00:01, 2026-11-12); marker `item14.done`.
  Branch `ar5d` is held by worktree `agent-a8a1ac8e80f8913cf` (tip fac79f9).
- Page v48 live; v49 after step 12. Step 12 fixes + key agent runs since 16:15, 13 Nov.

## Queue
Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0 to AR2, AR3 1 to 11, AR5 1 to 3, 5, 7 to 11, AR4 1): `STATE_HISTORY.md`.

| step | what | state |
|---|---|---|
| AR5 step 4 | `S10`, `S11`, `S11v2` on the held-out sample (task `ar5d`) | waits for the owner (page): `.claude/reports/20261111-2330-ar5d-phaseB-page.md` (main checkout) goes to `icenine_py/docs/autoresearch/findings_ar5_held.md` in worktree `agent-a8a1ac8e80f8913cf` |
| AR5 step 6 | `S11n`, `S11h` on `mid` and `grad` (task `ar5f`) | reviewed and fixed; waits for the owner (page) and for step 4's merge. Branch tip d643eff, held by the fixes worktree `agent-a9420a98460078e81`; the first worktree `agent-a89458addd73332d9` (detached, 4c18740) holds the run cache `ar5/ver/cache/`: keep both. Opus review (`.claude/reports/20261112-0510-review-ar5f.md`): every headline number recomputed, none differs; 11 wording fixes, no rerun (N3). Fixes round done (7f272ce, 18260d5, d643eff; report `.claude/reports/20261112-0517-ar5f-fixes.md` in the fixes worktree; hand-back check: only the expected format flag; no number changed, `predictions[4].met` now false; 14 tests pass). Numbers: on `mid` `S11n` 14,105,137 evaluations (12.3% fewer than `S11`, 6 seeds kept, 0 wrong), `S11h` 15,680,025 (2.5%, 1 kept); both pass; "plain does as well" not met on `mid` (20.4% of the saving; 62.2% on `grad`, in-sample); six of seven predictions met. PAGE: not in the repo (Write refused). Corrected text: `.claude/reports/20261112-ar5f-page-v2.md` in the fixes worktree; the owner places it at `icenine_py/docs/autoresearch/findings_ar5_ver.md` there (or allows it). Then: sync check, commit, merge after step 4 (`finish_task.sh` in the fixes worktree), close (pareto points, page, PR #43); then a decision on an opt-in key for the `N1` rule |
| items 14 + 3 | both speed keys on E0 and the pilot, timing A/B (task `item14`) | E0 done; timing waits for a quiet machine |
| C++ time | measured C++ time on the pilot region (MPI mode, own config) | D-20261112-N1, task `cpptime`: phase A done (agent and worktree `ad44bb90c441eca35`, tip 5d6cec2, report in that worktree's `.claude/reports/20261112-0012-cpptime.md`). Run, after item 14's timing and on a quiet machine: `caffeinate -i bash scripts/autoresearch/cpptime/run_cpp_time.sh` in its `icenine_py/` (background command); marker `benchmarks/autoresearch/cpptime/cache/cpptime.done`; then phase B by `SendMessage`. R4 limit 3.5 (own load 2.0: the server rank spins). `lldb` cannot run (developer mode off): evaluations missing unless the owner enables it or accepts a print in `Src/` |
| AR5 step 7 | opt-in key `BoundaryPassSkipHitRatio` (task `ar5g`) | MERGED at 50c3b9c (N4, N5, N6), see "Last finished step". Open: only its page, for the owner to place: text `.claude/reports/20261112-ar5g-page-v2.md` (main checkout; the merge removed the fixes worktree, the text was recovered exactly from the recorded write and a table sync) goes to `icenine_py/docs/autoresearch/findings_ar5_gate.md` on `feature/eval-reduction`, then a sync check and a commit. Keep worktree `agent-a331e7882919e74aa` (stored passes, `ar5/gate/`). Step 6's page text and reports are now also copied to the main `.claude/reports/` |
| AR5 step 8 | sweep order of the seed step (task `ar5h`) | MERGED at 949dc61 as a study; failed rule 5 (N9), no key; see "Last finished step". Open: only its page for the owner: `.claude/reports/20261112-ar5h-page.md` (main checkout) goes to `icenine_py/docs/autoresearch/findings_ar5_order.md`. Keep worktree `agent-abdf64302f4022c0a` (detached): the step 8 cache `ar5/order/cache` that step 9 links, and the MC traces. Reports also in the main `.claude/reports/`: `20261112-1419-ar5h.md`, `-corrected.md`, `20261112-2215-ar5h-fixes.md` |
| AR5 step 9 | rule `G1`, key `SeedRescueCMA` (task `ar5i`) | MERGED at 28b75fa. Open: its page for the owner: `.claude/reports/20261113-ar5i-page.md` goes to `icenine_py/docs/autoresearch/findings_ar5_rescue.md`. Keep worktrees `agent-a8d20b3821621fd41`, `agent-ab6ee3e8ef1058af1` (`ar5/rescue/cache`) |
| AR5 step 12 | rule PH(3000): stop a stalled seed-step MC refinement only if its best is at hit ratio 0.95 or more; new test sample `held2` (task `ar5l`) | D-20261113-N7, N8 (efe12c6: passes; 6 fixes; key `SeedIndexMCStallHit` approved). Review `.claude/reports/20261113-1600-review-ar5l.md`. ROADMAP heads-up posted (e646244). FIXES + PHASE B agent `ab6f021fc3d6bed02` (own worktree, report `.claude/reports/20261113-1630-ar5l-B.md`); at hand-back: verify, short Opus review of the key, copy reports to main, merge from its worktree, close, page v49. Next decision N9. OLD: GATE PASSED, tip bc6c650, freeze 9e315e2, report `.claude/reports/20261113-1100-ar5l.md` there, parent config set. RUNS DETACHED 11:02 on 13 Nov (out: scratchpad `stallhit_run_all.out`). Check: `cat <worktree>/icenine_py/benchmarks/autoresearch/ar5/stallhit/cache/logs/stallhit_run_all.done` (first line `exit=0`; absent while running; `pgrep -f autoresearch/ar5/`); if it died: `scripts/autoresearch/ar5/stallhit_run_all.sh --resume` detached from that `icenine_py`. Implementer id: `a0ffe667a1a6ee068`, worktree `.claude/worktrees/agent-a0ffe667a1a6ee068`. Phase A items 1 to 4 (freeze, gate replay) since about 10:35, 13 Nov; hands back "READY TO RUN". THEN: hand-back check; start `scripts/autoresearch/ar5/stallhit_run_all.sh` detached from the implementer's worktree `icenine_py/` (about 3.5 hours; marker `benchmarks/autoresearch/ar5/stallhit/cache/logs/stallhit_run_all.done`, first line `exit=0`). After: stats phase, hand-back check, Opus review, decision, fixes, key `SeedIndexMCStallHit` only if it passes (ROADMAP heads-up before the merge), merge, close, page v49 |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |
| tooling | `audit_numbers.py`: thousands separators; false misses in wide table rows | small |

## Last finished step

AR5 step 11 (`ar5k`): merged as a study at 2939db9, 2026-11-13; suite 1880 passed. Rule P(3000)
failed rule 5 on `mid` by one voxel (net 6, limit 5); no key. Held 5.7% fewer evaluations
than G, mid 10.6%. Lesson: 597 of 1,181 cut refinements had a best under 0.95 (see step 12).

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-acc78ac2a483ed86f` (AR5 step 2), `agent-ab349f0b6193c2ca1` (AR4 step 1): caches, links.
- `agent-a5a3601b746680640` (AR5 step 3): the `S11` caches and maps under `ar5/mix/`, links.
- `agent-af8961874e131bda2` (AR3 step 11): the sample `held`, its image stack, the baseline
  cache (`S8_held/`, `S8_held_cap120/`).
- AR5 step 4: `agent-a814d0121e137b01c` (caches on `held`), `agent-a8a1ac8e80f8913cf`. Step 5:
  `agent-afc3f195cb4fffb30`. Item 14: `agent-a012b8f461d1630a1`. C++: `agent-ad44bb90c441eca35`.
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

- Free the machine for the item 14 timing: `cpptools` (VS Code C/C++ extension, pid 18415)
  has used a whole core for 2 days; the timing preflight refuses while it does. Yours to close.
- Place five pages, or allow me to: AR5 steps 4, 6, 7, 8 and 9. Paths: their Queue rows.
- Your word on PR #43 ("merge PR 43" in the session): AR3 2 to 11, AR5 1 to 3, 5, 7 to 11, AR4 1.
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
  D-20261112-N1 to N10, D-20261113-N1 to N7. N15 adds `held` (new seed); N17, N18: run caps
  of 480 minutes. The gradient levels of step 7's sample are my choice. Step 5's result files
  add about 50 MB to git; say so if you want them pruned.

## Bookkeeping

- Next decision ID: D-20261113-N8 (first of a day is N1).
- Wake-up job: session cron id b645c248, every 3 hours at :23, created 2026-11-10, expires
  2026-11-17. Create it again on or after 2026-11-16 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`. Its plan approvals are the owner's; so is a "finisher" rename.
