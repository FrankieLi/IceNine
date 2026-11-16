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
| open PR | #43 (AR3 steps 2 to 11, AR5 steps 1 to 3, 5, 7 to 12, AR4 step 1, items 14 + 3, C++ time), https://github.com/FrankieLi/IceNine/pull/43 |

## How to restart

1. In the repository root: `git status`, branch `feature/eval-reduction` (or a task branch).
2. Do "Next action" below: every ready step in turn, independent jobs side by side (owner,
   2026-11-10). After each step: rewrite this file, commit, push, one shmem fact.
3. A detail that is not here: `shmem ask "<question>" | head -30`, then `grep` in the page.

## Next action
Runbook: "A step, from plan to merge", points 5 to 9. Merge the parent into a task branch first.
1. AR5 step 4 (`ar5d`; N16, N18, N20, N23): `S11v2` on `held` 37.4% fewer evaluations than `S8`,
   69 wrong; approved with three fixes; fix 3 done (fac79f9). BLOCKED on fix 1 (the page; agent
   writes refused): the owner places it or allows it. Then: sync check, commit, fix 2, merge, close.
2. AR5 step 6 (`ar5f`; N3): fixed, waits. Steps 9 to 12: merged (10, 11 failed studies, no key).
   Step 12 (key `SeedIndexMCStallHit`, off): merged at fe3d277; page waits. Step 13 (`ar5m`):
   BLOCKED FOR THE OWNER (Queue).
3. Items 14 + 3 (`item14`): MERGED at 8f8ca20 (smooth key fails rule R3; not recommended). Page
   waits: `.claude/reports/20261114-item14-findings_page_draft-v2.md` goes to
   `icenine_py/docs/autoresearch/findings_item14.md`.
   C++ time (`cpptime`): MERGED at 003cc5d (C1 1,236.87 s, C2 1,231.87 s, two busy processes).
   Not a speed comparison: the C++ run left 416 of 2,000 pilot voxels wrong, Python 24; 393
   refit voxels were never fitted again, all 1,921 local optimizations changed nothing.
   Flagged for the owner. Page waits: `.claude/reports/20261114-cpptime-findings_cpp_time-v2.md`
   goes to `icenine_py/docs/autoresearch/findings_cpp_time.md`. 16 pages wait (+ar5v, cppchk).
   NOW: (a) step 22 MERGED 566eea3f; v58 live (15:34, 16 Nov); v59 after C5, slices.
   (b) OWNER TASK part 1: fix 2 (973b57b5): C5 56 wrong (C4 311, C3 416); review af0fc74035f694d0f;
   part 2 `fine1`: slices run (CS, J5S, PS); then SendMessage a22702f36dfefa989 to project; caps.
   (c) restart distance, C++ config: done (reports 20261116-*); next: local-opt. schedule step.

## Running jobs
- Step 14 key WAITS (held2). Steps 15 to 21 merged (caches: `J5r` `agent-aeba15b7c3850b007`,
  step 20 `agent-a8b2d0066621fd0f8`, step 21 `agent-a8c1c16b70cf43e39`); v58 live.
  Step 22 (`ar5v`; cache in `agent-aba6e1bdb95dff72a`): MERGED 566eea3f, key `SeedIndexMCDefer`.
  Owner task: part 1 trace done (ae912cf9), `Reset` fix runs. Part 2 `fine7` frozen (ac76bdd6);
  slices restarted 15:38, 16 Nov with `--allow-contended` (first try refused: machine busy).

## Queue
Reading the table: one row per step that is not finished, in order. Finished steps (items 6,
7, 18, 22, 23, AR0 to AR2, AR3 1 to 11, AR5 1 to 3, 5, 7 to 12, AR4 1): `STATE_HISTORY.md`.

| step | what | state |
|---|---|---|
| OWNER TASK 16 Nov | (1) check the C++ pilot result for a bug (393 refits never fitted, 1,921 local optimizations returned their start, same inputs, scoring); (2) new finer layer, about 98,000 voxels: stock Python, stock C++, J5S; slice timing first, ask the owner if a run projects over 48 hours | FIRST after step 22, ahead of all else. Spec: PLAN, two sections 'Owner task of 16 Nov'; part 1 running; `shmem query "owner task C++ bug finer layer"` |
| AR5 step 4 | `S10`, `S11`, `S11v2` on the held-out sample (task `ar5d`) | waits for the owner (page): `.claude/reports/20261111-2330-ar5d-phaseB-page.md` (main checkout) goes to `icenine_py/docs/autoresearch/findings_ar5_held.md` in worktree `agent-a8a1ac8e80f8913cf` |
| AR5 step 6 | `S11n`, `S11h` on `mid` and `grad` (task `ar5f`) | reviewed and fixed; waits for the owner (page) and for step 4's merge. Branch tip d643eff, held by the fixes worktree `agent-a9420a98460078e81`; the first worktree `agent-a89458addd73332d9` (detached, 4c18740) holds the run cache `ar5/ver/cache/`: keep both. Opus review (`.claude/reports/20261112-0510-review-ar5f.md`): every headline number recomputed, none differs; 11 wording fixes, no rerun (N3). Fixes round done (7f272ce, 18260d5, d643eff; report `.claude/reports/20261112-0517-ar5f-fixes.md` in the fixes worktree; hand-back check: only the expected format flag; no number changed, `predictions[4].met` now false; 14 tests pass). Numbers: on `mid` `S11n` 14,105,137 evaluations (12.3% fewer than `S11`, 6 seeds kept, 0 wrong), `S11h` 15,680,025 (2.5%, 1 kept); both pass; "plain does as well" not met on `mid` (20.4% of the saving; 62.2% on `grad`, in-sample); six of seven predictions met. PAGE: not in the repo (Write refused). Corrected text: `.claude/reports/20261112-ar5f-page-v2.md` in the fixes worktree; the owner places it at `icenine_py/docs/autoresearch/findings_ar5_ver.md` there (or allows it). Then: sync check, commit, merge after step 4 (`finish_task.sh` in the fixes worktree), close (pareto points, page, PR #43); then a decision on an opt-in key for the `N1` rule |
| AR5 step 7 | opt-in key `BoundaryPassSkipHitRatio` (task `ar5g`) | MERGED at 50c3b9c (N4, N5, N6), see "Last finished step". Open: only its page, for the owner to place: text `.claude/reports/20261112-ar5g-page-v2.md` (main checkout; the merge removed the fixes worktree, the text was recovered exactly from the recorded write and a table sync) goes to `icenine_py/docs/autoresearch/findings_ar5_gate.md` on `feature/eval-reduction`, then a sync check and a commit. Keep worktree `agent-a331e7882919e74aa` (stored passes, `ar5/gate/`). Step 6's page text and reports are now also copied to the main `.claude/reports/` |
| AR5 step 8 | sweep order of the seed step (task `ar5h`) | MERGED at 949dc61 as a study; failed rule 5 (N9), no key; see "Last finished step". Open: only its page for the owner: `.claude/reports/20261112-ar5h-page.md` (main checkout) goes to `icenine_py/docs/autoresearch/findings_ar5_order.md`. Keep worktree `agent-abdf64302f4022c0a` (detached): the step 8 cache `ar5/order/cache` that step 9 links, and the MC traces. Reports also in the main `.claude/reports/`: `20261112-1419-ar5h.md`, `-corrected.md`, `20261112-2215-ar5h-fixes.md` |
| AR5 step 9 | rule `G1`, key `SeedRescueCMA` (task `ar5i`) | MERGED at 28b75fa. Open: its page for the owner: `.claude/reports/20261113-ar5i-page.md` goes to `icenine_py/docs/autoresearch/findings_ar5_rescue.md`. Keep worktrees `agent-a8d20b3821621fd41`, `agent-ab6ee3e8ef1058af1` (`ar5/rescue/cache`) |
| AR5 step 12 | rule PH(3000), key `SeedIndexMCStallHit` (task `ar5l`) | MERGED at fe3d277 (13-N7, N8, N10), see "Last finished step". Open: only its page for the owner: `.claude/reports/20261113-ar5l-page-v2.md` (main checkout) goes to `icenine_py/docs/autoresearch/findings_ar5_stallhit.md`. Keep worktrees `agent-a0ffe667a1a6ee068` (held2 sample, runs; `passes_held2/D4_H/tries.npz` stays unopened until step 13's freeze) and `agent-ab6f021fc3d6bed02` (LibH, Hl runs on `mid`). Reports in the main `.claude/reports/`: `20261113-1100-ar5l.md`, `-1600-review-ar5l.md`, `-1630-ar5l-B.md`, `-1730-review-ar5l-B.md`, `-1800-ar5l-fixes2.md` |
| AR5 step 13 | arm HK = H + pass skip key 0.95, best current recipe; clean test on held2 (task `ar5m`) | D-20261113-N9 (6d13c12), PLAN "AR5 step 13". BLOCKED FOR THE OWNER since 17:05, 13 Nov (agent `afef88a15cf24bade`, worktree `.claude/worktrees/agent-afef88a15cf24bade`, tip be8d2fc, report `.claude/reports/20261113-1700-ar5m.md` there; no freeze commit; held2 tries unopened). (1) The agent's run of `recipe_freeze.py` was refused by the permission classifier; the main session does NOT run it for the agent; the owner runs it or allows it. (2) Gate g2 cannot hold as written: the library pass uses one generator, H's stored passes used aligned draws per try (mid: same skipped set 13,782, same 15,675 kept tries and 956,175 evaluations, only 105 tries identical, wrong 37 live vs 35 replay). Recommended amendment (needs a decision before any held2 read): g2 compares a script-layer aligned-draw pass with the skip against the replay; the library pass is compared on skipped set, counts and own hit ratios. Phase B's precondition (step 12 merged) is now met. |
| item 5 | full-sample baseline run, `FusedCost` on | item 14 merged; needs planning |
| tooling | `audit_numbers.py`: thousands separators; false misses in wide table rows | small |

## Last finished step

C++ time (`cpptime`): merged at 003cc5d, 14 Nov; suite 2030 passed, 47 skipped, 1 deselected.
- Walls: C1 1,236.87 s, C2 1,231.87 s (two busy processes); not a speed comparison to Python.
- C++ left 416 of 2,000 pilot voxels wrong (Python 24); 393 refit voxels never fitted again.

## Untracked data that later steps need

Do not remove these worktrees (all under `.claude/worktrees/`). Link the data, do not commit it.

- `agent-acc78ac2a483ed86f` (AR5 step 2), `agent-ab349f0b6193c2ca1` (AR4 step 1): caches, links.
- `agent-a5a3601b746680640` (AR5 step 3): the `S11` caches and maps under `ar5/mix/`, links.
- `agent-af8961874e131bda2` (AR3 step 11): sample `held`, image stack, baseline cache.
- AR5 step 4: `agent-a814d0121e137b01c` (caches on `held`), `agent-a8a1ac8e80f8913cf`. Step 5:
  `agent-afc3f195cb4fffb30`. Item 14: `agent-a012b8f461d1630a1`. (C++ time: no worktree left.)
- Step 12: `agent-a0ffe667a1a6ee068` (`held2`, `ar5/stallhit/cache`), `agent-ab6f021fc3d6bed02`
  (`mid`). Step 1: `agent-a68538c58665829a1` (detached at 29e9b56): `ar5/cache/` and links.
- Step 10: `agent-a9920a0aeaf5d1276` (`ar5/cut/`), `agent-a6c7c67b4a03d1843` (`ar3/keys/cache/`),
  `agent-a3145f902bf0a61fb` (two reruns). Step 11: `agent-a69be86211b6c2797` (`ar5/patience/`).
- Step 8: `agent-abdf64302f4022c0a` (`ar5/order/cache`, MC traces; step 9 links it),
  `agent-a635b2f904279bdc9` (`ar3/reseed/cache/`). Step 9: `agent-a9c88b610fb10e1b6`
  (`ar3/start/cache/`, image stack and samples of `mid`).
- Older (steps 3 to 7: moved to `STATE_HISTORY.md`, "Worktrees kept, steps 3 to 7"): `grep` it.

## For the owner (open items only; details by `grep` in `STATE_HISTORY.md`, "For the owner")

- Place eight pages, or allow me to: AR5 steps 4, 6, 7, 8, 9, 12, item 14, C++ time (Queue
  rows; the last two: "Next action" 3).
- Your word on PR #43 ("merge PR 43" in the session): AR3 2 to 11, AR5 1 to 3, 5, 7 to 12, AR4 1,
  items 14 + 3, C++ time.
- A choice: which seed step the documentation of the keys recommends (step 10 added the keys,
  recommends none). Mine: `S8` settings (`SeedIndexCandidates 20 5`, `SeedIndexStopHitRatio
  0.95`): cheapest on all three samples, wrong counts not worse.
- Item 12: apply and commit `docs/autoresearch/item12_owner.patch` (commands in the history
  file). Agents cannot: it needs your `ALLOW_CLAUDE_CONFIG=1`.
- Three questions that block nothing (gradient shapes, your settings, voxels per grain).
- A note: C++ parallel mode keeps the higher confidence for a voxel that two clients fit; that
  rule did harm in the Python test of step 5. A reading of the source, not a test.
- A choice: try test `P2 m=2` as an opt-in key (no gain lost; cuts 17% to 39% of the pass).
- A choice: whether the documentation recommends the two `S11` keys (N21) and the skip key of
  step 7 (N6; `held` is not a clean test for it). All are off; none is recommended now.
- Decided under the mandate, yours to overrule: D-20261110-N4 to N8, D-20261111-N1 to N25,
  D-20261112-N1 to N10, D-20261113-N1 to N10, D-20261114-N1 to N3. N15 adds `held` (new seed);
  N17, N18: run caps of 480 minutes. The gradient levels of step 7's sample are my choice.
  Step 5's result files add about 50 MB to git; say so if you want them pruned.

## Bookkeeping
- Next decision ID: D-20261116-N6 (first of a day is N1).
- Wake-up job: session cron id 24092880, every 3 hours at :53, created 2026-11-16, expires
  2026-11-23. Create it again on or after 2026-11-22 and write the new id and date here.
- Other agent: Ilium (Track S, `origin/feature/python-speed-parity`, last read at 347a8a4);
  mailbox in `docs/ROADMAP.md`. Its plan approvals are the owner's; so is a "finisher" rename.
