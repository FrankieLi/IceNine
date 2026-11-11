# Autoresearch stage: saved state

Read this file first after any restart. It says where the work stands and what to do next. The plan
is `docs/autoresearch/PLAN.md`. The mandate is D-20261109-N1 in `docs/DECISIONS.md`.

## How to restart

1. In the repository root, check `git status` and the branch (`feature/eval-reduction`
   or one of its task branches).
2. Read "Next action" below. For more context: `shmem ask "What is the state of the IceNine
   autoresearch stage?"` (it works again since the evening of 2026-11-09, after the owner updated
   the key; fall back to `shmem query "IceNine autoresearch"` if it fails).
3. Before a new task: `scripts/dev/read_on_fetch.sh --save`, merge
   `origin/feature/python-speed-parity`, rerun with `--mark`, read Ilium's mailbox in
   `docs/ROADMAP.md`.
4. Do every step that is ready, not just one (the owner on 2026-11-10: "you can be a bit more
   aggressive in running"): hand-back, review, fixes, merge, the next plan and the next
   implementer in one wake, and independent jobs side by side within the process limit. After
   each step update this file, commit, push, and save one shmem fact
   (`shmem add --source llm "IceNine autoresearch <step>: <result>"`).

Rules that hold for every step: Python only through `uv run`; one fresh implementer agent per
task in its own worktree (brief from `scripts/dev/make_brief.sh --plan docs/autoresearch/PLAN.md
--task <task> --section "<heading>"`), `scripts/dev/verify_handback.sh`, then `/review-task`, fixes
applied, merge with `scripts/dev/finish_task.sh`; at most 6 heavy processes; no merge to `develop`
without the owner's go; no change to `.claude/` files; no p-values; no extra seeds; findings in
plain English with the plot before the table.

## The wake-up job

A session cron job (id b645c248, created 2026-11-10, every 3 hours at :23) sends the wake-up
prompt while the session is open and idle. Its prompt says to do every step that is ready (the
earlier job 6f6deed3 said one step per wake and was replaced). It expires on 2026-11-17: create
it again on or after 2026-11-16 and write the new date here. It dies with the session. After a
closed session the owner types `/resume` or "continue the autoresearch"; this file is enough to
go on.

## Status (2026-11-10)

- Review page: 24 items decided (D-20261109-N2 to N4). Table on
  `docs/owner_review_2026-10-30.md`, "Decisions of 2026-11-09".
- Plan approved under the mandate (D-20261109-N5).
- Item 6 is done: the review-fix branch was reviewed (APPROVE, 4 required fixes, all applied) and
  merged and pushed (d113694). D-20261109-N6 settles how a boundary crossing is counted.
- AR0 is done: reviewed, 12 fixes applied, merged and pushed (0decc3a). Results below.
- AR1 is done: reviewed (18 fixes applied), merged and pushed (5bd9ecc). Full suite on the task
  branch: 1463 passed, 41 skipped, 1 deselected. Results below.
- D-20261109-N7 settles the review's three questions: the 3 broken voxels count against the
  best-1 arm; AR2 uses the finisher on the best 2 and starts with a holdout; the symmetry
  operators are fixed now for non-cubic crystals.
- The plan's AR2 section now has two steps: the holdout of the finisher arms, then the key.
- AR2 step 1 (the holdout) is done: reviewed (10 fixes applied, D-20261109-N8), merged and pushed
  (86dd905). Full suite on the task branch: 1461 passed, 47 skipped, 1 deselected (6 skips are
  `test_coarse_proxy.py` caches not linked into that worktree). Results below.
- AR2 step 2 (the key `SeedByIndexing`) is done: reviewed (14 fixes applied, D-20261109-N9 and
  N10), merged and pushed (a8c99b7). Full suite at the merge: 1473 passed, 47 skipped, 1
  deselected. MIGRATION_HISTORY entry, ROADMAP row AR and the heads-up to Ilium are written.
  Results below.
- The pre-commit check is tightened on the owner's word "tighten the pre-commit check"
  (D-20261109-N12): reviewed (2 required fixes, applied), merged and pushed (f76f776). Full suite
  at the merge: 1500 passed, 37 skipped, 1 deselected. What it does: the Developer-tooling table
  of `icenine_py/README.md`. What it cannot do: "For the owner" below.
- The AR3 section of the plan is written in full and approved under the mandate (D-20261110-N1,
  2026-11-10). Two steps: painting on the pilot region as a study, then the key and the whole
  layer.
- **AR3 step 1 is done** (2026-11-10): implemented, reviewed (12 required fixes, all applied),
  merged and pushed (8e3098f). Full suite at the merge: 1505 passed, 44 skipped, 1 deselected.
  Numbers: "Results so far". The gate for step 2 is met.
- **PR #42 is merged** into `develop` (66e93fd, 2026-11-10, at branch head 0b3852d) on the
  owner's word: the comment "merge now" on the review page, then "added comments" in the session
  (D-20261110-N3). It holds the work up to AR3 step 1. The branch is kept.
- **The owner corrected the premise of AR3** (2026-11-10, D-20261110-N3): a strain-free grain
  does not have one orientation. Micro-gradients exist, so every voxel in a grain may have its
  own orientation. Painting is a starting approximation and a refinement of each voxel follows.
  "One orientation per grain" is true of the simulated test sample only. The plan's AR3 idea
  paragraph is corrected (93be764).
- **AR3 step 2 is done** (2026-11-10): implemented, reviewed (verdict "request changes", 12
  required fixes plus one added by the main session, all 13 applied by an Opus implementer),
  merged and pushed (741e8fc). Full suite after the fixes: 1515 passed, 44 skipped, 1
  deselected. Numbers: "Results so far". **PR #43** to `develop` is open and waits for the
  owner's word in the session.
- **AR3 step 3 is done** (2026-11-10; plan D-20261110-N4): paint as a start, then refine each
  voxel, on a sample whose grains hold gradients. Implemented, reviewed (verdict "request
  changes", 13 required fixes plus 4 added by the main session, all 17 applied by an Opus
  implementer from the stored records, no new run), merged and pushed (a442712). Full suite at
  the merge: 1520 passed, 44 skipped, 1 deselected. PR #43 now holds steps 2 and 3. Numbers:
  "Results so far".
- **AR3 step 4 is done** (plan 592937a, D-20261110-N5; task `ar3d`): a local fit whose step
  matches the error of its start. Implemented, reviewed (14 required fixes plus 4 of the main
  session, text and summary output only), all 18 applied by an Opus implementer with no new
  run (6e7e471, 60cafe1, ad72aab, 654a09f; 27c72ca by the main session), merged and pushed
  (5ccbd2e). Full suite at the merge: 1520 passed, 44 skipped, 1 deselected. PR #43 now holds
  steps 2, 3 and 4 (title and body updated). Page: `docs/autoresearch/findings_ar3_fine.md`.
- **AR3 step 5 is running** (plan 05577ac, D-20261110-N6; task `ar3e`): after a run, let
  neighbouring grains compete for the boundary voxels. A Sonnet implementer started on
  2026-11-10 at about 18:40 (see "Running jobs").

## Next action

1. AR0: done. Nothing left.
2. AR1: done. Nothing left. Its old worktree `.claude/worktrees/agent-af5119989737208f0` could not be
   removed (it holds two untracked links into the main checkout and 216K of run logs); it is
   harmless and can be removed with `git worktree remove --force`.
3. If an agent died with the session (a task branch exists but no report): look at
   `git log feature/eval-reduction..<task branch>` in its worktree, and start a fresh
   implementer with the same brief (regenerate it with `make_brief.sh`) and the note "continue from
   the commits on the branch".
4. AR2 step 1 (the holdout): done, merged and pushed (86dd905). Its worktree
   `.claude/worktrees/agent-a1bb4d7c15c5bb0d5` could not be removed (untracked link, cache and
   marker files); harmless, removable with `git worktree remove --force`.
5. AR2 step 2 (the key `SeedByIndexing`): done, merged and pushed (a8c99b7). Reports:
   `.claude/reports/20261109-1700-ar2b.md`, `20261109-review-ar2b.md`,
   `20261109-1838-ar2b-fixes.md`. Its worktree `.claude/worktrees/agent-aa30ebdcbce5ec948` could
   not be removed (untracked link and cache); harmless, removable with
   `git worktree remove --force`. Left open and recorded: `key_run_bg.sh` still writes "single
   run" into the label of future runs.
6. AR3 step 2: done, merged and pushed (741e8fc). Reports: `.claude/reports/20261110-0658-ar3b.md`,
   `20261110-0900-review-ar3b.md` (with the main session's additions at its end),
   `20261110-0827-ar3b-fixes.md`. The fixes worktree `.claude/worktrees/agent-a5f0af83ed49c9fc0`
   is locked and was not removed; harmless. PR #43 holds the step.
   **AR3 step 3 (task `ar3c`): done, merged and pushed (a442712).** Reports:
   `.claude/reports/20261110-0852-ar3c.md`, `20261110-1215-review-ar3c.md` (with the main
   session's additions at its end), `20261110-1221-ar3c-fixes.md`. The fixes worktree
   `.claude/worktrees/agent-a46fac4f35d66b8b7` is locked and was not removed; harmless. The
   worktree `.claude/worktrees/agent-ae68d32fb972a8ad4` (detached at a359b87) holds the
   untracked run files of step 3 that step 4 needs: P's `recon.mic` and `records.json`, the
   pairing cache, the image symlink. Do not remove it before step 4 is merged.
   **AR3 step 4 (task `ar3d`)**: done, merged and pushed (5ccbd2e). Reports:
   `.claude/reports/20261110-1313-ar3d.md` (implementer), `20261110-1555-review-ar3d.md`
   (review, fixes 1 to 18), `20261110-1603-ar3d-fixes.md` (fixes). The worktree
   `.claude/worktrees/agent-ac609f1c28c619ed0` holds step 4's untracked run files
   (`fine/runs/P4/maps.npz`, `fine/cache/`, the image symlink); keep it until step 5 is
   merged. From the fixes round, for step 5: among the fits the caller accepted in P4, every
   wrong voxel has an acceptance ratio below 0.95 (113 of 342 wrong at 0.90 to 0.92, 22 of
   614 at 0.92 to 0.95, none above); a description, no threshold picked.
   **AR3 step 5 (task `ar3e`)**: the plan is "Step 5" of the AR3 section (05577ac), approved
   under the mandate (D-20261110-N6). `read_on_fetch.sh` was run and marked (Ilium's branch
   at 347a8a4, nothing to merge). A fresh Sonnet implementer runs in the background in its
   own worktree on `feature/eval-reduction-ar3e` (brief from `scripts/dev/make_brief.sh
   --plan icenine_py/docs/autoresearch/PLAN.md --task ar3e --section "AR3: index the layer
   and paint grains"`, plus notes: step 5 only; the gate after part 0; where the untracked
   maps are; the acceptance ratio as a candidate feature of the skip rule; the 200-voxel
   check first; no truth in any rule; one hand-back at the end). Its report:
   `.claude/reports/20261110-1836-ar3e.md` in its worktree.
   - If the agent died with the session (a branch `feature/eval-reduction-ar3e` with no
     report): see point 3 above.
   - When it hands back: copy the report into the main `.claude/reports/`; place the page
     `docs/autoresearch/findings_ar3_boundary.md` if its write was refused (the text is then
     in the report between `=== FINDINGS PAGE ===` lines); `scripts/dev/verify_handback.sh
     feature/eval-reduction <report> --branch feature/eval-reduction-ar3e --no-tests`; then
     `/review-task` with the extra focus: is the pass compared on one footing, the same pass
     on both sides; does a keep rule or a skip rule use the truth anywhere; are part 2's keys
     listed and the seed's fit described from the code; what the C++ source does with an
     already fitted voxel. Then the fixes (a fresh Opus implementer), `finish_task.sh --yes
     --push`, `gh pr edit 43`, this file, the review page, one shmem fact.
   - After step 5: an opt-in key for the pass if it is a gain; then the seed cost (the grid
     search behind failed index seeds is 25.8 of D4's 35.1 million evaluations). Next
     decision ID: D-20261110-N7 or later.
7. Item 7: **PR #41 is merged** into `develop` (357a2b3, 2026-11-09, at branch head 16c72db), on
   the owner's word in the session (D-20261109-N11). The branch `feature/eval-reduction` is kept.
   **PR #42 is merged** into `develop` (66e93fd, 2026-11-10, at branch head 0b3852d) on the
   owner's word (D-20261110-N3): the item 22 hook change, the tightened pre-commit check, AR3
   step 1. **PR #43** (opened 2026-11-10) holds AR3 step 2 and the plan of step 3; it waits
   for the owner's word.
   Do not merge any PR without the owner's word in the session.
   The owner's review page (https://claude.ai/artifact/A63khAHf7nEx6NBZioSnGW) is updated after
   each step; its source is in the session scratchpad (republish from a read of the page if the
   scratchpad is gone).
8. Item 12 is **not applied**. Its two edits to three `.claude/` files are kept as a patch,
   `icenine_py/docs/autoresearch/item12_owner.patch`. The owner applies and commits it (commands
   under "For the owner"). The working tree is clean, so `finish_task.sh` runs. Do not apply
   the patch in the main checkout again: a dirty tree makes `finish_task.sh` refuse. Never
   `git add -A`, never `git stash`.
9. Merges and the tightened check (f76f776). A merge that git commits by itself and that brings
   a `CLAUDE.md`/`.claude/` change is refused by `.githooks/pre-merge-commit` and left staged.
   Run `git commit` to conclude it. That passes when a fetched remote branch brings the file
   unchanged (for example a merge of `origin/feature/python-speed-parity`). If `git commit`
   refuses too, stop and tell the owner; do not go around it (`git merge --abort` is fine).

Note on implementer agents: since 2026-11-09 the tool refuses a subagent's Write of `.md` files
("subagents must return findings as text"). The agent stops and returns the text. The main
session then writes the findings page and the hand-back report from that text and the committed
result files, and says so in both.

## Queue

Reading the table: one row per step, in order, with its state.

| step | what | state |
|---|---|---|
| item 6 | review and merge the review fixes | **done** 2026-11-09 (d113694) |
| AR0 | stop rules from the stored results | **done** 2026-11-09 (0decc3a) |
| AR1 | recall of the ray-indexing proposer | **done** 2026-11-09 (5bd9ecc) |
| item 7 | open the PR to `develop` (merge waits for the owner) | **done**: PR #41 merged 2026-11-09 (357a2b3) on the owner's word (D-20261109-N11) |
| items 12, 22, 23 | owner said yes to all three (D-20261109-N11) | 22 **done** (fd1721f, hook exception for merges); 23 **done** (rule recorded, Ilium told); 12 is a patch file, **apply and commit wait for the owner** |
| pre-commit tightening | owner: "tighten the pre-commit check" (D-20261109-N12) | **done** 2026-11-09 (f76f776) |
| item 14 + 3 | both keys on E0 and the pilot, timing A/B | not started; needs a quiet machine |
| item 18 | does C++'s clang fuse the clip's multiply-adds | **done** 2026-11-09: yes (below) |
| C++ time | measured C++ time on the pilot region (owner's question on finding 6) | not started |
| AR2 | indexing as the BFS seed search | step 1 (holdout) **done** 2026-11-09 (86dd905); step 2 (the key) **done** 2026-11-09 (a8c99b7) |
| AR3 | index the layer and paint grains | step 1 (paint on the pilot, a study) **done** 2026-11-10 (8e3098f), gate met; step 2 (the key `PaintGrains`, the whole layer) **done** 2026-11-10 (741e8fc), in PR #43; step 3 (paint as a start, then refine, on grains with gradients; stop-rule readings) **done** 2026-11-10 (a442712), in PR #43; step 4 (a local fit whose step matches the error, task `ar3d`) **done** 2026-11-10 (5ccbd2e), in PR #43; step 5 (let neighbouring grains compete for the boundary voxels, task `ar3e`): **running** since 2026-11-10 18:40 (plan 05577ac, D-20261110-N6) |
| AR4 | the learned proposer | waits for AR1 |
| AR5 | the stop rule | waits for AR1; first job: a per-level record of full searches on full-sample voxels (the pilot's 69), then the AR0 stop grid on it |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |

## Running jobs

Step 5 (task `ar3e`) was handed back on 2026-11-10 at about 22:05, status partial: parts 0
and 1 done, run `KE` done, run `K` stopped at its cap with no map. The report is copied to the
main `.claude/reports/20261110-1836-ar3e.md`. The main session placed the page
`docs/autoresearch/findings_ar3_boundary.md` on the task branch (commits 2eda7ec, 3813d69; the
implementer's write was refused). `verify_handback.sh`: 0 FAIL, 3 expected FLAGs. One agent is
running: the Opus reviewer (started at about 22:20). Its report will be
`.claude/reports/20261110-2220-review-ar3e.md`; check with `ls -la` on that path. The task
branch is in the worktree `.claude/worktrees/agent-a4fd6282107a76186` (it holds untracked
`boundary/cache/`, preflight and progress files and `runs/KE/records.json.gz`).

Unreviewed numbers of step 5 (from the implementer's report; do not quote before the review):
gate 412 of 412; wrong voxels before and after the pass with the `hit` keep rule: D4 225 to
74, P4R4 281 to 82, E4 305 to 0, D 203 to 88, PR 170 to 78, D2 153 to 14; the pass costs 3.0
to 3.7 million evaluations; 94% of the voxels are boundary voxels, so the pass is a second
look at nearly everything; skip rule (hit ratio over the best neighbour's at least 0.94)
skips 92% of the tries; KE 5.4 against 13.6 million evaluations; K at its cap had 38.1
million evaluations, 31.07 of them in 90 grid fallbacks. Open after the review, to decide
under the mandate (next ID D-20261110-N7 or later): keep or drop about 47 MB of `tries.npz`
in git; whether to rerun K with a longer cap; whether a sample with large grains is needed
before an opt-in key. Also for the fixes: the owner's term is "local optimization", not
"finisher" (2026-11-10).

Step 4 (merged, 5ccbd2e; page `docs/autoresearch/findings_ar3_fine.md`). The review
recomputed these numbers from the stored records. Scope: the simulated linear-gradient sample
of step 3 (run E4: step 2's sample); one run of each arm, one voxel order. Every loss in part
3 is with the trigger "the cost did not fall" (the stock fallback ran on 3,974 of 53,373
calls in D4), so the plan's fallback was not tested. All 134 newly wrong voxels of P4 plus R4
were wrong after the paint, 107 settled by a fit that passed the 0.9 test, 104 of them from a
grain seeded in another true grain; E4's 170 were accepted at the first neighbour fit (114)
or a revisit (56).

- Part 1, 19,013 painted voxels, share more than 0.1 degrees off: start 2,240; stock fit
  2,230; box 5/8 degree 937; box 5/16 degree 542; the existing CMA fit (first step 0.05
  degrees, at most 60 evaluations) 334. Median evaluations per fit: 662, 398, 1,245 and 61.
  The plan had guessed about 46 per fit; box 5/8 degree is about 1.7 times cheaper than
  the stock fit and box 5/16 costs more. The box 5/32 and the ladder ran on 600 voxels
  only (22,145 and 22,189 evaluations per fit, median). Wrong
  voxels: 87 in every arm. The pick on set A is the CMA fit; on set B it is a gain.
- Part 2: with a fit of 61 evaluations a skip or stop rule saves little (best: 195,627 of
  567,605 on B at about 1% hurt); a second fit helps 73 of 4,768.
- Part 3, whole layer, "CMA first, stock fit when the cost did not fall": D4 35,129,496
  evaluations against 118,350,225, wrong 225 against 203 (no difference); P4 plus R4
  50,892,696 against 90,985,770, wrong 281 against 170 (a loss); E4 13,608,181 against
  33,040,126, wrong 305 against 153 (a loss). The share above 0.1 degrees is a gain in D4
  and in P4 plus R4.

Old agent worktrees under `.claude/worktrees/` are
detached and can stay; `agent-a3fcfcc95d9db6295` holds the untracked run files of the pilot
paint run, `agent-a849170b7d4502001` those of the layer runs (pairing cache, the image symlink).

## Results so far

**AR3 step 3, paint then refine on grains with gradients (done, reviewed, fixes applied).**
Page: `docs/autoresearch/findings_ar3_gradient.md`. Whole layer, 24,570 voxels, 497 grains,
each grain turning linearly from its centre by its level (0, 0.1, 0.2, 0.5 or 1.0 degrees at the
edge); clean images, Q-max 8, FusedCost on, one worker, one run per arm, walls contended.

- Paint (P): 78,142,641 evaluations, 5,369 grain entries for 497 grains (21.9% of the voxels
  are seeds), 170 wrong of 24,570 (0.60 to 0.80%), median error 0.0342 degrees. 172 entries came
  from the grid-search fallback and are 74.8% of the evaluations (58,461,276; the page's 76.6%
  also counted the rejected grid searches). Step 2's paint on the sample without gradients cost
  1,316,345: 59 times less. Without the fallback the factor is 13.9, at 10.2 times the entries.
- Why so many entries: a flood from one orientation stops early. The centre orientation fails
  the 0.9 hit-ratio test on 21.8% of the voxels at level 0.1 and on 93.3% at level 1.0.
- Refinement (R, one local fit from each of the 19,013 painted starts): 12,843,129 evaluations;
  9 voxels came closer by more than 0.05 degrees, 2 were hurt, wrong count unchanged (87). The
  fit ended within 0.02 degrees of its start on 18,977. A second fit (R2) moved 9 of 4,768.
- BFS with indexed seeds (D): 118,350,225 evaluations, 203 wrong (0.72 to 0.95%). P plus R
  against D: 91 fixed, 58 broken, 112 wrong in both; intervals overlap. At level 1.0 paint is
  worse (95 wrong against 78).
- Stop rules read from the records sit at the edges of their grids (they save 98 to 100% of
  the fit evaluations at 0 to 3 of 9,305 voxels hurt). That says the fit moves nothing here; it
  is not a tested stop rule.
- **Why the fit moves nothing (from the stored records and the code; the test is step 4).**
  The records fit a step that is too coarse. It is not a flat cost. The true
  orientation costs less than the painted start on 16,591 of 19,013 voxels (87.3%), and on all
  2,240 that start more than 0.1 degrees off. The fit lowered the cost on 36. It searches a box
  of 5 degrees (LocalOrientationGridRadius 5, MinLocalResolution 0, as in every config file of
  the repository). Its trial steps have a median of 0.70, 1.41, 2.80 and 4.90 degrees in its
  four blocks, and it restarts 4.9 to 16.8 degrees away after a failed block. 18,500 of the
  19,013 fits are exactly four failed blocks (662 evaluations). The C++ source has the same
  logic; C++ was not run.
- The same holds in the BFS (D): 74,739 of its 90,192 fits (82.9%) do not lower the cost. A fit
  helps only when its start is far off: the cost fell on 6 of 8,973 fits that start under 0.05
  degrees off, and on 2,212 of 3,493 that start 0.5 to 1 degree off.
- A gain the page did not report: the share of voxels more than 0.1 degrees off is 9.9 to
  10.7% for paint plus refine against 14.0 to 14.9% for the BFS; in the level-0 grains 3 wrong
  against 17.
- Step 4 tests the remedy: a fit that starts with a small box (plan "Step 4", task `ar3d`).
- Deviations: P ran 123 minutes against a 90-minute stop rule, D 156 against 150.
- Total compute 4.95 hours against a budget of 4.
- Full suite at the merge: 1520 passed, 44 skipped, 1 deselected.

**AR3 step 2, the key `PaintGrains` on the whole layer (done, reviewed, fixes applied).** Page:
`docs/autoresearch/findings_ar3_layer.md`. Whole layer: 24,570 voxels, 497 grains, clean images,
Q-max 8, FusedCost on, one worker, one voxel order; every wall time contended. The truth map has
exactly one orientation per grain by construction, so every reading below holds for such a
sample and for nothing else (the owner's correction, D-20261110-N3).

- Checks: the pilot through the key gives step 1 again (156,766 evaluations, 0 of 2,000 voxels
  differ). The key-on BFS with the fit log is the stored run (2,816,897 evaluations, 0 differ).
- Evaluations: paint 1,316,345 against 33,040,126 for the BFS with indexed seeds (4.0%). In
  paint: seed paths 990,350 (528 paths), scoring 79,979, rejected local fits 246,016 (372 fits,
  none accepted), grid search 0. In the BFS, local fits are 96.6%.
- Wrong voxels (over 1 degree off): paint 14 of 24,570 (0.03 to 0.10%), BFS 153 (0.53 to
  0.73%). Paired: 145 fixed, 6 broken, 8 wrong in both. The intervals do not overlap.
- Direct test on the 363 contested voxels whose passing grains are over 1 degree apart: the
  lowest strict cost names the true grain on 350, the highest hit ratio on 360.
- Grain list: 528 entries for 497 grains; 32 near-duplicates (5.2% of the evaluations).
- Fit log of the BFS: of 48,639 local fits, 56.8% do not move at all, 31.9% move under 0.01
  degrees, 11.1% over 0.1 degrees. The review's recount by voxel: 25 wrong of the 4,747 voxels
  with a fit that moved over 0.1 degrees, 128 wrong of the 19,491 with none.
- Seed record: the lowest-strict-cost candidate is the finisher's end in 513 of 528 seeds; the
  second finisher run is 41.6% of the seed-path evaluations (the review's measurement). On the
  other 15 seeds that candidate's own finisher run ends 14.3 to 55.8 degrees from the truth (0
  of 15 within 1 degree), so "finisher on the lowest-cost candidate only" is not safe as it is.
- Side arm (a fit on each of the 363 contested voxels): 0 fixed, 0 broken, but every kept fit
  ended within 0.01 degrees of its start, so the test could not fail.
- Wall, contended: paint 194.3 s, BFS 2,263.8 s. A reading, not a speed-up.
- The review found no leak of the truth and no error in the main numbers. It found: claims
  scoped too widely (the premise), fit-log counts per fit instead of per voxel, per-voxel
  evaluation records of the key that do not add up to the total (page numbers unaffected), and
  a side arm that could not fail because its fits never moved.
- Full suite after the fixes: 1515 passed, 44 skipped, 1 deselected. Run C's committed
  `records.json` predates the fix of the per-voxel records (it sums to 1,070,857); the page says so.

**AR2 step 2, the key (done, reviewed).** Page: `docs/autoresearch/findings_ar2_key.md`. Pilot
region, clean images, Q-max 8, one run per arm, every time contended (the preflight refused).

- The stored key-off numbers were recomputed from the stored files (no rerun): 37 wrong,
  25,547,935 evaluations, 69 full searches.
- Evaluations: 2,819,906 with the key on against 25,547,935 (about 9.1 times fewer). 64 indexing
  paths at 1,909 each against 69 grid searches at 331,636. Local fits: 2,697,708 on, 2,665,070 off.
- The key's own path (N9: the chain on the seed voxel itself): 0 of 64 indexed seeds wrong (0.0 to
  5.7%) against 2 of 62 for the grid search (0.9 to 11.0%); not shown different. 64 of 64 passed
  the seed test (grid: 7 of 69 failed). 0 fallbacks, so the fallback ran on no real voxel.
- Whole run: 22 wrong of 2,000 against 37 (22 fixed, 7 broken, 15 both). With `FusedCost` on in
  both: 22 against 24 (14 fixed, 12 broken, 10 both). Two key-off runs alone differ by 20 fixed
  and 7 broken. So no gain in wrong voxels is claimed (N9).
- Every flip is an edge-sharing boundary voxel that takes a touching grain's orientation. Three
  fixes are on the key's own path and count for it (N10): 5876, 24334, and 16476 in the
  `FusedCost` pair. The rest count neither way. The trace is an orientation match within 1 degree;
  how an orientation reached a voxel is inferred, since no parent is recorded.
- The key changes which voxels are seeds, so the BFS route changes. Not tested on its own.
- The two key-on runs (`FusedCost` off and on) differ in stored orientation on 2 of 2,000 voxels.
  The review's 637 came from an angle function with a floor.
- Wall time, `FusedCost` on, one worker, contended, one run per arm: 213 s key on, 1633 s key off.
  Pairing 17.6 s once, indexing 0.186 s per voxel. A first reading only; item 14 stays open.
- Voxel 17361: both finisher stages ran and accepted nothing. One `RayGaussNewton` step reaches
  0.1283 (exact truth 0.0949). The strict cost is 0.1556 at 0.02 degrees from the truth and the
  finisher's first step is 0.132 degrees: consistent with a minimum narrower than the first step,
  not tested with a smaller step. A named case for AR5.
- Not tested: the fallback on real voxels; key-off bit-identity on real images (shown with fakes
  and the golden grid); the library's hull and crystal build against the scripts'; the whole
  layer (the pilot pairs rays with the hull of the 2,000-voxel region); few-ray voxels;
  realistic images.

**AR0 (done, reviewed).** Page: `docs/autoresearch/findings_ar0_stopping.md`. E0 is 200 hard
voxels; the pilot is the 2,000-voxel BFS key-off rerun. Lower cost is better.

- After a search, the final cost tells right from wrong with no mistake on E0: 159 right voxels all
  at 0.390 or less, 41 wrong voxels all at 0.673 or more. Read on the same 200 voxels.
- During a search, the best cost after a level does not. On E0 a rule "stop after level k when the
  best cost is at most c" saves at most 1.26% of the evaluations with no voxel broken (finisher
  skipped), 0.37% if the finisher still runs.
- That is E0 only. Level-0 screen and re-cost: 90.7% of the evaluations on E0, 23.7% on full-sample
  voxels (332,046 per search), where 39.2% come after level 0. A stop after level 0 could save far
  more there. Untested: the pilot's `records.json` has no per-level record. First job of AR5.
- On full-sample voxels the screen is 13.2% of the evaluations but 0.9% of the wall time (contended
  timing). So evaluations and wall time do not rank the stages the same way.
- Pilot: at a hit-ratio cut of 0.93, all 1,948 answers at or above are right; all 37 wrong ones
  are below, with 15 right local fits. Picked and read on the same run.

**AR1 (done, reviewed).** Page:
`docs/autoresearch/findings_ar1_ray_indexing.md`. Module: `icenine/ray_indexing.py`. The idea: a
spot on detector 0 paired with its spot on detector 1 gives a ray, a ray gives a scattering vector,
and two vectors of different families give an orientation. So candidates come from votes, not from
a grid.

- Pilot, the 69 voxels that took a full search (images with all grains' spots): the truth is among
  the first 5 candidates for 69 of 69, and the first candidate for 41 of 69. In the other 28 the
  first candidate is the real orientation of a neighbouring grain.
- Evaluations for one full search on those voxels: median 1,138 (strict cost on 5 candidates, then
  the finisher on the best one) against 332,154 for the stored grid search.
- Wrong voxels (more than 1 degree off): 3 of 69 against 4 of 69 stored. It fixes the 4 and breaks
  3. The review traced the 3: the arm picks one candidate before refining, and the neighbour grain
  is ahead only on the unrefined cost; after the finisher the true grain is lower in all 3. So the
  breaks count against this arm. Two of the 4 "fixed" are against an unresolved stored fit (7 of the
  69 are rejected seeds, not plain full searches). With the finisher on the best 2 or the best 3:
  0 of 69, at a mean of 1,926 or 2,546 evaluations; both arms were added after the first result was
  seen, so they need the holdout of AR2 step 1.
- Pairing and indexing time are outside the evaluation count. Indexing one voxel takes about as
  long as 400 strict evaluations (contended).
- In 27 of the 28 voxels whose first candidate is another grain, that grain touches the voxel.
- E0 (the easy case: its images hold only each voxel's own spots): 200 of 200 at every K; 0 wrong
  against 41 stored; median 530 evaluations against 48,724.
- Indexing one voxel: median 0.18 s; one-time pairing of all spots: 17 s. Both contended. No
  wall-time comparison yet.
- The gate for AR2 (recall at 20 of at least 90%) is met.

**AR2 step 1, the holdout (done, reviewed).** Page: `docs/autoresearch/findings_ar2_holdout.md`. 300 voxels outside the pilot region, fixed by
stride before the run: 150 boundary voxels (134 grains) and 150 interior voxels (116 grains).
Arms: strict cost on the first 5 indexing candidates, then the finisher on the best 1, 2 or 3.

- Wrong voxels (more than 1 degree off), both sets: best 1 has 5 of 300, best 2 has 1 of 300,
  best 3 has 1 of 300. All wrong voxels are boundary voxels; interior is 0 of 150 in every arm.
- The plan's rule: best 2 does not have more wrong voxels than best 3, so best 2 stays the main arm
  for the key.
- Mean evaluations per voxel: 1,302 (best 1), 2,046 (best 2), 2,634 (best 3).
- The truth is among the first 5 candidates for 300 of 300, and is the first for 207 of 300
  (boundary 78 of 150, interior 129 of 150).
- Four of best 1's wrong voxels (1687, 8343, 11568, 21535) are the AR1 pattern: the truth is
  second on the unrefined cost and wins after the finisher. They count against best 1.
- Voxel 17361 is wrong in every arm. The first reading (the strict cost prefers the touching
  grain) was wrong. The review found the strict cost at the exact truth is 0.0949, below the
  touching grain's 0.133; the finisher started 0.044 degrees from the truth and did not move
  (0.162 before and after). So it is a failure of the finisher and counts against all three arms
  (D-20261109-N8). Why the finisher did not move is not known: a named test case for AR2 step 2
  and AR5.
- The review also found: the error angle has a floor of a few thousandths of a degree (the
  function clipped; wrong counts are not affected; fixed in the library with atan2), and 74 of the 150
  boundary voxels touch another grain only at a corner (all 5 wrong voxels share an edge).
- Not covered: every holdout voxel has at least 30 explained rays, so few-ray voxels are untested
  (target of AR4). Timing was contended (3 processes). No stored grid search for these voxels.

Note on the `.md` refusal: the AR0 fix agent's Write of a new report file was allowed the same
afternoon. The refusal is not consistent; the rule in "Next action" still holds when it happens.

**Item 18.** The C++ build on this machine fuses the clip's multiply-add: Apple clang 21.0.0 (arm64) at
`-O3 -DNDEBUG` compiles `p0 + slope * (m - p0)` of `SutherlandHodgman.h` to one `fmadd` (single
precision); with `-ffp-contract=off` it is `fmul` then `fadd`. Test: `scripts/autoresearch/item18/clip_fma.cpp`.
Not checked: GCC, x86.

**AR3 step 1 (paint grains on the pilot; reviewed).** Page: `docs/autoresearch/findings_ar3_paint.md`.
Pilot region, 2,000 voxels, clean images, Q-max 8, FusedCost on, one fixed voxel order.
- Evaluations: paint 156,766; key-on BFS 2,819,906 (paint is 5.6%); key-off BFS 25,547,935
  (0.6%). By stage in paint: seed paths 119,033 (62 seeds), scoring 6,697, rejected local fits
  31,036 (48 fits, none accepted), grid search 0.
- Wrong voxels: paint 0 of 2,000 (0.0 to 0.2%); key on 22; key off 37. Against key on: 22
  fixed, 0 broken. Two key-off runs alone differ by 20 one way and 7 the other, so the count
  alone is not called a gain.
- The direct test carries the result: on the 44 contested voxels whose passing grains are over
  1 degree apart, the lowest strict cost names the true grain in 44 of 44 (92 to 100%). The
  other 82 contested voxels are passed only by two entries of one grain; the test cannot fail
  there.
- Of the 22 voxels wrong in key on: on 19 the wrong grain passes in paint too and loses on cost;
  on 3 it does not pass without a fit.
- Passes outside the true grain: 44, all within 1 voxel side. 62 entries for 60 grains.
- Accuracy without a fit: median 0.0169 degrees against 0.0159 for key on (62 distinct values,
  not 1,978 independent ones).
- Wall, contended, one worker: paint 37.9 s (17.2 s of it the one-time pairing); key on 224.3 s
  in the same session. A first reading of an unbatched loop.
- Review: no leak of the truth into any decision; same inputs and tests as the stored key-on
  run; evaluations counted on the same footing.
- Limits: the truth has exactly one orientation per grain, so the premise is true by
  construction. Not tested: the whole layer, realistic images, another voxel order, the
  fallback and the fit-settled branch (they ended no voxel).
- Found on the way, not traced: the key-on BFS map has only 64 distinct orientations (key off
  63), so the BFS local fits nearly always return their start.

## For the owner

- Still yours: item 12 (apply and commit, below), Ilium's plan approvals, and the note on what
  the pre-commit check cannot do (below).
- 2026-11-10: PR #42 was merged on your comment "merge now" on the review page together with
  your session message "added comments" (D-20261110-N3). If the session message was not meant
  as the go, say so and the merge is reverted by a PR.
- 2026-11-10: **PR #43** (https://github.com/FrankieLi/IceNine/pull/43) holds AR3 steps 2,
  3 and 4 and the plan of step 5. It is merged only on your word in the session ("merge PR 43").
- 2026-11-10, a question that blocks nothing: which local-fit settings do you use on real data?
  Every config file here has `LocalOrientationGridRadius 5` and `MinLocalResolution 0`, which
  gives first trial steps of about 0.7 degrees. Step 3 found that such a fit almost never
  improves a start that is under 0.1 degrees off. Your answer names the fair baseline of step 4.
- 2026-11-10, a question that blocks nothing: step 3 simulates grains whose orientation turns
  linearly with position, by 0 to 1 degree from the centre to the edge. Which sizes and shapes
  of gradient are typical in the samples you care about (smooth, or with sharp sub-grain
  boundaries)?
- 2026-11-10, decided under the mandate (D-20261110-N4), yours to overrule: step 2's whole-layer
  result files (about 44 MB raw) stay in git; later tasks commit only what their summary reads.
- 2026-11-09 (evening): the owner commented on the review page (merge PR #41; yes to items 12,
  22, 23; the embedder is no longer run) and then said in the session "Go ahead and execute
  according to the comments". That session message is the authority (D-20261109-N11); the page
  comments alone were not acted on. Done: PR #41 merged (357a2b3); item 22 implemented, reviewed
  (1 required fix, applied) and merged (fd1721f; suite 1489 passed, 37 skipped, 1 deselected);
  item 23 recorded and Ilium told (ROADMAP mailbox).
- Item 12, left for the owner. The edits are a patch file, not in the working tree (uncommitted
  `.claude/` files made `finish_task.sh` refuse every merge, so they were taken out again on
  2026-11-09). From the repository root, on `feature/eval-reduction`:
  `git apply icenine_py/docs/autoresearch/item12_owner.patch`, then
  `git add .claude/skills/implement-plan/SKILL.md .claude/skills/apply-fixes/SKILL.md
  .claude/hooks/run-tests-on-change.sh`, then `ALLOW_CLAUDE_CONFIG=1 git commit -m "Item 12:
  show make_brief notes; drop the extra -q in the test hook"`. The agents cannot do this: the
  Bash hook refuses an agent that sets an `ALLOW_*` override. An earlier offer on the review
  page to commit with the override was wrong and is withdrawn.
- 2026-11-09 (night): the owner said in the session "tighten the pre-commit check"
  (D-20261109-N12). Done and merged (f76f776). The three gaps noted before are closed: a merge
  that git commits by itself is now checked (`.githooks/pre-merge-commit`); a deletion, a rename
  away, a change of type and a submodule entry are refused; any letter case is guarded. The
  merge exception of item 22 is narrower: the merged commit must be on a remote-tracking branch
  and not already part of `HEAD`, and the merged branch must have changed the file since every
  merge base. The reviewer required 2 fixes (submodule entries; a merge head already in `HEAD`);
  both are applied. One further narrowing came from a reviewer note (a merge could bring back an
  older copy that the merged branch merely still had).
- What the pre-commit check cannot do. It is a git hook, so it runs only for `git commit` and
  for merges. The reviewer tried each route below and each one changes a guarded file with no
  check. None is new. An agent taking one would be breaking the rule "agents never bypass the
  hooks" (D-20261028-N3); the hook guards against accidents, not against that.
  (a) `git commit-tree` with `git update-ref`. (b) A fast-forward merge of a commit made that
  way. (c) `git cherry-pick`, `git revert`, `git rebase`, `git am` of such a commit.
  (d) `git merge --no-verify` (the Bash hook blocks `--no-verify` only for `git commit`).
  (e) `git -c core.hooksPath=... commit`. (f) Making a local commit look fetched with
  `git update-ref refs/remotes/...` or a remote that points at the clone itself.
  Routes (d) to (f) could be refused in `.claude/hooks/check-bash-command.sh`. That file is the
  owner's; no agent has changed it. A stronger backstop is on GitHub: a rule on `develop` that
  needs the owner's review for `CLAUDE.md` and `.claude/**` (a CODEOWNERS file).
- Small, known, not fixed: `precommit_check.py` stops with a traceback (so it refuses) when a
  submodule entry is staged anywhere outside `.claude/`. The repository has no submodules.
- The `shmem` warning about the local embedder (ollama, port 11434) is expected: the owner said
  on 2026-11-09 that the embedder is no longer run. Do not report it.
