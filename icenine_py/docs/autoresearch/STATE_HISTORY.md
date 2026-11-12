# Autoresearch stage: state history

**Archive. Not read on a wake.** This is the state file as it stood on 2026-11-10 at the merge
of AR3 step 5 (b69a66f), kept whole. The live state is `STATE.md`; the procedures are in
`RUNBOOK.md`. Look things up here with `grep -n`; do not read it whole. Notes of later
finished steps are added at the end, under a dated heading.

The text below is unchanged. Where it says "this file" it means the old state file.

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
- **AR3 step 5 is in its fixes round** (plan 05577ac, D-20261110-N6; task `ar3e`): after a
  run, let neighbouring grains compete for the boundary voxels. Handed back partial, reviewed
  (request changes, 15 fixes), an Opus implementer applies the fixes since 2026-11-10 at
  about 22:40 (D-20261110-N7; see "Running jobs").
- **AR3 step 6 is planned** (plan c88fdd8, D-20261110-N8; task `ar3f`): the seed cost. Why
  the index seed fails where the grid search succeeds, and what can replace the grid
  fallback. It starts after step 5 is merged.

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
| AR3 | index the layer and paint grains | step 1 (paint on the pilot, a study) **done** 2026-11-10 (8e3098f), gate met; step 2 (the key `PaintGrains`, the whole layer) **done** 2026-11-10 (741e8fc), in PR #43; step 3 (paint as a start, then refine, on grains with gradients; stop-rule readings) **done** 2026-11-10 (a442712), in PR #43; step 4 (a local fit whose step matches the error, task `ar3d`) **done** 2026-11-10 (5ccbd2e), in PR #43; step 5 (let neighbouring grains compete for the boundary voxels, task `ar3e`): handed back partial, reviewed, **fixes running** since 2026-11-10 22:40 (D-20261110-N7); step 6 (the seed cost, task `ar3f`): plan c88fdd8 approved (D-20261110-N8), starts after step 5 is merged |
| AR4 | the learned proposer | waits for AR1 |
| AR5 | the stop rule | waits for AR1; first job: a per-level record of full searches on full-sample voxels (the pilot's 69), then the AR0 stop grid on it |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |

## Running jobs

Step 5 (task `ar3e`): handed back partial on 2026-11-10 at about 22:05 (report
`.claude/reports/20261110-1836-ar3e.md`), page placed by the main session (2eda7ec, 3813d69),
reviewed by Opus (`.claude/reports/20261110-2220-review-ar3e.md`: request changes, 15 fixes,
the numbers reproduce). Decisions D-20261110-N7 (all 15 fixes and additions A to G with no
rerun; the 47 MB of `tries.npz` kept; no rerun of K before step 6; step 6 is the seed cost)
and D-20261110-N8 (step 6 plan, PLAN.md c88fdd8, task `ar3f`).

**One agent is running: the Opus fixes implementer** (started 2026-11-10 at about 22:40, in
its own worktree on `feature/eval-reduction-ar3e`, head before it 3813d69). Its report will
be `.claude/reports/20261110-2233-ar3e-fixes.md`; check with `ls -la` on that path and
`git log --oneline feature/eval-reduction..feature/eval-reduction-ar3e`. When it is back:

1. `scripts/dev/verify_handback.sh feature/eval-reduction <fixes report> --branch
   feature/eval-reduction-ar3e --fixes .claude/reports/20261110-2220-review-ar3e.md
   --no-tests`; read the changed page lead and the summary.json diff statement.
2. In the worktree that holds the branch: `scripts/dev/finish_task.sh --yes --push`
   (background, about 12 minutes; keep the main checkout clean while it runs).
3. `gh pr edit 43` (add step 5), this file, one shmem fact, the review page.
4. Step 6: `scripts/dev/read_on_fetch.sh --save` then `--mark`; brief from
   `scripts/dev/make_brief.sh --plan icenine_py/docs/autoresearch/PLAN.md --task ar3f
   --section "AR3: index the layer and paint grains"`; a Sonnet implementer in its own
   worktree. Notes for it: step 6 only; the stored maps of step 5 are in the old worktree
   `.claude/worktrees/agent-a4fd6282107a76186` (untracked `boundary/cache/`, the image
   symlink); results slim, under 5 MB; "local optimization" in the page; plot first.

If the fixes agent died with the session: a fresh Opus implementer with the same brief
(`make_brief.sh ... --task ar3e --fixes .claude/reports/20261110-2220-review-ar3e.md`) and the
additions A to G of D-20261110-N7, with the note "continue from the commits on the branch".

Reviewed numbers of step 5 (scope: the simulated gradient sample of step 3 and step 2's
sample; one run each; grains of about 50 voxels, so 93 to 94% of the voxels are boundary
voxels). Keep rule `hit`, picked on A by about 2 voxels over `cost`; `conf` breaks 4,445
right voxels on A. Wrong before and after the pass (A and B pooled): D4 225 to 74, P4R4 281
to 82, E4 305 to 0, D 203 to 88, PR 170 to 78, D2 153 to 14. The pass costs 3.0 to 3.7
million evaluations per map. With the pass on both sides: D4 74 against D 88 wrong at 38.8M
against 122.1M evaluations; P4R4 82 against PR 78 at 54.6M against 94.7M; E4 0 against D2 14
at 16.6M against 36.1M; KE 0 at 8.45M. A second pass fixes 17 to 23 more per gradient map.
The skip rule (hit ratio over the best neighbour's at least 0.94) skips about 92% of the
tries and loses 0 of 360 fixed on A, 1 of 313 on B. K at its cap: 38.1M evaluations, 31.07M
in 90 grid fallbacks, no map. C++ note for the owner: in parallel mode other clients fit
accepted voxels again and the server keeps the higher confidence (the `conf` rule); not
tested.

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

### AR3 step 6 (the seed cost, task `ar3f`), merged at 73992b8 on 2026-11-11

- Page `findings_ar3_seed.md`; results `benchmarks/autoresearch/ar3/seed/`; decisions
  D-20261110-N8 (plan), D-20261111-N1 (review: 11 fixes to text, tables and one plot label, no
  rerun; three optional suggestions of the reviewer left undone).
- Replay of 189 grid-fallback voxels: on 188 the indexer had a candidate within 1 degree among
  its first 5; on 175 it was costed and not among the 2 that were refined.
- Remedy 1b (cost 20 candidates, refine the best 5): passes the seed test on 82 of 91 (pick
  set) and 88 of 98 (test set); about 4,300 evaluations.
- Whole layer, gradient sample, run plus boundary pass: `S5` (1b alone, added after seeing
  `S`) 15 fallbacks, 20,629,760 evaluations, 62 wrong; `S` (neighbour tries, then 1b) 11,
  16,140,762, 100; `K2` (library CMA keys, tries, 1b) 53, 30,343,782, 96; `D4` 76, 38,807,186,
  74; `D` 142, 122,070,127, 88.
- `S5` against `S`: 57 voxels wrong only in `S`, 19 only in `S5`. Of the 602 voxels the tries
  settle, 59 are wrong before the pass and 40 after; in `D4` 27 and 24 of the same voxels.
- `S5` by stage: seed step 5,990,794; grid fallbacks 5,098,353; neighbour 2,786,526; revisit
  3,072,127; pass 3,681,960.
- The step's rows are on the cost-against-accuracy panels (review page version 31).

### AR3 step 7 (larger grains, task `ar3g`), merged at b98a14b on 2026-11-11

- Page `findings_ar3_big.md`; results `benchmarks/autoresearch/ar3/big/`; decisions
  D-20261111-N2 (plan), D-20261111-N3 (review: 10 fixes to wording and numbers, no rerun, no
  result file changed). Suite at the merge: 1517 passed, 47 skipped, 1 deselected.
- Sample `big`: 32 grains on the 24,570-voxel grid (about 770 voxels per grain; the first
  gradient sample has 497 grains, about 49 each). Gradient levels 0, 0.2, 0.5, 1, 2 and 4
  degrees on 5, 5, 1, 9, 6 and 6 grains. The levels are my choice.
- Run plus boundary pass on `big`: `D` 65,873,391 evaluations, 26 wrong; `D4` 11,943,692, 21;
  `S5` 12,999,852, 30; `P4R4` 19,960,941, 31. Evaluations per voxel against the first
  gradient sample: `D` 2,681 and 4,968; `D4` 486 and 1,579; `S5` 529 and 840; `P4R4` 812 and
  2,221. The page does not credit grain size alone for the drop.
- The plan's predictions. 1 does not agree: `S5` has 314 reseeds on 32 grains and 35 first
  seeds; its seed step is 7,469,574 of 12,999,852 evaluations (57.5%). 2 does not agree: 529.1
  evaluations per voxel against 225.1 for the local optimizations and the pass alone. 3
  agrees on the plan's two matched pairs of levels (1 to 9 grains per level).
- `S5` against `D4` on `big`: 8.8% more evaluations; 6 voxels fixed, 15 broken, 15 wrong in
  both; the intervals overlap. `S5`'s seed step costs 1.75 to 1.77 times `D4`'s per index run.
  `D4` fell back to the grid search on 8 of 329 seed steps (2.4%) here and 76 of 1,298 (5.9%)
  on the first gradient sample; `S5` breaks even at about 3.0% and 1.75%.
- Where `S5`'s 349 seed steps are: 267 in the 6 grains at 4 degrees, 54 in the 6 grains at 2
  degrees. Its 30 wrong voxels are in 10 grains, 15 of them in two grains at 4 degrees.
- The plan's rule: the first branch fails on cost; the second has two readings. The seed step
  for the opt-in key is the owner's choice.
- The fixes, one commit each: 8b3ca8a (1), 7810869 (2), e1e4410 (3), a00bb87 (4), 1aad88c
  (5), d742b39 (6), a2df75e (7), 5d8efb7 (8), 21d9529 (9), 409961b (10), then 9396c30 and
  073ac59 (follow-ups). The fixes report (`20261111-0355-ar3g-fixes.md`) was an untracked file
  in the fixes worktree, and that worktree was removed after the merge; this list replaces it.
- The suite at the merge skipped 47 tests, 3 more than at step 6 (1520 passed, 44 skipped);
  the same 1,564 tests ran or skipped.
- The step's eight rows are a new panel of the cost-against-accuracy plots (review page
  version 33). PR #43 now holds steps 2 to 7.

### AR3 step 8 (where the seed step's evaluations go, task `ar3h`), merged at ce9b106 on 2026-11-11

- Page `findings_ar3_reseed.md`; results `benchmarks/autoresearch/ar3/reseed/` (about 2 MB);
  decisions D-20261111-N4 (plan), D-20261111-N5 (review: 10 fixes to wording and to what the
  page discloses, no rerun, no result file changed). Suite at the merge: 1517 passed, 47 skipped, 1 deselected.
- One seed step, split. In `S5` on `big`: 349 seed steps, 7,469,574 evaluations, 20 of each
  step for costing candidates, the rest in the 5 refinements (99.9%). 62.2% of the seed-step
  evaluations (55.5% on the first gradient sample) run after a refinement that, judged with
  the truth, had already passed the seed test at a right orientation; in 279 of 349 steps the
  first refinement is that one. Only the first refinement is dearer on `big`.
- Classes of the 1,471 voxels in reseed clusters of `S5` on `big`: A (the fit ended within 1
  degree and was rejected by the relative test) 1,069, median error 0.294 degrees, median hit
  ratio 0.835 against a best of 1.000; B (own grain, not reached) 22; C (start from another
  grain) 380. Why class A is rejected is not tested: the fit may stop short, or the truth
  itself may score under 0.9 times the best.
- Remedies on stored cases (pick set: true grains with an even index; test set: odd). The
  plan's `E1` (stop at the first refinement that passes 0.80) failed: test sets 6 broken and 0
  fixed on `big`, 51 and 2 on the first gradient sample. A variant with a stop test of 0.95
  (grid 0.80 to 0.98, picked on the pick set, added after the 0.80 rows of one test set had
  been read) breaks none and cuts the seed step by 58.4% and 48.4%. `E2` (keep a rejected fit
  under a second test) failed on the pick set. `E3` was not tried (class B is 22 voxels).
- Run `S8` = `S5` with the variant, plus the boundary pass. `big`: 8,340,428 evaluations, 24
  wrong (`D4` 11,943,692 and 21; `S5` 12,999,852 and 30). First gradient sample: 18,000,230
  and 63 (`S5` 20,629,760 and 62; `D4` 38,807,186 and 74). The intervals overlap.
- The plan's predictions. 1 agrees (most rejected fits in the clusters at 2 and 4 degrees are
  of class A or B). 2 agrees (the cost is in the refinements). 3 is met by the variant only,
  not a confirmed prediction (D-20261111-N5); the page makes no recommendation.
- The fixes, one commit each: 79a5f51 (1), bb27a7f (2), 14c0b75 (3), d2fca4a (4), c47c3af (5),
  b092172 (6), 2fe1e59 (7), f46fd45 (8), 64ccd78 (9), 1f5e000 (10), then 36eaf21 (fix 5 in the
  history file). The fixes report is in the main `.claude/reports/`.
- The hand-back check stopped with exit 141 on branches with more than 8 commits (a pipe into
  `head` under `pipefail`); repaired in 9fdbbe9. On the fixes report it still reports one
  field missing ("left undone"): the field is there, in the middle of a line.
- Run `S8` is on the plots of cost against accuracy (review page version 35). PR #43 now
  holds steps 2 to 8.

### AR3 step 9 (a start test, the rejected fits that end right, a third sample; task `ar3i`), merged at 0b081fb on 2026-11-11

- Page `findings_ar3_start.md`; results `benchmarks/autoresearch/ar3/start/` (about 4 MB);
  decisions D-20261111-N6 (plan), D-20261111-N7 (review: 13 fixes to wording and to what the
  page discloses, no rerun; no earlier value of `summary.json` changed). Suite at the merge:
  1517 passed, 47 skipped, 1 deselected.
- Third sample `mid` (24,570 voxels, a 120-cell Voronoi grain map, not step 7's generator: a
  deviation the page states). Evaluations with the boundary pass and wrong voxels after it:
  `D4` 40,560,874 and 48; `S5` 39,093,613 and 34; `S8` 26,443,074 and 31; `S9` 22,672,792 and
  40. The intervals overlap. D-20261111-N7 counts `S8` as confirmed on one new sample.
- Run `S9` = `S8` plus a start test (one evaluation of the start; hit ratio at least 0.90
  keeps the start and skips the fit). With the pass: `big` 5,634,362 and 29 (`S8` 8,340,428
  and 24); first gradient sample 14,631,262 and 70 (`S8` 18,000,230 and 63). Cut against `S8`:
  32.4%, 18.7%, 14.3% (`mid`). Fixed and broken against `S8`: 10 and 15, 31 and 38, 9 and 18.
  The median error rises on all three (0.0316 to 0.0550, 0.0226 to 0.0360, 0.0297 to 0.0429
  degrees). Not a candidate for a key. The test at 0.95 was not run on a whole layer.
- Class A (a fit that ended right and was rejected). The truth passes the relative test in
  all of them, by construction of the samples. One more fit from the rejected end passes in
  807 of 1,069 (`big`) and 451 of 589 (first gradient sample).
- Remedy `F2` (one more fit before rejecting, at 0.2 degrees) failed: net cut of the
  seed-step evaluations of the reseed clusters 34.1% on `big`, -11.8% on the first gradient
  sample, 4.8% pooled; the plan asks a third. 137 of 672 clusters vanish.
- First gradient sample at 1 degree: 53 against 56 wrong of 5,539; 23 fixed, 26 broken.
- The plan's predictions: 1 agrees (by construction), 2 does not agree, 3 agrees (`S8` on
  `mid`), 4 agrees on `big`.
- The fixes, one commit each: afec122 (1), d78468a (2), 8ff31cb (3), cf8744e (4), 3cc10cb (5),
  75b4173 (6), 556c105 (7), f312e84 (8), 3e3c646 (9), 25f1360 (10), 94bc34f (11), 47cfcda
  (12), 8f71633 (13); then a152ac5, 33d958e, 47769e0, 21d88ac. The fixes report is in the
  main `.claude/reports/`.
- The hand-back check on the fixes report showed one FAIL: report fields missing. The cause
  is the main session's brief, which asked for other field names; the content is there.
- The findings page was placed by the main session from the implementer's text (its Write
  was refused; the standing exception the owner knows).
- The runs on `mid` and run `S9` are on the plots of cost against accuracy (review page
  version 36). PR #43 now holds steps 2 to 9.

### AR3 step 10 (opt-in keys for the indexed seed step and the boundary pass; task `ar3j`), merged at 17790e5 on 2026-11-11

- Page `findings_ar3_keys.md`; results `benchmarks/autoresearch/ar3/keys/` (about 1 MB);
  decisions D-20261111-N8 (plan), D-20261111-N9 (review: 8 fixes, a test of the reading about
  the pass inside the run, the repo page keeps the word "fit"). Suite at the merge:
  1576 passed, 47 skipped, 1 deselected.
- The first AR3 step since step 2 that changes the library: new `icenine/boundary_pass.py`;
  changes in `reconstructor.py` (`_local_fit`, `_index_seed`, `local_optimization`,
  `reconstruct_sample`, new `check_bfs_keys`), `orientation_search.py`, `config_file.py`,
  `ray_indexing.py`. Tests: `tests/test_bfs_ar3_keys.py`, 59 tests (57 written before the code).
- Keys, all off by default: `BFSNeighborCMAFirst 0|1`; `SeedIndexCandidates n_cost n_refine`
  (default `5 2`; a third token raises `ValueError`); `SeedIndexStopHitRatio t` (default 0);
  `BoundaryPass 0|1`. The stored runs as keys: `D4L` = `BFSNeighborCMAFirst 1`, `CMASigma0
  0.05`, `CMANeighborMaxEvals 60`, `BoundaryPass 1`; `S5L` adds `SeedIndexCandidates 20 5`;
  `S8L` adds `SeedIndexStopHitRatio 0.95`.
- Sample `big`, clean, Q-max 8, contended, one run each. Before the pass `D4L`, `S5L`, `S8L`
  equal the stored `D4`, `S5`, `S8` exactly: seeds / reseeds / grid fallbacks 38/291/8,
  35/314/0, 37/292/0; evaluations 10,811,715, 11,867,509, 7,201,863; wrong 68, 71, 67; every
  counter, the 24,570 records, the `.mic` bytes. `S8L`: 322 early stops in 329 seed steps.
- The pass on the stored (rounded) maps is exact: tries / evaluations 18,557 / 1,131,977,
  18,563 / 1,132,343, 18,665 / 1,138,565. In the run: `S5L` identical (30 wrong); `D4L` the
  same counts, 22 wrong against 21 (voxel 19614); `S8L` 18,664 tries / 1,138,504 evaluations,
  24 wrong in both, 5 fixed and 5 broken. Wilson 95% intervals after the pass: `D4L` 22
  (0.059% to 0.136%), `S5L` 30 (0.086% to 0.174%), `S8L` 24 (0.066% to 0.145%).
- The reading test (N9, fixes round): `S8L` and `D4L` rerun (886 s, 1,133 s, contended). For
  both: the rerun's in-run pass equals step 10's; the pass on the `.mic`-rounded map in the
  same process gives the stored tries exactly; a fresh process on the unrounded map gives the
  in-run pass try for try. So the rounding accounts for the difference, the process does not.
  `S5L` was not rerun. Rounded against unrounded in one process: every boundary voxel's
  float64 ends differ (11,087 and 11,054 voxels; largest rounding 8.8e-7); not a prediction.
- Evaluations not counted in the totals: the strict score of each boundary voxel's old end
  (11,087 / 11,049 / 11,054) and `S8L`'s 465 hit-ratio evaluations of the early stop.
- Off path: the 2,000-voxel pilot with every key absent gives byte-identical files on the
  parent commit (4ff3ce6) and on the branch. Not rerun after the fixes (they changed
  annotations, docstrings, an error message and the `ValueError` on a third token).
- No timing: the preflight refused twice (another program at about 100% of a CPU). The
  single-worker wall pair stays for "Machine-time runs".
- Commits: e289607 (tests first), 3fbfe2f (keys), b0716ed, 5b04c9f; fixes 07037d4, 7b992ab,
  6989d53, 159dcc2, 5c0487f, cd3a2b5. Reports in the main `.claude/reports/`:
  `20261111-0834-ar3j.md`, `20261111-0941-review-ar3j.md`, `20261111-1000-ar3j-fixes.md`.
- The findings page of the first round was placed by the main session from the implementer's
  text (its Write was refused; the standing exception the owner knows). In the fixes round the
  implementer's own edit of the page worked.
- Tooling trap found by this step: in a worktree `uv run pytest` could import the main
  checkout's library. Fixed in fb06e95 (the dev scripts use `uv run --extra dev pytest`;
  briefs tell implementers to `uv sync --extra dev` and check `uv run which pytest`).
- The heads-up to Ilium (protocol rule 3(a)) was posted late, after the task started and
  before the merge. Ilium's branch did not move (347a8a4) and its mailbox had no objection.
- The plots of cost against accuracy get no new point (review page version 37 says so).
  PR #43 now holds steps 2 to 10.

### AR5 step 1 (stop rules for the boundary pass, the grid fallbacks and the finished map; task `ar5a`), merged at e51e883 on 2026-11-11

- Page `findings_ar5_stop.md`; results `benchmarks/autoresearch/ar5/` (about 4.5 MB); scripts
  `scripts/autoresearch/ar5/stop_*.py`; tests `tests/test_ar5_stop_rules.py` (20). Decisions
  D-20261111-N10 (plan), N11 (review). Suite at the merge:
  1596 passed, 47 skipped, 1 deselected.
- Scripts only. A negative result: no rule passed its pass rule; run `T1` was not run; no new
  point on the plots. All on stored `S8` runs, clean images, one run each; picks on `big` and
  the first gradient sample, test on `mid`.
- Part 1, try test for the boundary pass. Pick `P2 m=1` (keep the try with the best start hit
  ratio): pick samples 4,821,806 to 2,164,026 evaluations (55.1% cut), net gain 181 of 190
  (needs 180.5). On `mid`: 62.7% cut, 9.4% of the gain lost (48 of 53; 0 fixed, 5 broken):
  fails. `P2 m=2` loses no gain on any sample, cuts 16.6%, 30.2%, 38.9% (under the 50% rule).
  `P1` and `P3` cut 95.7% to 98.2% and lose 17.9% to 60.0%.
- Part 2, the 31 grid-fallback cases (16 + 15). A near-truth candidate was costed and refined
  in all 31 (rank 1 in 27). `R1` costed every candidate (lists of 32 to 62) and repeated the
  seed step's refinements in 20 of 31; it passes the seed test in 0 (23 right ends at hit
  ratio 0.6561 to 0.7998, 8 wrong at 0.6964 to 0.7855). `R3` passes 17, all right. Ladder
  11,364,144 against grid 10,606,351 evaluations (107.1%). What-if with level 0 reused:
  8,536,861 (80.5%), above two thirds.
- Part 3, the finished map. Hit ratio under 0.95 flags 87 of 87 and 31 of 31 wrong voxels and
  8.8% and 15.1% of all voxels. `V1` on `mid`: 13 of 31 fixed, 0 broken, 27,882,533
  evaluations (105.4% of the run; 27,717,299 of them on right voxels). Vote count: weak; best
  cut moves from 1,975 to 3,035 between samples; cut picked on the pick samples (2,919)
  catches 12 of 31 wrong and 281 of 3,678 right on `mid` (a what-if).
- Review (fix 1): a first computation of `mid`'s start scores ran before the Part 1 pick was
  committed and was redone. The picks are mechanical; no plan setting passes on `mid` with
  hindsight. Not adopted (N11): `V1`, the vote count as a stop rule.
- Open for the owner: `P2 m=2` as an opt-in key.
- Commits: 53f8951 (tests first) to 29e9b56; fixes 35527a7 (parent merged), dfefe3a, 242f30c,
  0f5aa0b, 169abb5. Reports in the main `.claude/reports/`: `20261111-1008-ar5a.md`,
  `20261111-1215-review-ar5a.md`, `20261111-1134-ar5a-fixes.md`. The findings page of the
  first round was placed by the main session from the implementer's text (its Write was
  refused; the standing exception). Worktree with the cache: `agent-a68538c58665829a1`.
- What it led to: AR4 recast from proposer to verifier (N13), AR5 step 2 (N12).

### AR5 step 2 (two long local optimizations; task `ar5b`), merged at c688191 on 2026-11-11

- Page `findings_ar5_fits.md`; results `benchmarks/autoresearch/ar5/fit/`; scripts
  `scripts/autoresearch/ar5/fit_*.py`; tests `tests/test_ar5_fit_rules.py` (13). Decisions
  D-20261111-N12 (plan), N14 (review). Suite at the merge: 1609 passed, 47 skipped, 1
  deselected.
- Scripts only, clean images, one run each; picks on `big` and the first gradient sample, test
  on `mid`. The pick file was committed at 12:00:26, the first `mid` replay started 12:00:47.
- Base check: the script layer with both changes off gives stored `S8` on `big` exactly
  (7,201,863 evaluations, 67 wrong before the pass).
- Part 1, replay of 3,219 seed steps. Arm `M` (stored MC) matched the log on every step. No
  CMA arm met the 1% rule: `C200` loses 60 seed-test passes net on the 1,608 steps of the pick
  samples (64 with `M` only, 4 with the arm only; limit 16). Pick stays `M`. End errors on
  `mid`, median / 90th percentile / within 0.2 degree: `M` 0.020 / 0.076 / 97.8%, `C60` 0.160
  / 23.341 / 56.0%, `C200` 0.013 / 0.317 / 85.5%, `C0` 0.423 / 41.431 / 16.1%.
- MC fallback after CMA ended lower than the CMA end in 3 of 3,197, 3 of 3,903 and 45 of 4,568
  runs (`big`, first gradient sample, `mid`; the last is 0.985%).
- Part 2, run `S10` (`S8` with no MC fallback), with the pass: 5,898,258 / 15,167,998 /
  22,668,704 evaluations against 8,340,428 / 18,000,230 / 26,443,074 (29.3 / 15.7 / 14.3%
  fewer). Wrong after the pass 22 / 67 / 35 against 24 / 63 / 31; fixed and broken 11/9,
  29/33, 11/15. Passes the 20% rule on `big` only. `S10a` and `S10b` were not run (with pick
  `M` they equal `S8` and `S10`). Predictions 1, 3, 5 agree; 2 and 4 do not.
- Review's what-if from the stored replay (arithmetic, not a run): `C200` first, MC only when
  the CMA end hit ratio is below 0.95: passes 1,596 against 1,592 (pick samples), 1,600
  against 1,596 (`mid`); seed-step evaluations 2,501,991 against 5,865,469 and 5,643,374
  against 12,007,352. This became AR5 step 3 (N14). No key for `S10` alone (N14).
- The findings page was placed by the main session from the implementer's text (commit
  591ad30). Six `S10` rows in `docs/autoresearch/pareto/pareto_rows.json`.
- Commits: 71bbdee (tests first) to aeda9e4, 591ad30; fixes dbd4a9f (parent merged), 573a599,
  a186b34. Reports in `.claude/reports/`: `20261111-1134-ar5b.md`,
  `20261111-1310-review-ar5b.md`, `20261111-1317-ar5b-fixes.md`. Cache (untracked):
  worktree `agent-acc78ac2a483ed86f`, `ar5/fit/cache/`.

### AR4 step 1 (a learned verifier, offline; task `ar4a`), merged at 1ffee63 on 2026-11-11

- Page `findings_ar4_verifier.md`; results `benchmarks/autoresearch/ar4/`; scripts
  `scripts/autoresearch/ar4/verify_*.py`; tests `tests/test_ar4_verifier.py` (26). Decisions
  D-20261111-N13 (plan; AR4 recast from a learned proposer to a learned verifier), N15
  (review: request changes, 10 fixes, all applied; no new run). Suite at the merge: 1635
  passed, 47 skipped, 1 deselected.
- Scripts only, clean images. Models `B0` (hit ratio), `B1` (logistic regression), `N1` (small
  set network). Trained on `big` and the first gradient sample; thresholds committed at
  12:43:51 (669b9a3), first read of `mid` at 12:44:18.
- Use 1, flag of unfinished voxels on `mid` (31 wrong of 24,570, all at boundary voxels): `B0`
  flags 31 of 31 at 1,760 voxels, `N1` 23 of 31 at 190, `B1` 31 of 31 at 13,653. No model
  passes. `N1` scores 8 wrong voxels at 0.83 to 0.99 (it misses 4 of 87 in-sample).
- Use 2, accept a seed that failed the seed test (15 fallback cases on `mid`): `N1` accepts 7
  of 10 right ends and 0 of 5 wrong. Arithmetic on stored numbers, not a run: as an added rule
  it skips 2,325,301 of 5,002,750 grid evaluations on `mid` (3,744,760 of 5,449,501 on the
  first gradient sample, in-sample). In place of the seed test it refuses 20 right and 9 wrong
  passing seeds on `mid` (about 6,670,333 and 3,001,650 extra evaluations at the mean grid
  search of 333,517).
- Grid-only evaluations of the fallback cases equal the library counter of `S8`; AR5 step 1's
  figures were larger by the seed steps of the reseed cases (85,014 and 69,086).
- The share of peaks with any lit pixel does not separate right from wrong at boundary voxels
  (AUC 0.496 against 0.974 for the hit ratio; every peak lit for 77.1% of right and 77.4% of
  wrong). Predictions 1 and 3 agree; 2 and 4 do not. The vote-count input was skipped.
- Per-peak sums differ from the fused evaluation on 2.5% to 2.8% of pairs, the same for right
  and wrong; consistent with D-20261027-I2, not traced pair by pair.
- Next (N15): use 2 only as an added rule, tested on the unread sample `held` (AR3 step 11).
- Worktrees: `agent-ab349f0b6193c2ca1` (first round; `ar4/cache`), `agent-a5a51f7021e6cc8c1`
  (fixes round; a link to that cache). Reports `.claude/reports/20261111-1208-ar4a.md`,
  `20261111-1335-review-ar4a.md`, `20261111-1351-ar4a-fixes.md`.

### AR5 step 3 (a mixed seed refinement, run `S11`; task `ar5c`), merged at ab7e9fc on 2026-11-11

- Page `findings_ar5_mix.md`; results `benchmarks/autoresearch/ar5/mix/`; scripts
  `scripts/autoresearch/ar5/mix_*.py`; tests `tests/test_ar5_mix_rules.py` (14). Decisions
  D-20261111-N14 (plan), N16 (review: approve, 8 wording and traceability fixes, all applied;
  no new run). Suite at the merge: 1649 passed, 47 skipped, 1 deselected.
- Scripts only, clean images, one run per sample. The pick was committed at 13:57:26; the
  20-step check started 13:58:20 (it matched Part 1 on all 20 steps) and the runs at 14:05:51.
- Part 1, from the stored replay: pick `X200` (CMA with 200 evaluations first; MC from the same
  candidate only when the CMA end hit ratio is under 0.95; the lower cost of the two). Pick
  samples: 57.3% fewer seed-step evaluations than `M`, 1,596 seed-test passes against 1,592,
  1,589 ends within 1 degree against 1,583. `X60` meets both rules too but saves 25.2%. On
  `mid` 5,643,575 seed-step evaluations; the plan's what-if 5,643,374 is a variant that stops
  when either end reaches 0.95 and differs at replay step 442 only.
- Part 2, run `S11` (`S8` + `X200` + no MC fallback), with the pass: 4,015,342 / 12,307,994 /
  16,082,696 evaluations on `big`, the first gradient sample and `mid` (51.9 / 31.6 / 39.2%
  fewer than `S8`). Wrong 18 / 62 / 36 against 24 / 63 / 31. Broken minus fixed against `S8`:
  -6, -1, +5 (`mid`: 15 broken, 10 fixed; exactly at the limit). The pass rule is met on all
  three.
- Grid fallbacks 0 / 11 / 12 against 0 / 16 / 15. They are 30.0% and 9.6% of the saving against
  `S8` on the first gradient sample and `mid` (47.8% and 15.2% against `S10`); one fallback
  costs 333,517 to 340,594 evaluations. Why they drop was not tested.
- Predictions 1, 2 and 4 agree; 3 (fallbacks within 2 of `S8`) does not.
- The rule was found in a replay that included `mid`: no sample is held out for it. N16: no
  key before the test on `held` (`S10`, `S11`, `S11` plus verifier `N1`, by stage).
- Worktrees: `agent-a5a3601b746680640` (first round; the `S11` caches), `agent-a566790c79fbe9096`
  (fixes round; links). Reports `.claude/reports/20261111-1351-ar5c.md`,
  `20261111-1503-review-ar5c.md`, `20261111-1520-ar5c-fixes.md`.

### AR3 step 11 (a held-out sample with its `S8` baseline; task `ar3k`), merged at 504fadb on 2026-11-11

Scripts and results only. Suite after the merge: 1656 passed, 47 skipped, 1 deselected.
Decisions D-20261111-N15 (plan), N17 (second run with a cap of 480 minutes), N19 (review:
request changes, 7 document and script fixes, no rerun). Page:
`docs/autoresearch/findings_ar3_held.md`. Reports in `.claude/reports/`: `20261111-1351-ar3k.md`,
`20261111-1900-review-ar3k.md`, `20261111-1930-ar3k-fixes.md`.

- The sample `held`: 250 Voronoi grains (9 to 271 voxels, median 89.5), generator seed 1, the
  24,570 voxel centres and six gradient levels of `mid`; 21,131 boundary voxels (86.0%);
  10,897,893 lit pixels, 179,718 spots. The builder rebuilds `mid` and `held` exactly.
- The first `S8` run reached its cap of 120 minutes with no map written (48,312,402
  evaluations before the pass, 61 grid fallbacks, 3,237 seed steps; lower bounds).
- The second run (cap 480 minutes) finished: 54,733,887 evaluations before the pass,
  59,923,523 with it. By stage: seed step 21,257,244, grid fallback 25,915,118 (75 fallbacks),
  neighbour 3,004,901, revisit 4,556,624, pass 5,189,636. 3,567 seed steps, 3,262 reseeds.
- Wrong voxels: 95 before the pass (Wilson 0.32 to 0.47%), 70 after (0.23 to 0.36%). Median
  error 0.0288 and 0.0286 degree. The pass fixed 30 and broke 5. At gradient level 4 degrees it
  fixed 0 and broke 4 (49 wrong before, 53 after); the five lower levels went from 46 to 17.
- Predictions: 1 (20 to 40 million evaluations) missed, 2 (30 to 120 wrong) met, 3 (5 to 40
  grid fallbacks) missed. Why `held` has so many fallbacks, and why the pass harms level 4,
  was not tested.
- Base check on `big` exact (7,201,863 evaluations, 67 wrong). No other run on `held` in this
  step. Four heavy processes shared the machine in the second run (the plan allowed two); all
  wall times are contended.
- Data that stays untracked: worktree `agent-af8961874e131bda2` (the sample, the stack
  `held_clean.npy`, `ar3/held/cache/S8_held/` and `S8_held_cap120/`).

### AR5 step 5 (two opt-in library keys for run `S11`; task `ar5e`), merged at a494bcd on 2026-11-11

Library change (four files), tests, scripts, results. Suite after the merge:
1724 passed, 47 skipped, 1 deselected.
Decisions D-20261111-N21 (plan), N22 (review: request changes, six fixes to the page and the
tests, no rerun, no library change). Page: `docs/autoresearch/findings_ar5_keys.md`. Reports in
`.claude/reports/`: `20261111-2030-ar5e.md`, `20261111-2115-review-ar5e.md`,
`20261111-2200-ar5e-fixes.md`. Heads-up to Ilium in the mailbox of `docs/ROADMAP.md` (dc0e6a5).

- The keys, both off by default: `SeedIndexCMAFirst n` (mixed seed refinement: a CMA fit first
  for each indexed seed candidate; MC after it only when the CMA end misses the stop ratio;
  `S11` uses 200) and `BFSNeighborCMANoFallback 0|1` (a neighbour fit returns its CMA end).
- Parity before the pass, library run against stored script run, all exact (evaluations,
  seeds and reseeds, grid fallbacks, voxel records, map bytes): check 1 keys off on `big`
  7,201,863 (67 wrong); check 2 `S10` on `big` 4,767,135; check 3 `S11` on `big` 2,879,705;
  check 4 `S11` on the first gradient sample 8,622,191 (348 seeds, 986 reseeds, 11 grid
  fallbacks, 207 wrong). The reviewer recomputed checks 1 to 4 from the raw files.
- The boundary pass (check 5) is exact from the stored maps (totals 4,015,342 on `big`,
  12,307,994 on the gradient sample). Inside the `S11` library run it has the stored tries and
  evaluations but not the stored contents: costs differ in 4,372 of 18,617 tries on `big` and
  in 15,813 of 60,423 on the gradient sample; changed voxels 986 against 983 and 688 against
  681; wrong after the pass 19 against 18 and 62 against 62. The cause (the run reads the
  unrounded map) was shown only for the `S8` library run in AR3 step 10; untested for the others.
- The hit ratio of a CMA end is one more overlap evaluation per refined candidate; it is not
  counted (as in the script runs).
- Predictions: 1 to 3 met; 4 not met (62 new tests at the first hand-back, 68 after the fixes,
  against 25 to 50). No existing test changed.
- Not tested: timing; the keys on `mid` or `held`; no default changed; no setting recommended.
- Data that stays untracked: worktree `agent-afc3f195cb4fffb30` (`ar5/keys/cache/`: the library
  runs `S8L_big`, `S10L_big`, `S11L_big`, `S11L_grad` and the pass checks).

### AR5 step 7 (opt-in key `BoundaryPassSkipHitRatio`; task `ar5g`), merged at 50c3b9c on 2026-11-12

Library change (four files: `boundary_pass.py`, `config_file.py`, `orientation_search.py`,
`reconstructor.py`), tests, scripts, results. Suite after the merge: 1759 passed, 47 skipped,
1 deselected. Decisions D-20261112-N4 (plan), N5 (phase D: aligned draws, decided after the
first live result), N6 (review: request changes, twelve fixes, 1 to 11 applied; fix 12, placing
the page, waits for the owner). Reports: `20261112-1302-ar5g.md` (worktree
`agent-a331e7882919e74aa`), `20261112-1345-review-ar5g.md` (main checkout),
`20261112-1355-ar5g-fixes.md` and the page text `20261112-ar5g-page-v2.md` (worktree
`agent-a0db870adb7758a69`). The page `findings_ar5_gate.md` is not in the tree yet.

- The key, off by default: with `BoundaryPass 1`, a boundary voxel whose own hit ratio before
  the pass is at least `h` gets no try. A skipped try still consumes its random draw (aligned
  draws), so the kept tries are the same as in the full pass.
- Replay at 0.95 on stored passes: 0 fixed, 0 broken on `held` and on five disclosed passes;
  the pass is cut by 69% to 84% (`held` 72.3%, pass 1,436,977 evaluations).
- Variant 1 (the plan's pre-set live test, skipped tries do not draw): same evaluations, but
  fixed / broken against the full pass 4/7 `mid`, 15/12 `grad`, 2/6 `big`, 4/10 `held` (75
  wrong against 69). It fails the plan's rule on `held`. Runs only at commit 7b4fd0f.
- Variant 2 (aligned draws): every kept try identical to the stored full pass (15,592 `mid`,
  10,074 `grad`, 4,491 `big`, 23,557 `held`), 0 fixed, 0 broken, wrong 35, 59, 18, 69. `held`
  run total 33,753,877: 10.0% below the full pass, 43.7% below `S8`. The rule passes on all
  four samples, but `held` is not a clean test: the variant was decided after `held` was read.
- Lesson: one shifted draw stream alone moved 4 to 15 voxels each way per sample (one re-draw
  per sample). A live test of a change that only removes work uses aligned draws.

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
