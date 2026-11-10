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
4. Do one step. Then update this file, commit, push, and save one shmem fact
   (`shmem add --source llm "IceNine autoresearch <step>: <result>"`).

Rules that hold for every step: Python only through `uv run`; one fresh implementer agent per
task in its own worktree (brief from `scripts/dev/make_brief.sh --plan docs/autoresearch/PLAN.md
--task <task> --section "<heading>"`), `scripts/dev/verify_handback.sh`, then `/review-task`, fixes
applied, merge with `scripts/dev/finish_task.sh`; at most 6 heavy processes; no merge to `develop`
without the owner's go; no change to `.claude/` files; no p-values; no extra seeds; findings in
plain English with the plot before the table.

## The wake-up job

A session cron job (id 6f6deed3, created 2026-11-09, every 3 hours at :23) sends the wake-up
prompt while the session is open and idle. It expires on 2026-11-16: create it again on or after
2026-11-15 and write the new date here. It dies with the session. After a closed session the owner
types `/resume` or "continue the autoresearch"; this file is enough to go on.

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
  layer. No implementer has started yet.

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
6. **AR3 step 1 (task `ar3a`) is handed back and under review** (2026-11-10, about 04:30).
   - Branch `feature/eval-reduction-ar3a`, in the worktree
     `.claude/worktrees/agent-a3fcfcc95d9db6295`: 6af4cd2 (module `icenine/grain_paint.py` and
     tests), 83865d5 (runs, summary, figure), f7a1dcf (README, MIGRATION_HISTORY), d66b6bf (the
     findings page `docs/autoresearch/findings_ar3_paint.md`, written by the main session from
     the implementer's text, tables synced from `tables.md`). Not merged yet.
   - Report: `.claude/reports/20261110-0353-ar3a.md`. `verify_handback.sh`: 0 FAIL, 2 FLAG, both
     expected (untracked run files and the image symlink in the worktree; the new library module
     the plan asks for).
   - The implementer's numbers, not yet reviewed: paint 156,766 evaluations against 2,819,906
     for key on; 0 wrong of 2,000 against 22; the gate for step 2 is met. Do not quote them as
     results until the review is done.
   - A code-reviewer (model opus) is running; its report goes to
     `.claude/reports/20261110-0430-review-ar3a.md`. When it returns: if there are required
     fixes, a fresh Opus implementer applies them (`make_brief.sh --plan ... --task ar3a
     --section "AR3: index the layer and paint grains" --fixes <review report>`; it must work on
     the same branch, so no agent may still be running there), then `verify_handback.sh ...
     --fixes`, then `finish_task.sh --yes --push` from the task branch (run it in the worktree
     that holds the branch).
   - If the session is new and no review report exists, start the reviewer again with the
     `/review-task` brief.
   - The gate for step 2 is fixed in the plan: paint breaks at most 20 voxels against the key-on
     run and needs at most half of its evaluations.
   After step 1 is merged: write step 2's details into the plan (the key `PaintGrains`, the whole
   layer), log the approval, and go on the same way. Open the new PR to `develop` then (point 7).
7. Item 7: **PR #41 is merged** into `develop` (357a2b3, 2026-11-09, at branch head 16c72db), on
   the owner's word in the session (D-20261109-N11). The branch `feature/eval-reduction` is kept.
   There is no open PR now. Open a new PR to `develop` when AR3 is merged on the branch (title
   and text like PR #41's, listing what landed since 357a2b3: the item 22 hook change fd1721f,
   the tightened pre-commit check f76f776 and AR3). Do not merge that PR without the owner's word in the session.
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
| AR3 | index the layer and paint grains | plan written and approved 2026-11-10 (D-20261110-N1); **next**: step 1, the implementer |
| AR4 | the learned proposer | waits for AR1 |
| AR5 | the stop rule | waits for AR1; first job: a per-level record of full searches on full-sample voxels (the pilot's 69), then the AR0 stop grid on it |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |

## Running jobs

One agent is running (2026-11-10, started about 04:30): the code-reviewer of AR3 step 1 (task
`ar3a`). It only reads. Check: does `.claude/reports/20261110-0430-review-ar3a.md` exist? The
session is told when it finishes. The implementer has finished (see Next action 6).

## Results so far

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

## For the owner

- Still yours: item 12 (apply and commit, below), Ilium's plan approvals, and the note on what
  the pre-commit check cannot do (below).
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
