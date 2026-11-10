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

## Status (2026-11-09, evening)

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
5. AR2 step 2 (the key `SeedByIndexing`, main arm: finisher on the best 2): the implementer
   handed back on 2026-11-09 at 18:25. Branch `feature/eval-reduction-ar2b`, worktree
   `.claude/worktrees/agent-aa30ebdcbce5ec948`, 6 commits (head 7789887). Its write of the
   findings page was refused, so the main session placed `docs/autoresearch/findings_ar2_key.md`
   from its text, synced the tables and ran the number audit (commit 7789887).
   `verify_handback.sh` gave 0 FAIL. Its report is in the worktree:
   `.claude/reports/20261109-1700-ar2b.md` (copy it to the main checkout's `.claude/reports/`).
   The review is done: REQUEST CHANGES, 14 required fixes, all on the page, the summary and
   tables, or the figure; the code is correct and every number recomputed
   (`.claude/reports/20261109-review-ar2b.md` in the worktree). Its two questions are decided in
   D-20261109-N9 (938a342). **Fixes, state at 19:20:** the Opus fix implementer
   committed 1cee7261d (`key_summary.py`, summary, tables, figure: fixes 1, 4, 5, 6, 7, 9, 10,
   12). Its write of the revised page was refused, so the page is still as at 7789887 and fixes
   2, 3, 8, 11, 13 and the page side of the others are not in the repo yet. The agent was
   resumed to return the page text; the main session places it over
   `docs/autoresearch/findings_ar2_key.md`, runs `sync_doc_tables.py` (then `--check`) and
   `audit_numbers.py`, fills the audit count in the provenance sentence, and commits. If the
   agent died with the session: start a fresh Opus implementer with the same brief plus "return
   the page text, do not write the page". Its finding on fix 1: the two key-on runs differ in
   orientation on 2 of 2,000 voxels, not 637 (the review used the angle function with a floor).
   Its new case (seeds in both runs, only the indexed result right: voxels 5876, 24334; 16476 in
   the FusedCost pair) is decided in D-20261109-N10 (44c24f2): counts for the key. Report:
   `.claude/reports/20261109-1838-ar2b-fixes.md` in the worktree. After the page is placed: `scripts/dev/verify_handback.sh feature/eval-reduction <report>
   --branch feature/eval-reduction-ar2b --agent N` from the worktree, check that each fix number
   is applied, then merge with `../../../scripts/dev/finish_task.sh --yes --push` from the
   worktree. If the fix agent died with the session, regenerate the brief
   (`scripts/dev/make_brief.sh --plan icenine_py/docs/autoresearch/PLAN.md --task ar2b --section
   "AR2: indexing as the BFS seed search" --fixes <review report>`), add "work in the existing
   worktree", "N9 binds fixes 3, 4, 8, 9" and "fix 14 is the main session's", and start a fresh
   Opus implementer. At the merge: fix 14 (the MIGRATION_HISTORY entry), a heads-up to Ilium in the ROADMAP mailbox (the branch changes the
   shared files `icenine/orientation_search.py` and `icenine/reconstructor.py`, opt-in key only),
   ROADMAP row AR, MIGRATION_HISTORY, a shmem fact, and the owner's review page
   (https://claude.ai/artifact/A63khAHf7nEx6NBZioSnGW; source in the session scratchpad,
   republish from a read of the page if the scratchpad is gone). Unreviewed numbers from the
   hand-back are in "Results so far" below. Then AR3.
6. Item 7: PR #41 is open and has the test count. Do not merge it; the owner does.

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
| item 7 | open the PR to `develop` (merge waits for the owner) | **open**: PR #41 (2026-11-09). Not merged. Test count added to its body 2026-11-09 |
| item 14 + 3 | both keys on E0 and the pilot, timing A/B | not started; needs a quiet machine |
| item 18 | does C++'s clang fuse the clip's multiply-adds | **done** 2026-11-09: yes (below) |
| C++ time | measured C++ time on the pilot region (owner's question on finding 6) | not started |
| AR2 | indexing as the BFS seed search | step 1 (holdout) **done** 2026-11-09 (86dd905); step 2 (the key) handed back 2026-11-09, in review, branch `feature/eval-reduction-ar2b` |
| AR3 | index the layer and paint grains | after AR2's first result |
| AR4 | the learned proposer | waits for AR1 |
| AR5 | the stop rule | waits for AR1; first job: a per-level record of full searches on full-sample voxels (the pilot's 69), then the AR0 stop grid on it |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |

## Running jobs

Agents of the session of 2026-11-09 (they die with the session; `git worktree list` shows their
worktrees, and their branches keep the commits): the fix implementer of AR2 step 2 (branch
`feature/eval-reduction-ar2b`, started 2026-11-09 about 19:00). No compute job is running.

## Results so far

AR2 step 2, **as handed back, not yet reviewed** (2026-11-09; page `findings_ar2_key.md` on the
task branch). Pilot region, clean images, Q-max 8, every time contended (the preflight refused).
- The stored key-off baseline reproduced: 37 wrong, 25,547,935 evaluations, 69 full searches.
- Key on: 22 wrong of 2,000 against 37. Paired: 22 fixed, 7 broken, 15 wrong in both. All 7
  breaks are edge-sharing boundary voxels that take a touching grain's orientation, carried by
  local fits from an indexed seed that is right for its own voxel. 0 of 64 indexed seeds and
  reseeds are wrong. 0 fell back, so the fallback path ran on no real voxel.
- The key-off wrong count is not stable: the same key-off config with `FusedCost` on gives 24
  wrong, not 37 (20 fixed, 7 broken). Against that run the key gives 14 fixed, 12 broken. So the
  reading is "no worse than the grid search here, not shown better". The cause is not known.
- Evaluations: 2,819,906 against 25,547,935 (about 9.1 times fewer); 64 full searches at 1,909
  each against 69 at 331,636. The rest (2.70 M) is the local fits, unchanged.
- Wall time, `FusedCost` on, one worker, contended, one run per arm: 213 s key on, 1633 s key off.
  Pairing 17.6 s once, indexing 0.186 s per voxel, outside the evaluation count.
- Voxel 17361: both finisher stages ran and accepted nothing. One `RayGaussNewton` step reaches
  0.128 (exact truth 0.0949). The strict cost is 0.156 at 0.02 degrees from the truth and the
  finisher's first step is 0.132 degrees: consistent with a minimum narrower than the first step,
  not tested with a smaller step.
- Open questions it raises: the BFS pilot pairs rays with the hull of the 2,000-voxel region, not
  of the whole sample as AR1 did; key-off bit-identity was shown with fakes and the golden grid,
  not on real images.

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

- Still yours: the merge to `develop` (item 7), items 12, 22 and 23, and Ilium's plan approvals.
- `shmem ask` works again (owner updated the key, 2026-11-09). It still prints one warning: the
  local embedder (ollama on port 11434) is not running, so new facts are not embedded and a
  search uses what is already embedded plus keywords. Answers were correct in a test.
