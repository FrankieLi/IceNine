# Autoresearch stage: saved state

Read this file first after any restart. It says where the work stands and what to do next. The plan
is `docs/autoresearch/PLAN.md`. The mandate is D-20261109-N1 in `docs/DECISIONS.md`.

## How to restart

1. In the repository root, check `git status` and the branch (`feature/eval-reduction`
   or one of its task branches).
2. Read "Next action" below. For more context: `shmem query "IceNine autoresearch"` (`shmem ask`
   was broken on 2026-11-09; try it, fall back to `query`).
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
- AR1 is built and measured (branch `feature/eval-reduction-ar1`, 4 commits up to 3182d40). The
  hand-back check passes (0 FAIL). The review is running. Not merged yet. Numbers below are
  not reviewed yet.

## Next action

1. AR0: done. Nothing left.
2. AR1: the review report will be `.claude/reports/20261109-review-ar1.md`. If it is not there, run
   `/review-task` again from the worktree `.claude/worktrees/agent-af5119989737208f0`. Apply its
   required fixes with a fresh Opus implementer in that worktree. Merge from the worktree with
   `../../../scripts/dev/finish_task.sh --yes --push` (with tests, in the background). Do not commit
   the untracked symlink and cache directories of that worktree. Then write the reviewed result
   here, in MIGRATION_HISTORY, in the README (new module `icenine/ray_indexing.py`), and in the
   ROADMAP (row AR and a note to Ilium about the new module).
3. If an agent died with the session (a task branch exists but no report): look at
   `git log feature/eval-reduction..feature/eval-reduction-ar1` in its worktree, and start a fresh
   implementer with the same brief (regenerate it with `make_brief.sh`) and the note "continue from
   the commits on the branch".
4. After AR1 is merged: the gate is met (69 of 69, if the review confirms it). Start AR2 (opt-in
   key `SeedByIndexing`, pilot run paired against the key-off rerun) and AR3 (index the layer, paint
   grains). Open choice for AR2: the finisher on the best 1 candidate or on the best 2. Best 2 had
   0 of 69 wrong but was chosen after the fact; AR2's pilot run is the test on new voxels. AR4 (the
   learned proposer) aims at what voting misses: voxels with few rays.
5. Item 7: PR #41 is open. Do not merge it. On a quiet machine run `uv run pytest tests/` in
   `icenine_py` and put the count line in the PR body.

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
| AR1 | recall of the ray-indexing proposer | built and measured; in review; branch `feature/eval-reduction-ar1` |
| item 7 | open the PR to `develop` (merge waits for the owner) | **open**: PR #41 (2026-11-09). Not merged. Add a full-suite test count to its body before asking for the merge |
| item 14 + 3 | both keys on E0 and the pilot, timing A/B | not started; needs a quiet machine |
| item 18 | does C++'s clang fuse the clip's multiply-adds | **done** 2026-11-09: yes (below) |
| C++ time | measured C++ time on the pilot region (owner's question on finding 6) | not started |
| AR2 | indexing as the BFS seed search | next, after the AR1 merge |
| AR3 | index the layer and paint grains | after AR2's first result |
| AR4 | the learned proposer | waits for AR1 |
| AR5 | the stop rule | waits for AR1; first job: a per-level record of full searches on full-sample voxels (the pilot's 69), then the AR0 stop grid on it |
| item 5 | full-sample baseline run, `FusedCost` on | after item 14 |

## Running jobs

Agents of the session of 2026-11-09 (they die with the session; `git worktree list` shows their
worktrees, and their branches keep the commits): the AR1 reviewer (it writes
`.claude/reports/20261109-review-ar1.md`). No detached compute job is running.

## Results so far

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

**AR1 (measured, not reviewed yet).** Page on the task branch:
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
  3 boundary voxels. With the finisher on the best 2: 0 of 69, at a mean of 1,926 evaluations; this
  arm was added after the first result was seen.
- E0 (the easy case: its images hold only each voxel's own spots): 200 of 200 at every K; 0 wrong
  against 41 stored; median 530 evaluations against 48,724.
- Indexing one voxel: median 0.18 s; one-time pairing of all spots: 17 s. Both contended. No
  wall-time comparison yet.
- The gate for AR2 (recall at 20 of at least 90%) is met.

Note on the `.md` refusal: the AR0 fix agent's Write of a new report file was allowed the same
afternoon. The refusal is not consistent; the rule in "Next action" still holds when it happens.

**Item 18.** The C++ build on this machine fuses the clip's multiply-add: Apple clang 21.0.0 (arm64) at
`-O3 -DNDEBUG` compiles `p0 + slope * (m - p0)` of `SutherlandHodgman.h` to one `fmadd` (single
precision); with `-ffp-contract=off` it is `fmul` then `fadd`. Test: `scripts/autoresearch/item18/clip_fma.cpp`.
Not checked: GCC, x86.

## For the owner

- Still yours: the merge to `develop` (item 7), items 12, 22 and 23, and Ilium's plan approvals.
- `shmem ask` fails here: the API key in `~/.shmem/.env` is refused and the local embedder is
  not running.
