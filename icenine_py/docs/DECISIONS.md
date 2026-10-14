# Decision log (owner decisions; shared by Ilium Agent and Newt Agent)

One log for both agents (ROADMAP protocol rule 6, owner 2026-10-14). Rules:

- **What goes here:** every decision the owner makes: plan approvals, go / no-go calls, default
  flips or "stays off", order and priority, scope cuts, keep / drop of data, protocol rules.
  An agent's own implementation choices are not decisions; they go in reports and MIGRATION_HISTORY.
- **When:** in the session where the owner decides, before acting on it. Commit and push it
  (ROADMAP push timing (d)).
- **Append only.** New entries go at the end. An entry is never edited or deleted; a changed decision
  gets a new entry that names the one it supersedes. Git merges this file with the union driver
  (`.gitattributes`), so appends on both branches merge without conflicts.
- **ID:** `D-YYYYMMDD-<agent><n>`, agent `I` (Ilium) or `N` (Newt), n counting from 1 per day and
  agent, so IDs never collide.
- **shmem:** save every entry you record as one atomic fact,
  `shmem add --source user "Owner decision <ID> (<date>, IceNine <track/row>): <decision>"`. On read
  on fetch, save the other agent's new entries to your own shmem too (each machine has its own store).
- **Entry format:**

  ```
  ### <ID> · <track/row> · <one-line decision>
  - <date>, owner. Recorded by <agent>. Supersedes: <ID or none>.
  - <the decision, and the reason when the owner gave one>. Where: <ROADMAP row, plan or doc>.
  ```

The entries up to D-20261014-I3 were backfilled on 2026-10-14 from the owner decisions marked in
ROADMAP.md (rows and "Decisions waiting on the owner"). Older decisions are in MIGRATION_HISTORY.md.

## Entries

### D-20261006-I1 · Q-F1 · F1 and F1b stay off by default
- 2026-10-06, owner; confirmed 2026-10-13. Recorded by Ilium Agent (backfill). Supersedes: none.
- The CSL-relative FindOptimal fixes F1 and F1b stay opt-in, off by default. Where: ROADMAP Q-F1;
  `todo_findoptimal_wrong_candidates.md`.

### D-20261012-I1 · all tasks · Rules for every task
- 2026-10-12, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- New switches are opt-in and off by default; no default flips and no golden regeneration inside a
  task (list them for the owner); quote evaluations and wrong counts side by side; no p-values; clean
  full-sample images and the E0 200-voxel set until the noise model is redone (Q8). Where: ROADMAP
  "How to use this file".

### D-20261012-I2 · Q9 · Ray rule pixel counting in cases b and e deferred
- 2026-10-12, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- Deferred. Where: ROADMAP Q9; finisher findings §7.9.

### D-20261013-I1 · protocol · Between-agents protocol, rules 1-4
- 2026-10-13, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- The agents are named Ilium Agent and Newt Agent and talk only through ROADMAP.md: own mailboxes,
  read on fetch, push timing, own rows only. Where: ROADMAP "Between agents".

### D-20261013-I2 · Tracks E and S · Per-task review rule
- 2026-10-13, owner (Track E; Track S adopted it the same day). Recorded by Ilium Agent (backfill).
  Supersedes: none.
- The main session plans and the owner approves the plan; an `implementer` builds it on a task branch;
  `/review-task` (Opus) reviews it when done; required fixes go to a fresh implementer with no owner
  review in between. Where: ROADMAP Tracks E and S.

### D-20261013-I3 · protocol · Rule 5: clear context between roadmap items
- 2026-10-13, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- When an item is finished (merged, pushed, mailbox posted): `/handoff`, then the owner runs `/clear`
  and `/resume`. Where: ROADMAP protocol rule 5.

### D-20261013-I4 · A-exact · `CppEdgeRules` stays opt-in; A-exact before Task B
- 2026-10-13, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- Where: ROADMAP Track E (commit 948adf8).

### D-20261013-I5 · S2, S3 · A GPU path is wanted for the search cost
- 2026-10-13, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- S3 (batched torch) goes before S2 (fused C); S2 stays open. Where: ROADMAP S2, S3.

### D-20261013-I6 · Q3 · `CppOverlapRule` default stays off
- 2026-10-13, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- Where: ROADMAP Q3.

### D-20261013-I7 · Q8 · Remove the random spot drop; Ilium Agent does Q8
- 2026-10-13, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- Remove the unphysical 10% random spot drop (keep jitter, hot pixels, blobs), re-render and re-run the
  realistic comparisons; Ilium Agent does the Q8 re-render. Where: ROADMAP Q8.

### D-20261013-I8 · X1 · The 500-grain runs come later
- 2026-10-13, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- Track X stays deferred. Where: ROADMAP Track X.

### D-20261013-I9 · S1, S3 · S3 first; S3 plan approved with draw option A
- 2026-10-13, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- S3 goes before S1 (it does not wait for Newt Agent's A-exact); the S3 plan is approved with draw
  option A (exact). The draw option is superseded by D-20261014-I1. Where: ROADMAP S1, S3;
  `docs/todo_batched_quick_mc.md`.

### D-20261013-I10 · Q8 · Q8 plan approved
- 2026-10-13, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- Re-runs on the 200-voxel set only, after S3. Where: ROADMAP Q8; `docs/todo_q8_noise_nodrop.md`.

### D-20261014-I1 · S3 · Draw option B; Steps 1-5 approved
- 2026-10-14, owner. Recorded by Ilium Agent (backfill). Supersedes: D-20261013-I9 (draw option).
- S3 Step 0 found option A a no-go: speculation costs 21-67x the counted evaluations. Option B gives
  each candidate a fixed slice of uniforms, so results under the key are compared by wrong counts.
  Steps 1-5 approved. Where: ROADMAP S3; `docs/todo_batched_quick_mc.md`.

### D-20261014-I2 · S3 · `BatchedQuickMC` stays off by default
- 2026-10-14, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- The default flip is decided by a later large-scale test that includes it. Where: ROADMAP S3.

### D-20261014-I3 · S3 · The Step 0 traces stay committed
- 2026-10-14, owner. Recorded by Ilium Agent (backfill). Supersedes: none.
- `icenine_py/benchmarks/s3_step0/` (about 940 KB of JSON) stays in the repository. Where: ROADMAP S3.

### D-20261014-I4 · protocol · Rule 6: one shared decision log, saved to shmem
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- Every owner decision is logged in this file, shared by both agents, and saved to shmem. Where:
  ROADMAP protocol rule 6; this file.

### D-20261012-N1 · ray rule (finisher study) · Ray-rule windows: half-widths, left as set; levels 0-2 only
- 2026-10-12, owner. Recorded by Newt Agent (backfill, 2026-10-14). Supersedes: none.
- Each term of the arm A window is a half-width (the intended reading); both windows (R<n> pixels, arm A
  formula, d_omega 2 deg) stay as set, untuned; R2 and A go to step 3. The ray rule is used at levels 0-2
  only; level 3, FindOptimal, VarianceMinimizing and the final cost stay strict C++ (revisit only if step 3
  loses the truth at level 3). Case f1 dropped; case f always forgiven. Do not re-tune the windows without
  asking. Where: `docs/findings_ray_overlap_rule_2026-10.md`; code comments "WINDOW DECISION"
  (`icenine/ray_overlap.py`, `icenine/_rasterize.c`).

### D-20261012-N2 · all tasks · No extra seeds to firm up comparisons
- 2026-10-12, owner. Recorded by Newt Agent (backfill, 2026-10-14). Supersedes: none.
- Do not propose re-running with more random seeds to firm up a comparison: the search uses randomness to
  explore, not as independent samples. Where: this file (no ROADMAP row).

### D-20261013-N1 · A · `BatchedScreen` on by default
- 2026-10-13, owner. Recorded by Newt Agent (backfill, 2026-10-14). Supersedes: none.
- `BatchedScreen` defaults to 1 (config key and `SearchParameters`; `screen_chunk` 32; `BatchedScreen 0`
  restores the per-orientation screen), on condition that no pinned test or golden moves. The first flip
  moved a profiling test and was reverted; the flip was made once the profiling instrument counted
  `ScreenEvaluator` calls. `CppEdgeRules` and `CppOmegaInterp` stay off (D-20261013-I4). Where: ROADMAP
  Track E row A; MIGRATION_HISTORY "BatchedScreen as the default: made".

### D-20261014-N1 · Track E · Review on the task branch before the merge
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- `/review-task` runs on the task branch before `finish_task`; the required fixes go to a fresh implementer
  on the same branch (no separate fixes branch), then the task is merged. Refines the order of
  D-20261013-I2. Where: ROADMAP Track E.

### D-20261014-N2 · B · Task B report accepted
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- Task B (anatomy of the cost function's edges) accepted after the review fixes, as merged in 6561d26.
  Where: ROADMAP Track E row B; `docs/findings_cost_edges_2026-10.md`.

### D-20261014-N3 · C · Merge Task C before its code review; status page
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none (a one-off exception to D-20261014-N1
  for Task C only).
- Task C merged into `feature/eval-reduction` (e609c74) before its `/review-task`; the review runs on
  the merged commits e5a2cbb..e609c74. The owner also asked
  for a plain-English status page of Track E, kept in the branch:
  `docs/status_eval_reduction_2026-10-14.md`. Where: ROADMAP Track E row C; that status page.
