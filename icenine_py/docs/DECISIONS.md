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
### D-20261012-I3 · ray rule (finisher study) · The ray rule may become the default
- 2026-10-12, owner. Recorded by Ilium Agent (backfill, 2026-10-14). Supersedes: none.
- "The ray rule as the default is fine": it ties the lenient rule on clean images and costs 6% more
  evaluations than the C++ rule on full-sample voxels, against 46% for the lenient rule. The code
  default was not changed (still off); the flip is a follow-up. Where:
  `docs/findings_ray_overlap_rule_2026-10.md` (decisions table).

### D-20261014-I5 · Q8 · Studies stay on clean images
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- The "clean images only" rule is not lifted after Q8: the noisy images are not needed for the
  studies. Where: ROADMAP "How to use this file"; Q8 row.

### D-20261014-I6 · Q8 · The realistic noise default becomes no drop
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: the "no default flips" rule of the Q8 plan.
- The realistic noise models default to no random spot drop (`p_miss 0`; jitter, hot pixels and blobs
  kept), in Phase D (`NoiseParams`) and E0 (`RealismConfig`). Done on the Q8 task branch. Where:
  ROADMAP Q8; `docs/todo_q8_noise_nodrop.md`.

### D-20261014-I7 · Q8, X1 · No Phase D runs on the noisy images for now
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- The Phase D pilot and reseed are not re-run on `full/realistic_nodrop`; noisy-image runs are
  skipped for now. Where: ROADMAP Q8, X1.

### D-20261014-I8 · Q8 · Check why noisy images are easier to reconstruct
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- Owner's experience: noisy images give more candidates and a wider cost function, so the search is
  easier. Check it where the data does not already show it (every rule has fewer wrong voxels on the
  no-drop E0 images than on clean). Where: ROADMAP Q8.

### D-20261014-I9 · protocol, Q8 · The clean-images rule loses its "until Q8"
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- The shared ROADMAP lines that limited studies to clean images "until Q8" are reworded to state the
  rule without an end (D-20261014-I5). Where: ROADMAP "How to use this file"; Track E left-overs.

### D-20261014-I10 · Q3 · The Q3 headline is restated like for like
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- With `CppOverlapRule` on, E0 wrong goes 35 to 61 clean and 27 to 42 realistic of 200 (centroid
  screen on both sides), not the first-quoted 40 to 61 and 25 to 42 (rule-off runs with the vertex-0
  screen). Where: `docs/findings_finisher_mc_study_2026-10.md` §3(b) and §7.3.

### D-20261014-I11 · Q8 · E0 default readers stay on the old variant
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- A default E0 run writes `realistic_nodrop` results while the stored-result readers (`C.VARIANTS`)
  read `clean` and `all` by name; accepted as is. Where: `scripts/findoptimal_robustness/common.py`.

### D-20261014-I12 · Q8 · Other study scripts and configs keep the drop variant
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- The five study scripts with their own `VARIANTS = ["clean", "all"]` and the six
  `*_realistic*.config` Phase D arms are not moved to no drop. Where: MIGRATION_HISTORY "Q8".

### D-20261014-I13 · Q14 · Noise against neighbour spots: open question, not run now
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- Whether the detector noise or the neighbour/twin spots make the no-drop E0 images easier (Q8 Step 7)
  is logged as open question Q14; the `noise_nodrop` split is not run now. Where: ROADMAP Q14.
### D-20261014-N4 · C, G · Task C's filters are not tried in Task G for now
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- Task G starts from the baseline (`CandidateDedup` off). The two filters that looked harmless on
  full-sample voxels but lose E0 voxels (`clique=cpp:1`, `pre=0.5`) are not tried in Task G's BFS
  pilot for now. Where: ROADMAP Track E row C; `docs/status_eval_reduction_2026-10-14.md`.

### D-20261014-N5 · A · Delete the untracked Task A memory-run log folder
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- The untracked `benchmarks/eval_reduction/A/memory/cache/` (one run log and a done marker) was
  deleted. Where: this entry.

### D-20261014-I14 · S1 · S1 is the next Track S item
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- After Q8 and the Task C merge, the next Track S item is S1 (remove the per-call Python loops in
  `evaluate`; bit-identical). It starts with a plan, `docs/todo_s1_evaluate_loops.md`, which the
  owner approves before an implementer builds it. Where: ROADMAP Track S row S1.

### D-20261014-I15 · S1, S3 · S1 plan approved, with Step 4
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- The S1 plan (`docs/todo_s1_evaluate_loops.md`, Steps 0-6) is approved as written, including Step 4:
  `CppEdgeRules` becomes supported in the batched local cost (S3's open follow-up), so runs with both
  opt-in keys `BatchedQuickMC` and `CppEdgeRules` use the lockstep quick MC instead of the
  sequential fallback; no default changes. Where: ROADMAP Track S rows S1 and S3.

### D-20261014-I16 · S1 · Repeat the S1 timings on a quiet machine
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- The contended S1 timings are not accepted as final: the review-fix implementer re-runs the S1
  timing script (`run_s1_timing.sh`) only if the timing preflight passes; otherwise the repeat is
  listed as left undone and the docs keep the "contended" label. Where: ROADMAP row S1;
  MIGRATION_HISTORY "S1: per-call loops out of evaluate".

### D-20261014-I17 · S1, Track E scripts · Keep S1's edits to three Track E scripts
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- S1's edits to `scripts/eval_reduction/cost_edges.py`, `parity_common.py`, `screen_timing.py` and
  one check in `tests/test_cost_edges.py` (they also wrap the new `scattering_omegas_fixed_beam`,
  because `evaluate` no longer calls `get_scattering_omegas_torch`) are kept; Newt Agent is told in
  the merge heads-up. Where: Ilium mailbox; MIGRATION_HISTORY "S1".

### D-20261014-I18 · S1 · The evaluate fixture is strict only on its own platform
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- `tests/data/s1_evaluate_ref.npz` is compared exactly on the platform and torch version that made
  it (arm64, torch 2.8.0); on any other the test skips with a message naming the difference.
  Where: `tests/test_s1_evaluate.py`.

### D-20261014-N6 · D · Task D Step 4 (geometric-feature classifier) is skipped
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- The small classifier on pair-geometry features (plan Task D Step 4) is not built. The plan's trigger
  (single-source recall below 0.95) is undefined because almost no spot has one source; one-grain
  recall is 0.990. Where: ROADMAP Track E row D; `docs/findings_spot_pairing_2026-10.md`.

### D-20261014-N7 · D · Task D read as promising: spot-to-grain association is enough
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- Associating spots to grains (not to voxels) is sufficient: (1) the grain's spots give a starting
  orientation, and the area around the grain can be sampled with it; (2) from that start, run the usual
  FindOptimal search, or other searches outward from the orientation; (3) the orientations defined by
  the combinations of peaks and grains give a starting candidate list that may be much smaller. These
  are follow-ups to plan; none is started. Where: ROADMAP Track E row D.

### D-20261014-N8 · R · Learned spot-association model is parked as a research idea
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- A different classifier from Task D Step 4, needing much more planning: simulate images as usual with
  each peak labelled by its originating grain; build or fine-tune a ConvNet or ViT that looks at the
  pair of detector images and predicts the association given the second image. The model should learn
  (1) the origin-independent affine transformation of peak shapes between the detectors and (2) the
  association of peak shapes. Not started; needs its own plan. Where: ROADMAP Track R.

### D-20261014-N9 · D · Task D review fixes all applied; the indexing proposal is kept, corrected
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- All 16 required fixes of the Task D review are applied by a fresh Opus implementer, then the branch
  is reviewed again and fixed before it is merged. The "Proposal, not tested" block of the findings doc
  is kept with a corrected error budget, and the main session's notes (how to test the indexing
  proposal; notes for planning the learned model of N8) are added. Where: ROADMAP Track E row D;
  `docs/findings_spot_pairing_2026-10.md`.

### D-20261014-N10 · D · The indexing error budget leads with the measured ray-direction error
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- In the findings doc's "Proposal, not tested" block, the error budget leads with the measured
  ray-direction error on clean images (median 0.0044 degrees, p90 0.0155), with one pixel kept as the
  tolerance and upper bound; on clean images the 1 degree omega frame term dominates. The second
  review's 10 required fixes are applied, then checked again before merge. Where: ROADMAP Track E
  row D; `docs/findings_spot_pairing_2026-10.md`.

### D-20261014-N11 · D · Planning the grain-seeding follow-up (Task D2) is deferred
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- No plan is written now for the N7 follow-up (seed orientations by indexing each grain's matched
  spots, then FindOptimal from the seeds). It is deferred; the test design stays in the findings doc,
  "How to test it". Where: ROADMAP Track E row D; `docs/findings_spot_pairing_2026-10.md`.

### D-20261014-N12 · R · The learned spot-association model (N8) is listed in ROADMAP Track R
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- The parked ConvNet/ViT association model of D-20261014-N8 is added to the owner's Track R row, as a
  parked research idea that needs its own plan. Where: ROADMAP Track R.

### D-20261014-N13 · E (all tasks) · Broken-voxel rule: only breaks traced to the arm count
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: the plan's rule "no broken voxels beyond the
  baseline's count" (eval-reduction plan, Task C decision point; open question from the Task C review).
- An arm that changes the candidate list or the random draws flips E0 voxels by the change of draw alone
  (14-23 each way between seeds, `s3_ab_paired`), so no such arm could meet the old rule. From now on each
  E0 voxel an arm breaks is traced, and the arm fails only if a break is caused by the arm itself, not by
  the change of random draw. Reports give traced-to-arm, traced-to-draw and untraced breaks separately.
  Newt's default, open to the owner: an untraced break counts against the arm until traced. Task C's
  outcome is unchanged (N4). Where: plan "Ground rules for every task"; ROADMAP Track E row C.
### D-20261014-I19 · S4 · S4 is next, widened to a correctness check and a loop-ordering study
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none (replaces the S4 plan's Step 0
  go/no-go gate, which assumed a per-voxel screen at the seeds only).
- S4 is the next Track S item. Two questions lead it: (1) correctness: a translated voxel's diffraction
  and lit pixels must be the original voxel's pixels moved by the per-peak transform of the voxel
  translation; (2) speed: which loop ordering gives the best throughput for the level-0 grid search
  when it runs one orientation at a time over every voxel position (or over detector positions),
  instead of today's per-voxel screen followed by per-voxel optimisation. Investigate and plan; the
  implementer is Opus 5.5 at high effort, for the investigative nature. Where: ROADMAP Track S row S4;
  `docs/todo_translation_screen_prototype.md`.

### D-20261014-I20 · S4 · Fewer voxels per test, and run the work in parallel
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none (amends the S4 plan of I19).
- To cut S4's total run time: test fewer voxels per layer. That means 20 Phase D reference voxels instead of 50,
  10 timed `BatchedScreen` voxels instead of 20, and orderings timed on a 2,000-voxel patch with a linear-scaling check and
  a projection to the full layer, which is run only for the best ordering. Exploit parallelism: the
  reference, the correctness checks and the equality runs run concurrently in the background (up to 6
  processes); timings run alone. Where: `docs/todo_translation_screen_prototype.md`.

### D-20261014-I21 · S4 · S4 plan approved; implement now
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- The S4 plan (`docs/todo_translation_screen_prototype.md`, Steps 0-5, with I19 and I20) is approved
  as written; an Opus 5.5 implementer at high effort builds it on
  `feature/python-speed-parity-s4-translation-screen`. Where: ROADMAP Track S row S4.

### D-20261016-I22 · S4 · Hybrid library path parked; S4 stays unmerged for now
- 2026-10-16, owner. Recorded by Ilium Agent. Supersedes: none.
- The S4 hybrid library path (BFS as today, level-0 candidates from one global screen and a stored
  layer-wide re-cost table for seeds and failed voxels) is parked. Nothing else is started for now:
  S4 (review fixes applied, 37b49a7) stays on `feature/python-speed-parity-s4-translation-screen`,
  unmerged, and the projections are not redone under CppOverlapRule. Where: ROADMAP Track S row S4;
  `.claude/reports/20261015-0056-s4_steps0-5.md`.
### D-20261014-N14 · E · Task E implementer is Opus, not Sonnet
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none (an exception, for Task E only, to the plan's
  "implementers are fresh Sonnet agents" and the standing feedback rule).
- Task E (smoothing the cost) is implemented by a fresh background Opus implementer. Review fixes stay with
  fresh Opus implementers as before. Where: eval-reduction plan "Coordination"; ROADMAP Track E row E.

### D-20261018-I1 · H12 (O1) · Decision IDs keep the written rule: n restarts per day and agent
- 2026-10-18, owner. Recorded by Ilium Agent. Supersedes: none.
- IDs stay `D-YYYYMMDD-<agent><n>` with n counting from 1 per day and agent, as the header says
  ("keep current rule works"). D-20261016-I22, which used a running counter, stays as written (append
  only). `log_decision.sh` follows this rule. Where: DECISIONS.md header;
  `docs/todo_ilium_workflow_tooling.md` O1.

### D-20261018-I2 · H12 (O2, O5-O8) · Workflow tooling plan approved; implement now in three tasks
- 2026-10-18, owner. Recorded by Ilium Agent. Supersedes: none.
- The plan `docs/todo_ilium_workflow_tooling.md` is approved with its grouping: (a) `log_decision.sh`
  and `roadmap_row.py`, (b) `verify_handback.sh` and `read_on_fetch.sh`, (c) `make_brief.sh` and the
  `/decide`, `/implement-plan`, `/apply-fixes` skills; implement now, each task reviewed by
  `/review-task`. The other open choices as the plan proposes ("rest of the decisions are fine"): O5
  two timing-label fields (`workers`, `machine`); O6 task branches off `feature/python-speed-parity`;
  O7 `roadmap_row.py` edits Track S rows only; O8 `read_on_fetch.sh` prints the shmem lines and runs
  them only with `--save`. Where: ROADMAP Track H row H12; the plan.

### D-20261018-I3 · H12 (O3) · ROADMAP row H12 placed as proposed
- 2026-10-18, owner. Recorded by Ilium Agent. Supersedes: none.
- The owner agreed to the proposed Track H row H12 for the workflow tooling. Ilium Agent places it
  on `feature/python-speed-parity` with the tooling (an owner exception to protocol rule 4, since
  `feature/housekeeping` is not started) and keeps it up to date. Where: ROADMAP Track H row H12.

### D-20261018-I4 · all tasks (O4) · At most 6 heavy processes on Ilium's 16 GB machine; exceptions allowed
- 2026-10-18, owner. Recorded by Ilium Agent. Supersedes: none.
- The limit of at most 6 heavy processes at once on this 16 GB machine stands as the standing rule
  ("this process limit works fine"); a plan or an owner decision may grant an exception. It is kept
  per machine (`git config --local icenine.maxHeavyProcs 6` here), so Newt Agent's machine sets its
  own. Where: `docs/todo_ilium_workflow_tooling.md` O4 and O9; study-conventions (owner edit).

### D-20261018-I5 · H12 (O9) · Tool instructions live in tracked files shared with Newt Agent
- 2026-10-18, owner. Recorded by Ilium Agent. Supersedes: none.
- The instructions to use the tools go into tracked files (study-conventions, the skills, one sentence
  each in protocol rules 2 and 6, the README table), not into either machine's memory; the tools are
  agent-neutral (`icenine.agent`); the implementer drafts the study-conventions and protocol text and
  the owner commits it. Where: `docs/todo_ilium_workflow_tooling.md` "Instruction changes and sharing
  with Newt Agent".
### D-20261018-N15 · E, G · Task E's remaining 24 full-sample voxels wait for Task G
- 2026-10-18, owner. Recorded by Newt Agent. Supersedes: none.
- Task E tested the smooth ranking (`SmoothCostRank`) on the six interior full-sample voxels only (cut from 30 for
  time); all six are right with the strict search, so they can show harm but not a gain. The other 24 voxels of the
  30-voxel set, which include voxels the strict search gets wrong, are run in Task G (carried arms E3+E4 tau 0.25 and
  E1+E2+E3+E4, lenient and C++ rule, paired, breaks traced per N13), before the smooth key is chosen. Task E merges
  without them. Where: eval-reduction plan Task G step 2b; `docs/findings_cost_smoothing_2026-10.md` "Limits".

### D-20261018-N16 · E · Task E merges without a quiet-machine re-timing
- 2026-10-18, owner. Recorded by Newt Agent. Supersedes: none.
- The Task E review (`.claude/reports/20261018-2151-review-E.md`, fix 1) asks for whole-voxel times, key off and on.
  The fix implementer re-times only if the timing preflight passes; otherwise the docs keep the "contended" label and
  the repeat is listed as not done. The merge does not wait for it. Where: ROADMAP Track E row E; MIGRATION_HISTORY
  "Task E".

### D-20261018-N17 · E, G · The per-image-set SmoothCost cache is built in Task G
- 2026-10-18, owner. Recorded by Newt Agent. Supersedes: none.
- `SmoothCost` is rebuilt on every `reconstruct_voxel` call. Caching it per image set is done in Task G, if the smooth
  ranking key is used there. Task E's docs lead with the as-shipped cost (rebuild per voxel) and give the
  one-build-per-image-set figure as a projection. Where: eval-reduction plan Task G step 2c.

### D-20261021-I1 · H12 · Commit the workflow skills and the owner text; do the small tooling follow-ups
- 2026-10-21, owner. Recorded by Ilium Agent. Supersedes: none.
- The owner approved committing the skills /decide, /implement-plan and /apply-fixes and the
  drafted owner text (implementer.md, study-conventions.md, the caffeinate wrapper in the Bash
  hook, one sentence each for ROADMAP protocol rules 2 and 6), committed by Ilium Agent at the
  owner's request. Also approved: fix read_on_fetch.sh listing the caller's own merged commits as
  the other agent's, let make_brief.sh accept a plan outside the repository, and update the S4
  row. Where: ROADMAP H12; .claude/reports/20261018-owner-text-workflow-tooling.md (Ilium's
  machine).

### D-20261021-I2 · S2 · S2 is the next Track S item; plan it
- 2026-10-21, owner. Recorded by Ilium Agent. Supersedes: none.
- Go for S2, the fused C evaluation behind an opt-in key. Ilium Agent writes the plan for owner
  approval before any implementer runs. Where: ROADMAP Track S row S2; plan to come in
  docs/todo_s2_fused_c_eval.md.
### D-20261021-N1 · Track E row F/G · Compare every Task F carry option side by side before choosing
- 2026-10-21, owner. Recorded by Newt Agent. Supersedes: none.
- Before choosing what Task G carries, run every option in parallel and summarize capture, timing
  and evaluation counts so the owner can compare: F4 (1-frame window), F4 (5-frame window), F5, F5
  then F4, the policy 'F4 then F5+F4 on rejection' with each window, and CMA-ES (sigma0 1 and 2
  deg). The runs use voxels Task F did not use, so the policy and window numbers are out of
  sample. The choice of carry (open from Task F) waits for this summary. Where: ROADMAP Track E
  row F; docs/findings_outward_search_2026-10.md.

### D-20261021-I3 · S4 · Merge S4 as a finished study; library path stays parked
- 2026-10-21, owner. Recorded by Ilium Agent. Supersedes: none.
- Option A: merge feature/python-speed-parity-s4-translation-screen (scripts, results, findings,
  review fixes 37b49a7; no library change) into feature/python-speed-parity as a finished study.
  The hybrid library path stays parked (D-20261016-I22); re-projecting under CppOverlapRule
  (option B) or reviving the library path (option D) waits until Task G settles the overlap rule
  and the carry options. Where: ROADMAP Track S row S4; docs/todo_translation_screen_prototype.md.

### D-20261021-I4 · S2 · S2 plan approved; Opus at high effort implements it
- 2026-10-21, owner. Recorded by Ilium Agent. Supersedes: none.
- The S2 plan docs/todo_s2_fused_c_eval.md is approved as written: Step 0 (go/no-go prototype,
  scripts only) is its own implementer run, and Steps 1-6 continue without a new approval on a go
  (median <= 30 us per call on the triangle arm and every difference traced); 30-50 us goes to the
  owner first. The implementer, including fix implementers, is Opus at high effort. Where: ROADMAP
  Track S row S2; docs/todo_s2_fused_c_eval.md.

### D-20261027-I1 · S2 · Rest of S2 delegated to Ilium Agent; owner reviews in full at the end
- 2026-10-27, owner. Recorded by Ilium Agent. Supersedes: none.
- From Step 1 on, Ilium Agent goes through the whole S2 plan with its own recommendations, records
  each decision here, and has every task reviewed and merged step by step (subagents, in parallel
  where possible). The owner does a full review with Ilium Agent when the plan is finished.
  Decisions I2-I5 and later S2 entries marked 'delegated' are Ilium Agent's under this delegation.
  Where: ROADMAP Track S row S2; docs/todo_s2_fused_c_eval.md.

### D-20261027-I2 · S2 · Parity rule: pixel-edge flips accepted (bound 0.01 px), real gates are C++ and E0
- 2026-10-27, owner. Recorded by Ilium Agent. Supersedes: none.
- Owner: a difference of 0.005 px does not matter much. A fused-versus-torch difference is
  accepted when it is a pixel-edge flip, with the torch or fused pixel coordinate within 0.01 px
  of an integer (a bug guard, not a precision target); frame, eta, observability or t>0 flips
  fail. The gates that matter are Task A agreement with C++ (at least torch's) and the paired E0
  comparison. Replaces the 8-ulp rule of the plan. Where: ROADMAP Track S row S2;
  docs/todo_s2_fused_c_eval.md.

### D-20261027-I3 · S2 · Delegated: portable C, no Accelerate
- 2026-10-27, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261027-I1: the fused kernel stays portable C. Accelerate's sgemm made the
  reflection vectors bit-identical to torch on macOS (129 to 44 differing orientations of 4,000)
  but does not remove the flips, ties the build to macOS and torch's BLAS, and the parity rule of
  I2 does not need it. Where: ROADMAP Track S row S2; docs/todo_s2_fused_c_eval.md.

### D-20261027-I4 · S2 · Delegated: Task A check with the prototype runs in parallel with Steps 1-2 and gates their merge
- 2026-10-27, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261027-I1, owner asked for parallel subagents: the Step 0 prototype is
  compared with C++ on Task A's 2,230 orientations now (scripts only), in parallel with the Steps
  1-2 library task; if the prototype agrees with C++ on fewer orientations than torch, Steps 1-2
  are not merged and S2 stops for the owner. Step 3 follows Steps 1-2; Steps 4-6 follow Step 3.
  Where: ROADMAP Track S row S2; docs/todo_s2_fused_c_eval.md.

### D-20261027-I5 · S2 · Delegated: plain operation order for g_lab tried in Step 1
- 2026-10-27, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261027-I1: Step 1 tries the plain order for g_lab (Step 0: matches torch
  on 0.929 of elements against 0.804 for the fmaf chain) and keeps it only if it gives fewer flips
  under the parity rule. Where: ROADMAP Track S row S2; docs/todo_s2_fused_c_eval.md.

### D-20261027-I6 · S2 · Delegated: Task A inputs accepted on aggregates; no C++ toolchain install
- 2026-10-27, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261027-I1: Task A kept neither its input files nor the C++ output, and the
  C++ toolchain (cmake, boost, eigen, open-mpi) is not installed. The regenerated inputs are
  accepted because labels and indices match on 2,230 of 2,230, the stored C++ product statistics
  are reproduced on 2,000 of 2,000 pair records, and torch on them gives the stored Task A Python
  aggregates on all 4 arms for every orientation. Installing the toolchain is a system change
  outside S2; it waits for the owner's full review. Where: ROADMAP Track S row S2;
  benchmarks/s2_taska; MIGRATION_HISTORY S2.

### D-20261027-I7 · S2 · Delegated: plain g_lab order kept; Step 4 compares the library path with C++ directly
- 2026-10-27, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261027-I1 and I5: Step 1 measured fewer flips with the plain order (93 of
  4,000 orientations against 129 with the fmaf chain, all pixel-edge flips), so the library uses
  it. The Task A C++ check ran on the prototype (fmaf order), so it does not carry over bit for
  bit; Step 4 reruns the Task A comparison on the library path with the order it ships, and the
  gate (agreement with C++ at least torch's) applies there again. Where: ROADMAP Track S row S2;
  docs/todo_s2_fused_c_eval.md Steps 4-6; benchmarks/s2_core.

### D-20261027-I8 · S2 · Delegated: Steps 1-2 review choices (header, fp-contract flag, detector contract, merge before the library C++ check)
- 2026-10-27, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261027-I1, on the Steps 1-2 review (APPROVE, 5 fixes): (1) merge with the
  plain g_lab order now; the C++ check on the library path is Step 4 (I7); (2) accept a small
  internal header for the shared stage-D declarations, added to the plan's library changes; (3)
  document the detector-tensor contract instead of per-call _version checks (the context's
  geometry check covers moves, per-call checks would cost every call); (4) put -ffp-contract=off
  in the extension's compile flags so GCC also keeps the float order, and state results are
  validated on macOS clang only; key-off identity is rechecked after the flag. Where:
  docs/todo_s2_fused_c_eval.md Steps 1-2; .claude/reports/20261027-2232-review-s2-core.md (Ilium's
  machine).

### D-20261027-I9 · S2 · Delegated: fp-contract off only in the fused file, not the whole extension
- 2026-10-27, owner. Recorded by Ilium Agent. Supersedes: D-20261027-I8.
- Ilium Agent under D-20261027-I1, correcting part (4) of I8; parts (1)-(3) stand. A global
  -ffp-contract=off also changes the key-off _rasterize.c: clang stops fusing the multiply-adds of
  the Sutherland-Hodgman clip, and triangle_overlap counts change on 106 of about 400,000
  border-crossing fuzz triangles. Key off must stay bit-identical (rules for every task,
  D-20261012-I1), so the global flag is removed; _fused_eval.c alone keeps fp contraction off (its
  pragma, plus a GCC-honoured per-file pragma if one is portable), and results are stated as
  validated on macOS clang only. Whether C++'s clang build fuses the same clip expression is
  recorded as an open question for the owner's review. Where: docs/todo_s2_fused_c_eval.md Steps
  1-2; benchmarks/s2_core.

### D-20261027-I10 · S2 · Delegated: Steps 1-2 merge; GCC stays unvalidated and documented
- 2026-10-27, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261027-I1: Steps 1-2 merge after review (APPROVE) and fixes, key off
  bit-identical again after I9 (s1_parity 4 x 15 modes; triangle_overlap fuzz 0 of 600,000
  differing). No GCC is installed here, so the guarded GCC pragma in _fused_eval.c is uncompiled;
  FusedCost is documented as validated on macOS clang only, and a GCC or Linux check is listed for
  the owner's review rather than installing a compiler. Where: ROADMAP Track S row S2; README
  FusedCost row; benchmarks/s2_core.

### D-20261028-I1 · S2 · Delegated: Step 3 design accepted; two wording items and a GPU check deferred to Steps 4-6
- 2026-10-28, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261027-I1, after the Step 3 review (APPROVE, 4 record fixes): accept the
  batched fused dispatch in BatchedCostEvaluator.evaluate_many (not inside _forward_compact, since
  the fused entry includes stage D), the new variant="fused", HAS_FUSED staying 1 with a guarded
  import for the batched entry, and the one declaration line added to _stage_d.h outside the
  listed scope (the header is where shared entry points live). Deferred to Steps 4-6: the fallback
  warning text in cost_functions.py (says torch evaluate; the batched fallback is the torch
  compact variant), README rows 86 and 91, and a Step 5 check of whether GPU_MIN_BATCH (64,
  measured against the torch compact variant) should still send batches to the GPU when FusedCost
  is on. Where: ROADMAP Track S row S2; plan todo_s2_fused_c_eval.md Steps 4-6;
  .claude/reports/20261028-0930-review-s2-batched.md.

### D-20261028-I2 · S2 · Delegated: parity checker keeps the I2 wording; GPU_MIN_BATCH and BatchedQuickMC-with-FusedCost go to the owner's review unchanged
- 2026-10-28, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261027-I1, after the Steps 4-6 review (APPROVE, 8 record fixes): (1)
  s2_lib_parity.py drops pixel_negative from its accepted causes, matching s2_core_parity.py and
  the D-20261027-I2 wording (d = 0 is pixel coordinate 0.5, the centre of pixel 0, not an edge);
  no pixel_negative case occurred, so no result changes. (2) GPU_MIN_BATCH stays 64: with
  FusedCost on, MPS was 2.4-6.8x slower than the fused CPU entry at B 64-512 and CUDA is
  unmeasured; whether a GPU SearchDevice should yield to the fused CPU entry is listed for the
  owner's full review, not changed (no default flip, D-20261012-I1). (3) BatchedQuickMC with
  FusedCost on gives no per-evaluation gain (E0 20 vs 20 us); no recommended-configuration change,
  listed for the owner's review. Where: ROADMAP Track S row S2; MIGRATION_HISTORY S2;
  .claude/reports/20261028-2100-review-s2-results.md.

### D-20261028-I3 · H12 · While the owner reviews S2, Ilium does the S2 tooling follow-ups; Opus implements
- 2026-10-28, owner. Recorded by Ilium Agent. Supersedes: none.
- Owner: 'continue through the rest of the plan while I review' and 'use opus to implement'. The
  S2 plan has no steps left, so Ilium Agent takes the tooling gaps S2 found
  (todo_ilium_workflow_tooling.md, 'Follow-ups found in S2'): task (d) make_brief --section
  without another section's markers and verify_handback's ancestor-hash flag and test summary;
  task (e) sync_doc_tables --section and the 21 over-long lines in the S4 scripts. The two tasks
  share no files and run in parallel (3 heavy processes each). Opus implements; /review-task
  (Opus) reviews; no library change. Not started: S5 and S6 (need a plan and the owner's pick),
  the quiet-machine timing repeats of S1 and S3, and FusedCost's default (owner's review). Where:
  ROADMAP Track H row H12; docs/todo_ilium_workflow_tooling.md.

### D-20261028-I4 · H12 · Delegated: task (e) fix round also fixes audit_numbers' code-fence headings; keep --tables with several files; s2_e0_trace gets a marker
- 2026-10-28, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261028-I3, after the task (e) review (APPROVE, fixes): (1) keep commit
  1732831 (--tables takes several files; S2 Steps 1-6 draws tables from two files); (2) the fix
  round also makes audit_numbers.find_section skip headings inside code fences, sharing the
  heading scan with sync_doc_tables (today it cuts 6 MIGRATION_HISTORY and 9 README sections short
  without saying so; no S2 or S4 section affected), outside task (e)'s listed items but the same
  family and about 10 lines; (3) the S2 table s2_e0_trace, generated but referenced by no marker,
  gets a marker in MIGRATION_HISTORY S2, written by the main session at merge with the tooling
  completion entry. Where: ROADMAP Track H row H12; docs/todo_ilium_workflow_tooling.md task (e);
  .claude/reports/20261029-0000-review-tooling-e.md.

### D-20261028-I5 · H12 · Delegated: task (d) fix round; section-scoped markers in plans without # Steps; .claude/ items left for the owner
- 2026-10-28, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent under D-20261028-I3, after the task (d) review (APPROVE, 2 fixes): (1) the stop rule
  keeps falling back to the whole plan (taking a stop rule from elsewhere only tightens a brief);
  (2) in a plan with no '# Steps' heading, make_brief --section treats every '## ' section as its
  own task's: a library-changes or process-limit marker inside another '## ' section is never
  used, only markers in the named section or outside every '## ' section, with a stderr note
  naming any marker it skipped (stricter than review fix 1, which asked only for the note); (3)
  verify_handback: no count line with exit 0 is a FLAG (review fix 2); (4)
  scripts/dev/finish_task.sh drops its extra -q so its pytest output keeps the count line. Left
  for the owner's review, because .claude/ files are committed only by the owner: whether the
  implement-plan and apply-fixes skills show make_brief's stderr notes, and the
  run-tests-on-change.sh hook's extra -q. Where: ROADMAP Track H row H12;
  docs/todo_ilium_workflow_tooling.md task (d); .claude/reports/20261029-0010-review-tooling-d.md.
### D-20261028-N1 · Track E (all remaining work) · Opus implements every remaining task; one Fable review at the end
- 2026-10-28, owner. Recorded by Newt Agent. Supersedes: none.
- From now on every implementer for Track E is Opus (fresh agent per task, fixes included); the
  Sonnet default of the plan no longer applies. The per-task /review-task loop is dropped for the
  remaining work: once everything is done, one review by a Fable reviewer covers the carry-option
  comparison (Fcompare) and Task G, and replaces the plan's Review stage. Where:
  ~/.claude/plans/eval-reduction.md Task G (owner note) and Review stage.

### D-20261028-N2 · Track E row G · Task G plan approved; Newt runs the rest of the plan without sign-off between sub-tasks
- 2026-10-28, owner. Recorded by Newt Agent. Supersedes: none.
- The Task G plan in ~/.claude/plans/eval-reduction.md (sub-tasks G1, G2, G3, with its Library
  changes and Stop rule paragraphs) is approved; implement now. Newt Agent runs the rest of the
  plan that Ilium Agent has not done, merging each sub-task without asking; decision points take
  the plan's stated default, each recorded as a Delegated entry. Where:
  ~/.claude/plans/eval-reduction.md Task G; ROADMAP Track E row G.

### D-20261028-N3 · Track E (merge of python-speed-parity) · One-time pre-commit bypass for the merge of Ilium's branch
- 2026-10-28, owner. Recorded by Newt Agent. Supersedes: none.
- The merge of origin/feature/python-speed-parity through de3e29d into feature/eval-reduction
  (5e15bcd) was committed by the owner with --no-verify: the pre-commit hook refused the
  owner-approved .claude/ files (D-20261021-I1) and lines over 100 characters in Ilium's S2/S4
  scripts. One time only; agents still never bypass the hooks. Where: merge commit 5e15bcd;
  scripts/dev/precommit_check.py.

### D-20261028-N4 · Track E rows F, G · Delegated: the outward search carried into Task G is the F fallback policy, 1-frame window
- 2026-10-28, owner. Recorded by Newt Agent. Supersedes: none.
- Delegated (owner 2026-10-28: Newt decides the plan's decision points). From the carry-option
  comparison (D-20261021-N1, 24 fresh voxels per set): Task G uses the policy F4, then F5 then F4
  where F4's acceptance test rejects, with F4's frame_window at 1 (library default). It reached
  480/480 within 1 deg at every offset on both sets, and the 1-frame window ends within 0.05 deg
  more often on the full sample (85-86% against 77-78%). CMA-ES is not carried (34-57% from 5
  deg). The owner can override. Where: docs/findings_outward_search_2026-10.md 'Comparing the
  carry options'; ~/.claude/plans/eval-reduction.md Task G.

### D-20261029-N1 · Track E row G (G1) · Delegated: G2 runs with SmoothCostRank E3+E4 tau 0.25
- 2026-10-29, owner. Recorded by Newt Agent. Supersedes: none.
- Delegated (D-20261028-N2: the plan's default at the decision point). On Task E's other 24
  full-sample voxels the strict search is wrong on 1 voxel (1729, lenient) and none under the C++
  rule; E3+E4 tau 0.25 fixes 1729 (traced to the arm) and breaks none under either rule,
  E1+E2+E3+E4 fixes none. By the plan's rule (more fixed minus broken across both rules) G2
  carries SmoothCostRank E3+E4 tau 0.25, with the per-image-set SmoothCost cache on whenever the
  key is on. The choice rests on one voxel; the owner can override. Where:
  ~/.claude/plans/eval-reduction.md Task G1; MIGRATION_HISTORY Task G1;
  .claude/reports/20261028-2338-G1.md.

### D-20261029-N2 · Track E row G (G2) · Delegated: OutwardSearch stays off; the pilot baseline is today's key-off rerun
- 2026-10-29, owner. Recorded by Newt Agent. Supersedes: none.
- Delegated (owner 2026-10-28: Newt takes the plan's decision points). On the 2,000-voxel pilot
  region OutwardSearch (Task F fallback policy for BFS neighbours and reseeds) gave 36 wrong
  against 37 key off (18 fixed, 17 broken, 8 traced to it) at about 2.4x the evaluations; Task F's
  acceptance test (strict cost at most 0.8) accepts a neighbouring grain at boundaries, and all 9
  reseeds on its path ended 33-55 deg wrong. The key stays opt-in and off; a revisit would use F4
  alone for neighbours and no policy for reseeds. The stored 2026-10-10 pilot no longer reproduces
  on today's code (same 37 wrong, 14 fixed / 14 broken by voxel), so G compares against a key-off
  rerun. The overlap rule stayed lenient (no decision names one for Task G). Where:
  MIGRATION_HISTORY Task G2; .claude/reports/20261029-0212-G2.md.

### D-20261029-N3 · Track E row G (G3) · Delegated: G3 adds one BFS pilot arm with SmoothCostRank alone
- 2026-10-29, owner. Recorded by Newt Agent. Supersedes: none.
- Delegated (owner 2026-10-28: Newt takes the plan's decision points). The plan's G3 allows a new
  ablation run when two levers interact visibly: on the pilot region the integrated configuration
  (SmoothCostRank E3+E4 plus OutwardSearch) has 42 wrong and OutwardSearch alone 36, so
  SmoothCostRank may cost accuracy in BFS while it fixes 34 of 41 on E0. G3 runs one pilot arm
  with SmoothCostRank E3+E4 tau 0.25 (with the cache) and OutwardSearch off, against the key-off
  rerun, paired, breaks traced (D-20261014-N13), before the total is written. Where:
  ~/.claude/plans/eval-reduction.md Task G3.

### D-20261031-I1 · H12 · Delegated: make_brief supersede fix (Newt's request); a Supersedes field that begins with 'none' supersedes nothing
- 2026-10-31, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent, on Newt Agent's request of 2026-10-29 and the review of the fix (APPROVE, 3 fixes):
  make_brief.sh refuses an approval only when a later entry's own Supersedes field (below the
  header line) names the approval's ID as a whole ID (I1 is not I10). A field that begins with
  'none' supersedes nothing even when it cites an ID, for example 'Supersedes: none (a one-off
  exception to D-20261014-N1)': the entry says it replaces nothing. A withdrawal that names the
  plan still refuses as before. On the real DECISIONS.md the fix changes one brief:
  todo_batched_quick_mc.md, which the old check wrongly refused. Where: ROADMAP Track H row H12;
  scripts/dev/make_brief.sh; .claude/reports/20261031-review-brief-supersede.md.

### D-20261031-I2 · H12 · Delegated: correction to D-20261031-I1; no brief needed for the S3 plan, withdrawal check left as is
- 2026-10-31, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium Agent: D-20261031-I1 names the S3 plan's file and uses the word for taking a plan back, so
  make_brief.sh's withdrawal check, which this fix left unchanged, now reads I1 as taking back the
  S3 plan, and I1's sentence that the fix changes one brief no longer holds. S3 finished on
  2026-10-14 and needs no brief, so no action: the log stays append-only and the S3 plan gets a
  fresh owner entry only if it is ever reopened. That check is a keyword test on any later entry
  that names a plan; narrowing it is listed for the owner's review, not started. Lesson for
  entries: describe a plan by its row (S3) when the text also mentions taking plans back. Where:
  ROADMAP Track H row H12; scripts/dev/make_brief.sh; docs/owner_review_2026-10-30.md item 16.

### D-20261101-I1 · S3 · BatchedQuickMC: the E0 difference is probably fine; a large-scale test comes next, and its results must be statistically the same as key off
- 2026-11-01, owner. Recorded by Ilium Agent. Supersedes: none.
- Owner, on the review page (comment on S3's E0 result, 107 of 600 wrong against 111): 'This is
  probably fine. I want to do large scale test next. We expect minor difference with different
  stochasticity. We need answers to be statistically the same.' BatchedQuickMC stays off by
  default until that test (D-20261014-I2). The large-scale test is the next Track S item; Ilium
  Agent drafts its plan for the owner's approval, and the plan sets the scale and how
  'statistically the same' is judged (against the key-off arm's own seed-to-seed spread; no
  p-values). Where: ROADMAP Track S row S3; docs/owner_review_2026-10-30.md (review page, S3
  finding).

### D-20261101-I2 · S4 (later stage) · Ordering of the layer-wide screen and the BFS goes to a later autoresearch stage, with experiments
- 2026-11-01, owner. Recorded by Ilium Agent. Supersedes: none.
- Owner, on the review page (S4 thread, after Ilium's answer that the S4 hybrid keeps the BFS):
  add to a todo for a later autoresearch stage the question of which order of operations keeps the
  advantage of the layer-wide screen and of the BFS, for example screening in batches, growing the
  BFS from every voxel whose candidates take, and alternating back to screening; with experiments.
  Recorded as idea A1 in docs/todo_autoresearch.md (new file for such ideas). Not planned; waits
  for that stage and the owner's pick. Where: docs/todo_autoresearch.md A1; ROADMAP Track S row
  S7.

### D-20261109-N1 · Both tracks (review page) and Track A (autoresearch) · Owner delegates the 24 review items to Newt and sets five goals for an autoresearch stage that restarts by itself
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- Owner, 2026-11-09, to Newt Agent: 'go through the doc and make the best possible decision based
  on what you understand and answers from shmem.' Goals, in the owner's words: (1) 'Find ways to
  reduce both the total number of operations and the total wall-clock time of performing the HEDM
  reconstruction for a strain-free sample using the current physical assumptions.' (2) 'Look for
  ways to accelerate with deep learning methods', in particular an approximation that accelerates
  the grid search, given the repetitive nature of the problem. (3) 'We are allowed to discard this
  assumption of voxel-by-voxel reconstruction.' (4) A lower-quality approximation followed by a
  final optimization is fine. (5) 'General early stopping algorithms are desired': a way to
  predict whether further optimization will improve a voxel. Also: 'Create a mechanism so that you
  will restart automatically and save states accordingly', using shmem. Newt reads this as the go
  for an autoresearch stage (Track A) that runs without a per-plan approval: its plan is
  docs/autoresearch/PLAN.md and its saved state docs/autoresearch/STATE.md. Limits Newt keeps: no
  merge to develop without the owner's go, no change to .claude/ files, and no loosening of a
  guard that constrains the agents. Where: docs/autoresearch/PLAN.md; docs/autoresearch/STATE.md;
  docs/owner_review_2026-10-30.md, section "Decisions of 2026-11-09".

### D-20261109-N2 · Track E rows G, X1 (review items 1-8) · Delegated: review items 1-8 (Track E)
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- Newt Agent, under D-20261109-N1. Item 1: the G3 arm's 3.2 h run is accepted after the fact;
  D-20261029-N3 approved the arm and set no time. From now on a run that passes its plan's stop
  line is stopped, or gets a decision entry before it goes on. Item 2: the key-off rerun on
  today's code is the baseline of the BFS pilot. Item 3: the timing A/B is run as part of item 14,
  with FusedCost on in both arms, SmoothCostRank on against off, one worker each, on a quiet
  machine; no separate run with FusedCost off. Item 4: SmoothCostRank stays opt-in and off by
  default; no golden moves. Item 5: the full-sample BFS run goes ahead once, after item 14, with
  FusedCost on and SmoothCostRank off; it is the wall-clock and accuracy baseline of Track A and
  the key-off arm of the BatchedQuickMC large-scale test (D-20261101-I1). A SmoothCostRank arm at
  full scale is run only if item 14 shows a gain larger than the pilot's 3 voxels. Item 6:
  /review-task on the review-fix branch, then merge it. Item 7: Newt opens the PR to develop after
  item 6; the merge waits for the owner's go. Item 8: the ScreenEvaluator follow-up is folded into
  the screen port (item 21) and is not scheduled by itself. Where:
  docs/owner_review_2026-10-30.md, section "Decisions of 2026-11-09".

### D-20261109-N3 · Both tracks (review items 14-16) · Delegated: review items 14-16 (both tracks)
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- Newt Agent, under D-20261109-N1. Item 14: FusedCost and SmoothCostRank are run together on E0
  (200 voxels) and on the 2,000-voxel BFS pilot; Newt runs it; the arms compare wrong counts, not
  records, and the timing arms are single-worker on a quiet machine. Item 15: Ilium runs S6 on
  Ilium's machine, where the Phase D images now are, after the screen port (item 21), with
  FusedCost on and off. Item 16: the withdrawal check of make_brief.sh stays as it is; it is
  narrowed only when it refuses a plan that is in use. Where: docs/owner_review_2026-10-30.md,
  section "Decisions of 2026-11-09".

### D-20261109-N4 · Track S rows S1-S6, Q8, H12 (review items 9-13, 17-24) · Delegated: review items 9-13 and 17-24 (Track S); items 12, 22 and 23 stay with the owner
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- Newt Agent, under D-20261109-N1. Item 9: S2 is accepted on its record (the agent reviews, E0 35
  of 200 wrong with and without the key, agreement with C++ on Task A); the owner's own reading
  stays open and nothing waits for it. Item 10: the fused CPU path takes priority over a GPU
  SearchDevice for the strict cost; no GPU work on it until CUDA hardware is in use. Item 11:
  BatchedQuickMC is not part of the recommended configuration when FusedCost is on, because it
  gives no gain there; its large-scale test (D-20261101-I1) stays, and Newt recommends that its
  accuracy arms run with FusedCost on, which is cheaper and does not change the random-draw
  question; Ilium's plan sets the scale. Item 13: Ilium reruns the S1 and S3 timing scripts
  single-worker under the preflight and relabels. Item 17: FusedCost is validated on GCC and Linux
  before its first use there; no such machine is in use, so nothing runs now; Ilium makes the key
  warn on a compiler other than the validated one if it does not already. Item 18: Newt checks on
  Newt's machine, which has the C++ build, whether clang fuses the clip's multiply-adds. Item 19:
  the C++ toolchain is not installed on Ilium's machine; C++ checks go to Newt through the
  mailbox. Item 20: FusedCost stays opt-in and off by default; every timing baseline of new work
  uses FusedCost on. Item 21: porting the batched screen to the fused C path is the next Track S
  speed item; Ilium drafts the plan, designed layer-wide (one orientation across all voxels at a
  time, the owner's point of 2026-11-09 and idea A1); the plan's approval stays with the owner.
  Item 24: no reseed (owner rule of 2026-10-12: no extra seeds); the Phase D pilot on no-drop
  images is dropped as a separate item; studies stay on clean images; scripts with their own
  VARIANTS use the no-drop images at their next run and are not rerun for that alone. Not decided
  by Newt, because they change what the agents themselves may do: item 12 (two .claude/ edits;
  Newt recommends yes to both), item 22 (let a merge commit bring .claude/ files; Newt recommends
  yes when the files are identical to the merged branch's committed files) and item 23 (rerunning
  a refused command in the form the refusal message suggests; Newt recommends confirming). Where:
  docs/owner_review_2026-10-30.md, section "Decisions of 2026-11-09".

### D-20261109-N5 · Track E row AR (autoresearch) · Delegated: the autoresearch plan docs/autoresearch/PLAN.md is approved under the owner's mandate of 2026-11-09
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- Newt Agent, under D-20261109-N1 (the owner asked for an autoresearch stage that runs and
  restarts by itself): docs/autoresearch/PLAN.md is approved for implementation, steps AR0 to AR5,
  in the plan's order and with its gates (AR2, AR3 and AR4 start only after AR1). The owner has
  not read the plan; any word from the owner on it replaces this entry. Limits as in N1: no merge
  to develop without the owner's go, no change to .claude/ files, no loosening of a guard on the
  agents. Where: docs/autoresearch/PLAN.md; docs/autoresearch/STATE.md.

### D-20261109-N6 · Track E row G (G3, broken-voxel rule) · Delegated: a boundary crossing counts neither for nor against a lever under D-20261014-N13
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- Newt Agent, under D-20261109-N1, on the review of the review-fix branch
  (.claude/reports/20261109-review-review-fixes.md). A broken voxel whose trace is a boundary
  crossing (a boundary voxel that takes a touching grain's orientation) counts neither for nor
  against the lever under the broken-voxel rule D-20261014-N13; it is reported in its own column.
  Reason: the owner's comment of 2026-11-01 on the review page, on the finding that every wrong
  BFS voxel is a boundary voxel carrying a touching grain's orientation: 'this is acceptable and
  expected'. No default depends on this; SmoothCostRank stays opt-in (D-20261109-N2). Where:
  docs/findings_eval_reduction_2026-10.md (Terms and the N13 paragraph); MIGRATION_HISTORY Task
  G3.

### D-20261109-N7 · Track E row AR (AR1 review, AR2 arm) · Delegated: AR1's 3 broken voxels count against the best-1 arm; AR2 uses the finisher on the best 2, set before its first run; non-cubic operators fixed now
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- Newt Agent, under D-20261109-N1, on the review of AR1 (.claude/reports/20261109-review-ar1.md,
  'Decisions for the owner' 1-3). (1) The 3 pilot voxels that the arm 'strict cost on 5
  candidates, finisher on the best 1' breaks (7141, 8848, 23942) are traced to the arm's own rule:
  it picks one candidate before refining, and after the finisher the true grain has the lower
  strict cost in all 3. Under D-20261014-N13 a break caused by the arm counts against it, so they
  count against best-1. They are also shown in a boundary-crossing column (D-20261109-N6) if the
  grain taken touches the voxel; that column does not excuse them here. (2) AR2 uses the finisher
  on the best 2 candidates as its main arm, with best 1 and best 3 as side arms. Best 2 and best 3
  were both chosen after the best-1 result was seen, so their 0 of 69 is not a test. The test is a
  holdout: voxels outside the pilot region, picked by a fixed stride (no random draw), all three
  arms fixed before the run, wrong counts against the truth. It runs as the first step of AR2. (3)
  The plan says symmetry comes from the crystal's own operators, so ray_indexing.py is fixed now
  to use Cartesian operators, with a hexagonal unit test. The findings page still says the
  measurements are on one cubic copper sample. Reason for (2): at a boundary voxel two real grains
  are among the candidates, and the unrefined cost does not rank them reliably; refining two costs
  about 800 more evaluations per full search (mean 1,926 against 1,125), small against 332,000.
  Where: docs/autoresearch/PLAN.md (AR1, AR2); docs/autoresearch/findings_ar1_ray_indexing.md.

### D-20261109-N8 · AR / AR2 step 1 · Holdout review: voxel 17361 counts against all arms; fix the angle function in the library; split boundary voxels by edge and corner
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- Taken by Newt under the owner mandate D-20261109-N1, on the three questions of the AR2 step 1
  review. (1) Voxel 17361 counts against all three finisher arms. The review found the strict cost
  at the exact truth is 0.0949, lower than the touching grain's 0.133; the finisher started 0.044
  degrees from the truth and did not move (0.162 before and after). So the cause is the path's own
  finisher, not the cost. By the rule of D-20261014-N13 and item 1 of D-20261109-N7 a break traced
  to the arm counts against it. It is still shown in the boundary-crossing column (the voxel took
  a touching grain's orientation), which does not excuse it. This holds if the fix run reproduces
  the truth cost; if not, the question is open again. The best-2 choice for the key is unchanged
  (1 of 300 in best 2 and best 3). The finisher failure inside 0.05 degrees becomes a named test
  case for AR2 step 2 and AR5. (2) misorientation_deg in icenine/ray_indexing.py is fixed now, not
  left as a TODO: a robust angle (atan2 of the sine and cosine parts, or an orthonormalised input)
  with a unit test for small angles and slightly non-orthonormal input. The AR2 errors are
  regenerated. The AR1 page gets a note that its small distances have a floor of a few thousandths
  of a degree and that its wrong counts are not affected; AR1 medians are recomputed only if the
  stored files allow it without a rerun. (3) From now on a report with a boundary set gives the
  split into edge-sharing and corner-only voxels. Where:
  icenine_py/docs/autoresearch/findings_ar2_holdout.md; review report
  .claude/reports/20261109-review-ar2a.md; PLAN.md AR2.

### D-20261109-N9 · AR / AR2 step 2 · AR2 step 2 review: route-caused boundary crossings are not the key's own path, on both sides; the contended timing pair is a first reading only
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- Taken by Newt Agent under the owner's mandate D-20261109-N1, after the review of AR2 step 2
  (review report 20261109-review-ar2b, 14 required fixes, all to be applied). (1) The key's own
  path, in the sense of N7 and N8, is the chain on the seed voxel itself: index, strict cost on
  the first 5 candidates, finisher on the best 2, seed test. A boundary voxel that takes a
  touching grain's orientation through local fits from a seed that is right for its own voxel is a
  boundary crossing under N6. It counts neither for nor against the key. This holds on both sides:
  the 7 (12) broken voxels and the fixed voxels are treated alike, and the fixed voxels are traced
  in the key-off run the same way. The key does change which voxels are seeds, so it changes the
  route of the BFS; that is an effect of the key on the run, it is reported in full (whole-run
  wrong counts and the paired split for both pairs), and it is named as not tested. The claim on
  accuracy is therefore limited to: no wrong result among the 64 indexed seeds, and a whole-run
  wrong count that is not shown different from the grid search, given that two key-off runs that
  differ only in FusedCost differ by 27 voxels. No gain in wrong voxels is claimed. (2) The
  wall-time pair (213 s against 1633 s, FusedCost on) was taken after a refused preflight on a
  contended machine, one run per arm. It is kept as a first contended reading. Review item 14
  (both keys on E0 and the pilot, timing on a quiet machine) stays open and is not replaced by it.
  The owner can overrule either point. Where: icenine_py/docs/autoresearch/findings_ar2_key.md;
  docs/autoresearch/STATE.md.

### D-20261109-N10 · AR / AR2 step 2 · AR2 step 2: a voxel that is a seed in both runs and that only the indexed search gets right counts for the key (the mirror of N7 and N8)
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- Taken by Newt Agent under the owner's mandate D-20261109-N1. It fills a gap in D-20261109-N9
  that the fix pass found. In the stored pair two fixed voxels (5876 and 24334), and in the
  FusedCost-on pair one (16476), are seeds in both runs: the grid search's result was rejected by
  the seed test and the indexed result is right. These fixes are on the key's own path as N9
  defines it (the chain on the seed voxel itself). N7 and N8 count a break on the lever's own path
  against it, so a fix on the lever's own path counts for it. They are reported as counts, with
  the voxel numbers, beside the count of breaks on the key's own path (0). They are not turned
  into a rate and they do not change the limit set in N9: the whole-run wrong count is not shown
  different from the grid search, and every other fixed or broken voxel is a boundary crossing
  that counts neither way. The owner can overrule this. Where:
  icenine_py/docs/autoresearch/findings_ar2_key.md.

### D-20261109-N11 · Track E (PR #41, review items 12, 22, 23) · Owner: merge PR #41; yes to review items 12, 22 and 23
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- The owner's word in the session of 2026-11-09 ('I have left comment in the artifact. Go ahead
  and execute according to the comments'), on the owner's own comments on the review page. (1) PR
  #41 (feature/eval-reduction to develop) is merged with a merge commit (357a2b3, at branch head
  16c72db); the branch is kept, later autoresearch steps go in a new PR. (2) Item 12: yes to both
  edits under .claude/ (show the 'note:' lines of make_brief.sh in /implement-plan and
  /apply-fixes; drop the extra -q in run-tests-on-change.sh). Newt applies them to the working
  tree; the commit of .claude/ files stays with the owner. (3) Item 22: yes. precommit_check.py
  lets a merge commit bring a .claude/ or CLAUDE.md file only when the staged file is identical to
  that file on the branch being merged; every other .claude/ change is still refused. (4) Item 23:
  yes. When the worktree guard refuses a combined shell command, an agent may rerun it as plain
  commands, as the refusal message suggests. Doing a refused action another way (such as writing a
  file through the shell after Edit or Write was refused) stays forbidden. Also noted: the shmem
  embedder is no longer run, so its warning is expected. A page comment alone was not acted on;
  the session message is the authority. Where: docs/owner_review_2026-10-30.md items 12, 22, 23;
  PR #41; scripts/dev/precommit_check.py; ROADMAP Newt mailbox.

### D-20261109-N12 · Shared tooling (pre-commit check) · Owner: tighten the pre-commit check on CLAUDE.md/.claude files
- 2026-11-09, owner. Recorded by Newt Agent. Supersedes: none.
- The owner's word in the session of 2026-11-09 ('tighten the pre-commit check'), on the gaps the
  review of the item 22 change found. (1) A new hook .githooks/pre-merge-commit runs the same
  check for a merge that git commits by itself. git has not recorded the merged branch at that
  point, so a merge that brings CLAUDE.md/.claude files is left staged and is concluded with 'git
  commit', where the item 22 exception applies. (2) A staged deletion, a rename away and a change
  of type of a guarded file are refused like any other change. (3) Guarded names are matched in
  any letter case. (4) The item 22 exception (D-20261109-N11) is narrowed so that closing (1)
  holds: the merged commit must be on a remote-tracking branch (a stash or another local-only
  commit gets no exception), and a deletion passes only when the merged branch lacks the file and
  the merge base has it. The override ALLOW_CLAUDE_CONFIG=1 stays with the owner. Where:
  scripts/dev/precommit_check.py; .githooks/pre-merge-commit; icenine_py/README.md, Developer
  tooling.

### D-20261110-N1 · Track E / AR3 · Delegated: the AR3 plan (paint grains on the pilot, then the key and the whole layer) is approved under the owner's mandate
- 2026-11-10, owner. Recorded by Newt Agent. Supersedes: none.
- Taken by Newt Agent under the owner's mandate D-20261109-N1; not the owner's own word. The
  section 'AR3: index the layer and paint grains' of docs/autoresearch/PLAN.md is approved for
  implementation. Step 1 is a study on the 2,000-voxel pilot region: a new module
  icenine/grain_paint.py lists grains by indexing unpainted voxels one at a time, scores each
  grain's orientation on the voxels around it with one strict evaluation and no fit, and gives
  each voxel the passing grain with the lowest strict cost. It is compared with the stored key-on
  run (22 wrong, 2,819,906 evaluations) and the stored key-off run (37 wrong, 25,547,935
  evaluations). The gate for step 2 is fixed before the run: paint breaks at most 20 voxels
  against the key-on run and needs at most half of its evaluations. A wrong-count difference no
  larger than that between two runs of the same arm (20 one way, 7 the other) is not called a gain
  or a loss. A margin rule for contested voxels is only picked in step 1 and is tested in step 2
  on the voxels outside the pilot region. Step 2 (an opt-in key PaintGrains, off by default, and
  the whole layer) has its details fixed after step 1 and its review. The layer-wide vote over all
  rays of the first draft is not built. Where: docs/autoresearch/PLAN.md, section AR3;
  docs/autoresearch/STATE.md.

### D-20261110-N2 · Track A / AR3 step 2 · AR3 step 2 plan approved under the owner's mandate (key PaintGrains, whole layer)
- 2026-11-10, owner. Recorded by Newt Agent. Supersedes: none.
- Approved by Newt Agent under D-20261109-N1 (delegated): the details of AR3 step 2 in PLAN.md.
  Opt-in key PaintGrains (default off, needs SeedByIndexing 1, loop of grain_paint.py unchanged).
  Checks first: the pilot through the key gives step 1 again; the key-on BFS with the fit log is
  the stored run. Then the whole layer (24,570 voxels), paint against the key-on BFS, counts split
  pilot / rest. Fixed here: a fit-settled voxel keeps its fit in the key, both rules are reported;
  a paired count is called a gain only when the error-rate intervals do not overlap; no extra
  seeds. The fit log and the seed record are readings for AR5, not tests. Budget 3 hours. Task
  ar3b. No merge to develop without the owner's word. Where: icenine_py/docs/autoresearch/PLAN.md,
  section AR3, Step 2, the details.

### D-20261110-N3 · Track E / PR #42 · PR #42 merged into develop; painting is a starting approximation, grains hold orientation gradients
- 2026-11-10, owner. Recorded by Newt Agent. Supersedes: none.
- Owner, 2026-11-10: comment 'merge now' on the PR #42 card of the review page, then 'added
  comments' in the session. PR #42 merged with a merge commit (66e93fd, at branch head 0b3852d);
  branch feature/eval-reduction kept. Two more owner comments, same day: (1) 'one cost evaluation
  and no fit per voxel is fine as a starting approximation, and then we can refine it'; (2) it is
  NOT true that a strain-free grain has one orientation: there can be micro-gradients, so every
  voxel in a grain may have a different orientation. Consequence: 'one orientation per grain' is a
  property of the simulated test sample only; painting gives the start, a per-voxel refinement
  follows; the next step needs a test sample with gradients inside grains. Also in the session:
  the agent may run more per wake (more token limit left). Where: PR #42;
  icenine_py/docs/autoresearch/STATE.md; review page comments of 2026-11-10.

### D-20261110-N4 · Track A / AR3 step 3 · AR3 step 3 plan approved under the owner's mandate (paint as a start, refine, on grains with gradients); step 2 review items settled
- 2026-11-10, owner. Recorded by Newt Agent. Supersedes: none.
- Taken by Newt Agent under the owner's mandate D-20261109-N1; the owner may overrule any part.
  (1) AR3 step 3 (task ar3c) is approved as written in PLAN.md at 93be764: a gradient sample (same
  voxels and grains, a linear rotation field inside each grain, levels 0, 0.1, 0.2, 0.5 and 1.0
  degrees drawn per grain with seed 0), clean images of it, runs P (paint), R (one local fit from
  each painted start, with a trace), R2 (a second fit on every fourth voxel) and D (BFS with
  indexed seeds), and stop-rule readings for goal 5 picked on even-id grains and reported on
  odd-id grains. The gradient model is the agent's choice, not measured data; the owner is asked
  which sizes and shapes are typical, and the task does not wait for the answer. It starts after
  ar3b is merged. (2) Step 2 review (verdict: request changes, 12 fixes): all 12 go to a fresh
  Opus implementer, plus fix 13, a dated correction note on step 1's findings page. (3) The
  reviewer's three owner items: the plan's premise line is already corrected (D-20261110-N3); the
  whole-layer result files of step 2 stay in git as committed and history is not rewritten, and
  from step 3 on a task commits only the result files its summary reads; the off-path claim cites
  check B, with no stored output from the base commit. Where:
  icenine_py/docs/autoresearch/PLAN.md, AR3 Step 3; review report 20261110-0900-review-ar3b.md
  under .claude/reports.

### D-20261110-N5 · Track A / AR3 step 4 · AR3 step 4 plan approved under the mandate: a local fit whose step matches the error (task ar3d)
- 2026-11-10, owner. Recorded by Newt Agent. Supersedes: none.
- Approved under the mandate D-20261109-N1 by Newt Agent; the owner can overrule it. Step 3's
  review found that the BFS local fit moves almost no painted voxel because its step is too
  coarse: a box of 5 degrees, smallest trial step about 0.7 degrees, against a start that is 0.04
  degrees off (median), while the true orientation costs less than the start on 16,591 of 19,013
  painted voxels. Step 4 (PLAN.md 'Step 4', commit 592937a) tests a fit that starts with a small
  box. Part 1: boxes of 5/8, 5/16 and 5/32 degree, a ladder and the existing CMA fit, from step
  3's stored painted starts; picked on set A, reported on set B; stop if no arm is a gain over the
  start. Part 2: the stop readings again on the picked fit. Part 3: whole-layer runs with a script
  wrapper 'fine first, stock as fallback' (BFS with indexed seeds and paint plus refine on the
  gradient sample, and the BFS on step 2's sample), each with a coded time cap. The box is set by
  the script, not by MinLocalResolution. No library change, no default changed, C++ not run. An
  opt-in key comes after, only if the whole-layer runs show a gain. Where:
  icenine_py/docs/autoresearch/PLAN.md, AR3 'Step 4'.

### D-20261110-N6 · Track A / AR3 step 5 · AR3 step 4 review decisions and step 5 plan approved under the mandate: let neighbouring grains compete for the boundary voxels (task ar3e)
- 2026-11-10, owner. Recorded by Newt Agent. Supersedes: none.
- Approved under the mandate D-20261109-N1 by Newt Agent; the owner can overrule it. (1) Step 4's
  review (verdict request changes, 14 required fixes, text and summary output only) asked two
  questions. Decided: no rerun of step 4 and no wrapper that reads the caller's frame; the plan's
  fallback on the caller's 0.9 test stays untested and the page says so. Reason: the review's
  recomputation from the stored records shows the extra wrong voxels sit in fits the caller
  accepts (paint run: all 107 fit-settled newly wrong voxels passed the 0.9 test, 104 from a grain
  seeded in another true grain; step 2's sample: 114 accepted at the first neighbour fit, 56 at a
  revisit), so a fallback for rejected voxels would not act on them. (2) Step 5 (PLAN.md 'Step 5',
  commit 05577ac): after a run, each voxel with a neighbour whose orientation differs by more than
  0.5 degrees gets one CMA fit (0.05 degrees, 61 evaluations) from that neighbour's orientation
  and keeps the better end. Part 0 reads from stored maps and single scores whether the pass can
  work (a gate). Part 1 runs the pass on six stored maps with three keep rules (cost, hit ratio,
  confidence), picked on set A, reported on set B, with a skip rule for when a try is needed. Part
  2 makes two whole-layer runs with the existing library keys (LocalOptimizer cma, CMASigma0,
  CMANeighborMaxEvals, CMARetrySigma0) and no wrapper, each with a coded 60-minute cap. No library
  change, no default changed, C++ source read but not run. An opt-in key for the pass comes after,
  only if the pass is a gain. Where: icenine_py/docs/autoresearch/PLAN.md, 'Step 5' of the AR3
  section.

### D-20261110-N7 · Track A / AR3 step 5 review, step 6 · AR3 step 5 review decisions approved under the mandate: all 15 fixes with no rerun, per-try arrays kept, no K rerun yet, step 6 is the seed cost
- 2026-11-10, owner. Recorded by Newt Agent. Supersedes: none.
- Approved under the mandate D-20261109-N1 by Newt Agent; the owner can overrule it. Step 5's
  review (verdict request changes, 15 required fixes, the numbers reproduce) left four choices.
  (1) All 15 fixes are applied, to text, the summary script and the tables only; no pass and no
  run is repeated. Added to the fixes: a plot first (wrong voxels against evaluations per
  pipeline, pass on both sides), the old-end score check cited, 'fit-settled' and 'grid fallback'
  defined, a reading line before the ratio table, wall times labelled as contended, the run order
  of boundary_summary.py and boundary_ratio.py documented, and 'local optimization' in place of
  'finisher' in the page and the step 5 scripts. The lead result is the one-footing comparison:
  with the boundary pass on both sides the CMA-first run shows no detectable difference in wrong
  voxels on the gradient sample at about a third to 58% of the evaluations (D4+pass 74 against
  D+pass 88 wrong at 38.8M against 122.1M; P4R4+pass 82 against PR+pass 78 at 54.6M against
  94.7M), and 0 against 14 on step 2's sample at 16.6M against 36.1M. (2) The 14 tries.npz files
  (about 47 MB) stay as committed: D-20261110-N4 allows them because the summary reads every one,
  and removing them now would not shrink the history. From the next task on per-try arrays are
  stored slim (only the fields the summary reads). The owner can ask for a prune or Git LFS. (3)
  Run K (cma keys with the paint) is not repeated now: at its 60-minute cap it had already used
  38.1M evaluations, more than D4's 35.1M, and 31.07M of them were grid fallbacks behind failed
  index seeds; it also started under 8-process contention. It is repeated after step 6. (4) Step 6
  is the seed cost: grid fallbacks behind failed index seeds are the largest remaining cost (25.8M
  of D4's 35.1M). First read the stored records (does the grid search succeed where the index seed
  failed, and why does the index miss), then test remedies. The opt-in library key for the
  boundary pass waits for a check on a sample with larger grains, because here 93 to 94% of the
  voxels are boundary voxels. The review's C++ note (in parallel mode several clients fit accepted
  voxels again and the server keeps the higher confidence, the keep rule that breaks thousands of
  right voxels here) goes to the owner as a question; it is not tested. Where:
  icenine_py/docs/autoresearch/PLAN.md, AR3 section; review report
  .claude/reports/20261110-2220-review-ar3e.md; findings_ar3_boundary.md after the fixes.

### D-20261110-N8 · Track A / AR3 step 6 · AR3 step 6 plan approved under the mandate: why the index seed fails and what replaces the grid fallback (task ar3f)
- 2026-11-10, owner. Recorded by Newt Agent. Supersedes: none.
- Approved under the mandate D-20261109-N1 by Newt Agent; the owner can overrule it. Step 6
  (PLAN.md 'Step 6', commit c88fdd8) follows D-20261110-N7. Reason: in run D4 the grid search
  behind 76 failed index seeds took 25.8M of 35.1M evaluations, about 340,000 each; 74 of the 76
  are at reseeds, 36 of them for a cluster of one voxel; and after the grid search 73 of the 76
  voxels have a hit ratio of 0.88 or more, so a passing orientation exists that the index seed did
  not reach (main session's reading of the stored records, to be derived again by the task). Part
  0 replays the seed step on the fallback voxels and on 200 control seeds and puts each failure in
  one class (the truth fails the seed test; no candidate near the truth; candidate not costed; not
  refined; local optimization ends away). Part 1 tests six remedies by replay, with no whole-layer
  run: more candidates, the cma-mode seed step, a wide CMA fit, tries from the neighbours, no seed
  for clusters of one or two voxels, the existing outward search. Picked on set A, reported on B,
  counts with intervals. Part 2 makes two whole-layer runs with the pick through a script wrapper
  (S, and K2 which repeats step 5's unfinished run K), each with a coded 60-minute cap, then the
  boundary pass; it is skipped if no remedy removes half of the grid fallbacks on A. No library
  change, no default changed. Committed results stay under 5 MB and hold only what the summary
  reads. Where: icenine_py/docs/autoresearch/PLAN.md, 'Step 6' of the AR3 section.

### D-20261111-N1 · Track A / AR3 step 6 review · AR3 step 6 review decisions approved under the mandate: all 11 fixes with no rerun, additions 7a, 7b and S5 kept as marked, next key direction left to the owner with a default
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Approved under the mandate D-20261109-N1 by Newt Agent; the owner can overrule it. Step 6's
  review (verdict request changes, 11 required fixes, the code, the evaluation counts and the
  truth checks are sound) left three choices. (1) All 11 fixes are applied, to the findings page,
  the summary script, its tables and the plot label only; no replay and no whole-layer run is
  repeated. The summary script may be run again on the stored records; existing values in the
  summary files stay unchanged and new keys (the intervals) are added. (2) The remedies 7a and 7b
  and the whole-layer run S5 were added by the implementer beyond the plan after it saw the
  results. They are kept, marked 'added after seeing results', and are not part of the plan's
  pick. (3) The plan's pick (tries from the neighbours, then refine the best 5 of 20 candidates)
  cut the grid fallbacks from 76 to 11 and the evaluations from 38.8M to 16.1M (run S against D4,
  both with the boundary pass), but left 100 against 74 wrong voxels (intervals overlap), and the
  voxels settled by a neighbour's try were wrong on 59 against 27 of 602 before the pass (40 of
  the 59 still wrong after it). Run S5 (refine the best 5 of 20, no neighbour tries) gave 15
  fallbacks, 20.6M evaluations and 62 wrong; against S that is 62 against 100 with 57 fixed and 19
  broken; against D4 no detectable difference in wrong voxels at about half the evaluations. The
  page presents the choice of the next opt-in key as the owner's. The default taken until the
  owner says otherwise: the next work uses the seed step of S5 (best 5 of 20, no neighbour tries),
  checks it on a sample with larger grains before any opt-in key is added, and does not add the
  neighbour tries. One sample, one run each; no extra seeds. Where:
  icenine_py/docs/autoresearch/findings_ar3_seed.md (task branch feature/eval-reduction-ar3f);
  review report .claude/reports/20261111-0330-review-ar3f.md.

### D-20261111-N2 · Track A / AR3 step 7 · AR3 step 7 plan approved under the mandate: the same pipelines on a sample with larger grains (task ar3g)
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Approved under the mandate D-20261109-N1 by Newt Agent; the owner can overrule it. The plan
  (PLAN.md 'Step 7', commit 96b6902) follows D-20261111-N1. Reason: in step 6's run S5 (20.6M
  evaluations with the boundary pass, 24,570 voxels) the seed step took 6.0M and 15 grid fallbacks
  5.1M; both happen once per reseed, and the gradient sample has 497 grains of about 49 voxels.
  The step builds one sample 'big': the 32 orientations of random_cubic_32grains.mic on the
  gradient sample's voxel grid (about 770 voxels per grain), gradients by the step 3 function with
  seed 0 and spreads of 0, 0.2, 0.5, 1, 2 and 4 degrees per grain (0.4 to 4 degrees gives the same
  voxel-to-voxel change as 0.1 to 1 degree on the gradient sample), clean images. Four whole-layer
  runs with coded caps, each followed by the boundary pass: D (stock BFS with grid search, 180
  min), D4 (step 5), S5 (step 6's seed step, best 5 of 20), P4R4 (paint then CMA local
  optimization), 90 min each. Three predictions are written in the plan before the runs. The
  result decides which seed step the opt-in keys step uses, and names the next cost target from
  the split by stage. No library change, no default change, no extra seeds; one sample, one run
  each; the levels are a choice of the agent and the real gradient shape stays an open question
  for the owner. Where: icenine_py/docs/autoresearch/PLAN.md, section 'Step 7: the same pipelines
  on larger grains'.

### D-20261111-N3 · Track A / AR3 step 7 · AR3 step 7 review: apply all 10 required fixes, no rerun
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- The review of task ar3g (report .claude/reports/20261111-0500-review-ar3g.md) asks for changes:
  the data and code are sound, the findings page has wording and number errors. Apply all 10
  required fixes in one fixes round. No whole-layer run is repeated; existing values in the
  summary files stay unchanged; nothing under icenine_py/icenine/ changes. Any number a fix adds
  is recomputed from the stored records, not copied from the review. The choice of the seed step
  for the opt-in keys stays with the owner: on the sample with larger grains the seed step of run
  S5 costs more than that of run D4, so the plan's rule does not pick one. Taken under the mandate
  D-20261109-N1; the owner can overrule. Where: icenine_py/docs/autoresearch/PLAN.md, Step 7;
  findings_ar3_big.md.

### D-20261111-N4 · Track A / AR3 step 8 · AR3 step 8 plan approved: where the seed step's evaluations go in grains with strong gradients, remedies on stored cases first
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Approved under the mandate D-20261109-N1. Step 7 left the seed step as the largest stage on the
  sample with larger grains (57.5% of run S5), with 267 of 349 seed steps in the 6 grains at 4
  degrees and one seed step costing about 4.6 times what it costs on the first gradient sample.
  Step 8 (task ar3h, scripts and results only, no key, no default change) measures where those
  evaluations go (the split of one seed step; why the BFS rejected the fits that led to reseeds),
  tries three listed remedies on stored cases with a pick set and a test set split by grain, and
  spends whole-layer runs only on a remedy that cuts its target by a third and breaks no more
  cases than it fixes. Three predictions are written in the plan before the runs. The choice of
  the seed step for the opt-in key stays the owner's. The opt-in keys step waits for this result.
  Where: icenine_py/docs/autoresearch/PLAN.md, 'Step 8: why the seed step costs so much in grains
  with strong gradients'.

### D-20261111-N5 · Track A / AR3 step 8 · AR3 step 8 review: apply all 10 required fixes, no rerun; the plan's remedy E1 did not pass, run S8 is a variant chosen after results
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- The review of task ar3h (report .claude/reports/20261111-0545-review-ar3h.md) asks for changes:
  every number it recomputed matches, the page's wording and what it discloses do not. Apply all
  10 required fixes in one fixes round. No whole-layer run is repeated; existing values in the
  result files stay unchanged; nothing under icenine_py/icenine/ changes. On the reviewer's
  question: by the plan's own rule the listed remedy E1 (stop at the seed test, hit ratio 0.80)
  failed Part 2, so no listed remedy passed. Run S8 uses a variant (stop at hit ratio 0.95) whose
  test was added after the 0.80 rows of one test set had been read. Its whole-layer numbers stand
  as measured, and they meet the three parts of prediction 3, but the page reports prediction 3 as
  met by that variant only, not as a confirmed prediction of the plan. The variant is not
  confirmed until it is an arm of a whole-layer run on a sample it was not chosen on; the next
  step that makes such a run includes it. The hit ratio at the true orientation of the class A
  voxels (the reviewer's cheap check) is the first part of the next step, not of the fixes round.
  The choice of the seed step for the opt-in key stays with the owner; the key takes the stop test
  as a setting that is off by default. Taken under the mandate D-20261109-N1; the owner can
  overrule. Where: icenine_py/docs/autoresearch/PLAN.md, Step 8; findings_ar3_reseed.md.

### D-20261111-N6 · Track A / AR3 step 9 · AR3 step 9 plan approved: a start test for the neighbour fits, the rejected fits that end right, and a third sample that confirms or refutes run S8
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Step 9 (task ar3i) follows from step 8 and its review. Part 1 reads the stored runs: the hit
  ratio and error at the start and end of every neighbour fit and revisit, and for class A the hit
  ratio at the true orientation (the reviewer's check). Part 2 tries two remedies on stored cases
  with step 8's pick and test split: F1, a start test that skips the fit (threshold 0.90, 0.95 or
  0.98), and F2, one more fit from a rejected fit's end (first step 0.05, 0.1 or 0.2 degrees). The
  settings and the pick rule are fixed in the plan; nothing is added after a result, and the pick
  is committed before the test set is read. Part 3 builds a third sample, mid (about 120 grains,
  step 7's generator and seed value), and runs D4, S5 and S8 with the boundary pass on it: the
  confirmation that D-20261111-N5 asks for. Run S9 (S8 plus passed remedies) runs on all three
  samples only if a remedy passes. Four predictions are written in the plan. No library change, no
  key, no default change. The opt-in keys wait for this step, since the settings they would hold
  are still moving. Taken under the mandate D-20261109-N1; the owner can overrule. Where:
  icenine_py/docs/autoresearch/PLAN.md, Step 9.

### D-20261111-N7 · Track A / AR3 step 9 · AR3 step 9 review: apply all 13 required fixes, no rerun; S8 counts as confirmed on one new sample; the start test F1 at 0.90 is not a candidate for a key; F2's fail stands
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Review of task ar3i (report .claude/reports/20261111-0810-review-ar3i.md): changes requested, 13
  required fixes, all to wording and to what the page discloses. All 13 are applied; no
  whole-layer run is repeated; fixes 4 and 7 add numbers to start_summary.py and regenerate the
  summary from stored files. The reviewer's three open questions, decided under D-20261109-N1 (the
  owner can overrule each): (1) The sample mid counts as the confirmation that D-20261111-N5 asked
  for, worded as one run per arm on one new sample whose grain map (a 120-cell Voronoi map)
  departs from the plan. No setting was chosen on mid; the reviewer checked that both pilots ran
  on big and that mid's first result came after the pick commit. S8's seed step is therefore the
  recommended seed step for the opt-in key; what the key recommends stays the owner's choice. (2)
  The start test F1 at 0.90 passed on stored cases but in the whole layer run S9 broke more voxels
  than it fixed on all three samples (10 fixed and 15 broken on big, 31 and 38 on the first
  gradient sample, 9 and 18 on mid) and raised the median error. It is not a candidate for a key.
  F1 at 0.95 was a listed setting, kept no wrong start on stored cases, and was not run on a whole
  layer; a whole-layer test of it would be a variant chosen after results and needs its own plan.
  (3) F2's verdict (fails) stands: the pooled rule was written in the scripts before any F2
  result, and with F2's own cost set to zero the pooled cut is still under a third (21.5% and
  26.2%). The page gives F2 by sample (big alone 34.1%) and says that the pooled rule decides.
  Where: icenine_py/docs/autoresearch/findings_ar3_start.md;
  .claude/reports/20261111-0810-review-ar3i.md.

### D-20261111-N8 · Track A / AR3 step 10 · AR3 step 10 plan approved: opt-in keys in the library for the indexed seed step, its early stop, the CMA-first neighbour fit and the boundary pass; all off by default, no recommendation
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Plan in icenine_py/docs/autoresearch/PLAN.md, section 'Step 10: opt-in keys for the indexed seed
  step and the boundary pass' (task ar3j), approved under D-20261109-N1; it starts after step 9 is
  merged. First AR3 step that changes the library. Four keys, all off by default, off path
  bit-identical: BFSNeighborCMAFirst (step 4's neighbour fit: CMA first, the stock fit if the cost
  does not fall), SeedIndexCandidates (S5 uses 20 5), SeedIndexStopHitRatio (S8 uses 0.95),
  BoundaryPass. Tests are written first. Validation: the library runs with the keys must reproduce
  the stored script runs D4, S5 and S8 on the sample big exactly before the pass (evaluations,
  seeds, reseeds, wrong counts); a difference stops the step. Decided with the plan, each the
  owner's to overrule: (1) the third key BFSNeighborCMAFirst is in scope, because without it no
  stored run can be reproduced; (2) icenine/boundary_pass.py is added to the plan's list of new
  modules; (3) the pass is exact on the stored map and a prediction in the run, the library does
  not round its map; (4) the pass uses its own fixed generator; (5) never-fitted voxels are
  skipped by the pass; (6) the fixes round of this step may change the library files the step
  touched, with its tests and bit-identity checks run again (an exception to the runbook's fixes
  rule, this step only). The step recommends no seed step: the choice stays the owner's. No key
  for the start test (D-20261111-N7). No default flip, no golden change, no new dependency. Where:
  icenine_py/docs/autoresearch/PLAN.md (Step 10).

### D-20261111-N9 · Track A / AR3 step 10 · Review of AR3 step 10 (task ar3j): apply the 8 required fixes, test the rounded-map reading, keep the word fit on the repo page
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Under the owner mandate D-20261109-N1, yours to overrule. The review (Opus,
  .claude/reports/20261111-0941-review-ar3j.md) asks for changes in documentation and one test; it
  found the library code, the tests and the validation runs sound. (1) All 8 required fixes are
  applied in one fixes round (Opus), and the cheap suggestions with them (the both-wrong column,
  the uncounted evaluations given as numbers, the type annotations, the README rows, the
  docstrings, the note on the stored S8 label). (2) The reviewer left open whether to test the
  reading that the boundary pass inside the run differs from the stored run because its input map
  is not rounded. It is tested: the library pass is applied in a fresh process to the unrounded
  map of the library run for D4L and S8L if the cache holds that map bit for bit, otherwise by one
  BFS rerun of S8L. Equal tries to the in-run pass means the rounding accounts for the difference
  and the process does not; anything else is reported as it is. No setting is changed and no key
  default changes. (3) The page keeps the word fit, like the other repo pages; a repo-wide rename
  to local optimization waits for the owner's word. The review page for the owner says local
  optimization. No setting is recommended; the choice of the seed step stays the owner's. Where:
  icenine_py/docs/autoresearch/PLAN.md, Step 10; review report 20261111-0941-review-ar3j.md.

### D-20261111-N10 · Track A / AR5 step 1 · AR5 step 1 plan approved: stop rules for the boundary pass, the grid fallbacks and the finished map; scripts only
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Under the owner mandate D-20261109-N1, yours to overrule. Plan in
  icenine_py/docs/autoresearch/PLAN.md, section AR5 step 1 (commit a633bfd), task ar5a. No library
  file, key or default changes. Four parts: (1) a test before each try of the boundary pass, one
  evaluation of the try's start, picked on big and the first gradient sample and tested on mid;
  (2) a ladder of cheaper searches before a grid fallback (100 candidates costed; grid level 0
  then a local optimization), on the 31 stored fallback cases of run S8; (3) a flag for voxels
  that are not finished: the hit ratio after the pass, the vote count of unexplained rays (H3),
  and a check that reruns the seed step at flagged voxels; (4) one whole-layer run T1 per sample
  with the rules that passed. Pass rules and six predictions are written in the plan before any
  run. The original AR5 line (AR0 thresholds with the H3 vote count) is kept as part 3; parts 1
  and 2 are added because the stage table of S8 shows the boundary pass at 12 to 20 percent and
  the grid fallbacks at 19 and 30 percent of the evaluations. It starts while the fixes round of
  AR3 step 10 runs: it does not depend on step 10. Where: icenine_py/docs/autoresearch/PLAN.md,
  AR5 step 1.

### D-20261111-N11 · Track A / AR5 step 1 (task ar5a) · AR5 step 1 review: apply the 12 required fixes; V1 and the vote count are not adopted
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Decided by Newt Agent under the owner mandate D-20261109-N1; the owner may overrule. Review
  report .claude/reports/20261111-1215-review-ar5a.md: verdict REQUEST CHANGES, documents only, no
  result changes, the negative results stand at their stated scope. (1) All 12 required fixes are
  applied in a fixes round (Opus), with the parent merged into the task branch first (fix 2). The
  suggestions that cost little are applied too: remove the unused held parameter, offset the plot
  labels, draw the vote cut, give the share of boundary voxels, Wilson intervals for the Part 2
  counts. The audit_numbers tool artifact goes to the tooling queue. (2) The implementer's two
  questions: the check V1 after the 0.95 flag is not adopted (13 of 31 wrong voxels fixed on mid
  for 105.4% of the run's evaluations; 99.4% of them spent on right voxels); the vote count is
  dropped as a stop rule (its cut was chosen after the result and moves between samples). (3) The
  run T1 is not run. No new point on the plots of cost against accuracy. (4) The reviewer's three
  next-step options: option 1 (try test P2 m=2 as an opt-in key) is left to the owner; options 2
  and 3 (a second check beside the hit ratio) go into AR4 step 1. Where: PLAN.md AR5 step 1;
  findings_ar5_stop.md.

### D-20261111-N12 · Track A / AR5 step 2 (task ar5b) · AR5 step 2 plan approved: seed refinement by CMA and no MC fallback after CMA, replay then run S10
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Decided by Newt Agent under the owner mandate D-20261109-N1; the owner may overrule. Plan in
  PLAN.md, 'AR5 step 2: two long local optimizations that rarely change the result' (commit
  2030e9e). Scripts only, no library file, no key, no default. Why: in run S8 the stock MC
  refinement of indexed candidates and the MC fallback after CMA in neighbour fits are about 62%,
  31% and 58% of the evaluations (big, first gradient sample, mid), and the MC fallback ended
  lower than CMA in 3 to 6 of 3,000 to 4,000 runs where that was counted. Part 0 recomputes these
  from the stored runs. Part 1 replays the stored seed steps with arms M, C60, C200, C0; pick on
  big and the first gradient sample, test on mid, no script runs on mid before the pick commit.
  Part 2 runs S10 on the three samples with the pass, after a base check that reproduces S8
  exactly on big. Pass rules and five predictions are written in the plan before any run.
  Implementer: Sonnet, own worktree, at most 3 heavy processes. Where: PLAN.md AR5 step 2.

### D-20261111-N13 · Track A / AR4 step 1 (task ar4a) · AR4 recast from a learned proposer to a learned verifier; step 1 plan approved (offline)
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Decided by Newt Agent under the owner mandate D-20261109-N1; the owner may overrule, and the
  recast in particular is the owner's to overrule. Plan in PLAN.md, 'AR4: learning where the
  search is uncertain' and 'AR4 step 1: a learned verifier, measured offline against the hit
  ratio' (commit 2030e9e). Why the recast: the premise of H4 (a network proposes where the voting
  fails) does not hold on these samples: AR1 had the truth among the first 5 candidates for 69 of
  69 and 200 of 200 voxels, and in all 31 grid-fallback cases of S8 a candidate near the truth was
  already costed and refined (AR5 step 1). What fails is the test of a result: right and wrong
  ends overlap in hit ratio, a 0.95 flag marks 8.8% and 15.1% of voxels. Step 1 trains a logistic
  regression and a small set network on per-peak overlap shares (train on big and the first
  gradient sample, test on mid, thresholds committed before mid is scored) and compares them with
  the hit ratio for two uses: the flag of unfinished voxels, and accepting a seed without a grid
  search. Offline only; no BFS change. H4 comes back if ray indexing fails on realistic images
  (not tested). Starts after AR5 step 1 is merged. Implementer: Sonnet, own worktree, at most 2
  heavy processes. Where: PLAN.md AR4, AR4 step 1.

### D-20261111-N14 · Track A / AR5 step 2 · AR5 step 2 review: 4 fixes to the documents; no key for 'no MC fallback' alone; next is a whole-layer test of a mixed seed refinement
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Review verdict APPROVE with 4 required fixes (report
  .claude/reports/20261111-1310-review-ar5b.md); all four are applied in a fixes round, no new
  run. On the implementer's two questions, decided under the mandate D-20261109-N1 and the owner's
  to overrule: (1) No library key now for 'no MC fallback after CMA' alone. Run S10 keeps accuracy
  on all three samples (broken minus fixed -2, 4, 4; limit 5) but meets the plan's 20% rule on 1
  of 3 (29.3, 15.7, 14.3% fewer evaluations than S8 with the pass). The plan made a key
  conditional on S10 passing, and I do not waive a rule written before the runs. (2) The next
  step, AR5 step 3, tests a mixed seed refinement on whole layers: CMA with at most 200
  evaluations first, the stock MC optimizer from the same start only when the CMA end hit ratio is
  below 0.95, together with no MC fallback (run S11). The reviewer's what-if from the stored
  replay (not a run): the same seed-test passes as MC (1,596 against 1,592 on the pick samples,
  1,600 against 1,596 on mid) with 57.3% and 53.0% fewer seed-step evaluations. The sample mid is
  no longer held out for this idea; the step says so and the whole-layer runs, with predictions
  written first, are the test. A key for both changes together is considered only if S11 passes
  its rule. Where: icenine_py/docs/autoresearch/PLAN.md (AR5 step 2; AR5 step 3 to be added),
  .claude/reports/20261111-1310-review-ar5b.md.

### D-20261111-N15 · Track A / AR4 step 1 · AR4 step 1 review: 10 fixes; the verifier's seed-acceptance use gets a whole-layer test only on a sample no choice has read; build that sample (held) with its S8 baseline
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Review verdict REQUEST CHANGES with 10 required fixes, 2 of them critical and both in the
  findings page and its derived numbers (report .claude/reports/20261111-1335-review-ar4a.md); no
  leakage and no code defect found. All 10 are applied in a fixes round, no new run. Result as
  reviewed: no model passed the plan's rule for use 1 (flag of unfinished voxels): on mid the hit
  ratio flags 31 of 31 wrong voxels at 1,760 voxels, the set network N1 23 of 31 at 190. Use 2
  (accept a seed that fails the hit-ratio test when N1 scores it high) on the 15 fallback cases of
  mid: 7 of 10 right ends accepted, 0 of 5 wrong; corrected what-if 2,325,301 of 5,002,750 grid
  evaluations. It is usable only as an added rule (accept when the hit ratio or N1 says so); as a
  replacement of the seed test it refuses 20 right seeds. Decided under the mandate D-20261109-N1,
  the owner's to overrule: (1) AR4 step 1 closes as a script-only result after the fixes. (2) Use
  2 as an added rule gets a whole-layer test with N1 frozen as trained, but only on a sample that
  no choice has read, together with the mixed seed refinement of AR5 step 3 if run S11 passes. (3)
  Approved now: AR3 step 11 (task ar3k) builds that sample, held (250 Voronoi grains with
  gradients, built like mid, generator seed 1), and runs its S8 baseline with the boundary pass;
  no other run reads held until the plans of the steps that use it are committed with their
  predictions. The new generator seed makes a different sample for a carry-over test; it is not a
  repeat to firm up a comparison (owner rule of 2026-10-12). Where:
  icenine_py/docs/autoresearch/PLAN.md (AR4 step 1; AR3 step 11),
  .claude/reports/20261111-1335-review-ar4a.md.

### D-20261111-N16 · Track A / AR5 step 3 · AR5 step 3 review: approve, 8 document fixes; no key for the mixed seed refinement before the test on the sample held; that test splits evaluations by stage
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Review of AR5 step 3 (task ar5c): APPROVE with 8 required fixes, all wording or traceability, no
  new run. All 8 are applied in a fixes round (Opus), then the step is merged as a study with
  scripts only. Reviewed numbers: pick X200; run S11 with the boundary pass needs 51.9%, 31.6% and
  39.2% fewer evaluations than S8 on big, the first gradient sample and mid and meets the pass
  rule on all three; on mid broken minus fixed is 5 (15 broken, 10 fixed), exactly at the limit,
  with 36 wrong voxels against 31. The mixed rule was found in a replay that included mid, so no
  sample is held out for it. Therefore no library key is added now. The test on the sample held
  (AR3 step 11) comes first: runs S10 and S11, and S11 plus the verifier N1 as an added
  seed-accept rule, against the S8 baseline, with evaluations by stage and the grid-fallback cases
  paired by seed step, so that the fallback part of the saving (30% and 10% of it against S8,
  untested why) is seen separately. Its plan and predictions are written after held exists and
  before any run on it. Decided under the mandate D-20261109-N1; the owner may overrule. Where:
  icenine_py/docs/autoresearch/PLAN.md (AR5 step 3); review
  .claude/reports/20261111-1503-review-ar5c.md; page
  icenine_py/docs/autoresearch/findings_ar5_mix.md.

### D-20261111-N17 · Track A / AR3 step 11 · AR3 step 11 hand-back (partial): rerun the S8 baseline on the sample held with a cap of 480 minutes
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- The S8 run on held reached its 120-minute cap: all 24,570 voxels had a record, 48,312,402
  evaluations before the pass and 61 grid fallbacks (both lower bounds), and no map or log was
  written. Predictions 1 and 3 of the plan missed. Decision under the mandate: run the same S8 run
  on held again with a cap of 480 minutes and nothing else changed, started by the main session
  from the task branch at e45bb0b; the counters of the capped run are kept as S8_held_cap120. The
  120-minute cap was a budget taken from prediction 1; it is not part of the method. No choice of
  method is made from the capped counters. After the run: a completion round (statistics, the
  pass, pareto rows, the page), then one review of the whole step. The owner can overrule. Where:
  icenine_py/docs/autoresearch/PLAN.md, AR3 step 11; report .claude/reports/20261111-1351-ar3k.md.

### D-20261111-N18 · Track A / AR5 step 4 · AR5 step 4 plan approved with two amendments (cap 480 minutes; starts next to the baseline rerun)
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- The plan of AR5 step 4 (runs S10, S11 and S11v = S11 plus verifier N1 as an added seed-accept
  rule, on the sample held, against S8) was committed at 17a27c1 with six predictions before any
  number of held was read. Approved under the mandate with two amendments made after the hand-back
  of AR3 step 11: the cap of a run is 480 minutes instead of 150, and the step starts now on a
  task branch that merges the branch of step 11 first, with its runs next to the baseline rerun;
  it merges after step 11 and reads no comparison before the baseline has finished. The
  predictions, the pass rule, the frozen X200 rule and the frozen N1 model and cut were not
  changed. The main session starts the long runs so that no agent has to wait on them. The owner
  can overrule. Where: icenine_py/docs/autoresearch/PLAN.md, AR5 step 4.

### D-20261111-N19 · Track A / AR3 step 11 · AR3 step 11 review: request changes, 7 document and script fixes, no rerun; then merge
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- The Opus review recomputed the finished S8 run on held from the cache (54,733,887 evaluations
  before the pass, 59,923,523 with it, 75 grid fallbacks, wrong 95 then 70, the pass fixed 30 and
  broke 5, boundary share 86.0%) and rebuilt mid and held from the builder; no code defect in the
  run path. Decision under the mandate: apply all 7 required fixes of
  .claude/reports/20261111-1900-review-ar3k.md in one Opus fixes round (fixed and broken by
  gradient level, the scope of the held-out statement with the AR5 step 4 runs of N18 named,
  contended labels, the process-count deviation, derived values in summary.json, the history note,
  the cap as an argument of held_chain.sh), no new run, then merge into feature/eval-reduction.
  The owner can overrule. Where: .claude/reports/20261111-1900-review-ar3k.md;
  icenine_py/docs/autoresearch/findings_ar3_held.md.

### D-20261111-N20 · Track A / AR5 step 4 · AR5 step 4 review: request changes; S10 and S11 stand; S11v was a defective arm, repair the hook and rerun it once as S11v2
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Opus review of task ar5d (report .claude/reports/20261111-2000-review-ar5d.md): every number
  recomputed from the raw run files matches. S10 and S11 on held stand: with the pass S11 needs
  23.3% fewer evaluations than S8 and meets the pass rule; S10 needs 3.6% fewer and fails part
  (a). S11v did not test the planned rule: all 24 seeds that N1 accepted were refused by the
  library's seed check afterwards, so no grain grew from them; the run only skipped grid
  fallbacks. Decision under the mandate: (1) one Opus fixes round with the 10 required fixes; the
  hook is repaired in the script layer only, checked on the stored fallback cases of mid, not on
  held; (2) S11v runs once more on held as S11v2 with the frozen model, cut, rule and pass rule,
  cap 480 minutes; the first run stays reported as a defective arm, its pareto rows are removed;
  the rule part of prediction 5 and prediction 6 are judged on S11v2 and the page says the rerun
  came after the defective numbers were seen; this is a bug repair, not tuning and not an extra
  random seed; (3) because S11 passed on held, the library step for one opt-in key pair is planned
  next under its own decision. The owner may overrule. Where: icenine_py/docs/autoresearch/PLAN.md
  (AR5 step 4, amendment N20).

### D-20261111-N21 · Track A / AR5 step 5 · AR5 step 5 plan approved: two opt-in library keys for run S11 (SeedIndexCMAFirst, BFSNeighborCMANoFallback), off by default, parity with the stored script runs
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- S11 met the pass rule on the held-out sample (23.3% fewer evaluations than S8 with the pass,
  reviewed, D-20261111-N20), after passing on the three earlier samples; N14 and N16 made the key
  wait for that test. Decision under the mandate: add two opt-in keys, both off by default, with
  no default change and no recommended setting: SeedIndexCMAFirst n (the mixed seed refinement;
  first step from CMASigma0, stop ratio from SeedIndexStopHitRatio) and BFSNeighborCMANoFallback
  0|1. Tests first. Parity checks against the stored script runs S10 and S11 on big and S11 on the
  first gradient sample (counts, records, map bytes); keys off is bit-identical. Four predictions
  are in the plan, committed before any code. The verifier gets no key. Task ar5e, Sonnet
  implementer; it may share the machine with the one rerun of AR5 step 4. The owner may overrule;
  whether the documentation recommends the keys is the owner's choice. Where:
  icenine_py/docs/autoresearch/PLAN.md (AR5 step 5).

### D-20261111-N22 · Track A / AR5 step 5 · AR5 step 5 review: request changes; six fixes to the page and the tests, no rerun, no library change
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- The Opus review (.claude/reports/20261111-2115-review-ar5e.md) recomputed parity checks 1 to 4
  from the raw files and found them exact; with both keys off the library path is untouched. It
  requests changes: (1) the page must say that the boundary pass inside the S11 library run has
  the stored try counts and evaluations but not the stored contents (one try with another
  candidate on big, costs differ in 4,372 of 18,617 tries on big and 15,813 of 60,423 on the
  gradient sample, ends up to 1.2 degrees apart, changed sets 986 against 983 and 688 against
  681); (2) the unrounded-map cause is shown only for S8L (AR3 step 10) and is untested for S10L
  and S11L, and the page must say so; (3) refusal tests that only the new code refuses; (4) README
  and page state that the hit ratio of the CMA end is one more overlap evaluation per candidate
  that is not counted; (5) more terms defined; (6) a tracked source for the suite counts. Decided
  under the mandate: apply all six in one Opus fixes round on the task branch, no new run, no
  library behaviour change, then merge into feature/eval-reduction. The keys stay off by default
  with no recommended setting. Where: icenine_py/docs/autoresearch/PLAN.md (AR5 step 5);
  .claude/reports/20261111-2115-review-ar5e.md.

### D-20261111-N23 · Track A / AR5 step 4 · AR5 step 4 phase B re-review: approve; three small fixes; the verifier gets one whole-layer step next
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- The Opus re-review of phase B (report .claude/reports/20261111-2215-review-ar5d-phaseB.md)
  recomputed every S11v2 number from the stored files and approves. Three required fixes, no rerun
  and no library change: (1) place the page text of
  .claude/reports/20261111-2330-ar5d-phaseB-page.md as findings_ar5_held.md and confirm the table
  sync check; this waits for the owner, because the implementer's write and the main session's
  copy were both refused; (2) MIGRATION_HISTORY names the commit that places the page; (3) a note
  in s11v2_validity and summary.json that the hook's fitted-voxel fields miss reseed growth and
  that the log count 16 is the figure. Then merge the task into feature/eval-reduction. On the
  question the reviewer left: S11v2 passed both rules on held (37.4% fewer evaluations than S8,
  18.4% fewer than S11), but broken minus fixed against S11 is 4 of a limit of 5 from one run on
  one sample, so no key and no recommendation follow from it. The next plan tests the same frozen
  rule (N1, cut 0.74914) on samples it has not been run on with the whole layer, with predictions
  written first; no extra random seeds. The owner can overrule. Where:
  icenine_py/docs/autoresearch/PLAN.md (AR5 step 4);
  .claude/reports/20261111-2215-review-ar5d-phaseB.md.

### D-20261111-N24 · Track A / AR5 step 6 · AR5 step 6 plan approved: the frozen N1 seed-accept rule on mid and grad, next to a plain hit-ratio rule
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Plan in PLAN.md, section AR5 step 6 (2ed6618), task ar5f, scripts only, no library change. Arms,
  one run each with the boundary pass: S11n (S11 plus N1 at the frozen cut 0.74914 through the
  repaired hook) and S11h (S11 plus a plain rule: keep a failed seed whose hit ratio is at least
  h_c; h_c is fixed in phase A from a pool without mid, as the lowest cut that accepts no wrong
  end). Test sample mid; grad runs too but N1 was trained on it, so its scores are labelled
  in-sample. No sample is fully clean and none is generated: a new sample would be close to a
  repeat with a new random seed, which the owner ruled out. The C++ converged-seed rule is not a
  separate arm because it never fires at these costs. Pass rule unchanged; the reading plain does
  as well (S11h passes on mid and saves at least 80% of what S11n saves against S11) and seven
  numbered predictions are written before any run. Known limit: about 12 failed seeds per sample,
  so counts are small. Cost at most 56.8 million evaluations, four runs side by side. The
  implementer starts only after AR5 step 4 is merged, which waits for the owner to place its page.
  The owner can overrule. Where: icenine_py/docs/autoresearch/PLAN.md (AR5 step 6).

### D-20261111-N25 · Both tracks / review items 14 and 3 · Items 14 and 3 plan approved: both speed keys together on E0 and the pilot, timing A/B on a quiet machine
- 2026-11-11, owner. Recorded by Newt Agent. Supersedes: none.
- Plan in PLAN.md, section Items 14 and 3, task item14; scripts and one config only, no library
  change. Stored arms are read (E0 off 41 wrong, smooth key 7; pilot off 37, smooth 34, FusedCost
  24; all stored pilot times are contended). New runs: E0 with FusedCost alone and with both keys
  (accuracy only); on the pilot the timing arms A (FusedCost) and B (both keys), one worker, alone
  on the machine, order A1, B1, A2, B2 after a warm-up. Repeats of a timing are not extra random
  seeds: each repeat must give byte-identical results and A1 must equal the stored run, or the
  step stops. Decision rules R1 to R4 and seven predictions are written before any run; a run that
  is not quiet reports no time. The measured C++ time is left as a separate step (its BFS runs
  only in MPI mode and needs its own config). AR5 step 6 runs and any merge wait while the timing
  arms run. The owner can overrule. Where: icenine_py/docs/autoresearch/PLAN.md (Items 14 and 3).

### D-20261112-N1 · Track A / C++ time on the pilot (review page item 3) · C++ time plan approved: one client on the pilot region, stored binary, beside the Python timing A/B
- 2026-11-12, owner. Recorded by Newt Agent. Supersedes: none.
- Plan approved under D-20261109-N1 (owner mandate) and N2 (item 3 delegated). The stored Release
  binary runs mode p with one server and one client (mpirun -n 2) on the 2,000-voxel pilot region,
  config generated from key_off_fused.config minus the three Python-only lines. Two timed runs
  under the quiet rule R4 of item 14, plus one untimed run that reads the evaluation counter with
  lldb. No change under Src/ or icenine_py/icenine/. The page must say the two sides do not do the
  same work (C++ one-client has no revisit and no reseed; predicted 740 to 1,300 wrong voxels
  against 24 for Python arm A) and compares unit costs (seconds per global search, per local
  optimization, microseconds per evaluation). A print added in Src/ plus a rebuild, and a Python
  arm with C++ rules, are left as open choices for the owner. Timed runs start only after item
  14's timing arms, never beside them. The owner may overrule. Where:
  icenine_py/docs/autoresearch/PLAN.md, section 'C++ time on the pilot: one client, beside the
  Python timing A/B' (af6514a).

### D-20261112-N2 · Track A / AR5 step 6 (order of work) · AR5 step 6 starts now from the tip of the step 4 branch; its merge waits for step 4's
- 2026-11-12, owner. Recorded by Newt Agent. Supersedes: none.
- The plan (D-20261111-N24) said step 6 starts after step 4 (task ar5d) merges. Step 4 is
  re-reviewed and approved (N23); only the placing of its page text waits for the owner, and that
  does not touch the scripts step 6 reuses. The machine is idle until the owner frees it for the
  timing runs. So the step 6 task branch (feature/eval-reduction-ar5f) is cut from the tip of
  feature/eval-reduction-ar5d (fac79f9), with feature/eval-reduction recorded as its parent. Its
  review reads the diff against fac79f9. It merges only after step 4 has merged. Its four runs are
  not timing results and never run beside a timing run of item 14 or of the C++ step. Nothing else
  in the plan changes. The owner may overrule. Where: icenine_py/docs/autoresearch/PLAN.md,
  section 'AR5 step 6: the frozen seed-accept rule on two more whole layers, against a plain
  rule'.

### D-20261112-N3 · Track A / AR5 step 6 (review) · AR5 step 6 review: request changes, all 11 fixes accepted; wording only, no rerun; prediction 5 counts as not met in full
- 2026-11-12, owner. Recorded by Newt Agent. Supersedes: none.
- The Opus review of task ar5f (diff fac79f9..4c18740) recomputed every headline number from the
  raw run files and found none that differs; frozen.json (b451086, 04:00) was committed before the
  runs started (04:18). Verdict: request changes, 11 required fixes, all to wording on the page
  text, in ver_tables.py and in the MIGRATION_HISTORY entry. All 11 are accepted. Prediction 5
  names no arm in the plan, so it is quoted as written and reported as met for S11n (median growth
  1 and 1) and not met for S11h on grad (4): six of seven predictions met, not seven. The page
  states that nothing was tuned in this step but that N1 was trained and its cut picked on big and
  grad, that the h_c pool holds 16 grad ends, and every earlier use of mid. The claim that the
  plain rule does not replace the verifier is dropped: the page reports 20.4 percent of the saving
  on mid (rule not met) and 62.2 percent on grad (in-sample). No rerun; no number in a summary
  file changes; nothing under icenine_py/icenine/ changes. The fixes go to a fresh Opus
  implementer. The page file itself is still for the owner to place (its Write was refused); the
  fixes are made in the page text. The merge waits for step 4's merge (D-20261112-N2). The owner
  may overrule. Where: .claude/reports/20261112-0510-review-ar5f.md (review report, 'Required
  fixes').

### D-20261112-N4 · Track A / AR5 step 7 (plan) · AR5 step 7 approved: an opt-in key that skips the boundary pass at voxels whose own hit ratio is at least 0.95, tested on held
- 2026-11-12, owner. Recorded by Newt Agent. Supersedes: none.
- With the owner's two pages and the quiet machine still open, the next step is planned and
  started now. A table of the stored runs shows what is left in S11n with the pass: on mid the
  seed step 39.2 percent, the boundary pass 21.7, the grid fallback 16.6, revisits 12.0, neighbour
  local optimizations 10.6; on grad the pass is the largest stage at 38.6 percent. The step chosen
  is a skip rule for the boundary pass: a boundary voxel whose own hit ratio before the pass is at
  least 0.95 gets no try. A replay on five stored passes (mid, grad, big; written with the plan, a
  what-if, not a result) removes 69 to 84 percent of the pass with 0 fixed and 0 broken against
  the full pass; the highest hit ratio of a wrong voxel before the pass is 0.944, so the margin is
  small. The cut 0.95 is frozen in the plan; held is the test sample and its pass files were not
  opened. The library gains one opt-in key, BoundaryPassSkipHitRatio, default 0 (off), off path
  bit-identical. Not chosen now: a library key for the N1 seed-accept rule (it makes the measured
  12.3 percent usable but adds no new saving; it stays an open choice after step 6 closes), and
  the MC part of the seed step (the next target; its log lacks the MC end). The task branch ar5g
  starts from feature/eval-reduction and does not wait for steps 4 and 6. One implementer session
  runs all phases; no run needs a quiet machine and none runs beside a timing run. The owner may
  overrule. Where: icenine_py/docs/autoresearch/PLAN.md, section 'AR5 step 7: skip the boundary
  pass at voxels that already fit' (2a4db6b).

### D-20261112-N5 · Track A / AR5 step 7 (first result, added phase D) · AR5 step 7: the skip fails the plan's rule on held live (4 fixed, 10 broken); the key stays off; a phase D tests whether the order of random draws is the cause
- 2026-11-12, owner. Recorded by Newt Agent. Supersedes: none.
- First result, recorded as it is: the replay of the frozen skip (own hit ratio at least 0.95)
  gives 0 fixed and 0 broken on held and on the five disclosed passes and cuts the pass by 69 to
  84 percent (72.3 percent on held, the run 10.0 percent cheaper). The live pass with the key
  spends exactly the replay's evaluations but its wrong voxels differ from the full pass: fixed
  and broken 4 and 7 on mid, 15 and 12 on grad, 2 and 6 on big, 4 and 10 on held (75 wrong against
  69). The plan's rule allows broken minus fixed of 2 on held, so the skip fails there; prediction
  4 fails on all four samples. The cut is not tuned and the key stays off and not recommended.
  Every flipped voxel is a kept voxel, none is skipped. Two causes are possible and untested: the
  pass uses one random generator, so skipped tries shift the draws of the kept tries; or kept
  tries start from neighbour ends that the full pass had changed. Added phase D, decided after
  seeing this result and to be labelled so everywhere: with the key on, a skipped try consumes the
  draw it would have made, so kept tries get the draws of the full pass; the off path stays
  bit-identical. The four live passes run again with the same seed (no new random seed). If every
  kept try then equals the stored full pass, the draw order was the whole cause, the key keeps
  this behaviour, and the page reports both variants and the size of the draw noise of one pass
  (the flips above). If not, the step stops and reports the first difference. held has now been
  read once live; that is disclosed. The key for the N1 rule and the MC part of the seed step stay
  as in N4. The owner may overrule. Where: report .claude/reports/20261112-1302-ar5g.md in
  worktree agent-a331e7882919e74aa; PLAN.md section 'AR5 step 7: skip the boundary pass at voxels
  that already fit'.

### D-20261112-N6 · Track A / AR5 step 7 · AR5 step 7 review: fixes 1 to 11 accepted; the key keeps aligned draws and stays off; the page placement (fix 12) waits for the owner
- 2026-11-12, owner. Recorded by Newt Agent. Supersedes: none.
- The Opus review (.claude/reports/20261112-1345-review-ar5g.md) reproduced every headline number
  of both variants from the raw files and asked for changes in labels, documents and one test.
  Fixes 1 to 11 are accepted and go to a fresh Opus implementer. Fix 12 (place the page at
  docs/autoresearch/findings_ar5_gate.md) is not done by an agent: the write of that file was
  refused by a guard, so the owner places it or allows it; the fixed page text stays in the task's
  report folder and MIGRATION_HISTORY says the page waits for the owner. The key
  BoundaryPassSkipHitRatio keeps the aligned draws (variant 2: a skipped try consumes its random
  draw) as its only behaviour and stays off by default. The record is: variant 1, the plan's
  pre-set test, failed the rule on held (4 fixed, 10 broken); variant 2 passes on all four samples
  with every kept try identical to the full pass, but held is not a clean test for it, because it
  was decided after held was read live. Whether the documentation recommends the key is the
  owner's choice. For later plans: a live test of a change that only removes work is compared with
  aligned draws, because one shifted draw stream alone moved 4 to 15 voxels each way per sample
  (one re-draw per sample). Where: icenine_py/docs/autoresearch/PLAN.md (AR5 step 7); review
  report .claude/reports/20261112-1345-review-ar5g.md; branch feature/eval-reduction-ar5g.

### D-20261112-N7 · Track A / AR5 step 8 (plan) · AR5 step 8 approved: CMA on every seed candidate before any MC refinement (sweep order), tested on held with aligned draws
- 2026-11-12, owner. Recorded by Newt Agent. Supersedes: none.
- The plan section 'AR5 step 8: fit every candidate with CMA before any MC refinement' (PLAN.md,
  commit 334d2e2) is approved. It changes only the order of the seed step: CMA on every candidate
  first, the MC refinement only if no CMA end reaches hit ratio 0.95. No new threshold, budget or
  model. Reason: in the stored S11n runs 48.0% (mid) and 64.3% (grad) of the MC evaluations buy an
  end that is then discarded, and a replay on stored seed steps cuts seed-step evaluations by
  75.4% (big), 16.8% (grad), 42.1% (mid) with 0 fixed and 0 broken; those three samples are
  disclosed, held is the test sample. Both orders are rerun as whole layers with aligned random
  draws (lesson of step 7, D-20261112-N6), at most 4 heavy processes. Phase A has no library
  change and ends at READY FOR KEY; phase B (the opt-in key SeedIndexCMASweep) needs a further
  recorded decision. Phase A also records the MC traces for a later early-stop step. Not chosen
  now: an early stop inside the MC run (needs a fitted predictor and a replay), a skip of MC from
  the CMA end (needs a new cut). The timing runs (items 14 + 3, C++ time) keep priority if the
  machine becomes quiet: they wait until these runs end, or the runs are started after them.
  Where: icenine_py/docs/autoresearch/PLAN.md (AR5 step 8).

### D-20261112-N8 · Track A / AR5 step 8 (plan) · AR5 step 8 arms run on S11 (merged keys), not S11n, so the step needs no unmerged code
- 2026-11-12, owner. Recorded by Newt Agent. Supersedes: none.
- Amends D-20261112-N7. Both arms of AR5 step 8 (interleaved order and sweep order) are the run
  S11 with the merged keys SeedIndexCMAFirst 200 and BFSNeighborCMANoFallback 1, without the N1
  seed-accept rule. Reason: the N1 code is on the unmerged branches of steps 4 and 6, which wait
  for the owner to place their pages; stacking a third step on them would block it too. N1 acts
  only on failed seeds after the seed step, so it does not touch the order under test. Predictions
  and decision rules are unchanged. The sweep order together with N1 and the boundary-pass skip is
  a later combined run, after those branches merge. Where: icenine_py/docs/autoresearch/PLAN.md
  (AR5 step 8, paragraph Arms).

### D-20261112-N9 · A / AR5 step 8 · AR5 step 8 fails its frozen rule: no key for the sweep order; review fixes 1 to 7 accepted; merge as a study
- 2026-11-12, owner. Recorded by Newt Agent. Supersedes: none.
- The review reproduced every headline number from the raw run files. The sweep order (arm Y: CMA
  on every seed candidate before any MC refinement) passes rules 1 to 4 against the present order
  (arm A) on all four samples and saves 22.6% (held), 38.8% (mid), 16.4% (grad) and 60.1% (big) of
  the seed-step evaluations, 7.7%, 13.4%, 2.1% and 10.6% in total. It fails rule 5, the house rule
  against S8, on held (broken minus fixed 8) and mid (7), limit 5; arm A passes there. By the rule
  frozen before the runs the step fails: no key SeedIndexCMASweep, no phase B, nothing tuned. For
  the record: the nets and the margin of rule 5 are within the draw noise seen before (step 7: 4
  to 15 voxels; Yl against Y on big: 5 fixed, 12 broken), and the comparisons with S8 are not
  draw-aligned; the page says so and does not claim the order loses accuracy. The planning replay
  (0 fixed, 0 broken on seed steps) did not predict the whole-layer result, because a different
  taken end changes which voxels are grown and seeded afterwards; later plans treat a seed-step
  replay as a cost estimate only. Fixes 1 to 7 of the review are accepted; fix 6: the implementer
  writes the findings page, and if its write is refused the text stays in .claude/reports and
  waits for the owner. The step merges as a study with scripts only. Decided under the mandate
  D-20261109-N1; the owner's to overrule. Where: review
  .claude/reports/20261112-2200-review-ar5h.md; branch feature/eval-reduction-ar5h tip e601888.

### D-20261112-N10 · A / AR5 step 9 (plan) · AR5 step 9 approved: rescue a failed seed with a wide CMA local optimization before the grid fallback (rule G1), with a cheap gate; tested on held
- 2026-11-12, owner. Recorded by Newt Agent. Supersedes: none.
- The plan in PLAN.md (4fad072) is approved. When a seed fails the seed test, rule G1 runs a CMA
  local optimization (initial step 0.25 degree, 1,000 evaluations) from each kept seed end in
  order of cost and takes the first end with hit ratio at least 0.95; if none, the grid fallback
  runs as today. The reason: in arm A of step 8, in 22 of 23 fallbacks on grad and mid a kept seed
  end lies within 0.6 degree of the end the grid search finds, and all 23 grid ends are right; the
  grid fallback is the largest stage on held. The what-if saving (21.2% of the run on mid, 27.7%
  on grad) is an upper estimate that assumes the CMA escapes where the seed refinements stopped; a
  replay on the 23 disclosed cases is the gate (at least 12 rescued, no wrong rescued end) before
  any whole-layer run. Values frozen before any run; held is the test sample and none of its seed
  logs or fallback ends was opened; arm A is the stored step 8 run, arm G runs with aligned draws.
  Phase A changes no library file; the key SeedRescueCMA needs a further decision. Not chosen now:
  an early stop of an MC refinement after 2,000 evaluations without improvement (drafted; saves
  12.6% mid, 8.0% big, 0.4% grad) and a library key for the N1 seed-accept rule (its code waits on
  two unmerged branches). Correction to N9: the each-way draw noise of step 7 is 2 to 15 voxels
  per sample, not 4 to 15, and 14 voxels flip on held and mid together (6 before the pass, 8 in
  it), not 13. The step 8 fixes are applied (tip 493f9ec); its findings page waits for the owner
  (the write was refused); it merges as a study. Decided under the mandate D-20261109-N1; the
  owner's to overrule. Where: PLAN.md section AR5 step 9; scratchpad drafts plan_step9_draft.md
  and plan_step9b_draft.md.

### D-20261113-N1 · Track A / AR5 step 9 · AR5 step 9 passes its frozen rules; 8 review fixes accepted; phase B (opt-in key SeedRescueCMA, default off) approved
- 2026-11-13, owner. Recorded by Newt Agent. Supersedes: none.
- Arm G (rule G1: a wide CMA local optimization, 0.25 degree, 1,000 evaluations, at most 5 tries,
  from the kept seed ends after a failed seed test, before the grid fallback) passes rules (1) to
  (5) on held, mid, grad and big. The review (Opus) reproduced every number from the stored runs.
  Totals A to G: held 46,516,396 to 30,768,867 (33.9% fewer), mid 15,645,472 to 13,004,457
  (16.9%), grad 12,189,777 to 10,196,692 (16.4%), big identical. Grid fallbacks 65 to 19, 12 to 4,
  11 to 5. Rescued 46 of 65, 8 of 12, 6 of 11, no wrong rescued seed. Fixed and broken against A
  after the pass: 4 and 0 on held, 0 and 0 elsewhere. Predictions 2 (2.87 tries a case, at most
  2.5 predicted) and 5 (held net -4, range -3 to 5) not met. Close margins: held is 22.4% under S8
  (limit 20%); mid is 1 voxel inside the net limit. Accepted: all 8 required fixes of the review
  (docs, tables, tests; no number changes). Approved: phase B as the plan writes it: one opt-in
  key, SeedRescueCMA sigma evals, default off, off path bit-identical, no golden change, key-on
  library run on mid equal to script arm Gl. The key is not recommended as a default; that choice
  is the owner's. Shared library files change, so a ROADMAP heads-up goes out before the merge.
  Made by Newt Agent under D-20261109-N1; the owner can overrule. Where:
  icenine_py/docs/autoresearch/PLAN.md (AR5 step 9); review
  .claude/reports/20261113-0400-review-ar5i.md; branch feature/eval-reduction-ar5i a30c57d.

### D-20261113-N2 · Track A / AR5 step 9 · AR5 step 9 phase B: the stop rule fired on the in-run boundary pass and was reported; the key SeedRescueCMA is accepted on the before-the-pass equality; 6 review fixes accepted
- 2026-11-13, owner. Recorded by Newt Agent. Supersedes: none.
- The key-on library run on mid equals script arm Gl before the boundary pass in every compared
  item (map byte for byte, 10,064,185 evaluations, 13 rescues tried, 10 rescued). The pass made
  inside the library run differs from the script's pass (3,051,586 against 3,051,647 evaluations;
  8 voxels over 1 degree apart; 37 wrong voxels in both, 4 fixed and 4 broken). Read strictly, the
  plan's stop rule (the key-on library run differs from the script arm) fires. The run stopped and
  reported. The second review (Opus) found: the gap does not depend on the key (a key-off control
  shows it; boundary_pass.py is untouched; the same gap is on record in AR3 step 10 and AR5 step
  5); the in-run pass starts from the unrounded map in memory and the script from the rounded
  stored map; the library's pass on the stored map matches the script exactly; the in-run pass is
  the faithful one and this is not a library bug. Decision: accept the key on the same ground as
  in AR3 step 10 and AR5 step 5, and say openly in the docs and the page that the stop rule fired.
  Later stop rules say 'before the pass'. Accepted: all 6 required fixes of the second review.
  Then merge task ar5i into feature/eval-reduction (key default off, off path bit-identical, suite
  1842 passed with no golden change). The key is not recommended as a default; that and this
  reading of the stop rule are the owner's to overrule. Made by Newt Agent under D-20261109-N1.
  Where: review .claude/reports/20261113-0520-review-ar5i-B.md; branch feature/eval-reduction-ar5i
  107db3c.

### D-20261113-N3 · Track A / AR5 step 10 · AR5 step 10 approved: a level-0 cut of the grid search (quick MC on the first 1,000 listed candidates only), with a cheap gate; tested on held
- 2026-11-13, owner. Recorded by Newt Agent. Supersedes: none.
- After rule G1 the grid search still takes 21.3% of the run on held (19 searches), 10.3% on mid,
  16.6% on grad. The per-level record of 23 stored searches on mid and grad shows a candidate
  within 5 degrees of the end at place 165 or better of about 10,200 in all 23. Step 10 runs level
  0's quick MC only on the first 1,000 candidates of the list sorted by the discrete score (arm C
  = arm G plus the cut). N = 1,000 is frozen, chosen on mid and grad (disclosed); held is the test
  sample. Gate before any whole-layer run: a replay on the 23 cases must be right in at least 22,
  reach hit ratio 0.95 in at least 20, with at least 50% fewer evaluations. Predictions, rules,
  stop rule (library against script compared before the pass) as in the plan section. Phase B
  (opt-in key GridLevel0Top, default off) only after a recorded decision. Not chosen now: an early
  stop of MC refinements, a learned ranking of level-0 candidates, a stronger local optimization
  from G1's near ends. Made by Newt Agent under D-20261109-N1; the owner can overrule. Where:
  icenine_py/docs/autoresearch/PLAN.md (AR5 step 10).

### D-20261113-N4 · Track A / AR5 step 10 · AR5 step 10 fails its replay gate: no whole-layer run, no key; 6 review fixes accepted; merge as a study with scripts only
- 2026-11-13, owner. Recorded by Newt Agent. Supersedes: none.
- The level-0 cut (quick MC on the first 1,000 listed candidates) was replayed on the 23 stored
  grid searches of mid and grad. Right in 20 of 23 (gate: at least 22); hit ratio at least 0.95 in
  15 (gate: at least 20; the full searches reach 21); 2,499,899 evaluations against 7,733,964,
  67.7% fewer (gate met). Three cut searches end 24 to 50 degrees off where the full search was
  right. The review (Opus) reproduced the numbers and found no replay bug; no held file was read.
  By the rule fixed before the run: stop, no tuning, no other N, no whole-layer run, no key. Read
  after the gate (post hoc, not a test): the near candidate is inside the first 1,000 and survives
  the kept quarter in 7 of the 8 searches under 0.95; it is lost after level 0. The cut also
  shrinks levels 1 to 3 (250 starts against about 2,540), which the plan's reading of the record
  did not cover. Lesson: a claim about where the right candidate sits at one level does not show
  that later levels keep it. Accepted: all 6 required fixes (docs, one docstring, the shared
  Wilson function, two tests). No findings page for this step: the MIGRATION_HISTORY section is
  its record. Then merge task ar5j as a study with scripts only. Made by Newt Agent under
  D-20261109-N1; the owner can overrule. Where: review
  .claude/reports/20261113-0625-review-ar5j.md; branch feature/eval-reduction-ar5j 9b4c3ad.

### D-20261113-N5 · Track A / AR5 step 11 · AR5 step 11 approved: stop the seed step's MC refinement after 3,000 evaluations with no improvement (rule P(3000)), with a whole-path gate; tested on held
- 2026-11-13, owner. Recorded by Newt Agent. Supersedes: none.
- The seed step is arm G's largest stage on held (46.0%) and mid (41.7%); about 76% of its MC
  refinements' evaluations on mid come after the last improvement. Step 11 tests rule P(3000): an
  MC refinement in the seed step stops after 3,000 evaluations with no strict gain of the best
  cost (arm S = arm G plus the rule). On stored traces (a cost estimate only) it saves 11.4% of
  the run on mid, 6.1% on big, 0.2% on grad, and changes 24 of 1,607 taken seed ends on mid. k =
  3000 is frozen, chosen on mid and grad (disclosed); held is the test sample. No cost-trace
  number predicts at a stall whether a refinement will improve again (ranking score 0.38 to 0.71,
  chance 0.5), so no learned predictor is used. Gate before any whole-layer run: a replay of every
  seed step with an MC refinement along its whole changed path (the lesson of step 10).
  Predictions, rules and stop rule as in the plan section. Said now: on mid arm G is at net 4 of
  limit 5 against S8, so the house rule there is about even odds; the limit is not moved and a
  fail on mid alone means no key. Phase B (opt-in key SeedIndexMCPatience, default off) only after
  a recorded decision. Made by Newt Agent under D-20261109-N1; the owner can overrule. Where:
  icenine_py/docs/autoresearch/PLAN.md (AR5 step 11).

### D-20261113-N6 · Track A / AR5 step 11 · AR5 step 11 fails rule 5 on mid by one voxel: no key; 6 fixes; merge as a study
- 2026-11-13, owner. Recorded by Newt Agent. Supersedes: none.
- Arm S (arm G with rule P(3000)) passes rules 1 to 5 on held, grad and big and rules 1 to 4 on
  mid, and fails rule 5 on mid: against S8 it has 8 fixed and 14 broken (net 6, limit 5); G has
  net 4, so the clause for a failing G does not apply. The plan said before the run that a fail on
  mid alone means no key. So: no SeedIndexMCPatience key, no findings page, no tuning, no other k,
  no rerun. Evaluations against G: held 5.71% fewer (seed step 14.92%), mid 10.60%, grad 0.14%,
  big 5.01%; wrong voxels S/G 69/69, 37/35, 55/56, 21/21. The review reproduced the numbers and
  the freeze (rule committed 07:08:21, before any held run). Apply the 6 fixes of the review
  report (fix 6 included) and merge the branch as a study. The 597 of 1,181 cut refinements with a
  best under 0.95 are kept as data for a later stop that reads the hit ratio. The owner can
  overrule. Where: .claude/reports/20261113-1000-review-ar5k.md; PLAN.md AR5 step 11.

### D-20261113-N7 · Track A / AR5 step 12 · AR5 step 12 approved: stop a stalled seed-step MC refinement only when its best is already at hit ratio 0.95 (rule PH(3000)); new test sample held2
- 2026-11-13, owner. Recorded by Newt Agent. Supersedes: none.
- The plan of AR5 step 12 (PLAN.md, task ar5l) is approved. Rule PH(3000): at a stall of 3,000
  evaluations, a seed-step MC refinement stops only if its best is at a hit ratio of 0.95 or more;
  else it goes on. From the stored record it saves 7.5% to 8.3% on mid with 7 to 11 changed taken
  seed ends (step 11's rule: 20), 4.5% on big, 0.16% on grad. Gate first: a replay of the whole
  changed path. held was run and read in step 11 and this rule's cuts are a subset of step 11's,
  so held is disclosed; the test sample is a new sample held2, built exactly as held with
  generator seed 2. This is a new sample for a carry-over test, as held was in AR3 step 11, not a
  repeat of a comparison with another random seed; the owner may see this differently and can
  overrule. Cost about 181 million evaluations, about 3.5 hours wall with 4 processes. Said before
  the run: rule 5 on mid stays a risk (G is at net 4 of limit 5), judged at about one in three;
  the limit is not moved; a fail on mid alone means no key. Not chosen now: a stronger local
  optimization from the rescue's near ends (9 cases, no record of what reaches 0.95). Where:
  icenine_py/docs/autoresearch/PLAN.md, AR5 step 12.

### D-20261113-N8 · Track A / AR5 step 12 · AR5 step 12 passes its frozen rules on all five samples: 6 fixes; opt-in key SeedIndexMCStallHit approved (phase B)
- 2026-11-13, owner. Recorded by Newt Agent. Supersedes: none.
- Arm H (arm G with rule PH(3000)) passes rules 1 to 5 on held2 (the new test sample), held, mid,
  grad and big. Against G with the pass it needs 4.33%, 3.84%, 8.11%, 0.16% and 4.53% fewer
  evaluations (seed step on held2 9.21% fewer); wrong voxels are G's on every sample (67, 69, 35,
  56, 21; 0 fixed, 0 broken); against S8 broken minus fixed is -4, -1, +4, -7, -3 (limit 5; S8 is
  not draw-aligned). The review (APPROVE) recomputed the numbers from the raw maps and logs,
  confirmed the freeze (9e315e2, before held2 existed) and that held2 is built as held except
  generator seed 2. Apply the 6 fixes of the review. Phase B is approved: one opt-in key
  SeedIndexMCStallHit k (default 0, off) in the seed step's MC refinement; off path bit-identical;
  no default or golden change; key-on library run on mid equal to a script arm with the library's
  draws before the pass, evaluation counts included. The key is not recommended as a default; that
  is the owner's choice. The owner can overrule. Where:
  .claude/reports/20261113-1600-review-ar5l.md; PLAN.md AR5 step 12.

### D-20261113-N9 · Track A / AR5 step 13 · AR5 step 13 approved: arm H with the boundary-pass skip key (arm HK) as the best current recipe, and the skip's clean test on held2
- 2026-11-13, owner. Recorded by Newt Agent. Supersedes: none.
- The plan of AR5 step 13 (PLAN.md, task ar5m) is approved. Arm HK is arm H (rescue rule G1 plus
  rule PH(3000)) with the step 7 key BoundaryPassSkipHitRatio 0.95. An exact replay on H's stored
  passes gives 17.5% (mid), 29.9% (grad), 23.9% (big) and 12.8% (held) fewer evaluations than H
  with 0 fixed and 0 broken voxels, and 57% to 67% fewer than S8. held2's pass tries were never
  opened and 0.95 was fixed in step 7 before held2 existed, so held2 is an untouched test of the
  skip; step 7's test on held was not clean. Gate first: the live pass on mid must equal the
  replay try by try. No library change; about 15 million evaluations. Phase A may start before
  step 12's merge (it needs only step 12's cache); phase B (one recipe run on mid with all three
  keys) waits for that merge. Not chosen now: a local optimization from the rescue's ends to avoid
  grid searches (needs a probe replay first; a candidate for the step after). Where:
  icenine_py/docs/autoresearch/PLAN.md, AR5 step 13.
