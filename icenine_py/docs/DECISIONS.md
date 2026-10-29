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
