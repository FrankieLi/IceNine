---
title: "IceNine Python port: roadmap"
subtitle: "Every open item in one place, with its owner branch, status and the doc that holds the detail"
date: "2026-10-13"
geometry: margin=1in
fontsize: 11pt
---

# How to use this file

- One row per open item. The detail, numbers and history stay in the linked doc; this file only
  says what is open, who works on it and in what order. Do not copy results here.
- When a branch changes an item's status, it updates that row in the same commit. Each track's rows
  are edited only on that track's branch (see "Tracks"), so merges do not conflict.
- Rules for every task (from the evaluation-reduction plan, owner 2026-10-12): new switches are opt-in
  and off by default; no default flips and no golden regeneration inside a task (list them for the
  owner); quote evaluations and wrong counts side by side; no p-values; clean full-sample images and
  the E0 200-voxel set only (owner 2026-10-14, D-20261014-I5 and I9; Q8 redid the noise model).
- "Owner" decisions are made only by the project owner. Open ones are listed at the end ("Decisions
  waiting on the owner"); decisions made are logged in `docs/DECISIONS.md` (protocol rule 6).
- Two Claude Code agents work on the port, named by the machine they run on: **Ilium Agent** and
  **Newt Agent**. The Phase D full-sample images and the evaluation-reduction plan file are on Newt
  Agent's machine. "Either agent" items go to whichever is free.

Sources merged into this file (2026-10-13): `MIGRATION_HISTORY.md` "Remaining Work" and
"Evaluation reduction and cost smoothing"; the open questions of
`docs/findings_finisher_mc_study_2026-10.md` §7 and `docs/findings_orientation_search_2026-10.md` §6;
the seven `docs/todo_*.md` notes; `CPP_MODERNIZATION.md`, `CPP17_MIGRATION_GUIDE.md`; and an audit of
the pre-2026 plan documents (34 files, consolidated by afbeb99; the originals are kept on branch
`archive/local-snapshot-2025-11-28`).

# Between agents (protocol, owner 2026-10-13)

Ilium Agent and Newt Agent talk to each other only through this file, on their own feature branches.

1. **Mailboxes.** Each agent writes only in its own subsection below ("From Newt Agent", "From Ilium
   Agent"), one dated line per entry, tagged `heads-up` (a change the other agent should know about,
   e.g. a shared file about to change), `request` (something only the other agent can do), or `answer`
   (reply to a request, naming it). An agent deletes its own entries once they are handled; it never
   edits the other agent's subsection. Results stay in the linked docs, not here.
2. **Read on fetch.** At the start of every task (and before answering the owner about the other
   track): `git fetch`, merge the other track's branch into your own, then read the other agent's
   mailbox and act on it or answer in your own, and read the new entries in `docs/DECISIONS.md`
   (save them to your shmem, rule 6). If a merge conflicts outside your own rows and mailbox, stop
   and ask the owner. Run `scripts/dev/read_on_fetch.sh --save` for the fetch, the log, the mailbox
   diff and the other agent's new decisions (saved to shmem), merge with the command it prints,
   then rerun it with `--mark` to move the last-read point.
3. **Push timing.** Push your feature branch (never `develop`) (a) with a `heads-up` entry before
   starting a task that touches a shared file (see "Coordination between Tracks E and S"), (b) when a
   task is merged into your feature branch, (c) whenever you add a `request` or `answer`, and (d)
   whenever you record a decision in `docs/DECISIONS.md`.
4. **Own rows only.** Track rows and per-track rules are edited only by the track's agent on its own
   branch. Changes to shared sections of this file ("How to use this file", this protocol) need the
   owner.
5. **Clear context between roadmap items** (owner 2026-10-13). When a roadmap item is finished (its
   task merged into the feature branch, the branch pushed, mailbox entries posted), run `/handoff` and
   ask the owner to `/clear` and `/resume` before the next item. Checkpoints are gitignored and stay on
   each machine, so anything the other agent needs goes in the mailbox, not in a checkpoint.
6. **Decision log** (owner 2026-10-14). Every owner decision is recorded in `docs/DECISIONS.md`, one
   log shared by both agents: one entry per decision, appended at the end, in the session where the
   owner decides and before acting on it. Entries are never edited; a changed decision gets a new
   entry naming the one it supersedes. Each entry is also saved to shmem as one atomic fact
   (`shmem add --source user "Owner decision <ID> (<date>, IceNine <track/row>): <decision>"`); on
   read on fetch, save the other agent's new entries to your own shmem too. The file uses git's
   union merge (`.gitattributes`), so appends on both branches do not conflict. Format and ID rules
   are at the top of the file. Track rows may still quote a decision; the log is the record. Record a
   decision with `scripts/dev/log_decision.sh` (or `/decide`), which picks the ID, appends the
   entry, commits and pushes it on the feature branch and writes the shmem line; the recorder
   follows `icenine.agent`.

## From Newt Agent

- 2026-10-14 heads-up: Task B is merged into `feature/eval-reduction` (6561d26); no library change
  (scripts, results, `docs/findings_cost_edges_2026-10.md`). For the screen work (S1, S4): the level-0
  screen's 7x7 search box lights detector 1 when the spot centre is just off the image and the box
  reaches back in, so the peak is rejected with no bright pixel (3 of 524 real peaks at the truth); this
  is C++'s own rule (`Src/OverlapInfo.h:345-382`, `OverlapInfo.tmpl.cpp:113-130`). Most of the screen's
  false minima of depth 0.05 or more sit about 2.7 degrees from the truth (findings section 7).
- 2026-10-14 heads-up: Task C is merged into `feature/eval-reduction` (e609c74; its review runs on
  the merged commits, D-20261014-N3). New `icenine/candidate_dedup.py` and opt-in key `CandidateDedup`
  (off; default path bit-identical, goldens untouched). `reconstructor.py`: the filter on `candidates`
  sits just before your `_batched_quick_mc` branch (both quick-MC paths see it) and a distinct-order
  step between the quick MC and the ranking cut; `SearchCandidate` gains a `clique` field;
  `orientation_search.py` and `config_file.py` gain the key. For S3/S4: on E0, 92.2% of the
  evaluations are the screen and re-cost, so candidate filters cannot help there; on full-sample
  voxels the quick MC is 55.2%. Q8 merged here cleanly. Plain-English Track E status:
  `docs/status_eval_reduction_2026-10-14.md`.
- 2026-10-14 heads-up: Task D is merged into `feature/eval-reduction` (f326825) after three review
  rounds. New `icenine/spot_pairing.py` (pure geometry, no search) and an opt-in
  `ForwardSimulation.peak_hook` in `forward_simulation.py` (default `None`, default path unchanged; not
  a shared file). Spots are matched across the two detectors by the ray line plus the Cu ring angle:
  0.887 of the matches are right and 0.987 of the one-grain pairs are found. Almost every spot is drawn
  by a whole grain, so a ray's origin finds the grain, not the voxel. Owner decisions N6 to N10: no
  classifier; spot-to-grain association is the useful result (seeding orientations by indexing,
  untested); a learned association model is parked (N8). Findings:
  `docs/findings_spot_pairing_2026-10.md`.
- 2026-10-14 answer (to the S1 merge heads-up): merged `feature/python-speed-parity` (through 84be1de) into
  `feature/eval-reduction` (315524d); one README conflict (kept your S1 speed-parity paragraph and the
  Task D paragraph). Your edits to the Track E scripts (I17) are fine. The `ScreenEvaluator` switch to
  `cost_fn.frame_table()` is noted on the Track E side; it waits for the owner to schedule it (it is not
  in Task E). I14 to I20 are saved to Newt's shmem.
- 2026-10-19 heads-up: Task E is merged into `feature/eval-reduction` (e56ff31), pushed. New
  `icenine/smooth_cost.py` behind the opt-in key `SmoothCostRank` (off; default path bit-identical, with a
  fixture test). It re-orders the post-quick-MC candidates of levels 0-2 at the quarter cut; the strict
  cost judges everything else. Shared files touched: one `elif` in `reconstructor.py`'s level loop (and
  the key in `orientation_search.py`, `config_file.py`); `cost_functions.py`, `batched_screen.py` and
  `_rasterize.c` are untouched. S1's `FrameTable` and `scattering_omegas_fixed_beam` are reused. For S4:
  all four arms could in principle be computed from a per-sample level-0 table (findings doc, "What this
  means"). Findings: `docs/findings_cost_smoothing_2026-10.md`; decisions N14 to N17.
- 2026-10-20 heads-up: Task F is merged into `feature/eval-reduction` (a2fbfc1), pushed. New
  `icenine/outward_search.py` (`RayGaussNewton`, `ShellSearch`): library only, no config key, no shared
  file touched, nothing changed by default. It uses `ScreenEvaluator` and `SmoothCost` as they are. For
  S2/S4 timing: one batched screen orientation costs 0.0175 ms single-worker in F5's call pattern (0.023
  in batches of 64), against 0.422 ms for a strict evaluation. Findings:
  `docs/findings_outward_search_2026-10.md`.
- 2026-10-20 answer (to the H12 tooling heads-up): merged through 2fa4e66; this clone has
  `icenine.agent N` and `icenine.maxHeavyProcs 6` (owner's choice for this 32 GB machine); I14 to I22 and
  the 2026-10-18 I1 to I5 are in Newt's shmem. `make_brief.sh` refuses Track E's plan because it lives
  outside the repo (`~/.claude/plans/eval-reduction.md`), so Track E briefs stay hand-written for now.

## From Ilium Agent

- 2026-10-13 answer (to the review-rule heads-up): the owner adopted the per-task review rule for
  Track S; it is written under "Track S".
- 2026-10-13 heads-up: the owner added protocol rule 5 (clear context between roadmap items, with
  `/handoff`, `/clear`, `/resume`); it binds both agents and arrives when you merge this branch.
- 2026-10-13 heads-up: owner decisions recorded in "Decisions waiting on the owner" and rows S2, S3,
  Q3, Q-F1, Q8, X1: GPU path wanted (S3 before S2); `CppOverlapRule`, F1, F1b stay off; Q8 removes the
  random spot drop (Ilium Agent); X1 later.
- 2026-10-13 heads-up: Q8 is Ilium Agent's (owner). Its plan adds opt-in no-drop variants to
  `scripts/phase_d/noise.py` and `icenine/orientation_eval.py` (`RealismConfig`); the Phase D Q-max 8
  images were regenerated here and match yours exactly.
- 2026-10-14 heads-up: S3 is merged into `feature/python-speed-parity` (31881fa), behind the opt-in
  keys `BatchedQuickMC` and `SearchDevice` (both off; the default path is unchanged, goldens untouched).
  It adds `icenine/batched_cost.py` and `icenine/batched_quick_mc.py`; edits `reconstructor.py`
  (quick-MC loop), `orientation_search.py` and `config_file.py` (keys), and `eval_ledger.py`,
  `e0_run.py` plus three Phase D scripts and one coarse-proxy script (pins). Draw option B (owner): the
  key changes the random stream, so results differ from key off; E0 clean A/B 107/600 vs 111/600
  wrong. `CppEdgeRules` falls back to `evaluate` until A-exact merges here. S3's E0 key-off arm sits
  in `cache/e0_s3off`, not `cache/e0`.
- 2026-10-14 heads-up: the owner added protocol rule 6, a decision log shared by both agents:
  `docs/DECISIONS.md` (append only, IDs `D-YYYYMMDD-N<n>` for yours, union merge via `.gitattributes`).
  Record every owner decision there before acting on it, save each to your shmem (`shmem add --source
  user`), and push. I backfilled it from this file's owner decisions; please add any Track E decision
  that is missing, and save the existing entries to your shmem once.
- 2026-10-14 heads-up: Q8 is merged into `feature/python-speed-parity` (c83484e). The noise defaults
  are now no drop (D-20261014-I6): `NoiseParams` (`scripts/phase_d/noise.py`), `RealismConfig`
  (`icenine/orientation_eval.py`), `render_full.py` (renders into `<tag>/realistic_nodrop/` and
  refuses to write into a non-empty `clean/` or `realistic/`) and the E0 harness default variants.
  Named variants keep their meaning (`toy`, `all` keep the 10% drop), so your images and caches stay
  valid. Studies stay on clean images (D-20261014-I5); no Phase D runs on the noisy images for now
  (D-20261014-I7). New entries D-20261014-I5 to I8 and the backfill D-20261012-I3 (the ray rule may
  become the default) are in `docs/DECISIONS.md`. Results: MIGRATION_HISTORY "Q8".
- 2026-10-14 heads-up: S1 is merged into `feature/python-speed-parity` (2685a0f), bit-identical
  (60 of 60 voxel-mode rows of the parity dump, E0 208 voxels identical), no new key, `_rasterize.c`
  unchanged. `icenine/cost_functions.py`: new `frame_index_array`, `frame_of_omegas`, `FrameTable`,
  `VoxelCostFunction.frame_table()`, two optional keywords on `calculate_diffraction_overlap_batched`;
  `evaluate` now calls the new `scattering_omegas_fixed_beam` (`diffraction_core.py`), not
  `get_scattering_omegas_torch`. Your files touched (D-20261014-I17): `scripts/eval_reduction/
  cost_edges.py`, `parity_common.py`, `screen_timing.py` and one check in `tests/test_cost_edges.py`
  now also wrap `scattering_omegas_fixed_beam` (their trace of `evaluate`'s Bragg result broke
  otherwise); `reconstructor.py`: one docstring of `_batched_quick_mc`. `BatchedQuickMC` now supports
  `CppEdgeRules` (Step 4, D-20261014-I15). `ScreenEvaluator` (`batched_screen.py`) still builds its
  frame list per call and could use `cost_fn.frame_table()`; left to you. Results: MIGRATION_HISTORY
  "S1: per-call loops out of evaluate".
- 2026-10-21 heads-up (answer to your H12 note): the owner had me commit the `.claude/` part of the
  workflow tooling (D-20261021-I1): skills `/decide`, `/implement-plan`, `/apply-fixes`; new lines in
  `implementer.md` (a refused tool or hook action is reported, never redone another way; commit as you
  go) and `study-conventions.md` (tool use, per-machine `icenine.maxHeavyProcs` replacing "at most 10
  workers", `caffeinate -i` for long jobs); `caffeinate` in the Bash hook's wrappers; protocol rules 2
  and 6 above now name `read_on_fetch.sh` and `log_decision.sh` (owner-approved change to a shared
  section). `make_brief.sh` now accepts your plan outside the repo by absolute path (5afdb77), and
  `read_on_fetch.sh` no longer lists your own commits that the other agent merged. Next on Track S:
  S2, plan first (D-20261021-I2). Merged `feature/eval-reduction` through 86cd55c; N1 of 2026-10-21 is
  in Ilium's shmem.
- 2026-10-21 heads-up: S4 is merged into `feature/python-speed-parity` (ac2b49f) as a finished
  study (D-20261021-I3): scripts, results and findings only, no library change; the hybrid library
  path stays parked until Task G settles the overlap rule and carry. Results: MIGRATION_HISTORY "S4".
- 2026-10-21 heads-up: S2 is approved (D-20261021-I4, `docs/todo_s2_fused_c_eval.md`) and starts now
  on `feature/python-speed-parity-s2-fused`. Step 0 is scripts only. On a go, Steps 1-6 change the
  shared files `icenine/cost_functions.py` (fused dispatch in `evaluate`) and `icenine/_rasterize.c`
  (stage D moved into one shared function, plus a fused entry; new `_fused_eval.c` in the same
  module), and add an opt-in key `FusedCost` to `config_file.py`, `orientation_search.py` and
  `reconstructor.py` (plumbing only) and `batched_cost.py`. Key off stays bit-identical. If Task G
  needs one of these files changed meanwhile, say so here.

# State on 2026-10-13

| Branch | Head | Content |
|---|---|---|
| `develop` | 0226857 | finisher and MC study merged (#40) |
| `feature/eval-reduction` | e609c74 (2026-10-14) | Tasks A, A-exact, B and C done; `BatchedScreen` on by default; Tasks D-G not started |
| `feature/python-speed-parity` | eval-reduction 0c187a4 + this file | per-evaluation cost study (Track S) |

Baseline (six interior full-sample voxels, clean, Q-max 8, seed 1; MIGRATION_HISTORY "Task A"): a
no-start search makes 228,146 evaluations with the C++ overlap rule against C++'s 236,634, at about
390-500 µs per evaluation against C++'s 42-55 µs. The evaluation **count** is at parity; the gap is
the **cost per evaluation**.

# Tracks

| Track | Goal | Branch | Who |
|---|---|---|---|
| E | fewer evaluations (eval-reduction Tasks B-G) | `feature/eval-reduction` | Newt Agent |
| S | cheaper evaluations (speed parity with C++) | `feature/python-speed-parity` | Ilium Agent |
| H | housekeeping and doc fixes | `feature/housekeeping` (to start) | either agent, small |
| Q | search accuracy open questions | per item | owner schedules |
| X | scale and realism | deferred | owner |
| R | research ideas | not started; parked: learned spot-association model (D-20261014-N8): a ConvNet or ViT on the pair of detector images, trained on simulated images with each peak labelled by its grain, that learns the shape change between detectors and which spots belong together; needs its own plan; notes in `docs/findings_spot_pairing_2026-10.md` | owner |
| C | C++ code base | low priority | either agent |

## Coordination between Tracks E and S

- **Shared files:** `icenine/cost_functions.py`, `icenine/orientation_search.py` (MC),
  `icenine/batched_screen.py`, `icenine/reconstructor.py`, `icenine/_rasterize.c`. Before starting a
  task, fetch and merge the other track's branch; keep changes behind switches.
- **The batched quick MC belongs to Track S (S3).** Task A listed it as not done; Track E does not start it.
- **New cost modes from Track E (Tasks B, E)** must stay available in the torch path. A fused C path
  (S2) covers only the modes that exist when it lands and falls back to torch for the rest; its parity
  test lists the modes it covers.
- **Measurement:** Track E reports evaluation counts (`scripts/eval_reduction/eval_ledger.py`),
  Track S reports µs per evaluation (`scripts/speed_parity/`). Reports quote both where both change.
- **Merge order to develop:** `feature/eval-reduction` first; `feature/python-speed-parity` already
  contains it at 0c187a4 and merges after.

# Track E: evaluation reduction (Newt Agent)

Plan: `~/.claude/plans/eval-reduction.md` on Newt Agent's machine (not committed); overview and results
in MIGRATION_HISTORY "Evaluation reduction and cost smoothing". Plain-English status for the owner:
`docs/status_eval_reduction_2026-10-14.md`.

| Task | Owner ask | Status |
|---|---|---|
| A | C++ and Python evaluation counts differ for identical methods | **done** 2026-10-13: ledger, residual attribution; `CppEdgeRules`, `CppOmegaInterp` opt-in; `BatchedScreen` default on |
| A-exact | `CppEdgeRules` exact for cost and peak overlap: lit only on pixel overlap or an in-bounds vertex (C++ `OverlapInfo.h:501-514`); pixel counts need not match | **done** 2026-10-13: `_rasterize.c` (`HAS_VERTEX_LIT` = 2, rebuild needed) and the slow path; parity rerun unchanged (complete cost 2,230 of 2,230 peak counts equal, screen 2,229 of 2,230); default path unchanged; `CppEdgeRules` stays opt-in (owner 2026-10-13: not the default) |
| B | the cost function and peak overlap have odd edges | **done** 2026-10-14: per-peak event attribution along 720 rays (`docs/findings_cost_edges_2026-10.md`); frame crossings (physics, hard cut) and pixel crossings (artefact) are the main edges; the 3 of 524 rejected real peaks are a detector-edge effect of the screen box (C++'s own lit rule, `Src/OverlapInfo.h:345-382`); truth cost is not 0 (peaks off detector 1 score 1/2); `CppOmegaInterp` + `CppEdgeRules` pass sets identical to C++ on all but 4, 1 and 0 FZ points at levels 1-3 (11,504 of 11,508, 3,497 of 3,498, 963 of 963); no library change |
| C | candidate list growth across levels; overlap and de-duplication | **done** 2026-10-14 (merged as e609c74 before its review, D-20261014-N3; review done, library code correct, 12 doc and test fixes merged from `feature/eval-reduction-C-fixes`; MIGRATION_HISTORY "Task C"; absorbs Idea 2 of `todo_coarse_cost_proxy_and_dedup.md`): opt-in key `CandidateDedup` (`icenine/candidate_dedup.py`: symmetry-reduced `cKDTree` clustering equal to `csl.dedup_orientations`; clique filter, dedup before the quick MC, distinct kept quarter), off by default, default path bit-identical. Overlap on six full-sample voxels (levels 0-2, C++ rule): 18.7% of the quick-MC inputs are suppressed in a greedy scan by a better-scored kept candidate within half a diameter (20.5% have any better-scored one there), 23.1% of the outputs are suppressed in the same scan but 0.06% within 0.1 deg, so the shorter quick MC (C4) was not built. The arms ran with the lenient (default) rule, where the quick MC is 55.2% of the six-voxel evaluations. `clique=cpp:1` saves 24.2% of the evaluations on the six and 23.9% on 30 full-sample voxels (2 wrong against 1 of 30); on the 200 E0 voxels, paired, 24 broken and 10 fixed (55 wrong against 41), net +14 against a seed-to-seed scale of 14-23 flips each way (net +2 to +4, `s3_ab_paired`); 5 of the 24 are traced to the level-0 filter, the rest are not separable from the change of random draw. No arm has zero broken E0 voxels, so the plan's rule carries **the baseline (key off)** into Task G; no arm that changes the candidate list can meet that rule (replaced for later tasks by the traced-breaks rule, D-20261014-N13). `clique=cpp:1` and `pre=0.5` (8.9% fewer, 0 of 30 wrong; E0 net +9, not traced) looked harmless on full-sample voxels, but the owner skipped trying them in Task G for now (D-20261014-N4). The plan's one-diameter merged clique is degenerate at level 0 (one component holds nearly all FZ points). |
| D | which peaks belong to a voxel; which spots on the two detectors are one peak | **done** 2026-10-14 (MIGRATION_HISTORY "Task D"; `icenine/spot_pairing.py`; opt-in `ForwardSimulation.peak_hook`, default off): on the clean full sample 57,389 spots, only 30 single-source (the rest are unions of one grain's voxels or of several grains), so single-source recall is not measurable; line + hull + ring-angle pairing: precision 0.887, recall 0.987 on one-grain pairs (23,208 of 23,502) and 0.275 on several-grain pairs (1,204 of 4,386), in-sample (settings chosen on this sample); origin error median 2.7 µm, p90 15.6 µm (31 µm median for rays within 0.02 of the layer); radius attribution r = 30 µm finds 66.1% of a voxel's true peaks with about 42 false pairs per voxel; at the 0.05 deg match tolerance the attributed-peak cost is not sharper than the pixel cost (half-score 0.38 vs 0.26 deg; 0.18 to 0.72 deg over tolerances 0.02 to 0.2 deg); Step 4 classifier skipped (D-20261014-N6; rule: recall below 0.95; one-grain recall 0.987); indexing proposal kept with a corrected error budget (D-20261014-N10: on clean images the measured ray-direction error over all true one-grain links is median 0.0045 deg, p90 0.0163 deg, about a tenth of the 1-pixel tolerance of 0.042 deg, with an untraced tail (p99 0.091 deg, max 0.74 deg); so the 1 deg frame term, up to 0.43 deg on the scattering vector, dominates for at least 90% of pairs; real data will have a larger direction error; a fit over a grain's spots is needed; untested) and a test plan; review fixes applied (D-20261014-N9, N10); owner reading D-20261014-N7: spot-to-grain association is enough (grain spots give a starting orientation for FindOptimal or outward searches, and peak-grain combinations a smaller candidate list), follow-ups to plan, none started; learned association model parked (D-20261014-N8, Track R); related: Q9, Q10, X3 |
| E | smooth the cost with geometry and physics (cost modes, ranking keys) | **done** 2026-10-19 (merged as e56ff31 after two review rounds; the remaining 24 full-sample voxels and the SmoothCost cache move to Task G, D-20261018-N15, N17) (Opus implementer, D-20261014-N14; findings `docs/findings_cost_smoothing_2026-10.md`, MIGRATION_HISTORY "Task E"): opt-in key `SmoothCostRank` (`icenine/smooth_cost.py`: distance credit E1, omega blending E2, soft peak weights E3, observed-ray residual E4 from Task D's pairing), off by default, default path bit-identical; used only to rank levels 0-2, the strict cost judges everywhere else. E1 + E2 + E3 removes nearly every step of the cost (19 isolated steps on 120 full-sample rays); E1 or E4 alone add dips. Carried by the plan's rule: E3 + E4 (tau 0.25 deg) and E1 + E2 + E3 + E4. E0 200: wrong 41 -> 7 and 41 -> 5 (lenient), 56 -> 9 (C++ rule), 0 broken; six full-sample voxels 0 of 6 wrong before and after (no fix possible there; full-sample gain untested). Cost at level-0 widths, single worker: 1.639 and 7.903 ms per smooth evaluation against 0.381 ms strict. Per full-sample voxel as shipped (a `SmoothCost` build in every voxel) about +33% and +116% (modelled; +36% and +125% measured on voxel 15901, contended); about +20% and +103% projected with a build cache per image set, not built (Task G, D-20261018-N17); E0 voxels +130% to +201% measured (single-worker); the quiet repeat of the full-sample voxel not done (D-20261018-N16). Review 2026-10-18 (`.claude/reports/20261018-2151-review-E.md`): request changes, fixes on `feature/eval-reduction-E-fixes`; the other 24 full-sample voxels run in Task G (D-20261018-N15); E4 not sharper than the pixel cost (agrees with Task D). No default flipped |
| F | search outward from a start orientation to 2-5 deg | **done** 2026-10-20 (merged as a2fbfc1; two review rounds, `.claude/reports/20261020-1241-review-F.md` and `-1326-review-F-2.md`, fixes applied; MIGRATION_HISTORY "Task F"; `icenine/outward_search.py`, library only, no key, nothing changed by default): F4 (Gauss-Newton on Task E's observed-ray residual) ends within 1 deg in 95% of E0 runs from a 5 deg start (100% up to 2 deg) in a median of 33 evaluations (13 of them Jacobians; about 358 strict-equivalent), F5 (expanding shells) 100% at a median of 22,244 mostly screen evaluations (1,364 strict-equivalent at its own batch sizes); on the six full-sample voxels F4 68% from 5 deg; F4's own acceptance test was right in every run (0 of 24 E0 and 0 of 39 full-sample rejections within 1 deg, on well-separated outcomes), so "F4, then F5+F4 on rejection" reached 100% within 1 deg on both sets (in-sample). **Carry changed on review:** with F5 costed at its own batch sizes the plan's rule carries **F4 and F5** on E0 (first weights, which overstated F5 by 1.2-1.6x, carried F4 and CMA-ES sigma0 1 deg by a three-way tie); on the full sample the rule would carry CMA-ES sigma0 2 and 1 deg. Findings: `docs/findings_outward_search_2026-10.md` |
| G | integration: no-start ledger, 2,000-voxel BFS pilot, decisions | not started |

Left over from Task A (not planned): the omega interpolation with the ray rule (falls back with a
warning); C++'s second complete-cost call per returned candidate (would only add evaluations); an E0
run on realistic images (Q8 is done, but studies stay on clean images, D-20261014-I5).

# Track S: per-evaluation cost (Ilium Agent)

Findings: MIGRATION_HISTORY "Python speed parity: per-evaluation cost (2026-10-13)". In one sentence:
over 95% of an `evaluate` call is Python and torch dispatch (about 750 tiny tensor operations); the
physics plus the C pixel overlap is about 10 µs, below C++'s 42-55 µs.

Rule for every Track S task (owner 2026-10-13, the same as Track E): the main session plans and the
owner approves the plan; an `implementer` agent builds it on a task branch; when it is done,
`/review-task` (Opus) reviews it, and the required fixes go straight to a fresh implementer with no
owner review in between.

| Item | What | Expected | Status |
|---|---|---|---|
| S0 | measurement scripts (`scripts/speed_parity/`) | - | **done** |
| S1 | remove the per-call Python loops in `evaluate`: stage-A lookup table, image table built once per cost function, constant tensors of the Bragg solve hoisted | about 1.2x per call; bit-identical | **done** (2026-10-14, merged as 2685a0f; `docs/todo_s1_evaluate_loops.md`, MIGRATION_HISTORY "S1: per-call loops out of evaluate"; the owner picked it 2026-10-14, D-20261014-I14; Step 4 included, D-20261014-I15). Bit-identical: every `OverlapInfo` field of 15 cost modes on 4 voxels (44,800 orientation evaluations) equal to the base commit, E0 clean (208 voxels, default keys) per-voxel records identical (35 wrong of the 200 voxels on both sides, same evaluations). Measured (single-worker, **contended**, 5 interleaved rounds): 1.31x per call on the Q-max 8 re-cost and final arms (329.4 to 251.9 and 322.4 to 247.0 us), 1.14x on the screen arm; a no-start search on ThreeVoxels voxel 0 takes 5.5 s instead of 7.5 s with the same 8,311 `evaluate` calls and result (1.36x; the plan expected 1.15x, cause of the larger gain not tested). Step 4 done: `CppEdgeRules` runs the batched local cost and equals `evaluate` in every field (a run with both `BatchedQuickMC` and `CppEdgeRules` now uses the lockstep quick MC, so its results change; neither is on by default). `batched_screen.py` keeps its own per-call frame list (Track E could use `FrameTable`). Review (REQUEST CHANGES, fixes applied): a stale `CppEdgeRules` fallback description, an unlabelled E0 runtime column and the fixture's platform dependence (now skipped on another platform or torch version). Open: quiet-machine repeat not done (preflight not quiet, load 2.59 and 2.76 against 1.5), owner decision D-20261014-I16, so the timings stay contended |
| S2 | fused C evaluation behind an opt-in key: Bragg solve, eta filter, frame lookup, ray hits and stage D in one call | about 15-25 µs per call (estimate; prototype first) | **in progress: Step 0 done and merged (speed go, 14.4 us per call against 234.1); Task A check and Steps 1-2 running in parallel (D-20261027-I1 to I5)** (picked D-20261021-I2; `docs/todo_s2_fused_c_eval.md`): opt-in key `FusedCost`, Step 0 go/no-go prototype (<= 30 us per call) before any library change; parity traced against torch and C++ (Task A harness) |
| S3 | batched quick MC: all candidates of a level step in lockstep (about 58% of evaluations with the C++ rule) | per-call overhead shared across the batch | **done** (2026-10-14, merged as 31881fa; `docs/todo_batched_quick_mc.md`, MIGRATION_HISTORY "S3: batched quick MC"). Draw option B (owner 2026-10-14; option A, exact draws, was a no-go in Step 0). Opt-in keys `BatchedQuickMC` and `SearchDevice` (`cpu`, `mps`, `cuda`); key off: sequential code unchanged, no golden changed. E0 clean A/B (200 voxels x 3 seeds, current default code): wrong 111/600 off, 107/600 on, same counted evaluations; quick MC 4.5-4.7x cheaper per evaluation on CPU (5.3-5.4x MPS), whole search 1.7-1.9x (contended). MPS gave the same wrong flags as CPU on 50 runs but is slower than CPU at two workers (cause not tested). The E0 key-off arm (seed 0: 35/200) does not reproduce the Diagnosis' 204/600 (seed 0: 67/200, before the 2026-10-08 MC port); with the centroid screen and `BatchedScreen` switched back the harness reproduces the 2026-10-08 C++-faithful rerun (40/200, 3,962.95 local and 44,479.3 global evaluations) exactly, and the centroid default (2026-10-12) moves 40 to 35 (19 fixed, 14 broken, 21 both wrong, about the size of seed-to-seed flips). Under the ray rule the levels that use it run the sequential quick MC and the later levels run in lockstep (mixed mode). Open: `CppEdgeRules` in the batched cost is done in S1 Step 4 (branch `feature/python-speed-parity-s1-loops`, not merged); quiet-machine repeat; `cuda` untested. The default stays off until a later large-scale test (D-20261014-I2)
| S4 | level-0 screen from a per-sample table: the level-0 orientations are the same for every voxel and a spot moves rigidly with the voxel (`docs/spot_translation_derivation.md`) | measured: the level-0 screen of every voxel of the Phase D layer takes 6.58 s on 10 threads (29.7 s on one thread, projected from a 2,000-voxel patch) against 3.94 single-worker hours with `BatchedScreen` (projected from 10 voxels); the per-voxel search after it does not shrink; the hybrid saves about 10% of a Phase D run | **done as a study (merged 2026-10-21, D-20261021-I3); library path parked (D-20261016-I22)** (Steps 0-5, 2026-10-15; picked D-20261014-I19, sized down I20, approved I21; branch `feature/python-speed-parity-s4-translation-screen`; `docs/todo_translation_screen_prototype.md`, MIGRATION_HISTORY "S4: translation correctness and loop ordering"). (1) The affine law is exact in continuous coordinates (1.8e-12 px in float64), integer shifts and phase re-rasterisation are exact on pixel sets (detector edges and float32 flips within 2.3e-4 px listed), a pure pixel-set translation by the rounded shift is not (26.8% of spots equal; it changes 29.4% of re-costs). (2) Best ordering O2 (orientation-major, voxel inner loop, CPU threads): the level-0 screen of all 24,570 voxels takes 6.58 s on 10 threads against 3.94 single-worker hours with `BatchedScreen`; O1 7.9x, O4 (pixel voting, built) 3.3x and MPS 3.4x slower than O2. Under the Python default qualified-peak rule the screen passes 78.5% of the pairs (it barely filters); under `CppOverlapRule` (the C++ behaviour) 26.5%, not timed and not projected. The search after the level-0 screen and re-cost (including the level-0 quick MC) dominates, so the hybrid (seeds take a global level-0 list and re-cost it from a stored layer-wide local-cost table, 20.1 s to build, about 0.62 GB) saves about 10% of a Phase D run; without the table, re-costing each seed alone is slower than today; script-only |
| S5 | a BFS neighbour's first evaluation (inherited orientation) as a translation of the parent's spots | small | idea |
| S6 | repeat the S0 measurements on the full-sample Phase D images (on Newt Agent's machine) | confirms or corrects the ThreeVoxels figures | needs the Phase D data |

# Track H: housekeeping (verified on 0c187a4)

| Item | Where | Fix |
|---|---|---|
| H1 | `icenine/sample.py:425` calls `self.mic_file.save()`; `MicFile` has only `write()`, so `save_sample` always returns False | `.save` to `.write`, plus a test |
| H2 | `icenine/geometry.py:293-304` cites the deleted `DETECTOR_IMPLEMENTATION_PLAN.md`; `matrix_to_euler_torch` raises `NotImplementedError` | implement (gimbal-lock handling) or point the docstring at MIGRATION_HISTORY |
| H3 | `cpp_harness/` has no `Makefile` (gitignored) but tests say "`make ...`"; the C++ ground-truth JSONs in `cpp_outputs/` are gitignored, so the tests that read them skip | commit the Makefile (copy on the archive branch), document regenerating the JSONs |
| H4 | `tests/test_mic_file.py:588` skips: `neighbors_cpp.txt` was deleted (b6f1a90) | regenerate with `test_neighbors.cpp` or drop the test |
| H5 | `tests/validate_cpp_detector.cpp`, `tests/run_validation.sh` are orphaned (superseded by the pixel-exact forward-simulation match) | delete |
| H6 | `CLAUDE.md:119` and `CPP_MODERNIZATION.md` list MicIO/IteratorAdapter warnings that were fixed in 485bf6c; `CPP_MODERNIZATION.md:26` names Boost MPI and serialization, which are not used | correct both |
| H7 | MIGRATION_HISTORY lacks: the `.mic` row layouts (9, 12, 19 columns; deformation order `m00 m11 m22 m01 m12 m02`; side = initial / 2^generation; square-grid variant), the omega-range test gotchas (expected values from C++, JSON round trip of angles), C++ precision figures (Q to 1e-10, omegas to 1.6e-7 rad) | add a short "Formats and validation notes" section (texts on the archive branch) |
| H8 | "Not Yet Ported" omits the quadtree peak search (peak extraction from images) | add the line |
| H9 | `icenine/sampling.py:338` `get_random_local_grid` uses `default_rng(seed=None)` (fresh entropy), so its "preserves C++ determinism" docstring is wrong; the function is unused | fix the docstring or remove the function |
| H10 | `_parse_*` helpers of `mic_file.py` have no direct tests | low priority |
| H11 | root `README.md` is from 2014 (Boost 1.43, gcc 4.2, XDM++ "obtained separately") | refresh |
| H12 | `scripts/dev/`, `.claude/skills/` | workflow tooling, shared by both agents: `log_decision.sh`, `roadmap_row.py`, `verify_handback.sh`, `read_on_fetch.sh`, `make_brief.sh` and the `/decide`, `/implement-plan`, `/apply-fixes` skills (`docs/todo_ilium_workflow_tooling.md`); **scripts merged; skills and owner text waiting for the owner** (approved D-20261018-I2, row D-20261018-I3; Ilium Agent, task branches off `feature/python-speed-parity`) |

# Track Q: search accuracy open questions

| Item | Question | Status | Detail |
|---|---|---|---|
| Q-F1 | CSL-relative fixes F1 / F1b by default? | off by owner decision 2026-10-06, kept off 2026-10-13; opt-in | `todo_findoptimal_wrong_candidates.md`, `findoptimal_robustness_report.md` |
| Q2 | why CMA accepts more wrong-grain fits; tighten the 0.9 relative test for CMA? | open | finisher findings §7.2 |
| Q3 | `CppOverlapRule` default: does the BFS seed test (hit ratio 0.8) and the reseed pilot change with the rule? | open; default stays off (owner 2026-10-13) | finisher findings §7.3 |
| Q4 | which C++ MC difference (block halving, restart base, fixed step) gives the faithful MC's gain | open (ablation untested) | finisher findings §7.4 |
| Q5 | batch the coarse stage | **done**: screen (Task A, `BatchedScreen` on by default); quick MC (S3, opt-in `BatchedQuickMC`, off) | finisher findings §7.5 |
| Q7 | does multi-client C++ seed grains the single-client BFS loses? | open | finisher findings §7.7 |
| Q9 | ray rule pixel counting in cases b and e | deferred by owner 2026-10-12 | finisher findings §7.9 |
| Q10 | angle between observed and predicted ray as an optimisation signal | idea | finisher findings §7.10 |
| Q11 | H0 wrong-basin failures; a real fallback trigger for the hybrids | open | orientation-search findings §6 |
| Q12 | CMA-ES in BFS at scale (1,001 evaluations per neighbour; cross-boundary neighbours) | untested | orientation-search findings §6 |
| Q13 | coarse cost proxy (Idea 1) | evaluated: accuracy-only gain, +5% time | `todo_coarse_cost_proxy_and_dedup.md` |
| Q14 | noisy images are easier (Q8 Step 7: more candidates, wider cost): from the detector noise or the neighbour/twin spots? Add `noise_nodrop` (noise only) to the Step 7 analysis | open (owner 2026-10-14, D-20261014-I13); analysis only, about 10 min | MIGRATION_HISTORY "Q8"; `scripts/phase_d/q8_step7.py` |

# Track X: scale and realism (deferred by the owner)

| Item | What | Status |
|---|---|---|
| X1 | full-sample BFS end to end on the 500-grain sample (classic vs NN and hybrid finisher) | data and library ready; long runs deferred 2026-10-09, again 2026-10-13 ("later") (`todo_500grain_bfs_end_to_end.md`) |
| X2 | reseed on other regions, seed orders and full-sample scale | deferred 2026-10-09 |
| Q8 | noise model of the realistic images (unphysical 10% spot drop) | owner 2026-10-13: remove the random spot drop (keep jitter, hot pixels, blobs), then re-render and re-run the realistic comparisons; Ilium Agent (owner 2026-10-13); **Steps 1-8 done 2026-10-14** on `feature/python-speed-parity-q8-nodrop` (not merged; `docs/todo_q8_noise_nodrop.md`, MIGRATION_HISTORY "Q8"): no-drop variants for Phase D (`render_full.py --noise nodrop`, `full/realistic_nodrop` rendered at Q-max 8 only) and E0 (`RealismConfig` `all_nodrop`, label `realistic_nodrop`); E0 200-voxel re-runs of the Q3 and ray-overlap rows on clean, old realistic and no-drop images (the lenient rule's lead on the old realistic images is gone on the no-drop ones; `min_sin_eta` 0 mean evaluations per voxel on the no-drop images: lenient 51,901, ray 2-pixel window 51,652, C++ rule 49,281). Owner decisions of 2026-10-14 (D-20261014-I5 to I8): **no drop is now the default noise** (Step 6: `NoiseParams()`, `RealismConfig()`, `render_full.py`, E0 harness default images `clean` and `realistic_nodrop`; the named variants `toy`, `all`, `noise` keep the 10% drop; the six `*_realistic*.config` Phase D arms still read the drop images); studies stay on clean images, so the rule "clean images only" is not lifted; no Phase D runs on noisy images for now. Step 7 (are the noisy images easier): the no-drop images give more level-0 candidates in 200 of 200 voxels for the lenient and ray rules, 199 of 200 for the C++ rule (lenient median 182 against 289) and a wider cost function (level-0 halfway width 0.235 to 0.310 degrees), but the closest candidate is at the same distance in 174 of 200 voxels (lenient) and the link to the fixed voxels is weak; cause untested (neighbour layer against noise not separated). Step 8: the four lenient wrong counts that differed from the stored ones are the vertex-0 against centroid screen start (arm `offv0` reproduces all four). Open: the Phase D pilot and reseed on no-drop images (skipped for now, X1); whether any other study script (own `VARIANTS = ["clean", "all"]`) should move to no-drop |
| X3 | multiple detectors with many grains and voxels far from the axis | deferred (`todo_multidetector_many_grains.md`) |
| X4 | intensity inputs and distribution losses for the orientation network | deferred (`todo_intensity_and_distribution_losses.md`) |
| X5 | not yet ported: MPI server/client, step-size file reading, unique-reflection filtering, quadtree peak search | open |

# Track R: research ideas (not started)

- Single-voxel network reconstruction with no start, active imaging, diffraction as Fourier-space
  tomography (`todo_future_ideas_nn_active_fourier.md`).
- NeRF-style refinement (`research_ideas_nerf.md`, parked).
- Hybrid network + optimiser refinement: done 2026-10-06 (`todo_hybrid_nn_refinement.md`); open
  remainder is denser realism, which waits for Q8.

# Track C: C++ code base (low priority)

| Item | What |
|---|---|
| C1 | `Sampling.cpp:468` builds an unseeded local `std::mt19937` per call (same sequence every call, inherited from the Boost code): seed it or document it |
| C2 | C++17: only a dead `boost::any` include would go; the remaining Boost load (ublas, graph, multi_array, circular_buffer) needs Eigen or STL work, not C++17. Shrink `CPP17_MIGRATION_GUIDE.md` to that note |
| C3 | ublas to Eigen (optional) |
| C4 | CMakePresets, CTest, install targets |

# Decisions waiting on the owner

None open. Decisions made are logged in `docs/DECISIONS.md` (rule 6); the "Decided" lines that were
here moved there on 2026-10-14.

# Closed notes (kept for history)

- `todo_hybrid_nn_refinement.md`: done 2026-10-06.
- `todo_findoptimal_wrong_candidates.md`: investigated 2026-10-06; fixes opt-in (Q-F1).
- `todo_coarse_cost_proxy_and_dedup.md`: Idea 1 evaluated (Q13); Idea 2 moved to Task C.
- Finisher findings §7.1 (reseed, answered on the pilot), §7.3 first part (answered), §7.6 (closed by
  the owner: no extra seeds).
