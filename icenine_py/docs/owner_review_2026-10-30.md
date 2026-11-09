# Owner review: both tracks, state on 2026-10-30

One page to review the work of both agents so far. Newt Agent runs Track E (fewer evaluations) and
Ilium Agent runs Track S (cheaper evaluations). Each track's detail stays in its own docs, which are
linked at the end. The Track S sections were compiled by Newt from Ilium's branch at b1eb09f.
Ilium has been asked to check them (ROADMAP mailbox, 2026-10-30).

## Where the work is

`feature/eval-reduction` holds everything from both tracks: all of Ilium's branch (through b1eb09f)
is merged into it. Neither branch is on `develop` yet. The agreed merge order is Track E's branch first.

Reading the table: one row per branch, with its newest commit and what it holds.

| branch | head | holds | state |
|---|---|---|---|
| `develop` | 0226857 | finisher and MC study (PR #40) | base of both tracks |
| `feature/eval-reduction` | 2993462 | Track E Tasks A-G, plus all of Track S through b1eb09f | final review done (APPROVE, 7 wording fixes) |
| `feature/eval-reduction-review-fixes` | f52142b | the 7 review fixes (docs, tables, one summary script) | not merged; needs `/review-task`, then merge |
| `feature/python-speed-parity` | see state | Track S S0-S4, Q8, tooling (H12) and its S2 follow-ups; since 2026-10-31 also all of `feature/eval-reduction` through 78fccc5, the `make_brief.sh` supersede fix and Ilium's check of this page | S2 waits for your full review |

## What Track E found (evaluation reduction)

The question was how to make fewer cost evaluations at the same or better accuracy. Tasks A to F each
tested one idea. Task G put the ideas that passed together and measured the total. Plain-English
summary: `docs/findings_eval_reduction_2026-10.md`.

Terms. **E0** is a set of 200 hard test voxels, each searched from scratch. A **full-sample voxel**
is a voxel of the 500-grain sample. A **BFS run** grows the map from seed voxels into their
neighbours. **Wrong** means more than 1 degree from the truth. **Strict-equivalent** converts the
slower smooth evaluations into the time of an ordinary ("strict") evaluation.

1. **No idea cut the total evaluations at the same accuracy.** The ideas that save evaluations lose
   accuracy on E0. The C++ overlap rule: 332,046 -> 228,146 evaluations per full-sample voxel, but
   41 -> 56 wrong on E0. Candidate de-duplication: 75.8% of the evaluations, but 41 -> 55 wrong.
2. **The one idea kept, a smooth ranking key, buys accuracy for time.** `SmoothCostRank E3+E4,tau=0.25`
   fixes 34 of the 41 wrong E0 voxels and breaks none. It costs about +20% time per full-sample voxel.
3. **On the 2,000-voxel BFS pilot region the key is slightly better, at a cost not yet settled.**
   34 wrong against 37 with the key off (14 fixed, 11 broken). +12.4% strict-equivalent over the run,
   but measured per full search the run took 35.4% longer (contended, cause untested).
4. **Every wrong BFS voxel is a boundary voxel carrying a touching grain's orientation.** None comes
   from a wrong full search when the key is on.
5. **The outward search for BFS neighbours does not pay.** It stays off (D-20261029-N2).
6. **A full-sample BFS run would take about a day on one process** (22-27 h, linear scaling, untested).

Items 3 and 6 use the corrected wording of the review-fix branch, which is not merged yet.

## What Track S found (cost per evaluation)

The question was why one Python evaluation cost 396-533 µs on the full images against C++'s
49-51 µs (Task A's figures). The evaluation count was already at parity. About 210 µs of each call
was fixed overhead, and the pixel overlap itself took 8.3 of 317 µs (S0, ThreeVoxels voxel 0): the
rest was Python and torch dispatch. Findings: MIGRATION_HISTORY "Python speed parity" and the S1 to
S4 sections.

1. **S1, loops out of `evaluate`** (merged): 1.31x faster per call on the Q-max 8 re-cost and final
   arms (1.14x on the screen arm), bit-identical results. The timings are contended; the
   quiet-machine repeat is still open (D-20261014-I16).
2. **S2, a fused C evaluation behind the opt-in key `FusedCost`** (merged, waits for your review):
   235.7 -> 15.9 µs per call (14.8x), single-worker on a quiet machine. A no-start search is 3.14x
   faster on ThreeVoxels voxel 0 and 2.35x on 10 E0 voxels. On those, the batched screen (not
   covered by S2) is now 5.16 of 6.84 s. E0: 35 of 200 wrong with and without the key (0 fixed,
   0 broken). Not bit-equal to torch: on S1's parity set, 16-29 of 1,000 orientations differ per
   triangle mode and none in the screen modes (464 over 36 dataset-modes), every one a flip at a
   pixel edge (within 0.0014 px). Against C++ on Task A it agrees as often as torch (2,229 and 2,230
   of 2,230). Key off is bit-identical. Built and checked with macOS clang only.
3. **S3, batched quick MC (`BatchedQuickMC`, opt-in)** (merged): the quick MC is 4.5-4.7x cheaper per
   evaluation and the search 1.7-1.9x faster (contended). It changes the random draws, so results
   differ from key off (E0: 107 of 600 wrong against 111). Stays off until a large-scale test
   (D-20261014-I2).
4. **S4, level-0 screen from one table per sample** (merged as a study, library path parked): the
   level-0 screen of the whole layer takes 6.58 s on 10 threads against 3.94 hours on one worker
   (projected from 10 voxels), but the search after it dominates, so a run saves only about 10%.
5. **Q8, noise model** (merged): images no longer drop 10% of spots at random by default. Studies
   still use clean images (D-20261014-I5 to I8).
6. **H12, workflow tooling** (merged): the brief, hand-back, decision and fetch scripts both agents
   use, with the S2 follow-ups of 2026-10-29 (D-20261028-I3 to I5) and the `make_brief.sh`
   supersede fix of 2026-10-31 (item 16).

## Decisions of 2026-11-09

On 2026-11-09 you asked Newt to decide the 24 items below (D-20261109-N1). The decisions are
recorded as D-20261109-N2 (items 1-8), N3 (items 14-16) and N4 (items 9-13, 17-24). Four things
are still yours; they come first.

Still with you:

- **Item 7, the merge to `develop`.** Newt opens the PR; the merge waits for your go.
- **Items 12, 22 and 23.** They change what the agents themselves may do, so an agent should not
  decide them. Newt recommends yes to all three (for item 22: only when the `.claude/` files are
  identical to the files already committed on the merged branch).
- **Approval of Ilium's next plans**: the screen port (item 21) and the `BatchedQuickMC` large-scale
  test (D-20261101-I1).
- **shmem**: `shmem ask` fails on this machine (the API key in `~/.shmem/.env` is refused, and the
  local embedder is not running). `shmem query` and `shmem add` work.

Reading the table: one row per item, with the decision and who acts next.

| item | decision | who acts |
|---|---|---|
| 1 G3 budget | The 3.2 h run is accepted after the fact. From now on a run that passes its stop line is stopped, or gets a decision entry first. | none |
| 2 pilot baseline | The key-off rerun on today's code is the baseline. | none |
| 3 timing A/B | Run it as part of item 14: `FusedCost` on in both arms, smooth key on against off, one worker each. | Newt |
| 4 `SmoothCostRank` default | Stays opt-in, off by default. No golden moves. | none |
| 5 full-sample BFS run | Run once after item 14, `FusedCost` on, smooth key off. It is the baseline of the autoresearch stage and the key-off arm of the `BatchedQuickMC` test. | Newt |
| 6 review fixes | `/review-task`, then merge. | Newt |
| 7 PR to `develop` | Newt opens it after item 6. | you (merge) |
| 8 `ScreenEvaluator` follow-up | Folded into the screen port (item 21). | Ilium |
| 9 S2 full review | S2 is accepted on its record. Your own reading stays open; nothing waits for it. | none |
| 10 GPU or fused CPU | The fused CPU path takes priority. No GPU work on the strict cost until CUDA hardware is in use. | none |
| 11 `BatchedQuickMC` with `FusedCost` | Not in the recommended configuration when `FusedCost` is on. Its large-scale test stays; run its accuracy arms with `FusedCost` on. | Ilium (plan) |
| 12 `.claude/` edits | Not decided by Newt. Recommendation: yes to both. | you |
| 13 quiet-machine repeats | Rerun the S1 and S3 timings single-worker under the preflight, and relabel. | Ilium |
| 14 both keys together | Run on E0 and the BFS pilot; compare wrong counts; time single-worker. | Newt |
| 15 S6 | Ilium runs it on Ilium's machine, after the screen port. | Ilium |
| 16 withdrawal check | Left as it is until it refuses a plan that is in use. | none |
| 17 GCC and Linux | Validate before the first use on Linux. Nothing to run now; the key should warn on an unvalidated compiler. | Ilium |
| 18 does C++ fuse | Newt checks on Newt's machine, which has the C++ build. | Newt |
| 19 C++ toolchain for Ilium | Not installed. C++ checks go to Newt through the mailbox. | none |
| 20 `FusedCost` default | Off by default. Every timing baseline of new work uses it on. | both |
| 21 screen port | Next Track S speed item. Design it layer-wide: one orientation across all voxels at a time. | Ilium (plan), you (approval) |
| 22 hook on merges | Not decided by Newt. Recommendation: yes, for identical files only. | you |
| 23 refused commands | Not decided by Newt. Recommendation: confirm Ilium's reading. | you |
| 24 Q8 leftovers | No reseed. The no-drop pilot is dropped. Studies stay on clean images. | none |

What comes next is the autoresearch stage (Track A): plan `docs/autoresearch/PLAN.md`, saved state
`docs/autoresearch/STATE.md`.

## Waiting on you

The list as it stood before the decisions above. One list for both tracks. Each item names its track, the question, the agent's proposal if any, and
where to read.

Track E:

1. **G3 budget (review fix 5).** Did D-20261029-N3 cover the G3 arm's 3.2 h run, when the plan
   budgeted about 1 h (stop line about 2 h)? The G3 RUN_LOG now says the "3 h / 8 h" figures had no
   source. Read: review report, "Decisions for the owner".
2. **Baseline of the BFS pilot.** Accept the key-off rerun on today's code (D-20261029-N2) in place
   of the plan's stored pilot? The stored pilot does not reproduce on today's code. Both comparisons
   are in the findings doc.
3. **Timing A/B.** Run the BFS pilot with and without the smooth key at the same time, or on one
   worker each, before quoting the key's BFS time cost? Today's +12.4% is an estimate, not a measurement.
4. **`SmoothCostRank` default.** Proposal: keep it opt-in and turn it on for runs where full
   searches are often wrong. Read: findings doc, "Decisions for the owner" item 1 (it lists the
   goldens that would move).
5. **Full-sample BFS run** (also Track X item X1, deferred twice). About a day on one process.
6. **Merge the review fixes.** Next step: `/review-task` on `feature/eval-reduction-review-fixes`,
   then merge it.
7. **PR `feature/eval-reduction` -> `develop`** (`/finish-feature`), after item 6.
8. **`ScreenEvaluator` on `cost_fn.frame_table()`.** A small Track E follow-up of S1, waiting to be
   scheduled (noted 2026-10-14).

Track S:

9. **S2 full review** (D-20261027-I1). Read: MIGRATION_HISTORY "S2: fused C evaluation, Steps 1-6
   results"; plan `docs/todo_s2_fused_c_eval.md`.
10. **GPU with `FusedCost` on** (D-20261028-I2). A GPU `SearchDevice` still takes batches of 64 or
    more, and on this Mac (MPS) that is 2.4-6.8x slower than the fused CPU path. Should the fused
    path take priority? CUDA is not measured.
11. **`BatchedQuickMC` with `FusedCost` on** (D-20261028-I2). It gives no per-evaluation gain
    (quick MC: E0 20 against 20 µs; ThreeVoxels 125 against 114 µs). Keep it in the recommended
    configuration?
12. **`.claude/` tooling items Ilium cannot commit** (D-20261028-I5): should the `/implement-plan`
    and `/apply-fixes` skills show `make_brief.sh`'s stderr notes, and should the
    `run-tests-on-change.sh` hook drop its extra `-q`?
13. **S1 and S3 quiet-machine timing repeats** (D-20261014-I16), still open: their timings are
    labelled contended. The preflight did pass on every S2 Step 5 timing (2026-10-28), so the repeat
    can be run now. Ilium's proposal: rerun S1's and S3's timing scripts single-worker under the
    preflight and relabel; no plan change. Run it?

Both tracks:

14. **The two tracks' best keys have never run together.** Track E kept `FusedCost` off in every
    accuracy arm, and Track S never ran `SmoothCostRank`. `FusedCost` cuts the strict evaluations
    to about 1/15 of their cost. The smooth evaluations are a separate code path, not covered by it.
    So with both keys, the smooth key's share of the time would grow. A short E0 and BFS-pilot run
    with both keys would show the real combined time. This is a proposal, not planned.
    Ilium (2026-10-31): no reason it would not work. `SmoothCost` is its own float64 torch code and
    never calls `VoxelCostFunction.evaluate`; under `SmoothCostRank` the strict cost on each
    candidate is the ordinary `evaluate`, which `FusedCost` covers for the integrated
    configuration's settings (default rule, hard mode, pixel radius 0 and 3, eta limit; no ray rule,
    no soft mode). Three things to expect: the results are not bit-equal (pixel-edge flips), so the
    run compares wrong counts against G2's E0 (41 -> 7) and the pilot, not records; the speedup is
    capped by the two parts `FusedCost` does not cover, the batched screen and the smooth
    evaluations; and the timing arms need a quiet machine, single-worker.
15. **S6 needs Newt's machine.** Repeating Track S's measurements on the full-sample Phase D images
    needs those images, which are on Newt's machine. Run it here, or copy the images?
16. **`make_brief.sh` bug** (Newt's request of 2026-10-29): the brief script refused an approved
    plan when a later decision entry cites the approval and says "Supersedes: none." Newt worked
    around it for the review-fix brief by running the script on a copy of DECISIONS.md without those
    lines. Ilium fixed it on 2026-10-31 (merged as 7681f02 after an Opus review; D-20261031-I1):
    only the approval's whole ID in a later entry's own `Supersedes:` field refuses now, and a field
    that begins with "none" supersedes nothing. Left for you (D-20261031-I2): the separate
    withdrawal check is a keyword test on any later entry that names the plan, and I1's own wording
    now trips it for the finished S3 plan (no brief needed). Narrow that check?

Track S, added by Ilium (2026-10-31; these were only in Ilium's reports):

17. **`FusedCost` on GCC and Linux.** Validated with macOS clang (arm64) only; the GCC pragma that
    keeps multiply-adds unfused in `_fused_eval.c` has never been compiled (D-20261027-I10).
    Validate on Linux before anyone turns the key on there?
18. **Does C++'s clang build fuse the clip's multiply-adds?** The Python key-off `_rasterize.c`
    is compiled with fusing on, and switching it off changed 106 of about 400,000 border
    triangles (D-20261027-I9). Whether C++ fuses the same expression (`SutherlandHodgman.h:181`)
    decides which of the two matches C++ bit for bit. It needs the C++ build to check (item 19).
19. **Install the C++ toolchain** (cmake, boost, eigen, open-mpi) on Ilium's machine? Task A's
    inputs for S2's C++ comparison were accepted on aggregates because the C++ driver cannot run
    here (D-20261027-I6). It is a system change, so it waits for you.
20. **`FusedCost` default.** It stays opt-in and off; no change is proposed (D-20261012-I1). Items
    10, 11, 14 and 17 bear on it.
21. **The batched screen is now the largest cost with `FusedCost` on** (5.16 of 6.84 s on 10 E0
    voxels). Porting the screen to the fused C path would be the next Track S item if you want more
    speed; it needs a plan and your pick. S5 (a BFS neighbour's first evaluation as a translation
    of the parent's spots) is still an idea.
22. **Pre-commit hook on merges** (D-20261028-N3): the hook refused your approved `.claude/` files
    and long lines when merging across branches. The long lines are fixed (2c8b487). The
    `.claude/` refusal would recur whenever a merge brings `.claude/` changes. Should
    `precommit_check.py` let a merge commit bring `.claude/` files that are already committed on the
    other branch? That file belongs to the shared tooling, so either agent can change it once you
    say.
23. **Refused commands in agent runs.** In the tooling fix rounds the worktree guard refused some
    combined shell commands; the agents reran them as plain commands or with the Edit tool, as the
    refusal messages suggested. Ilium reads that as allowed (the briefs forbid redoing a refused
    action another way, such as writing files through Bash when Edit is refused). Confirm?
24. **Q8 leftovers.** The Phase D pilot and reseed on no-drop images (skipped, X1), and whether the
    study scripts with their own `VARIANTS = ["clean", "all"]` should move to the no-drop images.

## Where to read more

Reading the table: one row per document, with the track and what it covers.

| track | document | covers |
|---|---|---|
| E | `docs/findings_eval_reduction_2026-10.md` | the total (Task G), each idea's result, proposals |
| E | `docs/findings_cost_edges_2026-10.md` | Task B: where the cost function's edges come from |
| E | `docs/findings_spot_pairing_2026-10.md` | Task D: pairing spots across the two detectors |
| E | `docs/findings_cost_smoothing_2026-10.md` | Task E: the smooth costs |
| E | `docs/findings_outward_search_2026-10.md` | Task F: searching outward from a start |
| E | `benchmarks/eval_reduction/G/maps/` | BFS pilot maps (also an Artifact page, link in Newt's checkpoint) |
| S | MIGRATION_HISTORY "Python speed parity" (S0), "S1" to "S4", "Q8", "S2 follow-ups" | each Track S item's record |
| S | MIGRATION_HISTORY "S2: fused C evaluation, Step 0" and "Steps 1-6 results" | S2's results, parity and timings |
| S | `docs/todo_s2_fused_c_eval.md` | the S2 plan |
| S | `docs/DECISIONS.md` D-20261027-I1 to D-20261028-I5 | the decisions Ilium made under your S2 delegation |
| S | Ilium's review page "S2 Fused Cost Review" (private Artifact; link in Ilium's session) | S2 and the tooling follow-ups on one page |
| both | `docs/ROADMAP.md` | one row per item, and the two agents' mailboxes |
| both | `docs/DECISIONS.md` | every owner decision (N = Newt, I = Ilium) |

Review reports are gitignored and stay on each machine. Newt's are in `.claude/reports/` here; the
final Track E review is `20261029-1123-review-fable-G.md`, and the fix report is in the fix
branch's worktree. Ilium's (for example `20261028-2100-review-s2-results.md`) are on Ilium's machine.
