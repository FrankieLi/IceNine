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
| `feature/python-speed-parity` | b1eb09f | Track S S0-S4, Q8, tooling (H12) | S2 waits for your full review |

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

The question was why one Python evaluation cost about 390-500 µs against C++'s 42-55 µs (Task A's
figures). Over 95% of it was Python and torch overhead. Findings: MIGRATION_HISTORY "Python speed
parity" and the S1 to S4 sections.

1. **S1, loops out of `evaluate`** (merged): 1.31x faster per call, bit-identical results. The timings
   are contended; the quiet-machine repeat is still open (D-20261014-I16).
2. **S2, a fused C evaluation behind the opt-in key `FusedCost`** (merged, waits for your review):
   235.7 -> 15.9 µs per call (14.9x). A no-start search is 3.14x faster on ThreeVoxels voxel 0 and
   2.35x on 10 E0 voxels. On those, the batched screen (not covered by S2) is now 5.16 of 6.84 s.
   E0: 35 of 200 wrong with and without the key (0 fixed, 0 broken). Not bit-equal to torch: 93 of
   4,000 orientations differ, every one at a pixel edge (within 0.0014 px). Key off is bit-identical.
   Built and checked with macOS clang only.
3. **S3, batched quick MC (`BatchedQuickMC`, opt-in)** (merged): the quick MC is 4.5-4.7x cheaper per
   evaluation and the search 1.7-1.9x faster (contended). It changes the random draws, so results
   differ from key off (E0: 107 of 600 wrong against 111). Stays off until a large-scale test
   (D-20261014-I2).
4. **S4, level-0 screen from one table per sample** (merged as a study, library path parked): the
   level-0 screen of the whole layer takes 6.58 s on 10 threads against 3.94 hours on one worker,
   but the search after it dominates, so a run saves only about 10%.
5. **Q8, noise model** (merged): images no longer drop 10% of spots at random by default. Studies
   still use clean images (D-20261014-I5).
6. **H12, workflow tooling** (merged): the brief, hand-back, decision and fetch scripts both agents
   use.

## Waiting on you

One list for both tracks. Each item names its track, the question, the agent's proposal if any, and
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
    (E0 20 against 20 µs). Keep it in the recommended configuration?
12. **`.claude/` tooling items Ilium cannot commit** (D-20261028-I5): should the `/implement-plan`
    and `/apply-fixes` skills show `make_brief.sh`'s stderr notes, and should the
    `run-tests-on-change.sh` hook drop its extra `-q`?
13. **S1 quiet-machine timing repeat** (D-20261014-I16), still open: the machine was never quiet.

Both tracks:

14. **The two tracks' best keys have never run together.** Track E kept `FusedCost` off in every
    accuracy arm, and Track S never ran `SmoothCostRank`. `FusedCost` cuts the strict evaluations
    to about 1/15 of their cost. The smooth evaluations are a separate code path, not covered by it.
    So with both keys, the smooth key's share of the time would grow. A short E0 and BFS-pilot run
    with both keys would show the real combined time. This is a proposal, not planned.
15. **S6 needs Newt's machine.** Repeating Track S's measurements on the full-sample Phase D images
    needs those images, which are on Newt's machine. Run it here, or copy the images?
16. **`make_brief.sh` bug** (Newt's request of 2026-10-29, not yet answered by Ilium): the brief
    script refuses an approved plan when a later decision entry cites the approval and says
    "Supersedes: none." Newt worked around it for the review-fix brief by running the script on a
    copy of DECISIONS.md without those lines.

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
| S | MIGRATION_HISTORY "S1" to "S4", "Q8", "S2 follow-ups" | each Track S item's record |
| S | `docs/todo_s2_fused_c_eval.md` | the S2 plan |
| both | `docs/ROADMAP.md` | one row per item, and the two agents' mailboxes |
| both | `docs/DECISIONS.md` | every owner decision (N = Newt, I = Ilium) |

Review reports are gitignored and stay on each machine. Newt's are in `.claude/reports/` here; the
final Track E review is `20261029-1123-review-fable-G.md`, and the fix report is in the fix
branch's worktree. Ilium's (for example `20261028-2100-review-s2-results.md`) are on Ilium's machine.
