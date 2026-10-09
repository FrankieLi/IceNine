---
title: "Finisher and MC study, and the first full-sample BFS test: findings of 2026-10-07 to 2026-10-09"
subtitle: "Cost sensitivity, MC behaviour, two port discrepancies against C++, CMA-ES, and a limit of the BFS seed rule seen in the first multi-grain test"
date: "2026-10-09"
geometry: margin=1in
fontsize: 11pt
---

# 1. Purpose and setup

The earlier report, `docs/findings_orientation_search_2026-10.md`, ended with these open questions:

- Why does FindOptimal (the finisher, `refine_from_candidates`) stop short of the truth?
- Is the cause the sensitivity of the cost function?
- Why does MC run out of restarts, and is MC the right tool?
- Was the April 2026 MC-vs-optimizer comparison adopted too soon?

The owner then asked for an end-to-end test: a BFS reconstruction of the 500-grain sample with new orientations,
comparing the classic optimizer with the new options for timing and accuracy.

This report gathers the answers. It also covers two findings the work turned up on the way:

- two places where the Python port did not follow the C++ algorithm;
- a limit of IceNine's BFS seed rule (a seed is drawn only from unvisited voxels), first seen in a 2,000-voxel
  multi-grain region; whether it appears only at larger scale was not tested.

**Where the details are.** The detailed records live in `icenine_py/MIGRATION_HISTORY.md` (MH), section "Finisher
and MC study (2026-10-07)". The subsections are:

| Phase | MH subsection |
|---|---|
| A | "Phase A results" |
| B1/B2 | "Phase B1/B2 results" |
| B3 | "Phase B3 results" |
| C | "Phase C results" |
| D | "Phase D library prerequisites", "Phase D data", "Phase D seed-cost diagnosis", "C++ vs Python variance stage", "C++-faithful MC: reruns", "Phase D pilot results" |

Result files are under `icenine_py/benchmarks/`. Section 6 lists which earlier numbers are superseded.

**Conventions.**

- "Wrong" means an error over 1° (symmetry-reduced misorientation).
- Cases within a voxel or grain are not independent. Where a p-value is quoted, it comes from a voxel- or
  grain-clustered test.
- Timings are contended unless marked single-worker.

# 2. Short answers

1. **Is the stopping short caused by cost-function sensitivity? No.**
   - The cost still slopes clearly downhill toward the truth well below the finisher's error.
   - It is flat only below about 0.0005–0.001°.
   - Precision near 0.003° is reachable.
2. **Why did MC stop short and run out of restarts? Mostly because the Python port was not C++.**
   - Two code differences, not config differences, made the Python "MC" behave unlike IceNine's C++ algorithm: the
     main MC loop and the VarianceMinimizing restart.
   - Both are now ported faithfully and checked against C++ debug traces.
   - With the C++-faithful code, the classic finisher is better and cheaper than recorded before.
3. **Is MC the right tool? For coarse candidate ranking, yes. For the final polish, CMA-ES is better.**
   - CMA-ES is about 5× more precise at a similar cost.
   - It gives no fewer wrong answers.
   - In BFS it accepts more neighbours with a wrong-grain orientation.
4. **Was the April comparison adopted too soon? It did not test what was deployed.**
   - It scored "success" as ending under 1° from starts 1–5° away. Its successful runs ended about 0.5° from the
     truth.
   - It used MC settings and a Python MC loop that differ from the deployed C++ algorithm.
   - So it neither supports nor contradicts the deployed configuration at the 0.01–0.03° scale.
5. **End-to-end BFS: no verdict yet.**
   - The 2,000-voxel pilot found that the BFS seed rule leaves 25–27 of 60 grains without a full search (a 26th to
     28th lost grain had its one seed rejected). The C++ server uses the same rule. The seed-count oracle had
     predicted that about 28% of grains would be swallowed; the pilot rate is higher.
   - The voxels of those grains end wrong whatever the optimizer, about a third of all voxels.
   - This must be fixed before the classic-vs-new comparison means anything.

# 3. Findings

## (a) The cost function is not the limit (Phase A)

On the 300 T5 cases, Phase A sampled the cost at 400 random directions × 11 radii around the truth, on both clean
and realistic images.

- **Steep.**
  - A 0.01° rotation changes the cost by a median 0.025 (realistic) to 0.037 (clean).
  - Within each case, cost tracks the angle from the truth (Spearman median 0.97).
  - Along individual directions, the cost is non-decreasing from 0.005° outward in 99% or more of rays.
- **Flat only at the bottom.** The region with cost at or below the truth's has a median radius of 0.0005° (q75
  0.001°).
- **The realistic minimum.** Its offset from the truth is not resolved below the 0.0005° sampling grid.
- **Pixel and frame quantisation.** These give a centroid-quantisation scale of about 0.0125° per voxel.
  - Gauss-Newton on clean windows reaches it.
  - It is a typical scale, not a lower bound: searches that use the exact pixel sets go below it.
- **Conclusion.** The finisher's former 0.023° median error was an optimizer limit, not an information limit.
  - Phase A's continuation results (A3, T5) were measured with the old, non-faithful MC and VarianceMinimizing.
    Section 6 has the caveat.

## (b) Two port discrepancies against C++

Both were found by tracing the same voxel through C++ (a scratch debug build) and Python. Both are code differences,
which no config change could fix.

1. **The VarianceMinimizing restart.**
   - C++ restarts a failed subregion run within ±the current (doubled) radius, about the **initial** orientation.
   - Python restarted within ±half the box, about the **current best**.
   - On full-sample images, Python's version stayed on the steep slope, kept the cost variance high and kept
     extending its own budget.

   | Variance-stage evaluations per seed (mean of 6 voxels) | Value |
   |---|---|
   | Python before | 707,732 |
   | Python after | 1,609 |
   | C++ | 1,295 |

2. **The main MC loop** (`MCOptimizer.optimize`, against C++ `RandomRestartZeroTemp`). The quick MC and FindOptimal
   both use it.

   | | C++ | Python port (before) |
   |---|---|---|
   | Step halving | once per fixed block of `nMinErgodicSteps` (31 here), when the block beats the global best | on every improvement |
   | Block length | fixed | recomputed (×8 per halving) |
   | Restart trigger | after a failed block | after N steps in a row without improvement |
   | Restart point | ±tan(box)/√48 about the **initial** orientation | ±box/2 about the **best** (about 3.5× wider) |
   | Stop | more than `SuccessiveRestarts` failures in a row | total restarts |

   After the port, a debug trace shows zero violations of the C++ block rules. The quick MC uses 12 evaluations per
   call in both.

**Remaining gap to C++.**

- A Python no-start search still makes 1.31–1.57× C++'s evaluations: 332k against 237k per seed.
- The excess comes from the coarse stage: C++ passes 1.24–1.53× fewer candidates out of level 0. The cause is
  uninvestigated; C++'s "Num Cliques" step is a candidate.
- Per evaluation, Python is about 10× slower (about 500 µs against 50 µs).

**Tests.**

- `tests/test_variance_stage_parity.py` and `tests/test_mc_cpp_parity.py` pin the C++ semantics.
- Goldens that pinned a random draw were re-recorded, with the old values kept in comments.
- No check against C++ ground truth (`cpp_outputs/`) changed.

## (c) What the faithful MC changes

These results are from reruns on the same cases, images, starts and seeds as before.

**Finisher precision, T5 net-started cases (realistic).**

| | Before | After |
|---|---|---|
| Default finisher, median error | 0.0229° | 0.0168° |
| Default finisher, share under 0.02° | 44% | 56% |
| Default finisher, median evaluations | 2,629 | 560 |
| Deployed MC alone, median error | 0.0382° | 0.0175° |
| Sweep cases (1,000), deployed MC | 0.1609° | 0.0249° |

**Wrong answers, 200-voxel search from scratch (E0).** The gain comes from the quick MC in the coarse stage: the
CMA-ES arm, whose finisher is unchanged so that the coarse-stage quick MC is its only changed component, fell by a
similar amount (66 to 40 clean, 44 to 24 realistic).

| Images | Old MC (seeds 0–2) | Faithful MC (seed 0) |
|---|---|---|
| Clean | 66–71 wrong of 200 | 40 |
| Realistic | 45–52 wrong of 200 | 25 |

**The B1 mechanism described the old loop.**

- The step no longer collapses to 0.0003°.
- 87% of improving runs restart after their first improvement (before: 0%).

**The quarter-box VarianceMinimizing continuation is worse than recorded**, especially from far starts:

- T5 H0 at 10,000 evaluations: 0.0096° became 0.1378°.
- Wrong voxels: 20 became 29 of 100.
- 1 voxel improved and 13 got worse, of 16 (p = 0.0018).

**Unchanged.** BFS neighbour refinement under `mc` keeps the inherited start: 47/50 cases return the start
unchanged. Its 5° box is far too coarse for starts about 0.07° off.

## (d) Local finishers head to head (B3, faithful MC)

Equal evaluation budgets, T5 H3 realistic:

| Finisher | Median error | Share under 0.02° | Evaluations |
|---|---|---|---|
| CMA-ES (σ0 0.2°) | 0.0031° | 94% | 250 |
| Nelder–Mead | 0.0031° | 94% | 250 |
| CMA-ES (σ0 0.2°) | 0.0022° | 96% | 1,000 |
| Deployed MC alone | 0.0175° | 55% | 193 |
| Default finisher | 0.0168° | 56% | 560 |

The ratio of median errors, CMA-ES at 250 evaluations against the faithful MC, is:

| Comparison | T5 H3 | Sweep |
|---|---|---|
| Against the deployed MC | 5.5× | 6.5× |
| Against the whole finisher | 5.3× | 4.8× |

Other results from B3:

- **Far starts.** Starts more than 1° off stay mostly wrong for every local method.
  - CMA-ES repairs about half of the 33 T5 starts at 1.5–3°.
  - Nelder–Mead repairs none.
- **Gauss-Newton** is precise on clean windows (about 0.013°) but not on realistic ones.
- **Hybrid Riemannian Adam** does not help.
- **Success-rate step rule.** The (1+1)-ES run in B3, with one untuned setting, shrinks its step too early and
  stalls.

## (e) CMA-ES in the library (Phase C)

- **The switch.** Opt-in `LocalOptimizer cma`; `cma` is now a core dependency.
  - It replaces the finisher's MC and VarianceMinimizing stages.
  - It also replaces BFS neighbour refinement (`CMANeighborMaxEvals`, default 250).
  - The default (`mc`) path is bit-identical to before.
- **Library check.** The library version reproduces B3's CMA-ES orientations exactly (50/50).
- **E0 with the faithful MC.**
  - `cma` and `mc` differ in the wrong count by 1 voxel (p = 1.0).
  - `cma` stays about 4× more precise among right answers: 0.0021° against 0.0087° (realistic).
  - It costs +2–7% evaluations.
  - The Phase C claim "CMA fixes 8 voxels" was an artefact of the old MC.

## (f) Full-sample data and cost (Phase D data, seed-cost diagnosis)

- **Sample.** The 500-grain Example2 sample keeps its geometry, with a new random orientation per grain (seed 0).
  - It has 497 grains and 24,570 voxels.
  - Every new orientation is at least 1.29° from every old one, and at least 3.9° from adjacent grains' new
    orientations.
  - Reconstruction configs read a grid-only copy, so the truth is used only for scoring.
- **Images.** No Example2 images existed before; an earlier claim that they did was wrong. The new ones are Python
  forward renders, both detectors, 180 frames each.

  | Image set | Lit pixels | Render time (10 workers) |
  |---|---|---|
  | Clean, plus realistic (detector noise added) | 7.41 M clean, 6.59 M realistic | 48 s for both |
  | Realistic Q-max 16 | 58.8 M clean (about 8× more) | 246 s |

  - C++ reads the images and gets the same quality on a test voxel (0.94697).
- **Memory.** A shared read-only image loader (`ExperimentalData.from_binary_memmap`) cuts a process from about
  7.7 GB to 0.37 GB, bit-identical to the dense loader. Many runs can now go at once on a 32 GB machine.
- **Seed cost.**
  - A no-start search takes about 174 s single-worker in Python (332k evaluations), (contended), against about 18 s in C++ (contended; 10–21 s).
  - Before the variance-stage fix it took about 9 minutes.
  - Searches on a voxel's isolated image are cheap only because they prune hard: on those images, 3 of 6 seeds
    ended wrong.
  - So per-voxel timings from earlier studies do not carry over to full samples.

## (g) The BFS seed rule (Phase D pilot)

**Setup.**

- 2,000 voxels covering 60 grains.
- Five single-process BFS runs (MC, CMA and CMA without retry on clean images; MC and CMA on realistic Q-max 16).
- The revisit option was on in every run.
- Each run took about 2.1–2.4 h (contended).

| Run | Wrong of 2,000 | Median error of right answers | Unresolved |
|---|---|---|---|
| MC, clean | 740 | 0.0119° | 710 |
| CMA, clean | 742 | 0.0061° | 681 |
| CMA without retry, clean | 734 | 0.0055° | 677 |
| MC, Q-max 16 | 709 | 0.0124° | 691 |
| CMA, Q-max 16 | 617 | 0.0141° | 580 |

**Why so many are wrong.**

- In every run, 26–28 of the 60 grains were lost (no right voxel): 25–27 never got a full search, and in each run
  one more got a full search that was rejected (grain 394 in four runs, grain 148 in one). All five runs start
  from the same seed order in one region, so this is one draw, not five replications (28 of 60 grains, Wilson
  34.6–59.1%, treating grains as independent). Truncation by the region does not explain it: 15 of 33 grains
  wholly inside the disc were lost against 13 of 27 truncated ones.
- The expansions of neighbouring grains reach all their voxels first and reject them (REFIT). A REFIT voxel is
  not given a full search: a seed is drawn only from unvisited voxels.
- The revisit only offers neighbouring grains' orientations.
- So those voxels end a median 41–44° off. They account for 98% or more of the unresolved voxels.
- This is IceNine's own rule. The C++ server alone chooses seeds, and its `Pop` skips every voxel that is not
  unvisited in the server grid (REFIT included). A single-client C++ run is therefore expected to lose the same
  grains. In a multi-client run another client's expansion can refit such a voxel only from that client's
  orientation (what the revisit does); the server might also dispatch a seed into the grain before the neighbours'
  results arrive. Neither effect has been tested with a multi-client C++ run.
- The seed-count oracle (MH, "Phase D seed-cost diagnosis") had predicted this: 497 pieces under the BFS radius
  (not 612) and about 140 of 497 grains (28%) swallowed before any seed is drawn. The pilot's 42–45% (43–47% with
  the rejected seed) is higher; the oracle ignores wrong cross-boundary acceptances, and the pilot is one region.

**What the pilot does measure.**

- **Seeds are all right** (0 wrong).
- **Every accepted wrong voxel carries a neighbouring grain's orientation.**
- **CMA accepts about twice as many of these as MC:** 61 against 32 on clean images, 37 against 19 on Q-max 16
  (grain-clustered p 0.0001 and 0.0007). The cause is untested.
- **The CMA retry is not worth keeping.** 2 of 1,895 retries were accepted, both wrong, for about 60% more
  neighbour evaluations.
- **Seeds dominate the cost:** 84–92% of the wall time, at about 3.5 minutes each (contended).
- **Projected full-sample time** with each of the 497 grain pieces seeded: about 30–34 h per run. Runs can go concurrently at
  about one core each.

# 4. Recommendations

| Topic | Recommendation |
|---|---|
| Port parity | Keep the C++-faithful MC and VarianceMinimizing as the default classic path. |
| BFS | Before any end-to-end comparison, add an opt-in reseed. When no unvisited voxels remain, give one voxel of each connected cluster of still-REFIT voxels a full search, then expand with revisit. This departs from the C++ seed rule, so label it as an option. One reseed may fail (a rejected seed recurred in the same grain). Measure its cost and accuracy on the pilot region. |
| Finisher | Offer CMA-ES (σ0 0.2°, 250–1,000 evaluations) where precision matters. Report its higher rate of wrong-grain acceptances in BFS, and investigate it. |
| CMA retry | Drop the wider-sigma retry (`CMARetrySigma0 0`). |
| Seed cost | The lever is the cost per seed: the level-0 candidate count, and about 10× per evaluation against C++. The neighbour budget is not. |
| Scale | Re-pilot with the reseed (about 3 h). Then choose between several 2,000-voxel regions (hours) and the full sample (about 30–34 h per run with every piece seeded, projected, contended; several seeds multiply the runs). |

# 5. Decisions on record

- Port the C++ algorithm faithfully rather than tune around the port's behaviour. This was the owner's decision on
  2026-10-08, after the variance-stage finding.
- `cma` is a core dependency. CMA-ES is opt-in and off by default.
- BFS revisit (`BFSRevisitRefit`) approximates the multi-client C++ refits of REFIT voxels within one grid (C++
  refits through other clients' independent grids, and also refits fitted voxels, keeping the higher confidence). The
  C++ restart-only Refit path is not used for this.
- Phase D uses truth-free grid configs and the shared image loader.

# 6. Superseded earlier claims

| Earlier claim | Status |
|---|---|
| "MC as deployed stops far above the truth (0.038° / 0.161°)" (B3) | Old Python loop. Now 0.0175° / 0.0249°. |
| B1 mechanism: step collapse to 0.0003°, no restart after the first improvement | Old Python loop. Not C++ behaviour. |
| "Default finisher 0.0229°, about 2,600 evaluations" (T5, B3, earlier report) | Now 0.0168°, 560 evaluations. |
| "Quarter-box VarianceMinimizing closes 93% of the gap, reaching 0.006°" (T5, A3) | Old variance stage. The faithful stage is worse from far starts. |
| "CMA-ES fixes 8 E0 voxels MC got wrong" (Phase C) | Old MC. Now no difference in wrong count. |
| "No-start search about 20 s" (per-voxel images) | About 174 s on full-sample images. Per-voxel pruning is not representative. |
| "Full-sample BFS about 1.4 h per pipeline on 10 workers" | The Python BFS is single-process. About 30–34 h per run with all 497 pieces seeded (as piloted, 18–30 h). |
| "ManyGrains images exist" | False. They were first rendered on 2026-10-07. |

These still stand (they do not depend on the optimizer):

- the Phase A landscape results (A1, A2).

Not re-run; their baselines changed: the hybrid/network and proxy-rerank findings of the earlier report, the HP
sweep, the finisher diagnosis and the FindOptimal-robustness results. They are paired against, or run through,
the old Python quick MC and old variance stage; with the faithful quick MC the E0 wrong count fell from 67 to 40 of
200 (clean) and 52 to 25 (realistic), so their wrong rates are not shown to be independent of the MC.

The earlier report's precision and timing numbers for the default finisher carry the caveat above.

# 7. Open questions

1. **Reseed.** What does the reseed of still-REFIT clusters cost, and how much does it fix?
2. **Wrong-grain acceptances.** Why does CMA accept more of them? Should the 0.9 relative test be tightened for
   CMA fits?
3. **Coarse candidates.** Why does C++ pass fewer level-0 candidates? Closing that gap would cut about 30% of
   Python's seed evaluations.
4. **Which difference matters.** Which of the C++ differences (block halving, restart base, fixed step) produces
   the faithful MC's gain? An ablation is untested.
5. **Per-evaluation cost.** Can the 10× per-evaluation gap be cut, for example with batched evaluation in the
   coarse stage?
6. **Seeds 1–2 for E0.** Run them with the faithful MC, to firm up "no difference in wrong count".
7. **Multi-client C++ reference.** Does C++ with several clients seed any of the grains that the single-client BFS
   loses? By the C++ code only timing (asynchrony) could do it; untested.
