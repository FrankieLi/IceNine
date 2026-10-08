---
title: "Orientation search for NF-HEDM: findings of 2026-10-01 to 2026-10-07"
subtitle: "Learned Gauss-Newton network, FindOptimal, hybrids, a cost proxy, run-time profiling and the follow-ups"
date: "2026-10-07"
geometry: margin=1in
fontsize: 11pt
---

# 1. Purpose and setup

IceNine reconstructs crystal orientations from near-field (NF) HEDM detector images with a forward model. For a
candidate orientation, the simulator predicts where the spots land on the detector frames, and a cost (1 minus the
pixel-overlap quality) scores the match. The Python port searches in three stages: a multi-level discrete search, a
quick Monte Carlo (MC) step per level, and a finisher (`refine_from_candidates`, "FindOptimal", followed by a
VarianceMinimizing pass).

This report collects one week of work on finding and refining a voxel's orientation:

- a learned network;
- the network combined with FindOptimal;
- a cheaper key for pruning candidates;
- where the run time goes;
- why the full search sometimes returns a wrong answer.

**Data.**

- **Samples.** The Cu (FCC) ThreeVoxels and ManyGrains examples. Most single-voxel studies use random voxels of the
  500-grain ManyGrains sample.
- **"Realistic" data.** Exact synthetic detector windows plus neighbour/twin overlap plus detector noise (missing
  spots, edge jitter, spurious blobs). The labels stay exact.
- **Images.** All studies render per-voxel images, with at most 3 distractor sources, not full-sample renders.

**Conventions.**

- "Wrong" means an error > 1° (cubic-reduced).
- Intervals are 95% Wilson; paired tests are exact McNemar.
- Section names refer to `icenine_py/MIGRATION_HISTORY.md` (MH); benchmark paths are under `icenine_py/benchmarks/`.

# 2. Findings

## (a) Toy orientation network (learned Gauss-Newton layer)

**Question.** Can a network predict a voxel's orientation offset from its windows in one pass? Does building the
physics into the network help?

**Answer.** Yes.

- **The network.** GNLayerNet has a shared per-peak encoder that gives weights and measurement corrections, followed
  by pooled Gauss-Newton normal equations with a closed-form solve.
- **Clean multi-voxel data.** It removes the 3x gap of a pooled set network. On seen voxels it matches plain GN's
  accuracy; on unseen voxels it is 1.1–1.7x GN's error.
- **Realistic data.**
  - Plain GN degrades about 20x: median 0.22–0.33°, against 0.012–0.014° on clean data.
  - Huber-robust GN removes about a third of that.
  - A network trained on realistic windows reaches 0.06–0.08°, 2–3x better than robust GN.
  - A network trained only on clean data is barely better than plain GN.
- **Pixel noise alone.** Robust GN is as good as the network.
- **Detector pairing** changes little.

**Details.** MH "Toy Orientation NN" sections (2026-09-28 to 2026-10-01); `docs/orientation_nn_design.md`.

## (b) Single-voxel perturbation sweep: network vs optimizers vs FindOptimal

**Question.** Starting from the truth rotated by r°, how does each method's error grow?

**Answer.**

- **Clean data.** Iterating the network (3 passes, re-centred) holds the error at its floor of about 0.012° for every
  r up to 5°.
- **Realistic data.** The realism-trained network stays at 0.06–0.08° up to r = 3. It was trained in a 1° ball, so
  r > 1 is extrapolation.
- **FindOptimal** needs no guess of r, but from far starts it ends in the wrong basin.
- **The classic local optimizers do worse near the truth.**
  - MC, told r through its search box, ends at about 0.55 r.
  - Adam on the differentiable cost ends at about 0.13°. Its gradient is zero inside binary blobs.
  - Centroid Gauss-Newton on the network's windows reaches 0.0127° on clean data.

| r (deg) | FindOptimal median, realistic | FindOptimal fraction < 0.1° |
|---|---|---|
| 0.05 | 0.0189 | 0.990 |
| 1 | 0.0905 | 0.533 |
| 2 | 0.869 | 0.172 |
| 3 | 2.98 | 0.012 |

**Details.** MH "Toy Orientation NN — Perturbation sweep (2026-10-04)" and "Comparison with existing optimizers";
`benchmarks/toy_orientation_sweep/`.

## (c) Why full reconstruction returns wrong answers

**Question.** With no start, `reconstruct_voxel` is wrong for about a third of voxels. Where is the truth lost?

**Answer.** The search loses it, not the cost function.

- **Rates.** Over 200 voxels × 3 seeds, 204/600 = 34.0% of clean runs are wrong (CI 30.3–37.9), and 143/600 = 23.8%
  of realistic runs (20.6–27.4).
- **The cost would have chosen right.** In 347/347 wrong runs the cost at the truth is below the cost of the returned
  answer.
- **Where it happens.** The truth basin is among the level-0 candidates in 98–99% of runs. In 88% (clean) and 83%
  (realistic) of wrong runs it is lost at pruning: present, but not in the best quarter by post-quick-MC cost.
- **Why the ranking fails.** After quick MC a basin candidate is still 1–2° off, and the cost is sharp to 0.2–0.5°.
  Its cost is therefore 0.94–0.99, and it loses a ranking among near-1 costs.
- **What the wrong answers are.** Mostly near-degenerate alternatives that explain many of the same spots. In cubic
  Cu these are CSL relatives: Σ3 twins (62–66%), then Σ5 and Σ7.

**Fixes** (wrong rate clean / realistic):

| Fix | Wrong rate | Notes |
|---|---|---|
| Baseline | 33.5% / 26.0% | seed 0 |
| F1: CSL-relative check after the search, cubic-specific | 3.5% / 3.3% | 3 seeds; +11% evaluations; no right runs broken |
| F1b: relatives at every level | 0% / 1.5% | +26–28% evaluations |
| F2 (keep more candidates), F3 (full Q_max at coarse levels) | less improvement | |
| Smaller coarse pixel tolerance | worse | |

**Details.** MH "FindOptimal robustness — wrong candidates and CSL traps (2026-10-06)";
`docs/findoptimal_robustness_report.md`; `benchmarks/findoptimal_robustness/`.

## (d) Hybrid network → FindOptimal (H0 / H3 / HG)

**Question.** Does handing a network estimate to FindOptimal beat either one alone?

**Answer.** It depends on how far the start is from the truth.

- **The pipelines.**
  - H0: FindOptimal alone, from the nominal.
  - H3: 3 network passes, then FindOptimal.
  - HG: Huber GN, then 3 network passes, then FindOptimal.
- **Wrong rates on realistic data:**

  | Start r | H0 | net ×3 alone | H3 | HG |
  |---|---|---|---|---|
  | 3° | 94.6% | 6.7% | 4.1% | 0.0% |
  | 5° | 100% | 52.9% | 43.6% | 2.31% |

- **Precision.** H3's median error is about 0.020° for r ≤ 3, against 0.062–0.080° for the network alone. On clean
  data FindOptimal adds nothing after the network (0.011°).
- **Dropped variants.**
  - A covariance-sized FindOptimal box.
  - A network → MC finisher.
  - One network pass (H1) equals three only up to r = 1.5.

**Details.** MH "Task 1 results"; `benchmarks/nn_hybrid/`.

## (e) Q_max-8 cost proxy and rerank

**Question.** Can a cheap low-Q model replace the raw cost as the pruning key?

**Answer.** Yes, for accuracy.

- **The model.** A gradient-boosted model on low-Q (|q| ≤ 5) features, plus the Q8 cost, which is already computed
  at that point.
- **Pruning recall** (keep 1/4, levels 0–2, pooled over held-out grains):

  | Pruning key | Recall |
  |---|---|
  | Raw cost | 0.900 |
  | Untrained hit-rate key | 0.953 |
  | This proxy | 0.977 |
  | Full E2 classifier | 0.985 |

- **End to end** (200 voxels × 2 variants), the proxy rerank cuts the wrong rate from 33.5% / 26% to 6.5% / 7.5%.
  The E2 rerank gives 4.0% / 4.5%. The difference between them is not significant (McNemar p = 0.125 / 0.238).
- **The rerank is symmetry-agnostic.** Combined with the cubic-specific F1, it is wrong in 0% / 1.5% of runs over 3
  seeds.
- **How it runs.** The low-Q pass now runs batched; see `docs/batched_lowq_pass.md`.

**Details.** MH "Task 2 results"; `benchmarks/coarse_proxy/`.

## (f) Run-time profiling (U0 / U1 / U2) and re-timing

**Question.** Does the network or the proxy reduce total run time?

**Answer.** Only when a start orientation is 1.5° or more off.

- **No start (U0, single worker).** Every alternative to the baseline costs between +5% and +39% time. The reranks are
  accuracy gains, not speed-ups.
- **Start within 0.1° (U1).** Use H0. H1 and H3 are not faster.
- **Start 0.5–3° off (U2).** The hybrids help only at r ≥ 1.5°.
  - Pooled over r = 0.5–3, H3 is 3.9x (clean) and 3.3x (realistic) faster than H0 with an ideal fallback.
  - At r = 5 only HG helps.
- **Where the time goes.**
  - 98.3% of a baseline run is in `VoxelCostFunction.evaluate`.
  - The network costs 0.15–0.38 s of a 1.4–2.1 s hybrid case.
  - MPS does not matter per case.
- **Contention.** 10 workers vs 1 gives a factor of 1.25 (U0) and 1.21 (U1/U2).
- **Re-timing.** 32 tasks that had started under CPU contention were re-timed. No verdict changed.

**Details.** MH "Task 3 results" and its re-timing note; `benchmarks/profiling/`.

## (g) Follow-ups T1–T5

| Task | Result |
|---|---|
| T1 locked environment | A fresh Python 3.12 install broke the golden bit-identity tests. The cause is float-level drift (6e-9 in matrix entries), not a different search path. `uv sync --extra dev` from `uv.lock` (Python 3.9) is the documented setup. On a version mismatch the golden tests fail rather than skip. |
| T2 shared statistics | The study scripts use `scripts/common/stats.py`. Regenerated summaries are byte-identical except four Wilson upper bounds (1.0000000000000002 → 1.0). |
| T3 batched low-Q pass | Exact feature equality. 11.6x faster per candidate at batch 200 (0.26 evaluation equivalents). |
| T4 keep 1/6, 1/8 | Pooled wrong rate 41/400 (0.102) and 66/400 (0.165), against 28/400 (0.070) at keep 1/4 (McNemar p = 0.011 and 2.6e-7), for 2.5% and 3.3% of single-worker time saved. Not adopted. |
| T5 finisher diagnosis | See below. |

**T5: what the finisher does at the end.** Re-running the finisher on realistic sweep cases (200 H3, 100 H0)
reproduces the earlier results exactly.

- **It stops above the truth's cost.** The result's cost is higher than the truth's in 191/200 H3 cases (96%, CI
  92–98) and in 100/100 H0 cases.
- **How far.** For H3 the result is a median 0.023° from the truth, with a cost gap of 0.028, whatever the network's
  starting error.
- **No barrier in the way.** On the straight path from the result to the truth, a cost rise larger than the gap
  occurs in only 11/200 cases.
- **The cost is rough at this scale.** A 0.01° rotation at the truth changes the cost by a median 0.025, about the
  size of the gap.
- **The result is not a local minimum** on the 0.001–0.05° scale in 175/200 (H3) and 80/100 (H0) sampled cases.
- **How the MC stops.** The default finisher MC stops on its 200-step budget, or on exhausted restarts in exactly the
  runs that never accept a move.
- **Continuations from the result lower the cost.**
  - A step/10 MC pass improves 194/200 H3 cases.
  - A quarter-box VarianceMinimizing pass closes a median 93% of the gap. It ends 0.006° from the truth, but hit its
    50,000-step cap in 198/200 runs (about 19x the default evaluations).
- **The H0 failures at 2–3°** are a different problem: a wrong basin.

The effect of these continuations on reconstruction success has not been tested.

# 3. Recommendations

| Situation | Use | Reason |
|---|---|---|
| No start (U0) | Proxy or E2 rerank | Wrong rate 33.5% / 26% → 6.5% / 7.5% (proxy) or 4.0% / 4.5% (E2), at +5–7% time. An accuracy gain, not a speed-up. |
| Start within 0.1° | H0 | H1 and H3 are not faster. |
| Start possibly ≥ 1.5° off | HG (or H3 for r ≤ 3) | H0 is 15.9% wrong at 1.5° and 46.8% at 2° (realistic). At r = 5 only HG works. |
| Speed work | `VoxelCostFunction.evaluate` | 98.3% of a baseline run. |
| Proxy runs | `--batched` | Identical features, about 12x cheaper scoring at batch 200. |

# 4. Decisions on record

- **F1/F1b are off by default**, for parity with the C++ reconstruction. Both are opt-in, and both are cubic-specific.
- **Keep 1/6 and 1/8 are not adopted** (pooled verdict). The keep-1/4 proxy rerank stays, for accuracy.
- **The Wilson clamp is accepted.**
- **The environment is locked:** `uv sync --extra dev` from `uv.lock`. The golden tests fail rather than skip under
  other versions.

# 5. Caveats

- **Images:** per-voxel images with at most 3 distractor sources, not full-sample renders.
- **Material:** Cu FCC only.
- **Timing:**
  - Speed claims rest on single-worker timing (one thread, quiet machine, interleaved arms).
  - 10-worker wall times are contended and support no speed claim.
  - "H0 + fallback" assumes an ideal oracle trigger.
- **Seeds:** most wrong-rate comparisons use one seed (200 voxels per variant).

# 6. Open questions and next steps

- **Why does the finisher stop short of the truth's cost?** In particular, how much of this is the cost function's
  sensitivity: a binary overlap that is flat or rough at 0.01–0.02°, and a realistic-data minimum that is not at the
  truth in 59/200 H3 cases?
  - *Phase A (finisher/MC study, `benchmarks/cost_sensitivity/`):* no evidence that cost-function sensitivity limits the finisher at its 0.023° error scale. The landscape still has clear downhill directions there: the within-case Spearman correlation of cost with angle has median 0.97, 99.1% or more of the rays are monotone from 0.005° outward (84.8% of realistic H3 rays through all radii), a 0.01° rotation changes the cost by a median 0.025 (realistic), and the cost is flat only below about 0.0005–0.001°. The offset of the realistic minimum is not resolved below the 0.0005° sampling grid (q75 0.0020°). A quarter-box VarianceMinimizing pass reaches 0.0060° (at ~19× the default evaluations; step cap hit in 198/200), in 83% of cases. The finisher's error is a median 1.86× the centroid-quantisation scale (0.0125°; a scale, not a lower bound). So the optimizer stops early on a landscape resolvable at that scale; why it stops is Phase B. See MIGRATION_HISTORY "Phase A results"
- **Why does MC run out of restarts, and is MC the right local optimizer at all?** The April HP sweep compared methods
  by success from 1–5° starts, at 3500 MC steps. It did not measure final precision, nor FindOptimal's deployed
  configuration: 200 steps, a 0.33° box.
  - *Phase B1 (`benchmarks/mc_mechanism/`):* the restart rule holds as written. No restart fires after a run's first
    improvement (0/122); every run that exhausts its restarts never improved, and every never-improving run exhausts
    them (78/78 both ways in H3 realistic; their start is 0.0232° from the truth against a 0.1317° step); 5/122
    improving runs restarted before their first improvement.
    One improvement raises `min_ergodic` from 31 to 250, above the 200-step budget. Each improvement also halves the
    step, and an improvement is cheap at any step up to the distance to the truth (probability 0.27–0.47). The
    remaining travel falls below the distance to the truth at some improvement in 93% of the improving runs (77/122
    by the third), and the
    MC output is a median 0.0705° from the truth with probability 0 of an improving proposal at its own step in 89%
    (0.28 at the best grid step, a median 0.0075°). The step that maximises expected cost progress is about the
    distance to the truth. *Consistent with* this being why MC stops short; a step rule that follows the acceptance
    rate is the hypothesis B3 tests.
  - *Phase B2 (`benchmarks/sweep_audit/`):* the sweep's recorded 96% is "final error below 1.0°" from a 1° start, not
    0.5°. Its successes end at a median 0.47° (Adam) and 0.53° (MC); 0/100 end under 0.02°. MC at 100 steps equals MC
    at 3500 in that protocol, there is no detectable r_perp dependence (Spearman, n = 100 voxels) inside 75–563 µm, and the hybrid benchmark shows no
    difference from MC (9 runs only hybrid, 8 only MC, McNemar p = 1.0). So the comparison neither tested nor
    contradicts the deployed configuration at the 0.01–0.03° scale.
  - *Phase B3 (`benchmarks/finisher_bench/`):* a head-to-head of local finishers on the same cost function with counted
    evaluations (300 T5 cases and 1000 sweep cases, clean and realistic). MC is not the right tool for the finishing
    stage. At 250 evaluations CMA-ES (sigma0 0.2°) ends a median 0.0031° (T5 H3) and 0.0038° (sweep) from the truth, 94% and
    92% under 0.02°; the deployed MC at its 201 evaluations ends 0.0382° and 0.1609° (26% and 11% under 0.02°), and the
    whole default finisher 0.0229° and 0.0357° at 2629 and 1776 evaluations (44% and 30%). Nelder-Mead on the rotation
    vector does nearly as well at 250 evaluations (0.0031° and 0.0034°). The deployed MC ends early exactly when it
    never improved (74/200 T5 H3 realistic, 37.0%), and in those and in the improving runs the simplex and CMA methods
    end 0.0021°–0.0028° away at 1000 evaluations. Neither MC with local restarts (0.0362°) nor the success-rate rule (0.0240°)
    closes the gap: the success-rate rule as run (one untuned setting) collapses its step and stalls, so it does not
    support the B1 step-rule hypothesis; MC with local restarts keeps the halving and so tests only the restart part.
    The CMA and Nelder-Mead results reach the truth's cost (median gap 0.0000), so the remaining 0.002–0.003° is where
    the cost itself has its minimum, not early stopping. Starts 2–3° away often stay wrong (T5 H0: start 33/100 wrong;
    CMA-ES 15/100). Centroid Huber GN is good on clean windows (0.0132°) and poor on realistic (0.2628°);
    the hybrid Adam hardly moves. Single-worker time per cost evaluation is 0.55 ms for every method, so time follows
    evaluations (CMA-ES at 250: 0.143 s; default finisher 1.071 s). Phase C candidates: CMA-ES (sigma0 0.2°) and
    Nelder-Mead. See MIGRATION_HISTORY "Phase B3 results".
  - *Phase C (`benchmarks/cma_finisher/`):* CMA-ES (sigma0 0.2°, 1000 evaluations) is now an opt-in local refinement in
    the library (`SearchParameters.local_optimizer = "cma"`, config key `LocalOptimizer cma`), replacing the refinement MC
    in `refine_from_candidates` (seed voxels) and `local_optimization` (BFS neighbours); the default is bit-identical. It
    reproduces B3's `cma_02` exactly on 50 T5 H3 realistic cases (median 0.0025° at 1000 evaluations, 48/50 under
    0.02°). On the 200-voxel E0 set (no start, per-voxel images, seed 0) the right answers end a median 0.0019° (clean)
    and 0.0021° (realistic) from the truth against 0.0301° and 0.0278° under MC; wrong (> 1°) counts are 66/200 against
    67/200 (clean) and 44/200 against 52/200 (realistic, 8 fixed and 0 broken, exact McNemar p = 0.0078), at +0.19% and
    +3.5% mean cost evaluations and +1.7% and +4.8% single-worker time per no-start `reconstruct_voxel` (20 voxels per variant). A BFS neighbour
    (`local_optimization`) costs 1001 evaluations under CMA against a median 662 under MC (about +51%); MC returned the inherited
    start unchanged in 47/50 seeded cases (5° box, strictly-lower-cost acceptance). Not shown:
    the cost of 1001 evaluations per BFS neighbour at scale (the MC call used a median 662 in a 50-case check),
    behaviour for neighbours across a grain boundary, and seeds 1-2. To switch it on for a BFS run set the key (or
    `search_params.local_optimizer = "cma"`). See MIGRATION_HISTORY "Phase C results".
- **A full-sample, end-to-end BFS reconstruction of the 500-grain sample with new orientations,** comparing classic
  BFS (C++-parity optimizers) with BFS using the network and hybrid finisher, on timing and accuracy.
- **Other open items:**
  - seeds 1–2 for the keep rows;
  - H0 wrong-basin failures;
  - a real fallback trigger;
  - candidate de-duplication;
  - the three ideas in `docs/todo_future_ideas_nn_active_fourier.md`.

The plan for the first three items is the MIGRATION_HISTORY section "Finisher and MC study (2026-10-07)".

# 7. Symmetry framing

The failure mechanism in (c) is general. Ranking among near-degenerate alternatives loses the right basin whenever
other orientations explain much of the same spots.

- **Cubic crystals.** The alternatives are the CSL (Σ) relatives, mostly Σ3 twins.
- **Other symmetries.** Other near-degenerate orientations play that role.

The symmetry-agnostic results are the rerank, the hybrids, the finisher diagnosis and the profiling. F1, F1b and
proxy + F1 are cubic-specific: they enumerate Σ relatives.
