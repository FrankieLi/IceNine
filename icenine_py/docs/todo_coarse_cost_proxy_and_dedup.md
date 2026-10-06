---
title: "TODO: a high-Q_max cost proxy and cheap candidate de-duplication for the coarse search"
subtitle: "Two ideas raised after the FindOptimal robustness study (2026-10-06); not started"
date: "2026-10-06"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Ideas, not started.** Raised by the project owner after the FindOptimal robustness study
(`docs/findoptimal_robustness_report.md`, `MIGRATION_HISTORY.md` "FindOptimal robustness — wrong
candidates and CSL traps"). Both target the stage where that study found the truth is lost: the
coarse levels of `AdaptiveVoxelReconstructor.reconstruct_voxel`, where 83–88% of wrong runs had a
truth-basin candidate that was then pruned (S2).

# Idea 1: a proxy that predicts the Q_max = 8 score from a cheap low-Q_max evaluation

**What.** Train a model that maps a cheap cost evaluation at Q_max = 3, 4 or 5 (few reflections,
fast) to the score the candidate would get at Q_max = 8 (all reflections the final cost uses), and
rank/prune the coarse candidates by the predicted Q_max = 8 score.

**Why it might help.**
- The coarse levels start at n_q_max = 5 and add one per level; the final local cost uses
  Q_max = 8. A CSL relative shares a large subset of the low-|q| reflections with the truth (for
  Σ3 about 60–76% of all reflections up to Q_max 8 are shared), so at low Q_max a trap can look as
  good as the truth; the reflections that tell them apart are mostly at higher |q|.
- Simply evaluating at Q_max = 8 throughout (fix F3a) was expensive (+30–45% evaluations) and broke
  some right runs, so a cheap proxy is the interesting middle ground.
- The E2 classifier already showed that per-|q|-family hit rates from one pass carry the signal
  (AUC 0.997 vs 0.963 for the cost); a proxy restricted to low-|q| inputs tests whether that signal
  survives without evaluating the high-|q| reflections at all.

**Inputs to try.** Per-|q|-family hit rates (exact and ±1/±3 px) at Q_max 3/4/5, number of predicted
spots, the low-Q_max cost itself; optionally the candidate's orientation relative to the sample
(which reflections are near the detector edges or the eta limit).

**Target.** The Q_max = 8 local cost (regression), or directly "would this candidate survive
pruning / is it in the truth's basin" (classification). Regression keeps the reconstructor's
current ranking logic and is easier to sanity-check.

**Experiments.**
1. Offline: on the E0 candidate dataset (`scripts/findoptimal_robustness/e2_dataset.py`, about
   360k candidates with features and errors), fit low-Q_max → Q_max-8 score per Q_max in {3, 4, 5};
   report rank correlation with the true Q_max-8 cost, AUC basin vs trap, and pruning recall,
   against (a) the raw low-Q_max cost and (b) the current coarse ranking.
2. End to end: use the proxy as the pruning key via the existing `rank_key` hook; wrong rate with
   95% CI and evaluations per run vs baseline, F1 and F1b, on the same 200 voxels × seeds.
3. Cost accounting: time per proxy evaluation including the low-Q_max forward pass; the gain only
   counts if it is cheaper than evaluating at Q_max 8.

**Risks.** The mapping depends on structure (copper FCC here), Q_max, detector geometry and realism
(noise, overlap); the E2 classifier lost final-choice precision when trained on clean data and
tested on realistic data (0.83). Train with realistic data and hold out by grain.

# Idea 2: low-cost clustering of N candidate orientations to remove duplicates

**What.** Before pruning (and before FindOptimal), collapse candidates that are the same
orientation up to cubic symmetry (and a small tolerance) into one, keeping the best-scored
representative, so the kept slots go to distinct orientations.

**Why it might help.**
- Each level keeps only the top quarter (median 172 → 45 → 11 → 2 candidates on clean data), and
  FindOptimal receives 2–4 candidates. If several kept candidates are copies of the same trap
  (e.g. the same Σ3 relative reached from different grid points after quick MC), they crowd out the
  truth's basin candidate. Whether this happens is **not yet measured**.
- Duplicates also waste quick-MC and FindOptimal evaluations.

**Cheap method to try.**
- Map each candidate to the fundamental zone (unique representative under the 24 cubic operators),
  as a unit quaternion or Rodrigues vector; build a `scipy.spatial.cKDTree` (project rule: never
  `KDTree`) and do a radius query at the tolerance; handle the FZ boundary by also inserting the
  symmetric copies that fall within the tolerance of a boundary.
- Then greedy non-maximum suppression: sort by score, keep a candidate, drop its neighbours within
  the tolerance. O(N log N) instead of the existing O(N² · 24) helper
  (`scripts/findoptimal_robustness/csl.py::dedup_orientations`), which is the correctness reference.
- Tolerance: tie it to the level's resolution (local grid diameter / quick-MC box), e.g. 0.5–1×
  the box at that level.

**Experiments.**
1. Measure first: on the recorded E0 levels, how many kept candidates are within 0.5°/1°/2°
   (symmetry-reduced) of a better-ranked kept candidate, per level, in wrong vs right runs. If
   duplicates are rare among kept candidates, this idea is about cost, not robustness.
2. Correctness test: the fast clustering returns the same groups as `dedup_orientations` on random
   and adversarial sets (orientations near FZ boundaries, symmetry-equivalent copies).
3. End to end: dedup before each pruning step (keep count unchanged, so freed slots go to the next
   distinct candidates); wrong rate and evaluations vs baseline, alone and combined with F1/F1b.

# Relation to other work

- Complements F1/F1b (CSL-relative checks), which fix most failures today; these ideas aim to keep
  the truth from being pruned in the first place, possibly at lower cost than F1b (+28%).
- The `rank_key` hook and the instrumented E0 recorder in `icenine/reconstructor.py` and
  `scripts/findoptimal_robustness/` make both ideas testable without further reconstructor changes.
