---
title: "TODO: hybrid network + optimizer refinement"
subtitle: "Raised after the perturbation sweep vs existing optimizers (2026-10-05); deferred until that work is finished"
date: "2026-10-05"
geometry: margin=1in
fontsize: 11pt
---

# Status

**In progress (2026-10-06): implemented as Task 1 of `feature/nn-hybrid-proxy-profiling`** (`scripts/nn_hybrid/`, results in `benchmarks/nn_hybrid/`, see `MIGRATION_HISTORY.md`, "Hybrid NN finisher, Q_max-8 cost proxy, run-time profiling"). Originally raised after the perturbation-sweep work (`MIGRATION_HISTORY.md`, "Toy Orientation NN — Perturbation sweep (2026-10-04)"; `benchmarks/toy_orientation_sweep/`). **Task 1 done (2026-10-06)**: net x3 -> seeded FindOptimal reaches a median of about 0.020 deg for r <= 3 (net alone 0.062-0.080, FindOptimal alone 0.037-2.98), with a net-failure tail at r = 3 (4.1% wrong); Huber GN -> net -> FindOptimal removes it (0.00% wrong for r <= 3, 2.3% at r = 5); the covariance-sized box and the net -> MC finisher did not help. Results: `MIGRATION_HISTORY.md`, "Task 1 results". Open: run-time cost on one worker (Task 3) and denser realism (full-sample renders).

# Why

On the same 1000 cases per radius (50 voxels of the 500-grain sample, realistic data), the
methods fail in complementary places:

- The realism-trained network (iterated x3) is the best method from about r = 0.75° to 3°
  (median 0.061–0.079°; FindOptimal 0.074° at 0.75°, 0.22° at 1.5°, 2.98° at 3°; it wins 58–95%
  of paired cases there). It has a floor near 0.06° set by the distractors: started closer than
  that (r = 0.05°) it makes the estimate worse in more than half of the cases, and it collapses
  at r = 5° (median 1.97°).
- Near the truth the existing optimizers are better on realistic data. FindOptimal (the seeded
  final stage, `refine_from_candidates`) has median 0.019–0.055° for r ≤ 0.5° against the net's
  0.063–0.065°, and the net wins only 24–49% of paired cases there. MC (told r through its
  search box) is best at r = 0.05° (0.028°), ties the net at r = 0.1°, and then stops at about
  0.5–0.6 r.
- FindOptimal's search box is fixed (0.33°, independent of r), so it cannot undo more than about
  1–1.5° of error (median 0.22° at r = 1.5°, 0.87° at 2°, 2.98° at 3°).
- Huber Gauss–Newton is robust at large r (0.23° at r = 5°) but about 3x worse than the network
  inside 3°.

# Idea

Use the network to bring the estimate to within about 0.5–0.75° of the truth (where it is
reliable and FindOptimal's box can reach the truth), then hand off to a finisher that is good
close to the truth:

1. network (iterated) → seeded FindOptimal (`AdaptiveVoxelReconstructor.refine_from_candidates`
   with the network estimate as the single candidate; FindOptimal + VarianceMinimizing + final
   evaluation) as the finisher, which reaches about 0.02° near the truth; MC with a small search
   box (≈ 0.1–0.2°) is the cheaper alternative;
2. optionally Huber GN first for very large starts (r ≳ 3°), then the network, then the finisher;
3. use the network's predicted covariance to set the finisher's search box per case
   (its Mahalanobis² is ≈ 3 on realistic data, i.e. calibrated).

# Questions to answer

- Does net → MC reach MC's near-truth accuracy (≈ 0.03°) from starts up to 3°, and at what cost
  per voxel vs MC alone and vs the full multi-level reconstruction?
- Is the covariance-sized box better than a fixed box?
- Does the gain survive denser realism (full-sample renders with every grain lighting the
  detector, rather than ≤ 3 distractor sources)?
