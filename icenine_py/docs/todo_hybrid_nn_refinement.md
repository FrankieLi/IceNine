---
title: "TODO: hybrid network + optimizer refinement"
subtitle: "Raised after the perturbation sweep vs existing optimizers (2026-10-05); deferred until that work is finished"
date: "2026-10-05"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Deferred, not started.** The project owner asked to revisit this after the perturbation-sweep
work (single-voxel sweep, comparison with MC/Adam/GN and the multi-level FindOptimal
reconstruction) is finished. Results it builds on: `MIGRATION_HISTORY.md`, "Toy Orientation NN
— Perturbation sweep (2026-10-04)" and its "Comparison with existing optimizers" subsection;
`benchmarks/toy_orientation_sweep/`.

# Why

On the same 1000 cases per radius (50 voxels of the 500-grain sample, realistic data), the
methods fail in complementary places:

- The realism-trained network (iterated x3) is the best method from r ≈ 0.1° to 3°
  (median ≈ 0.065–0.08°), but it has a floor near 0.06° set by the distractors: started closer
  than that (r = 0.05°) it makes the estimate worse in more than half of the cases, and it
  collapses at r = 5° (median 1.97°).
- MC (told r through its search box) is best very close to the truth (0.028° at r = 0.05°) but
  stops at about 0.5–0.6 r for larger starts.
- Huber Gauss–Newton is robust at large r (0.23° at r = 5°) but about 3x worse than the network
  inside 3°.

# Idea

Use the network to bring the estimate within about 0.1° (where it is reliable), then hand off to
a finisher that is good close to the truth:

1. network (iterated) → MC with a small search box (≈ 0.1–0.2°) or FindOptimal/VarianceMinimizing
   seeded at the network estimate;
2. optionally Huber GN first for very large starts (r ≳ 3°), then the network, then the finisher;
3. use the network's predicted covariance to set the finisher's search box per case
   (its Mahalanobis² is ≈ 3 on realistic data, i.e. calibrated).

# Questions to answer

- Does net → MC reach MC's near-truth accuracy (≈ 0.03°) from starts up to 3°, and at what cost
  per voxel vs MC alone and vs the full multi-level reconstruction?
- Is the covariance-sized box better than a fixed box?
- Does the gain survive denser realism (full-sample renders with every grain lighting the
  detector, rather than ≤ 3 distractor sources)?
