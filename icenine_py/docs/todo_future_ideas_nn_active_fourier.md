---
title: "TODO: single-voxel NN reconstruction, active imaging, and diffraction as Fourier-space tomography"
subtitle: "Three ideas raised by the project owner during the hybrid / proxy / profiling work (2026-10-06); not started"
date: "2026-10-06"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Ideas, not started.** Raised by the project owner on 2026-10-06, while the run-time profiling
of `feature/nn-hybrid-proxy-profiling` was running (`MIGRATION_HISTORY.md`, "Hybrid NN finisher,
Q_max-8 cost proxy, run-time profiling"). The notes below each idea are a first framing for
discussion, not a plan.

# Idea 1: single-voxel NN reconstruction straight from the detector images

**What.** Have a network predict a voxel's orientation from its detector images with **no
starting orientation**. This would replace or seed the coarse search, which is where the no-start
case spends its time and where the truth is lost (S2 pruning).

**Why it might help.**
- In the no-start case, `reconstruct_voxel` spends about 44.5k global evaluations on the coarse
  discrete search, against 5–8k local ones (about 25 s per voxel on 10 workers). Task 3 of the
  current feature measures this split on a single worker.
- The current net (`GNLayerNet`) is only a local refiner. It needs a nominal orientation to
  predict the ROI windows it reads, so it cannot start from nothing.

**Framing.**
- **Input.** Without an orientation, windows cannot be cut, so the input must be one of:
  - everything in the voxel's projected footprint: at each ω, the voxel position confines its
    spots to a known band on the detector;
  - a spot list, (ω, detector x, y, detector distance) per spot, fed to a permutation-invariant
    set network.
- **Output.**
  - Cubic symmetry is 24-fold: predict in the fundamental zone or use a symmetry-invariant loss.
  - The answer can be multi-modal (Σ3 twins and other CSL relatives look alike), so output the
    top-k orientations or a mixture rather than a single point.
- **Use.** As a candidate generator: hand its top-k to `refine_from_candidates` (FindOptimal),
  optionally with the F1 CSL-relative check. Using the network only to propose candidates, and
  keeping the physics-based cost as the final judge, is also what made Tasks 1 and 2 work.

**Experiments (sketch).**
1. Top-1 and top-k basin hit rates (within 1° and 3°, cubic-reduced) on the 200-voxel E0 set,
   clean and realistic, with whole grains held out of training.
2. End to end: net top-k → FindOptimal (→ F1). Measure the wrong rate with Wilson CI, the median
   error and the wall time, against the baseline, F1, F1b and proxy + F1.
3. Cost: net time vs the coarse search time from Task 3. It only helps if it is much cheaper than
   about 44k evaluations at the same wrong rate.

**Risks.** Generalising across structures, detector geometries and realism (overlap and noise).
Training data comes from the Python forward model, which is validated against C++.

# Idea 2: active imaging, reconstructing from incomplete data and refining as images arrive

**What.**
- Reconstruct from a subset of the images: some ω steps, or one detector distance.
- Keep the candidate list, not just the best orientation.
- Add images progressively, re-scoring and pruning the kept candidates with each new image.
- Choose the next images to acquire by where the remaining candidates disagree most.

**Why it might help.**
- It saves beam time and enables in-situ or time-resolved experiments, where full ω scans are
  costly.
- The FindOptimal study showed that the hard cases are the truth versus a few CSL relatives. Their
  predicted spots differ on specific reflections and ω ranges, so a well-chosen image could settle
  a case that many uninformative images would not.
- The reconstructor already holds a candidate list per level (recorder, `extra_candidates`,
  `rank_key` hooks).

**Framing.**
- **Cost on a subset of images.** Restrict `VoxelCostFunction` to an ω subset or detector subset.
  Check how this interacts with the per-level `n_q_max` schedule.
- **Candidate-list persistence.** Keep the top-N distinct candidates, deduplicated as in
  `todo_coarse_cost_proxy_and_dedup.md` Idea 2, between rounds, rather than collapsing to one
  answer.
- **Acquisition rule.** Choose the ω (and detector distance) that maximise the expected
  disagreement between the remaining candidates' predicted spots: an information-gain or
  query-by-committee rule.
- **Stopping rule.** Stop once the cost margin between the best candidate and the runner-up, or a
  posterior margin, exceeds a threshold.

**Experiments (sketch).**
1. Wrong rate and error against the fraction of images used (10%, 25%, 50%, 100%), with the
   images chosen randomly, uniformly in ω, or adaptively.
2. How fast the candidate list shrinks, and whether the truth survives, compared with the current
   one-shot pruning.
3. Whole sample: BFS neighbours inherit candidate lists. When does a voxel's answer stop changing?

# Idea 3: diffraction imaging as tomography in Fourier space, and what it means for resolution

**What.** Make precise the view that NF-HEDM is a form of tomography, and work out the resolution
consequences.

- **Real space.** Each diffraction spot is a projection of the diffracting volume along the
  diffracted beam. The spots across ω and detector distances form a set of projections, as in
  tomography.
- **Reciprocal (Fourier) space.** Each reflection samples the crystal's reciprocal lattice at one
  scattering vector **G**. The set of measured **G** (bounded by Q_max, the η limit and the ω
  range) is a sampling of Fourier space.

**Questions.**
1. **Spatial resolution.** How does it depend on the pixel size (1.48 µm here), the number of
   usable reflections per voxel, and the projection geometry?
   - Does incomplete ω coverage or the η limit leave a "missing wedge" that makes the resolution
     anisotropic, as in limited-angle tomography?
2. **Orientation resolution.** It should improve with |G|: larger |q| gives a larger spot
   displacement for a given rotation.
   - Relate this to the measured sharpness of the cost function (0.2–0.5°), the Q_max = 8 choice,
     and the finding that low-Q reflections alone cannot separate CSL relatives (Task 2).
3. **Coupling.** Position and orientation errors trade off: a small rotation and a small shift can
   move a spot by the same amount. Which combinations of detector distances and reflections remove
   that degeneracy?
4. **Sampling.** Is there a Nyquist-like rule for how many reflections and ω steps are needed per
   voxel for a given spatial and orientation resolution? It bears directly on Idea 2: how few images
   are enough.

**Check against DCT.**
- **Diffraction Contrast Tomography (DCT)**, confirmed by the project owner as the intended
  comparison: grain mapping from diffraction spots treated as projections, reconstructed with
  ART/SIRT-type algorithms (Ludwig et al., J. Appl. Cryst. 2008 and later).
- Secondary: **discrete tomography** (e.g. DART, Batenburg & Sijbers, IEEE Trans. Image Process.
  2011). Grain maps are piecewise constant, which is exactly the prior that discrete
  tomography uses to reconstruct from few projections. That prior may give the resolution versus
  number of images trade-off for Idea 2.
- Read both for their resolution analyses (spatial resolution vs pixel size and number of spots;
  orientation sensitivity) and compare with the adaptive forward-model approach (Li & Suter 2013).
- The references above are from memory and still to be checked.

# Relation to other work

- Idea 1 extends the hybrid NN study (`todo_hybrid_nn_refinement.md`) from "a start exists" to
  "no start".
- Idea 2 reuses the candidate-list machinery and the dedup idea
  (`todo_coarse_cost_proxy_and_dedup.md` Idea 2), and targets the CSL traps from
  `findoptimal_robustness_report.md`.
- Idea 3 is the theory behind Q_max, pixel size and image-count choices, and would set the limits
  that Ideas 1 and 2 are judged against.
