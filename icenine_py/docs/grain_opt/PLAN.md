# AR6 (`grainopt`): full-grain fitting research program

Owner's request (Frankie, 2026-11-16), in the owner's words where possible. Branch
`feature/grain-opt`, worktree `.claude/worktrees/grain-opt`. Separate from Newt's autoresearch
(`docs/autoresearch/`): never edit those files or the main checkout.

## Goals

- **G1. Find an algorithm that reaches the optimum of the objectives we have.**
  - **G1a.** For a fixed objective, sketch the Pareto frontier of cost (evaluations) against
    quality, across algorithms and settings.
  - **G1b.** Robustness against the kind of gradient: from linear gradients (easy) to large,
    fairly random voxel-to-voxel misorientations (0 < d < 5 deg; plastically deformed samples
    with geometrically necessary dislocations).
  - **G1c.** Sensitivity to elastic strain (lattice deformation): can it be captured and
    recovered?
- **G2. Better objectives to guide the reconstruction.**
  - **G2a.** Different ways to use intensity.
  - **G2b.** More sophisticated diffraction simulation (rocking curve alpha > 0, partial-volume
    coverage, strain in the forward model, detector response).
  - **G2c.** More sophisticated models, including neural networks (reuse
    `docs/orientation_nn_design.md`, `docs/nn_inverse_problem_formulation.md`,
    `docs/research_ideas_nerf.md`, `docs/todo_intensity_and_distribution_losses.md`,
    `docs/todo_hybrid_nn_refinement.md`, and shmem).

## Standing judging rule (applies to every step unless its section overrides it)

Written before any run of the step.

- **Benchmark G20**: the 20 `grad` grains picked by `run_grains.pick_true_grains(T, 20)` (4 per
  gradient level 0, 0.1, 0.2, 0.5, 1.0 deg; sizes at quantiles 0.1 to 0.95 of the level). The
  fitted group of a true grain is the start-map flood-fill group (2 deg) holding most of its
  voxels. Start map: the painted map (AR3 step 3, `runs/P/maps.npz` "P").
- **Quality**: per-voxel symmetry-reduced misorientation to the voxel's own truth: median, 95th
  percentile, share above 0.1 deg, wrong (above 1 deg), with Wilson 95% intervals. Boundary
  voxels (a neighbour outside the group) reported separately. The strict binary cost
  (`VoxelCostFunction`, local, hard) is the referee and is reported per voxel.
- **Cost**: objective evaluations per grain (one evaluation renders every voxel of the group
  once) and wall seconds (labelled single-worker or contended by `preflight_label`). For a
  per-voxel method, strict evaluations per voxel. The Pareto x axis is "renders per voxel":
  grain-objective evaluations, or per-voxel strict evaluations.
- **Reference points**: start (0), AR3 step 4 CMA per-voxel fit `C` (61 strict evaluations per
  voxel; on the painted voxels 334 of 19,013 above 0.1 deg, 1.6 to 2.0%), on the same voxels.
- **Comparisons**: a gain or loss only when Wilson intervals do not overlap. No p-values, no
  extra seeds unless the step says so. One run per arm.
- **Acceptance inside a fit**: a grain update is kept only if the objective and the mean strict
  cost of its voxels do not get worse.
- **Stop rule**: a run stops at 2x its time budget; a step at 2 h of compute unless a decision
  entry raises it. Check `uptime` before runs (Newt's autoresearch shares the machine; at most
  4 heavy processes for this track).

## Steps

Each step gets its own section below before it runs: question, arms, judging rule overrides,
budget, what is committed. Ideas live in `IDEAS.md`; a step takes one or more ideas.

### GO1: objective bias (ideas I1, I2)

Question: where does each objective's minimum sit relative to the truth? Level-0 grains of
`grad` (one true orientation each): grains 381 and 302. CMA over the grain rotation only
(3 numbers), started AT the truth, 150 evaluations, sigma 0.03 deg. Arms: unbalanced Sinkhorn
(reach 3 px) at blur 1 and 0.5 px, balanced at blur 1, with the pixel-box dilation of the
prediction at blur 1 and 0.5, intensity-weighted, least squares. Reported: median offset of the
optimum from the truth (deg), objective and mean strict cost at truth and optimum. Rule: an
objective is "unbiased enough" if its offset is below 0.01 deg on both grains (the start's median
error on level-0 grains is 0.016 deg, on all painted voxels 0.041 deg). If no arm qualifies, the
next step studies the cause (I3, I4) before any optimizer study.

### GO2 onwards

Planned order (each gets a section when it starts): GO2 optimizer Pareto on G20 with the GO1
objective (I5 to I10); GO3 sample generator for random and correlated voxel-to-voxel
misorientation (I11, I12); GO4 robustness on those samples (I13); GO5 strain in the forward
model and the strained samples (I14 to I16); GO6 objective variants (I17 to I22); GO7 richer
simulation (I23, I24); GO8 learned models (I25 to I27).
