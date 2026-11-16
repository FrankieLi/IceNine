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
    recovered? Reference to look up (owner, 2026-11-16): HEXRD far-field HEDM strain methods,
    https://hexrd.github.io (local repo `~/Research/github/hexrd`, `fitgrains.py`). Far field
    puts the detector far away, so a grain covers only a few pixels; strain resolution depends
    only on the ratio of detector pixel size to sample-detector distance.
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
- **Quality**: per-voxel symmetry-reduced misorientation to the voxel's own ground truth: median, 95th
  percentile, share above 0.1 deg, wrong (above 1 deg), with Wilson 95% intervals. Boundary
  voxels (a neighbour outside the group) reported separately. The strict binary cost
  (`VoxelCostFunction`, local, hard) is the referee and is reported per voxel.
- **Cost**: objective evaluations per grain (one evaluation renders every voxel of the group
  once) and wall seconds (labelled single-worker or contended by `preflight_label`). For a
  per-voxel method, strict evaluations per voxel. The Pareto x axis is "objective evaluations
  per voxel": one evaluation of a grain objective scores every voxel of the grain once (counts 1
  per voxel); a per-voxel fit counts its strict-cost calls per voxel; chained stages add up.
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

Question: where does each objective's minimum sit relative to the ground truth? Level-0 grains of
`grad` (one true orientation each): grains 381 and 302. CMA over the grain rotation only
(3 numbers), started AT the ground truth, 150 evaluations, sigma 0.03 deg. Arms: unbalanced Sinkhorn
(reach 3 px) at blur 1 and 0.5 px, balanced at blur 1, with the pixel-box dilation of the
prediction at blur 1 and 0.5, intensity-weighted, least squares. Reported: median offset of the
optimum from the ground truth (deg), objective and mean strict cost at ground truth and optimum. Rule: an
objective is "unbiased enough" if its offset is below 0.01 deg on both grains (the start's median
error on level-0 grains is 0.016 deg, on all painted voxels 0.041 deg). If no arm qualifies, the
next step studies the cause (I3, I4) before any optimizer study.

### GO2: optimizer Pareto on G20 (ideas I5, I6, I9, I10, I28, I29)

Question: for each objective, how many objective evaluations per voxel buy how much accuracy, and which chain
of stages reaches the strict-cost optimum cheapest? GO1 found every Sinkhorn variant flat within
0.03 to 0.04 deg of the ground truth, and the strict cost (grain mean) within 0.006 to 0.008 deg. So
Sinkhorn is a basin finder at best, and the finish must use the strict cost (or a sharper smooth
cost, I3).

Benchmark G20 from the painted start, start-map flood-fill groups. Arms (pipelines of
`run_pipeline.py`; A = grain-field CMA, B = per-voxel Adam with smoothness, V = per-voxel CMA on
the strict cost):

1. Field CMA on the strict cost: `Astrict/50`, `/100`, `/200`, `/400` (I28).
2. The same with the quadratic field: `Astrict/200/o2`.
3. Field then per-voxel polish: `Astrict/100,V/10`, `Astrict/100,V/30`, `Astrict/100,V/61` (I10).
4. Per-voxel only, from the start: `V/61` (AR3 C re-run on the same voxels, I9) and `V/20`.
5. Sinkhorn first: `Asink/150,Astrict/100,V/20`; and `Asink/150,Bsink/40,Astrict/50` (I6).
6. Least squares first: `Als/150,Astrict/100,V/20` (I21 as basin finder).

Judging: the standing rule. Pareto axes: objective evaluations per voxel (x) against share above 0.1 deg and
median error (y), on the voxels of the 20 grains with a C record; also wall seconds per voxel.
The frontier is the set of arms no other arm beats on both axes (with interval overlap counted as
a tie). Budget: strict-cost arms take seconds per grain; Sinkhorn arms about 1 min per grain;
whole step under 1 h of compute. Committed: run.json and voxels.json.gz per arm (each under
1 MB), tables.md, Pareto rows.

### GO3: deformed samples (ideas I11, I12) -- built

`rand` (white-noise voxel rotations) and `corr3` (Gaussian-smoothed, length 3 voxel sides), RMS
levels 0.5, 1, 2, 5 deg drawn per grain (seed 0). Rendered clean with the stock forward model;
coverage mass on lit pixels at the ground truth 97 to 99% on 10 grains of each. In `rand`, 69% (2 deg)
and 97% (5 deg) of neighbour pairs exceed 2 deg, so flood-fill grouping fragments those grains:
GO4 uses true grain ids (`--group true`) and the grain-centre start (`--start centre`).

### GO4: robustness across kinds of orientation variation (ideas I13, I6)

Question: does the GO2 frontier (grain field on the strict cost, then a per-voxel polish) hold
when the voxel-to-voxel variation is random or spatially correlated, from 0.5 to 5 deg RMS?
Samples `rand` (white noise) and `corr3` (correlation length 3 voxel sides), GO3. Start: each
voxel at its grain's orientation before deformation (`--start centre`, a perfect paint that knows
nothing of the variation); membership: true grain ids (`--group true`; flood fill fails at 2 and
5 deg). Benchmark: 20 grains per sample by `pick_true_grains` (5 per level 0.5, 1, 2, 5 deg).

Arms (per sample): `Astrict/100`; `Astrict/100,V/61`; `Astrict/100,V/61/0.3` (per-voxel step
0.3 deg); `V/61/0.3`; `Als/150,Astrict/100,V/61/0.3` (least squares as basin finder).
Judging: the standing rule, reported per level (median, p95, above 0.1 deg, wrong above 1 deg,
Wilson intervals); no C reference exists on these samples, so all voxels of the groups count.
Expected: the field cannot follow white noise; the question is which polish recovers which level.
Budget: strict arms seconds per grain, the LS arm about 10 s per grain; under 1 h in total.

### GO5 onwards

Planned order (each gets a section when it starts): GO5 strain in the forward
model and the strained samples (I14 to I16); GO6 objective variants (I17 to I22); GO7 richer
simulation (I23, I24); GO8 learned models (I25 to I27).
